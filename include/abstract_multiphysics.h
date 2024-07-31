//
// Created by xlluo on 24-7-31.
//

#ifndef CRACKS_ABSTRACT_MULTIPHYSICS_H
#define CRACKS_ABSTRACT_MULTIPHYSICS_H

#include "abstract_field.h"
#include "controller.h"
#include "dealii_includes.h"
#include "elasticity.h"
#include "parameters.h"
#include "phase_field.h"
#include "utils.h"
#include <fstream>
#include <iostream>
using namespace dealii;

template <int dim> class AbstractMultiphysics {
public:
  explicit AbstractMultiphysics(Parameters::AllParameters &prms);

  void run();
  Controller<dim> ctl;

private:
  virtual void setup_system() { AssertThrow(false, ExcNotImplemented()); };
  virtual void refine_grid() { AssertThrow(false, ExcNotImplemented()); };
  virtual void record_old_solution() {
    AssertThrow(false, ExcNotImplemented());
  };
  virtual void return_old_solution() {
    AssertThrow(false, ExcNotImplemented());
  };
  virtual double staggered_scheme() {
    AssertThrow(false, ExcNotImplemented());
  };
  virtual void respective_output_results(DataOut<dim> &data_out) {
    AssertThrow(false, ExcNotImplemented());
  };

  void setup_mesh();
  void output_results();
};

template <int dim>
AbstractMultiphysics<dim>::AbstractMultiphysics(Parameters::AllParameters &prms)
    : ctl(prms) {}

template <int dim> void AbstractMultiphysics<dim>::run() {
  ctl.dcout << "Project: " << ctl.params.project_name << std::endl;
  ctl.dcout << "Mesh from: " << ctl.params.mesh_from << std::endl;
  ctl.dcout << "Load sequence from: " << ctl.params.load_sequence_from
            << std::endl;
  ctl.dcout << "Output directory: " << ctl.params.output_dir << std::endl;
  ctl.dcout << "Solving " << ctl.params.dim << " dimensional PFM problem"
            << std::endl;
  ctl.dcout << "Running on " << Utilities::MPI::n_mpi_processes(ctl.mpi_com)
            << " MPI rank(s)" << std::endl;
  ctl.dcout << "Number of threads " << MultithreadInfo::n_threads()
            << std::endl;
  ctl.dcout << "Number of cores " << MultithreadInfo::n_cores() << std::endl;

  ctl.timer.enter_subsection("Set mesh");
  setup_mesh();
  ctl.timer.leave_subsection("Set mesh");

  ctl.timer.enter_subsection("Initialize system");
  setup_system();
  ctl.timer.leave_subsection("Initialize system");

  //  if (ctl.params.enable_phase_field) {
  //    enforce_phase_field_limitation();
  //  }

  ctl.timer.enter_subsection("Solve Newton system");
  unsigned int refinement_cycle = 0;
  double finishing_timestep_loop = 0;
  double tmp_timestep = 0.0;

  ctl.current_timestep = ctl.params.timestep;
  // Initialize old and old_old timestep sizes
  ctl.old_timestep = ctl.current_timestep;

  do {
    double newton_reduction = 1.0;
    if (ctl.timestep_number > ctl.params.switch_timestep &&
        ctl.params.switch_timestep > 0)
      ctl.current_timestep = ctl.params.timestep_size_2;

    double tmp_current_timestep = ctl.current_timestep;
    ctl.old_timestep = ctl.current_timestep;

  mesh_refine_checkpoint:
    ctl.pcout << std::endl;
    ctl.pcout << "\n=============================="
              << "=========================================" << std::endl;
    ctl.pcout << "Time " << ctl.timestep_number << ": " << ctl.time << " ("
              << ctl.current_timestep << ")" << "   "
              << "Cells: " << ctl.triangulation.n_global_active_cells();
    ctl.pcout << "\n--------------------------------"
              << "---------------------------------------" << std::endl;
    ctl.pcout << std::endl;

    ctl.time += ctl.current_timestep;

    do {
      // The Newton method can either stagnate or the linear solver
      // might not converge. To not abort the program we catch the
      // exception and retry with a smaller step.
      //          use_old_timestep_pf = false;

      record_old_solution();
      try {
        newton_reduction = staggered_scheme();
        while (newton_reduction > ctl.params.upper_newton_rho) {
          //              use_old_timestep_pf = true;
          ctl.time -= ctl.current_timestep;
          ctl.current_timestep = ctl.current_timestep / 10.0;
          ctl.time += ctl.current_timestep;
          return_old_solution();
          staggered_scheme();

          if (ctl.current_timestep < 1.0e-9) {
            ctl.pcout << "Step size too small - keeping the step size"
                      << std::endl;
            break;
          }
        }

        break;

      } catch (SolverControl::NoConvergence &e) {
        ctl.pcout << "Solver did not converge! Adjusting time step."
                  << std::endl;
        ctl.time -= ctl.current_timestep;
        return_old_solution();
        ctl.current_timestep = ctl.current_timestep / 10.0;
        ctl.time += ctl.current_timestep;
      }
    } while (true);

    //    LA::MPI::Vector distributed_solution(elasticity.locally_owned_dofs,
    //    ctl.mpi_com); distributed_solution = elasolution;
    //    elasticity.distribute_hanging_node_constraints(distributed_solution,
    //    ctl);

    // Refine mesh and return to the beginning if mesh is changed.
    //    if (ctl.params.refine) {
    //      bool changed = refine_grid();
    //      if (changed) {
    //        // redo the current time step
    //        ctl.pcout << "Mesh changed! Re-do the current time step" <<
    //        std::endl; ctl.time -= ctl.current_timestep; solution =
    //        old_solution; goto mesh_refine_checkpoint; continue;
    //      }
    //    }

    // Recover time step
    ctl.current_timestep = tmp_current_timestep;
    ctl.timer.leave_subsection("Solve Newton system");
    ctl.timer.enter_subsection("Calculate outputs");
    ctl.computing_timer.enter_subsection("Calculate outputs");
    output_results();
    ctl.computing_timer.leave_subsection("Calculate outputs");
    ctl.timer.leave_subsection("Calculate outputs");
    ctl.timer.enter_subsection("Solve Newton system");
    ++ctl.timestep_number;

    ctl.computing_timer.print_summary();
    ctl.computing_timer.reset();
    ctl.pcout << std::endl;
    refine_grid();
  } while (ctl.timestep_number <= ctl.params.max_no_timesteps);
  ctl.timer.leave_subsection("Solve Newton system");
  ctl.timer.manual_print_summary(ctl.dcout.fout);
}

template <int dim> void AbstractMultiphysics<dim>::setup_mesh() {
  GridIn<dim> grid_in;
  /**
   * similar to normal use of GridIn.
   */
  grid_in.attach_triangulation(ctl.triangulation);
  if (!checkFileExsit(ctl.params.mesh_from)) {
    throw std::runtime_error("Mesh file does not exist");
  }
  std::filebuf fb;
  if (fb.open(ctl.params.mesh_from, std::ios::in)) {
    std::istream is(&fb);
    grid_in.read_abaqus(is);
    fb.close();
  }
  //  GridGenerator::hyper_cube(ctl.triangulation);
  //  ctl.triangulation.refine_global(5);

  if (dim == 2) {
    std::ofstream out(ctl.params.output_dir + "initial_grid.svg");
    GridOut grid_out;
    grid_out.write_svg(ctl.triangulation, out);
  }

  //  std::tuple<std::vector<Point<dim>>, std::vector<CellData<dim>>,
  //  SubCellData> info; info =
  //  GridTools::get_coarse_mesh_description(ctl.triangulation);
  ctl.dcout << "Find " << ctl.triangulation.n_global_active_cells()
            << " elements" << std::endl;
}

template <int dim> void AbstractMultiphysics<dim>::output_results() {
  DataOut<dim> data_out;
  data_out.attach_triangulation(ctl.triangulation);

  Vector<float> subdomain(ctl.triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = ctl.triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");

  // Record statistics
  ctl.statistics.add_value("Step", ctl.timestep_number);
  ctl.statistics.set_precision("Step", 1);
  ctl.statistics.set_scientific("Step", false);
  ctl.statistics.add_value("Time", ctl.time);
  ctl.statistics.set_precision("Time", 8);
  ctl.statistics.set_scientific("Time", true);

  respective_output_results(data_out);

  data_out.build_patches();

  data_out.write_vtu_with_pvtu_record(ctl.params.output_dir, "solution",
                                      ctl.timestep_number, ctl.mpi_com, 2, 8);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
    std::ofstream stat_file(
        (ctl.params.output_dir + "/statistics.txt").c_str());
    ctl.statistics.write_text(stat_file);
    stat_file.close();
  }
}

#endif // CRACKS_ABSTRACT_MULTIPHYSICS_H
