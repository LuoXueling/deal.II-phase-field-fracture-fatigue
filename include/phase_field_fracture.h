/**
 * Xueling Luo @ Shanghai Jiao Tong University, 2022
 * This code is for multiscale phase field fracture.
 **/

#ifndef PHASE_FIELD_FRACTURE_H
#define PHASE_FIELD_FRACTURE_H

#include "abstract_field.h"
#include "abstract_multiphysics.h"
#include "controller.h"
#include "dealii_includes.h"
#include "elasticity.h"
#include "parameters.h"
#include "utils.h"
#include <fstream>
#include <iostream>
using namespace dealii;

template <int dim> class PhaseFieldFracture : public AbstractMultiphysics<dim> {
public:
  explicit PhaseFieldFracture(Parameters::AllParameters &prms);

private:
  void setup_system() override;
  void refine_grid() override;
  void record_old_solution() override;
  void return_old_solution() override;
  double staggered_scheme() override;
  void respective_output_results(DataOut<dim> &data_out) override;

  Elasticity<dim> elasticity;
  PhaseField<dim> phasefield;
};

template <int dim>
PhaseFieldFracture<dim>::PhaseFieldFracture(Parameters::AllParameters &prms)
    : AbstractMultiphysics<dim>(prms),
      elasticity(dim, (this->ctl).params.boundary_from, "newton", this->ctl),
      phasefield(prms.phase_field_scheme, this->ctl) {}

template <int dim> void PhaseFieldFracture<dim>::setup_system() {
  this->ctl.debug_dcout << "Initialize system - elasticity" << std::endl;
  elasticity.setup_system(this->ctl);
  if ((this->ctl).params.enable_phase_field) {
    this->ctl.debug_dcout << "Initialize system - phase field" << std::endl;
    phasefield.setup_system(this->ctl);
  }
}

template <int dim> void PhaseFieldFracture<dim>::record_old_solution() {
  elasticity.record_old_solution(this->ctl);
  if ((this->ctl).params.enable_phase_field) {
    phasefield.record_old_solution(this->ctl);
  }
}

template <int dim> void PhaseFieldFracture<dim>::return_old_solution() {
  elasticity.return_old_solution(this->ctl);
  if ((this->ctl).params.enable_phase_field) {
    phasefield.return_old_solution(this->ctl);
  }
}

template <int dim> double PhaseFieldFracture<dim>::staggered_scheme() {
  if ((this->ctl).params.enable_phase_field) {
    (this->ctl).dcout << "Staggered scheme - Solving phase field" << std::endl;
    (this->ctl).computing_timer.enter_subsection("Solve phase field");
    double newton_reduction_phasefield = phasefield.update(this->ctl);
    (this->ctl).debug_dcout << "Staggered scheme - Solving phase field - point_history" << std::endl;
    (this->ctl).finalize_point_history();
    (this->ctl).debug_dcout << "Staggered scheme - Solving phase field - phase field limitation" << std::endl;
    phasefield.enforce_phase_field_limitation(this->ctl);
    (this->ctl).computing_timer.leave_subsection("Solve phase field");

    (this->ctl).dcout << "Staggered scheme - Solving elasticity" << std::endl;
    (this->ctl).computing_timer.enter_subsection("Solve elasticity");
    double newton_reduction_elasticity = elasticity.update(this->ctl);
    (this->ctl).debug_dcout << "Staggered scheme - Solving elasticity - point_history" << std::endl;
    (this->ctl).finalize_point_history();
    (this->ctl).computing_timer.leave_subsection("Solve elasticity");

    return std::max(newton_reduction_elasticity, newton_reduction_phasefield);
  } else {
    (this->ctl).dcout
        << "Solve Newton system - staggered scheme - Solving elasticity"
        << std::endl;
    (this->ctl).computing_timer.enter_subsection("Solve elasticity");
    double newton_reduction_elasticity = elasticity.update(this->ctl);
    (this->ctl).finalize_point_history();
    (this->ctl).computing_timer.leave_subsection("Solve elasticity");
    return newton_reduction_elasticity;
  }
}

template <int dim>
void PhaseFieldFracture<dim>::respective_output_results(
    DataOut<dim> &data_out) {
  (this->ctl).dcout << "Computing output - elasticity" << std::endl;
  elasticity.output_results(data_out, this->ctl);
  if ((this->ctl).params.enable_phase_field) {
    (this->ctl).dcout << "Computing output - phase field" << std::endl;
    phasefield.output_results(data_out, this->ctl);
  }
}

template <int dim> void PhaseFieldFracture<dim>::refine_grid() {
  typename DoFHandler<dim>::active_cell_iterator
      cell = phasefield.dof_handler.begin_active(),
      endc = phasefield.dof_handler.end();

  FEValues<dim> fe_values(phasefield.fe, (this->ctl).quadrature_formula,
                          update_gradients);

  unsigned int n_q_points = (this->ctl).quadrature_formula.size();
  std::vector<Tensor<1, dim>> phasefield_grads(n_q_points);

  // Define refinement criterion and mark cells to refine
  unsigned int will_refine = 0;
  double a1 = (this->ctl).params.refine_influence_initial;
  double a2 = (this->ctl).params.refine_influence_final;
  double phi_ref = std::exp(-a2) / std::exp(-a1);
  for (; cell != endc; ++cell) {
    if (cell->is_locally_owned()) {
      if (cell->diameter() < (this->ctl).params.l_phi *
                                 (this->ctl).params.refine_minimum_size_ratio) {
        cell->clear_refine_flag();
        continue;
      }
      fe_values.reinit(cell);
      fe_values.get_function_gradients((phasefield.solution), phasefield_grads);
      double max_grad = 0;
      for (unsigned int q = 0; q < n_q_points; ++q) {
        double prod = std::sqrt(phasefield_grads[q]*phasefield_grads[q]);
        max_grad = std::max(max_grad, prod);
      }
      if (max_grad > 1 / (this->ctl).params.l_phi * phi_ref * exp(-a1)) {
        cell->set_refine_flag();
        will_refine = 1;
      }
    }
  }
  (this->ctl).debug_dcout << "Refine - finish marking" << std::endl;
  double will_refine_global =
      Utilities::MPI::sum(will_refine, (this->ctl).mpi_com);
  if (!static_cast<bool>(will_refine_global)) {
    (this->ctl).dcout << "No cell to refine" << std::endl;
  } else {
    (this->ctl).debug_dcout << "Refine - prepare" << std::endl;
    // Prepare transferring of point history
    parallel::distributed::ContinuousQuadratureDataTransfer<dim, PointHistory>
        point_history_transfer(FE_Q<dim>((this->ctl).params.poly_degree),
                               QGauss<dim>((this->ctl).params.poly_degree + 1),
                               QGauss<dim>((this->ctl).params.poly_degree + 1));
    point_history_transfer.prepare_for_coarsening_and_refinement(
        (this->ctl).triangulation, (this->ctl).quadrature_point_history);

    // Prepare transferring of fields
    parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector>
        soltrans_elasticity = elasticity.prepare_refine();
    parallel::distributed::SolutionTransfer<dim, LA::MPI::Vector>
        soltrans_phasefield = phasefield.prepare_refine();

    (this->ctl).debug_dcout << "Refine - start refinement" << std::endl;
    // Execute refinement
    (this->ctl).triangulation.execute_coarsening_and_refinement();
    setup_system();

    (this->ctl).debug_dcout << "Refine - after refinement - point history"
                            << std::endl;
    // Finalize transferring of point history
    (this->ctl).initialize_point_history();
    point_history_transfer.interpolate();
    (this->ctl).debug_dcout << "Refine - after refinement - transfer fields"
                            << std::endl;
    // Finalize transferring of fields
    elasticity.post_refine(soltrans_elasticity, this->ctl);
    phasefield.post_refine(soltrans_phasefield, this->ctl);
    (this->ctl).debug_dcout << "Refine - done" << std::endl;
  }
}

#endif