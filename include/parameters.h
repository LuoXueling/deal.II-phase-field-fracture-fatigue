/**
 * Xueling Luo @ Shanghai Jiao Tong University, 2022
 * This code is for multiscale phase field fracture.
 **/

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "dealii_includes.h"
#include <ctime>

namespace Parameters {
    struct Project {
        std::string mesh_from;
        std::string boundary_from;
        std::string project_name;
        std::string output_dir_top;
        std::string output_dir_sub;
        std::string load_sequence_from;
        bool enable_phase_field;
        bool enable_fatigue;
        bool debug_output;

        static void subsection_declare_parameters(ParameterHandler &prm);

        void subsection_parse_parameters(ParameterHandler &prm);
    };

    void Project::subsection_declare_parameters(ParameterHandler &prm) {
        prm.enter_subsection("Project"); {
            prm.declare_entry("Mesh from", "script",
                              Patterns::FileName(Patterns::FileName::FileType::input));
            prm.declare_entry("Boundary from", "none",
                              Patterns::FileName(Patterns::FileName::FileType::input));
            prm.declare_entry("Project name", "Default project",
                              Patterns::FileName(Patterns::FileName::FileType::output));
            prm.declare_entry("Output directory", "../output/",
                              Patterns::FileName(Patterns::FileName::FileType::input));
            prm.declare_entry("Output sub-directory", "",
                              Patterns::FileName(Patterns::FileName::FileType::input));
            prm.declare_entry("Load sequence from", "script",
                              Patterns::FileName(Patterns::FileName::FileType::input));
            prm.declare_entry("Enable phase field", "true", Patterns::Bool());
            prm.declare_entry("Enable fatigue", "false", Patterns::Bool());

            prm.declare_entry("Debug output", "false", Patterns::Bool());
        }
        prm.leave_subsection();
    }

    void Project::subsection_parse_parameters(ParameterHandler &prm) {
        prm.enter_subsection("Project"); {
            mesh_from = prm.get("Mesh from");
            boundary_from = prm.get("Boundary from");
            project_name = prm.get("Project name");
            output_dir_top = prm.get("Output directory");
            output_dir_sub = prm.get("Output sub-directory");
            load_sequence_from = prm.get("Load sequence from");
            enable_phase_field = prm.get_bool("Enable phase field");
            enable_fatigue = prm.get_bool("Enable fatigue");
            debug_output = prm.get_bool("Debug output");
        }
        prm.leave_subsection();
    }

    struct Runtime {
        unsigned int max_no_timesteps;
        std::string adaptive_timestep;
        std::string adaptive_timestep_parameters;
        double timestep;
        double timestep_size_2;
        unsigned int switch_timestep;
        std::string norm_type;
        bool direct_solver;
        std::string direct_solver_type;
        double lower_bound_newton_residual;
        unsigned int max_no_newton_steps;
        bool skip_first_iter;
        double upper_newton_rho;
        std::string adjustment_method;
        std::string adjustment_method_elasticity;
        std::string linesearch_parameters;
        std::string modified_newton_parameters;
        unsigned int max_adjustment_steps;
        bool is_monolithic;
        bool multipass_staggered;
        unsigned int max_multipass;
        double multipass_residual_tol;
        bool quit_multipass_if_increase;
        bool throw_if_multipass_increase;
        std::string phase_field_scheme;
        std::string decomposition;
        double constant_k;
        unsigned int save_vtk_per_step;
        double max_crack_length;
        double fix_phasefield_near_boundary_distance;

        static void subsection_declare_parameters(ParameterHandler &prm);

        void subsection_parse_parameters(ParameterHandler &prm);
    };

    void Runtime::subsection_declare_parameters(ParameterHandler &prm) {
        prm.enter_subsection("Runtime"); {
            prm.declare_entry("Max No of timesteps", "1", Patterns::Integer(0));

            prm.declare_entry(
                "Adaptive timestep", "exponential",
                Patterns::Selection(
                    "exponential|constant|KristensenCLA|"
                    "CojocaruCycleJump|LiCycleJump|JonasCycleJump|YangCycleJump|JacconCycleJump"));
            prm.declare_entry("Adaptive timestep parameters", "", Patterns::Anything());
            prm.declare_entry("Timestep size", "1.0", Patterns::Double(0));

            prm.declare_entry("Timestep size to switch to", "1.0", Patterns::Double(0));

            prm.declare_entry("Switch timestep after steps", "0", Patterns::Integer(0));
            prm.declare_entry("Use Direct Inner Solver", "false", Patterns::Bool());
            prm.declare_entry("Use Direct Inner Solver type", "Amesos_Klu", Patterns::Anything());

            prm.declare_entry("Norm type", "linfty",
                              Patterns::Selection("linfty|l2|l1"));

            prm.declare_entry("Newton lower bound", "1.0e-10", Patterns::Double(0));

            prm.declare_entry("Newton maximum steps", "10", Patterns::Integer(0));

            prm.declare_entry("Upper Newton rho", "0.999", Patterns::Double(0));
            prm.declare_entry("Allow skip first Newton iteration", "true",
                              Patterns::Bool());

            prm.declare_entry(
                "Adjustment method", "linesearch",
                Patterns::Selection("none|linesearch|KristensenModifiedNewton"));

            prm.declare_entry(
                "Adjustment method for elasticity", "linesearch",
                Patterns::Selection(
                    "none|linesearch|AndersonNewton|KristensenModifiedNewton"));

            prm.declare_entry("Parameters of line search", "0.1", Patterns::Anything());
            prm.declare_entry("Parameters of modified newton", "",
                              Patterns::Anything());

            prm.declare_entry("Maximum number of adjustment steps of Newton solution",
                              "5", Patterns::Integer(0));

            prm.declare_entry("Line search damping", "0.5", Patterns::Double(0));

            prm.declare_entry("Use monolithic", "false", Patterns::Bool());

            prm.declare_entry("Use multipass staggered", "false", Patterns::Bool());
            prm.declare_entry("Maximum number of multipass steps", "5",
                              Patterns::Integer(0));
            prm.declare_entry("Residual tolerance of multipass", "1e-8",
                              Patterns::Double(0));
            prm.declare_entry("Quit multipass if residual increasing", "true",
                              Patterns::Bool());
            prm.declare_entry("Throw if multipass residual increasing", "false",
                              Patterns::Bool());

            prm.declare_entry("Phase field update", "newton",
                              Patterns::Selection("newton|linear"));

            prm.declare_entry(
                "Decomposition", "hybrid",
                Patterns::Selection("none|hybrid|sphere|eigen|hybridnotension"));

            prm.declare_entry("Constant small quantity k", "1.0e-6",
                              Patterns::Double(0));
            prm.declare_entry("Save vtk per step", "1", Patterns::Integer(0));

            prm.declare_entry("Maximum crack length", "1e8", Patterns::Double(0));

            // Fix the phase field to 0 (intact) at nodes within this distance of
            // any boundary face carrying a Dirichlet/Neumann BC. Negative = disabled.
            prm.declare_entry("Fix phase field near boundary", "-1",
                              Patterns::Double());
        }
        prm.leave_subsection();
    }

    void Runtime::subsection_parse_parameters(ParameterHandler &prm) {
        prm.enter_subsection("Runtime"); {
            max_no_timesteps = prm.get_integer("Max No of timesteps");
            adaptive_timestep = prm.get("Adaptive timestep");
            adaptive_timestep_parameters = prm.get("Adaptive timestep parameters");
            timestep = prm.get_double("Timestep size");
            timestep_size_2 = prm.get_double("Timestep size to switch to");
            switch_timestep = prm.get_integer("Switch timestep after steps");

            direct_solver = prm.get_bool("Use Direct Inner Solver");
            direct_solver_type = prm.get("Use Direct Inner Solver type");

            norm_type = prm.get("Norm type");

            // Newton tolerances and maximum steps
            lower_bound_newton_residual = prm.get_double("Newton lower bound");
            max_no_newton_steps = prm.get_integer("Newton maximum steps");

            skip_first_iter = prm.get_bool("Allow skip first Newton iteration");

            // Criterion when time step should be cut
            // Higher number means: almost never
            // only used for simple penalization
            upper_newton_rho = prm.get_double("Upper Newton rho");

            adjustment_method = prm.get("Adjustment method");
            adjustment_method_elasticity = prm.get("Adjustment method for elasticity");
            max_adjustment_steps = prm.get_integer(
                "Maximum number of adjustment steps of Newton solution");
            linesearch_parameters = prm.get("Parameters of line search");
            modified_newton_parameters = prm.get("Parameters of modified newton");

            is_monolithic = prm.get_bool("Use monolithic");

            multipass_staggered = prm.get_bool("Use multipass staggered");
            max_multipass = prm.get_integer("Maximum number of multipass steps");
            multipass_residual_tol = prm.get_double("Residual tolerance of multipass");
            quit_multipass_if_increase =
                    prm.get_bool("Quit multipass if residual increasing");
            throw_if_multipass_increase =
                    prm.get_bool("Throw if multipass residual increasing");

            phase_field_scheme = prm.get("Phase field update");
            decomposition = prm.get("Decomposition");

            constant_k = prm.get_double("Constant small quantity k");

            save_vtk_per_step = prm.get_integer("Save vtk per step");
            max_crack_length = prm.get_double("Maximum crack length");
            fix_phasefield_near_boundary_distance =
                    prm.get_double("Fix phase field near boundary");
        }
        prm.leave_subsection();
    }

    struct Material {
        double E;
        double v;
        double Gc;
        double l_phi;
        double lambda;
        double mu;
        double lame_coefficient_mu;
        double lame_coefficient_lambda;
        double density;
        std::string plane_state;
        std::string phasefield_model;
        std::string degradation;
        std::string fatigue_degradation;
        std::string fatigue_degradation_parameters;
        std::string fatigue_accumulation;
        std::string fatigue_accumulation_parameters;
        std::string fatigue_increment;
        std::string fatigue_increment_parameters;
        std::string fatigue_alpha_t;

        static void subsection_declare_parameters(ParameterHandler &prm);

        void subsection_parse_parameters(ParameterHandler &prm);
    };

    void Material::subsection_declare_parameters(ParameterHandler &prm) {
        prm.enter_subsection("Material"); {
            prm.declare_entry("Young's modulus", "1000", Patterns::Double(0));
            prm.declare_entry("Poisson's ratio", "0.3", Patterns::Double(0, 0.5));
            prm.declare_entry("Critical energy release rate", "1", Patterns::Double(0));
            prm.declare_entry("Phase field length scale", "0.01", Patterns::Double(0));
            prm.declare_entry("Density (for self weight)", "0.0", Patterns::Double(0));
            prm.declare_entry("Plane state", "stress",
                              Patterns::Selection("stress|strain"));
            prm.declare_entry("Phase field model", "AT2",
                              Patterns::Selection("AT2|AT1"));
            prm.declare_entry("Degradation", "quadratic",
                              Patterns::Selection("quadratic|cubic"));
            prm.declare_entry(
                "Fatigue degradation", "CarraraAsymptotic",
                Patterns::Selection("CarraraAsymptotic|CarraraLogarithmic|"
                    "KristensenAsymptotic|CojocaruAsymptotic"));
            prm.declare_entry("Fatigue degradation parameters", "",
                              Patterns::Anything());
            prm.declare_entry(
                "Fatigue accumulation", "CarraraNoMeanEffect",
                Patterns::Selection(
                    "CarraraNoMeanEffect|CarraraMeanEffect|Kristensen|"
                    "KristensenCLA|Cojocaru|CojocaruCLA|Li|LiCLA|Jonas|JonasCLA|JonasNodegrade|"
                    "Yang|Jaccon|JacconNodegrade"));
            prm.declare_entry("Fatigue accumulation parameters", "",
                              Patterns::Anything());
            // Which per-cycle increment law the resolved-cycle branch of an
            // acceleration algorithm uses. "Auto" keeps each algorithm's
            // built-in choice, so existing parameter files are unaffected.
            // CarraraNoMeanEffect/Kristensen/CarraraMeanEffect accumulations
            // are fixed to their own law and reject anything else.
            prm.declare_entry(
                "Fatigue increment", "Auto",
                Patterns::Selection("Auto|CarraraNoMeanEffect|Kristensen|"
                    "CarraraMeanEffect"));
            // Parameters for the increment law selected above, when it needs
            // any. Only CarraraMeanEffect reads it, taking the leading entry as
            // alpha_n. Empty (the default) falls back to "Fatigue alpha_t" if
            // set, else the shared Gc/(12*l_phi) default -- so existing
            // parameter files are unaffected. This exists because the host
            // acceleration algorithm already owns "Fatigue accumulation
            // parameters", leaving the increment law no slot of its own.
            prm.declare_entry("Fatigue increment parameters", "",
                              Patterns::Anything());
            // Global fatigue threshold. Empty (the default) means every
            // consumer keeps its own hardcoded formulation, so existing
            // parameter files are unaffected. When set, it overrides that
            // default everywhere alpha_t/alpha_n would otherwise be derived
            // from Gc/l_phi: the Carrara/Kristensen asymptotic degradations,
            // the CarraraMeanEffect accumulation's alpha_n, and JonasCycleJump.
            // A scheme's own parameter string still wins over this if present.
            prm.declare_entry("Fatigue alpha_t", "", Patterns::Anything());
        }
        prm.leave_subsection();
    }

    void Material::subsection_parse_parameters(ParameterHandler &prm) {
        prm.enter_subsection("Material"); {
            E = prm.get_double("Young's modulus");
            v = prm.get_double("Poisson's ratio");
            Gc = prm.get_double("Critical energy release rate");
            l_phi = prm.get_double("Phase field length scale");
            density = prm.get_double("Density (for self weight)");
            plane_state = prm.get("Plane state");
            phasefield_model = prm.get("Phase field model");
            degradation = prm.get("Degradation");
            fatigue_degradation = prm.get("Fatigue degradation");
            fatigue_degradation_parameters = prm.get("Fatigue degradation parameters");
            fatigue_accumulation = prm.get("Fatigue accumulation");
            fatigue_accumulation_parameters =
                    prm.get("Fatigue accumulation parameters");
            fatigue_increment = prm.get("Fatigue increment");
            fatigue_increment_parameters =
                    prm.get("Fatigue increment parameters");
            fatigue_alpha_t = prm.get("Fatigue alpha_t");
        }
        prm.leave_subsection();
        lame_coefficient_mu = E / (2.0 * (1 + v));
        lame_coefficient_lambda = (2 * v * lame_coefficient_mu) / (1.0 - 2 * v);
    }

    struct FESystemInfo {
        unsigned int dim;
        unsigned int poly_degree;
        unsigned int quad_order;
        // "hex" (quadrilateral in 2D, hexahedron in 3D) or "tet" (tetrahedron).
        // Selects the reference cell, and with it the finite element,
        // quadrature, mapping, triangulation type and mesh reader.
        std::string element_type;
        bool refine;
        double refine_influence_initial;
        double refine_influence_final;
        double refine_minimum_size_ratio;

        static void subsection_declare_parameters(ParameterHandler &prm);

        void subsection_parse_parameters(ParameterHandler &prm);
    };

    void FESystemInfo::subsection_declare_parameters(ParameterHandler &prm) {
        prm.enter_subsection("Finite element system"); {
            prm.declare_entry("Physical dimension", "2", Patterns::Integer(0),
                              "Physical dimension");
            prm.declare_entry("Polynomial degree", "2", Patterns::Integer(0),
                              "Displacement system polynomial order");

            prm.declare_entry("Quadrature order", "3", Patterns::Integer(0),
                              "Gauss quadrature order");
            prm.declare_entry("Element type", "hex",
                              Patterns::Selection("hex|tet"),
                              "Reference cell of the mesh. 'hex' (the default, "
                              "quadrilateral in 2D / hexahedron in 3D) is read "
                              "from an Abaqus .inp file and supports adaptive "
                              "refinement. 'tet' (3D linear tetrahedra) is read "
                              "from a gmsh .msh file (convert an Abaqus C3D4 "
                              ".inp with meshes/inp2msh.py) and does not "
                              "support adaptive refinement.");
            prm.declare_entry("Refine", "false", Patterns::Bool());

            prm.declare_entry("Phase field initial influential ratio (for refinement)",
                              "2", Patterns::Double(0));
            prm.declare_entry("Phase field final influential ratio (for refinement)",
                              "3", Patterns::Double(0));
            prm.declare_entry("Minimum relative size of refined cells w.r.t l_phi",
                              "0.2", Patterns::Double(0));
        }
        prm.leave_subsection();
    }

    void FESystemInfo::subsection_parse_parameters(ParameterHandler &prm) {
        prm.enter_subsection("Finite element system"); {
            dim = prm.get_integer("Physical dimension");
            poly_degree = prm.get_integer("Polynomial degree");
            quad_order = prm.get_integer("Quadrature order");
            element_type = prm.get("Element type");
            refine = prm.get_bool("Refine");
            refine_influence_final =
                    prm.get_double("Phase field final influential ratio (for refinement)");
            refine_influence_initial = prm.get_double(
                "Phase field initial influential ratio (for refinement)");
            refine_minimum_size_ratio =
                    prm.get_double("Minimum relative size of refined cells w.r.t l_phi");

            // Fail loudly on a combination we cannot honour, rather than
            // silently ignoring the user's intent later on.
            if (element_type == "tet") {
                AssertThrow(dim == 3,
                            ExcMessage("'Element type = tet' is only supported "
                                "for 'Physical dimension = 3'."));
                AssertThrow(poly_degree >= 1 && poly_degree <= 2,
                            ExcMessage("'Element type = tet' requires "
                                "'Polynomial degree' of 1 or 2; FE_SimplexP is "
                                "not implemented for higher degrees."));
                AssertThrow(!refine,
                            ExcMessage("'Refine = true' is not supported for "
                                "'Element type = tet': deal.II cannot adaptively "
                                "refine simplex meshes. Set 'Refine = false'."));
            }
        }
        prm.leave_subsection();
    }

    struct AllParameters : public FESystemInfo,
                           public Project,
                           public Runtime,
                           public Material {
        AllParameters() = default;

        AllParameters(const std::string &input_file);

        static void declare_parameters(ParameterHandler &prm);

        void set_parameters(const std::string &input_file);

        void parse_parameters(ParameterHandler &prm);

        std::string output_dir;
        std::string param_dir;
    };

    AllParameters::AllParameters(const std::string &input_file) {
        set_parameters(input_file);
    }

    void AllParameters::set_parameters(const std::string &input_file) {
        ParameterHandler prm;
        declare_parameters(prm);
        prm.parse_input(input_file);
        parse_parameters(prm);

        // set output directory
        std::string stime;
        if (output_dir_sub.length() == 0) {
            std::time_t currenttime = std::time(0);
            char tAll[255];
            std::strftime(tAll, sizeof(tAll), "%Y-%m-%d-%H-%M-%S",
                          std::localtime(&currenttime));
            std::stringstream strtime;
            strtime << tAll;
            stime = strtime.str();
            output_dir = output_dir_top + this->project_name + "-" + stime + "/";
        } else {
            output_dir = output_dir_top + output_dir_sub + "/";
        }

        param_dir = input_file;
    }

    void AllParameters::declare_parameters(ParameterHandler &prm) {
        FESystemInfo::subsection_declare_parameters(prm);
        Project::subsection_declare_parameters(prm);
        Runtime::subsection_declare_parameters(prm);
        Material::subsection_declare_parameters(prm);
    }

    void AllParameters::parse_parameters(ParameterHandler &prm) {
        FESystemInfo::subsection_parse_parameters(prm);
        Project::subsection_parse_parameters(prm);
        Runtime::subsection_parse_parameters(prm);
        Material::subsection_parse_parameters(prm);
    }
} // namespace Parameters

#endif
