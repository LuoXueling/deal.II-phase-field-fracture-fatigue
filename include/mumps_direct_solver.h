//
// Faithful drop-in for TrilinosWrappers::SolverDirect restricted to
// Amesos_Mumps, with the single addition of an Amesos SetParameters() call that
// raises MUMPS ICNTL(14) (percentage increase in the estimated working space).
//
// The stock TrilinosWrappers::SolverDirect (deal.II 9.4-9.6) never calls
// SetParameters and exposes no way to set ICNTL, so it occasionally fails the
// numeric factorization with MUMPS INFOG(1)=-9 ("main internal real workarray S
// too small") on ill-balanced 3D problems. The initialize()/solve() bodies are
// copied from deal.II 9.5.2 (source/lac/trilinos_solver.cc) with two deviations,
// both marked below:
//   1. the Amesos_Mumps solver is constructed directly instead of through the
//      Amesos factory. Amesos::Create("Amesos_Mumps", problem) internally does
//      exactly `new Amesos_Mumps(problem)` (Trilinos packages/amesos/src/
//      Amesos.cpp), so this is functionally identical for a MUMPS build, but it
//      avoids including the heavy <Amesos.h> factory header (which pulls in every
//      enabled Amesos backend) into every translation unit that includes this.
//   2. the SetParameters() block that raises ICNTL(14). The "mumps" sublist with
//      integer "ICNTL(i)" entries is the officially supported path (Trilinos
//      packages/amesos/src/Amesos_Mumps.cpp).
//
// Guarding mirrors the Amesos factory exactly. <Amesos_Mumps.h> #includes
// dmumps_c.h unconditionally, so it (and the real solver) cannot be compiled
// unless Amesos was built with MUMPS+MPI. We therefore compile the real class
// only under #if defined(HAVE_AMESOS_MUMPS) && defined(HAVE_MPI) -- the same
// condition Amesos::Create/Query use -- and otherwise compile a stub that, like
// the factory's else-branch, raises an error at runtime when MUMPS is selected
// but not configured.
//

#ifndef CRACKS_MUMPS_DIRECT_SOLVER_H
#define CRACKS_MUMPS_DIRECT_SOLVER_H

#include "dealii_includes.h"

#ifdef DEAL_II_WITH_TRILINOS

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/trilinos_solver.h> // TrilinosWrappers::SolverDirect::ExcTrilinosError
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <Amesos_ConfigDefs.h> // HAVE_AMESOS_MUMPS, HAVE_MPI

using namespace dealii;

#if defined(HAVE_AMESOS_MUMPS) && defined(HAVE_MPI)

// Real implementation -- compiled only when Amesos provides MUMPS. These headers
// (notably <Amesos_Mumps.h>, which pulls in dmumps_c.h) only exist/compile on a
// MUMPS-enabled build, which is exactly why they live behind this guard.
#include <Amesos_BaseSolver.h>
#include <Amesos_Mumps.h>
#include <Epetra_LinearProblem.h>
#include <Teuchos_ParameterList.hpp>

#include <iostream>
#include <memory>

class MumpsDirectSolver {
public:
  MumpsDirectSolver(SolverControl &cn, int icntl14 = 20,
                    bool output_solver_details = false)
    : solver_control(cn), icntl14(icntl14),
      output_solver_details(output_solver_details) {}

  // TrilinosWrappers::SolverDirect::initialize(const SparseMatrix &A), with the
  // two marked deviations (direct construction, SetParameters).
  void initialize(const TrilinosWrappers::SparseMatrix &A) {
    // We need an Epetra_LinearProblem object to let the Amesos solver know
    // about the matrix and vectors.
    linear_problem = std::make_unique<Epetra_LinearProblem>();

    // Assign the matrix operator to the Epetra_LinearProblem object
    linear_problem->SetOperator(
      const_cast<Epetra_CrsMatrix *>(&A.trilinos_matrix()));

    // Fetch return value of Amesos Solver functions
    int ierr;

    // First set whether we want to print the solver information to screen or
    // not.
    ConditionalOStream verbose_cout(std::cout, output_solver_details);

    // ---- deviation 1: construct Amesos_Mumps directly (== what the Amesos
    // factory does for "Amesos_Mumps"), avoiding the heavy <Amesos.h>. ----
    solver = std::make_unique<Amesos_Mumps>(*linear_problem);

    // ---- deviation 2: raise the MUMPS working-space margin via the documented
    // "mumps" sublist to avoid intermittent INFOG(1)=-9. ----
    Teuchos::ParameterList params;
    params.sublist("mumps").set("ICNTL(14)", icntl14);
    ierr = solver->SetParameters(params);
    AssertThrow(ierr == 0,
                TrilinosWrappers::SolverDirect::ExcTrilinosError(ierr));
    // ----------------------------------------------------------------------

    verbose_cout << "Starting symbolic factorization" << std::endl;
    ierr = solver->SymbolicFactorization();
    AssertThrow(ierr == 0,
                TrilinosWrappers::SolverDirect::ExcTrilinosError(ierr));

    verbose_cout << "Starting numeric factorization" << std::endl;
    ierr = solver->NumericFactorization();
    AssertThrow(ierr == 0,
                TrilinosWrappers::SolverDirect::ExcTrilinosError(ierr));
  }

  // Verbatim TrilinosWrappers::SolverDirect::solve(MPI::Vector &x,
  //                                                const MPI::Vector &b).
  void solve(TrilinosWrappers::MPI::Vector &x,
             const TrilinosWrappers::MPI::Vector &b) {
    // Assign the empty LHS vector to the Epetra_LinearProblem object
    linear_problem->SetLHS(&x.trilinos_vector());

    // Assign the RHS vector to the Epetra_LinearProblem object
    linear_problem->SetRHS(
      const_cast<Epetra_MultiVector *>(&b.trilinos_vector()));

    // First set whether we want to print the solver information to screen or
    // not.
    ConditionalOStream verbose_cout(std::cout, output_solver_details);

    verbose_cout << "Starting solve" << std::endl;

    // Fetch return value of Amesos Solver functions
    int ierr = solver->Solve();
    AssertThrow(ierr == 0,
                TrilinosWrappers::SolverDirect::ExcTrilinosError(ierr));

    // Finally, force the SolverControl object to report convergence
    solver_control.check(0, 0);
  }

  SolverControl &control() const { return solver_control; }

private:
  SolverControl &solver_control;
  int icntl14;
  bool output_solver_details;
  std::unique_ptr<Epetra_LinearProblem> linear_problem;
  std::unique_ptr<Amesos_BaseSolver> solver;
};

#else // !(HAVE_AMESOS_MUMPS && HAVE_MPI)

// Stub -- no Amesos_Mumps dependency, so it compiles on a build whose Amesos
// lacks MUMPS. Mirrors the Amesos factory's else-branch (which yields no usable
// solver for "Amesos_Mumps"): selecting this solver raises the same error
// deal.II's SolverDirect raises (via Amesos::Query) for an unconfigured backend.
// The throw is in the constructor so the failure surfaces immediately at the
// point of selection, exactly like the factory's Create.
class MumpsDirectSolver {
public:
  MumpsDirectSolver(SolverControl &cn, int = 20, bool = false)
    : solver_control(cn) {
    AssertThrow(false,
                ExcMessage(
                  "You tried to select the solver type <Amesos_Mumps> but this "
                  "solver is not supported by Trilinos either because it does "
                  "not exist, or because Trilinos was not configured for its "
                  "use."));
  }

  void initialize(const TrilinosWrappers::SparseMatrix &) {}

  void solve(TrilinosWrappers::MPI::Vector &,
             const TrilinosWrappers::MPI::Vector &) {}

  SolverControl &control() const { return solver_control; }

private:
  SolverControl &solver_control;
};

#endif // HAVE_AMESOS_MUMPS && HAVE_MPI

#endif // DEAL_II_WITH_TRILINOS

#endif // CRACKS_MUMPS_DIRECT_SOLVER_H
