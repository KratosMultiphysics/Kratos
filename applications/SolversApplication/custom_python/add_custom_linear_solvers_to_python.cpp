//
//   Project Name:        KratosSolversApplication $
//   Created by:          $Author:     JMCarbonell $
//   Last modified by:    $Co-Author:              $
//   Date:                $Date:      January 2019 $
//   Revision:            $Revision:           0.0 $
//
//

// System includes

// External includes

// Project includes
#include "spaces/ublas_space.h"
#include "custom_python/add_custom_linear_solvers_to_python.h"

//linear solvers
#ifdef INCLUDE_SUPERLU_MT
  #include "linear_system/linear_solvers/superlu_mt_direct_solver.hpp"
#else
  #include "linear_system/linear_solvers/superlu_direct_solver.hpp"
#endif

#ifdef INCLUDE_FEAST
  #include "linear_system/linear_solvers/feast_solver.hpp"
#endif

// #include "solvers_application.h"
// #include "includes/standard_linear_solver_factory.h"

namespace Kratos
{

namespace Python
{
namespace py = pybind11;

void AddCustomLinearSolversToPython(pybind11::module& m)
{
  typedef UblasSpace<double, CompressedMatrix, Vector>                       SparseSpaceType;
  typedef UblasSpace<double, Matrix, Vector>                                  LocalSpaceType;

#ifdef INCLUDE_SUPERLU_MT
  typedef SuperLUmtDirectSolver<SparseSpaceType, LocalSpaceType>   SuperLUmtDirectSolverType;

  py::class_<SuperLUmtDirectSolverType, typename SuperLUmtDirectSolverType::Pointer, DirectSolverType>
      (m, "SuperLU_MT_DirectSolver")
      .def(py::init<>() )
      .def(py::init<Parameters>());
#else
  typedef DirectSolver<SparseSpaceType, LocalSpaceType>                     DirectSolverType;
  typedef SuperLUDirectSolver<SparseSpaceType, LocalSpaceType>       SuperLUDirectSolverType;
  //typedef SuperLUIterativeSolver<SparseSpaceType, LocalSpaceType> SuperLUIterativeSolverType;

  py::class_<SuperLUDirectSolverType, typename SuperLUDirectSolverType::Pointer, DirectSolverType>
      (m, "SuperLU_DirectSolver")
      .def(py::init<>() )
      .def(py::init<Parameters>());

  // py::class_<SuperLUIterativeSolverType, typename SuperLUIterativeSolverType::Pointer, SuperLUDirectSolverType>
  //     (m, "SuperLU_IterativeSolver")
  //     .def(py::init<>() )
  //     .def(py::init<Parameters>());
#endif

#ifdef INCLUDE_FEAST
  typedef FEASTSolver<SparseSpaceType, LocalSpaceType>                          FEASTEigenValueSolverType;
  typedef LinearSolver<SparseSpaceType, LocalSpaceType>                                  LinearSolverType;

  typedef UblasSpace<std::complex<double>, ComplexCompressedMatrix, ComplexVector> ComplexSparseSpaceType;
  typedef UblasSpace<std::complex<double>, ComplexMatrix, ComplexVector>            ComplexLocalSpaceType;
  typedef LinearSolver<ComplexSparseSpaceType, ComplexLocalSpaceType>             ComplexLinearSolverType;

  py::class_<FEASTEigenValueSolverType, FEASTEigenValueSolverType::Pointer, LinearSolverType>
      (m, "FEAST_EigenValueSolver")
      .def(py::init<Parameters>() )
      .def(py::init<Parameters, ComplexLinearSolverType::Pointer>() )
      ;
#endif


}

}  // namespace Python.

// Register located here: (to avoid multiple defined symbols when including the external C libraries)
// SolversApplicationRegisterLinearSolvers::SolversApplicationRegisterLinearSolvers()
// {
//   typedef UblasSpace<double, CompressedMatrix, Vector>                       SparseSpaceType;
//   typedef UblasSpace<double, Matrix, Vector>                                  LocalSpaceType;

// #ifdef INCLUDE_SUPERLU_MT
//   typedef SuperLUmtDirectSolver<SparseSpaceType, LocalSpaceType>   SuperLUmtDirectSolverType;

//   const StandardLinearSolverFactory<SparseSpaceType, LocalSpaceType, SuperLUmtDirectSolverType> mSuperLUmtDirectSolverFactory;
//   KRATOS_REGISTER_LINEAR_SOLVER("SuperLU_MT_DirectSolver", mSuperLUmtDirectSolverFactory);
// #else
//   typedef DirectSolver<SparseSpaceType, LocalSpaceType>                     DirectSolverType;
//   typedef SuperLUDirectSolver<SparseSpaceType, LocalSpaceType>       SuperLUDirectSolverType;
//   //typedef SuperLUIterativeSolver<SparseSpaceType, LocalSpaceType> SuperLUIterativeSolverType;

//   const StandardLinearSolverFactory<SparseSpaceType, LocalSpaceType, SuperLUDirectSolverType> mSuperLUDirectSolverFactory;
//   //const StandardLinearSolverFactory<SparseSpaceType, LocalSpaceType, SuperLUIterativeSolverType> mSuperLUIterativeSolverFactory;

//   KRATOS_REGISTER_LINEAR_SOLVER("SuperLU_DirectSolver", mSuperLUDirectSolverFactory);
//   //KRATOS_REGISTER_LINEAR_SOLVER("SuperLU_IterativeSolver", mSuperLUIterativeSolverFactory);
// #endif

// #ifdef INCLUDE_FEAST
//   typedef FEASTSolver<SparseSpaceType, LocalSpaceType> FEASTEigenValueSolverType;
//   const StandardLinearSolverFactory<SparseSpaceType, LocalSpaceType, FEASTEigenValueSolverType> mFEASTEigenValueSolverFactory;
//   KRATOS_REGISTER_LINEAR_SOLVER("FEAST_EigenValueSolver", mFEASTEigenValueSolverFactory);
// #endif
// }


} // Namespace Kratos
