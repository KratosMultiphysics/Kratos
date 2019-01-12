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
//#include "linear_system/linear_solvers/superlu_direct_solver.hpp"
#include "linear_system/linear_solvers/superlu_mt_direct_solver.hpp"

namespace Kratos
{

namespace Python
{
namespace py = pybind11;

void AddCustomLinearSolversToPython(pybind11::module& m)
{
  typedef UblasSpace<double, CompressedMatrix, Vector>                       SparseSpaceType;
  typedef UblasSpace<double, Matrix, Vector>                                  LocalSpaceType;
  typedef DirectSolver<SparseSpaceType, LocalSpaceType>                     DirectSolverType;
  //typedef SuperLUDirectSolver<SparseSpaceType, LocalSpaceType>       SuperLUDirectSolverType;
  typedef SuperLUmtDirectSolver<SparseSpaceType, LocalSpaceType>   SuperLUmtDirectSolverType;
  //typedef SuperLUIterativeSolver<SparseSpaceType, LocalSpaceType> SuperLUIterativeSolverType;

  // py::class_<SuperLUDirectSolverType, typename SuperLUDirectSolverType::Pointer, DirectSolverType>
  //     (m, "SuperLU_DirectSolver")
  //     .def(py::init<>() )
  //     .def(py::init<Parameters>());

  py::class_<SuperLUmtDirectSolverType, typename SuperLUmtDirectSolverType::Pointer, DirectSolverType>
      (m, "SuperLU_MT_DirectSolver")
      .def(py::init<>() )
      .def(py::init<Parameters>());

  // py::class_<SuperLUIterativeSolverType, typename SuperLUIterativeSolverType::Pointer, SuperLUDirectSolverType>
  //     (m, "SuperLU_IterativeSolver")
  //     .def(py::init<>() )
  //     .def(py::init<Parameters>());

}

}  // namespace Python.

} // Namespace Kratos
