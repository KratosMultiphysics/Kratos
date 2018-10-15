//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// System includes

// External includes
#include "pybind11/pybind11.h"


// Project includes
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "linear_solvers/linear_solver.h"

#include "custom_utilities/turbulence_statistics_container.h"

namespace Kratos
{

namespace Python
{


  void  AddCustomUtilitiesToPython(pybind11::module& m)
  {
    using namespace pybind11;

    //typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    //typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    //typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    class_<Variable<TurbulenceStatisticsContainer::Pointer>,VariableData>(m, "TurbulenceStatisticsContainerVariable")
    .def("__str__", KRATOS_DEF_PYTHON_STR(Variable<TurbulenceStatisticsContainer::Pointer>))
    ;
  }

}  // namespace Python.

} // Namespace Kratos
