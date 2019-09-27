//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    YOUR_NAME_HERE
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "costum_utilities/calculate_mean_temp.h"
#include "costum_utilities/get_rom_residuals.h"

namespace Kratos
{

namespace Python
{
using namespace pybind11;

  void  AddCustomUtilitiesToPython(pybind11::module& m)
  {
    class_<CalculateMeanTemperature, typename CalculateMeanTemperature::Pointer>(m, "CalculateMeanTemperature")
    .def(init<ModelPart&>()) // The input parameters is a model part 
    .def("Execute",&CalculateMeanTemperature::Calculate) // When we call "Execute" in python, Calculate is called in C++. Notice we don't write the input parameters here. NO, IT DOESNT LIKE THE TYPE OF THE SCHEME 
    ;

    class_<GetRomResiduals, typename GetRomResiduals::Pointer>(m, "GetRomResiduals")
    .def( Parameters, typename TSchemeType::Pointer ,init<ModelPart&  >() ) //  Will it recognize the type of the scheme?
    .def("Execute",&GetRomResiduals::Calculate) // 

  }

}  // namespace Python.

} // Namespace Kratos


//   ORIGINAL, DOES NOT WORK :(
// //    |  /           |
// //    ' /   __| _` | __|  _ \   __|
// //    . \  |   (   | |   (   |\__ `
// //   _|\_\_|  \__,_|\__|\___/ ____/
// //                   Multi-Physics
// //
// //  License:		 BSD License
// //					 Kratos default license: kratos/license.txt
// //
// //  Main authors:    Author1 Fullname
// //                   Author2 Fullname
// //


// // System includes

// // External includes
// #include <pybind11/pybind11.h>


// // Project includes
// #include "includes/define.h"
// #include "custom_python/add_custom_utilities_to_python.h"

// #include "spaces/ublas_space.h"
// #include "linear_solvers/linear_solver.h"


// namespace Kratos {
// namespace Python {

// void AddCustomUtilitiesToPython(pybind11::module& m)
// {
//     namespace py = pybind11;

//     typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
//     typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
//     typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

// }

// } // namespace Python.
// } // Namespace Kratos