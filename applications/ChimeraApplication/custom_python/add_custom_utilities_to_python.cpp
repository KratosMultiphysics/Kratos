//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					        Navaneeth K Narayanan
//					        Rishith Ellath Meethal
//
// External includes

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/solver_settings.h"
#include "custom_utilities/fractional_step_settings_for_chimera.h"

//Processes
namespace Kratos
{

namespace Python
{

    // Utility function to transfer the solution step data between the modelparts.
    // Used in the python file chimera_modelpart_import.py
    void TransferSolutionStepData(ModelPart& rFromModelPart, ModelPart& rToModelPart)
    {
      rToModelPart.SetNodalSolutionStepVariablesList( rFromModelPart.pGetNodalSolutionStepVariablesList() );
    }

    namespace py = pybind11;
    void  AddCustomUtilitiesToPython(pybind11::module& m)
    {

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    // Base settings
    typedef SolverSettings<SparseSpaceType,LocalSpaceType,LinearSolverType> BaseSettingsType;
    typedef FractionalStepSettingsForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType> FractionalStepSettingsForChimeraType;

    typedef void (FractionalStepSettingsForChimeraType::*SetStrategyByParamsWithBSType)(SolverSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel const&,
                                                                                    LinearSolverType::Pointer,
                                                                                    const double,
                                                                                    const unsigned int);
    SetStrategyByParamsWithBSType ThisSetStrategyOverload = &FractionalStepSettingsForChimeraType::SetStrategy;

    py::class_< FractionalStepSettingsForChimeraType, BaseSettingsType>
        (m, "FractionalStepSettings")
        .def(py::init<ModelPart&,unsigned int,unsigned int,bool,bool,bool>())
        .def("SetStrategy",ThisSetStrategyOverload)
        .def("GetStrategy",&FractionalStepSettingsForChimeraType::pGetStrategy)
        .def("SetEchoLevel",&FractionalStepSettingsForChimeraType::SetEchoLevel)
        ;

    // Utility function to transfer the solution step data between the modelparts.
    // Used in the python file chimera_modelpart_import.py
    m.def("TransferSolutionStepData", &TransferSolutionStepData, "Utility function to transfer the solution step data between the modelparts.");
  }

}  // namespace Python.

} // Namespace Kratos
