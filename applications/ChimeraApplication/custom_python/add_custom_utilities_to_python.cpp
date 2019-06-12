//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//                   Navaneeth K Narayanan
//
// ==============================================================================
// System includes

// External includes


// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/solver_settings.h"
#include "custom_utilities/fractional_step_settings_for_chimera.h"

//Processes
namespace Kratos
{

namespace Python
{
    using namespace pybind11;
    void  AddCustomUtilitiesToPython(pybind11::module& m)
    {

 // Base settings
    typedef SolverSettings<SparseSpaceType,LocalSpaceType,LinearSolverType> BaseSettingsType;

    class_ < BaseSettingsType >(m, "BaseSettingsType");

    // Fractional step settings
    enum_<FractionalStepSettingsForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel>(m,"ChimeraStrategyLabel")
        .value("Velocity",FractionalStepSettingsForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::Velocity)
        .value("Pressure",FractionalStepSettingsForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::Pressure)
        //.value("EddyViscosity",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::EddyViscosity)
        ;
    ;

    /* enum_<FractionalStepSettingsForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel>("TurbulenceModelLabel")
        .value("SpalartAllmaras",FractionalStepSettingsForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::SpalartAllmaras)
    ;
 */
    typedef void (FractionalStepSettingsForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::*SetStrategyByParamsType)(FractionalStepSettingsForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel const&,LinearSolverType::Pointer,const double,const std::size_t);
    SetStrategyByParamsType ThisSetStrategyOverload = &FractionalStepSettingsForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy;

    class_< FractionalStepSettingsForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>, BaseSettingsType>
        (m, "FractionalStepSettingsForChimera")
        .def(init<ModelPart&,unsigned int,unsigned int,bool,bool,bool>())
        .def("SetStrategy",ThisSetStrategyOverload)
        //.def("SetTurbulenceModel",SetTurbModel_Build)
        //.def("SetTurbulenceModel",SetTurbModel_Pass)
        .def("GetStrategy",&FractionalStepSettingsForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
        .def("SetEchoLevel",&FractionalStepSettingsForChimera<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
    ;
  }


    }
}  // namespace Python.

} // Namespace Kratos
