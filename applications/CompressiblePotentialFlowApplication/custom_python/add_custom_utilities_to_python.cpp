//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"

//Utilities
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto mod_potential_flow_utilities = m.def_submodule("PotentialFlowUtilities");

    mod_potential_flow_utilities.def("CheckIfWakeConditionsAreFulfilled2D",&PotentialFlowUtilities::CheckIfWakeConditionsAreFulfilled<2>);
    mod_potential_flow_utilities.def("CheckIfWakeConditionsAreFulfilled3D",&PotentialFlowUtilities::CheckIfWakeConditionsAreFulfilled<3>);
    mod_potential_flow_utilities.def("CalculateArea",&PotentialFlowUtilities::CalculateArea<ModelPart::ElementsContainerType>);
    mod_potential_flow_utilities.def("CalculateArea",&PotentialFlowUtilities::CalculateArea<ModelPart::ConditionsContainerType>);
    mod_potential_flow_utilities.def("ComputePotentialJump2D",&PotentialFlowUtilities::ComputePotentialJump<2,3>);
    mod_potential_flow_utilities.def("ComputePotentialJump3D",&PotentialFlowUtilities::ComputePotentialJump<3,4>);
}

}  // namespace Python.
} // Namespace Kratos
