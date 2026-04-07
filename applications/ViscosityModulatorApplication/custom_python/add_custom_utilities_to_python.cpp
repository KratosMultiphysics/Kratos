// KRATOS
//  _   _            _   _            _
// | | | (_)        |  \/  |         | |
// | | | |_ ___  ___| .  . | ___   __| |
// | | | | / __|/ __| |\/| |/ _ \ / _` |
// \ \_/ / \__ \ (__| |  | | (_) | (_| |
//  \___/|_|___/\___\_|  |_/\___/ \__,_|  APPLICATION
//                                      
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/viscosity_modulator_coupling_utilities.h"



namespace Kratos
{
namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // py::class_<Variable<ViscosityModulatorSettings::Pointer>, VariableData>(m, "ViscosityModulatorSettingsVariable")
    // .def("__str__", PrintObject<Variable<ViscosityModulatorSettings::Pointer>>)
    // ;

    py::class_<ViscosityModulatorCouplingUtilities, ViscosityModulatorCouplingUtilities::Pointer>(m, "ViscosityModulatorCouplingUtilities")
    .def(py::init<>())
    .def_static("ComputeRelativeResidual", &ViscosityModulatorCouplingUtilities::ComputeRelativeResidual)
    .def_static("ApplyRelaxation", &ViscosityModulatorCouplingUtilities::ApplyRelaxation)
    .def_static("ComputeQuasiNewtonUpdateVectors", &ViscosityModulatorCouplingUtilities::ComputeQuasiNewtonUpdateVectors)
    .def_static("UpdateConvergenceVariables", &ViscosityModulatorCouplingUtilities::UpdateConvergenceVariables)
    ;

    // py::class_< ViscosityModulatorSettings, ViscosityModulatorSettings::Pointer >	(m,"ViscosityModulatorSettings")
    // .def(py::init<	>() )
    // .def("SetDensityVariable",&ViscosityModulatorSettings::SetDensityVariable)
    // .def("SetDiffusionVariable",&ViscosityModulatorSettings::SetDiffusionVariable)
    // .def("SetUnknownVariable",&ViscosityModulatorSettings::SetUnknownVariable)
    // .def("SetVolumeSourceVariable",&ViscosityModulatorSettings::SetVolumeSourceVariable)
    // .def("SetSurfaceSourceVariable",&ViscosityModulatorSettings::SetSurfaceSourceVariable)
    // .def("SetProjectionVariable",&ViscosityModulatorSettings::SetProjectionVariable)
    // .def("SetMeshVelocityVariable",&ViscosityModulatorSettings::SetMeshVelocityVariable)
    // .def("SetConvectionVariable",&ViscosityModulatorSettings::SetConvectionVariable)
    // .def("SetGradientVariable",&ViscosityModulatorSettings::SetGradientVariable)
    // .def("SetTransferCoefficientVariable",&ViscosityModulatorSettings::SetTransferCoefficientVariable)
    // .def("SetSpecificHeatVariable",&ViscosityModulatorSettings::SetSpecificHeatVariable)
    // .def("SetVelocityVariable",&ViscosityModulatorSettings::SetVelocityVariable)
    // .def("SetReactionVariable",&ViscosityModulatorSettings::SetReactionVariable)
    // .def("SetReactionGradientVariable",&ViscosityModulatorSettings::SetReactionGradientVariable)
    //     .def("GetDensityVariable",&ViscosityModulatorSettings::GetDensityVariable, py::return_value_policy::reference_internal )
    // .def("GetDiffusionVariable",&ViscosityModulatorSettings::GetDiffusionVariable, py::return_value_policy::reference_internal )
    // .def("GetUnknownVariable",&ViscosityModulatorSettings::GetUnknownVariable, py::return_value_policy::reference_internal )
    // .def("GetVolumeSourceVariable",&ViscosityModulatorSettings::GetVolumeSourceVariable, py::return_value_policy::reference_internal )
    // .def("GetSurfaceSourceVariable",&ViscosityModulatorSettings::GetSurfaceSourceVariable, py::return_value_policy::reference_internal )
    // .def("GetProjectionVariable",&ViscosityModulatorSettings::GetProjectionVariable, py::return_value_policy::reference_internal )
    // .def("GetMeshVelocityVariable",&ViscosityModulatorSettings::GetMeshVelocityVariable, py::return_value_policy::reference_internal )
    // .def("GetConvectionVariable",&ViscosityModulatorSettings::GetConvectionVariable, py::return_value_policy::reference_internal )
    // .def("GetGradientVariable",&ViscosityModulatorSettings::GetGradientVariable, py::return_value_policy::reference_internal )
    // .def("GetTransferCoefficientVariable",&ViscosityModulatorSettings::GetTransferCoefficientVariable, py::return_value_policy::reference_internal)
    // .def("GetSpecificHeatVariable",&ViscosityModulatorSettings::GetSpecificHeatVariable, py::return_value_policy::reference_internal )
    // .def("GetVelocityVariable",&ViscosityModulatorSettings::GetVelocityVariable, py::return_value_policy::reference_internal )
    // .def("GetReactionVariable",&ViscosityModulatorSettings::GetReactionVariable, py::return_value_policy::reference_internal )
    // .def("GetReactionGradientVariable",&ViscosityModulatorSettings::GetReactionGradientVariable, py::return_value_policy::reference_internal )
    //     .def("IsDefinedDensityVariable",&ViscosityModulatorSettings::IsDefinedDensityVariable)
    // .def("IsDefinedDiffusionVariable",&ViscosityModulatorSettings::IsDefinedDiffusionVariable)
    // .def("IsDefinedUnknownVariable",&ViscosityModulatorSettings::IsDefinedUnknownVariable)
    // .def("IsDefinedVolumeSourceVariable",&ViscosityModulatorSettings::IsDefinedVolumeSourceVariable)
    // .def("IsDefinedSurfaceSourceVariable",&ViscosityModulatorSettings::IsDefinedSurfaceSourceVariable)
    // .def("IsDefinedProjectionVariable",&ViscosityModulatorSettings::IsDefinedProjectionVariable)
    // .def("IsDefinedMeshVelocityVariable",&ViscosityModulatorSettings::IsDefinedMeshVelocityVariable)
    // .def("IsDefinedConvectionVariable",&ViscosityModulatorSettings::IsDefinedConvectionVariable)
    // .def("IsDefinedGradientVariable",&ViscosityModulatorSettings::IsDefinedGradientVariable)
    // .def("IsDefinedSpecificHeatVariable",&ViscosityModulatorSettings::IsDefinedSpecificHeatVariable)
    // .def("IsDefinedVelocityVariable",&ViscosityModulatorSettings::IsDefinedVelocityVariable)
    // .def("IsDefinedTransferCoefficientVariable",&ViscosityModulatorSettings::IsDefinedTransferCoefficientVariable)
    // .def("IsDefinedReactionVariable",&ViscosityModulatorSettings::IsDefinedReactionVariable)
    // .def("IsDefinedReactionGradientVariable",&ViscosityModulatorSettings::IsDefinedReactionGradientVariable)
    // ;

}

}  // namespace Python.
} // Namespace Kratos
