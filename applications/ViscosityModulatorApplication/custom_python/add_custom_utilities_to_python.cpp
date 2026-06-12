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
#include "custom_utilities/vm_compute_flux_vector.h"



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

    py::class_<ViscosityModulatorCouplingUtilities, ViscosityModulatorCouplingUtilities::Pointer>(m, "BoussinesqCouplingUtilities")
    .def(py::init<>())
    .def_static("ComputeRelativeResidual", &ViscosityModulatorCouplingUtilities::ComputeRelativeResidual)
    .def_static("ApplyRelaxation", &ViscosityModulatorCouplingUtilities::ApplyRelaxation)
    .def_static("ComputeQuasiNewtonUpdateVectors", &ViscosityModulatorCouplingUtilities::ComputeQuasiNewtonUpdateVectors)
    .def_static("UpdateConvergenceVariables", &ViscosityModulatorCouplingUtilities::UpdateConvergenceVariables)
    ;

    py::class_<VmComputeFluxUtility, VmComputeFluxUtility::Pointer>(m, "ComputeFluxUtility")
    .def(py::init<>())
    .def_static("ComputeVectorialFlux", &VmComputeFluxUtility::ComputeVectorialFlux)
    ;

}

}  // namespace Python.
} // Namespace Kratos
