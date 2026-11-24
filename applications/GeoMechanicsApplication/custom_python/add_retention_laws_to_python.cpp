// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Richard Faasse
//

#include "custom_python/add_retention_laws_to_python.h"

#include "custom_retention/retention_law.h"
#include "custom_retention/retention_law_factory.h"
#include "custom_retention/saturated_below_phreatic_level_law.h"
#include "custom_retention/saturated_law.h"
#include "custom_retention/van_genuchten_law.h"

namespace Kratos::Python
{

void AddRetentionLawsToPython(const pybind11::module& rModule)
{
    pybind11::class_<RetentionLaw, RetentionLaw::Pointer>(rModule, "RetentionLaw", pybind11::module_local())
        .def("CalculateSaturation", &RetentionLaw::CalculateSaturation)
        .def("CalculateEffectiveSaturation", &RetentionLaw::CalculateEffectiveSaturation)
        .def("CalculateDerivativeOfSaturation", &RetentionLaw::CalculateDerivativeOfSaturation)
        .def("CalculateRelativePermeability", &RetentionLaw::CalculateRelativePermeability);

    pybind11::class_<SaturatedLaw, SaturatedLaw::Pointer>(rModule, "SaturatedLaw", pybind11::module_local())
        .def(pybind11::init<>())
        .def("CalculateSaturation", &SaturatedLaw::CalculateSaturation)
        .def("CalculateEffectiveSaturation", &SaturatedLaw::CalculateEffectiveSaturation)
        .def("CalculateDerivativeOfSaturation", &SaturatedLaw::CalculateDerivativeOfSaturation)
        .def("CalculateRelativePermeability", &SaturatedLaw::CalculateRelativePermeability);

    pybind11::class_<SaturatedBelowPhreaticLevelLaw, SaturatedBelowPhreaticLevelLaw::Pointer>(
        rModule, "SaturatedBelowPhreaticLevelLaw", pybind11::module_local())
        .def(pybind11::init<>())
        .def("CalculateSaturation", &SaturatedBelowPhreaticLevelLaw::CalculateSaturation)
        .def("CalculateEffectiveSaturation", &SaturatedBelowPhreaticLevelLaw::CalculateEffectiveSaturation)
        .def("CalculateDerivativeOfSaturation", &SaturatedBelowPhreaticLevelLaw::CalculateDerivativeOfSaturation)
        .def("CalculateRelativePermeability", &SaturatedBelowPhreaticLevelLaw::CalculateRelativePermeability);

    pybind11::class_<VanGenuchtenLaw, VanGenuchtenLaw::Pointer>(rModule, "VanGenuchtenLaw",
                                                                pybind11::module_local())
        .def(pybind11::init<>())
        .def("CalculateSaturation", &VanGenuchtenLaw::CalculateSaturation)
        .def("CalculateEffectiveSaturation", &VanGenuchtenLaw::CalculateEffectiveSaturation)
        .def("CalculateDerivativeOfSaturation", &VanGenuchtenLaw::CalculateDerivativeOfSaturation)
        .def("CalculateRelativePermeability", &VanGenuchtenLaw::CalculateRelativePermeability);

    pybind11::class_<RetentionLaw::Parameters, RetentionLaw::Parameters::Pointer>(
        rModule, "RetentionLawParameters", pybind11::module_local())
        .def(pybind11::init<const Properties&>())
        .def("SetFluidPressure", &RetentionLaw::Parameters::SetFluidPressure)
        .def("GetFluidPressure", &RetentionLaw::Parameters::GetFluidPressure);
}

} // Namespace Kratos::Python.
