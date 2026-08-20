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

// External includes

// Project includes
#include "custom_python/add_retention_laws_to_python.h"

#include "../custom_retention/retention_law.h"
#include "../custom_retention/van_genuchten_law.h"
#include "custom_retention/retention_law.h"
#include "custom_retention/retention_law_factory.h"

namespace Kratos::Python
{

void AddRetentionLawsToPython(const pybind11::module& rModule)
{
    pybind11::class_<RetentionLawFactory, RetentionLawFactory::Pointer>(
        rModule, "RetentionLawFactory", pybind11::module_local())
        .def("Clone", &RetentionLawFactory::Clone);

    pybind11::class_<VanGenuchtenLaw, VanGenuchtenLaw::Pointer>(rModule, "VanGenuchtenLaw",
                                                                pybind11::module_local())
        .def(pybind11::init<>())
        .def("CalculateSaturation", &VanGenuchtenLaw::CalculateSaturation)
        .def("CalculateEffectiveSaturation", &VanGenuchtenLaw::CalculateSaturation)
        .def("CalculateDerivativeOfSaturation", &VanGenuchtenLaw::CalculateDerivativeOfSaturation)
        .def("CalculateRelativePermeability", &VanGenuchtenLaw::CalculateRelativePermeability);

    pybind11::class_<RetentionLaw::Parameters, RetentionLaw::Parameters::Pointer>(
        rModule, "RetentionLawParameters", pybind11::module_local())
        .def(pybind11::init<const Properties&>())
        .def("SetFluidPressure", &RetentionLaw::Parameters::SetFluidPressure)
        .def("GetFluidPressure", &RetentionLaw::Parameters::GetFluidPressure);
}

} // Namespace Kratos::Python.
