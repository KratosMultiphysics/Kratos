//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Max Friedrichs Dachale, Juan Ignacio Camarotti,
//                   Ricky Aristio, Alicia Knauer
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

// Constitutive laws
#include "custom_constitutive/bernoulli_beam_elastic_constitutive_law.h"

namespace Kratos::Python {

void AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_< BernoulliBeamElasticConstitutiveLaw, typename BernoulliBeamElasticConstitutiveLaw::Pointer, ConstitutiveLaw >
    (m, "BernoulliBeamElasticConstitutiveLaw").def(py::init<>() )
    ;
}

}  // namespace Kratos::Python.
