// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: IgaApplication/license.txt
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