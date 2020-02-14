// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

// Elastic laws
#include "custom_constitutive/truss_constitutive_law.h"
#include "custom_constitutive/beam_constitutive_law.h"
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/axisym_elastic_isotropic.h"
#include "custom_constitutive/linear_plane_stress.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/user_provided_linear_elastic_law.h"

namespace Kratos {
namespace Python {

void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_< TrussConstitutiveLaw, typename TrussConstitutiveLaw::Pointer, ConstitutiveLaw >
    (m, "TrussConstitutiveLaw").def(py::init<>() )
    ;

    py::class_< BeamConstitutiveLaw, typename BeamConstitutiveLaw::Pointer, ConstitutiveLaw >
    (m, "BeamConstitutiveLaw").def(py::init<>() )
    ;

    py::class_< LinearPlaneStress, typename LinearPlaneStress::Pointer, ConstitutiveLaw >
    (m, "LinearElasticPlaneStress2DLaw").def(py::init<>() )
    ;

    py::class_< LinearPlaneStrain, typename LinearPlaneStrain::Pointer, ConstitutiveLaw >
    (m, "LinearElasticPlaneStrain2DLaw").def(py::init<>() )
    ;

    py::class_< ElasticIsotropic3D, typename ElasticIsotropic3D::Pointer, ConstitutiveLaw >
    (m, "LinearElastic3DLaw").def(py::init<>() )
    ;

    py::class_< AxisymElasticIsotropic, typename AxisymElasticIsotropic::Pointer, ConstitutiveLaw >
    (m, "LinearElasticAxisym2DLaw").def(py::init<>() )
    ;

    py::class_< UserProvidedLinearElasticLaw<2>, typename UserProvidedLinearElasticLaw<2>::Pointer, ConstitutiveLaw >
    (m, "UserProvidedLinearElastic2DLaw").def(py::init<>() )
    ;

    py::class_< UserProvidedLinearElasticLaw<3>, typename UserProvidedLinearElasticLaw<3>::Pointer, ConstitutiveLaw >
    (m, "UserProvidedLinearElastic3DLaw").def(py::init<>() )
    ;
}

}  // namespace Python.
}  // namespace Kratos.
