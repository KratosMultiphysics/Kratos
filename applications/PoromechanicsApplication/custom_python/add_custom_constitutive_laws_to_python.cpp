//
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// System includes

// Project includes
#include "includes/constitutive_law.h"

//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//constitutive laws
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"
#include "custom_constitutive/bilinear_cohesive_2D_law.hpp"
#include "custom_constitutive/exponential_cohesive_3D_law.hpp"
#include "custom_constitutive/exponential_cohesive_2D_law.hpp"

#include "custom_constitutive/simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/simo_ju_local_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/simo_ju_local_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/simo_ju_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/simo_ju_nonlocal_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/modified_mises_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/modified_mises_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/modified_mises_nonlocal_damage_plane_stress_2D_law.hpp"
#include "custom_constitutive/history_linear_elastic_3D_law.hpp"
#include "custom_constitutive/history_linear_elastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/history_linear_elastic_plane_stress_2D_law.hpp"

#include "custom_constitutive/hyperelastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_stress_2D_law.hpp"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    py::class_< BilinearCohesive3DLaw, BilinearCohesive3DLaw::Pointer, ConstitutiveLaw >
    (m, "BilinearCohesive3DLaw")
    .def( py::init<>() );
    py::class_< BilinearCohesive2DLaw, BilinearCohesive2DLaw::Pointer, ConstitutiveLaw >
    (m, "BilinearCohesive2DLaw")
    .def( py::init<>() ) ;
    py::class_< ExponentialCohesive3DLaw, ExponentialCohesive3DLaw::Pointer, ConstitutiveLaw >
    (m, "ExponentialCohesive3DLaw")
    .def( py::init<>() );
    py::class_< ExponentialCohesive2DLaw, ExponentialCohesive2DLaw::Pointer, ConstitutiveLaw >
    (m, "ExponentialCohesive2DLaw")
    .def( py::init<>() );

    py::class_< LinearElastic3DLaw, LinearElastic3DLaw::Pointer, ConstitutiveLaw >
    (m, "LinearElasticSolid3DLaw")
    .def( py::init<>() );
    py::class_< LinearElasticPlaneStrain2DLaw, LinearElasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "LinearElasticPlaneStrainSolid2DLaw")
    .def(py::init<>() );
    py::class_< LinearElasticPlaneStress2DLaw, LinearElasticPlaneStress2DLaw::Pointer, ConstitutiveLaw >
    (m, "LinearElasticPlaneStressSolid2DLaw")
    .def(py::init<>() );
    py::class_< HyperElastic3DLaw, HyperElastic3DLaw::Pointer, ConstitutiveLaw >
    (m, "HyperElasticSolid3DLaw")
    .def( py::init<>() );

    py::class_< SimoJuLocalDamage3DLaw, SimoJuLocalDamage3DLaw::Pointer, ConstitutiveLaw >
    (m, "SimoJuLocalDamage3DLaw")
    .def( py::init<>() );
    py::class_< SimoJuLocalDamagePlaneStrain2DLaw, SimoJuLocalDamagePlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "SimoJuLocalDamagePlaneStrain2DLaw")
    .def( py::init<>() );
    py::class_< SimoJuLocalDamagePlaneStress2DLaw, SimoJuLocalDamagePlaneStress2DLaw::Pointer, ConstitutiveLaw >
    (m, "SimoJuLocalDamagePlaneStress2DLaw")
    .def( py::init<>() );

    py::class_< SimoJuNonlocalDamage3DLaw, SimoJuNonlocalDamage3DLaw::Pointer, ConstitutiveLaw >
    (m, "SimoJuNonlocalDamage3DLaw")
    .def( py::init<>() );
    py::class_< SimoJuNonlocalDamagePlaneStrain2DLaw, SimoJuNonlocalDamagePlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "SimoJuNonlocalDamagePlaneStrain2DLaw")
    .def( py::init<>() );
    py::class_< SimoJuNonlocalDamagePlaneStress2DLaw, SimoJuNonlocalDamagePlaneStress2DLaw::Pointer, ConstitutiveLaw >
    (m, "SimoJuNonlocalDamagePlaneStress2DLaw")
    .def( py::init<>() );

    py::class_< ModifiedMisesNonlocalDamage3DLaw, ModifiedMisesNonlocalDamage3DLaw::Pointer, ConstitutiveLaw >
    (m, "ModifiedMisesNonlocalDamage3DLaw")
    .def( py::init<>() );
    py::class_< ModifiedMisesNonlocalDamagePlaneStrain2DLaw, ModifiedMisesNonlocalDamagePlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "ModifiedMisesNonlocalDamagePlaneStrain2DLaw")
    .def( py::init<>() );
    py::class_< ModifiedMisesNonlocalDamagePlaneStress2DLaw, ModifiedMisesNonlocalDamagePlaneStress2DLaw::Pointer, ConstitutiveLaw >
    (m, "ModifiedMisesNonlocalDamagePlaneStress2DLaw")
    .def( py::init<>() );

    py::class_< HistoryLinearElastic3DLaw, HistoryLinearElastic3DLaw::Pointer, ConstitutiveLaw >
    (m, "HistoryLinearElastic3DLaw")
    .def( py::init<>() );
    py::class_< HistoryLinearElasticPlaneStrain2DLaw, HistoryLinearElasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "HistoryLinearElasticPlaneStrain2DLaw")
    .def( py::init<>() );
    py::class_< HistoryLinearElasticPlaneStress2DLaw, HistoryLinearElasticPlaneStress2DLaw::Pointer, ConstitutiveLaw >
    (m, "HistoryLinearElasticPlaneStress2DLaw")
    .def( py::init<>() );
}

}  // namespace Python.
}  // namespace Kratos.
