//
//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:               $
//   Revision:            $Revision:           $
//

// System includes

// Project includes
#include "includes/constitutive_law.h"

//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//constitutive laws
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress.hpp"

#include "custom_constitutive/linear_elastic_3D_law_nodal.hpp"
#include "custom_constitutive/linear_elastic_2D_plane_strain_nodal.hpp"
#include "custom_constitutive/linear_elastic_2D_plane_stress_nodal.hpp"

#include "custom_constitutive/thermal_linear_elastic_3D_law_nodal.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain_nodal.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress_nodal.hpp"

#include "custom_constitutive/thermal_simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/thermal_modified_mises_nonlocal_damage_plane_stress_2D_law.hpp"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    class_< ThermalLinearElastic3DLaw, ThermalLinearElastic3DLaw::Pointer, ConstitutiveLaw >
    (m, "ThermalLinearElastic3DLaw")
    .def(init<>());
    class_< ThermalLinearElastic2DPlaneStrain, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalLinearElastic2DPlaneStrain",init<>() );
    class_< ThermalLinearElastic2DPlaneStress, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalLinearElastic2DPlaneStress",init<>() );

    class_< LinearElastic3DLawNodal, bases< ConstitutiveLaw >, boost::noncopyable >( "LinearElastic3DLawNodal",init<>() );
    class_< LinearElastic2DPlaneStrainNodal, bases< ConstitutiveLaw >, boost::noncopyable >( "LinearElastic2DPlaneStrainNodal",init<>() );
    class_< LinearElastic2DPlaneStressNodal, bases< ConstitutiveLaw >, boost::noncopyable >( "LinearElastic2DPlaneStressNodal",init<>() );

    class_< ThermalLinearElastic3DLawNodal, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalLinearElastic3DLawNodal",init<>() );
    class_< ThermalLinearElastic2DPlaneStrainNodal, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalLinearElastic2DPlaneStrainNodal",init<>() );
    class_< ThermalLinearElastic2DPlaneStressNodal, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalLinearElastic2DPlaneStressNodal",init<>() );

    class_< ThermalSimoJuLocalDamage3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalSimoJuLocalDamage3DLaw",init<>() );
    class_< ThermalSimoJuLocalDamagePlaneStrain2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalSimoJuLocalDamagePlaneStrain2DLaw",init<>() );
    class_< ThermalSimoJuLocalDamagePlaneStress2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalSimoJuLocalDamagePlaneStress2DLaw",init<>() );

    class_< ThermalSimoJuNonlocalDamage3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalSimoJuNonlocalDamage3DLaw",init<>() );
    class_< ThermalSimoJuNonlocalDamagePlaneStrain2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalSimoJuNonlocalDamagePlaneStrain2DLaw",init<>() );
    class_< ThermalSimoJuNonlocalDamagePlaneStress2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalSimoJuNonlocalDamagePlaneStress2DLaw",init<>() );

    class_< ThermalModifiedMisesNonlocalDamage3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalModifiedMisesNonlocalDamage3DLaw",init<>() );
    class_< ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw",init<>() );
    class_< ThermalModifiedMisesNonlocalDamagePlaneStress2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalModifiedMisesNonlocalDamagePlaneStress2DLaw",init<>() );
}

}  // namespace Python.
}  // namespace Kratos.
