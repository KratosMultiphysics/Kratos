//
//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:                 $
//   Revision:            $Revision:          $
//

// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"

#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"
#include "python/add_mesh_to_python.h"

//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//constitutive laws
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress.hpp"

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

using namespace boost::python;

void  AddCustomConstitutiveLawsToPython()
{
    class_< ThermalLinearElastic3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalLinearElastic3DLaw",init<>() );

    class_< ThermalLinearElastic2DPlaneStrain, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalLinearElastic2DPlaneStrain",init<>() );

    class_< ThermalLinearElastic2DPlaneStress, bases< ConstitutiveLaw >, boost::noncopyable >( "ThermalLinearElastic2DPlaneStress",init<>() );

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
