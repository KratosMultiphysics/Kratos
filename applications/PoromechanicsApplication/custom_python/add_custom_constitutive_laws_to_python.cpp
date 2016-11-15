//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// System includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"

//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//constitutive laws
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"
#include "custom_constitutive/bilinear_cohesive_2D_law.hpp"

#include "custom_constitutive/simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/simo_ju_local_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/simo_ju_local_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/simo_ju_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/simo_ju_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/simo_ju_nonlocal_damage_plane_stress_2D_law.hpp"

#include "custom_constitutive/modified_mises_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/modified_mises_nonlocal_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/modified_mises_nonlocal_damage_plane_stress_2D_law.hpp"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void  AddCustomConstitutiveLawsToPython()
{        
    class_< BilinearCohesive3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "BilinearCohesive3DLaw",init<>() );
    class_< BilinearCohesive2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "BilinearCohesive2DLaw",init<>() );
    
    class_< SimoJuLocalDamage3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "SimoJuLocalDamage3DLaw",init<>() );
    class_< SimoJuLocalDamagePlaneStrain2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "SimoJuLocalDamagePlaneStrain2DLaw",init<>() );
    class_< SimoJuLocalDamagePlaneStress2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "SimoJuLocalDamagePlaneStress2DLaw",init<>() );

    class_< SimoJuNonlocalDamage3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "SimoJuNonlocalDamage3DLaw",init<>() );
    class_< SimoJuNonlocalDamagePlaneStrain2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "SimoJuNonlocalDamagePlaneStrain2DLaw",init<>() );
    class_< SimoJuNonlocalDamagePlaneStress2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "SimoJuNonlocalDamagePlaneStress2DLaw",init<>() );

    class_< ModifiedMisesNonlocalDamage3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ModifiedMisesNonlocalDamage3DLaw",init<>() );
    class_< ModifiedMisesNonlocalDamagePlaneStrain2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ModifiedMisesNonlocalDamagePlaneStrain2DLaw",init<>() );
    class_< ModifiedMisesNonlocalDamagePlaneStress2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "ModifiedMisesNonlocalDamagePlaneStress2DLaw",init<>() );
}

}  // namespace Python.
}  // namespace Kratos.
