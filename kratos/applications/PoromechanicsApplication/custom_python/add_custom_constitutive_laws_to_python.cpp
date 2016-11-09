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

#include "custom_constitutive/custom_flow_rules/restore_damage_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/restore_nonlocal_damage_flow_rule.hpp"

//constitutive laws
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"
#include "custom_constitutive/bilinear_cohesive_2D_law.hpp"
#include "custom_constitutive/restore_simo_ju_3D_law.hpp"
#include "custom_constitutive/restore_simo_ju_plane_strain_2D_law.hpp"
#include "custom_constitutive/restore_simo_ju_plane_stress_2D_law.hpp"
#include "custom_constitutive/restore_simo_ju_nonlocal_3D_law.hpp"
#include "custom_constitutive/restore_simo_ju_nonlocal_plane_strain_2D_law.hpp"
#include "custom_constitutive/restore_simo_ju_nonlocal_plane_stress_2D_law.hpp"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void  AddCustomConstitutiveLawsToPython()
{        
    class_< BilinearCohesive3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "BilinearCohesive3DLaw",init<>() );
    class_< BilinearCohesive2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "BilinearCohesive2DLaw",init<>() );
    
    class_< RestoreSimoJu3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "RestoreSimoJu3DLaw",init<>() );
    class_< RestoreSimoJuPlaneStrain2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "RestoreSimoJuPlaneStrain2DLaw",init<>() );
    class_< RestoreSimoJuPlaneStress2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "RestoreSimoJuPlaneStress2DLaw",init<>() );

    class_< RestoreSimoJuNonlocal3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "RestoreSimoJuNonlocal3DLaw",init<>() );
    class_< RestoreSimoJuNonlocalPlaneStrain2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "RestoreSimoJuNonlocalPlaneStrain2DLaw",init<>() );
    class_< RestoreSimoJuNonlocalPlaneStress2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "RestoreSimoJuNonlocalPlaneStress2DLaw",init<>() );
}

}  // namespace Python.
}  // namespace Kratos.
