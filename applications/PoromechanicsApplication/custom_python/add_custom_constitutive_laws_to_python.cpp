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

#include "custom_constitutive/restore_simo_ju_3D_law.hpp"
#include "custom_constitutive/restore_simo_ju_plane_strain_2D_law.hpp"
#include "custom_constitutive/restore_simo_ju_plane_stress_2D_law.hpp"

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
}

}  // namespace Python.
}  // namespace Kratos.
