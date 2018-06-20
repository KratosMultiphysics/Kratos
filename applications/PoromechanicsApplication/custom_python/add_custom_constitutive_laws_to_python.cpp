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

using namespace pybind11;

void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    class_< BilinearCohesive3DLaw, BilinearCohesive3DLaw::Pointer, ConstitutiveLaw >
    (m, "BilinearCohesive3DLaw")
    .def( init<>() );
    class_< BilinearCohesive2DLaw, BilinearCohesive2DLaw::Pointer, ConstitutiveLaw >
    (m, "BilinearCohesive2DLaw")
    .def( init<>() ) ;

    class_< SimoJuLocalDamage3DLaw, SimoJuLocalDamage3DLaw::Pointer, ConstitutiveLaw >
    (m, "SimoJuLocalDamage3DLaw")
    .def( init<>() );
    class_< SimoJuLocalDamagePlaneStrain2DLaw, SimoJuLocalDamagePlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "SimoJuLocalDamagePlaneStrain2DLaw")
    .def( init<>() );
    class_< SimoJuLocalDamagePlaneStress2DLaw, SimoJuLocalDamagePlaneStress2DLaw::Pointer, ConstitutiveLaw >
    (m, "SimoJuLocalDamagePlaneStress2DLaw")
    .def( init<>() );

    class_< SimoJuNonlocalDamage3DLaw, SimoJuNonlocalDamage3DLaw::Pointer, ConstitutiveLaw >
    (m, "SimoJuNonlocalDamage3DLaw")
    .def( init<>() );
    class_< SimoJuNonlocalDamagePlaneStrain2DLaw, SimoJuNonlocalDamagePlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "SimoJuNonlocalDamagePlaneStrain2DLaw")
    .def( init<>() );
    class_< SimoJuNonlocalDamagePlaneStress2DLaw, SimoJuNonlocalDamagePlaneStress2DLaw::Pointer, ConstitutiveLaw >
    (m, "SimoJuNonlocalDamagePlaneStress2DLaw")
    .def( init<>() );

    class_< ModifiedMisesNonlocalDamage3DLaw, ModifiedMisesNonlocalDamage3DLaw::Pointer, ConstitutiveLaw >
    (m, "ModifiedMisesNonlocalDamage3DLaw")
    .def( init<>() );
    class_< ModifiedMisesNonlocalDamagePlaneStrain2DLaw, ModifiedMisesNonlocalDamagePlaneStrain2DLaw::Pointer, ConstitutiveLaw >
    (m, "ModifiedMisesNonlocalDamagePlaneStrain2DLaw")
    .def( init<>() );
    class_< ModifiedMisesNonlocalDamagePlaneStress2DLaw, ModifiedMisesNonlocalDamagePlaneStress2DLaw::Pointer, ConstitutiveLaw >
    (m, "ModifiedMisesNonlocalDamagePlaneStress2DLaw")
    .def( init<>() );
}

}  // namespace Python.
}  // namespace Kratos.
