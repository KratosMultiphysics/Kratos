//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <pybind11/stl.h>

// External includes

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

// Constitutive laws
#include "custom_constitutive/hyperelastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_U_P_3D_law.hpp"
#include "custom_constitutive/hyperelastic_U_P_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_U_P_axisym_2D_law.hpp"

#include "custom_constitutive/linear_elastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_stress_2D_law.hpp"
#include "custom_constitutive/linear_elastic_axisym_2D_law.hpp"
#include "custom_constitutive/linear_elastic_orthotropic_3D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_J2_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_J2_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_U_P_J2_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_J2_axisym_2D_law.hpp"

#include "custom_constitutive/isotropic_damage_simo_ju_3D_law.hpp"
#include "custom_constitutive/isotropic_damage_simo_ju_plane_strain_2D_law.hpp"
#include "custom_constitutive/isotropic_damage_simo_ju_plane_stress_2D_law.hpp"

#include "custom_constitutive/isotropic_damage_modified_mises_3D_law.hpp"
#include "custom_constitutive/isotropic_damage_modified_mises_plane_strain_2D_law.hpp"
#include "custom_constitutive/isotropic_damage_modified_mises_plane_stress_2D_law.hpp"

namespace Kratos
{
namespace Python
{

using namespace pybind11;

void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
  // Linear Elastic laws

  class_< LinearElastic3DLaw, typename LinearElastic3DLaw::Pointer, ConstitutiveLaw >
      (m, "LinearElastic3DLaw").def(init<>() )
      ;

  class_< LinearElasticPlaneStrain2DLaw, typename LinearElasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "LinearElasticPlaneStrain2DLaw").def(init<>() )
      ;

  class_< LinearElasticPlaneStress2DLaw, typename LinearElasticPlaneStress2DLaw::Pointer, ConstitutiveLaw >
      (m, "LinearElasticPlaneStress2DLaw").def(init<>() )
      ;

  class_< LinearElasticAxisym2DLaw, typename LinearElasticAxisym2DLaw::Pointer, ConstitutiveLaw >
      (m, "LinearElasticAxisym2DLaw").def(init<>() )
      ;

  class_< LinearElasticOrthotropic3DLaw, typename LinearElasticOrthotropic3DLaw::Pointer, ConstitutiveLaw >
      (m, "LinearElasticOrthotropic3DLaw").def(init<>() )
      ;

  // Hyperelastic laws

  class_< HyperElastic3DLaw, typename HyperElastic3DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElastic3DLaw").def(init<>() )
      ;

  class_< HyperElasticPlaneStrain2DLaw, typename HyperElasticPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticPlaneStrain2DLaw").def(init<>() )
      ;

  class_< HyperElasticAxisym2DLaw, typename HyperElasticAxisym2DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticAxisym2DLaw").def(init<>() )
      ;


  // Hyperelastic laws U-P

  class_< HyperElasticUP3DLaw, typename HyperElasticUP3DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticUP3DLaw").def(init<>() )
      ;


  class_< HyperElasticUPPlaneStrain2DLaw, typename HyperElasticUPPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticUPPlaneStrain2DLaw").def(init<>() )
      ;

  class_< HyperElasticUPAxisym2DLaw, typename HyperElasticUPAxisym2DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticUPAxisym2DLaw").def(init<>() )
      ;


  // Hyperelastic Plastic J2 specilization laws

  class_<HyperElasticPlasticJ23DLaw, typename HyperElasticPlasticJ23DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticPlasticJ23DLaw").def(init<>() )
      ;

  class_<HyperElasticPlasticJ2PlaneStrain2DLaw, typename HyperElasticPlasticJ2PlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticPlasticJ2PlaneStrain2DLaw").def(init<>() )
      ;

  class_<HyperElasticPlasticJ2Axisym2DLaw, typename HyperElasticPlasticJ2Axisym2DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticPlasticJ2Axisym2DLaw").def(init<>() )
      ;

  // Hyperelastic Plastic J2 specilization laws U-P

  class_<HyperElasticPlasticUPJ23DLaw, typename HyperElasticPlasticUPJ23DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticPlasticUPJ23DLaw").def(init<>() )
      ;

  class_<HyperElasticPlasticUPJ2PlaneStrain2DLaw, typename HyperElasticPlasticUPJ2PlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticPlasticUPJ2PlaneStrain2DLaw").def(init<>() )
      ;

  class_<HyperElasticPlasticUPJ2Axisym2DLaw, typename HyperElasticPlasticUPJ2Axisym2DLaw::Pointer, ConstitutiveLaw >
      (m, "HyperElasticPlasticUPJ2Axisym2DLaw").def(init<>() )
      ;


  // Isotropic Damage laws

  class_<IsotropicDamageSimoJu3DLaw, typename IsotropicDamageSimoJu3DLaw::Pointer, ConstitutiveLaw >
      (m, "IsotropicDamageSimoJu3DLaw").def(init<>() )
      ;

  class_<IsotropicDamageSimoJuPlaneStrain2DLaw, typename IsotropicDamageSimoJuPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "IsotropicDamageSimoJuPlaneStrain2DLaw").def(init<>() )
      ;

  class_<IsotropicDamageSimoJuPlaneStress2DLaw, typename IsotropicDamageSimoJuPlaneStress2DLaw::Pointer, ConstitutiveLaw >
      (m, "IsotropicDamageSimoJuPlaneStress2DLaw").def(init<>() )
      ;

  class_<IsotropicDamageModifiedMises3DLaw, typename IsotropicDamageModifiedMises3DLaw::Pointer, ConstitutiveLaw >
      (m, "IsotropicDamageModifiedMises3DLaw").def(init<>() )
      ;

  class_<IsotropicDamageModifiedMisesPlaneStrain2DLaw, typename IsotropicDamageModifiedMisesPlaneStrain2DLaw::Pointer, ConstitutiveLaw >
      (m, "IsotropicDamageModifiedMisesPlaneStrain2DLaw").def(init<>() )
      ;

  class_<IsotropicDamageModifiedMisesPlaneStress2DLaw, typename IsotropicDamageModifiedMisesPlaneStress2DLaw::Pointer, ConstitutiveLaw >
      (m, "IsotropicDamageModifiedMisesPlaneStress2DLaw").def(init<>() )
      ;

}

}  // namespace Python.
}  // namespace Kratos.
