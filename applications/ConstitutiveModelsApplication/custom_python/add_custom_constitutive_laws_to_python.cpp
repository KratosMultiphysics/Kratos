//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <pybind11/stl.h>

// External includes

// Project includes
#include "includes/constitutive_law.h"
//#include "python/pointer_vector_set_python_interface.h"
//#include "python/variable_indexing_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"

// Outfitted python laws
//#include "custom_python/python_outfitted_constitutive_law.hpp"


// Constitutive laws

// Small strain laws
#include "custom_laws/small_strain_laws/small_strain_orthotropic_3D_law.hpp"
#include "custom_laws/small_strain_laws/small_strain_plane_strain_2D_law.hpp"
#include "custom_laws/small_strain_laws/small_strain_plane_stress_2D_law.hpp"
#include "custom_laws/small_strain_laws/small_strain_axisymmetric_2D_law.hpp"

// Large strain laws
#include "custom_laws/large_strain_laws/large_strain_plane_strain_2D_law.hpp"
#include "custom_laws/large_strain_laws/large_strain_axisymmetric_2D_law.hpp"

// Strain rate laws
#include "custom_laws/strain_rate_laws/strain_rate_plane_strain_2D_law.hpp"
#include "custom_laws/strain_rate_laws/newtonian_plane_strain_2D_law.hpp"

// Constitutive models

// Elasticity models
#include "custom_models/elasticity_models/linear_elastic_model.hpp"
#include "custom_models/elasticity_models/saint_venant_kirchhoff_model.hpp"
#include "custom_models/elasticity_models/neo_hookean_model.hpp"
#include "custom_models/elasticity_models/neo_hookean_lnJ_squared_model.hpp"
#include "custom_models/elasticity_models/neo_hookean_J_1_squared_model.hpp"
#include "custom_models/elasticity_models/isochoric_neo_hookean_model.hpp"
#include "custom_models/elasticity_models/isochoric_neo_hookean_lnJ_squared_model.hpp"
#include "custom_models/elasticity_models/incompressible_neo_hookean_model.hpp"
#include "custom_models/elasticity_models/borja_model.hpp"
#include "custom_models/elasticity_models/ogden_model.hpp"
#include "custom_models/elasticity_models/isochoric_ogden_model.hpp"
#include "custom_models/elasticity_models/incompressible_hypo_elastic_model.hpp"

// Plasticity models
#include "custom_models/plasticity_models/von_mises_linear_elastic_plasticity_model.hpp"
#include "custom_models/plasticity_models/von_mises_neo_hookean_plasticity_model.hpp"
#include "custom_models/plasticity_models/simo_J2_plasticity_model.hpp"
#include "custom_models/plasticity_models/simo_J2_thermo_plasticity_model.hpp"
#include "custom_models/plasticity_models/johnson_cook_J2_thermo_plasticity_model.hpp"
#include "custom_models/plasticity_models/baker_johnson_cook_J2_thermo_plasticity_model.hpp"
#include "custom_models/plasticity_models/cam_clay_model.hpp"
#include "custom_models/plasticity_models/simo_ju_exponential_damage_model.hpp"
#include "custom_models/plasticity_models/simo_ju_modified_exponential_damage_model.hpp"


namespace Kratos
{
namespace Python
{

using namespace pybind11;

typedef ConstitutiveLaw                                              ConstitutiveLawBaseType;
typedef ConstitutiveLaw::Pointer                                      ConstitutiveLawPointer;
typedef std::vector<ConstitutiveLaw::Pointer>                             MaterialsContainer;

typedef ConstitutiveModel                                          ConstitutiveModelBaseType;
typedef ConstitutiveModel::Pointer                                  ConstitutiveModelPointer;



void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{

  //outfitted python laws
  // class_< PythonOutfittedConstitutiveLaw, typename PythonOutfittedConstitutiveLaw::Pointer, ConstitutiveLawBaseType >
  //  	(m, "PythonOutfittedConstitutiveLaw")
  //  	.def( init<>() )
  //  	.def(init<PyObject* >())
  //  	;

  //general constitutive laws

  //small strain laws
  class_< SmallStrain3DLaw, typename SmallStrain3DLaw::Pointer, ConstitutiveLawBaseType >(m, "SmallStrain3DLaw")
      .def( init<>() )
      .def( init<ConstitutiveModelPointer>() )
      ;

  class_< SmallStrainOrthotropic3DLaw, typename SmallStrainOrthotropic3DLaw::Pointer, ConstitutiveLawBaseType >
      (m, "SmallStrainOrthotropic3DLaw")
      .def( init<>() )
      .def( init<ConstitutiveModelPointer>() )
      ;

  class_< SmallStrainPlaneStrain2DLaw, typename SmallStrainPlaneStrain2DLaw::Pointer, ConstitutiveLawBaseType >
      (m, "SmallStrainPlaneStrain2DLaw")
      .def( init<>() )
      .def( init<ConstitutiveModelPointer>() )
      ;

  class_< SmallStrainPlaneStress2DLaw, typename SmallStrainPlaneStress2DLaw::Pointer, ConstitutiveLawBaseType >
      (m, "SmallStrainPlaneStress2DLaw")
      .def( init<>() )
      .def( init<ConstitutiveModelPointer>() )
      ;

  class_< SmallStrainAxisymmetric2DLaw, typename SmallStrainAxisymmetric2DLaw::Pointer, ConstitutiveLawBaseType >
      (m, "SmallStrainAxisymmetric2DLaw")
      .def( init<>() )
      .def( init<ConstitutiveModelPointer>() )
      ;



  //large strain laws
  class_< LargeStrain3DLaw, typename LargeStrain3DLaw::Pointer, ConstitutiveLawBaseType >
      (m, "LargeStrain3DLaw")
      .def( init<ConstitutiveModelPointer>() )
      ;

  class_< LargeStrainPlaneStrain2DLaw, typename LargeStrainPlaneStrain2DLaw::Pointer, ConstitutiveLawBaseType >
      (m, "LargeStrainPlaneStrain2DLaw")
      .def( init<ConstitutiveModelPointer>() )
      ;

  class_< LargeStrainAxisymmetric2DLaw, typename LargeStrainAxisymmetric2DLaw::Pointer, ConstitutiveLawBaseType >
      (m, "LargeStrainAxisymmetric2DLaw")
      .def( init<ConstitutiveModelPointer>() )
      ;


  //strain rate laws
  class_< StrainRate3DLaw, typename StrainRate3DLaw::Pointer, ConstitutiveLawBaseType >
      (m, "StrainRate3DLaw")
      .def( init<ConstitutiveModelPointer>() )
      ;

  class_< StrainRatePlaneStrain2DLaw, typename StrainRatePlaneStrain2DLaw::Pointer, ConstitutiveLawBaseType >
      (m, "StrainRatePlaneStrain2DLaw")
      .def( init<ConstitutiveModelPointer>() )
      ;

  class_< Newtonian3DLaw, typename Newtonian3DLaw::Pointer, ConstitutiveLawBaseType >(m, "Newtonian3DLaw")
      .def( init<>() )
      ;

  class_< NewtonianPlaneStrain2DLaw, typename NewtonianPlaneStrain2DLaw::Pointer, ConstitutiveLawBaseType >(m, "NewtonianPlaneStrain2DLaw")
      .def( init<>() )
      ;

  //general constitutive models
  class_< ConstitutiveModelBaseType, ConstitutiveModelPointer>(m, "ConstitutiveModelModel")
      .def( init<>() )
      ;


  //elasticity models
  class_< LinearElasticModel, typename LinearElasticModel::Pointer, ConstitutiveModelBaseType >
      (m, "LinearElasticModel")
      .def( init<>() )
      ;

  class_< SaintVenantKirchhoffModel, typename SaintVenantKirchhoffModel::Pointer, ConstitutiveModelBaseType >
      (m, "SaintVenantKirchhoffModel")
      .def( init<>() )
      ;

  class_< NeoHookeanModel, typename NeoHookeanModel::Pointer, ConstitutiveModelBaseType >
      (m, "NeoHookeanModel")
      .def( init<>() )
      ;

  class_< NeoHookeanLnJSquaredModel, typename NeoHookeanLnJSquaredModel::Pointer, ConstitutiveModelBaseType >
      (m, "NeoHookeanLnJSquaredModel")
      .def( init<>() )
      ;

  class_< NeoHookeanJ_1SquaredModel, typename NeoHookeanJ_1SquaredModel::Pointer, ConstitutiveModelBaseType >
      (m, "NeoHookeanJ_1SquaredModel")
      .def( init<>() )
      ;

  class_< IsochoricNeoHookeanModel, typename IsochoricNeoHookeanModel::Pointer, ConstitutiveModelBaseType >
      (m, "IsochoricNeoHookeanModel")
      .def( init<>() )
      ;

  class_< IsochoricNeoHookeanLnJSquaredModel, typename IsochoricNeoHookeanLnJSquaredModel::Pointer, ConstitutiveModelBaseType >
      (m, "IsochoricNeoHookeanLnJSquaredModel")
      .def( init<>() )
      ;

  class_< IncompressibleNeoHookeanModel, typename IncompressibleNeoHookeanModel::Pointer, ConstitutiveModelBaseType >
      (m, "IncompressibleNeoHookeanModel")
      .def( init<>() )
      ;

  class_< BorjaModel, typename BorjaModel::Pointer, ConstitutiveModelBaseType >
      (m, "BorjaModel")
      .def( init<>() )
      ;

  class_< OgdenModel, typename OgdenModel::Pointer, ConstitutiveModelBaseType >
      (m, "OgdenModel")
      .def( init<>() )
      ;

  class_< IsochoricOgdenModel, typename IsochoricOgdenModel::Pointer, ConstitutiveModelBaseType >
      (m, "IsochoricOgdenModel")
      .def( init<>() )
      ;

  class_< HypoElasticModel, typename HypoElasticModel::Pointer, ConstitutiveModelBaseType >
      (m, "HypoElasticModel")
      .def( init<>() )
      ;

  class_< IsochoricHypoElasticModel, typename IsochoricHypoElasticModel::Pointer, ConstitutiveModelBaseType >
      (m, "IsochoricHypoElasticModel")
      .def( init<>() )
      ;

  class_< IncompressibleHypoElasticModel, typename IncompressibleHypoElasticModel::Pointer, ConstitutiveModelBaseType >
      (m, "IncompressibleHypoElasticModel")
      .def( init<>() )
      ;

  //plasticity models
  class_< VonMisesLinearElasticPlasticityModel, typename VonMisesLinearElasticPlasticityModel::Pointer, ConstitutiveModelBaseType >
      (m, "VonMisesLinearElasticPlasticityModel")
      .def( init<>() )
      ;

  class_< VonMisesNeoHookeanPlasticityModel, typename VonMisesNeoHookeanPlasticityModel::Pointer, ConstitutiveModelBaseType >
      (m, "VonMisesNeoHookeanPlasticityModel")
      .def( init<>() )
      ;

  class_< SimoJ2PlasticityModel, typename SimoJ2PlasticityModel::Pointer, ConstitutiveModelBaseType >
      (m, "SimoJ2PlasticityModel")
      .def( init<>() )
      ;

  class_< SimoJ2ThermoPlasticityModel, typename SimoJ2ThermoPlasticityModel::Pointer, ConstitutiveModelBaseType >
      (m, "SimoJ2ThermoPlasticityModel")
      .def( init<>() )
      ;

  class_< JohnsonCookJ2ThermoPlasticityModel, typename JohnsonCookJ2ThermoPlasticityModel::Pointer, ConstitutiveModelBaseType >
      (m, "JohnsonCookJ2ThermoPlasticityModel")
      .def( init<>() )
      ;

  class_< BakerJohnsonCookJ2ThermoPlasticityModel, typename BakerJohnsonCookJ2ThermoPlasticityModel::Pointer, ConstitutiveModelBaseType >
      (m, "BakerJohnsonCookJ2ThermoPlasticityModel")
      .def( init<>() )
      ;

  class_< CamClayModel, typename CamClayModel::Pointer, ConstitutiveModelBaseType >
      (m, "CamClayModel")
      .def( init<>() )
      ;

  class_< SimoJuExponentialDamageModel, typename SimoJuExponentialDamageModel::Pointer, ConstitutiveModelBaseType >
      (m, "SimoJuExponentialDamageModel")
      .def( init<>() )
      ;

  class_< SimoJuModifiedExponentialDamageModel, typename SimoJuModifiedExponentialDamageModel::Pointer, ConstitutiveModelBaseType >
      (m, "SimoJuModifiedExponentialDamageModel")
      .def( init<>() )
      ;
}

}  // namespace Python.
}  // namespace Kratos.
