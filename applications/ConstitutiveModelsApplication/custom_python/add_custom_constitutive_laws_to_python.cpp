//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "includes/properties.h"

#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"


//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//outfitted python laws
#include "custom_python/python_outfitted_constitutive_law.hpp"

//general constitutive laws

//small strain laws
#include "custom_laws/small_strain_laws/small_strain_orthotropic_3D_law.hpp"
#include "custom_laws/small_strain_laws/small_strain_plane_strain_2D_law.hpp"
#include "custom_laws/small_strain_laws/small_strain_plane_stress_2D_law.hpp"
#include "custom_laws/small_strain_laws/small_strain_axisymmetric_2D_law.hpp"

//large strain laws
#include "custom_laws/large_strain_laws/large_strain_plane_strain_2D_law.hpp"
#include "custom_laws/large_strain_laws/large_strain_axisymmetric_2D_law.hpp"

//general constitutive models

//elasticity models
#include "custom_models/elasticity_models/linear_elastic_model.hpp"
#include "custom_models/elasticity_models/saint_venant_kirchhoff_model.hpp"
#include "custom_models/elasticity_models/neo_hookean_model.hpp"
#include "custom_models/elasticity_models/neo_hookean_lnJ_squared_model.hpp"
#include "custom_models/elasticity_models/neo_hookean_J_1_squared_model.hpp"
#include "custom_models/elasticity_models/isochoric_neo_hookean_model.hpp"
#include "custom_models/elasticity_models/isochoric_neo_hookean_lnJ_squared_model.hpp"
#include "custom_models/elasticity_models/incompressible_neo_hookean_model.hpp"
#include "custom_models/elasticity_models/borja_model.hpp"

//plasticity models
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

    using namespace boost::python;

    typedef Properties::Pointer                                                PropertiesPointer;

    typedef ConstitutiveLaw                                              ConstitutiveLawBaseType;
    typedef ConstitutiveLaw::Pointer                                      ConstitutiveLawPointer;
    typedef std::vector<ConstitutiveLaw::Pointer>                             MaterialsContainer;
 
    typedef ConstitutiveModel                                          ConstitutiveModelBaseType;
    typedef ConstitutiveModel::Pointer                                  ConstitutiveModelPointer;
    
    
    void Push_Back_Constitutive_Laws( MaterialsContainer& ThisMaterialsContainer,
				      ConstitutiveLawPointer ThisConstitutiveLaw )
    {
      ThisMaterialsContainer.push_back( ThisConstitutiveLaw );
    }

    void  AddCustomConstitutiveLawsToPython()
    {
      class_< MaterialsContainer >( "MaterialsContainer", init<>() )
	.def( "PushBack", Push_Back_Constitutive_Laws )
	;

      //outfitted python laws
      class_< PythonOutfittedConstitutiveLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
       	( "PythonOutfittedConstitutiveLaw",
       	  init<>() )
       	.def(init<PyObject* >())
       	;
      
      //general constitutive laws
      
      //small strain laws
      class_< SmallStrain3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
       	( "SmallStrain3DLaw",
       	  init<>() )
      	.def( init<ConstitutiveModelPointer>() )
       	;

      class_< SmallStrainOrthotropic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
       	( "SmallStrainOrthotropic3DLaw",
       	  init<>() )
      	.def( init<ConstitutiveModelPointer>() )
      	;
      
      class_< SmallStrainPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "SmallStrainPlaneStrain2DLaw",
      	  init<>() )
      	.def( init<ConstitutiveModelPointer>() )
      	;

      class_< SmallStrainPlaneStress2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "SmallStrainPlaneStress2DLaw",
      	  init<>() )
      	.def( init<ConstitutiveModelPointer>() )
      	;

      class_< SmallStrainAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "SmallStrainAxisymmetric2DLaw",
      	  init<>() )
      	.def( init<ConstitutiveModelPointer>() )
      	;


      
      //large strain laws
      class_< LargeStrain3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "LargeStrain3DLaw",
      	  init<ConstitutiveModelPointer>() )
      	;
    
      class_< LargeStrainPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "LargeStrainPlaneStrain2DLaw",
      	  init<ConstitutiveModelPointer>() )
      	;

      class_< LargeStrainAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "LargeStrainAxisymmetric2DLaw",
      	  init<ConstitutiveModelPointer>() )
      	;

      //general constitutive models
      class_< ConstitutiveModelBaseType, ConstitutiveModelPointer, boost::noncopyable >
       	( "ConstitutiveModelModel",
       	  init<>() )
       	;

      
      //elasticity models      
      class_< LinearElasticModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
      	( "LinearElasticModel",
      	  init<>() )
       	;

      class_< SaintVenantKirchhoffModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
      	( "SaintVenantKirchhoffModel",
      	  init<>() )
       	;
      
      class_< NeoHookeanModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
      	( "NeoHookeanModel",
      	  init<>() )
       	;
      
      class_< NeoHookeanLnJSquaredModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
      	( "NeoHookeanLnJSquaredModel",
      	  init<>() )
       	;

      class_< NeoHookeanJ_1SquaredModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
      	( "NeoHookeanJ_1SquaredModel",
      	  init<>() )
       	;

      class_< IsochoricNeoHookeanModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
      	( "IsochoricNeoHookeanModel",
      	  init<>() )
       	;

      class_< IsochoricNeoHookeanLnJSquaredModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
      	( "IsochoricNeoHookeanLnJSquaredModel",
      	  init<>() )
       	;
      
      class_< IncompressibleNeoHookeanModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
      	( "IncompressibleNeoHookeanModel",
      	  init<>() )
       	;
      class_< BorjaModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
      	( "BorjaModel",
      	  init<>() )
       	;
      
      //plasticity models
      class_< VonMisesLinearElasticPlasticityModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
       	( "VonMisesLinearElasticPlasticityModel",
       	  init<>() )
       	;
      
      class_< VonMisesNeoHookeanPlasticityModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
       	( "VonMisesNeoHookeanPlasticityModel",
       	  init<>() )
       	;
      
      class_< SimoJ2PlasticityModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
       	( "SimoJ2PlasticityModel",
       	  init<>() )
       	;

      class_< SimoJ2ThermoPlasticityModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
       	( "SimoJ2ThermoPlasticityModel",
       	  init<>() )
       	;

      class_< JohnsonCookJ2ThermoPlasticityModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
       	( "JohnsonCookJ2ThermoPlasticityModel",
       	  init<>() )
       	;

      class_< BakerJohnsonCookJ2ThermoPlasticityModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
       	( "BakerJohnsonCookJ2ThermoPlasticityModel",
       	  init<>() )
       	;
      
      class_< CamClayModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
       	( "CamClayModel",
       	  init<>() )
       	;
      
      class_< SimoJuExponentialDamageModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
       	( "SimoJuExponentialDamageModel",
       	  init<>() )
       	;

      class_< SimoJuModifiedExponentialDamageModel, bases< ConstitutiveModelBaseType >, boost::noncopyable >
       	( "SimoJuModifiedExponentialDamageModel",
       	  init<>() )
       	;
    }

  }  // namespace Python.
}  // namespace Kratos.
