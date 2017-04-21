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

//elasticity laws

//isotropic linear elastic laws
#include "custom_laws/linear_elastic_laws/linear_elastic_plane_strain_2D_law.hpp"
#include "custom_laws/linear_elastic_laws/linear_elastic_plane_stress_2D_law.hpp"
#include "custom_laws/linear_elastic_laws/linear_elastic_axisymmetric_2D_law.hpp"

//orthotropic linear elastic laws
#include "custom_laws/linear_elastic_laws/linear_elastic_orthotropic_3D_law.hpp"

//isotropic hyperelastic laws
#include "custom_laws/hyperelastic_laws/hyperelastic_3D_law.hpp"
#include "custom_laws/hyperelastic_laws/hyperelastic_plane_strain_2D_law.hpp"
#include "custom_laws/hyperelastic_laws/hyperelastic_axisymmetric_2D_law.hpp"
#include "custom_laws/hyperelastic_laws/hyperelastic_U_P_plane_strain_2D_law.hpp"
#include "custom_laws/hyperelastic_laws/hyperelastic_U_P_axisymmetric_2D_law.hpp"

//specialized isotropic hyperelastic laws
#include "custom_laws/hyperelastic_laws/neo_hookean_3D_law.hpp"
#include "custom_laws/hyperelastic_laws/compressible_neo_hookean_3D_law.hpp"
#include "custom_laws/hyperelastic_laws/isochoric_neo_hookean_3D_law.hpp"
#include "custom_laws/hyperelastic_laws/saint_venant_kirchhoff_3D_law.hpp"

//plasticity laws

//isotropic linear elastic plasticity laws
//#include "custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_3D_law.hpp"
//#include "custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_plane_strain_2D_law.hpp"
//#include "custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_plane_stress_2D_law.hpp"
//#include "custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_axysimmetric_2D_law.hpp"

//isotropic hyperelastic plasticity laws
#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_3D_law.hpp"
#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_plane_strain_2D_law.hpp"
#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_axisymmetric_2D_law.hpp"

#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_U_P_3D_law.hpp"
#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_U_P_plane_strain_2D_law.hpp"
#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_U_P_axisymmetric_2D_law.hpp"

//specialized isotropic hyperelastic plastic laws
#include "custom_laws/hyperelastic_plastic_laws/von_mises_hyperelastic_plastic_3D_law.hpp"


//isotropic linear elastic damage laws
//#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_simo_ju_3D_law.hpp"
//#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_simo_ju_plane_strain_2D_law.hpp"
//#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_simo_ju_plane_stress_2D_law.hpp"

//#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_modified_mises_3D_law.hpp"
//#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_modified_mises_plane_strain_2D_law.hpp"
//#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_modified_mises_plane_stress_2D_law.hpp"

//elasticity models
#include "custom_models/elasticity_models/linear_elastic_model.hpp"

//hyperelastic models
#include "custom_models/elasticity_models/hyperelastic_models/saint_venant_kirchhoff_model.hpp"
#include "custom_models/elasticity_models/hyperelastic_models/neo_hookean_model.hpp"
#include "custom_models/elasticity_models/hyperelastic_models/compressible_neo_hookean_model.hpp"
#include "custom_models/elasticity_models/hyperelastic_models/isochoric_neo_hookean_model.hpp"
#include "custom_models/elasticity_models/hyperelastic_models/incompressible_neo_hookean_model.hpp"

//plasticity models
#include "custom_models/plasticity_models/von_mises_plasticity_model.hpp"
#include "custom_models/plasticity_models/von_mises_neo_hookean_plasticity_model.hpp"

//yield criteria
#include "custom_models/plasticity_models/yield_criteria/mises_huber_yield_criterion.hpp"
#include "custom_models/plasticity_models/yield_criteria/simo_ju_yield_criterion.hpp"
#include "custom_models/plasticity_models/yield_criteria/modified_mises_yield_criterion.hpp"

//hardening laws
#include "custom_models/plasticity_models/hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"
#include "custom_models/plasticity_models/hardening_laws/exponential_damage_hardening_law.hpp"
#include "custom_models/plasticity_models/hardening_laws/modified_exponential_damage_hardening_law.hpp"

namespace Kratos
{
  namespace Python
  {

    using namespace boost::python;

    typedef Properties::Pointer                                                PropertiesPointer;

    typedef ConstitutiveLaw                                              ConstitutiveLawBaseType;
    typedef ConstitutiveLaw::Pointer                                      ConstitutiveLawPointer;
    typedef std::vector<ConstitutiveLaw::Pointer>                             MaterialsContainer;
 
    typedef ElasticityModel                                                 ElasticModelBaseType;
    typedef ElasticityModel::Pointer                                         ElasticModelPointer;
    
    typedef HyperElasticModel                                          HyperElasticModelBaseType;
    typedef HyperElasticModel::Pointer                                  HyperElasticModelPointer;

    typedef HardeningLaw                                                    HardeningLawBaseType;
    typedef HardeningLaw::Pointer                                            HardeningLawPointer;

    typedef YieldCriterion<HardeningLawBaseType>                                   YieldCriterionBaseType;
    typedef typename YieldCriterionBaseType::Pointer                                YieldCriterionPointer;
    
    typedef PlasticityModel<ElasticModelBaseType,YieldCriterionBaseType>          PlasticityModelBaseType;  
    typedef typename PlasticityModelBaseType::Pointer                              PlasticityModelPointer;

    typedef PlasticityModel<HyperElasticModelBaseType,YieldCriterionBaseType>   HyperElasticPlasticModelBaseType;  
    typedef typename HyperElasticPlasticModelBaseType::Pointer                   HyperElasticPlasticModelPointer;

    
    // template<class THyperElasticModel>
    // void export_von_mises_plasticity_model(std::string name) {
      
    //   typedef typename THyperElasticModel::Pointer               THyperElasticModelPointer;
    //   typedef VonMisesPlasticityModel<THyperElasticModel>      VonMisesPlasticityModelType;
    //   typedef typename VonMisesPlasticityModelType::Pointer VonMisesPlasticityModelPointer;
			      
    //   class_< VonMisesPlasticityModel<THyperElasticModel>, VonMisesPlasticityModelPointer, boost::noncopyable >
    // 	(  name.c_str(),
    // 	   init<THyperElasticModelPointer>() )
    // 	;
    // }
    
    
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
      
      //elasticity laws
      
      //isotropic linear elastic plasticity laws
      class_< LinearElastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
       	( "LinearElastic3DLaw",
       	  init<>() )
       	;

      class_< LinearElasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "LinearElasticPlaneStrain2DLaw",
      	  init<>() )
      	;

      class_< LinearElasticPlaneStress2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "LinearElasticPlaneStress2DLaw",
      	  init<>() )
      	;

      class_< LinearElasticAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "LinearElasticAxisymmetric2DLaw",
      	  init<>() )
      	;

      //orthotropic linear elastic laws
      class_< LinearElasticOrthotropic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
       	( "LinearElasticOrthotropic3DLaw",
       	  init<>() )
      	;

      
      //isotropic hyperelastic laws
      class_< HyperElastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "HyperElastic3DLaw",
      	  init<>() )
      	.def( init<HyperElasticModelPointer>() )
      	;
    
      class_< HyperElasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "HyperElasticPlaneStrain2DLaw",
      	  init<>() )
      	.def( init<HyperElasticModelPointer>() )
      	;

      class_< HyperElasticAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "HyperElasticAxisymmetric2DLaw",
      	  init<>() )
      	.def( init<HyperElasticModelPointer>() )
      	;

      class_< HyperElasticUP3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "HyperElasticUP3DLaw",
      	  init<>() )
      	.def( init<HyperElasticModelPointer>() )
      	;
    
      class_< HyperElasticUPPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "HyperElasticUPPlaneStrain2DLaw",
      	  init<>() )
      	.def( init<HyperElasticModelPointer>() )
      	;

      class_< HyperElasticUPAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      	( "HyperElasticUPAxisymmetric2DLaw",
      	  init<>() )
      	.def( init<HyperElasticModelPointer>() )
      	;    

      //specialized isotropic hyperelastic laws
      class_< NeoHookean3DLaw, bases< HyperElastic3DLaw >, boost::noncopyable >
      	( "NeoHookean3DLaw",
      	  init<>() )
      	;

      class_< CompressibleNeoHookean3DLaw, bases< HyperElastic3DLaw >, boost::noncopyable >
      	( "CompressibleNeoHookean3DLaw",
      	  init<>() )
      	;
      
      class_< IsochoricNeoHookean3DLaw, bases< HyperElastic3DLaw >, boost::noncopyable >
      	( "IsochoricNeoHookean3DLaw",
      	  init<>() )
      	;      

      class_< SaintVenantKirchhoff3DLaw, bases< HyperElastic3DLaw >, boost::noncopyable >
      	( "SaintVenantKirchhoff3DLaw",
      	  init<>() )
      	;
      
      //plasticity laws
     

      //isotropic linear elastic plasticity laws
      // class_< LinearElasticPlastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // ( "LinearElasticPlastic3DLaw",
      //   init<>() )
      //   .def( init<HardeningLawPointer>() )
      // ;

      // class_< LinearElasticPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // ( "LinearElasticPlasticPlaneStrain2DLaw",
      //   init<>() )
      //   .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      // ;

      // class_< LinearElasticPlasticPlaneStress2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // ( "LinearElasticPlasticPlaneStress2DLaw",
      //   init<>() )
      //   .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      // ;
      
      // class_< LinearElasticPlasticAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // ( "LinearElasticPlasticAxisymmetric2DLaw",
      //   init<>() )
      //   .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      // ;

      //isotropic hyperelastic plasticity laws
      class_<HyperElasticPlastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      ( "HyperElasticPlastic3DLaw",
        init<>() )
        .def( init<HyperElasticPlasticModelPointer>() )
      ;

      class_<HyperElasticPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      ( "HyperElasticPlasticPlaneStrain2DLaw",
        init<>() )
        .def( init<HyperElasticPlasticModelPointer>() )
      ;

      class_<HyperElasticPlasticAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      ( "HyperElasticPlasticAxisym2DLaw",
        init<>() )
        .def( init<HyperElasticPlasticModelPointer>() )
      ;

      class_<HyperElasticPlasticUP3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      ( "HyperElasticPlasticUP3DLaw",
        init<>() )
        .def( init<HyperElasticPlasticModelPointer>() )
      ;

      class_<HyperElasticPlasticUPPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      ( "HyperElasticPlasticUPPlaneStrain2DLaw",
        init<>() )
        .def( init<HyperElasticPlasticModelPointer>() )
      ;

      class_<HyperElasticPlasticUPAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      ( "HyperElasticPlasticAxisym2DLaw",
        init<>() )
        .def( init<HyperElasticPlasticModelPointer>() )
      ;      

      
      //specialized isotropic hyperelastic plastic laws
      class_<VonMisesHyperElasticPlastic3DLaw, bases< HyperElasticPlastic3DLaw >, boost::noncopyable >
      ( "VonMisesHyperElasticPlastic3DLaw",
	 init<>())
      ;
      
      //isotropic linear elastic damage specilization laws
      // class_< IsotropicDamageSimoJu3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "IsotropicDamageSimoJu3DLaw",
      // 	  init<>() )
      // 	;

      // class_< IsotropicDamageSimoJuPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "IsotropicDamageSimoJuPlaneStrain2DLaw",
      // 	  init<>() )
      // 	;

      // class_< IsotropicDamageSimoJuPlaneStress2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "IsotropicDamageSimoJuPlaneStress2DLaw",
      // 	  init<>() )
      // 	;

      // class_< IsotropicDamageModifiedMises3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "IsotropicDamageModifiedMises3DLaw",
      // 	  init<>() )
      // 	;

      // class_< IsotropicDamageModifiedMisesPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "IsotropicDamageModifiedMisesPlaneStrain2DLaw",
      // 	  init<>() )
      // 	;

      // class_< IsotropicDamageModifiedMisesPlaneStress2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "IsotropicDamageModifiedMisesPlaneStress2DLaw",
      // 	  init<>() )
      // 	;


      //hyperelastic model
      class_< HyperElasticModel, HyperElasticModelPointer, boost::noncopyable >
       	( "HyperElasticModel",
       	  init<>() )
       	;

      class_< SaintVenantKirchhoffModel, bases< HyperElasticModelBaseType >, boost::noncopyable >
      	( "SaintVenantKirchhoffModel",
	  init<>() )
       	;
      
      class_< NeoHookeanModel, bases< HyperElasticModelBaseType >, boost::noncopyable >
      	( "NeoHookeanModel",
	  init<>() )
       	;
      
      class_< CompressibleNeoHookeanModel, bases< HyperElasticModelBaseType >, boost::noncopyable >
      	( "CompressibleNeoHookeanModel",
	  init<>() )
       	;

      class_< IncompressibleNeoHookeanModel, bases< HyperElasticModelBaseType >, boost::noncopyable >
      	( "IncompressibleNeoHookeanModel",
	  init<>() )
       	;

      
      //hardening laws
      class_< HardeningLawBaseType, HardeningLawPointer, boost::noncopyable >
       	( "HardeningLaw",
       	  init<>() )
       	;

      //yield criteria
      class_< YieldCriterionBaseType, YieldCriterionPointer, boost::noncopyable >
       	( "YieldCriterion",
       	  init<>() )
	.def( init<HardeningLawPointer>() )
       	;

      //plasticity models
      class_< PlasticityModelBaseType, PlasticityModelPointer, boost::noncopyable >
       	( "PlasticityModel",
       	  init<>() )
      	.def( init<ElasticModelPointer,YieldCriterionPointer>() )
       	;

      class_< HyperElasticPlasticModelBaseType, HyperElasticPlasticModelPointer, boost::noncopyable >
       	( "HyperElasticPlasticModel",
       	  init<>() )
      	.def( init<HyperElasticModelPointer,YieldCriterionPointer>() )
       	;
      
      //export_von_mises_plasticity_model<NeoHookeanModel>("VonMisesNeoHookeanModel");

      class_< VonMisesNeoHookeanPlasticityModel, bases< HyperElasticPlasticModelBaseType >, boost::noncopyable >
       	( "VonMisesNeoHookeanPlasticityModel",
       	  init<>() )
       	;
      
      
    }

  }  // namespace Python.
}  // namespace Kratos.
