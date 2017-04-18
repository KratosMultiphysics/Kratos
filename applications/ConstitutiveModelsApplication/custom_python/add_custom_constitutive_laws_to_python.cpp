//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
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
//#include "custom_python/python_outfitted_constitutive_law.hpp"

//general constitutive laws
#include "custom_laws/elastic_3D_law.hpp"

//isotropic linear elastic laws
#include "custom_laws/linear_elastic_laws/linear_elastic_3D_law.hpp"
//#include "custom_laws/linear_elastic_laws/linear_elastic_plane_strain_2D_law.hpp"
//#include "custom_laws/linear_elastic_laws/linear_elastic_plane_stress_2D_law.hpp"
//#include "custom_laws/linear_elastic_laws/linear_elastic_axisymmetric_2D_law.hpp"

//orthotropic linear elastic laws
//#include "custom_laws/linear_elastic_laws/linear_elastic_orthotropic_3D_law.hpp"

//isotropic hyperelastic laws
#include "custom_laws/hyperelastic_laws/hyperelastic_3D_law.hpp"
//#include "custom_laws/hyperelastic_laws/hyperelastic_plane_strain_2D_law.hpp"
//#include "custom_laws/hyperelastic_laws/hyperelastic_axisymmetric_2D_law.hpp"

//#include "custom_laws/hyperelastic_laws/hyperelastic_U_P_3D_law.hpp"
//#include "custom_laws/hyperelastic_laws/hyperelastic_U_P_plane_strain_2D_law.hpp"
//#include "custom_laws/hyperelastic_laws/hyperelastic_U_P_axisymmetric_2D_law.hpp"

//plasticity models

//isotropic linear elastic plasticity laws
//#include "custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_3D_law.hpp"
//#include "custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_plane_strain_2D_law.hpp"
//#include "custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_plane_stress_2D_law.hpp"
//#include "custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_axysimmetric_2D_law.hpp"

//isotropic hyperelastic plasticity laws
//#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_3D_law.hpp"
//#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_plane_strain_2D_law.hpp"
//#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_axisymmetric_2D_law.hpp"

//#include "custom_laws/hyperelastic_plastic_laws/incompressible_hyperelastic_plastic_3D_law.hpp"
//#include "custom_laws/hyperelastic_plastic_laws/incompressible_hyperelastic_plastic_plane_strain_2D_law.hpp"
//#include "custom_laws/hyperelastic_plastic_laws/incompressible_hyperelastic_plastic_axisymmetric_2D_law.hpp"

//#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_J2_3D_law.hpp"
//#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_J2_plane_strain_2D_law.hpp"
//#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_J2_axisymmetric_2D_law.hpp"

//#include "custom_laws/hyperelastic_plastic_laws/incompressible_hyperelastic_plastic_J2_3D_law.hpp"
//#include "custom_laws/hyperelastic_plastic_laws/incompressible_hyperelastic_plastic_J2_plane_strain_2D_law.hpp"
//#include "custom_laws/hyperelastic_plastic_laws/incompressible_hyperelastic_plastic_J2_axisymmetric_2D_law.hpp"

//isotropic linear elastic damage laws
//#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_simo_ju_3D_law.hpp"
//#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_simo_ju_plane_strain_2D_law.hpp"
//#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_simo_ju_plane_stress_2D_law.hpp"

//#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_modified_mises_3D_law.hpp"
//#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_modified_mises_plane_strain_2D_law.hpp"
//#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_modified_mises_plane_stress_2D_law.hpp"


//elasticity models
//#include "custom_models/elasticity_models/elasticity_model.hpp"

//hyperelastic models
#include "custom_models/elasticity_models/hyperelastic_models/neo_hookean_model.hpp"

//plastic models
#include "custom_models/plasticity_models/non_linear_associative_plastic_model.hpp"

//hardening laws
#include "custom_models/plasticity_models/hardening_laws/hardening_law.hpp"
//#include "custom_models/plasticity_models/hardening_laws/non_linear_isotropic_kinematic_hardening_law.hpp"
//#include "custom_models/plasticity_models/hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"
//#include "custom_models/plasticity_models/hardening_laws/exponential_damage_hardening_law.hpp"
//#include "custom_models/plasticity_models/hardening_laws/modified_exponential_damage_hardening_law.hpp"

//yield criteria
#include "custom_models/plasticity_models/yield_criteria/yield_criterion.hpp"
//#include "custom_models/plasticity_models/yield_criteria/mises_huber_yield_criterion.hpp"
//#include "custom_models/plasticity_models/yield_criteria/simo_ju_yield_criterion.hpp"
//#include "custom_models/plasticity_models/yield_criteria/modified_mises_yield_criterion.hpp"


namespace Kratos
{
  namespace Python
  {

    using namespace boost::python;

    // typedef YieldCriterion::Pointer            YieldCriterionPointer;
    // typedef HardeningLaw::Pointer                HardeningLawPointer;
    typedef Properties::Pointer                    PropertiesPointer;

    typedef ConstitutiveLaw                  ConstitutiveLawBaseType;
    typedef ConstitutiveLaw::Pointer          ConstitutiveLawPointer;
    typedef std::vector<ConstitutiveLaw::Pointer> MaterialsContainer;

    // typedef HyperElasticModel              HyperElasticModelBaseType;
    // typedef HyperElasticModel::Pointer      HyperElasticModelPointer;


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

      // //outfitted python laws
      // class_< PythonOutfittedConstitutiveLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      //  	( "PythonOutfittedConstitutiveLaw",
      //  	  init<>() )
      //  	.def(init<PyObject* >())
      //  	;
      
      // //elasticity models
      
      // //isotropic linear elastic plasticity laws
      // class_< LinearElastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      //  	( "LinearElastic3DLaw",
      //  	  init<>() )
      //  	;

      // class_< LinearElasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "LinearElasticPlaneStrain2DLaw",
      // 	  init<>() )
      // 	;

      // class_< LinearElasticPlaneStress2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "LinearElasticPlaneStress2DLaw",
      // 	  init<>() )
      // 	;

      // class_< LinearElasticAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "LinearElasticAxisymmetric2DLaw",
      // 	  init<>() )
      // 	;

      // //orthotropic linear elastic laws
      // // class_< LinearElasticOrthotropic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // //  	( "LinearElasticOrthotropic3DLaw",
      // //  	  init<>() )
      // // 	;

      // //hyperelastic model
      // class_< HyperElasticModel, HyperElasticModelPointer, boost::noncopyable >
      // 	( "HyperElasticModel",
      // 	  init<>() )
      // 	;
    
      // class_< NeoHookeanModel, bases< HyperElasticModelBaseType >, boost::noncopyable >
      // 	( "NeoHookeanModel",
      //  	  init<>() )
      // 	;

      
      // //isotropic hyperelastic laws
      // class_< HyperElastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "HyperElastic3DLaw",
      // 	  init<>() )
      // 	.def( init<HyperElasticModelPointer>() )
      // 	;
    
      // class_< HyperElasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "HyperElasticPlaneStrain2DLaw",
      // 	  init<>() )
      // 	.def( init<HyperElasticModelPointer>() )
      // 	;

      // class_< HyperElasticAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "HyperElasticAxisymmetric2DLaw",
      // 	  init<>() )
      // 	.def( init<HyperElasticModelPointer>() )
      // 	;

      // class_< HyperElasticUP3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "HyperElasticUP3DLaw",
      // 	  init<>() )
      // 	.def( init<HyperElasticModelPointer>() )
      // 	;
    
      // class_< HyperElasticUPPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "HyperElasticUPPlaneStrain2DLaw",
      // 	  init<>() )
      // 	.def( init<HyperElasticModelPointer>() )
      // 	;

      // class_< HyperElasticUPAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "HyperElasticUPAxisymmetric2DLaw",
      // 	  init<>() )
      // 	.def( init<HyperElasticModelPointer>() )
      // 	;    
      

      //plasticity models

      //isotropic linear elastic plasticity laws
      // class_< LinearElasticPlastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // ( "LinearElasticPlastic3DLaw",
      //   init<>() )
      //   .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
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
      // class_<HyperElasticPlastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // ( "HyperElasticPlastic3DLaw",
      //   init<>() )
      //   .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      // ;

      // class_<HyperElasticPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // ( "HyperElasticPlasticPlaneStrain2DLaw",
      //   init<>() )
      //   .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      // ;

      // class_<HyperElasticPlasticAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // ( "HyperElasticPlasticAxisym2DLaw",
      //   init<>() )
      //   .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      // ;

      // class_<IncompressibleHyperElasticPlastic3DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // ( "IncompressibleHyperElasticPlastic3DLaw",
      //   init<>() )
      //   .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      // ;

      // class_<IncompressibleHyperElasticPlasticPlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // ( "IncompressibleHyperElasticPlasticPlaneStrain2DLaw",
      //   init<>() )
      //   .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      // ;

      // class_<IncompressibleHyperElasticPlasticAxisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // ( "IncompressibleHyperElasticPlasticAxisymmetric2DLaw",
      //   init<>() )
      //   .def( init<FlowRulePointer, YieldCriterionPointer, HardeningLawPointer>() )
      // ;

      //plasticity J2 specilization laws 
      // class_<HyperElasticPlasticJ23DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "HyperElasticPlasticJ23DLaw",
      // 	  init<>() )
      // 	;

      // class_<HyperElasticPlasticJ2PlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "HyperElasticPlasticJ2PlaneStrain2DLaw",
      // 	  init<>() )
      // 	;

      // class_<HyperElasticPlasticJ2Axisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "HyperElasticPlasticJ2Axisymmetric2DLaw",
      // 	  init<>() )
      // 	;

      // class_<IncompressibleHyperElasticPlasticJ23DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "IncompressibleHyperElasticPlasticJ23DLaw",
      // 	  init<>() )
      // 	;

      // class_<IncompressibleHyperElasticPlasticJ2PlaneStrain2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "IncompressibleHyperElasticPlasticJ2PlaneStrain2DLaw",
      // 	  init<>() )
      // 	;

      // class_<IncompressibleHyperElasticPlasticJ2Axisymmetric2DLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
      // 	( "IncompressibleHyperElasticPlasticJ2Axisymmetric2DLaw",
      // 	  init<>() )
      // 	;


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

    }

  }  // namespace Python.
}  // namespace Kratos.
