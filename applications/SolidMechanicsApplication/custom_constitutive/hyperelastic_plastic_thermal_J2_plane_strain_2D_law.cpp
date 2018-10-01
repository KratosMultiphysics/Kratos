//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_yield_criteria/mises_huber_thermal_yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_thermal_hardening_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_thermal_J2_plane_strain_2D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticPlasticThermalJ2PlaneStrain2DLaw::HyperElasticPlasticThermalJ2PlaneStrain2DLaw()
  : HyperElasticPlasticPlaneStrain2DLaw()
  {
    mpHardeningLaw   = HardeningLaw::Pointer( new NonLinearIsotropicKinematicThermalHardeningLaw() );
    mpYieldCriterion = YieldCriterion::Pointer( new MisesHuberThermalYieldCriterion(mpHardeningLaw) );
    mpFlowRule       = FlowRule::Pointer( new NonLinearAssociativePlasticFlowRule(mpYieldCriterion) );
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticPlasticThermalJ2PlaneStrain2DLaw::HyperElasticPlasticThermalJ2PlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
  {
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  =  YieldCriterion::Pointer( new MisesHuberThermalYieldCriterion(mpHardeningLaw) );
    mpFlowRule        =  pFlowRule;
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  HyperElasticPlasticThermalJ2PlaneStrain2DLaw::HyperElasticPlasticThermalJ2PlaneStrain2DLaw(const HyperElasticPlasticThermalJ2PlaneStrain2DLaw& rOther)
  : HyperElasticPlasticPlaneStrain2DLaw(rOther)
  {

  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer HyperElasticPlasticThermalJ2PlaneStrain2DLaw::Clone() const
  {
    return Kratos::make_shared<HyperElasticPlasticThermalJ2PlaneStrain2DLaw>(*this);
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticPlasticThermalJ2PlaneStrain2DLaw::~HyperElasticPlasticThermalJ2PlaneStrain2DLaw()
  {
  }


  //******************************* COMPUTE DOMAIN TEMPERATURE  ************************
  //************************************************************************************


  double &  HyperElasticPlasticThermalJ2PlaneStrain2DLaw::CalculateDomainTemperature (const MaterialResponseVariables & rElasticVariables,
										      double & rTemperature)
  {

    //1.-Temperature from nodes
    const GeometryType& DomainGeometry = rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues = rElasticVariables.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();

    rTemperature=0;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
      {
     	rTemperature += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(TEMPERATURE);
      }


    return rTemperature;
  }

  //*******************************OPERATIONS FROM BASE CLASS***************************
  //************************************************************************************

  //***********************HAS : DOUBLE - VECTOR - MATRIX*******************************
  //************************************************************************************

  bool HyperElasticPlasticThermalJ2PlaneStrain2DLaw::Has( const Variable<double>& rThisVariable )
  {
    if(rThisVariable == DELTA_PLASTIC_DISSIPATION || rThisVariable == PLASTIC_DISSIPATION )
      return true;

    return false;
  }

  //***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  double& HyperElasticPlasticThermalJ2PlaneStrain2DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
  {
    if (rThisVariable==PLASTIC_STRAIN)
      {
	const FlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
	rValue=InternalVariables.EquivalentPlasticStrain;
      }

    if (rThisVariable==DELTA_PLASTIC_STRAIN)
      {
	const FlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
	rValue=InternalVariables.DeltaPlasticStrain;
      }


    if (rThisVariable==PLASTIC_DISSIPATION)
      {
	const FlowRule::ThermalVariables& ThermalVariables = mpFlowRule->GetThermalVariables();
	rValue=ThermalVariables.PlasticDissipation;
      }

    if (rThisVariable==DELTA_PLASTIC_DISSIPATION)
      {
	const FlowRule::ThermalVariables& ThermalVariables = mpFlowRule->GetThermalVariables();
	rValue=ThermalVariables.DeltaPlasticDissipation;
      }

    return( rValue );
  }


} // Namespace Kratos
