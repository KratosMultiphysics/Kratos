//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    JMCarbonell
//					 (adapted to Particle Mechanics by Peter Wilson)
//

// System includes

// External includes

// Project includes

#include "custom_constitutive/johnson_cook_thermal_plastic_3D_law.hpp"

#include "particle_mechanics_application_variables.h"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

    JohnsonCookThermalPlastic3DLaw::JohnsonCookThermalPlastic3DLaw()
  : HyperElastic3DLaw()
  {
    mpHardeningLaw   = ParticleHardeningLaw::Pointer( new JohnsonCookThermalHardeningLaw() );
    mpYieldCriterion = ParticleYieldCriterion::Pointer( new MisesHuberThermalYieldCriterion(mpHardeningLaw) );
    mpFlowRule       = ParticleFlowRule::Pointer( new NonLinearRateDependentPlasticFlowRule(mpYieldCriterion) );
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

    JohnsonCookThermalPlastic3DLaw::JohnsonCookThermalPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
  {
    mpHardeningLaw    =  pHardeningLaw;
    mpYieldCriterion  = ParticleYieldCriterion::Pointer( new MisesHuberThermalYieldCriterion(mpHardeningLaw) );
    mpFlowRule        =  pFlowRule;
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

    JohnsonCookThermalPlastic3DLaw::JohnsonCookThermalPlastic3DLaw(const JohnsonCookThermalPlastic3DLaw& rOther)
  : HyperElastic3DLaw(rOther)
  {

  }

  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer JohnsonCookThermalPlastic3DLaw::Clone() const
  {
    return Kratos::make_shared<JohnsonCookThermalPlastic3DLaw>(*this);
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  JohnsonCookThermalPlastic3DLaw::~JohnsonCookThermalPlastic3DLaw()
  {
  }

  //******************************* COMPUTE DOMAIN TEMPERATURE  ************************
  //************************************************************************************


  double & JohnsonCookThermalPlastic3DLaw::CalculateDomainTemperature (const MaterialResponseVariables & rElasticVariables,
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

  bool JohnsonCookThermalPlastic3DLaw::Has( const Variable<double>& rThisVariable )
  {
    if(rThisVariable == DELTA_PLASTIC_DISSIPATION || rThisVariable == PLASTIC_DISSIPATION )
      return true;

    return false;
  }

  //***********************GET VALUE: DOUBLE - VECTOR - MATRIX**************************
  //************************************************************************************

  double& JohnsonCookThermalPlastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
  {
    if (rThisVariable==PLASTIC_STRAIN)
      {
	const ParticleFlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
	rValue=InternalVariables.EquivalentPlasticStrain;
      }

    if (rThisVariable==DELTA_PLASTIC_STRAIN)
      {
	const ParticleFlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
	rValue=InternalVariables.DeltaPlasticStrain;
      }


    if (rThisVariable==PLASTIC_DISSIPATION)
      {
	const ParticleFlowRule::ThermalVariables& ThermalVariables = mpFlowRule->GetThermalVariables();
	rValue=ThermalVariables.PlasticDissipation;
      }

    if (rThisVariable==DELTA_PLASTIC_DISSIPATION)
      {
	const ParticleFlowRule::ThermalVariables& ThermalVariables = mpFlowRule->GetThermalVariables();
	rValue=ThermalVariables.DeltaPlasticDissipation;
      }

    return( rValue );
  }


} // Namespace Kratos
