//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "solid_mechanics_application_variables.h"
#include "custom_constitutive/custom_flow_rules/isotropic_damage_flow_rule.hpp"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

IsotropicDamageFlowRule::IsotropicDamageFlowRule() : FlowRule()
{
   
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

IsotropicDamageFlowRule::IsotropicDamageFlowRule(YieldCriterionPointer pYieldCriterion)
	:FlowRule(pYieldCriterion)
{
   
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

IsotropicDamageFlowRule& IsotropicDamageFlowRule::operator=(IsotropicDamageFlowRule const& rOther)
{
   FlowRule::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

IsotropicDamageFlowRule::IsotropicDamageFlowRule(IsotropicDamageFlowRule const& rOther)
	:FlowRule(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

FlowRule::Pointer IsotropicDamageFlowRule::Clone() const
{
  FlowRule::Pointer p_clone(new IsotropicDamageFlowRule(*this));
  return p_clone;
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

IsotropicDamageFlowRule::~IsotropicDamageFlowRule()
{
}

/// Operations.


//*************************** INITIALIZE MATERIAL ************************************
//************************************************************************************
void IsotropicDamageFlowRule::InitializeMaterial (YieldCriterionPointer& pYieldCriterion, HardeningLawPointer& pHardeningLaw, const Properties& rMaterialProperties)
{
    //set yield criterion
    mpYieldCriterion = pYieldCriterion;
    mpYieldCriterion->InitializeMaterial(pHardeningLaw, rMaterialProperties);	

    //initialize material variables
    mInternalVariables.clear();
    mThermalVariables.clear();
    
    //EquivalentPlasticStrain is the maximum historical equivalent strain, and DAMAGE_THRESHOLD is the minimum default value
    mInternalVariables.EquivalentPlasticStrain = rMaterialProperties[DAMAGE_THRESHOLD];
}


//*************************** CALCULATE RETURN MAPPING *******************************
//************************************************************************************

bool IsotropicDamageFlowRule::CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix )
{
    //Compute the new equivalent strain
    double NewEquivalentStrain;
    YieldCriterion::Parameters YieldCriterionParameters;
    YieldCriterionParameters.SetStrainMatrix(rReturnMappingVariables.StrainMatrix);
    YieldCriterionParameters.SetStressMatrix(rReturnMappingVariables.TrialIsoStressMatrix);
    mpYieldCriterion->CalculateYieldCondition(NewEquivalentStrain, YieldCriterionParameters);

    //Update maximum historical equivalent strain
    if(NewEquivalentStrain >= mInternalVariables.EquivalentPlasticStrain)
    {
        mInternalVariables.EquivalentPlasticStrain = NewEquivalentStrain;
        
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,true); //loading with growing damage
    }
    else
    {
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,false); //elastic loading or unloading
    }
    
    //Compute Damage variable from the internal historical variable (the maximum equivalent strain)
    YieldCriterionParameters.SetCharacteristicSize(rReturnMappingVariables.CharacteristicSize);
    YieldCriterionParameters.SetDeltaGamma(mInternalVariables.EquivalentPlasticStrain); //internal historical variable
    mpYieldCriterion->CalculateStateFunction(rReturnMappingVariables.TrialStateFunction, YieldCriterionParameters);

    //Compute Damaged Stresses
    noalias(rIsoStressMatrix) = (1-rReturnMappingVariables.TrialStateFunction)*rReturnMappingVariables.TrialIsoStressMatrix;
    
	return true;
}


//************************* COMPUTE TANGENT CONSTITUTIVE MATRIX **********************
//************************************************************************************

void IsotropicDamageFlowRule::ComputeElastoPlasticTangentMatrix( const RadialReturnVariables& rReturnMappingVariables, const Matrix& rElasticLeftCauchyGreen, 
                                                                 const double& rAlpha, Matrix& rElastoPlasticMatrix)
{   
    Vector StrainVector = MathUtils<double>::StrainTensorToVector( rReturnMappingVariables.StrainMatrix, 0 ); //For plane state or 3D cases
    Vector StressVector = MathUtils<double>::StressTensorToVector( rReturnMappingVariables.TrialIsoStressMatrix, 0 ); //For plane state or 3D cases
    
    Vector EquivalentStrainDerivative = ZeroVector(StrainVector.size());
    this->CalculateEquivalentStrainDerivative(EquivalentStrainDerivative, StrainVector, rElasticLeftCauchyGreen);
    
    double DamageDerivative;
    YieldCriterion::Parameters YieldCriterionParameters;
    YieldCriterionParameters.SetCharacteristicSize(rReturnMappingVariables.CharacteristicSize);
    YieldCriterionParameters.SetDeltaGamma(mInternalVariables.EquivalentPlasticStrain); //Maximum historical equivalent strain
    mpYieldCriterion->CalculateDeltaStateFunction(DamageDerivative, YieldCriterionParameters); //Compute damage derivative
    
    noalias(rElastoPlasticMatrix) = - rAlpha*DamageDerivative*outer_prod(StressVector,EquivalentStrainDerivative);
}


//*************************** UPDATE DAMAGE VARIABLE *********************************
//************************************************************************************

bool IsotropicDamageFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{
    //mInternalVariables.DeltaPlasticStrain is the damage variable that will be printed in post-process
    mInternalVariables.DeltaPlasticStrain = rReturnMappingVariables.TrialStateFunction;
    
	return true;
}


//************************ CALCULATE EQUIVALENT STRAIN DERIVATIVE ********************
//************************************************************************************


void IsotropicDamageFlowRule::CalculateEquivalentStrainDerivative(Vector& rEquivalentStrainDerivative, const Vector& rStrainVector, const Matrix& rLinearElasticMatrix)
{
    //The derivative of the equivalent strain with respect to the strain vector is obtained through the perturbation method
    
    Vector StrainVector_p = ZeroVector(rStrainVector.size());
    Vector StressVector_p = ZeroVector(rStrainVector.size());
    Matrix StrainMatrix_p, StressMatrix_p;
    double EquivalentStrain_f, EquivalentStrain_b;
    YieldCriterion::Parameters YieldCriterionParameters;
    
    //Compute the strains perturbations in each direction of the vector
    IsotropicDamageUtilities DamageUtilities;
    Vector PerturbationStrainVector = DamageUtilities.PerturbationVector(rStrainVector);

    for(unsigned int i=0;i<rStrainVector.size();i++)
    {
        //Forward perturbed equivalent strain
        noalias(StrainVector_p) = rStrainVector;
        StrainVector_p[i]       = StrainVector_p[i] + PerturbationStrainVector[i];
        noalias(StressVector_p) = prod(rLinearElasticMatrix, StrainVector_p);
        StrainMatrix_p = MathUtils<double>::StrainVectorToTensor(StrainVector_p);
        StressMatrix_p = MathUtils<double>::StressVectorToTensor(StressVector_p);
        YieldCriterionParameters.SetStrainMatrix(StrainMatrix_p);
        YieldCriterionParameters.SetStressMatrix(StressMatrix_p);
        mpYieldCriterion->CalculateYieldCondition(EquivalentStrain_f, YieldCriterionParameters);
        
        //Backward perturbed equivalent strain
        noalias(StrainVector_p) = rStrainVector;
        StrainVector_p[i]       = StrainVector_p[i] - PerturbationStrainVector[i];
        noalias(StressVector_p) = prod(rLinearElasticMatrix, StrainVector_p);
        StrainMatrix_p = MathUtils<double>::StrainVectorToTensor(StrainVector_p);
        StressMatrix_p = MathUtils<double>::StressVectorToTensor(StressVector_p);
        YieldCriterionParameters.SetStrainMatrix(StrainMatrix_p);
        YieldCriterionParameters.SetStressMatrix(StressMatrix_p);
        mpYieldCriterion->CalculateYieldCondition(EquivalentStrain_b, YieldCriterionParameters);
        
        rEquivalentStrainDerivative[i] = (EquivalentStrain_f - EquivalentStrain_b)/(2*PerturbationStrainVector[i]);
    }
}


//************************************************************************************
//************************************************************************************


void IsotropicDamageFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
}

void IsotropicDamageFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )
}


}  // namespace Kratos.
