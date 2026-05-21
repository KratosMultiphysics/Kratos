//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/continuum_laws/custom_flow_rules/nonlocal_damage_flow_rule.hpp"

namespace Kratos
{

//Default Constructor
NonlocalDamageFlowRule::NonlocalDamageFlowRule() : LocalDamageFlowRule() {}

//----------------------------------------------------------------------------------------

//Second Constructor
NonlocalDamageFlowRule::NonlocalDamageFlowRule(YieldCriterionPointer pYieldCriterion)
	:LocalDamageFlowRule(pYieldCriterion) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
NonlocalDamageFlowRule::NonlocalDamageFlowRule(NonlocalDamageFlowRule const& rOther)
	:LocalDamageFlowRule(rOther) {}

//----------------------------------------------------------------------------------------

//Assignment Operator
NonlocalDamageFlowRule& NonlocalDamageFlowRule::operator=(NonlocalDamageFlowRule const& rOther)
{
   LocalDamageFlowRule::operator=(rOther);
   return *this;
}

//----------------------------------------------------------------------------------------

//Destructor
NonlocalDamageFlowRule::~NonlocalDamageFlowRule() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

FlowRule::Pointer NonlocalDamageFlowRule::Clone() const
{
  FlowRule::Pointer p_clone(new NonlocalDamageFlowRule(*this));
  return p_clone;
}

//----------------------------------------------------------------------------------------

void NonlocalDamageFlowRule::InitializeMaterial (YieldCriterionPointer& pYieldCriterion, HardeningLawPointer& pHardeningLaw, const Properties& rMaterialProperties)
{
    IsotropicDamageFlowRule::InitializeMaterial(pYieldCriterion,pHardeningLaw,rMaterialProperties);
    
    mThermalVariables.PlasticDissipation = 0.0; // Initial value for the NonlocalEquivalentStrain
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bool NonlocalDamageFlowRule::CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, 
                                                        Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{
    //Compute LocalEquivalentStrain and Damage
    this->CalculateLocalInternalVariables( rReturnMappingVariables );
    
    //Compute Damaged Stresses
    noalias(rStressMatrix) = (1.0-rReturnMappingVariables.TrialStateFunction)*rReturnMappingVariables.TrialIsoStressMatrix; // S = (1-d)*Se
    
    return false;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bool NonlocalDamageFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{
    bool Restore;
    
    // Updates internal variables depending on the convergence of the solution
    if( rReturnMappingVariables.Options.IsNot(RETURN_MAPPING_COMPUTED) ) // Convergence was achieved
    {
        // Set New NonlocalEquivalentStrain (mThermalVariables.PlasticDissipation)
        mThermalVariables.PlasticDissipation = rReturnMappingVariables.NormIsochoricStress;
                
        //Update maximum historical equivalent strain
        if(mThermalVariables.PlasticDissipation >= mInternalVariables.EquivalentPlasticStrain)
        {
            mInternalVariables.EquivalentPlasticStrain = mThermalVariables.PlasticDissipation;
        }
        
        Restore = false;
    }
    else // There was no convergence
    {   
        Restore = true;
    }

    //Compute Damage variable from the internal historical variable (the maximum equivalent strain)
    YieldCriterion::Parameters YieldCriterionParameters;
    YieldCriterionParameters.SetCharacteristicSize(rReturnMappingVariables.CharacteristicSize);
    YieldCriterionParameters.SetDeltaGamma(mInternalVariables.EquivalentPlasticStrain); //internal historical variable
    mpYieldCriterion->CalculateStateFunction(rReturnMappingVariables.TrialStateFunction, YieldCriterionParameters);

    //mInternalVariables.DeltaPlasticStrain is the damage variable that will be printed in post-process
    mInternalVariables.DeltaPlasticStrain = rReturnMappingVariables.TrialStateFunction;
    
    return Restore;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void NonlocalDamageFlowRule::CalculateLocalInternalVariables(RadialReturnVariables& rReturnMappingVariables)
{
    //Compute the new local equivalent strain (mThermalVariables.DeltaPlasticDissipation)
    YieldCriterion::Parameters YieldCriterionParameters;
    YieldCriterionParameters.SetStrainMatrix(rReturnMappingVariables.StrainMatrix);
    YieldCriterionParameters.SetStressMatrix(rReturnMappingVariables.TrialIsoStressMatrix);
    mpYieldCriterion->CalculateYieldCondition(mThermalVariables.DeltaPlasticDissipation, YieldCriterionParameters);

    rReturnMappingVariables.Options.Set(PLASTIC_REGION,false); //elastic loading or unloading
        
    //Compute Damage variable from the internal historical variable (the maximum equivalent strain)
    YieldCriterionParameters.SetCharacteristicSize(rReturnMappingVariables.CharacteristicSize);
    YieldCriterionParameters.SetDeltaGamma(mInternalVariables.EquivalentPlasticStrain); //internal historical variable
    mpYieldCriterion->CalculateStateFunction(rReturnMappingVariables.TrialStateFunction, YieldCriterionParameters);
    
    //mInternalVariables.DeltaPlasticStrain is the damage variable that will be printed in post-process
    mInternalVariables.DeltaPlasticStrain = rReturnMappingVariables.TrialStateFunction;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bool NonlocalDamageFlowRule::CalculateInternalVariables(RadialReturnVariables& rReturnMappingVariables)
{
    // Set New NonlocalEquivalentStrain (mThermalVariables.PlasticDissipation)
    mThermalVariables.PlasticDissipation = rReturnMappingVariables.NormIsochoricStress;
    
    //Update maximum historical equivalent strain and set whether the tangent matrix needs to be computed
    bool Tangent;
    if(mThermalVariables.PlasticDissipation >= mInternalVariables.EquivalentPlasticStrain)
    {   
        //rReturnMappingVariables.Options.Set(PLASTIC_REGION,true); //loading with growing damage
        //Tangent = true;
        //TODO: For the moment we use the secant constitutive tensor
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,false); //elastic loading or unloading
        Tangent = false;
    }
    else
    {
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,false); //elastic loading or unloading
        Tangent = false;
    }
    
    //Compute Damage variable from the internal historical variable (the maximum equivalent strain)
    YieldCriterion::Parameters YieldCriterionParameters;    
    YieldCriterionParameters.SetCharacteristicSize(rReturnMappingVariables.CharacteristicSize);
    YieldCriterionParameters.SetDeltaGamma(mInternalVariables.EquivalentPlasticStrain); //internal historical variable
    mpYieldCriterion->CalculateStateFunction(rReturnMappingVariables.TrialStateFunction, YieldCriterionParameters);
    
    //mInternalVariables.DeltaPlasticStrain is the damage variable that will be printed in post-process
    mInternalVariables.DeltaPlasticStrain = rReturnMappingVariables.TrialStateFunction;
    
    return Tangent;
}

} // Namespace Kratos
