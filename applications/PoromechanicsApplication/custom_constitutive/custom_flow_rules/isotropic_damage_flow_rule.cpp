//
//   Project Name:        KratosPoromechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_flow_rules/isotropic_damage_flow_rule.hpp"

#include "poromechanics_application_variables.h"

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
  return Kratos::make_shared<IsotropicDamageFlowRule>(*this);
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
    FlowRule::InitializeMaterial(pYieldCriterion,pHardeningLaw,rMaterialProperties);

    //EquivalentPlasticStrain is the maximum historical equivalent strain, and DAMAGE_THRESHOLD is the minimum default value
    mInternalVariables.EquivalentPlasticStrain = rMaterialProperties[DAMAGE_THRESHOLD];
    mInternalVariables.EquivalentPlasticStrainOld = mInternalVariables.EquivalentPlasticStrain;
}


//*************************** CALCULATE RETURN MAPPING *******************************
//************************************************************************************

bool IsotropicDamageFlowRule::CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient,
                                                    Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{
    //Compute Damage variable
    bool Tangent = this->CalculateInternalVariables( rReturnMappingVariables );

    //Compute Damaged Stresses
    noalias(rStressMatrix) = (1.0-rReturnMappingVariables.TrialStateFunction)*rReturnMappingVariables.TrialIsoStressMatrix; // S = (1-d)*Se

    return Tangent;
}

bool IsotropicDamageFlowRule::CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix )
{
    //Compute Damage variable
    bool Tangent = this->CalculateInternalVariables( rReturnMappingVariables );

    //Compute Damaged Stresses
    noalias(rIsoStressMatrix) = (1.0-rReturnMappingVariables.TrialStateFunction)*rReturnMappingVariables.TrialIsoStressMatrix; // S = (1-d)*Se

	return Tangent;
}


//************************* COMPUTE TANGENT CONSTITUTIVE MATRIX **********************
//************************************************************************************

void IsotropicDamageFlowRule::ComputeElastoPlasticTangentMatrix( const RadialReturnVariables& rReturnMappingVariables, const Matrix& rElasticLeftCauchyGreen,
                                                                 const double& rAlpha, Matrix& rElastoPlasticMatrix)
{
    const unsigned int VoigtSize = rElasticLeftCauchyGreen.size1();

    Vector EffectiveStressVector(VoigtSize);
    noalias(EffectiveStressVector) = MathUtils<double>::StressTensorToVector( rReturnMappingVariables.TrialIsoStressMatrix, 0 );

    Vector EquivalentStrainDerivative(VoigtSize);
    this->CalculateEquivalentStrainDerivative(EquivalentStrainDerivative, rReturnMappingVariables, rElasticLeftCauchyGreen);

    double DamageDerivative;
    YieldCriterion::Parameters YieldCriterionParameters;
    YieldCriterionParameters.SetCharacteristicSize(rReturnMappingVariables.CharacteristicSize);
    YieldCriterionParameters.SetDeltaGamma(mInternalVariables.EquivalentPlasticStrain); //Maximum historical equivalent strain
    mpYieldCriterion->CalculateDeltaStateFunction(DamageDerivative, YieldCriterionParameters); //Compute damage derivative

    noalias(rElastoPlasticMatrix) += - rAlpha*DamageDerivative*outer_prod(EffectiveStressVector,EquivalentStrainDerivative);
}


//*************************** UPDATE INTERNAL VARIABLES ******************************
//************************************************************************************

bool IsotropicDamageFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{
    bool Restore;
    YieldCriterion::Parameters YieldCriterionParameters;

    // Updates internal variables depending on the convergence of the solution
    if( rReturnMappingVariables.Options.IsNot(RETURN_MAPPING_COMPUTED) ) // Convergence was achieved
    {
        //Compute the new equivalent strain
        double NewEquivalentStrain;
        YieldCriterionParameters.SetStrainMatrix(rReturnMappingVariables.StrainMatrix);
        YieldCriterionParameters.SetStressMatrix(rReturnMappingVariables.TrialIsoStressMatrix);
        mpYieldCriterion->CalculateYieldCondition(NewEquivalentStrain, YieldCriterionParameters);

        //Update maximum historical equivalent strain
        if(NewEquivalentStrain >= mInternalVariables.EquivalentPlasticStrain)
        {
            mInternalVariables.EquivalentPlasticStrain = NewEquivalentStrain;
        }

        Restore = false;
    }
    else // There was no convergence
    {
        Restore = true;
    }

    //Compute Damage variable from the internal historical variable (the maximum equivalent strain)
    YieldCriterionParameters.SetCharacteristicSize(rReturnMappingVariables.CharacteristicSize);
    YieldCriterionParameters.SetDeltaGamma(mInternalVariables.EquivalentPlasticStrain); //internal historical variable
    mpYieldCriterion->CalculateStateFunction(rReturnMappingVariables.TrialStateFunction, YieldCriterionParameters);

    //mInternalVariables.DeltaPlasticStrain is the damage variable that will be printed in post-process
    mInternalVariables.DeltaPlasticStrain = rReturnMappingVariables.TrialStateFunction;

    return Restore;
}


//************************** CALCULATE INTERNAL VARIABLES ****************************
//************************************************************************************

bool IsotropicDamageFlowRule::CalculateInternalVariables(RadialReturnVariables& rReturnMappingVariables)
{
    //Compute the new equivalent strain
    double NewEquivalentStrain;
    YieldCriterion::Parameters YieldCriterionParameters;
    YieldCriterionParameters.SetStrainMatrix(rReturnMappingVariables.StrainMatrix);
    YieldCriterionParameters.SetStressMatrix(rReturnMappingVariables.TrialIsoStressMatrix);
    mpYieldCriterion->CalculateYieldCondition(NewEquivalentStrain, YieldCriterionParameters);

    //Update maximum historical equivalent strain and set whether the tangent matrix needs to be computed
    bool Tangent;
    if(NewEquivalentStrain >= mInternalVariables.EquivalentPlasticStrain)
    {
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,true); //loading with growing damage
        Tangent = true;
    }
    else
    {
        rReturnMappingVariables.Options.Set(PLASTIC_REGION,false); //elastic loading or unloading
        Tangent = false;
    }

    //Compute Damage variable from the internal historical variable (the maximum equivalent strain)
    YieldCriterionParameters.SetCharacteristicSize(rReturnMappingVariables.CharacteristicSize);
    YieldCriterionParameters.SetDeltaGamma(mInternalVariables.EquivalentPlasticStrain); //internal historical variable
    mpYieldCriterion->CalculateStateFunction(rReturnMappingVariables.TrialStateFunction, YieldCriterionParameters);

    //mInternalVariables.DeltaPlasticStrain is the damage variable that will be printed in post-process
    mInternalVariables.DeltaPlasticStrain = rReturnMappingVariables.TrialStateFunction;

    return Tangent;
}


//************************ CALCULATE EQUIVALENT STRAIN DERIVATIVE ********************
//************************************************************************************

void IsotropicDamageFlowRule::CalculateEquivalentStrainDerivative(Vector& rEquivalentStrainDerivative, const RadialReturnVariables& ReturnMappingVariables,
                                                                    const Matrix& LinearElasticMatrix)
{
    //The derivative of the equivalent strain with respect to the strain vector is obtained through the perturbation method

    const unsigned int VoigtSize = LinearElasticMatrix.size1();
    unsigned int Dim = ReturnMappingVariables.StrainMatrix.size1();

    if(rEquivalentStrainDerivative.size() != VoigtSize)
        rEquivalentStrainDerivative.resize(VoigtSize,false);

    // Necessary variables
    Vector StrainVector(VoigtSize);
    noalias(StrainVector) = MathUtils<double>::StrainTensorToVector( ReturnMappingVariables.StrainMatrix, 0 );
    Vector StrainVector_p(VoigtSize);
    Vector StressVector_p(VoigtSize);
    Matrix StrainMatrix_p(Dim,Dim);
    Matrix StressMatrix_p(Dim,Dim);
    double EquivalentStrain_f, EquivalentStrain_b;
    YieldCriterion::Parameters YieldCriterionParameters;
    YieldCriterionParameters.SetStrainMatrix(StrainMatrix_p); // YieldCriterionParameters stores the direction of StrainMatrix_p (const Matrix* mpStrainMatrix)
    YieldCriterionParameters.SetStressMatrix(StressMatrix_p); // YieldCriterionParameters stores the direction of StressMatrix_p (const Matrix* mpStressMatrix)

    //Compute the strains perturbations in each direction of the vector
    Vector PerturbationStrainVector(VoigtSize);
    IsotropicDamageUtilities::ComputePerturbationVector(PerturbationStrainVector,StrainVector);

    for(unsigned int i = 0; i < VoigtSize; i++)
    {
        //Forward perturbed equivalent strain
        noalias(StrainVector_p) = StrainVector;
        StrainVector_p[i] = StrainVector_p[i] + PerturbationStrainVector[i];
        noalias(StrainMatrix_p) = MathUtils<double>::StrainVectorToTensor(StrainVector_p);
        noalias(StressVector_p) = prod(LinearElasticMatrix, StrainVector_p);
        noalias(StressMatrix_p) = MathUtils<double>::StressVectorToTensor(StressVector_p);
        mpYieldCriterion->CalculateYieldCondition(EquivalentStrain_f, YieldCriterionParameters);

        //Backward perturbed equivalent strain
        noalias(StrainVector_p) = StrainVector;
        StrainVector_p[i] = StrainVector_p[i] - PerturbationStrainVector[i];
        noalias(StrainMatrix_p) = MathUtils<double>::StrainVectorToTensor(StrainVector_p);
        noalias(StressVector_p) = prod(LinearElasticMatrix, StrainVector_p);
        noalias(StressMatrix_p) = MathUtils<double>::StressVectorToTensor(StressVector_p);
        mpYieldCriterion->CalculateYieldCondition(EquivalentStrain_b, YieldCriterionParameters);

        rEquivalentStrainDerivative[i] = (EquivalentStrain_f - EquivalentStrain_b)/(2.0*PerturbationStrainVector[i]);
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
