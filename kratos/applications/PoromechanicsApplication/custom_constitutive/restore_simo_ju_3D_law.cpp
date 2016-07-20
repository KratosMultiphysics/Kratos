//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/restore_simo_ju_3D_law.hpp"

namespace Kratos
{

//Default Constructor
RestoreSimoJu3DLaw::RestoreSimoJu3DLaw() : IsotropicDamageSimoJu3DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new RestoreDamageFlowRule(mpYieldCriterion) );
}

//----------------------------------------------------------------------------------------

//Second Constructor
RestoreSimoJu3DLaw::RestoreSimoJu3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : IsotropicDamageSimoJu3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
RestoreSimoJu3DLaw::RestoreSimoJu3DLaw(const RestoreSimoJu3DLaw& rOther) : IsotropicDamageSimoJu3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
RestoreSimoJu3DLaw::~RestoreSimoJu3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int RestoreSimoJu3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo)
{
    int ierr = IsotropicDamageSimoJu3DLaw::Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if ( IS_CONVERGED.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"IS_CONVERGED has Key zero! (check if the application is correctly registered", "" )
    
    return ierr;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer RestoreSimoJu3DLaw::Clone() const
{
    RestoreSimoJu3DLaw::Pointer p_clone(new RestoreSimoJu3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void RestoreSimoJu3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{    
    // Check
    rValues.CheckAllParameters();
    
    // Get values for the constitutive law
    Flags& Options = rValues.GetOptions();
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    Vector& rStrainVector = rValues.GetStrainVector();

    // Initialize main variables
    // ReturnMappingVariables
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize();
    // Strain and Stress matrices
    Matrix AuxMatrix = MathUtils<double>::StrainVectorToTensor(rStrainVector);
    const unsigned int Dim = AuxMatrix.size1();
    ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
    noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
    ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);
    // CharacteristicSize
    double CharacteristicSize = 1.0;
    this->CalculateCharacteristicSize(CharacteristicSize,rValues.GetElementGeometry());
    ReturnMappingVariables.CharacteristicSize = CharacteristicSize;
    
    // LinearElasticMatrix
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    const unsigned int VoigtSize = rStrainVector.size();
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);
    
    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {
        if(Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS))
        {
            // COMPUTE_CONSTITUTIVE_TENSOR
            Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            Vector EffectiveStressVector(VoigtSize);

            this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,EffectiveStressVector,LinearElasticMatrix,rStrainVector);

            this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
        }
        else
        {
            // COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
            Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            Vector& rStressVector = rValues.GetStressVector();
            
            this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
            
            this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
        }
    }
    else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
    {
        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();
        
        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
    }
}

//----------------------------------------------------------------------------------------

void RestoreSimoJu3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{    
    // Check
    rValues.CheckAllParameters();
    
    // Get values for the constitutive law
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    Vector& rStrainVector = rValues.GetStrainVector();
    Vector& rStressVector = rValues.GetStressVector();

    // Initialize main variables
    // ReturnMappingVariables
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize();
    // Strain and Stress matrices
    Matrix AuxMatrix = MathUtils<double>::StrainVectorToTensor(rStrainVector);
    const unsigned int Dim = AuxMatrix.size1();
    ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
    noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
    ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);
    // CharacteristicSize
    double CharacteristicSize = 1.0;
    this->CalculateCharacteristicSize(CharacteristicSize,rValues.GetElementGeometry());
    ReturnMappingVariables.CharacteristicSize = CharacteristicSize;
    
    // LinearElasticMatrix
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    const unsigned int VoigtSize = rStrainVector.size();
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);
    
    if(rValues.GetProcessInfo()[IS_CONVERGED]==true) //Convergence is achieved. Save equilibrium state variable
    {
        ReturnMappingVariables.Options.Set(FlowRule::RETURN_MAPPING_COMPUTED,false); // Restore sate variable = false
        
        this->UpdateInternalStateVariables(ReturnMappingVariables,rStressVector,LinearElasticMatrix,rStrainVector);
    }
    else // No convergence is achieved. Restore state variable to equilibrium
    {
        ReturnMappingVariables.Options.Set(FlowRule::RETURN_MAPPING_COMPUTED,true); // Restore sate variable = true
        
        this->UpdateInternalStateVariables(ReturnMappingVariables,rStressVector,LinearElasticMatrix,rStrainVector);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double& RestoreSimoJu3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if (rThisVariable==DAMAGE_VARIABLE)
    {
        const FlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
        rValue=InternalVariables.DeltaPlasticStrain;
    }
    else if(rThisVariable==STATE_VARIABLE)
    {
        const FlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
        rValue=InternalVariables.EquivalentPlasticStrain;
    }

    return( rValue );
}

//----------------------------------------------------------------------------------------

void RestoreSimoJu3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                        const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == STATE_VARIABLE)
    {
        FlowRule::RadialReturnVariables ReturnMappingVariables;
        ReturnMappingVariables.TrialStateFunction = rValue;
        
        FlowRule::PlasticFactors ScalingFactors;
        
        mpFlowRule->CalculateScalingFactors(ReturnMappingVariables,ScalingFactors);
    }
}

} // Namespace Kratos
