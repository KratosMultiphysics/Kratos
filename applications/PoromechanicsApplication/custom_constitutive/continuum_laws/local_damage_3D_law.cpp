//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/continuum_laws/local_damage_3D_law.hpp"

namespace Kratos
{

//Default Constructor
LocalDamage3DLaw::LocalDamage3DLaw() : LinearElasticPlastic3DLaw() {}

//----------------------------------------------------------------------------------------

//Second Constructor
LocalDamage3DLaw::LocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : LinearElasticPlastic3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
LocalDamage3DLaw::LocalDamage3DLaw(const LocalDamage3DLaw& rOther) : LinearElasticPlastic3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
LocalDamage3DLaw::~LocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int LocalDamage3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = HyperElasticPlastic3DLaw::Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if ( IS_CONVERGED.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"IS_CONVERGED Key is 0. Check if all applications were correctly registered.", "" )
    if ( STATE_VARIABLE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"STATE_VARIABLE Key is 0. Check if all applications were correctly registered.", "" )

    return ierr;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer LocalDamage3DLaw::Clone() const
{
    LocalDamage3DLaw::Pointer p_clone(new LocalDamage3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LocalDamage3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{    
    // Check
    rValues.CheckAllParameters();
    
    // Get values for the constitutive law
    Flags& Options = rValues.GetOptions();
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    Vector& rStrainVector = rValues.GetStrainVector();

    // Initialize main variables

    // LinearElasticMatrix
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    const unsigned int VoigtSize = rStrainVector.size();
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);
    
    // ReturnMappingVariables
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize();
    // Strain and Stress matrices
    const unsigned int Dim = this->WorkingSpaceDimension();
    Matrix AuxMatrix(Dim,Dim);
    noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(rStrainVector);
    ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
    noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
    ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);
    // CharacteristicSize
    double CharacteristicSize = 1.0;
    this->CalculateCharacteristicSize(CharacteristicSize,rValues.GetElementGeometry());
    ReturnMappingVariables.CharacteristicSize = CharacteristicSize;

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

void LocalDamage3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{    
    // Check
    rValues.CheckAllParameters();
    
    // Get values for the constitutive law
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    Vector& rStrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = rStrainVector.size();
    Vector EffectiveStressVector(VoigtSize);
    
    // Initialize main variables

    // LinearElasticMatrix
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);
    
    // ReturnMappingVariables
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize();
    // Strain and Stress matrices
    const unsigned int Dim = this->WorkingSpaceDimension();
    ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
    noalias(ReturnMappingVariables.StrainMatrix) = MathUtils<double>::StrainVectorToTensor(rStrainVector);
    ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);
    // CharacteristicSize
    double CharacteristicSize = 1.0;
    this->CalculateCharacteristicSize(CharacteristicSize,rValues.GetElementGeometry());
    ReturnMappingVariables.CharacteristicSize = CharacteristicSize;

    if(rValues.GetProcessInfo()[IS_CONVERGED]==true) //Convergence is achieved. Save equilibrium state variable
    {
        ReturnMappingVariables.Options.Set(FlowRule::RETURN_MAPPING_COMPUTED,false); // Restore sate variable = false
        
        this->UpdateInternalStateVariables(ReturnMappingVariables,EffectiveStressVector,LinearElasticMatrix,rStrainVector);
    }
    else // No convergence is achieved. Restore state variable to equilibrium
    {
        ReturnMappingVariables.Options.Set(FlowRule::RETURN_MAPPING_COMPUTED,true); // Restore sate variable = true
        
        this->UpdateInternalStateVariables(ReturnMappingVariables,EffectiveStressVector,LinearElasticMatrix,rStrainVector);
    }

    if(rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRESS))
    {
        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();

        this->UpdateStressVector(rStressVector,ReturnMappingVariables,EffectiveStressVector);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double& LocalDamage3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if(rThisVariable==STATE_VARIABLE)
    {
        const FlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
        rValue=InternalVariables.EquivalentPlasticStrain;
    }
    else
    {
        LinearElasticPlastic3DLaw::GetValue(rThisVariable,rValue);
    }

    return( rValue );
}

//----------------------------------------------------------------------------------------

void LocalDamage3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
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
