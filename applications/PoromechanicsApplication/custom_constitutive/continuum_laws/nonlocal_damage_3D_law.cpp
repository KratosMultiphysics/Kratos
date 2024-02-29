//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/continuum_laws/nonlocal_damage_3D_law.hpp"

namespace Kratos
{

//Default Constructor
NonlocalDamage3DLaw::NonlocalDamage3DLaw() : LocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------

//Second Constructor
NonlocalDamage3DLaw::NonlocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : LocalDamage3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
NonlocalDamage3DLaw::NonlocalDamage3DLaw(const NonlocalDamage3DLaw& rOther) : LocalDamage3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
NonlocalDamage3DLaw::~NonlocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int NonlocalDamage3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = LocalDamage3DLaw::Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if ( LOCAL_EQUIVALENT_STRAIN.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"LOCAL_EQUIVALENT_STRAIN Key is 0. Check if all applications were correctly registered.", "" )
    if ( NONLOCAL_EQUIVALENT_STRAIN.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"NONLOCAL_EQUIVALENT_STRAIN Key is 0. Check if all applications were correctly registered.", "" )

    return ierr;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer NonlocalDamage3DLaw::Clone() const
{
    NonlocalDamage3DLaw::Pointer p_clone(new NonlocalDamage3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------

void NonlocalDamage3DLaw::InitializeMaterial( const Properties& rMaterialProperties,
						   const GeometryType& rElementGeometry,
						   const Vector& rShapeFunctionsValues )
{
    HyperElasticPlastic3DLaw::InitializeMaterial(rMaterialProperties,rElementGeometry,rShapeFunctionsValues);
    
    mNonlocalEquivalentStrain = 0.0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void NonlocalDamage3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
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
    // CharacteristicSize (for nonlocal damage it must be 1.0)
    ReturnMappingVariables.CharacteristicSize = 1.0;

    if(Options.Is(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE)) // LOCAL QUANTITIES
    {
        if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
        {
            if(Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS))
            {
                // COMPUTE_CONSTITUTIVE_TENSOR
                Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                Vector EffectiveStressVector(VoigtSize);

                this->CalculateLocalReturnMapping(ReturnMappingVariables,AuxMatrix,EffectiveStressVector,LinearElasticMatrix,rStrainVector);

                this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
            }
            else
            {
                // COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
                Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                Vector& rStressVector = rValues.GetStressVector();
                
                this->CalculateLocalReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
                
                this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
            }
        }
        else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
        {
            // COMPUTE_STRESS
            Vector& rStressVector = rValues.GetStressVector();
            
            this->CalculateLocalReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
        }
    }
    else // NONLOCAL QUANTITIES
    {
        ReturnMappingVariables.NormIsochoricStress = mNonlocalEquivalentStrain;
        
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
}

//----------------------------------------------------------------------------------------

void NonlocalDamage3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
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
    ReturnMappingVariables.CharacteristicSize = 1.0;
    ReturnMappingVariables.NormIsochoricStress = mNonlocalEquivalentStrain;

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

double& NonlocalDamage3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if (rThisVariable==LOCAL_EQUIVALENT_STRAIN)
    {
        const FlowRule::ThermalVariables& ThermalVariables = mpFlowRule->GetThermalVariables();
        rValue=ThermalVariables.DeltaPlasticDissipation;
    }
    else
    {
        LocalDamage3DLaw::GetValue(rThisVariable,rValue);
    }

    return( rValue );
}

//----------------------------------------------------------------------------------------

void NonlocalDamage3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                        const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == NONLOCAL_EQUIVALENT_STRAIN)
    {
        mNonlocalEquivalentStrain = rValue;
    }
    else
    {
        LocalDamage3DLaw::SetValue(rThisVariable,rValue,rCurrentProcessInfo);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void NonlocalDamage3DLaw::CalculateLocalReturnMapping( FlowRule::RadialReturnVariables& rReturnMappingVariables, Matrix& rStressMatrix, 
                                                        Vector& rStressVector, const Matrix& LinearElasticMatrix, const Vector& StrainVector )
{    
    noalias(rStressVector) = prod(LinearElasticMatrix, StrainVector);
    noalias(rReturnMappingVariables.TrialIsoStressMatrix) = MathUtils<double>::StressVectorToTensor(rStressVector);
    
    Matrix Aux1 = Matrix();
    Matrix Aux2 = Matrix();
    
    mpFlowRule->CalculateReturnMapping( rReturnMappingVariables, Aux1, rStressMatrix, Aux2 );
    
    noalias(rStressVector) = MathUtils<double>::StressTensorToVector( rStressMatrix, StrainVector.size() );
}

} // Namespace Kratos