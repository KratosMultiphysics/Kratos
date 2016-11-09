//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/restore_simo_ju_nonlocal_3D_law.hpp"

namespace Kratos
{

//Default Constructor
RestoreSimoJuNonlocal3DLaw::RestoreSimoJuNonlocal3DLaw() : RestoreSimoJu3DLaw() 
{
  //mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
  //mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
  mpFlowRule = FlowRule::Pointer( new RestoreNonlocalDamageFlowRule(mpYieldCriterion) );
}

//----------------------------------------------------------------------------------------

//Second Constructor
RestoreSimoJuNonlocal3DLaw::RestoreSimoJuNonlocal3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : RestoreSimoJu3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
RestoreSimoJuNonlocal3DLaw::RestoreSimoJuNonlocal3DLaw(const RestoreSimoJuNonlocal3DLaw& rOther) : RestoreSimoJu3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
RestoreSimoJuNonlocal3DLaw::~RestoreSimoJuNonlocal3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int RestoreSimoJuNonlocal3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo)
{
    int ierr = RestoreSimoJu3DLaw::Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if ( LOCAL_EQUIVALENT_STRAIN.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"LOCAL_EQUIVALENT_STRAIN has Key zero! (check if the application is correctly registered", "" )
    if ( NONLOCAL_EQUIVALENT_STRAIN.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"NONLOCAL_EQUIVALENT_STRAIN has Key zero! (check if the application is correctly registered", "" )
    return ierr;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer RestoreSimoJuNonlocal3DLaw::Clone() const
{
    RestoreSimoJuNonlocal3DLaw::Pointer p_clone(new RestoreSimoJuNonlocal3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------

void RestoreSimoJuNonlocal3DLaw::InitializeMaterial( const Properties& rMaterialProperties,
						   const GeometryType& rElementGeometry,
						   const Vector& rShapeFunctionsValues )
{
    HyperElasticPlastic3DLaw::InitializeMaterial(rMaterialProperties,rElementGeometry,rShapeFunctionsValues);
    
    mNonlocalEquivalentStrain = 0.0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void RestoreSimoJuNonlocal3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
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
    const unsigned int Dim = this->WorkingSpaceDimension();
    Matrix AuxMatrix(Dim,Dim);
    noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(rStrainVector);
    ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
    noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
    ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);
    // CharacteristicSize (for nonlocal damage it must be 1.0)
    ReturnMappingVariables.CharacteristicSize = 1.0;
    
    // LinearElasticMatrix
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    const unsigned int VoigtSize = rStrainVector.size();
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);
    
    if(Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY)) // LOCAL QUANTITIES
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

void RestoreSimoJuNonlocal3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
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
    const unsigned int Dim = this->WorkingSpaceDimension();
    Matrix AuxMatrix(Dim,Dim);
    noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(rStrainVector);
    ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
    noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
    ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);
    // CharacteristicSize
    ReturnMappingVariables.CharacteristicSize = 1.0;
    ReturnMappingVariables.NormIsochoricStress = mNonlocalEquivalentStrain;

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

double& RestoreSimoJuNonlocal3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if (rThisVariable==LOCAL_EQUIVALENT_STRAIN)
    {
        const FlowRule::ThermalVariables& ThermalVariables = mpFlowRule->GetThermalVariables();
        rValue=ThermalVariables.DeltaPlasticDissipation;
    }
    else
    {
        RestoreSimoJu3DLaw::GetValue(rThisVariable,rValue);
    }

    return( rValue );
}

//----------------------------------------------------------------------------------------

void RestoreSimoJuNonlocal3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                        const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == NONLOCAL_EQUIVALENT_STRAIN)
    {
        mNonlocalEquivalentStrain = rValue;
    }
    else
    {
        RestoreSimoJu3DLaw::SetValue(rThisVariable,rValue,rCurrentProcessInfo);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void RestoreSimoJuNonlocal3DLaw::CalculateLocalReturnMapping( FlowRule::RadialReturnVariables& rReturnMappingVariables, Matrix& rStressMatrix, 
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