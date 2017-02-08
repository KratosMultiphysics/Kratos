//   
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/thermal_local_damage_3D_law.hpp"

namespace Kratos
{

//Default Constructor
ThermalLocalDamage3DLaw::ThermalLocalDamage3DLaw() : LocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------

//Second Constructor
ThermalLocalDamage3DLaw::ThermalLocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : LocalDamage3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalLocalDamage3DLaw::ThermalLocalDamage3DLaw(const ThermalLocalDamage3DLaw& rOther) : LocalDamage3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
ThermalLocalDamage3DLaw::~ThermalLocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int ThermalLocalDamage3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo)
{
    int ierr = LocalDamage3DLaw::Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if ( THERMAL_EXPANSION.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"THERMAL_EXPANSION Key is 0. Check if all applications were correctly registered.", "" )
    if ( REFERENCE_TEMPERATURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"REFERENCE_TEMPERATURE Key is 0. Check if all applications were correctly registered.", "" )
    if ( TEMPERATURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"TEMPERATURE Key is 0. Check if all applications were correctly registered.", "" )

    return ierr;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer ThermalLocalDamage3DLaw::Clone() const
{
    ThermalLocalDamage3DLaw::Pointer p_clone(new ThermalLocalDamage3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLocalDamage3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{    
    // Check
    rValues.CheckAllParameters();
    
    // Get values for the constitutive law
    Flags& Options = rValues.GetOptions();
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    Vector& rStrainVector = rValues.GetStrainVector();
    
    // Initialize main variables //

    // LinearElasticMatrix
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    const unsigned int VoigtSize = rStrainVector.size();
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);
    
    // MaterialResponseVariables (Thermal variables)
    HyperElastic3DLaw::MaterialResponseVariables ElasticVariables;
	ElasticVariables.SetShapeFunctionsValues(rValues.GetShapeFunctionsValues());
	ElasticVariables.SetElementGeometry(rValues.GetElementGeometry());
    ElasticVariables.LameMu = 1.0+PoissonCoefficient;
    ElasticVariables.ThermalExpansionCoefficient = MaterialProperties[THERMAL_EXPANSION]; 
    ElasticVariables.ReferenceTemperature = rValues.GetProcessInfo()[REFERENCE_TEMPERATURE];
    
    // ReturnMappingVariables
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize();
    // Strain and Stress matrices
    const unsigned int Dim = this->WorkingSpaceDimension();
    Matrix AuxMatrix(Dim,Dim);
    ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
    ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);
    // CharacteristicSize
    double CharacteristicSize = 1.0;
    this->CalculateCharacteristicSize(CharacteristicSize,rValues.GetElementGeometry());
    ReturnMappingVariables.CharacteristicSize = CharacteristicSize;

    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) //TOTAL STRESS
    {
        // Thermal strain
        Vector ThermalStrainVector(VoigtSize);
        this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables);
        // Mechanical strain
        noalias(rStrainVector) -= ThermalStrainVector;
        noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(rStrainVector);
        noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
    
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
    else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS)) //TOTAL STRESS
    {
        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();

        // Thermal strain
        Vector ThermalStrainVector(VoigtSize);
        this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables);
        // Mechanical strain
        noalias(rStrainVector) -= ThermalStrainVector;
        noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(rStrainVector);
        noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
        
        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
    }
    else if(Options.Is(ConstitutiveLaw::TOTAL_TENSOR)) //This should be COMPUTE_MECHANICAL_STRESS
    {
        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();
        
        // Total Strain
        noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(rStrainVector);
        noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
        
        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
    }
    else if(Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)) //This should be COMPUTE_THERMAL_STRESS
    {
        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();

        // Thermal strain
        this->CalculateThermalStrain(rStrainVector,ElasticVariables);
        noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(rStrainVector);
        noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
        
        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
    }
    else if(Options.Is(ConstitutiveLaw::COMPUTE_STRAIN)) //This should be COMPUTE_THERMAL_STRAIN
    {
        // COMPUTE_STRAIN
        
        // Thermal strain
        this->CalculateThermalStrain(rStrainVector,ElasticVariables);
    }
}

//----------------------------------------------------------------------------------------

void ThermalLocalDamage3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{    
    // Check
    rValues.CheckAllParameters();
    
    // Get values for the constitutive law
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    Vector& rStrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = rStrainVector.size();
    Vector EffectiveStressVector(VoigtSize);
    
    // Initialize main variables //

    // LinearElasticMatrix
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);
    
    // MaterialResponseVariables (Thermal variables)
    HyperElastic3DLaw::MaterialResponseVariables ElasticVariables;
	ElasticVariables.SetShapeFunctionsValues(rValues.GetShapeFunctionsValues());
	ElasticVariables.SetElementGeometry(rValues.GetElementGeometry());
    ElasticVariables.LameMu = 1.0+PoissonCoefficient;
    ElasticVariables.ThermalExpansionCoefficient = MaterialProperties[THERMAL_EXPANSION]; 
    ElasticVariables.ReferenceTemperature = rValues.GetProcessInfo()[REFERENCE_TEMPERATURE];
    // Compute Thermal strain
    Vector ThermalStrainVector(VoigtSize);
    this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables);
    // Compute Mechanical strain
    noalias(rStrainVector) -= ThermalStrainVector;
    
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

    if(rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRESS)) //TOTAL STRESS
    {
        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();
        
        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLocalDamage3DLaw::CalculateThermalStrain(Vector& rThermalStrainVector, const MaterialResponseVariables& ElasticVariables)
{
    KRATOS_TRY

    //1.-Temperature from nodes 
    const GeometryType& DomainGeometry = ElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues = ElasticVariables.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();
    
    double Temperature = 0.0;
    
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
      Temperature += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(TEMPERATURE);
    
    //Identity vector
    if(rThermalStrainVector.size()!=6)
        rThermalStrainVector.resize(6,false);
    rThermalStrainVector[0] = 1.0;
    rThermalStrainVector[1] = 1.0;
    rThermalStrainVector[2] = 1.0;
    rThermalStrainVector[3] = 0.0;
    rThermalStrainVector[4] = 0.0;
    rThermalStrainVector[5] = 0.0;

    // Delta T
    double DeltaTemperature = Temperature - ElasticVariables.ReferenceTemperature;

    //Thermal strain vector
    for(unsigned int i = 0; i < 6; i++)
        rThermalStrainVector[i] *= ElasticVariables.ThermalExpansionCoefficient * DeltaTemperature;
    
    KRATOS_CATCH( "" )
}

} // Namespace Kratos