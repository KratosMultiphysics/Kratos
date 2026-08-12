//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/thermal_local_damage_3D_law.hpp"
#include "custom_constitutive/thermal_output_utilities.hpp"

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

int ThermalLocalDamage3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = LocalDamage3DLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    return ierr;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer ThermalLocalDamage3DLaw::Clone() const
{
    ThermalLocalDamage3DLaw::Pointer p_clone(new ThermalLocalDamage3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Lifecycle
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bool ThermalLocalDamage3DLaw::RequiresInitializeMaterialResponse()
{
    // The damage state is initialized through InitializeMaterial (the flow-rule
    // history/threshold) and subsequently managed through the finalization
    // (commit/restore with IS_CONVERGED). No local-damage state is initialized
    // through the InitializeMaterialResponse hooks, so the (PK2/Kirchhoff)
    // initialization hooks are not required.
    return false;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLocalDamage3DLaw::ReinitializeMaterialProperties(const Properties& rMaterialProperties)
{
    // After a serialization/restart the transient Properties on the hardening
    // law are not restored by the Serializer; re-establish them so that the
    // flow/yield/hardening hierarchy can respond. The committed damage/history
    // state is untouched. This is idempotent and called automatically before
    // every entry into the damage hierarchy.
    if (mpHardeningLaw)
    {
        mpHardeningLaw->SetProperties(rMaterialProperties);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Material response entry points
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLocalDamage3DLaw::CalculateMaterialResponsePK2 (Parameters & rValues)
{
    this->CalculateThermalDamageResponse(rValues);
}

void ThermalLocalDamage3DLaw::CalculateMaterialResponseKirchhoff (Parameters & rValues)
{
    this->CalculateThermalDamageResponse(rValues);
}

void ThermalLocalDamage3DLaw::CalculateMaterialResponseCauchy (Parameters & rValues)
{
    this->CalculateThermalDamageResponse(rValues);
}

void ThermalLocalDamage3DLaw::FinalizeMaterialResponsePK2 (Parameters & rValues)
{
    this->FinalizeThermalDamageResponse(rValues);
}

void ThermalLocalDamage3DLaw::FinalizeMaterialResponseKirchhoff (Parameters & rValues)
{
    this->FinalizeThermalDamageResponse(rValues);
}

void ThermalLocalDamage3DLaw::FinalizeMaterialResponseCauchy (Parameters & rValues)
{
    this->FinalizeThermalDamageResponse(rValues);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Parameter validation
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLocalDamage3DLaw::CheckThermalDamageParameters(Parameters& rValues) const
{
    KRATOS_ERROR_IF_NOT(rValues.IsSetStrainVector()) << "ThermalLocalDamage3DLaw: StrainVector NOT SET" << std::endl;
    KRATOS_ERROR_IF_NOT(rValues.IsSetShapeFunctionsValues()) << "ThermalLocalDamage3DLaw: ShapeFunctionsValues NOT SET" << std::endl;
    KRATOS_ERROR_IF_NOT(rValues.IsSetElementGeometry()) << "ThermalLocalDamage3DLaw: ElementGeometry NOT SET" << std::endl;
    KRATOS_ERROR_IF_NOT(rValues.IsSetMaterialProperties()) << "ThermalLocalDamage3DLaw: MaterialProperties NOT SET" << std::endl;
    KRATOS_ERROR_IF_NOT(rValues.IsSetProcessInfo()) << "ThermalLocalDamage3DLaw: ProcessInfo NOT SET" << std::endl;

    const Flags& Options = rValues.GetOptions();
    if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
    {
        KRATOS_ERROR_IF_NOT(rValues.IsSetStressVector()) << "ThermalLocalDamage3DLaw: StressVector NOT SET" << std::endl;
    }
    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {
        KRATOS_ERROR_IF_NOT(rValues.IsSetConstitutiveMatrix()) << "ThermalLocalDamage3DLaw: ConstitutiveMatrix NOT SET" << std::endl;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Common response
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLocalDamage3DLaw::CalculateThermalDamageResponse(Parameters& rValues)
{
    // Check only the parameters genuinely consumed by this law.
    this->CheckThermalDamageParameters(rValues);

    // Automatic transient rebinding: the hardening-law Properties are not
    // serialized, so they must be re-established from the CURRENT parameters
    // before entering the damage hierarchy (required after restart).
    this->ReinitializeMaterialProperties(rValues.GetMaterialProperties());

    // Get values for the constitutive law
    Flags& Options = rValues.GetOptions();
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    const Vector& rTotalStrainVector = rValues.GetStrainVector();

    // Initialize main variables //

    // LinearElasticMatrix
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    const unsigned int VoigtSize = rTotalStrainVector.size();
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);

    // MaterialResponseVariables (Thermal variables)
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.SetShapeFunctionsValues(rValues.GetShapeFunctionsValues());
    ElasticVariables.SetElementGeometry(rValues.GetElementGeometry());
    ElasticVariables.LameMu = 1.0+PoissonCoefficient;
    ElasticVariables.ThermalExpansionCoefficient = MaterialProperties[THERMAL_EXPANSION];
    /* Calculate Nodal Reference Temperature */
    double NodalReferenceTemperature;
    this->CalculateNodalReferenceTemperature(ElasticVariables,NodalReferenceTemperature);

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

    // Mechanical strain is computed in a LOCAL vector so that the
    // element-provided total strain vector is never destructively modified.
    Vector MechanicalStrainVector(VoigtSize);
    noalias(MechanicalStrainVector) = rTotalStrainVector;

    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){ //TOTAL STRESS

      // Thermal strain
      Vector ThermalStrainVector(VoigtSize);
      this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables,NodalReferenceTemperature);
      // Mechanical strain
      noalias(MechanicalStrainVector) -= ThermalStrainVector;
      noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(MechanicalStrainVector);
      noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;

      if(Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS)){

	// COMPUTE_CONSTITUTIVE_TENSOR
	Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
	Vector EffectiveStressVector(VoigtSize);

	this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,EffectiveStressVector,LinearElasticMatrix,MechanicalStrainVector,rValues);

	this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
      }
      else{

	// COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
	Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
	Vector& rStressVector = rValues.GetStressVector();

	this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,MechanicalStrainVector,rValues);

	this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
      }
    }
    else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS)){ //TOTAL STRESS

      if(Options.Is(ConstitutiveLaw::MECHANICAL_RESPONSE_ONLY)){ //This should be COMPUTE_MECHANICAL_STRESS

        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();

        // Total Strain (used as mechanical strain)
        noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(MechanicalStrainVector);
        noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;

        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,MechanicalStrainVector,rValues);
      }
      else if(Options.Is(ConstitutiveLaw::THERMAL_RESPONSE_ONLY)){ //This should be COMPUTE_THERMAL_STRESS

	// COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();

        // Thermal strain
        Vector ThermalStrainVector(VoigtSize);
        this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables,NodalReferenceTemperature);
        noalias(MechanicalStrainVector) = ThermalStrainVector;
        noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(MechanicalStrainVector);
        noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;

        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,MechanicalStrainVector,rValues);
      }
      else{

	// COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();

        // Thermal strain
        Vector ThermalStrainVector(VoigtSize);
        this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables,NodalReferenceTemperature);
        // Mechanical strain
        noalias(MechanicalStrainVector) -= ThermalStrainVector;
        noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(MechanicalStrainVector);
        noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;

        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,MechanicalStrainVector,rValues);

      }

    }
    else if(Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)){ //This should be COMPUTE_THERMAL_STRAIN

      // VOLUMETRIC_TENSOR_ONLY
      if(Options.Is(ConstitutiveLaw::THERMAL_RESPONSE_ONLY)){

	// Thermal strain (the strain vector is the output channel of this query)
        Vector& rStrainVector = rValues.GetStrainVector();
        this->CalculateThermalStrain(rStrainVector,ElasticVariables,NodalReferenceTemperature);
      }
      //other strain: to implement

    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Common finalization
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLocalDamage3DLaw::FinalizeThermalDamageResponse(Parameters& rValues)
{
    // Check only the parameters genuinely consumed by this law.
    this->CheckThermalDamageParameters(rValues);

    // Automatic transient rebinding (idempotent; required after restart).
    this->ReinitializeMaterialProperties(rValues.GetMaterialProperties());

    // Get values for the constitutive law
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    const Vector& rTotalStrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = rTotalStrainVector.size();
    Vector EffectiveStressVector(VoigtSize);

    // Initialize main variables //

    // LinearElasticMatrix
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);

    // MaterialResponseVariables (Thermal variables)
    MaterialResponseVariables ElasticVariables;
	ElasticVariables.SetShapeFunctionsValues(rValues.GetShapeFunctionsValues());
	ElasticVariables.SetElementGeometry(rValues.GetElementGeometry());
    ElasticVariables.LameMu = 1.0+PoissonCoefficient;
    ElasticVariables.ThermalExpansionCoefficient = MaterialProperties[THERMAL_EXPANSION];
    /* Calculate Nodal Reference Temperature */
    double NodalReferenceTemperature;
    this->CalculateNodalReferenceTemperature(ElasticVariables,NodalReferenceTemperature);

    // Compute Thermal strain
    Vector ThermalStrainVector(VoigtSize);
    this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables,NodalReferenceTemperature);
    // Compute Mechanical strain in a LOCAL vector
    Vector MechanicalStrainVector(VoigtSize);
    noalias(MechanicalStrainVector) = rTotalStrainVector;
    noalias(MechanicalStrainVector) -= ThermalStrainVector;

    // ReturnMappingVariables
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize();
    // Strain and Stress matrices
    const unsigned int Dim = this->WorkingSpaceDimension();
    ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
    noalias(ReturnMappingVariables.StrainMatrix) = MathUtils<double>::StrainVectorToTensor(MechanicalStrainVector);
    ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);
    // CharacteristicSize
    double CharacteristicSize = 1.0;
    this->CalculateCharacteristicSize(CharacteristicSize,rValues.GetElementGeometry());
    ReturnMappingVariables.CharacteristicSize = CharacteristicSize;

    if(rValues.GetProcessInfo()[IS_CONVERGED]==true) //Convergence is achieved. Save equilibrium state variable
    {
        ReturnMappingVariables.Options.Set(FlowRule::RETURN_MAPPING_COMPUTED,false); // Restore sate variable = false

        this->UpdateInternalStateVariables(ReturnMappingVariables,EffectiveStressVector,LinearElasticMatrix,MechanicalStrainVector,rValues);
    }
    else // No convergence is achieved. Restore state variable to equilibrium
    {
        ReturnMappingVariables.Options.Set(FlowRule::RETURN_MAPPING_COMPUTED,true); // Restore sate variable = true

        this->UpdateInternalStateVariables(ReturnMappingVariables,EffectiveStressVector,LinearElasticMatrix,MechanicalStrainVector,rValues);
    }

    if(rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRESS)) //TOTAL STRESS
    {
        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();

        this->UpdateStressVector(rStressVector,ReturnMappingVariables,EffectiveStressVector);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double&  ThermalLocalDamage3DLaw::CalculateNodalReferenceTemperature (const MaterialResponseVariables & rElasticVariables, double & rNodalReferenceTemperature)
{
    KRATOS_TRY

    const GeometryType& DomainGeometry = rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues = rElasticVariables.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();

    rNodalReferenceTemperature = 0.0;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
      rNodalReferenceTemperature += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE);
    }

    return rNodalReferenceTemperature;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLocalDamage3DLaw::CalculateThermalStrain(Vector& rThermalStrainVector, const MaterialResponseVariables& ElasticVariables, double & rNodalReferenceTemperature)
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
    double DeltaTemperature = Temperature - rNodalReferenceTemperature;

    //Thermal strain vector
    for(unsigned int i = 0; i < 6; i++)
        rThermalStrainVector[i] *= ElasticVariables.ThermalExpansionCoefficient * DeltaTemperature;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Vector& ThermalLocalDamage3DLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    KRATOS_TRY

    if (rThisVariable == THERMAL_STRAIN_VECTOR ||
        rThisVariable == THERMAL_STRESS_VECTOR ||
        rThisVariable == MECHANICAL_STRESS_VECTOR) {

        // These outputs are computed from the current state carried by the
        // Parameters (element-provided total strain, shape functions, geometry),
        // so they reach the law through CalculateOnConstitutiveLaw. Has() is
        // intentionally left false so that the parameter-dependent CalculateValue
        // path is used instead of GetValue.
        const Properties& r_material_properties = rParameterValues.GetMaterialProperties();
        const double& r_young_modulus = r_material_properties[YOUNG_MODULUS];
        const double& r_poisson_ratio = r_material_properties[POISSON_RATIO];
        const Vector& r_total_strain = rParameterValues.GetStrainVector();
        const std::size_t voigt_size = r_total_strain.size();

        // Automatic transient rebinding before the damage-factor evaluation
        // (required after restart; the hardening-law Properties are transient).
        this->ReinitializeMaterialProperties(r_material_properties);

        // Constitutive matrix (dimensional specialization via virtual dispatch).
        Matrix constitutive_matrix(voigt_size, voigt_size);
        noalias(constitutive_matrix) = ZeroMatrix(voigt_size, voigt_size);
        this->CalculateLinearElasticMatrix(constitutive_matrix, r_young_modulus, r_poisson_ratio);

        // Thermal strain (dimensional specialization via virtual dispatch).
        MaterialResponseVariables elastic_variables;
        elastic_variables.SetShapeFunctionsValues(rParameterValues.GetShapeFunctionsValues());
        elastic_variables.SetElementGeometry(rParameterValues.GetElementGeometry());
        elastic_variables.LameMu = 1.0 + r_poisson_ratio;
        elastic_variables.ThermalExpansionCoefficient = r_material_properties[THERMAL_EXPANSION];

        double reference_temperature;
        this->CalculateNodalReferenceTemperature(elastic_variables, reference_temperature);

        Vector thermal_strain_vector(voigt_size);
        this->CalculateThermalStrain(thermal_strain_vector, elastic_variables, reference_temperature);

        // CURRENT damage factor (1 - d) from the committed maximum equivalent
        // strain, exactly as the trial response computes it for the current
        // mechanical strain. The evaluation is read-only and never commits.
        double characteristic_size = 1.0;
        this->CalculateCharacteristicSize(characteristic_size, rParameterValues.GetElementGeometry());
        const double damage_factor = ThermalOutputUtilities::CalculateCurrentDamageFactor(
            *mpFlowRule, *mpYieldCriterion, characteristic_size);

        // Assemble the three outputs from the SAME damage factor.
        Vector thermal_strain_out, thermal_stress_out, mechanical_stress_out;
        ThermalOutputUtilities::AssembleOutputs(thermal_strain_out, thermal_stress_out, mechanical_stress_out,
                                                r_total_strain, thermal_strain_vector, constitutive_matrix, damage_factor);

        if (rThisVariable == THERMAL_STRAIN_VECTOR) {
            if (rValue.size() != voigt_size)
                rValue.resize(voigt_size, false);
            noalias(rValue) = thermal_strain_out;
            return rValue;
        }
        if (rThisVariable == THERMAL_STRESS_VECTOR) {
            if (rValue.size() != voigt_size)
                rValue.resize(voigt_size, false);
            noalias(rValue) = thermal_stress_out;
            return rValue;
        }
        // MECHANICAL_STRESS_VECTOR
        if (rValue.size() != voigt_size)
            rValue.resize(voigt_size, false);
        noalias(rValue) = mechanical_stress_out;
        return rValue;
    }

    // Not one of the specialized outputs: keep the base behaviour.
    return rValue;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Matrix& ThermalLocalDamage3DLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue)
{
    KRATOS_TRY

    if (rThisVariable == THERMAL_STRAIN_TENSOR) {
        // Reuse the validated vector output as the single source of truth.
        Vector strain_vector = ZeroVector(this->GetStrainSize());
        this->CalculateValue(rParameterValues, THERMAL_STRAIN_VECTOR, strain_vector);
        ThermalOutputUtilities::AssignStrainTensor(rValue, strain_vector);
        return rValue;
    }

    if (rThisVariable == THERMAL_STRESS_TENSOR) {
        Vector stress_vector = ZeroVector(this->GetStrainSize());
        this->CalculateValue(rParameterValues, THERMAL_STRESS_VECTOR, stress_vector);
        ThermalOutputUtilities::AssignStressTensor(rValue, stress_vector);
        return rValue;
    }

    if (rThisVariable == MECHANICAL_STRESS_TENSOR) {
        Vector stress_vector = ZeroVector(this->GetStrainSize());
        this->CalculateValue(rParameterValues, MECHANICAL_STRESS_VECTOR, stress_vector);
        ThermalOutputUtilities::AssignStressTensor(rValue, stress_vector);
        return rValue;
    }

    // Not one of the specialized outputs: keep the base behaviour.
    return rValue;

    KRATOS_CATCH( "" )
}

} // Namespace Kratos
