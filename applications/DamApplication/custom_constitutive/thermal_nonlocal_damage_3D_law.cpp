//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/thermal_nonlocal_damage_3D_law.hpp"
#include "custom_constitutive/thermal_output_utilities.hpp"

namespace Kratos
{

//Default Constructor
ThermalNonlocalDamage3DLaw::ThermalNonlocalDamage3DLaw() : NonlocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------

//Second Constructor
ThermalNonlocalDamage3DLaw::ThermalNonlocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : NonlocalDamage3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalNonlocalDamage3DLaw::ThermalNonlocalDamage3DLaw(const ThermalNonlocalDamage3DLaw& rOther)
    : NonlocalDamage3DLaw(rOther)
{
    // The base NonlocalDamage3DLaw copy constructor does NOT preserve
    // mNonlocalEquivalentStrain; it is preserved here in the Dam derived copy
    // path (the member is protected and accessible in this class).
    mNonlocalEquivalentStrain = rOther.mNonlocalEquivalentStrain;
}

//----------------------------------------------------------------------------------------

//Destructor
ThermalNonlocalDamage3DLaw::~ThermalNonlocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int ThermalNonlocalDamage3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = NonlocalDamage3DLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    return ierr;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer ThermalNonlocalDamage3DLaw::Clone() const
{
    ThermalNonlocalDamage3DLaw::Pointer p_clone(new ThermalNonlocalDamage3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Lifecycle
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bool ThermalNonlocalDamage3DLaw::RequiresInitializeMaterialResponse()
{
    // The damage/history state is initialized through InitializeMaterial (flow
    // rule) and managed through the finalization. No local-damage state is
    // initialized through the InitializeMaterialResponse hooks.
    return false;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalNonlocalDamage3DLaw::ReinitializeMaterialProperties(const Properties& rMaterialProperties)
{
    // After a serialization/restart the transient Properties on the hardening
    // law are not restored by the Serializer; re-establish them so that a
    // subsequent material response works. The committed damage/history state is
    // untouched.
    if (mpHardeningLaw)
    {
        mpHardeningLaw->SetProperties(rMaterialProperties);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Parameter-aware scalar outputs
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double& ThermalNonlocalDamage3DLaw::CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue)
{
    if (rThisVariable == LOCAL_EQUIVALENT_STRAIN)
    {
        // Recompute the LOCAL driving quantity from the CURRENT kinematics
        // supplied by the element (the generic SMA CalculateOnConstitutiveLaw
        // path) and store it in the flow-rule state.
        this->CalculateLocalEquivalentStrain(rParameterValues);
        return this->GetValue(LOCAL_EQUIVALENT_STRAIN, rValue);
    }

    return HyperElasticPlastic3DLaw::CalculateValue(rParameterValues, rThisVariable, rValue);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Vector& ThermalNonlocalDamage3DLaw::CalculateValue(
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
        // strain driven by mNonlocalEquivalentStrain, exactly as the current
        // NONLOCAL response computes it. The evaluation is read-only: no LOCAL
        // recompute, no LOCAL/NONLOCAL change and no commit.
        const double damage_factor = ThermalOutputUtilities::CalculateCurrentDamageFactor(
            *mpFlowRule, *mpYieldCriterion, 1.0 /* nonlocal characteristic size */);

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

Matrix& ThermalNonlocalDamage3DLaw::CalculateValue(
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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Material response entry points
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalNonlocalDamage3DLaw::CalculateMaterialResponsePK2 (Parameters & rValues)
{
    this->CalculateThermalNonlocalDamageResponse(rValues);
}

void ThermalNonlocalDamage3DLaw::CalculateMaterialResponseKirchhoff (Parameters & rValues)
{
    this->CalculateThermalNonlocalDamageResponse(rValues);
}

void ThermalNonlocalDamage3DLaw::CalculateMaterialResponseCauchy (Parameters & rValues)
{
    this->CalculateThermalNonlocalDamageResponse(rValues);
}

void ThermalNonlocalDamage3DLaw::FinalizeMaterialResponsePK2 (Parameters & rValues)
{
    this->FinalizeThermalNonlocalDamageResponse(rValues);
}

void ThermalNonlocalDamage3DLaw::FinalizeMaterialResponseKirchhoff (Parameters & rValues)
{
    this->FinalizeThermalNonlocalDamageResponse(rValues);
}

void ThermalNonlocalDamage3DLaw::FinalizeMaterialResponseCauchy (Parameters & rValues)
{
    this->FinalizeThermalNonlocalDamageResponse(rValues);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Parameter validation
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalNonlocalDamage3DLaw::CheckThermalNonlocalDamageParameters(Parameters& rValues) const
{
    KRATOS_ERROR_IF_NOT(rValues.IsSetStrainVector()) << "ThermalNonlocalDamage3DLaw: StrainVector NOT SET" << std::endl;
    KRATOS_ERROR_IF_NOT(rValues.IsSetShapeFunctionsValues()) << "ThermalNonlocalDamage3DLaw: ShapeFunctionsValues NOT SET" << std::endl;
    KRATOS_ERROR_IF_NOT(rValues.IsSetElementGeometry()) << "ThermalNonlocalDamage3DLaw: ElementGeometry NOT SET" << std::endl;
    KRATOS_ERROR_IF_NOT(rValues.IsSetMaterialProperties()) << "ThermalNonlocalDamage3DLaw: MaterialProperties NOT SET" << std::endl;
    KRATOS_ERROR_IF_NOT(rValues.IsSetProcessInfo()) << "ThermalNonlocalDamage3DLaw: ProcessInfo NOT SET" << std::endl;

    const Flags& Options = rValues.GetOptions();
    if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
    {
        KRATOS_ERROR_IF_NOT(rValues.IsSetStressVector()) << "ThermalNonlocalDamage3DLaw: StressVector NOT SET" << std::endl;
    }
    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {
        KRATOS_ERROR_IF_NOT(rValues.IsSetConstitutiveMatrix()) << "ThermalNonlocalDamage3DLaw: ConstitutiveMatrix NOT SET" << std::endl;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Common LOCAL equivalent-strain calculation
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalNonlocalDamage3DLaw::CalculateLocalEquivalentStrain(
    Parameters& rValues,
    Vector* pStressVector,
    Matrix* pConstitutiveMatrix)
{
    // Validate only what is genuinely consumed (no shape-function derivatives,
    // no F/detF requirement).
    this->CheckThermalNonlocalDamageParameters(rValues);

    // Automatic transient rebinding (idempotent; required after restart).
    this->ReinitializeMaterialProperties(rValues.GetMaterialProperties());

    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    const Vector& rTotalStrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = rTotalStrainVector.size();

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

    // Thermal strain
    Vector ThermalStrainVector(VoigtSize);
    this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables,NodalReferenceTemperature);
    // Mechanical strain (local copy, the supplied total strain is not mutated)
    Vector MechanicalStrainVector(VoigtSize);
    noalias(MechanicalStrainVector) = rTotalStrainVector;
    noalias(MechanicalStrainVector) -= ThermalStrainVector;

    // ReturnMappingVariables
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize();
    const unsigned int Dim = this->WorkingSpaceDimension();
    Matrix AuxMatrix(Dim,Dim);
    noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(MechanicalStrainVector);
    ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
    noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
    ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);
    // CharacteristicSize (for nonlocal damage it must be 1.0)
    ReturnMappingVariables.CharacteristicSize = 1.0;

    Vector EffectiveStressVector(VoigtSize);
    this->CalculateLocalReturnMapping(ReturnMappingVariables,AuxMatrix,EffectiveStressVector,LinearElasticMatrix,MechanicalStrainVector,rValues);

    if (pConstitutiveMatrix != nullptr)
    {
        this->CalculateConstitutiveTensor(*pConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
    }

    if (pStressVector != nullptr)
    {
        noalias(*pStressVector) = EffectiveStressVector;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Common thermal NONLOCAL response
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalNonlocalDamage3DLaw::CalculateThermalNonlocalDamageResponse(Parameters& rValues)
{
    // Validate only the parameters genuinely consumed by this law.
    this->CheckThermalNonlocalDamageParameters(rValues);

    // Automatic transient rebinding (idempotent; required after restart).
    this->ReinitializeMaterialProperties(rValues.GetMaterialProperties());

    // Get values for the constitutive law
    Flags& Options = rValues.GetOptions();
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    const Vector& rTotalStrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = rTotalStrainVector.size();

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

    // ReturnMappingVariables
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize();
    const unsigned int Dim = this->WorkingSpaceDimension();
    Matrix AuxMatrix(Dim,Dim);
    ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
    ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);
    // CharacteristicSize (for nonlocal damage it must be 1.0)
    ReturnMappingVariables.CharacteristicSize = 1.0;

    if(Options.Is(ConstitutiveLaw::INITIALIZE_MATERIAL_RESPONSE)) // LOCAL QUANTITIES (temporary legacy compatibility)
    {
        // The legacy Dam thermo-mechanical element requests COMPUTE_STRESS (and
        // possibly COMPUTE_CONSTITUTIVE_TENSOR) together with this flag during
        // its nonlinear-iteration hooks; reproduce those outputs.
        Vector* pStress = Options.Is(ConstitutiveLaw::COMPUTE_STRESS) ? &rValues.GetStressVector() : nullptr;
        Matrix* pTangent = Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR) ? &rValues.GetConstitutiveMatrix() : nullptr;

        this->CalculateLocalEquivalentStrain(rValues, pStress, pTangent);
        return;
    }

    // NONLOCAL quantities
    ReturnMappingVariables.NormIsochoricStress = mNonlocalEquivalentStrain;

    // Mechanical strain in a LOCAL vector
    Vector MechanicalStrainVector(VoigtSize);
    noalias(MechanicalStrainVector) = rTotalStrainVector;

    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) // Compute constitutive tensor and total stress
    {
      // Thermal strain
      Vector ThermalStrainVector(VoigtSize);
      this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables,NodalReferenceTemperature);
      // Mechanical strain
      noalias(MechanicalStrainVector) -= ThermalStrainVector;
      noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(MechanicalStrainVector);
      noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;

      if(Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS))
      {
        // COMPUTE_CONSTITUTIVE_TENSOR
        Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
        Vector EffectiveStressVector(VoigtSize);

        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,EffectiveStressVector,LinearElasticMatrix,MechanicalStrainVector,rValues);

        this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
      }
      else
      {
        // COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
        Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
        Vector& rStressVector = rValues.GetStressVector();

        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,MechanicalStrainVector,rValues);

        this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
      }
    }
    else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
    {
      if(Options.Is(ConstitutiveLaw::MECHANICAL_RESPONSE_ONLY))
      {
        // COMPUTE_STRESS: MECHANICAL COMPONENT
        Vector& rStressVector = rValues.GetStressVector();

        // Total Strain
        noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(MechanicalStrainVector);
        noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;

        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,MechanicalStrainVector,rValues);
      }
      else if(Options.Is(ConstitutiveLaw::THERMAL_RESPONSE_ONLY))
      {
        // COMPUTE_STRESS : THEMAL COMPONENT
        Vector& rStressVector = rValues.GetStressVector();

        // Thermal strain
        Vector ThermalStrainVector(VoigtSize);
        this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables,NodalReferenceTemperature);
        noalias(MechanicalStrainVector) = ThermalStrainVector;
        noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(MechanicalStrainVector);
        noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;

        this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,MechanicalStrainVector,rValues);
      }
      else // Compute total stress
      {
        // COMPUTE_STRESS : ALL STRESS COMPONENTS
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
    else if(Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY))
    {
      // VOLUMETRIC_TENSOR_ONLY
      if(Options.Is(ConstitutiveLaw::THERMAL_RESPONSE_ONLY))
      {
        // Thermal strain (the strain vector is the output channel of this query)
        Vector& rStrainVector = rValues.GetStrainVector();
        this->CalculateThermalStrain(rStrainVector,ElasticVariables,NodalReferenceTemperature);
      }
      //other strain: to implement
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Common thermal NONLOCAL finalization
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalNonlocalDamage3DLaw::FinalizeThermalNonlocalDamageResponse(Parameters& rValues)
{
    // Validate only the parameters genuinely consumed by this law.
    this->CheckThermalNonlocalDamageParameters(rValues);

    // Automatic transient rebinding (idempotent; required after restart).
    this->ReinitializeMaterialProperties(rValues.GetMaterialProperties());

    // Get values for the constitutive law
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    const Vector& rTotalStrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = rTotalStrainVector.size();
    Vector EffectiveStressVector(VoigtSize);

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
    ReturnMappingVariables.CharacteristicSize = 1.0;
    ReturnMappingVariables.NormIsochoricStress = mNonlocalEquivalentStrain;

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

    if(rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRESS))
    {
        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();

        this->UpdateStressVector(rStressVector,ReturnMappingVariables,EffectiveStressVector);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double&  ThermalNonlocalDamage3DLaw::CalculateNodalReferenceTemperature (const MaterialResponseVariables & rElasticVariables, double & rNodalReferenceTemperature)
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

void ThermalNonlocalDamage3DLaw::CalculateThermalStrain(Vector& rThermalStrainVector, const MaterialResponseVariables& ElasticVariables, double & rNodalReferenceTemperature)
{
    KRATOS_TRY

    //1.-Temperature from nodes
    const GeometryType& DomainGeometry = ElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues = ElasticVariables.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();

    double Temperature = 0.0;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        Temperature += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(TEMPERATURE);
    }

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

} // Namespace Kratos
