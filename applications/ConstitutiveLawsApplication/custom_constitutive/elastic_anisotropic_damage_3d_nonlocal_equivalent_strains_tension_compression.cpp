// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Athira Vadakkekkara
//  Collaborators:
//
//ToDo: Transform the calculated dHdEpseq_tension and dHdEpseq_compression in to global coordinates.
// Changes from version 1: Strain equivalence, definition of m changed
// Project includes
#include "elastic_anisotropic_damage_3d_nonlocal_equivalent_strains_tension_compression.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "includes/checks.h"
#include "utilities/math_utils.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC()
    : ElasticIsotropic3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC(const ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC &rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::Clone() const
{
    return Kratos::make_shared<ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC>(ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::~ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC()
{
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

int ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS_TENSION));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS_COMPRESSION));
    KRATOS_CHECK(rMaterialProperties.Has(DAMAGE_MODEL_PARAMETER_BETA1_TENSION));
    KRATOS_CHECK(rMaterialProperties.Has(DAMAGE_MODEL_PARAMETER_BETA2_TENSION));
    KRATOS_CHECK(rMaterialProperties.Has(DAMAGE_MODEL_PARAMETER_BETA1_COMPRESSION));
    KRATOS_CHECK(rMaterialProperties.Has(DAMAGE_MODEL_PARAMETER_BETA2_COMPRESSION));
    return 0;
}

//************************************************************************************
//************************************************************************************

bool ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::Has(const Variable<Vector>& rThisVariable)
{
    if(rThisVariable == DAMAGE_VECTOR){
        // explicitly returning "false", so the element calls CalculateValue(...)s
        return true;
    } else if(rThisVariable == STRAIN){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    }
    return false;
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    mMaterialProperties = ZeroVector(6);
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double ft = rMaterialProperties[YIELD_STRESS_TENSION];
    const double nu = rMaterialProperties[POISSON_RATIO];
    mMaterialProperties[0] = (rMaterialProperties.Has(DAMAGE_THRESHOLD_TENSION)==true) ? rMaterialProperties[DAMAGE_THRESHOLD_TENSION] : ft/E;
    mMaterialProperties[1] = (rMaterialProperties.Has(DAMAGE_THRESHOLD_COMPRESSION)==true) ? rMaterialProperties[DAMAGE_THRESHOLD_COMPRESSION] : (1 + nu) * (2.15*ft)/(E);
    mMaterialProperties[2] = rMaterialProperties[DAMAGE_MODEL_PARAMETER_BETA1_TENSION];
    mMaterialProperties[3] = rMaterialProperties[DAMAGE_MODEL_PARAMETER_BETA2_TENSION];
    mMaterialProperties[4] = rMaterialProperties[DAMAGE_MODEL_PARAMETER_BETA1_COMPRESSION];
    mMaterialProperties[5] = rMaterialProperties[DAMAGE_MODEL_PARAMETER_BETA2_COMPRESSION];
}

//************************************************************************************
//************************************************************************************

Vector& ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValues
    )
{
    KRATOS_TRY
    if(rThisVariable == DAMAGE_VECTOR){
        rValues = mInternalDamageVariables.PrincipalDamageVector;
    }else if(rThisVariable == INTERNAL_VARIABLES){
        rValues = mInternalDamageVariables.StrainHistoryParameter;
    }
    return rValues;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValues,
    const ProcessInfo& rProcessInfo
    )
{
    //Why do I need this method?
    KRATOS_TRY
    if(rThisVariable == DAMAGE_VECTOR){
        mInternalDamageVariables.PrincipalDamageVector = rValues;
    }else if(rThisVariable == INTERNAL_VARIABLES){
        mInternalDamageVariables.StrainHistoryParameter = rValues;
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY

    InternalDamageVariables rInternalDamageParameters;
    this->CalculateStressResponse(rParametersValues, rInternalDamageParameters);
    mInternalDamageVariables.DamageMatrix = rInternalDamageParameters.DamageMatrix;
    mInternalDamageVariables.PrincipalDamageVector = rInternalDamageParameters.PrincipalDamageVector;
    mInternalDamageVariables.EquivalentStrainsTC = rInternalDamageParameters.EquivalentStrainsTC;
    mInternalDamageVariables.kappa = rInternalDamageParameters.kappa;
    mInternalDamageVariables.StrainHistoryParameter = rInternalDamageParameters.StrainHistoryParameter;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY

    InternalDamageVariables rInternalDamageVariables;
    this->CalculateStressResponse(rParametersValues, rInternalDamageVariables);
    //check whether internal parameters are set properly
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rParametersValues,
    InternalDamageVariables& rInternalDamageVariables)
{
    KRATOS_TRY
    Flags& r_constitutive_law_options = rParametersValues.GetOptions();
    Vector& r_strain_vector = rParametersValues.GetStrainVector();
    CalculateValue(rParametersValues, STRAIN, r_strain_vector);

    // If we compute the tangent moduli or the stress
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
        {
        Vector& r_stress_vector       = rParametersValues.GetStressVector();
        const Vector& r_strain_vector = rParametersValues.GetStrainVector();
        Matrix& r_constitutive_matrix = rParametersValues.GetConstitutiveMatrix();
        Matrix r_elasticity_matrix;
        r_elasticity_matrix.resize(6, 6, false);
        CalculateElasticMatrix(r_elasticity_matrix, rParametersValues);
        noalias(r_stress_vector) = prod(r_elasticity_matrix, r_strain_vector);
        ModelParameters rModelParameters(Dimension);
        BoundedMatrixVoigtType EffStiffnessMatrix = ZeroMatrix(6, 6);
        BoundedMatrixVoigtType damage_effect_tensor = ZeroMatrix(6, 6);
        BoundedMatrixVoigtType transformed_damage_effect_tensor = ZeroMatrix(6, 6);
        BoundedMatrixVoigtType Inv_M  = ZeroMatrix(6,6);
        BoundedVectorVoigtType H_uNL1 = ZeroVector(6);
        BoundedVectorVoigtType H_uNL2 = ZeroVector(6);
        BoundedVectorVoigtType H_NL1u = ZeroVector(6);
        BoundedVectorVoigtType H_NL2u = ZeroVector(6);
        const double H_NLNL = 0.0;
        double det_M = 0.0;
        Vector nonlocal_equivalent_strains = ZeroVector(2);
        const auto& N = rParametersValues.GetShapeFunctionsValues();
        for (SizeType i = 0; i < N.size(); ++i) {
            nonlocal_equivalent_strains[0] += N[i] * rParametersValues.GetElementGeometry()[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_TENSION, 0);
        }
        for (SizeType i = 0; i < N.size(); ++i) {
            nonlocal_equivalent_strains[1] += N[i] * rParametersValues.GetElementGeometry()[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_COMPRESSION, 0);
        }
        BoundedVectorType principal_strains = ZeroVector(3);
        GetEigenValues(principal_strains, STRAIN, r_strain_vector);
        BoundedVectorType principal_damage_increment_vector;
        BoundedMatrixType damage_increment_tensor;
        BoundedMatrixType total_damage_tensor = mInternalDamageVariables.DamageMatrix;
        BoundedVectorType principal_damage_vector = mInternalDamageVariables.PrincipalDamageVector;
        rModelParameters.kappa = mInternalDamageVariables.kappa;
        rModelParameters.strain_history_parameter = mInternalDamageVariables.StrainHistoryParameter;
        Vector local_equivalent_strains = ZeroVector(2);
        ComputePrincipalDamageIncrement(principal_damage_increment_vector, rParametersValues, rModelParameters, principal_strains, nonlocal_equivalent_strains);
        TransformPrincipalDamageIncrementToGlobal(damage_increment_tensor, principal_damage_increment_vector, r_strain_vector);
        total_damage_tensor += damage_increment_tensor;
        GetTotalPrincipalDamageVector(principal_damage_vector,total_damage_tensor);
        for(SizeType i =0; i < Dimension; ++i){
            if(principal_damage_vector[i] > 0.99){
                principal_damage_vector[i] = 0.99;
            }
        }
        if(principal_damage_increment_vector[0] > 0.0 || principal_damage_increment_vector[1] > 0.0 || principal_damage_increment_vector[2] > 0.0){
            GetDamageEffectTensor(damage_effect_tensor, principal_damage_vector);
            GetTransformedDamageEffectTensor(transformed_damage_effect_tensor, damage_effect_tensor, r_strain_vector);
            MathUtils<double>::InvertMatrix(transformed_damage_effect_tensor, Inv_M, det_M);
            const BoundedMatrixVoigtType temp =prod(Inv_M, r_elasticity_matrix);
            noalias(r_stress_vector)  = prod(temp, r_strain_vector);
            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                Calculate_tangent_Huu(r_elasticity_matrix, rParametersValues);
                Calculate_tangent_HuNL(H_uNL1, H_uNL2, rParametersValues, rModelParameters, principal_damage_vector, principal_strains);
                Calculate_tangent_HNLu(H_NL1u, H_NL2u, rParametersValues, principal_strains);
            }
        }
        AssembleConstitutiveMatrix(r_constitutive_matrix, r_elasticity_matrix, H_NL1u, H_NL2u, H_uNL1, H_uNL2, H_NLNL);
        rInternalDamageVariables.DamageMatrix = total_damage_tensor;
        rInternalDamageVariables.PrincipalDamageVector = principal_damage_vector;
        CalculateLocalEquivalentStrains(local_equivalent_strains, rParametersValues, principal_strains);
        rInternalDamageVariables.EquivalentStrainsTC = local_equivalent_strains;
        rInternalDamageVariables.kappa = rModelParameters.kappa ;
        rInternalDamageVariables.StrainHistoryParameter = rModelParameters.strain_history_parameter;
        }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::Calculate_tangent_Huu(
    Matrix& r_elasticity_matrix,
    ConstitutiveLaw::Parameters& rParametersValues)
{
    // Converged values to be storaged
    const Vector unperturbed_strain_vector_gp = Vector(rParametersValues.GetStrainVector());
    const Vector unperturbed_stress_vector_gp = Vector(rParametersValues.GetStressVector());

    // The number of components
    const SizeType num_components = unperturbed_stress_vector_gp.size();

    // The tangent tensor_Huu
    r_elasticity_matrix.clear();
    Matrix auxiliar_tensor = ZeroMatrix(num_components,num_components);

    // Calculate the perturbation
    double perturbation = PerturbationThreshold;
    if (rParametersValues.GetMaterialProperties().Has(PERTURBATION_SIZE)) {
        perturbation = rParametersValues.GetMaterialProperties()[PERTURBATION_SIZE];
    } else {
        for (IndexType i_component = 0; i_component < num_components; ++i_component) {
            double component_perturbation;
            CalculatePerturbation(unperturbed_strain_vector_gp, i_component, component_perturbation);
            perturbation = std::max(component_perturbation, perturbation);
        }
        // We check that the perturbation has a threshold value of PerturbationThreshold
        if (perturbation < PerturbationThreshold) perturbation = PerturbationThreshold;
    }

    // Loop over components of the strain
    Vector& r_perturbed_strain = rParametersValues.GetStrainVector();
    Vector& r_perturbed_integrated_stress = rParametersValues.GetStressVector();

    for (IndexType i_component = 0; i_component < num_components; ++i_component) {
        // Apply the perturbation (positive)
        PerturbStrainVector(r_perturbed_strain, unperturbed_strain_vector_gp, perturbation, i_component);

        // We continue with the calculations
        IntegratePerturbedStrain(rParametersValues, this, StressMeasure_PK2);

        // Compute stress (plus)
        const Vector strain_plus = r_perturbed_strain;
        const Vector stress_plus = r_perturbed_integrated_stress;

        // Reset the values to the initial ones
        noalias(r_perturbed_strain) = unperturbed_strain_vector_gp;
        noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;

        // Apply the perturbation (negative)
        PerturbStrainVector(r_perturbed_strain, unperturbed_strain_vector_gp, - perturbation, i_component);

        // We continue with the calculations
        IntegratePerturbedStrain(rParametersValues, this,StressMeasure_PK2);

        // Compute stress (minus)
        const Vector strain_minus = r_perturbed_strain;
        const Vector stress_minus = r_perturbed_integrated_stress;

        // Finally we compute the components
        CalculateComponentsToTangentTensorSecondOrder(auxiliar_tensor, strain_plus, strain_minus, stress_plus, stress_minus, i_component);

        // Reset the values to the initial ones
        noalias(r_perturbed_strain) = unperturbed_strain_vector_gp;
        noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;
    }

    noalias(r_elasticity_matrix) = 0.5*(auxiliar_tensor + trans(auxiliar_tensor));
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::Calculate_tangent_HuNL(
    BoundedVectorVoigtType& H_uNL1,
    BoundedVectorVoigtType& H_uNL2,
    ConstitutiveLaw::Parameters& rParametersValues,
    ModelParameters& rModelParameters,
    const Vector& Damage_Vector,
    const BoundedVectorType& Principal_Strains
)
{
    KRATOS_TRY
    Vector nonlocal_equivalent_strains = ZeroVector(2);
    const auto& N = rParametersValues.GetShapeFunctionsValues();
    for (SizeType i = 0; i < N.size(); ++i) {
        nonlocal_equivalent_strains[0] += N[i] * rParametersValues.GetElementGeometry()[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_TENSION, 0);
    }
    for (SizeType i = 0; i < N.size(); ++i) {
        nonlocal_equivalent_strains[1] += N[i] * rParametersValues.GetElementGeometry()[i].FastGetSolutionStepValue(NON_LOCAL_EQUIVALENT_STRAIN_COMPRESSION, 0);
    }
    const Vector unperturbed_stress_vector_gp = Vector(rParametersValues.GetStressVector());
    const Vector unperturbed_nonlocal_equivalent_strains_gp = nonlocal_equivalent_strains;

    // The tangent tensor_Huu
    H_uNL1.clear();
    H_uNL2.clear();

    double perturbation ;
    if (rParametersValues.GetMaterialProperties().Has(PERTURBATION_SIZE)) {
        perturbation = rParametersValues.GetMaterialProperties()[PERTURBATION_SIZE];
    } else {
        perturbation = PerturbationThreshold;
    }

    Vector r_perturbed_nonlocal_equivalent_strains = nonlocal_equivalent_strains;
    Vector r_perturbed_integrated_stress = Vector(rParametersValues.GetStressVector());

    //perturbing and integrating equivalent strain in tension alone
    PerturbNonlocalEquivalentStrains(r_perturbed_nonlocal_equivalent_strains[0], unperturbed_nonlocal_equivalent_strains_gp[0], perturbation);
    IntegratePerturbedNonlocalEquivalentStrains(r_perturbed_integrated_stress, r_perturbed_nonlocal_equivalent_strains, rParametersValues, rModelParameters);

    const Vector nonlocal_equivalent_strains_plus1 = r_perturbed_nonlocal_equivalent_strains;
    const Vector stress_plus1 = r_perturbed_integrated_stress;

    noalias(r_perturbed_nonlocal_equivalent_strains) = unperturbed_nonlocal_equivalent_strains_gp;
    noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;

    PerturbNonlocalEquivalentStrains(r_perturbed_nonlocal_equivalent_strains[0], unperturbed_nonlocal_equivalent_strains_gp[0], -perturbation);
    IntegratePerturbedNonlocalEquivalentStrains(r_perturbed_integrated_stress, r_perturbed_nonlocal_equivalent_strains, rParametersValues, rModelParameters);

    const Vector nonlocal_equivalent_strains_minus1 = r_perturbed_nonlocal_equivalent_strains;
    const Vector stress_minus1 = r_perturbed_integrated_stress;

    const double delta1 = nonlocal_equivalent_strains_plus1[0] - nonlocal_equivalent_strains_minus1[0];

    for (IndexType row = 0; row < VoigtSize; ++row) {
        H_uNL1[row] = (stress_plus1[row] - stress_minus1[row]) / delta1;
    }

    //perturbing and integrating equivalent strain in compresssion alone
    PerturbNonlocalEquivalentStrains(r_perturbed_nonlocal_equivalent_strains[1], unperturbed_nonlocal_equivalent_strains_gp[1], perturbation);
    IntegratePerturbedNonlocalEquivalentStrains(r_perturbed_integrated_stress, r_perturbed_nonlocal_equivalent_strains, rParametersValues, rModelParameters);

    const Vector nonlocal_equivalent_strains_plus2 = r_perturbed_nonlocal_equivalent_strains;
    const Vector stress_plus2 = r_perturbed_integrated_stress;

    noalias(r_perturbed_nonlocal_equivalent_strains) = unperturbed_nonlocal_equivalent_strains_gp;
    noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;

    PerturbNonlocalEquivalentStrains(r_perturbed_nonlocal_equivalent_strains[1], unperturbed_nonlocal_equivalent_strains_gp[1], -perturbation);
    IntegratePerturbedNonlocalEquivalentStrains(r_perturbed_integrated_stress, r_perturbed_nonlocal_equivalent_strains, rParametersValues, rModelParameters);

    const Vector nonlocal_equivalent_strains_minus2 = r_perturbed_nonlocal_equivalent_strains;
    const Vector stress_minus2 = r_perturbed_integrated_stress;
    const double delta2= nonlocal_equivalent_strains_plus2[1] - nonlocal_equivalent_strains_minus2[1];

    for (IndexType row = 0; row < VoigtSize; ++row) {
        H_uNL2[row] = (stress_plus2[row] - stress_minus2[row]) / delta2;
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::IntegratePerturbedNonlocalEquivalentStrains(
    Vector& r_perturbed_integrated_stress,
    const Vector& r_perturbed_nonlocal_equivalent_strains,
    ConstitutiveLaw::Parameters& rParametersValues,
    ModelParameters& rModelParameters)
{
    const Vector r_strain_vector = Vector(rParametersValues.GetStrainVector());
    BoundedVectorType principal_strains = ZeroVector(3);
    GetEigenValues(principal_strains, STRAIN, r_strain_vector);
    BoundedVectorType principal_damage_increment_vector;
    BoundedMatrixType damage_increment_tensor;
    BoundedMatrixType total_damage_tensor = mInternalDamageVariables.DamageMatrix;
    BoundedVectorType principal_damage_vector = mInternalDamageVariables.PrincipalDamageVector;
    BoundedMatrixVoigtType damage_effect_tensor = ZeroMatrix(6, 6);
    BoundedMatrixVoigtType transformed_damage_effect_tensor = ZeroMatrix(6, 6);
    BoundedMatrixVoigtType Inv_M  = ZeroMatrix(6,6);
    double det_M;
    Matrix r_elasticity_matrix;
    r_elasticity_matrix.resize(6, 6, false);
    CalculateElasticMatrix(r_elasticity_matrix, rParametersValues);
    ComputePrincipalDamageIncrement(principal_damage_increment_vector, rParametersValues, rModelParameters, principal_strains,r_perturbed_nonlocal_equivalent_strains);
    TransformPrincipalDamageIncrementToGlobal(damage_increment_tensor, principal_damage_increment_vector, r_strain_vector);
    total_damage_tensor += damage_increment_tensor;
    GetTotalPrincipalDamageVector(principal_damage_vector,total_damage_tensor);
    GetDamageEffectTensor(damage_effect_tensor, principal_damage_vector);
    GetTransformedDamageEffectTensor(transformed_damage_effect_tensor, damage_effect_tensor, r_strain_vector);
    MathUtils<double>::InvertMatrix(transformed_damage_effect_tensor, Inv_M, det_M);
    const BoundedMatrixVoigtType temp = prod(Inv_M, r_elasticity_matrix);
    r_perturbed_integrated_stress  = prod(temp, r_strain_vector);
}

//************************************************************************************
//************************************************************************************
void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::PerturbNonlocalEquivalentStrains(
    double& rPerturbedNonlocalEquivalentStrainComponent,
    const double& rUnperturbedNonlocalEquivalentStrainComponent,
    const double Perturbation
)
{
    rPerturbedNonlocalEquivalentStrainComponent = rUnperturbedNonlocalEquivalentStrainComponent;
    rPerturbedNonlocalEquivalentStrainComponent += Perturbation;

}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CalculatePerturbation(
    const Vector& rStrainVector,
    const IndexType Component,
    double& rPerturbation)
{
    double perturbation_1, perturbation_2;
    if (std::abs(rStrainVector[Component]) > tolerance) {
        perturbation_1 = PerturbationCoefficient1 * rStrainVector[Component];
    } else {
        double min_strain_component;
        GetMinAbsValue(rStrainVector, min_strain_component);
        perturbation_1 = PerturbationCoefficient1 * min_strain_component;
    }
    double max_strain_component;
    GetMaxAbsValue(rStrainVector, max_strain_component);
    perturbation_2 = PerturbationCoefficient2 * max_strain_component;
    rPerturbation = std::max(perturbation_1, perturbation_2);
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetMinAbsValue(
    const Vector& rArrayValues,
    double& rMinValue)
{
    const SizeType dimension = rArrayValues.size();

    IndexType counter = 0;
    double aux = std::numeric_limits<double>::max();
    for (IndexType i = 0; i < dimension; ++i) {
        if (std::abs(rArrayValues[i]) < aux) {
            aux = std::abs(rArrayValues[i]);
            ++counter;
        }
    }

    rMinValue = aux;
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetMaxAbsValue(
    const Vector& rArrayValues,
    double& rMaxValue)
{
    const SizeType dimension = rArrayValues.size();

    IndexType counter = 0;
    double aux = 0.0;
    for (IndexType i = 0; i < dimension; ++i) {
        if (std::abs(rArrayValues[i]) > aux) {
            aux = std::abs(rArrayValues[i]);
            ++counter;
        }
    }

    rMaxValue = aux;
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::PerturbStrainVector(
    Vector& rPerturbedStrainVector,
    const Vector& rStrainVectorGP,
    const double Perturbation,
    const IndexType Component
    )
{
    noalias(rPerturbedStrainVector) = rStrainVectorGP;
    rPerturbedStrainVector[Component] += Perturbation;
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::IntegratePerturbedStrain(
    ConstitutiveLaw::Parameters& rValues,
    ConstitutiveLaw* pConstitutiveLaw,
    const ConstitutiveLaw::StressMeasure& rStressMeasure)
{
    Flags& cl_options = rValues.GetOptions();

    // In order to avoid recursivity...
    const bool flag_back_up_1 = cl_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_back_up_2 = cl_options.Is(ConstitutiveLaw::COMPUTE_STRESS);

    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    pConstitutiveLaw->CalculateMaterialResponse(rValues, rStressMeasure);

    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_back_up_1);
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_back_up_2);
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CalculateComponentsToTangentTensorSecondOrder(
    Matrix& rTangentTensor,
    const Vector& rVectorStrainPlus,
    const Vector& rVectorStrainMinus,
    const Vector& rStressPlus,
    const Vector& rStressMinus,
    const IndexType Component
    )
{
    const double perturbation = (rVectorStrainPlus[Component] - rVectorStrainMinus[Component]);
    const SizeType voigt_size = rStressPlus.size();
    for (IndexType row = 0; row < voigt_size; ++row) {
        rTangentTensor(row, Component) = (rStressPlus[row] - rStressMinus[row]) / perturbation;
    }
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetEigenValues(
    BoundedVectorType& Pri_Values,
    const Variable<Vector>& rThisVariable,
    const Vector& VectorForm)
{
    KRATOS_TRY
    BoundedMatrixType MatrixForm = ZeroMatrix(3,3);
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    VectorToTensor(MatrixForm, VectorForm, rThisVariable);
    //prinicpal values, max and min
    MathUtils<double>::GaussSeidelEigenSystem(MatrixForm, EigenVectors, EigenValues);
    Pri_Values[0] = EigenValues(0,0);
    Pri_Values[1] = EigenValues(1,1);
    Pri_Values[2] = EigenValues(2,2);
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::ComputePrincipalDamageIncrement(
    BoundedVectorType& PrincipalDamageIncrement,
    ConstitutiveLaw::Parameters& rParametersValues,
    ModelParameters& rModelParameters,
    const BoundedVectorType& PrincipalStrains,
    const Vector& nonlocal_equivalent_strains
    )
{
    KRATOS_TRY
    BoundedVectorType principal_nonlocal_strains;
    BoundedVectorType principal_damage_n = ZeroVector(3);
    BoundedMatrixType MatrixForm = ZeroMatrix(3,3);
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    Vector local_equivalent_strains_tc = ZeroVector(2);
    MathUtils<double>::GaussSeidelEigenSystem(mInternalDamageVariables.DamageMatrix, EigenVectors, EigenValues);
    principal_damage_n[0] = EigenValues(0,0);
    principal_damage_n[1] = EigenValues(1,1);
    principal_damage_n[2] = EigenValues(2,2);
    CalculateLocalEquivalentStrains(local_equivalent_strains_tc, rParametersValues, PrincipalStrains);
    ScaleNonlocalEquivalentStrain(principal_nonlocal_strains, rParametersValues, PrincipalStrains, nonlocal_equivalent_strains, local_equivalent_strains_tc);
    SetModelParameters(rModelParameters, rParametersValues, PrincipalStrains, principal_nonlocal_strains);
    CheckDamageCriteria(PrincipalDamageIncrement, rModelParameters, principal_damage_n);
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CheckDamageCriteria(
    BoundedVectorType& PrincipalDamageIncrement,
    const ModelParameters& rModelParameters,
    const BoundedVectorType& PrincipalDamageVector_n
    )
{
    KRATOS_TRY
    for (SizeType i = 0; i < Dimension; ++i) {
        if (rModelParameters.kappa[i]>= 0 && rModelParameters.kappa[i]<=rModelParameters.kappa0[i]){
            PrincipalDamageIncrement[i]=0.0;
        }else if (rModelParameters.kappa[i]> rModelParameters.kappa0[i]){
            PrincipalDamageIncrement[i] = (1-PrincipalDamageVector_n[i]) * ((rModelParameters.Beta1[i]/rModelParameters.kappa[i]) + (rModelParameters.Beta2[i]/rModelParameters.kappa0[i])) * rModelParameters.del_kappa[i];
        }
        if(PrincipalDamageIncrement[i] < 0.0){
            PrincipalDamageIncrement[i] = 0.0;
        }

    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::SetModelParameters(
    ModelParameters& rModelParameters,
    ConstitutiveLaw::Parameters& rParametersValues,
    const BoundedVectorType& PrincipalStrains,
    const BoundedVectorType& PrincipalNonlocalStrains
    )
{
    KRATOS_TRY

    BoundedVectorType current_max_strain_history_parameter = ZeroVector(3);
    BoundedVectorType previous_kappa = rModelParameters.kappa;
    for(SizeType i = 0; i < Dimension; ++i){
        rModelParameters.Beta1[i] = (PrincipalStrains[i] > eps) ? mMaterialProperties[2] : mMaterialProperties[4];
        rModelParameters.Beta2[i] = (PrincipalStrains[i] > eps) ? mMaterialProperties[3] : mMaterialProperties[5];
        rModelParameters.kappa0[i] = (PrincipalStrains[i] > eps) ? mMaterialProperties[0] : mMaterialProperties[1];
        current_max_strain_history_parameter[i] = std::max(rModelParameters.strain_history_parameter[i], PrincipalNonlocalStrains[i]);
        rModelParameters.strain_history_parameter[i] = current_max_strain_history_parameter[i];
        rModelParameters.kappa[i] = std::max(current_max_strain_history_parameter[i], rModelParameters.kappa0[i]);
        rModelParameters.del_kappa[i] = rModelParameters.kappa[i] - previous_kappa[i];
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CalculateLocalEquivalentStrains(
    Vector& LocalEquivalentStrains,
    ConstitutiveLaw::Parameters& rParametersValues,
    const BoundedVectorType& PrincipalStrains
    )
{
    KRATOS_TRY
    Vector& r_strain_vector = rParametersValues.GetStrainVector();
    BoundedVectorType PrincipalDeviatoricStrains = ZeroVector(3);
    GetPrincipalDeviatoricStrains(PrincipalDeviatoricStrains, r_strain_vector, PrincipalStrains);
    LocalEquivalentStrains[0] = sqrt(pow(MacaulayBrackets(PrincipalStrains[0]),2) + pow(MacaulayBrackets(PrincipalStrains[1]),2) + pow(MacaulayBrackets(PrincipalStrains[2]),2));
    LocalEquivalentStrains[1] = sqrt(pow(MacaulayBrackets(PrincipalDeviatoricStrains[0]),2) + pow(MacaulayBrackets(PrincipalDeviatoricStrains[1]),2) + pow(MacaulayBrackets(PrincipalDeviatoricStrains[2]),2));
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetPrincipalDeviatoricStrains(
    BoundedVectorType& PrincipalDeviatoricStrains,
    const Vector& StrainVector,
    const BoundedVectorType& PrincipalStrains
    )
{
    KRATOS_TRY
    const double trace = StrainVector[0] + StrainVector[1] + StrainVector[2];
    const double hydrostatic_strain = trace/3;
    for (SizeType i = 0; i < Dimension; ++i){
        PrincipalDeviatoricStrains[i] = PrincipalStrains[i] - hydrostatic_strain;
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//***********************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::TransformPrincipalDamageIncrementToGlobal(
    BoundedMatrixType& DamageIncrementTensor,
    const BoundedVectorType& PrinicipalDamageIncrementVector,
    const Vector& StrainVector
    )
{
    KRATOS_TRY
    BoundedVectorType l = ZeroVector(3);
    BoundedVectorType m = ZeroVector(3);
    BoundedVectorType n = ZeroVector(3);
    BoundedMatrixType MatrixForm = ZeroMatrix(3,3);
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    VectorToTensor(MatrixForm, StrainVector, STRAIN);
    MathUtils<double>::GaussSeidelEigenSystem(MatrixForm, EigenVectors, EigenValues);
        for(SizeType i = 0; i < Dimension; ++i){
            l[i] = EigenVectors(i,0);
            m[i] = EigenVectors(i,1);
            n[i] = EigenVectors(i,2);
        }
    for(SizeType i = 0; i < Dimension; ++i){
        DamageIncrementTensor(i,i) = PrinicipalDamageIncrementVector[0] * l[i] * l[i] + PrinicipalDamageIncrementVector[1] * m[i] * m[i] + PrinicipalDamageIncrementVector[2] * n[i] * n[i];
    }
    DamageIncrementTensor(0,1) = PrinicipalDamageIncrementVector[0] * l[0] * l[1] + PrinicipalDamageIncrementVector[1] * m[0] * m[1] + PrinicipalDamageIncrementVector[2] * n[0] * n[1];
    DamageIncrementTensor(1,0) = PrinicipalDamageIncrementVector[0] * l[0] * l[1] + PrinicipalDamageIncrementVector[1] * m[0] * m[1] + PrinicipalDamageIncrementVector[2] * n[0] * n[1];
    DamageIncrementTensor(0,2) = PrinicipalDamageIncrementVector[0] * l[0] * l[2] + PrinicipalDamageIncrementVector[1] * m[0] * m[2] + PrinicipalDamageIncrementVector[2] * n[0] * n[2];
    DamageIncrementTensor(2,0) = PrinicipalDamageIncrementVector[0] * l[0] * l[2] + PrinicipalDamageIncrementVector[1] * m[0] * m[2] + PrinicipalDamageIncrementVector[2] * n[0] * n[2];
    DamageIncrementTensor(1,2) = PrinicipalDamageIncrementVector[0] * l[1] * l[2] + PrinicipalDamageIncrementVector[1] * m[1] * m[2] + PrinicipalDamageIncrementVector[2] * n[1] * n[2];
    DamageIncrementTensor(2,1) = PrinicipalDamageIncrementVector[0] * l[1] * l[2] + PrinicipalDamageIncrementVector[1] * m[1] * m[2] + PrinicipalDamageIncrementVector[2] * n[1] * n[2];
    KRATOS_CATCH("")
}

//************************************************************************************
//***********************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetTotalPrincipalDamageVector(
    BoundedVectorType& PrincipaldamageVector,
    const BoundedMatrixType& TotalDamageTensor
    )
{
    KRATOS_TRY
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    MathUtils<double>::GaussSeidelEigenSystem(TotalDamageTensor, EigenVectors, EigenValues);
    PrincipaldamageVector[0] = EigenValues(0,0);
    PrincipaldamageVector[1] = EigenValues(1,1);
    PrincipaldamageVector[2] = EigenValues(2,2);
    KRATOS_CATCH("")
}

//************************************************************************************
//***********************************************************************************

Vector& ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    KRATOS_TRY
    //Explicitly having STRAIN and INITAL_STRAIN_VECTOR calculated in base class
    if (rThisVariable == EQUIVALENT_STRAINS_TC) {
        InternalDamageVariables rInternalDamageVariables;
        this->CalculateStressResponse( rParametersValues, rInternalDamageVariables);
        rValue = rInternalDamageVariables.EquivalentStrainsTC;
    }else if(rThisVariable == DAMAGE_VECTOR){
        InternalDamageVariables rInternalDamageVariables;
        this->CalculateStressResponse(rParametersValues,rInternalDamageVariables);
        rValue = rInternalDamageVariables.PrincipalDamageVector;
    }else{
        ElasticIsotropic3D::CalculateValue(rParametersValues, rThisVariable, rValue);
    }
    return(rValue);
    KRATOS_CATCH("")
}

//************************************************************************************
//***********************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetDamageEffectTensor(
    BoundedMatrixVoigtType& DamageEffectTensor,
    const BoundedVectorType& DamageVector
)
{
    KRATOS_TRY
    const double D1 =  DamageVector[0];
    const double D2 =  DamageVector[1];
    const double D3 =  DamageVector[2];
    DamageEffectTensor(0,0) = pow((1-D1),(-1));
    DamageEffectTensor(1,1) = pow((1-D2),(-1));
    DamageEffectTensor(2,2) = pow((1-D3),(-1));
    DamageEffectTensor(3,3) = pow((sqrt ((1-D1)*(1-D2))),(-1));
    DamageEffectTensor(4,4) = pow((sqrt ((1-D2)*(1-D3))),(-1));
    DamageEffectTensor(5,5) = pow((sqrt ((1-D1)*(1-D3))),(-1));
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetTransformedDamageEffectTensor(
    BoundedMatrixVoigtType& TransformedDamageEffectTensor,
    const BoundedMatrixVoigtType& DamageEffectTensor,
    const Vector& StrainVector
)
{
    KRATOS_TRY

    BoundedVectorVoigtType W = ZeroVector(6);
    BoundedVectorType l = ZeroVector(3);
    BoundedVectorType m = ZeroVector(3);
    BoundedVectorType n = ZeroVector(3);
    BoundedMatrixType MatrixForm = ZeroMatrix(3,3);
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    for(SizeType i = 0; i < VoigtSize; ++i){
        W[i] = 1.0/DamageEffectTensor(i,i);
    }
    VectorToTensor(MatrixForm, StrainVector, STRAIN);
    MathUtils<double>::GaussSeidelEigenSystem(MatrixForm, EigenVectors, EigenValues);
    for(SizeType i = 0; i < Dimension; ++i){
        l[i] = EigenVectors(i,0);
        m[i] = EigenVectors(i,1);
        n[i] = EigenVectors(i,2);
    }
    for(SizeType i = 0; i < Dimension; ++i) {
        for(SizeType j = 0; j < Dimension; ++j) {
            TransformedDamageEffectTensor(i,j)  = ( (l[i] * l[i]) * (l[j] * l[j]) ) / W[0];
            TransformedDamageEffectTensor(i,j) += ( (m[i] * m[i]) * (m[j] * m[j]) ) / W[1] ;
            TransformedDamageEffectTensor(i,j) += ( (n[i] * n[i]) * (n[j] * n[j]) ) / W[2] ;
            TransformedDamageEffectTensor(i,j) += ( 2 * m[i] * n[i] * m[j] * n[j] ) / W[4] ;
            TransformedDamageEffectTensor(i,j) += ( 2 * n[i] * l[i] * n[j] * l[j] ) / W[5] ;
            TransformedDamageEffectTensor(i,j) += ( 2 * l[i] * m[i] * l[j] * m[j] ) / W[3] ;
        }
    }

    for(SizeType i = 3; i < 5; ++i ){
        TransformedDamageEffectTensor(i,i)  = ( 2.0 * l[2] * l[2] * l[4-i] * l[4-i] ) / W[0] ;
        TransformedDamageEffectTensor(i,i) += ( 2.0 * m[2] * m[2] * m[4-i] * m[4-i] ) / W[1] ;
        TransformedDamageEffectTensor(i,i) += ( 2.0 * n[2] * n[2] * n[4-i] * n[4-i] ) / W[2] ;
        TransformedDamageEffectTensor(i,i) += std::pow(( m[2] * n[4-i] + n[2] * m[4-i] ), 2.0) / W[4] ;
        TransformedDamageEffectTensor(i,i) += std::pow(( n[2] * l[4-i] + l[2] * n[4-i] ), 2.0) / W[5] ;
        TransformedDamageEffectTensor(i,i) += std::pow(( l[2] * m[4-i] + m[2] * l[4-i] ), 2.0) / W[3] ;
    }

    TransformedDamageEffectTensor(5,5)  =  ( 2.0 * l[0] * l[0] * l[1] * l[1] ) / W[0] ;
    TransformedDamageEffectTensor(5,5) +=  ( 2.0 * m[0] * m[0] * m[1] * m[1] ) / W[1] ;
    TransformedDamageEffectTensor(5,5) +=  ( 2.0 * n[0] * n[0] * n[1] * n[1] ) / W[2] ;
    TransformedDamageEffectTensor(5,5) +=  std::pow(( m[0] * n[1] + n[0] * m[1] ), 2.0) / W[4] ;
    TransformedDamageEffectTensor(5,5) +=  std::pow(( n[0] * l[1] + l[0] * n[1] ), 2.0) / W[5] ;
    TransformedDamageEffectTensor(5,5) +=  std::pow(( l[0] * m[1] + m[0] * l[1] ), 2.0) / W[3] ;

    for(SizeType i = 3; i < 5; ++i) {
        for(SizeType j = 0; j < Dimension; ++j) {
            TransformedDamageEffectTensor(i,j)  = ( l[2] * l[4-i] * l[j] * l[j] ) / W[0] ;
            TransformedDamageEffectTensor(i,j) += ( m[2] * m[4-i] * m[j] * m[j] ) / W[1] ;
            TransformedDamageEffectTensor(i,j) += ( n[2] * n[4-i] * n[j] * n[j] ) / W[2] ;
            TransformedDamageEffectTensor(i,j) += ( m[j] * n[j] * ( m[2] * n[4-i] +  n[2] * m[4-i] ) ) / W[4] ;
            TransformedDamageEffectTensor(i,j) += ( n[j] * l[j] * ( n[2] * l[4-i] +  l[2] * n[4-i] ) ) / W[5] ;
            TransformedDamageEffectTensor(i,j) += ( l[j] * m[j] * ( l[2] * m[4-i] +  m[2] * l[4-i] ) ) / W[3] ;
        }
    }

    TransformedDamageEffectTensor(4,3)  = ( 2 * l[0] * l[1] * l[2] * l[2] ) / W[0] ;
    TransformedDamageEffectTensor(4,3) += ( 2 * m[0] * m[1] * m[2] * m[2] ) / W[1] ;
    TransformedDamageEffectTensor(4,3) += ( 2 * n[0] * n[1] * n[2] * n[2] ) / W[2] ;
    TransformedDamageEffectTensor(4,3) += ( ( m[1] * n[2] + m[2] * n[1] ) * ( m[0] * n[2] + m[2] * n[0]) ) / W[4] ;
    TransformedDamageEffectTensor(4,3) += ( ( n[1] * l[2] + n[2] * l[1] ) * ( n[2] * l[0] + l[2] * n[0]) ) / W[5] ;
    TransformedDamageEffectTensor(4,3) += ( ( l[1] * m[2] + l[2] * m[1] ) * ( l[2] * m[0] + m[2] * l[0]) ) / W[3] ;

    for(SizeType j = 0; j < Dimension; ++j){
       TransformedDamageEffectTensor(5,j) =  ( l[0] * l[1] * l[j] * l[j] ) / W[0] ;
       TransformedDamageEffectTensor(5,j) += ( m[0] * m[1] * m[j] * m[j] ) / W[1] ;
       TransformedDamageEffectTensor(5,j) += ( n[0] * n[1] * n[j] * n[j] ) / W[2] ;
       TransformedDamageEffectTensor(5,j) += ( m[j] * n[j] * ( m[0] * n[1] + m[1] * n[0] ) ) / W[4] ;
       TransformedDamageEffectTensor(5,j) += ( n[j] * l[j] * ( n[0] * l[1] + l[0] * n[1] ) ) / W[5] ;
       TransformedDamageEffectTensor(5,j) += ( l[j] * m[j] * ( l[0] * m[1] + m[0] * l[1] ) ) / W[3] ;
    }

    for(SizeType j = 3; j < 5; ++j){
       TransformedDamageEffectTensor(5,j) =  ( 2 * l[0] * l[1] * l[2] * l[4-j] ) / W[0] ;
       TransformedDamageEffectTensor(5,j) += ( 2 * m[0] * m[1] * m[2] * m[4-j] ) / W[1] ;
       TransformedDamageEffectTensor(5,j) += ( 2 * n[0] * n[1] * n[2] * n[4-j] ) / W[2] ;
       TransformedDamageEffectTensor(5,j) += ( ( m[2] * n[4-j] + m[4-j] * n[2] ) * ( m[0] * n[1] + n[1] * m[0] ) ) / W[4] ;
       TransformedDamageEffectTensor(5,j) += ( ( n[2] * l[4-j] + n[4-j] * l[2])  * ( n[0] * l[1] + l[0] * n[1] ) ) / W[5] ;
       TransformedDamageEffectTensor(5,j) += ( ( l[2] * m[4-j] + l[4-j] * m[2] ) * ( l[0] * m[1] + m[0] * l[1] ) ) / W[3] ;
    }

    for(SizeType i = 3; i < VoigtSize; ++i){
        for(SizeType j = 0; i < Dimension; ++j){
            TransformedDamageEffectTensor(j,i) = 2 * TransformedDamageEffectTensor(i,j);
        }
    }

    for(SizeType i = 4; i < VoigtSize; ++i){
        for(SizeType j = 3; j < 5; ++j){
            TransformedDamageEffectTensor(j,i) = TransformedDamageEffectTensor(i,j);
        }
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::VectorToTensor(
    BoundedMatrixType& TensorForm,
    const Vector& VectorForm,
    const Variable<Vector>& rThisVariable
)
{
    KRATOS_TRY
    if(rThisVariable == STRESSES){
        TensorForm(0,1)= TensorForm(1,0)= VectorForm[3];
        TensorForm(0,2)= TensorForm(2,0)= VectorForm[5];
        TensorForm(2,1)= TensorForm(1,2)= VectorForm[4];
    }else if(rThisVariable == STRAIN){
        TensorForm(0,1)= TensorForm(1,0)= 0.5 * VectorForm[3];
        TensorForm(0,2)= TensorForm(2,0)= 0.5 * VectorForm[5];
        TensorForm(2,1)= TensorForm(1,2)= 0.5 * VectorForm[4];
    }
    for(SizeType i = 0; i < Dimension; ++i){
        TensorForm(i,i) = VectorForm[i];
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::ScaleNonlocalEquivalentStrain(
    BoundedVectorType& Principal_Nonlocal_Strains,
    ConstitutiveLaw::Parameters& rParametersValues,
    const BoundedVectorType& Principal_Strains,
    const Vector& Nonlocal_Equivalent_Strains,
    Vector& Local_Equivalent_Strains
)
{
    KRATOS_TRY

    const Vector& r_strain_vector = rParametersValues.GetStrainVector();
    BoundedVectorType PrincipalDeviatoricStrains = ZeroVector(3);
    GetPrincipalDeviatoricStrains(PrincipalDeviatoricStrains, r_strain_vector, Principal_Strains);
    if(Local_Equivalent_Strains[0] == 0){
        Local_Equivalent_Strains[0] = 2 * eps;
    }
    if(Local_Equivalent_Strains[1] == 0){
        Local_Equivalent_Strains[1] = 2 * eps;
    }
    for (SizeType i = 0; i < Dimension; ++i){
        if (Principal_Strains[i]>0.0){
            Principal_Nonlocal_Strains[i] = Principal_Strains[i] * Nonlocal_Equivalent_Strains[0]/(Local_Equivalent_Strains[0]);
        }else if(Principal_Strains[i]<0.0){
            Principal_Nonlocal_Strains[i] = fabs(PrincipalDeviatoricStrains[i]) * Nonlocal_Equivalent_Strains[1]/(Local_Equivalent_Strains[1]);
        }else{
            Principal_Nonlocal_Strains[i] = 0.0;
        }
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::Calculate_tangent_HNLu(
    BoundedVectorVoigtType& H_NL1u,
    BoundedVectorVoigtType& H_NL2u,
    ConstitutiveLaw::Parameters& rParametersValues,
    const BoundedVectorType& Principal_Strains
)
{
    KRATOS_TRY
    BoundedMatrix3x6Type derivatives_of_eigen_values;
    const Vector& r_strain_vector = rParametersValues.GetStrainVector();
    Vector local_equivalent_strains = ZeroVector(2);
    CalculateLocalEquivalentStrains(local_equivalent_strains, rParametersValues, Principal_Strains);
    double Local_Equivalent_Strain_tension = local_equivalent_strains[0];
    double Local_Equivalent_Strain_compression = local_equivalent_strains[1];
    CalculateDerivativesofEigenvalues(derivatives_of_eigen_values, Principal_Strains, r_strain_vector, STRAIN);
    BoundedVectorType PrincipalDeviatoricStrains = ZeroVector(3);
    BoundedVectorVoigtType dHydrostaticdEpsilon = ZeroVector(6);
    dHydrostaticdEpsilon[0] = 1/3;
    dHydrostaticdEpsilon[1] = 1/3;
    dHydrostaticdEpsilon[2] = 1/3;
    GetPrincipalDeviatoricStrains(PrincipalDeviatoricStrains, r_strain_vector, Principal_Strains);
    for(SizeType i = 0; i < Dimension; ++i){
        for (SizeType j = 0; j < VoigtSize; ++j){
            H_NL1u[j]   +=  MacaulayBrackets(Principal_Strains[i]) * derivatives_of_eigen_values(i,j);
            H_NL2u[j]   +=  MacaulayBrackets(PrincipalDeviatoricStrains[i]) * (derivatives_of_eigen_values(i,j) - dHydrostaticdEpsilon[j]);
        }
    }
    H_NL1u /= (Local_Equivalent_Strain_tension +eps) ;
    H_NL2u /= (Local_Equivalent_Strain_compression +eps) ;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::TensorProduct(
    BoundedMatrixVoigtType& dHdNL1,
    BoundedMatrixVoigtType& dHdNL2,
    const array_1d<BoundedMatrix<double, 6, 6>, 3>& dHdD,
    const BoundedVectorType& Vector1,
    const BoundedVectorType& Vector2
    )
{
    KRATOS_TRY
    for (SizeType i = 0; i < VoigtSize; ++i) {
        for (SizeType j = 0; j < VoigtSize; ++j) {
            for (SizeType k = 0; k < Dimension; ++k) {
                dHdNL1(i,j) += dHdD[k](i,j)* Vector1[k];
                dHdNL2(i,j) += dHdD[k](i,j)* Vector2[k];
            }
        }
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::TensorProduct6(
    Matrix& rOutput,
    const Vector& rVector1,
    const Vector& rVector2)
{
    KRATOS_TRY
    KRATOS_DEBUG_ERROR_IF(rVector1.size()!=6 || rVector2.size()!=6) << "check the size of vectors  ";

    if (rOutput.size1() != 6 || rOutput.size1() != 6) {
        rOutput.resize(6, 6,false);
    }
    for(int i=0; i<6; ++i){
        for(int j=0; j<6; ++j){
            rOutput(i,j)= rVector1[i] * rVector2[j];
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::CalculateDerivativesofEigenvalues(
     BoundedMatrix3x6Type &DerivativesofEigenvalues,
     const BoundedVectorType &EigenvaluesVector,
     const BoundedVectorVoigtType &Voigtform,
     const Variable<Vector>& rThisVariable
    )
{
    KRATOS_TRY
    BoundedMatrixType Matrixform;
    VectorToTensor(Matrixform, Voigtform, rThisVariable);
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    MathUtils<double>::GaussSeidelEigenSystem(Matrixform, EigenVectors, EigenValues);
    BoundedVectorType Eigen_vector_column = ZeroVector (3) ;
    BoundedMatrixType DerivativesMatrix = ZeroMatrix(3,3);
    for(SizeType i = 0; i < Dimension; ++i){
        Eigen_vector_column = column(EigenVectors,i);
        DerivativesMatrix = outer_prod(Eigen_vector_column, Eigen_vector_column);
        DerivativesofEigenvalues(i,0) = DerivativesMatrix(0,0);
        DerivativesofEigenvalues(i,1) = DerivativesMatrix(1,1);
        DerivativesofEigenvalues(i,2) = DerivativesMatrix(2,2);
        DerivativesofEigenvalues(i,3) = DerivativesMatrix(0,1);
        DerivativesofEigenvalues(i,4) = DerivativesMatrix(1,2);
        DerivativesofEigenvalues(i,5) = DerivativesMatrix(2,0);
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::AssembleConstitutiveMatrix(
    Matrix& ConstitutiveMatrix,
    const Matrix& H_uu,
    const Vector& H_NL1u,
    const Vector& H_NL2u,
    const Vector& H_uNL1,
    const Vector& H_uNL2,
    const double& H_NLNL
    )
{
    KRATOS_TRY
    const SizeType num_rows_C = H_uu.size1();
    const SizeType num_cols_C = H_uu.size2();
    const SizeType size_CuNL  = H_uNL1.size();
    const SizeType size_CNLu  = H_NL1u.size();

    for (SizeType i = 0; i < num_rows_C ; ++i) {
        for (SizeType j = 0; j < num_cols_C ; ++j) {
            ConstitutiveMatrix(i, j) = H_uu(i,j);
        }
    }
    for (SizeType i = 0; i < size_CuNL; ++i) {
        ConstitutiveMatrix(i, num_cols_C) = H_uNL1[i] ;
        ConstitutiveMatrix(i, num_cols_C + 1) = H_uNL2[i] ;
    }
    for (SizeType i = 0; i < size_CNLu; ++i) {
        ConstitutiveMatrix(num_rows_C, i) = H_NL1u[i] ;
        ConstitutiveMatrix(num_rows_C +1, i) = H_NL2u[i] ;
    }
    ConstitutiveMatrix(num_rows_C, num_cols_C) = H_NLNL;
    ConstitutiveMatrix(num_rows_C+1, num_cols_C) = H_NLNL;
    ConstitutiveMatrix(num_rows_C, num_cols_C+1) = H_NLNL;
    ConstitutiveMatrix(num_rows_C+1, num_cols_C+1) = H_NLNL;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = this->GetStrainSize();
    rFeatures.mSpaceDimension = this->WorkingSpaceDimension();
}

//************************************************************************************
//************************************************************************************



void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("DamageMatrix", mInternalDamageVariables.DamageMatrix);
    rSerializer.save("PrincipalDamageVector", mInternalDamageVariables.PrincipalDamageVector);
    rSerializer.save("EquivalentStrainsTC", mInternalDamageVariables.EquivalentStrainsTC);
    rSerializer.save("kappa", mInternalDamageVariables.kappa);
    rSerializer.save("StrainHistoryParameter", mInternalDamageVariables.StrainHistoryParameter);
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("DamageMatrix", mInternalDamageVariables.DamageMatrix);
    rSerializer.load("PrincipalDamageVector", mInternalDamageVariables.PrincipalDamageVector);
    rSerializer.load("EquivalentStrainsTC", mInternalDamageVariables.EquivalentStrainsTC);
    rSerializer.load("kappa", mInternalDamageVariables.kappa);
    rSerializer.load("StrainHistoryParameter", mInternalDamageVariables.StrainHistoryParameter);
}



} /* namespace Kratos.*/
