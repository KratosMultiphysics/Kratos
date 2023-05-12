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
//  Main authors:    Marcelo Raschi
//  Collaborators:
//

// Project includes
#include "small_strain_isotropic_damage_implex_3d.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "constitutive_laws_application_variables.h"
#include "includes/checks.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainIsotropicDamageImplex3D::SmallStrainIsotropicDamageImplex3D()
    : ElasticIsotropic3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

SmallStrainIsotropicDamageImplex3D::SmallStrainIsotropicDamageImplex3D(
    const SmallStrainIsotropicDamageImplex3D &rOther) = default;

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainIsotropicDamageImplex3D::Clone() const
{
    return Kratos::make_shared<SmallStrainIsotropicDamageImplex3D>(
        SmallStrainIsotropicDamageImplex3D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SmallStrainIsotropicDamageImplex3D::~SmallStrainIsotropicDamageImplex3D() = default;

//************************************************************************************
//************************************************************************************

bool SmallStrainIsotropicDamageImplex3D::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == SCALE_FACTOR) {
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    } else if (rThisVariable == STRAIN_ENERGY) {
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    } else if (rThisVariable == DAMAGE_VARIABLE) {
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

bool SmallStrainIsotropicDamageImplex3D::Has(const Variable<Vector>& rThisVariable)
{
    if(rThisVariable == INTERNAL_VARIABLES){
        return true;
    } else if(rThisVariable == STRAIN){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

Vector& SmallStrainIsotropicDamageImplex3D::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        rValue.resize(2);
        rValue[0] = mStrainVariable;
        rValue[1] = mStrainVariablePrevious;
    }

    return rValue;
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageImplex3D::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rProcessInfo
    )
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        mStrainVariable = rValue[0];
        mStrainVariablePrevious = rValue[1];
    }
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageImplex3D::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    const double yield_stress = rMaterialProperties[STRESS_LIMITS](0);
    const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
    mStrainVariable = yield_stress / std::sqrt(young_modulus);
    mStrainVariablePrevious = mStrainVariable;
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageImplex3D::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    Vector internal_variables(2);
    this->CalculateStressResponse(rParametersValues, internal_variables);
    mStrainVariable = internal_variables[0];
    mStrainVariablePrevious = internal_variables[1];
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageImplex3D::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    Vector internal_variables(2);
    CalculateStressResponse(rParametersValues, internal_variables);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageImplex3D::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rParametersValues,
    Vector& rInternalVariables)
{
    const Properties& r_material_properties = rParametersValues.GetMaterialProperties();
    Flags& r_constitutive_law_options = rParametersValues.GetOptions();
    Vector& r_strain_vector = rParametersValues.GetStrainVector();
    CalculateValue(rParametersValues, STRAIN, r_strain_vector);

    // Implex
    const double dt_n1 = rParametersValues.GetProcessInfo()[DELTA_TIME];
    double dt_n0 = rParametersValues.GetProcessInfo().GetPreviousTimeStepInfo(1)[DELTA_TIME];
    if (dt_n0 <= 0) {
          dt_n0 = dt_n1; //To avoid division by zero in 1st increment
    }
    double r_implex = mStrainVariable + (dt_n1 / dt_n0) * (mStrainVariable - mStrainVariablePrevious);
    const double q_implex = EvaluateHardeningLaw(r_implex, r_material_properties);
    const double d_implex = 1. - q_implex / r_implex;

    double r = mStrainVariable;

    // If we compute the tangent moduli or the stress
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ||
        r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR )) {
        Vector& r_stress_vector = rParametersValues.GetStressVector();
        Matrix& r_constitutive_matrix = rParametersValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rParametersValues);
        noalias(r_stress_vector) = prod(r_constitutive_matrix, r_strain_vector);

        Vector r_stress_vector_pos = r_stress_vector;
        ComputePositiveStressVector(r_stress_vector_pos, r_stress_vector);

        double energy = inner_prod(r_stress_vector_pos, r_strain_vector);
        // energy may be a small negative due to machine precision error, forcing zero
        if (energy < 0) {
            energy = 0;
        }

        const double strain_norm = std::sqrt(energy);

        // inelastic case
        if (strain_norm > mStrainVariable) {
            r = strain_norm;
        }

        r_constitutive_matrix *= (1 - d_implex);
        r_stress_vector *= (1 - d_implex);

    }

    rInternalVariables[0] = r;
    rInternalVariables[1] = mStrainVariable;
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageImplex3D::ComputePositiveStressVector(
            Vector& rStressVectorPos, Vector& rStressVector)
{
    // Auxiliary stress vector to allow derived models (e.g. traction-only damage)
    // to set the value of r_stress_vector_pos with the ComputePositiveStressVector
    // function.

    // In this symmetric damage model, ComputePositiveStressVector function does
    // nothing, as rStressVectorPos = rStressVector already.
}

//************************************************************************************
//************************************************************************************

double SmallStrainIsotropicDamageImplex3D::EvaluateHardeningModulus(
        double r,
        const Properties &rMaterialProperties)
{
    if (rMaterialProperties[HARDENING_CURVE] == 0) {

        // 0: exponential hardening

        const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
        const double yield_stress = rMaterialProperties[STRESS_LIMITS](0);
        const double inf_yield_stress = rMaterialProperties[STRESS_LIMITS](1);
        const double h0 = rMaterialProperties[HARDENING_PARAMETERS](0);
        const double r0 = yield_stress / std::sqrt(young_modulus);
        const double q1 = inf_yield_stress / std::sqrt(young_modulus);

        if (r < r0)
            return 0.;
        else // >= r0:
            return h0 * (q1/r0 - 1) * std::exp(h0 * (1 - r/r0));

    } else {

        // 1: multilinear hardening

        const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
        const double yield_stress = rMaterialProperties[STRESS_LIMITS](0);
        const double r0 = yield_stress / std::sqrt(young_modulus);
        const double h0 = rMaterialProperties[HARDENING_PARAMETERS](0);

        if (r < r0)
            return  0.;

        switch(rMaterialProperties[HARDENING_PARAMETERS].size())
        {
            case 1:  // linear
                return h0;
                break;

            case 2:  // bilinear
            {
                const double q0 = r0;
                const double s1 = rMaterialProperties[STRESS_LIMITS](1);
                const double q1 = s1 / std::sqrt(young_modulus);
                const double r1 = r0 + (q1 - q0) / h0;
                const double h1 = rMaterialProperties[HARDENING_PARAMETERS](1);

                if (r >= r0 && r < r1)
                    return h0;
                else  // r >= r1:
                    return h1;
                break;
            }

            case 3:   // trilinear
            {
                const double q0 = r0;
                const double s1 = rMaterialProperties[STRESS_LIMITS](1);
                const double q1 = s1 / std::sqrt(young_modulus);
                const double r1 = r0 + (q1 - q0) / h0;
                const double h1 = rMaterialProperties[HARDENING_PARAMETERS](1);
                const double s2 = rMaterialProperties[STRESS_LIMITS](2);
                const double q2 = s2 / std::sqrt(young_modulus);
                const double r2 = r1 + (q2 - q1) / h1;
                const double h2 = rMaterialProperties[HARDENING_PARAMETERS](2);

                if (r >= r0 && r < r1)
                    return h0;
                else if (r >= r1 && r < r2)
                    return h1;
                else  // r >= r2
                    return h2;
                break;
            }

            default:  // for other cases, implement as needed
                KRATOS_ERROR << "Multilinear hardening model of more than 3 " <<
                "regions not yet implemented." << std::endl;
                break;
        }
    }
}

//************************************************************************************
//************************************************************************************

double SmallStrainIsotropicDamageImplex3D::EvaluateHardeningLaw(
    double r,
    const Properties &rMaterialProperties)
{
    if (rMaterialProperties[HARDENING_CURVE] == 0) {

        // 0: exponential hardening

        const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
        const double yield_stress = rMaterialProperties[STRESS_LIMITS](0);
        const double inf_yield_stress = rMaterialProperties[STRESS_LIMITS](1);
        const double r0 = yield_stress / std::sqrt(young_modulus);
        const double h0 = EvaluateHardeningModulus(r0, rMaterialProperties);
        const double q0 = r0;  // strain_variable_init
        const double q1 = inf_yield_stress / std::sqrt(young_modulus);  // stress_variable_inf

        if (r < r0)
            return q0;
        else  // r >= r0:
            return q1 - (q1 - r0) * std::exp(h0 * (1 - r/r0));  //  for H0 > 0

    } else {

        // 1: bilinear hardening

        const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
        const double s0 = rMaterialProperties[STRESS_LIMITS](0);
        const double r0 = s0 / std::sqrt(young_modulus);
        const double h0 = EvaluateHardeningModulus(r0, rMaterialProperties);
        const double q0 = r0;

        if (r < r0)
            return q0;

        switch(rMaterialProperties[HARDENING_PARAMETERS].size())
        {
            case 1:  // linear
                return q0 + h0 * (r - r0);
                break;

            case 2:  // bilinear
            {
                const double s1 = rMaterialProperties[STRESS_LIMITS](1);
                const double q1 = s1 / std::sqrt(young_modulus);
                const double r1 = r0 + (q1 - q0) / h0;
                const double h1 = EvaluateHardeningModulus(r1, rMaterialProperties);

                if (r >= r0 && r < r1)
                    return q0 + h0 * (r - r0);
                else  // r >= r1:
                    return q1 + h1 * (r - r1);
                break;
            }

            case 3:  // trilinear
            {
                const double s1 = rMaterialProperties[STRESS_LIMITS](1);
                const double q1 = s1 / std::sqrt(young_modulus);
                const double r1 = r0 + (q1 - q0) / h0;
                const double h1 = EvaluateHardeningModulus(r1, rMaterialProperties);
                const double s2 = rMaterialProperties[STRESS_LIMITS](2);
                const double q2 = s2 / std::sqrt(young_modulus);
                const double r2 = r1 + (q2 - q1) / h1;
                const double h2 = EvaluateHardeningModulus(r2, rMaterialProperties);

                if (r >= r0 && r < r1)
                    return q0 + h0 * (r - r0);
                else if (r >= r1 && r < r2)
                    return q1 + h1 * (r - r1);
                else  // r >= r2
                    return q2 + h2 * (r - r2);
                break;
            }

            default:  // for other cases, implement as needed
                KRATOS_ERROR << "Multilinear hardening model of more than 3 " <<
                "regions not yet implemented." << std::endl;
                break;
        }
    }
}

//************************************************************************************
//************************************************************************************

double& SmallStrainIsotropicDamageImplex3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == SCALE_FACTOR) {
        const Properties& r_material_properties = rParametersValues.GetMaterialProperties();
        const double stress_like_variable = EvaluateHardeningLaw(
                                                mStrainVariable,
                                                r_material_properties);
        const double hardening_modulus = EvaluateHardeningModulus(
                                                mStrainVariable,
                                                r_material_properties);

        rValue = ((stress_like_variable - hardening_modulus * mStrainVariable) / (mStrainVariable*mStrainVariable))*(mStrainVariable - mStrainVariablePrevious);

    } else if (rThisVariable == STRAIN_ENERGY) {
        Vector& r_strain_vector = rParametersValues.GetStrainVector();
        this->CalculateValue(rParametersValues, STRAIN, r_strain_vector);
        const Properties& r_material_properties = rParametersValues.GetMaterialProperties();
        Matrix constitutive_matrix;
        CalculateElasticMatrix(constitutive_matrix, rParametersValues);
        const double stress_like_variable = EvaluateHardeningLaw(
                                                mStrainVariable,
                                                r_material_properties);
        const double damage_variable = 1. - stress_like_variable / mStrainVariable;

        rValue = 0.5 * ((1. - damage_variable) * inner_prod(r_strain_vector,
                                        prod(constitutive_matrix, r_strain_vector)));

    } else if (rThisVariable == DAMAGE_VARIABLE){
        const Properties& r_material_properties = rParametersValues.GetMaterialProperties();
        const double stress_like_variable = EvaluateHardeningLaw(
                                                mStrainVariable,
                                                r_material_properties);

        rValue = 1. - stress_like_variable / mStrainVariable;

    } else {
        ElasticIsotropic3D::CalculateValue(rParametersValues, rThisVariable, rValue);

    }

    return(rValue);
}

//************************************************************************************
//************************************************************************************

Vector& SmallStrainIsotropicDamageImplex3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    //Explicitly having STRAIN and INITAL_STRAIN_VECTOR calculated in base class
    ElasticIsotropic3D::CalculateValue(rParametersValues, rThisVariable, rValue);
    return(rValue);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageImplex3D::GetLawFeatures(Features& rFeatures)
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

int SmallStrainIsotropicDamageImplex3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    const double tolerance = std::numeric_limits<double>::epsilon();

    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(HARDENING_CURVE));
    KRATOS_CHECK(rMaterialProperties.Has(STRESS_LIMITS));
    KRATOS_CHECK(rMaterialProperties.Has(HARDENING_PARAMETERS));

    // current supported hardening models:
    KRATOS_CHECK(rMaterialProperties[HARDENING_CURVE] == 0 || rMaterialProperties[HARDENING_CURVE] == 1);

    // checks specific for exponential hardening
    if (rMaterialProperties[HARDENING_CURVE] == 0){
        KRATOS_CHECK_GREATER_EQUAL(rMaterialProperties[STRESS_LIMITS].size(), 2);
        KRATOS_CHECK_GREATER(rMaterialProperties[STRESS_LIMITS](0), tolerance);
        KRATOS_CHECK_GREATER(rMaterialProperties[STRESS_LIMITS](1), rMaterialProperties[STRESS_LIMITS](0));
        KRATOS_CHECK_GREATER_EQUAL(rMaterialProperties[HARDENING_PARAMETERS](0), 0.);
    }

    // checks specific for multilinear hardening
    if (rMaterialProperties[HARDENING_CURVE] == 1){
        KRATOS_CHECK_EQUAL(rMaterialProperties[HARDENING_PARAMETERS].size(),
                           rMaterialProperties[STRESS_LIMITS].size());
        for (const auto& h : rMaterialProperties[HARDENING_PARAMETERS]){
            KRATOS_CHECK_GREATER_EQUAL(h, 0.);
            KRATOS_CHECK_LESS_EQUAL(h, 1.);
        }
        for (const auto& h : rMaterialProperties[STRESS_LIMITS]){
            KRATOS_CHECK_GREATER(h, tolerance);
        }
    }

    return 0;
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageImplex3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ElasticIsotropic3D);
    rSerializer.save("mStrainVariable", mStrainVariable);
    rSerializer.save("mStrainVariablePrevious", mStrainVariablePrevious);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageImplex3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ElasticIsotropic3D);
    rSerializer.save("mStrainVariable", mStrainVariable);
    rSerializer.save("mStrainVariablePrevious", mStrainVariablePrevious);
}

} /* namespace Kratos.*/
