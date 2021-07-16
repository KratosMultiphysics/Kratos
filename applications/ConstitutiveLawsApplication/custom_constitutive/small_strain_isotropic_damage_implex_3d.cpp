// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License: BSD License (structural_mechanics_application/license.txt)
//
//  Main authors:    Marcelo Raschi
//  Collaborators:
//

// Project includes
#include "custom_constitutive/small_strain_isotropic_damage_implex_3d.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "includes/checks.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainIsotropicDamageImplex3D::SmallStrainIsotropicDamageImplex3D()
    : SmallStrainIsotropicDamage3D()
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
    if(rThisVariable == SCALE_FACTOR){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;

    } else {
        BaseType::Has(rThisVariable);

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
    if(rThisVariable == INTERNAL_VARIABLES){
        rValue.resize(2);
        rValue[0] = mStrainVariable;
        rValue[1] = mStrainVariablePrevious;

    } else {
        BaseType::GetValue(rThisVariable, rValue);

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

    } else {
        BaseType::SetValue(rThisVariable, rValue, rProcessInfo);

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
    CalculateStressResponse(rParametersValues, internal_variables);
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
    BaseType::CalculateValue(rParametersValues, STRAIN, r_strain_vector);

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

        if (strain_norm > mStrainVariable) {

            // inelastic case

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

    } else {
        BaseType::CalculateValue(rParametersValues, rThisVariable, rValue);


    }

    return(rValue);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageImplex3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    rSerializer.save("mStrainVariable", mStrainVariable);
    rSerializer.save("mStrainVariablePrevious", mStrainVariablePrevious);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageImplex3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    rSerializer.save("mStrainVariable", mStrainVariable);
    rSerializer.save("mStrainVariablePrevious", mStrainVariablePrevious);
}

} /* namespace Kratos.*/
