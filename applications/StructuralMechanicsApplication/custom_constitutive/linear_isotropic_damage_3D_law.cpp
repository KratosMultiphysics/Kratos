// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi
//  Collaborators:
//

// Project includes
#include "linear_isotropic_damage_3D_law.h"
#include "structural_mechanics_application_variables.h"
#include "includes/checks.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearIsotropicDamage3D::LinearIsotropicDamage3D()
    : ConstitutiveLaw()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

LinearIsotropicDamage3D::LinearIsotropicDamage3D(const LinearIsotropicDamage3D &rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearIsotropicDamage3D::Clone() const
{
    LinearIsotropicDamage3D::Pointer p_clone(new LinearIsotropicDamage3D(*this));
    return p_clone;
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

LinearIsotropicDamage3D::~LinearIsotropicDamage3D()
{
}

//************************************************************************************
//************************************************************************************

bool LinearIsotropicDamage3D::Has(const Variable<bool>& rThisVariable)
{
    if(rThisVariable == INELASTIC_FLAG){
        return true;
    }
    return false;
}

//************************************************************************************
//************************************************************************************

bool LinearIsotropicDamage3D::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        return true;
    }
    return false;
}

//************************************************************************************
//************************************************************************************

bool& LinearIsotropicDamage3D::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    if(rThisVariable == INELASTIC_FLAG){
        rValue = mInelasticFlag;
    }

    return rValue;
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
    mDamageThresholdOld = yield_stress / std::sqrt(young_modulus);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    // update of damage threshold
    mDamageThresholdOld = mDamageThreshold;
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::CalculateMaterialResponsePK1(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::CalculateMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const Flags &Options = rValues.GetOptions();

    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    Vector& strain_vector = rValues.GetStrainVector();
    Vector& stress_vector = rValues.GetStressVector();

    if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateConstitutiveMatrix(constitutive_matrix, rMaterialProperties);
    }

    if (Options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
        }
        Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();

        CalculateConstitutiveMatrix(constitutive_matrix, rMaterialProperties);
        noalias(stress_vector) = prod(constitutive_matrix, strain_vector);

        const double strain_norm = std::sqrt(inner_prod(stress_vector, strain_vector));
        if (strain_norm <= mDamageThresholdOld)
        {
            // ELASTIC
            mInelasticFlag = false;
            mDamageThreshold = mDamageThresholdOld;
            const double q = EvaluateHardeningLaw(mDamageThreshold, rMaterialProperties);
            const double d = 1. - q / mDamageThreshold;
            constitutive_matrix *= (1 - d);
            stress_vector *= (1 - d);
        }
        else
        {
            // INELASTIC
            mInelasticFlag = true;
            mDamageThreshold = strain_norm;
            const double q = EvaluateHardeningLaw(mDamageThreshold, rMaterialProperties);
            const double d = 1. - q / mDamageThreshold;
            const double H = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
            const double dpointcoeff = (q - H * mDamageThreshold) / (mDamageThreshold * mDamageThreshold * mDamageThreshold);
            constitutive_matrix *= (1. - d);
            constitutive_matrix -= dpointcoeff * outer_prod(stress_vector, stress_vector);
            stress_vector *= (1. - d);
        }
    }
}

//************************************************************************************
//************************************************************************************

double& LinearIsotropicDamage3D::CalculateValue(
    Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == STRAIN_ENERGY){
        Vector& strain_vector = rParameterValues.GetStrainVector();
        if (rParameterValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(strain_vector) += rParameterValues.GetProcessInfo()[INITIAL_STRAIN];
        }
        const Properties& rMaterialProperties = rParameterValues.GetMaterialProperties();
        Matrix& constitutive_matrix = rParameterValues.GetConstitutiveMatrix();
        CalculateConstitutiveMatrix(constitutive_matrix, rMaterialProperties);
        const double q = EvaluateHardeningLaw(mDamageThresholdOld, rMaterialProperties);
        const double d = 1. - q / mDamageThresholdOld;

        rValue = 0.5 * ((1. - d) * inner_prod(strain_vector,
                                              prod(constitutive_matrix, strain_vector)));
    }
    return(rValue);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::FinalizeMaterialResponsePK1(Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::FinalizeMaterialResponsePK2(Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

double LinearIsotropicDamage3D::EvaluateHardeningLaw(
        double DamageThreshold,
        const Properties &rMaterialProperties
)
{
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double inf_yield_stress = rMaterialProperties[INFINITY_YIELD_STRESS];
    const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
    const double H = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double damage_threshold_original = yield_stress / std::sqrt(young_modulus);
    const double q_inf = inf_yield_stress / std::sqrt(young_modulus);
    double q;

    if (DamageThreshold < damage_threshold_original)
        return DamageThreshold;
    q = damage_threshold_original + H * (DamageThreshold - damage_threshold_original);
    if ((H > 0 && q > q_inf) || (H < 0 && q < q_inf))
        q = q_inf;
    return q;
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::CalculateConstitutiveMatrix(
    Matrix &rConstitTensor,
    const Properties &rMaterialProperties
    )
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double nu = rMaterialProperties[POISSON_RATIO];

    if (rConstitTensor.size1() != 6 || rConstitTensor.size2() != 6)
        rConstitTensor.resize(6, 6, false);
    rConstitTensor.clear();

    rConstitTensor(0, 0) = (E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)));
    rConstitTensor(1, 1) = rConstitTensor(0, 0);
    rConstitTensor(2, 2) = rConstitTensor(0, 0);
    rConstitTensor(3, 3) = rConstitTensor(0, 0) * (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));
    rConstitTensor(4, 4) = rConstitTensor(3, 3);
    rConstitTensor(5, 5) = rConstitTensor(3, 3);
    rConstitTensor(0, 1) = rConstitTensor(0, 0) * nu / (1.0 - nu);
    rConstitTensor(1, 0) = rConstitTensor(0, 1);
    rConstitTensor(0, 2) = rConstitTensor(0, 1);
    rConstitTensor(2, 0) = rConstitTensor(0, 1);
    rConstitTensor(1, 2) = rConstitTensor(0, 1);
    rConstitTensor(2, 1) = rConstitTensor(0, 1);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = 6;
    rFeatures.mSpaceDimension = 3;
}

//************************************************************************************
//************************************************************************************

int LinearIsotropicDamage3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(INFINITY_YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(ISOTROPIC_HARDENING_MODULUS));
    KRATOS_CHECK_GREATER(rMaterialProperties[YIELD_STRESS], 0.);
    KRATOS_CHECK_GREATER(rMaterialProperties[INFINITY_YIELD_STRESS], 0.);
    KRATOS_CHECK_LESS(rMaterialProperties[ISOTROPIC_HARDENING_MODULUS], 1.);
    KRATOS_CHECK_NOT_EQUAL(rMaterialProperties[ISOTROPIC_HARDENING_MODULUS], 0.);

    if (rMaterialProperties[ISOTROPIC_HARDENING_MODULUS] > 0 &&
        rMaterialProperties[INFINITY_YIELD_STRESS] <= rMaterialProperties[YIELD_STRESS])
        KRATOS_ERROR << "If ISOTROPIC_HARDENING_MODULUS is positive, "
            "INFINITY_YIELD_STRESS must be greater than YIELD_STRESS" << std::endl;

    if (rMaterialProperties[ISOTROPIC_HARDENING_MODULUS] < 0 &&
        rMaterialProperties[INFINITY_YIELD_STRESS] >= rMaterialProperties[YIELD_STRESS])
        KRATOS_ERROR << "If ISOTROPIC_HARDENING_MODULUS is negative, "
                "INFINITY_YIELD_STRESS must be lesser than YIELD_STRESS" << std::endl;

    return 0;
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mInelasticFlag", mInelasticFlag);
    rSerializer.save("mDamageThreshold", mDamageThreshold);
    rSerializer.save("mDamageThresholdOld", mDamageThresholdOld);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mInelasticFlag", mInelasticFlag);
    rSerializer.save("mDamageThreshold", mDamageThreshold);
    rSerializer.save("mDamageThresholdOld", mDamageThresholdOld);
}

} /* namespace Kratos.*/
