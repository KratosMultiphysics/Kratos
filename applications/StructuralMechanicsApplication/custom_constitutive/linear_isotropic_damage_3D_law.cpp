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
    return Kratos::make_shared<LinearIsotropicDamage3D>(LinearIsotropicDamage3D(*this));
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
    if(rThisVariable == DAMAGE_VARIABLE){
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
    mStrainVariableOld = yield_stress / std::sqrt(young_modulus);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    mStrainVariableOld = mStrainVariable;
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
    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    Vector& strain_vector = rValues.GetStrainVector();
    Vector& stress_vector = rValues.GetStressVector();

    if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
        noalias(strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
    }
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();

    CalculateConstitutiveMatrix(constitutive_matrix, rMaterialProperties);
    noalias(stress_vector) = prod(constitutive_matrix, strain_vector);

    const double strain_norm = std::sqrt(inner_prod(stress_vector, strain_vector));
    if (strain_norm <= mStrainVariableOld)
    {
        // ELASTIC
        mInelasticFlag = false;
        mStrainVariable = mStrainVariableOld;
        const double stress_variable = EvaluateHardeningLaw(mStrainVariable, rMaterialProperties);
        const double damage_variable = 1. - stress_variable / mStrainVariable;
        constitutive_matrix *= (1 - damage_variable);
        stress_vector *= (1 - damage_variable);
    }
    else
    {
        // INELASTIC
        mInelasticFlag = true;
        mStrainVariable = strain_norm;
        const double stress_variable = EvaluateHardeningLaw(mStrainVariable, rMaterialProperties);
        const double damage_variable = 1. - stress_variable / mStrainVariable;
        const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
        const double damage_rate = (stress_variable - hardening_modulus * mStrainVariable)
                                   / (mStrainVariable * mStrainVariable * mStrainVariable);
        constitutive_matrix *= (1. - damage_variable);
        constitutive_matrix -= damage_rate * outer_prod(stress_vector, stress_vector);
        stress_vector *= (1. - damage_variable);
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
        const double stress_like_variable = EvaluateHardeningLaw(mStrainVariableOld, rMaterialProperties);
        const double damage_variable = 1. - stress_like_variable / mStrainVariableOld;

        rValue = 0.5 * ((1. - damage_variable) * inner_prod(strain_vector,
                                              prod(constitutive_matrix, strain_vector)));
    }

    if (rThisVariable == DAMAGE_VARIABLE){
        const Properties& rMaterialProperties = rParameterValues.GetMaterialProperties();
        const double stress_like_variable = EvaluateHardeningLaw(mStrainVariableOld, rMaterialProperties);

        rValue = 1. - stress_like_variable / mStrainVariableOld;
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
        double StrainVariable,
        const Properties &rMaterialProperties
)
{
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double inf_yield_stress = rMaterialProperties[INFINITY_YIELD_STRESS];
    const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double strain_variable_init = yield_stress / std::sqrt(young_modulus);
    const double stress_variable_inf = inf_yield_stress / std::sqrt(young_modulus);
    double stress_variable;
    const double tolerance = std::numeric_limits<double>::epsilon();

    if (StrainVariable < strain_variable_init)
        return StrainVariable;
    stress_variable = strain_variable_init + hardening_modulus * (StrainVariable - strain_variable_init);
    if ((hardening_modulus > tolerance && stress_variable > stress_variable_inf) ||
        (hardening_modulus < tolerance && stress_variable < stress_variable_inf))
        stress_variable = stress_variable_inf;
    return stress_variable;
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
    const double tolerance = std::numeric_limits<double>::epsilon();

    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(INFINITY_YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(ISOTROPIC_HARDENING_MODULUS));
    KRATOS_CHECK_GREATER(rMaterialProperties[YIELD_STRESS], tolerance);
    KRATOS_CHECK_GREATER(rMaterialProperties[INFINITY_YIELD_STRESS], tolerance);
    KRATOS_CHECK_LESS(rMaterialProperties[ISOTROPIC_HARDENING_MODULUS], 1.);
    KRATOS_CHECK_NOT_EQUAL(rMaterialProperties[ISOTROPIC_HARDENING_MODULUS], tolerance);

    if (rMaterialProperties[ISOTROPIC_HARDENING_MODULUS] > tolerance &&
        rMaterialProperties[INFINITY_YIELD_STRESS] <= rMaterialProperties[YIELD_STRESS])
        KRATOS_ERROR << "If ISOTROPIC_HARDENING_MODULUS is positive, "
            "INFINITY_YIELD_STRESS must be greater than YIELD_STRESS" << std::endl;

    if (rMaterialProperties[ISOTROPIC_HARDENING_MODULUS] < tolerance &&
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
    rSerializer.save("mStrainVariable", mStrainVariable);
    rSerializer.save("mStrainVariableOld", mStrainVariableOld);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamage3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mInelasticFlag", mInelasticFlag);
    rSerializer.save("mStrainVariable", mStrainVariable);
    rSerializer.save("mStrainVariableOld", mStrainVariableOld);
}

} /* namespace Kratos.*/
