// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi

#include "linear_isotropic_damage_plane_strain_2d.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearIsotropicDamagePlaneStrain2D::LinearIsotropicDamagePlaneStrain2D()
    : LinearIsotropicDamage3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************
LinearIsotropicDamagePlaneStrain2D::LinearIsotropicDamagePlaneStrain2D(
        const LinearIsotropicDamagePlaneStrain2D &rOther)
        : LinearIsotropicDamage3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearIsotropicDamagePlaneStrain2D::Clone() const
{
   LinearIsotropicDamagePlaneStrain2D::Pointer pclone(new LinearIsotropicDamagePlaneStrain2D(*this));
    return pclone;
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

LinearIsotropicDamagePlaneStrain2D::~LinearIsotropicDamagePlaneStrain2D()
{
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::CalculateMaterialResponsePK1(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::CalculateMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    Vector& strain_vector = rValues.GetStrainVector();
    Vector& stress_vector = rValues.GetStressVector();
    Vector stress_vector_pos;
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
    double sigma_xx, sigma_yy, sigma_xy;
    double hyp, sigma_1, sigma_2, angle, cos_a, sin_a;
    //bool TRACTION_ONLY = rMaterialProperties[FLOW_RULE_IS_TRACTION_ONLY];

    if (rValues.GetProcessInfo().Has(INITIAL_STRAIN))
    {
        noalias(strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
    }

    CalculateConstitutiveMatrix(constitutive_matrix, rMaterialProperties);
    stress_vector = prod(constitutive_matrix, strain_vector);
    stress_vector_pos = prod(constitutive_matrix, strain_vector);
    //TODO: for use with strain energy computation (see below)
    Matrix constitutive_elastic_matrix = constitutive_matrix;

    // for tension-only fluency law:
    // originally sigma and sigma_positive are the same (as it is in the
    // symmetrical case), this block modifies sigma_positive
    if (false)
    //if (TRACTION_ONLY)
    {
        sigma_xx = stress_vector(0);
        sigma_yy = stress_vector(1);
        sigma_xy = stress_vector(2);
        hyp = std::hypot(0.5 * (sigma_xx - sigma_yy), sigma_xy);
        sigma_1 = 0.5 * (sigma_xx + sigma_yy) + hyp;
        sigma_2 = 0.5 * (sigma_xx + sigma_yy) - hyp;
        angle = 0.5 * std::atan2(2.0 * sigma_xy, sigma_xx - sigma_yy);
        cos_a = std::cos(angle);
        sin_a = std::sin(angle);
        stress_vector_pos(0) = 0.0;
        stress_vector_pos(1) = 0.0;
        stress_vector_pos(2) = 0.0;
        if (sigma_1 > 0)
        {
            stress_vector_pos(0) += sigma_1 * cos_a * cos_a;
            stress_vector_pos(1) += sigma_1 * sin_a * sin_a;
            stress_vector_pos(2) += sigma_1 * sin_a * cos_a;
        }
        if (sigma_2 > 0)
        {
            stress_vector_pos(0) += sigma_2 * sin_a * sin_a;
            stress_vector_pos(1) += sigma_2 * cos_a * cos_a;
            stress_vector_pos(2) -= sigma_2 * sin_a * cos_a;
        }
    }

    const double strain_norm = std::sqrt(inner_prod(stress_vector_pos, strain_vector));

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
        const double damage_rate = (stress_variable - hardening_modulus * mStrainVariable) / (mStrainVariable * mStrainVariable * mStrainVariable);
        constitutive_matrix *= (1. - damage_variable);
        constitutive_matrix -= damage_rate * outer_prod(stress_vector_pos, stress_vector);
        stress_vector *= (1. - damage_variable);
    }
    //TODO: Add e.g. COMPUTE_STRAIN_ENERGY flag
    //mStrainEnergy = 0.5 * ((1. - damage_variable) * inner_prod(strain_vector, prod(constitutive_elastic_matrix, strain_vector)));
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

double LinearIsotropicDamagePlaneStrain2D::EvaluateHardeningLaw(double r, const Properties &material_prop)
{
    double H = material_prop[ISOTROPIC_HARDENING_MODULUS];
    double r0 = material_prop[YIELD_STRESS] / std::sqrt(material_prop[YOUNG_MODULUS]);
    double q_inf = material_prop[INFINITY_YIELD_STRESS] /
                   std::sqrt(material_prop[YOUNG_MODULUS]);
    double q;

    if (r < r0)
        return r;
    q = r0 + H * (r - r0);
    if ((H > 0 && q > q_inf) || (H < 0 && q < q_inf))
        q = q_inf;
    return q;
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::CalculateConstitutiveMatrix(Matrix &D, const Properties &props)
{
    double E = props[YOUNG_MODULUS];
    double nu = props[POISSON_RATIO];
    double Ebar = E / (1. - nu * nu);
    double nubar = nu / (1. - nu);

    D.clear();

    D(0, 0) = 1;
    D(0, 1) = nubar;
    D(0, 2) = 0;
    D(1, 0) = nubar;
    D(1, 1) = 1;
    D(1, 2) = 0;
    D(2, 0) = 0;
    D(2, 1) = 0;
    D(2, 2) = 0.5 * (1 - nubar);

    D *= Ebar / (1. - nubar * nubar);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = 3;
    rFeatures.mSpaceDimension = 2;
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LinearIsotropicDamage3D);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LinearIsotropicDamage3D);
}

} /* namespace Kratos.*/
