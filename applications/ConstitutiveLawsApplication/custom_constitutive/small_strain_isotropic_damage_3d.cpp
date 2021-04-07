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
#include "small_strain_isotropic_damage_3d.h"
#include "structural_mechanics_application_variables.h"
#include "includes/checks.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainIsotropicDamage3D::SmallStrainIsotropicDamage3D()
    : ElasticIsotropic3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

SmallStrainIsotropicDamage3D::SmallStrainIsotropicDamage3D(const SmallStrainIsotropicDamage3D &rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainIsotropicDamage3D::Clone() const
{
    return Kratos::make_shared<SmallStrainIsotropicDamage3D>(SmallStrainIsotropicDamage3D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SmallStrainIsotropicDamage3D::~SmallStrainIsotropicDamage3D()
{
}

//************************************************************************************
//************************************************************************************

bool SmallStrainIsotropicDamage3D::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        // explicitly returning "false", so we know we must call CalculateValue(...)
        return false;
    }
    if(rThisVariable == DAMAGE_VARIABLE){
        // explicitly returning "false", so we know we must call CalculateValue(...)
        return false;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

bool SmallStrainIsotropicDamage3D::Has(const Variable<Vector>& rThisVariable)
{
    if(rThisVariable == INTERNAL_VARIABLES){
        return true;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

Vector& SmallStrainIsotropicDamage3D::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if(rThisVariable == INTERNAL_VARIABLES){
        rValue.resize(1);
        rValue[0] = mStrainVariable;
    }

    return rValue;
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3D::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rProcessInfo
    )
{
    if(rThisVariable == INTERNAL_VARIABLES){
        mStrainVariable = rValue[0];
    }
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3D::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
    mStrainVariable = yield_stress / std::sqrt(young_modulus);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3D::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3D::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    Vector internal_variables(1);
    this->CalculateStressResponse(rValues, internal_variables);
    mStrainVariable = internal_variables[0];
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3D::CalculateMaterialResponsePK2(Parameters& rValues)
{
    Vector internal_variables(1);
    this->CalculateStressResponse(rValues, internal_variables);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3D::CalculateStressResponse(
        Parameters &rValues,
        Vector& rInternalVariables)
{
    double strain_variable = mStrainVariable;
    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    Flags & r_constitutive_law_options = rValues.GetOptions();
    Vector& r_strain_vector = rValues.GetStrainVector();
    if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
        noalias(r_strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
    }

    if( r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        //this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // If we compute the tangent moduli or the stress
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ||
        r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR )) {
        Vector& r_stress_vector = rValues.GetStressVector();
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rValues);
        noalias(r_stress_vector) = prod(r_constitutive_matrix, r_strain_vector);

        // Auxiliary stress vector to allow derived models (e.g. traction-only damage) to set
        // a different value of r_stress_vector_pos with ComputePositiveStressVector function
        // In the symmetric model, ComputePositiveStressVector does not do anything.
        Vector r_stress_vector_pos = r_stress_vector;
        ComputePositiveStressVector(r_stress_vector_pos, r_stress_vector);

        // energy may be a small negative due to machine precision error, forcing zero
        const double product = inner_prod(r_stress_vector_pos, r_strain_vector);
        const double strain_norm = (product >=0 ) ? std::sqrt(product) : 0;
        if (strain_norm <= mStrainVariable)
        {
            // ELASTIC
            strain_variable = mStrainVariable;
            const double stress_variable = EvaluateHardeningLaw(strain_variable, rMaterialProperties);
            const double damage_variable = 1. - stress_variable / strain_variable;
            r_constitutive_matrix *= (1 - damage_variable);
            r_stress_vector *= (1 - damage_variable);
        }
        else
        {
            // INELASTIC
            strain_variable = strain_norm;
            const double stress_variable = EvaluateHardeningLaw(strain_variable, rMaterialProperties);
            const double damage_variable = 1. - stress_variable / strain_variable;
            const double hardening_modulus = EvaluateHardeningModulus(strain_variable, rMaterialProperties);
            const double damage_rate = (stress_variable - hardening_modulus * strain_variable)
                                       / (strain_variable * strain_variable * strain_variable);
            r_constitutive_matrix *= (1. - damage_variable);
            r_constitutive_matrix -= damage_rate * outer_prod(r_stress_vector, r_stress_vector_pos);
            // Computing: real stress = (1-d) * effective stress
            r_stress_vector *= (1. - damage_variable);
        }
    }
    rInternalVariables[0] = strain_variable;
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3D::ComputePositiveStressVector(
            Vector& rStressVectorPos, Vector& rStressVector)
{
}

//************************************************************************************
//************************************************************************************

double SmallStrainIsotropicDamage3D::EvaluateHardeningModulus(
        double r,
        const Properties &rMaterialProperties)
{
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double inf_yield_stress = rMaterialProperties[INFINITY_YIELD_STRESS];
    const double young_modulus = rMaterialProperties[YOUNG_MODULUS];

    const double H0 = rMaterialProperties[HARDENING_MODULI_VECTOR][0];
    const double H1 = rMaterialProperties[HARDENING_MODULI_VECTOR][1];

    const double r0 = yield_stress / std::sqrt(young_modulus);
    const double q0 = r0;  // strain_variable_init
    const double q1 = inf_yield_stress / std::sqrt(young_modulus);  // stress_variable_inf
    const double r1 = r0 + (q1 - q0) / H0;

     if (r < r0)
        return 0.;
    if (r >= r0 && r < r1)
        return H0;
    //Case r >= r1:
    return H1;
}

//************************************************************************************
//************************************************************************************

double SmallStrainIsotropicDamage3D::EvaluateHardeningLaw(
        double r,
        const Properties &rMaterialProperties)
{
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double inf_yield_stress = rMaterialProperties[INFINITY_YIELD_STRESS];
    const double young_modulus = rMaterialProperties[YOUNG_MODULUS];

    const double r0 = yield_stress / std::sqrt(young_modulus);
    const double H0 = EvaluateHardeningModulus(r0, rMaterialProperties);
    const double q0 = r0;  // strain_variable_init
    const double q1 = inf_yield_stress / std::sqrt(young_modulus);  // stress_variable_inf
    const double r1 = r0 + (q1 - q0) / H0;
    const double H1 = EvaluateHardeningModulus(r1, rMaterialProperties);

    if (r < r0)
        return q0;
    if (r >= r0 && r < r1)
        return q0 + H0 * (r - r0);
    // Case r >= r1:
    return q1 + H1 * (r - r1);
}

//************************************************************************************
//************************************************************************************

double& SmallStrainIsotropicDamage3D::CalculateValue(
    Parameters& rValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == STRAIN_ENERGY){
        Vector& r_strain_vector = rValues.GetStrainVector();
        if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(r_strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
        }
        const Properties& rMaterialProperties = rValues.GetMaterialProperties();
        Matrix constitutive_matrix;
        CalculateElasticMatrix(constitutive_matrix, rValues);
        const double stress_like_variable = EvaluateHardeningLaw(mStrainVariable, rMaterialProperties);
        const double damage_variable = 1. - stress_like_variable / mStrainVariable;

        rValue = 0.5 * ((1. - damage_variable) * inner_prod(r_strain_vector,
                                              prod(constitutive_matrix, r_strain_vector)));
    }

    if (rThisVariable == DAMAGE_VARIABLE){
        const Properties& rMaterialProperties = rValues.GetMaterialProperties();
        const double stress_like_variable = EvaluateHardeningLaw(mStrainVariable, rMaterialProperties);

        rValue = 1. - stress_like_variable / mStrainVariable;
    }

    return(rValue);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3D::GetLawFeatures(Features& rFeatures)
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

int SmallStrainIsotropicDamage3D::Check(
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
    KRATOS_CHECK(rMaterialProperties.Has(HARDENING_MODULI_VECTOR));
    KRATOS_CHECK_GREATER(rMaterialProperties[YIELD_STRESS], tolerance);
    KRATOS_CHECK_GREATER(rMaterialProperties[INFINITY_YIELD_STRESS], tolerance);
    KRATOS_CHECK_LESS_EQUAL(rMaterialProperties[HARDENING_MODULI_VECTOR](0), 1.);
    KRATOS_CHECK_GREATER_EQUAL(rMaterialProperties[HARDENING_MODULI_VECTOR](0), 0.);
    KRATOS_CHECK_LESS_EQUAL(rMaterialProperties[HARDENING_MODULI_VECTOR](1), 1.);
    KRATOS_CHECK_GREATER_EQUAL(rMaterialProperties[HARDENING_MODULI_VECTOR](1), 0.);

    return 0;
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mStrainVariable", mStrainVariable);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamage3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mStrainVariable", mStrainVariable);
}

} /* namespace Kratos.*/
