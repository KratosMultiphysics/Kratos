// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//
// System includes

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_advanced_constitutive/truss_plasticity_constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TrussPlasticityConstitutiveLaw::TrussPlasticityConstitutiveLaw()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

TrussPlasticityConstitutiveLaw::TrussPlasticityConstitutiveLaw(const TrussPlasticityConstitutiveLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer TrussPlasticityConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<TrussPlasticityConstitutiveLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

TrussPlasticityConstitutiveLaw::~TrussPlasticityConstitutiveLaw()
{
    // TODO: Add if necessary
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void TrussPlasticityConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = 1;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}

//************************************************************************************
//************************************************************************************

double& TrussPlasticityConstitutiveLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == PLASTIC_STRAIN) rValue = mAccumulatedPlasticStrain;
    else if(rThisVariable == PLASTIC_ALPHA) rValue = mPlasticAlpha;
    else KRATOS_ERROR << "Can't get the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

array_1d<double, 3 > & TrussPlasticityConstitutiveLaw::GetValue(
    const Variable<array_1d<double, 3 > >& rThisVariable,
    array_1d<double, 3 > & rValue)
{
    KRATOS_ERROR << "Can't get the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

void TrussPlasticityConstitutiveLaw::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rThisVariable == PLASTIC_STRAIN) mAccumulatedPlasticStrain = rValue;
    else if(rThisVariable == PLASTIC_ALPHA) mPlasticAlpha = rValue;
    else KRATOS_ERROR << "Can't set the specified value" << std::endl;
}

//************************************************************************************
//************************************************************************************

double& TrussPlasticityConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == TANGENT_MODULUS)
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();
        const Properties& r_material_properties = rParameterValues.GetMaterialProperties();
        const double hardening_modulus = r_material_properties[HARDENING_MODULUS_1D];
        const double youngs_modulus = r_material_properties[YOUNG_MODULUS];
        if (mCurrentInElasticFlag)
        {
            KRATOS_DEBUG_ERROR_IF((hardening_modulus+youngs_modulus)<numerical_limit)
             << "Dividing by 0 when calculating the plastic tangent modulus" << std::endl;
            rValue = (hardening_modulus*youngs_modulus)/(hardening_modulus+youngs_modulus);
        }
        else rValue = youngs_modulus;

    }
    else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;

    return rValue;
}


//************************************************************************************
//************************************************************************************
void TrussPlasticityConstitutiveLaw::CalculateMaterialResponsePK2Custom(Parameters& rValues, double& rCurrentAccumulatedPlasticStrain, double& rCurrentPlasticAlpha)
{
    KRATOS_ERROR_IF_NOT(rValues.IsSetStrainVector()) << "Strain vector not set" << std::endl;
    KRATOS_DEBUG_ERROR_IF_NOT(rValues.IsSetStressVector()) << "Stress vector not set" << std::endl;

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const double prestress = r_material_properties[TRUSS_PRESTRESS_PK2];


    const double axial_strain (rValues.GetStrainVector()[0]);
    rCurrentAccumulatedPlasticStrain = mAccumulatedPlasticStrain;
    rCurrentPlasticAlpha = mPlasticAlpha;
    const double elastic_trial_strain = axial_strain-rCurrentAccumulatedPlasticStrain;
    const double youngs_modulus = r_material_properties[YOUNG_MODULUS];
    double current_stress  = youngs_modulus*elastic_trial_strain;
    // consider pre-stress in plastic algorithm
    current_stress += prestress;
    double temp_stress(current_stress);

    // start plastic algorithm
    mCurrentInElasticFlag = CheckIfIsPlasticRegime(rValues,current_stress);
    if (mCurrentInElasticFlag)
    {
        const double hardening_modulus = r_material_properties[HARDENING_MODULUS_1D];
        const double youngs_modulus = r_material_properties[YOUNG_MODULUS];
        const double trial_yield_function = TrialYieldFunction(r_material_properties,current_stress);

        const double delta_gamma = trial_yield_function/(youngs_modulus+hardening_modulus);
        current_stress = 1.00 - ((delta_gamma*youngs_modulus)/std::abs(temp_stress));
        current_stress = current_stress * temp_stress;

        rCurrentAccumulatedPlasticStrain += delta_gamma*MathUtils<double>::Sign(temp_stress);
        rCurrentPlasticAlpha += delta_gamma;
    }

    // return only material response
    Vector& stress_vector = rValues.GetStressVector();
    stress_vector = ZeroVector(1);
    stress_vector[0] = current_stress-prestress;
}


void TrussPlasticityConstitutiveLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    double temp_accumulated_plastic_strain;
    double temp_plastic_alpha;
    CalculateMaterialResponsePK2Custom(rValues,temp_accumulated_plastic_strain,temp_plastic_alpha);
}


//************************************************************************************
//************************************************************************************

double TrussPlasticityConstitutiveLaw::TrialYieldFunction(const Properties& rMaterialProperties,
     const double& rCurrentStress)
{
    const double yield_stress_limit = rMaterialProperties[YIELD_STRESS];
    const double hardening_modulus = rMaterialProperties[HARDENING_MODULUS_1D];

    double trial_yield_function =  std::abs(rCurrentStress);
    trial_yield_function -= yield_stress_limit + (hardening_modulus*mPlasticAlpha);

    return trial_yield_function;
}


//************************************************************************************
//************************************************************************************

bool TrussPlasticityConstitutiveLaw::CheckIfIsPlasticRegime(Parameters& rValues,const double& rCurrentStress)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    bool is_in_plastic_regime = false;
    const double trial_yield_function = TrialYieldFunction(rValues.GetMaterialProperties(),rCurrentStress);
    if (trial_yield_function > numerical_limit) is_in_plastic_regime = true;

    const double check_limits = (rCurrentStress/rValues.GetMaterialProperties()[YOUNG_MODULUS])+mAccumulatedPlasticStrain;
    if (std::abs(check_limits)<numerical_limit) is_in_plastic_regime=false;

    return is_in_plastic_regime;
}

//************************************************************************************
//************************************************************************************

void TrussPlasticityConstitutiveLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    double temp_accumulated_plastic_strain;
    double temp_plastic_alpha;
    CalculateMaterialResponsePK2Custom(rValues,temp_accumulated_plastic_strain,temp_plastic_alpha);
    mAccumulatedPlasticStrain = temp_accumulated_plastic_strain;
    mPlasticAlpha = temp_plastic_alpha;
}

int TrussPlasticityConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));

    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK(rMaterialProperties.Has(DENSITY));

    KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS);
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS));

    KRATOS_CHECK_VARIABLE_KEY(HARDENING_MODULUS_1D);
    KRATOS_CHECK(rMaterialProperties.Has(HARDENING_MODULUS_1D));
    return 0;
}
} // Namespace Kratos
