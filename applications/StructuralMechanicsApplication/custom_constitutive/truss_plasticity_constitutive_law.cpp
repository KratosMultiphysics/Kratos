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
#include "custom_constitutive/truss_plasticity_constitutive_law.h"
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

bool& TrussPlasticityConstitutiveLaw::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    if(rThisVariable == INELASTIC_FLAG) rValue = this->mInElasticFlagVector[0];
    else KRATOS_ERROR << "Can't get the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

double& TrussPlasticityConstitutiveLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == PLASTIC_STRAIN) rValue = this->mAccumulatedPlasticStrainVector[0];
    else if(rThisVariable == PLASTIC_ALPHA) rValue = this->mPlasticAlphaVector[0];
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
    if(rThisVariable == PLASTIC_STRAIN) this->mAccumulatedPlasticStrainVector[0] = rValue;
    else if(rThisVariable == PLASTIC_ALPHA) this->mPlasticAlphaVector[0] = rValue;
    else KRATOS_ERROR << "Can't set the specified value" << std::endl;
}

//************************************************************************************
//************************************************************************************

void TrussPlasticityConstitutiveLaw::SetValue(const Variable<bool>& rThisVariable,
                          const bool& rValue,
                          const ProcessInfo& rCurrentProcessInfo)
{
    if(rThisVariable == INELASTIC_FLAG) this->mInElasticFlagVector[0] = rValue;
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
        if (this->mInElasticFlagVector[0])
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

Vector& TrussPlasticityConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{

    if(rThisVariable == NORMAL_STRESS)
    {
        constexpr SizeType dofs = 6;
        const Properties& r_material_properties = rParameterValues.GetMaterialProperties();

        rValue = ZeroVector(dofs);
        double current_stress = this->mStressState;


        if (this->mInElasticFlagVector[0])
        {
            const double hardening_modulus = r_material_properties[HARDENING_MODULUS_1D];
            const double youngs_modulus = r_material_properties[YOUNG_MODULUS];
            const double trial_yield_function = this->TrialYieldFunction(r_material_properties,current_stress);


            const double delta_gamma = trial_yield_function/(youngs_modulus+hardening_modulus);
            current_stress = 1.00 - ((delta_gamma*youngs_modulus)/std::abs(this->mStressState));
            current_stress = current_stress * this->mStressState;

            this->mAccumulatedPlasticStrainVector[0] = this->mAccumulatedPlasticStrainVector[0] + (delta_gamma*MathUtils<double>::Sign(this->mStressState));
            this->mPlasticAlphaVector[0] = this->mPlasticAlphaVector[0] + delta_gamma;
        }

        rValue[0] = -1.0 * current_stress;
        rValue[3] = 1.0 * current_stress;

        this->mStressState = current_stress;

    }
    else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

array_1d<double, 3 > & TrussPlasticityConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<array_1d<double, 3 > >& rVariable,
	array_1d<double, 3 > & rValue)
    {
        if (rVariable == FORCE)
        {
            constexpr SizeType dimension = 3;
            rValue = ZeroVector(dimension);
            rValue[0] = this->mStressState;
            rValue[1] = 0.0;
            rValue[2] = 0.0;
        }
        else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
        return rValue;
    }


//************************************************************************************
//************************************************************************************
// calculate the trial stress !
void TrussPlasticityConstitutiveLaw::CalculateMaterialResponse(
    const Vector& rStrainVector,const Matrix& rDeformationGradient,
    Vector& rStressVector,Matrix& rAlgorithmicTangent,
    const ProcessInfo& rCurrentProcessInfo,const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues,
    bool CalculateStresses,int CalculateTangent,bool SaveInternalVariables)
{
    const double axial_strain = rStrainVector[0];
    double current_plastic_strain = 0.00;
    this->GetValue(PLASTIC_STRAIN,current_plastic_strain);

    const double elastic_trial_strain = axial_strain-current_plastic_strain;
    const double youngs_modulus = rMaterialProperties[YOUNG_MODULUS];

    if (rStressVector.size() != 1) rStressVector.resize(1);
    rStressVector[0] = youngs_modulus*elastic_trial_strain;

    if (SaveInternalVariables)
    {
        this->mStressState = rStressVector[0];
        this->mInElasticFlagVector[0] = this->CheckIfIsPlasticRegime(rMaterialProperties,rStressVector[0]);

        if (this->CheckPlasticIterationHistory())
        {
            this->mAccumulatedPlasticStrainVector[0] = this->mAccumulatedPlasticStrainVector[1];
            this->mPlasticAlphaVector[0] = this->mPlasticAlphaVector[1];
        }
    }

}


//************************************************************************************
//************************************************************************************

double TrussPlasticityConstitutiveLaw::TrialYieldFunction(const Properties& rMaterialProperties,
     const double& rCurrentStress)
{
    const double yield_stress_limit = rMaterialProperties[YIELD_STRESS];
    const double hardening_modulus = rMaterialProperties[HARDENING_MODULUS_1D];

    double trial_yield_function =  std::abs(rCurrentStress);
    trial_yield_function -= yield_stress_limit + (hardening_modulus*this->mPlasticAlphaVector[0]);

    return trial_yield_function;
}


//************************************************************************************
//************************************************************************************

bool TrussPlasticityConstitutiveLaw::CheckIfIsPlasticRegime(const Properties& rMaterialProperties,
     const double& rCurrentStress)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();

    bool is_in_plastic_regime = false;
    const double trial_yield_function = this->TrialYieldFunction(rMaterialProperties,rCurrentStress);
    if (trial_yield_function > 0.00) is_in_plastic_regime = true;


    const double check_limits = (rCurrentStress/rMaterialProperties[YOUNG_MODULUS])+this->mAccumulatedPlasticStrainVector[0];
    if (std::abs(check_limits)<numerical_limit) is_in_plastic_regime=false;

    return is_in_plastic_regime;
}

//************************************************************************************
//************************************************************************************




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


void TrussPlasticityConstitutiveLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                    const GeometryType& rElementGeometry,
                    const Vector& rShapeFunctionsValues,
                    const ProcessInfo& rCurrentProcessInfo)
{
    this->mInElasticFlagVector[1] = this->mInElasticFlagVector[0];
}

void TrussPlasticityConstitutiveLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
                    const GeometryType& rElementGeometry,
                    const Vector& rShapeFunctionsValues,
                    const ProcessInfo& rCurrentProcessInfo)
{
    this->mAccumulatedPlasticStrainVector[1] = this->mAccumulatedPlasticStrainVector[0];
    this->mPlasticAlphaVector[1] = this->mPlasticAlphaVector[0];
    this->mInElasticFlagVector[0] = false;
    this->mInElasticFlagVector[1] = false;
}

} // Namespace Kratos
