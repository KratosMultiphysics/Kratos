// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi
//                   Manuel Caicedo
//                   Alfredo Huespe

#include "linear_j2_plasticity_plane_strain_2d.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearJ2PlasticityPlaneStrain2D::LinearJ2PlasticityPlaneStrain2D()
    : ConstitutiveLaw()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

LinearJ2PlasticityPlaneStrain2D::LinearJ2PlasticityPlaneStrain2D(const LinearJ2PlasticityPlaneStrain2D &rOther)
            : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearJ2PlasticityPlaneStrain2D::Clone() const
{
    return ConstitutiveLaw::Pointer(new LinearJ2PlasticityPlaneStrain2D());
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

LinearJ2PlasticityPlaneStrain2D::~LinearJ2PlasticityPlaneStrain2D()
{
}

//************************************************************************************
//************************************************************************************

bool LinearJ2PlasticityPlaneStrain2D::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        return true;
    }
    if(rThisVariable == INELASTIC_FLAG){
        return true;
    }
    return false;
}

//************************************************************************************
//************************************************************************************
/*
LinearJ2PlasticityPlaneStrain2D::SizeType LinearJ2PlasticityPlaneStrain2D::WorkingSpaceDimension()
{
    return 2;
}

LinearJ2PlasticityPlaneStrain2D::SizeType LinearJ2PlasticityPlaneStrain2D::GetStrainSize()
{
    return 4;
}
*/

double& LinearJ2PlasticityPlaneStrain2D::GetValue(const Variable<double>& rThisVariable,
                                                         double& rValue)
{
    if(rThisVariable == INELASTIC_FLAG){
        rValue = mInelasticFlag;
    }
    return rValue;
}
/*
LinearJ2PlasticityPlaneStrain2D::StrainMeasure LinearJ2PlasticityPlaneStrain2D::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

LinearJ2PlasticityPlaneStrain2D::StressMeasure LinearJ2PlasticityPlaneStrain2D::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}
*/

void LinearJ2PlasticityPlaneStrain2D::InitializeMaterial(const Properties& material_prop,
                                                                const GeometryType& rElementGeometry,
                                                                const Vector& rShapeFunctionsValues)
{
    mPlasticStrainOld = ZeroVector(this->GetStrainSize());
    mAccumulatedPlasticStrainOld = 0.;
}

void LinearJ2PlasticityPlaneStrain2D::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    mPlasticStrainOld = mPlasticStrain;
    mAccumulatedPlasticStrainOld = mAccumulatedPlasticStrain;
}

void LinearJ2PlasticityPlaneStrain2D::CalculateMaterialResponsePK1(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearJ2PlasticityPlaneStrain2D::CalculateMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearJ2PlasticityPlaneStrain2D::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearJ2PlasticityPlaneStrain2D::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const Properties& matprops = rValues.GetMaterialProperties();
    Vector& epsilon = rValues.GetStrainVector();
    Vector& StressVector = rValues.GetStressVector();
    Matrix ElasticityTensor;
    Matrix& TangentTensor =
        rValues.GetConstitutiveMatrix(); // TODO find proper getter
    double hardening_modulus = matprops[ISOTROPIC_HARDENING_MODULUS];
    double delta_k = matprops[INFINITY_HARDENING_MODULUS];
    double hardening_exponent = matprops[HARDENING_EXPONENT];
    double trial_yield_function;
    double E = matprops[YOUNG_MODULUS];
    double poisson_ratio = matprops[POISSON_RATIO];
    double mu = E / (2. + 2. * poisson_ratio);
    double volumetric_modulus = E / (3 * (1 - 2 * poisson_ratio));

    if (rValues.GetProcessInfo().Has(INITIAL_STRAIN))
    {
        noalias(epsilon) += rValues.GetProcessInfo()[INITIAL_STRAIN];
    }

    mPlasticStrain = mPlasticStrainOld;
    mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld;

    ElasticityTensor.resize(4, 4);
    CalculateElasticityTensor(matprops, ElasticityTensor);
    Vector sigma_trial;
    sigma_trial.resize(4);
    sigma_trial = prod(ElasticityTensor, epsilon - mPlasticStrainOld);

    // StressTrialDev = sigma - 1/3 tr(sigma) * I
    Vector StressTrialDev;
    StressTrialDev.resize(4);
    StressTrialDev = sigma_trial;

    double trace =
        0.33333333333333333 * (sigma_trial(0) + sigma_trial(1) + sigma_trial(2));
    StressTrialDev(0) -= trace;
    StressTrialDev(1) -= trace;
    StressTrialDev(2) -= trace;
    double norm_dev_stress = std::sqrt(StressTrialDev(0) * StressTrialDev(0) +
                                       StressTrialDev(1) * StressTrialDev(1) +
                                       StressTrialDev(2) * StressTrialDev(2) +
                                       2. * StressTrialDev(3) * StressTrialDev(3));
    trial_yield_function = this->yieldFunction(norm_dev_stress, matprops);
    double dgamma = 0;

    //KRATOS_WATCH(epsilon)

    if (trial_yield_function <= 0.)
    {
        // ELASTIC
        StressVector = sigma_trial;
        TangentTensor = ElasticityTensor;
    }
    else
    {
        // INELASTIC
        mInelasticFlag = 1;

        Vector YieldFunctionNormalVector = StressTrialDev / norm_dev_stress;
        if (delta_k != 0.0 && hardening_exponent != 0.0)
        {
            // Exponential softening
            dgamma = GetDeltaGamma(norm_dev_stress, matprops);
        }
        else
        {
            // Linear softening
            dgamma = trial_yield_function /
                     (2. * mu * (1. + (hardening_modulus / (3. * mu))));
        }

        StressVector(0) =
            volumetric_modulus * (epsilon(0) + epsilon(1) + epsilon(2)) +
            StressTrialDev(0) - 2. * mu * dgamma * YieldFunctionNormalVector(0);
        StressVector(1) =
            volumetric_modulus * (epsilon(0) + epsilon(1) + epsilon(2)) +
            StressTrialDev(1) - 2. * mu * dgamma * YieldFunctionNormalVector(1);
        StressVector(2) =
            volumetric_modulus * (epsilon(0) + epsilon(1) + epsilon(2)) +
            StressTrialDev(2) - 2. * mu * dgamma * YieldFunctionNormalVector(2);
        StressVector(3) =
            StressTrialDev(3) - 2. * mu * dgamma * YieldFunctionNormalVector(3);

        mPlasticStrain(0) = mPlasticStrainOld(0) + dgamma * YieldFunctionNormalVector(0);
        mPlasticStrain(1) = mPlasticStrainOld(1) + dgamma * YieldFunctionNormalVector(1);
        mPlasticStrain(2) = mPlasticStrainOld(2) + dgamma * YieldFunctionNormalVector(2);
        mPlasticStrain(3) =
            mPlasticStrainOld(3) + dgamma * YieldFunctionNormalVector(3) * 2;

        mAccumulatedPlasticStrain =
            mAccumulatedPlasticStrainOld + 0.8164965809277260 * dgamma;

        // Update derivative of the hardening-softening modulus

        CalculateTangentTensor(dgamma, norm_dev_stress, YieldFunctionNormalVector,
                               matprops, TangentTensor);


    }

    // Linear + exponential hardening
    mStrainEnergy = 0.5 * inner_prod(epsilon - mPlasticStrain, prod(ElasticityTensor, epsilon - mPlasticStrain))
                    + GetPlasticPotential(matprops);
}

double& LinearJ2PlasticityPlaneStrain2D::CalculateValue(Parameters& rParameterValues,
        const Variable<double>& rThisVariable, double& rValue) {

    //const Properties& MaterialProperties = rParameterValues.GetMaterialProperties();
    //Vector& StrainVector = rParameterValues.GetStrainVector();
    //Vector& StressVector = rParameterValues.GetStressVector();
    //const double& E = MaterialProperties[YOUNG_MODULUS];
    //const double& NU = MaterialProperties[POISSON_RATIO];

    if(rThisVariable == STRAIN_ENERGY){
        //TODO (marcelo): calculate STRAIN_ENERGY here and remove mStrainEnergy
        //CalculateCauchyGreenStrain(rParameterValues, StrainVector);
        //CalculatePK2Stress( StrainVector, StressVector, E, NU );
        //rValue = 0.5 * inner_prod(StrainVector,StressVector); // Strain energy = 0.5*E:C:E
        rValue = mStrainEnergy;
    }
    return(rValue);
}

void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
}

double LinearJ2PlasticityPlaneStrain2D::GetSaturationHardening(const Properties& rMaterialProperties)
{
    double yield_stress = rMaterialProperties[YIELD_STRESS];
    double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];

    double k_new = yield_stress + (theta * hardening_modulus * mAccumulatedPlasticStrain) +
               delta_k * (1. - std::exp(-hardening_exponent * mAccumulatedPlasticStrain));
    return k_new;
}

double LinearJ2PlasticityPlaneStrain2D::GetPlasticPotential(const Properties& rMaterialProperties)
{
   double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
   double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
   double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
   double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];

   double Wp_new = 0.5*(theta * hardening_modulus * std::pow(mAccumulatedPlasticStrain, 2.0)) +
                   delta_k * (mAccumulatedPlasticStrain - (1/hardening_exponent) * (1- std::exp(-hardening_exponent * mAccumulatedPlasticStrain)));
   return Wp_new;
}

double LinearJ2PlasticityPlaneStrain2D::GetDeltaGamma(double norm_s_trial,
                                                             const Properties& rMaterialProperties)
{
    double E = rMaterialProperties[YOUNG_MODULUS];
    double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    double yield_stress = rMaterialProperties[YIELD_STRESS];
    double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    double tolerance = 1e-6 * yield_stress;
    double mu = E / (2. * (1. + poisson_ratio));
    double dgamma = 0.0;
    double norm_yieldfunction = 1.0;

    while (norm_yieldfunction > tolerance)
    {
        double k_new = GetSaturationHardening(rMaterialProperties);
        double kp_new =
            theta * hardening_modulus +
            delta_k * (hardening_exponent *
                       std::exp(-hardening_exponent * mAccumulatedPlasticStrain));
        double yieldfunction = -0.8164965809277260 * k_new + norm_s_trial - 2. * mu * dgamma;
        double derivative_yieldfunction = -2. * mu * (1. + kp_new / (3. * mu));
        dgamma = dgamma - yieldfunction / derivative_yieldfunction;
        mAccumulatedPlasticStrain =
            mAccumulatedPlasticStrainOld + 0.8164965809277260 * dgamma;
        norm_yieldfunction = std::abs(yieldfunction);
    }
    // TODO (marcelo): handle the case when no convergence is achieved.
    return dgamma;
}

double LinearJ2PlasticityPlaneStrain2D::yieldFunction(const double norm_dev_stress,
                                                             const Properties& rMaterialProperties)
{
    double yield_stress = rMaterialProperties[YIELD_STRESS];
    double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    double k_old =
        yield_stress + (theta * hardening_modulus * mAccumulatedPlasticStrainOld) +
        (delta_k) * (1. - std::exp(-hardening_exponent * mAccumulatedPlasticStrainOld));

    return norm_dev_stress - k_old * 0.8164965809277260; // sqrt(2/3)
}

void LinearJ2PlasticityPlaneStrain2D::CalculateElasticityTensor(const Properties& props,
                                                                       Matrix& D)
{
    double E = props[YOUNG_MODULUS];
    double poisson_ratio = props[POISSON_RATIO];
    double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    double mu = E / (2. + 2. * poisson_ratio);

    D.clear();

    D(0, 0) = lambda + 2. * mu;
    D(0, 1) = lambda;
    D(0, 2) = lambda;
    D(0, 3) = 0.;
    D(1, 0) = lambda;
    D(1, 1) = lambda + 2. * mu;
    D(1, 2) = lambda;
    D(1, 3) = 0.;
    D(2, 0) = lambda;
    D(2, 1) = lambda;
    D(2, 2) = lambda + 2. * mu;
    D(2, 3) = 0.;
    D(3, 0) = 0.;
    D(3, 1) = 0.;
    D(3, 2) = 0.;
    D(3, 3) = mu;
}

void LinearJ2PlasticityPlaneStrain2D::CalculateTangentTensor(
    double dgamma, double norm_s_trial, const Vector& N_new, const Properties& props, Matrix& D)
{
    double hardening_modulus = props[ISOTROPIC_HARDENING_MODULUS];
    double theta = props[REFERENCE_HARDENING_MODULUS];
    double delta_k = props[INFINITY_HARDENING_MODULUS];
    double hardening_exponent = props[HARDENING_EXPONENT];
    double E = props[YOUNG_MODULUS];
    double poisson_ratio = props[POISSON_RATIO];
    double mu = E / (2. + 2. * poisson_ratio);
    double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));

    double kp_new = (theta * hardening_modulus) +
                    delta_k * (hardening_exponent *
                               std::exp(-hardening_exponent * mAccumulatedPlasticStrain));

    double theta_new = 1 - (2. * mu * dgamma) / norm_s_trial;
    double theta_new_b = 1. / (1. + kp_new / (3. * mu)) - (1. - theta_new);

    D(0, 0) = volumetric_modulus + (2 * mu * theta_new * 2. / 3.) -
              (2 * mu * theta_new_b * (N_new(0) * N_new(0)));
    D(0, 1) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
              (2 * mu * theta_new_b * (N_new(0) * N_new(1)));
    D(0, 2) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
              (2 * mu * theta_new_b * (N_new(0) * N_new(2)));
    D(0, 3) = -(2 * mu * theta_new_b * (N_new(0) * N_new(3)));

    D(1, 0) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
              (2 * mu * theta_new_b * (N_new(1) * N_new(0)));
    D(1, 1) = volumetric_modulus + (2 * mu * theta_new * 2. / 3.) -
              (2 * mu * theta_new_b * (N_new(1) * N_new(1)));
    D(1, 2) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
              (2 * mu * theta_new_b * (N_new(1) * N_new(2)));
    D(1, 3) = -(2 * mu * theta_new_b * (N_new(1) * N_new(3)));

    D(2, 0) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
              (2 * mu * theta_new_b * (N_new(2) * N_new(0)));
    D(2, 1) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
              (2 * mu * theta_new_b * (N_new(2) * N_new(1)));
    D(2, 2) = volumetric_modulus + (2 * mu * theta_new * 2. / 3.) -
              (2 * mu * theta_new_b * (N_new(2) * N_new(2)));
    D(2, 3) = -(2 * mu * theta_new_b * (N_new(2) * N_new(3)));

    D(3, 0) = -(2 * mu * theta_new_b * (N_new(3) * N_new(0)));
    D(3, 1) = -(2 * mu * theta_new_b * (N_new(3) * N_new(1)));
    D(3, 2) = -(2 * mu * theta_new_b * (N_new(3) * N_new(2)));
    D(3, 3) = mu * theta_new - (2 * mu * theta_new_b * (N_new(3) * N_new(3)));
}

void LinearJ2PlasticityPlaneStrain2D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = 4;
    rFeatures.mSpaceDimension = 2;
}

int LinearJ2PlasticityPlaneStrain2D::Check(const Properties& rMaterialProperties,
                                                  const GeometryType& rElementGeometry,
                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
    KRATOS_CHECK_VARIABLE_KEY(POISSON_RATIO);
    KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS);
    KRATOS_CHECK_VARIABLE_KEY(REFERENCE_HARDENING_MODULUS);
    KRATOS_CHECK_VARIABLE_KEY(ISOTROPIC_HARDENING_MODULUS);
    KRATOS_CHECK_VARIABLE_KEY(INFINITY_HARDENING_MODULUS);
    KRATOS_CHECK_VARIABLE_KEY(HARDENING_EXPONENT);

    return 0;
}

} /* namespace Kratos.*/
