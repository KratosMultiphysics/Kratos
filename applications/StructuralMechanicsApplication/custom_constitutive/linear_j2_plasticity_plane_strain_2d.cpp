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
    LinearJ2PlasticityPlaneStrain2D::Pointer pclone(new LinearJ2PlasticityPlaneStrain2D(*this));
    return pclone;
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

double& LinearJ2PlasticityPlaneStrain2D::GetValue(const Variable<double>& rThisVariable,
                                                  double& rValue)
{
    if(rThisVariable == INELASTIC_FLAG){
        rValue = mInelasticFlag;
    }
    return rValue;
}

void LinearJ2PlasticityPlaneStrain2D::InitializeMaterial(const Properties& material_prop,
                                                         const GeometryType& rElementGeometry,
                                                         const Vector& rShapeFunctionsValues)
{
    mPlasticStrainOld = ZeroVector(this->GetStrainSize());
    mAccumulatedPlasticStrainOld = 0.0;
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
    Flags &Options = rValues.GetOptions();

    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    Vector& strain_vector = rValues.GetStrainVector();
    Vector& stress_vector = rValues.GetStressVector();

    if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(ConstitutiveMatrix, rMaterialProperties);
    }

    if (Options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
        }
        Matrix elastic_tensor;
        Matrix& tangent_tensor = rValues.GetConstitutiveMatrix();
        const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
        const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
        const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
        const double E = rMaterialProperties[YOUNG_MODULUS];
        const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
        const double mu = E / (2. + 2. * poisson_ratio);
        const double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));
        const double sqrt_two_thirds = std::sqrt(2.0 / 3.0); // =0.8164965809277260
        double trial_yield_function;

        mPlasticStrain = mPlasticStrainOld;
        mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld;

        elastic_tensor.resize(4, 4, false);
        CalculateElasticMatrix(elastic_tensor, rMaterialProperties);
        Vector sigma_trial(4);
        noalias(sigma_trial) = prod(elastic_tensor, strain_vector - mPlasticStrainOld);

        // stress_trial_dev = sigma - 1/3 tr(sigma) * I
        Vector stress_trial_dev = sigma_trial;

        const double trace = 1.0 / 3.0 * (sigma_trial(0) + sigma_trial(1) + sigma_trial(2));
        stress_trial_dev(0) -= trace;
        stress_trial_dev(1) -= trace;
        stress_trial_dev(2) -= trace;
        const double norm_dev_stress = std::sqrt(stress_trial_dev(0) * stress_trial_dev(0) +
                                           stress_trial_dev(1) * stress_trial_dev(1) +
                                           stress_trial_dev(2) * stress_trial_dev(2) +
                                           2. * stress_trial_dev(3) * stress_trial_dev(3));
        trial_yield_function = this->yieldFunction(norm_dev_stress, rMaterialProperties);

        if (trial_yield_function <= 0.) {
            // ELASTIC
            mInelasticFlag = false;
            stress_vector = sigma_trial;
            tangent_tensor = elastic_tensor;
        }
        else {
            // INELASTIC
            mInelasticFlag = true;
            double dgamma = 0;
            Vector yield_function_normal_vector = stress_trial_dev / norm_dev_stress;
            if (delta_k != 0.0 && hardening_exponent != 0.0) {
                // Exponential softening
                dgamma = GetDeltaGamma(norm_dev_stress, rMaterialProperties);
            }
            else {
                // Linear softening
                dgamma = trial_yield_function /
                         (2. * mu * (1. + (hardening_modulus / (3. * mu))));
            }

            stress_vector(0) =
                volumetric_modulus * (strain_vector(0) + strain_vector(1) + strain_vector(2)) +
                stress_trial_dev(0) - 2. * mu * dgamma * yield_function_normal_vector(0);
            stress_vector(1) =
                volumetric_modulus * (strain_vector(0) + strain_vector(1) + strain_vector(2)) +
                stress_trial_dev(1) - 2. * mu * dgamma * yield_function_normal_vector(1);
            stress_vector(2) =
                volumetric_modulus * (strain_vector(0) + strain_vector(1) + strain_vector(2)) +
                stress_trial_dev(2) - 2. * mu * dgamma * yield_function_normal_vector(2);
            stress_vector(3) =
                stress_trial_dev(3) - 2. * mu * dgamma * yield_function_normal_vector(3);

            mPlasticStrain(0) = mPlasticStrainOld(0) + dgamma * yield_function_normal_vector(0);
            mPlasticStrain(1) = mPlasticStrainOld(1) + dgamma * yield_function_normal_vector(1);
            mPlasticStrain(2) = mPlasticStrainOld(2) + dgamma * yield_function_normal_vector(2);
            mPlasticStrain(3) = mPlasticStrainOld(3) + dgamma * yield_function_normal_vector(3) * 2;
            mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld + sqrt_two_thirds * dgamma;

            // Update derivative of the hardening-softening modulus

            CalculateTangentTensor(dgamma, norm_dev_stress, yield_function_normal_vector,
                                   rMaterialProperties, tangent_tensor);
        }

    // Linear + exponential hardening
    mStrainEnergy = 0.5 * inner_prod(strain_vector - mPlasticStrain, prod(elastic_tensor, strain_vector - mPlasticStrain))
                    + GetPlasticPotential(rMaterialProperties);
    }
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
}

void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponsePK2(Parameters& rValues)
{
}

void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
}

void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
}

double LinearJ2PlasticityPlaneStrain2D::GetSaturationHardening(const Properties& rMaterialProperties)
{
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];

    const double k_new = yield_stress + (theta * hardening_modulus * mAccumulatedPlasticStrain) +
               delta_k * (1. - std::exp(-hardening_exponent * mAccumulatedPlasticStrain));
    return k_new;
}

double LinearJ2PlasticityPlaneStrain2D::GetPlasticPotential(const Properties& rMaterialProperties)
{
   const double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
   const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
   const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
   const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];

   const double wp_new = 0.5*(theta * hardening_modulus * std::pow(mAccumulatedPlasticStrain, 2.0)) +
                   delta_k * (mAccumulatedPlasticStrain -
                   (1/hardening_exponent) * (1- std::exp(-hardening_exponent * mAccumulatedPlasticStrain)));
   return wp_new;
}

double LinearJ2PlasticityPlaneStrain2D::GetDeltaGamma(double NormSTrial,
                                                      const Properties& rMaterialProperties)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    const double tolerance = 1e-6 * yield_stress;
    const double mu = E / (2. * (1. + poisson_ratio));
    const double sqrt_two_thirds = std::sqrt(2.0 / 3.0); // =0.8164965809277260
    double dgamma = 0.0;
    double norm_yieldfunction = 1.0;

    while (norm_yieldfunction > tolerance)
    {
        const double k_new = GetSaturationHardening(rMaterialProperties);
        const double kp_new = theta * hardening_modulus +
            delta_k * (hardening_exponent * std::exp(-hardening_exponent * mAccumulatedPlasticStrain));
        const double yieldfunction = - sqrt_two_thirds * k_new + NormSTrial - 2. * mu * dgamma;
        const double derivative_yieldfunction = -2. * mu * (1. + kp_new / (3. * mu));
        dgamma = dgamma - yieldfunction / derivative_yieldfunction;
        mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld + sqrt_two_thirds * dgamma;
        norm_yieldfunction = std::abs(yieldfunction);
    }
    // TODO (marcelo): handle the case when no convergence is achieved.
    return dgamma;
}

double LinearJ2PlasticityPlaneStrain2D::yieldFunction(const double norm_dev_stress,
                                                      const Properties& rMaterialProperties)
{
    const double sqrt_two_thirds = std::sqrt(2.0 / 3.0);
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    const double k_old =
        yield_stress + (theta * hardening_modulus * mAccumulatedPlasticStrainOld) +
        (delta_k) * (1. - std::exp(-hardening_exponent * mAccumulatedPlasticStrainOld));

    return norm_dev_stress - k_old * sqrt_two_thirds;
}

void LinearJ2PlasticityPlaneStrain2D::CalculateElasticMatrix(Matrix &D, const Properties &props)
{
    const double E = props[YOUNG_MODULUS];
    const double poisson_ratio = props[POISSON_RATIO];
    const double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    const double mu = E / (2. + 2. * poisson_ratio);

    if (D.size1() != 4 || D.size2() != 4)
        D.resize(4, 4, false);
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
    double dgamma, double NormSTrial, const Vector& N_new, const Properties& props, Matrix& D)
{
    const double hardening_modulus = props[ISOTROPIC_HARDENING_MODULUS];
    const double theta = props[REFERENCE_HARDENING_MODULUS];
    const double delta_k = props[INFINITY_HARDENING_MODULUS];
    const double hardening_exponent = props[HARDENING_EXPONENT];
    const double E = props[YOUNG_MODULUS];
    const double poisson_ratio = props[POISSON_RATIO];
    const double mu = E / (2. + 2. * poisson_ratio);
    const double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));

    const double kp_new = (theta * hardening_modulus) +
                    delta_k * (hardening_exponent *
                               std::exp(-hardening_exponent * mAccumulatedPlasticStrain));

    const double theta_new = 1 - (2. * mu * dgamma) / NormSTrial;
    const double theta_new_b = 1. / (1. + kp_new / (3. * mu)) - (1. - theta_new);

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
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(REFERENCE_HARDENING_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(ISOTROPIC_HARDENING_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(INFINITY_HARDENING_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(HARDENING_EXPONENT));

    return 0;
}

void LinearJ2PlasticityPlaneStrain2D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mInelasticFlag", mInelasticFlag);
    rSerializer.save("mStrainEnergy", mStrainEnergy);
    rSerializer.save("mPlasticStrain", mPlasticStrain);
    rSerializer.save("mPlasticStrainOld", mPlasticStrainOld);
    rSerializer.save("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
    rSerializer.save("mAccumulatedPlasticStrainOld", mAccumulatedPlasticStrainOld);
}

void LinearJ2PlasticityPlaneStrain2D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mInelasticFlag", mInelasticFlag);
    rSerializer.load("mStrainEnergy", mStrainEnergy);
    rSerializer.load("mPlasticStrain", mPlasticStrain);
    rSerializer.load("mPlasticStrainOld", mPlasticStrainOld);
    rSerializer.load("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
    rSerializer.load("mAccumulatedPlasticStrainOld", mAccumulatedPlasticStrainOld);
}

} /* namespace Kratos.*/
