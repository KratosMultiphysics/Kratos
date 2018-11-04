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
//  Collaborator:    Vicente Mataix Ferrandiz
//

#include "linear_j2_plasticity_plane_strain_2d.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearJ2PlasticityPlaneStrain2D::LinearJ2PlasticityPlaneStrain2D()
    : LinearJ2Plasticity3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

LinearJ2PlasticityPlaneStrain2D::LinearJ2PlasticityPlaneStrain2D(const LinearJ2PlasticityPlaneStrain2D &rOther)
    : LinearJ2Plasticity3D(rOther)
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

void LinearJ2PlasticityPlaneStrain2D::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
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
        Vector yield_tensionrial(4);
        noalias(yield_tensionrial) = prod(elastic_tensor, strain_vector - mPlasticStrainOld);

        // stress_trial_dev = sigma - 1/3 tr(sigma) * I
        Vector stress_trial_dev = yield_tensionrial;

        const double trace = 1.0 / 3.0 * (yield_tensionrial(0) + yield_tensionrial(1) + yield_tensionrial(2));
        stress_trial_dev(0) -= trace;
        stress_trial_dev(1) -= trace;
        stress_trial_dev(2) -= trace;
        const double norm_dev_stress = std::sqrt(stress_trial_dev(0) * stress_trial_dev(0) +
                                           stress_trial_dev(1) * stress_trial_dev(1) +
                                           stress_trial_dev(2) * stress_trial_dev(2) +
                                           2. * stress_trial_dev(3) * stress_trial_dev(3));
        trial_yield_function = this->YieldFunction(norm_dev_stress, rMaterialProperties);

        if (trial_yield_function <= 0.) {
            // ELASTIC
            mInelasticFlag = false;
            stress_vector = yield_tensionrial;
            tangent_tensor = elastic_tensor;
        } else {
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
    }
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************


void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************


void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::CalculateElasticMatrix(
    Matrix &rElasticityTensor,
    const Properties &rMaterialProperties
    )
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    const double mu = E / (2. + 2. * poisson_ratio);

    if (rElasticityTensor.size1() != 4 || rElasticityTensor.size2() != 4)
        rElasticityTensor.resize(4, 4, false);
    rElasticityTensor.clear();

    rElasticityTensor(0, 0) = lambda + 2. * mu;
    rElasticityTensor(0, 1) = lambda;
    rElasticityTensor(0, 2) = lambda;
    rElasticityTensor(0, 3) = 0.;
    rElasticityTensor(1, 0) = lambda;
    rElasticityTensor(1, 1) = lambda + 2. * mu;
    rElasticityTensor(1, 2) = lambda;
    rElasticityTensor(1, 3) = 0.;
    rElasticityTensor(2, 0) = lambda;
    rElasticityTensor(2, 1) = lambda;
    rElasticityTensor(2, 2) = lambda + 2. * mu;
    rElasticityTensor(2, 3) = 0.;
    rElasticityTensor(3, 0) = 0.;
    rElasticityTensor(3, 1) = 0.;
    rElasticityTensor(3, 2) = 0.;
    rElasticityTensor(3, 3) = mu;
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::CalculateTangentTensor(
    const double DeltaGamma,
    const double NormStressTrial,
    const Vector& YieldFunctionNormalVector,
    const Properties& rMaterialProperties,
    Matrix& rElasticityTensor
    )
{
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double mu = E / (2. + 2. * poisson_ratio);
    const double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));

    const double kp_new = (theta * hardening_modulus) +
                    delta_k * (hardening_exponent *
                               std::exp(-hardening_exponent * mAccumulatedPlasticStrain));

    const double theta_new = 1 - (2. * mu * DeltaGamma) / NormStressTrial;
    const double theta_new_b = 1. / (1. + kp_new / (3. * mu)) - (1. - theta_new);

    rElasticityTensor(0, 0) = volumetric_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(0)));
    rElasticityTensor(0, 1) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(1)));
    rElasticityTensor(0, 2) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(2)));
    rElasticityTensor(0, 3) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(3)));

    rElasticityTensor(1, 0) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(0)));
    rElasticityTensor(1, 1) = volumetric_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(1)));
    rElasticityTensor(1, 2) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(2)));
    rElasticityTensor(1, 3) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(3)));

    rElasticityTensor(2, 0) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(0)));
    rElasticityTensor(2, 1) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(1)));
    rElasticityTensor(2, 2) = volumetric_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(2)));
    rElasticityTensor(2, 3) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(3)));

    rElasticityTensor(3, 0) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(0)));
    rElasticityTensor(3, 1) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(1)));
    rElasticityTensor(3, 2) = -(2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(2)));
    rElasticityTensor(3, 3) = mu * theta_new - (2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(3)));
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = 4;
    rFeatures.mSpaceDimension = 2;
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LinearJ2Plasticity3D);
}

//************************************************************************************
//************************************************************************************

void LinearJ2PlasticityPlaneStrain2D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LinearJ2Plasticity3D);
}

} /* namespace Kratos.*/
