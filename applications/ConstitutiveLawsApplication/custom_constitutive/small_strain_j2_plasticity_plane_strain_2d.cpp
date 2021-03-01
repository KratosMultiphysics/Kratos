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

// System includes

// External includes

// Project includes
#include "small_strain_j2_plasticity_plane_strain_2d.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainJ2PlasticityPlaneStrain2D::SmallStrainJ2PlasticityPlaneStrain2D()
    : SmallStrainJ2Plasticity3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

SmallStrainJ2PlasticityPlaneStrain2D::SmallStrainJ2PlasticityPlaneStrain2D(const SmallStrainJ2PlasticityPlaneStrain2D &rOther)
    : SmallStrainJ2Plasticity3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainJ2PlasticityPlaneStrain2D::Clone() const
{
    return Kratos::make_shared<SmallStrainJ2PlasticityPlaneStrain2D>(SmallStrainJ2PlasticityPlaneStrain2D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SmallStrainJ2PlasticityPlaneStrain2D::~SmallStrainJ2PlasticityPlaneStrain2D()
{
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2PlasticityPlaneStrain2D::CalculateStressResponse(
        ConstitutiveLaw::Parameters& rValues,
        Vector& rPlasticStrain,
        double& rAccumulatedPlasticStrain)
{
    rPlasticStrain.resize(4, false);
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    Flags& r_constitutive_law_options = rValues.GetOptions();
    Vector& r_strain_vector = rValues.GetStrainVector();
    if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
        noalias(r_strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
    }

    if( r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        //this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    if( r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS) ||
        r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Vector& r_stress_vector = rValues.GetStressVector();
        Matrix& tangent_tensor = rValues.GetConstitutiveMatrix();
        const double hardening_modulus = r_material_properties[ISOTROPIC_HARDENING_MODULUS];
        const double delta_k = r_material_properties[EXPONENTIAL_SATURATION_YIELD_STRESS] - r_material_properties[YIELD_STRESS];
        const double hardening_exponent = r_material_properties[HARDENING_EXPONENT];
        const double E = r_material_properties[YOUNG_MODULUS];
        const double poisson_ratio = r_material_properties[POISSON_RATIO];
        const double mu = E / (2. + 2. * poisson_ratio);
        const double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));
        const double sqrt_two_thirds = std::sqrt(2. / 3.); // = 0.8164965809277260
        double trial_yield_function;

        rPlasticStrain = mPlasticStrain;
        rAccumulatedPlasticStrain = mAccumulatedPlasticStrain;

        Matrix elastic_tensor;
        elastic_tensor.resize(4, 4, false);
        CalculateElasticMatrix(r_material_properties, elastic_tensor);
        Vector yield_tension(4);
        noalias(yield_tension) = prod(elastic_tensor, r_strain_vector - mPlasticStrain);

        // stress_trial_dev = sigma - 1/3 tr(sigma) * I
        Vector stress_trial_dev = yield_tension;

        const double trace = 1. / 3. * (yield_tension(0) + yield_tension(1) + yield_tension(2));
        stress_trial_dev(0) -= trace;
        stress_trial_dev(1) -= trace;
        stress_trial_dev(2) -= trace;
        const double norm_dev_stress = std::sqrt(stress_trial_dev(0) * stress_trial_dev(0) +
                                           stress_trial_dev(1) * stress_trial_dev(1) +
                                           stress_trial_dev(2) * stress_trial_dev(2) +
                                           2. * stress_trial_dev(3) * stress_trial_dev(3));
        trial_yield_function = this->YieldFunction(norm_dev_stress, r_material_properties, mAccumulatedPlasticStrain);

        if (trial_yield_function <= 0.) {
            // ELASTIC
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
                r_stress_vector = yield_tension;
            }
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
                tangent_tensor = elastic_tensor;
            }
        } else {
            // INELASTIC
            double dgamma = 0;
            Vector yield_function_normal_vector = stress_trial_dev / norm_dev_stress;
            if (delta_k != 0.0 && hardening_exponent != 0.0) {
                // Exponential hardening
                dgamma = GetAccumPlasticStrainRate(norm_dev_stress, r_material_properties, mAccumulatedPlasticStrain);
            }
            else {
                // Linear hardening
                dgamma = trial_yield_function /
                         (2. * mu * (1. + (hardening_modulus / (3. * mu))));
            }

            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
                r_stress_vector(0) =
                    volumetric_modulus * (r_strain_vector(0) + r_strain_vector(1) + r_strain_vector(2)) +
                    stress_trial_dev(0) - 2. * mu * dgamma * yield_function_normal_vector(0);
                r_stress_vector(1) =
                    volumetric_modulus * (r_strain_vector(0) + r_strain_vector(1) + r_strain_vector(2)) +
                    stress_trial_dev(1) - 2. * mu * dgamma * yield_function_normal_vector(1);
                r_stress_vector(2) =
                    volumetric_modulus * (r_strain_vector(0) + r_strain_vector(1) + r_strain_vector(2)) +
                    stress_trial_dev(2) - 2. * mu * dgamma * yield_function_normal_vector(2);
                r_stress_vector(3) =
                        stress_trial_dev(3) - 2. * mu * dgamma * yield_function_normal_vector(3);
                }

            rPlasticStrain(0) += dgamma * yield_function_normal_vector(0);
            rPlasticStrain(1) += dgamma * yield_function_normal_vector(1);
            rPlasticStrain(2) += dgamma * yield_function_normal_vector(2);
            rPlasticStrain(3) += dgamma * yield_function_normal_vector(3) * 2;
            rAccumulatedPlasticStrain += sqrt_two_thirds * dgamma;

            // We update the tangent tensor
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
                CalculateTangentMatrix(dgamma, norm_dev_stress, yield_function_normal_vector,
                                       r_material_properties, rAccumulatedPlasticStrain, tangent_tensor);
            }
        }
    }
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2PlasticityPlaneStrain2D::CalculateElasticMatrix(const Properties &rMaterialProperties, Matrix &rElasticityTensor)
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

void SmallStrainJ2PlasticityPlaneStrain2D::CalculateTangentMatrix(const double DeltaGamma, const double NormStressTrial,
                                                             const Vector &YieldFunctionNormalVector,
                                                             const Properties &rMaterialProperties,
                                                             const double AccumulatedPlasticStrain,
                                                             Matrix &rElasticityTensor)
{
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[EXPONENTIAL_SATURATION_YIELD_STRESS] - rMaterialProperties[YIELD_STRESS];
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double mu = E / (2. + 2. * poisson_ratio);
    const double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));

    const double kp_new = hardening_modulus +
                    delta_k * (hardening_exponent *
                               std::exp(-hardening_exponent * AccumulatedPlasticStrain));

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

void SmallStrainJ2PlasticityPlaneStrain2D::GetLawFeatures(Features& rFeatures)
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

void SmallStrainJ2PlasticityPlaneStrain2D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainJ2Plasticity3D);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2PlasticityPlaneStrain2D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainJ2Plasticity3D);
}

} /* namespace Kratos.*/
