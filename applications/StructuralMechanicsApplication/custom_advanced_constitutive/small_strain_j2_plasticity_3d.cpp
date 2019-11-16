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
#include "small_strain_j2_plasticity_3d.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "includes/checks.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainJ2Plasticity3D::SmallStrainJ2Plasticity3D()
    : ConstitutiveLaw()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

SmallStrainJ2Plasticity3D::SmallStrainJ2Plasticity3D(const SmallStrainJ2Plasticity3D &rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainJ2Plasticity3D::Clone() const
{
    return Kratos::make_shared<SmallStrainJ2Plasticity3D>(SmallStrainJ2Plasticity3D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SmallStrainJ2Plasticity3D::~SmallStrainJ2Plasticity3D()
{
}

//************************************************************************************
//************************************************************************************

bool SmallStrainJ2Plasticity3D::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == ACCUMULATED_PLASTIC_STRAIN){
        return true;
    }
    return false;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rThisVariable == ACCUMULATED_PLASTIC_STRAIN){
        mAccumulatedPlasticStrain = rValue;
    }
}

//************************************************************************************
//************************************************************************************

double& SmallStrainJ2Plasticity3D::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == ACCUMULATED_PLASTIC_STRAIN){
        rValue = mAccumulatedPlasticStrain;
    }

    return rValue;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    mPlasticStrain = ZeroVector(this->GetStrainSize());
    mAccumulatedPlasticStrain = 0.0;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::InitializeMaterialResponseCauchy(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::InitializeMaterialResponsePK2(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::InitializeMaterialResponsePK1(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::InitializeMaterialResponseKirchhoff(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    InitializeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::FinalizeMaterialResponseCauchy(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    Vector plastic_strain;
    double accumulated_plastic_strain;
    this->CalculateStressResponse(rValues, plastic_strain, accumulated_plastic_strain);
    mPlasticStrain = plastic_strain;
    mAccumulatedPlasticStrain = accumulated_plastic_strain;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::FinalizeMaterialResponsePK2(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::FinalizeMaterialResponsePK1(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::FinalizeMaterialResponseKirchhoff(
    Kratos::ConstitutiveLaw::Parameters &rValues)
{
    // In small deformation is the same as compute Cauchy
    FinalizeMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // In small deformation is the same as compute Cauchy
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // In small deformation is the same as compute Cauchy
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    Vector plastic_strain;
    double accumulated_plastic_strain;
    this->CalculateStressResponse(rValues, plastic_strain, accumulated_plastic_strain);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rPlasticStrain,
    double& rAccumulatedPlasticStrain)
{
    rPlasticStrain.resize(6, false);
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    Flags& r_constitutive_law_options = rValues.GetOptions();
    Vector& r_strain_vector = rValues.GetStrainVector();
    if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
        noalias(r_strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
    }

    if( r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        //this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // If we compute the tangent moduli or the stress
    if( r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS) ||
        r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Vector& r_stress_vector = rValues.GetStressVector();
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        const double hardening_modulus = r_material_properties[ISOTROPIC_HARDENING_MODULUS];
        const double delta_k = r_material_properties[EXPONENTIAL_SATURATION_YIELD_STRESS] - r_material_properties[YIELD_STRESS];
        const double hardening_exponent = r_material_properties[HARDENING_EXPONENT];
        const double E = r_material_properties[YOUNG_MODULUS];
        const double poisson_ratio = r_material_properties[POISSON_RATIO];
        const double mu = E / (2. + 2. * poisson_ratio);
        const double bulk_modulus = E / (3. * (1. - 2. * poisson_ratio));
        const double sqrt_two_thirds = std::sqrt(2. / 3.); // = 0.8164965809277260
        double trial_yield_function;

        rPlasticStrain = mPlasticStrain;
        rAccumulatedPlasticStrain = mAccumulatedPlasticStrain;

        Matrix elastic_tensor;
        elastic_tensor.resize(6, 6, false);
        CalculateElasticMatrix(r_material_properties, elastic_tensor);
        Vector trial_stress(6);
        noalias(trial_stress) = prod(elastic_tensor, r_strain_vector - mPlasticStrain);

        // trial_stress_dev = sigma - 1/3 tr(sigma) * I
        Vector trial_stress_dev = trial_stress;
        const double trace = 1. / 3. * (trial_stress(0) + trial_stress(1) + trial_stress(2));
        trial_stress_dev(0) -= trace;
        trial_stress_dev(1) -= trace;
        trial_stress_dev(2) -= trace;
        const double norm_trial_stress_dev = std::sqrt(trial_stress_dev(0) * trial_stress_dev(0) +
                                        trial_stress_dev(1) * trial_stress_dev(1) +
                                        trial_stress_dev(2) * trial_stress_dev(2) +
                                        2. * trial_stress_dev(3) * trial_stress_dev(3) +
                                        2. * trial_stress_dev(4) * trial_stress_dev(4) +
                                        2. * trial_stress_dev(5) * trial_stress_dev(5));
        trial_yield_function = this->YieldFunction(norm_trial_stress_dev, r_material_properties, mAccumulatedPlasticStrain);

        if (trial_yield_function <= 0.) {
            // ELASTIC
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
                r_stress_vector = trial_stress;
            }
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
                r_constitutive_matrix = elastic_tensor;
            }
        } else {
            // INELASTIC
            double accum_plastic_strain_rate = 0;
            Vector yield_function_normal_vector = trial_stress_dev / norm_trial_stress_dev;
            if (delta_k != 0.0 && hardening_exponent != 0.0) {
                // Exponential hardening
                accum_plastic_strain_rate = GetAccumPlasticStrainRate(norm_trial_stress_dev, r_material_properties,
                                                                      mAccumulatedPlasticStrain);
            }
            else {
                // Linear hardening
                accum_plastic_strain_rate = trial_yield_function /
                        (2. * mu * (1. + (hardening_modulus / (3. * mu))));
            }
            // We update the stress
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
                r_stress_vector(0) =
                    bulk_modulus * (r_strain_vector(0) + r_strain_vector(1) + r_strain_vector(2)) +
                    trial_stress_dev(0) - 2. * mu * accum_plastic_strain_rate * yield_function_normal_vector(0);
                r_stress_vector(1) =
                    bulk_modulus * (r_strain_vector(0) + r_strain_vector(1) + r_strain_vector(2)) +
                    trial_stress_dev(1) - 2. * mu * accum_plastic_strain_rate * yield_function_normal_vector(1);
                r_stress_vector(2) =
                    bulk_modulus * (r_strain_vector(0) + r_strain_vector(1) + r_strain_vector(2)) +
                    trial_stress_dev(2) - 2. * mu * accum_plastic_strain_rate * yield_function_normal_vector(2);
                r_stress_vector(3) =
                    trial_stress_dev(3) - 2. * mu * accum_plastic_strain_rate * yield_function_normal_vector(3);
                r_stress_vector(4) =
                    trial_stress_dev(4) - 2. * mu * accum_plastic_strain_rate * yield_function_normal_vector(4);
                r_stress_vector(5) =
                    trial_stress_dev(5) - 2. * mu * accum_plastic_strain_rate * yield_function_normal_vector(5);
            }

            rPlasticStrain(0) += accum_plastic_strain_rate * yield_function_normal_vector(0);
            rPlasticStrain(1) += accum_plastic_strain_rate * yield_function_normal_vector(1);
            rPlasticStrain(2) += accum_plastic_strain_rate * yield_function_normal_vector(2);
            rPlasticStrain(3) += accum_plastic_strain_rate * yield_function_normal_vector(3) * 2;
            rPlasticStrain(4) += accum_plastic_strain_rate * yield_function_normal_vector(4) * 2;
            rPlasticStrain(5) += accum_plastic_strain_rate * yield_function_normal_vector(5) * 2;
            rAccumulatedPlasticStrain += sqrt_two_thirds * accum_plastic_strain_rate;

            // We update the tangent tensor
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
                CalculateTangentMatrix(accum_plastic_strain_rate, norm_trial_stress_dev, yield_function_normal_vector,
                                       r_material_properties, rAccumulatedPlasticStrain, r_constitutive_matrix);
            }
        }
    }
}

//************************************************************************************
//************************************************************************************

double& SmallStrainJ2Plasticity3D::CalculateValue(
    ConstitutiveLaw::Parameters& rValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == STRAIN_ENERGY){
        Vector& r_strain_vector = rValues.GetStrainVector();
        if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(r_strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
        }
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        Matrix elastic_tensor;
        CalculateElasticMatrix(r_material_properties, elastic_tensor);

        rValue = 0.5 * inner_prod(r_strain_vector - mPlasticStrain,
                                  prod(elastic_tensor, r_strain_vector - mPlasticStrain))
                 + GetPlasticPotential(r_material_properties, mAccumulatedPlasticStrain);
    }

    return(rValue);
}

//************************************************************************************
//************************************************************************************

Vector& SmallStrainJ2Plasticity3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == STRAIN) {
        const SizeType space_dimension = this->WorkingSpaceDimension();

        // Compute total deformation gradient
        const Matrix& r_F = rParameterValues.GetDeformationGradientF();
        KRATOS_DEBUG_ERROR_IF(r_F.size1()!= space_dimension || r_F.size2() != space_dimension)
            << "expected size of F " << space_dimension << "x" << space_dimension
            << ", got " << r_F.size1() << "x" << r_F.size2() << std::endl;

        const Matrix C_tensor = prod(trans(r_F),r_F);
        ConstitutiveLawUtilities<6>::CalculateGreenLagrangianStrain(C_tensor, rValue);
    }

    return(rValue);
}

//************************************************************************************
//************************************************************************************

double SmallStrainJ2Plasticity3D::GetSaturationHardening(const Properties& rMaterialProperties,
    const double accumulated_plastic_strain)
{
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[EXPONENTIAL_SATURATION_YIELD_STRESS] - yield_stress;
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];

    const double k_new = yield_stress + (hardening_modulus * accumulated_plastic_strain) +
                delta_k * (1. - std::exp(-hardening_exponent * accumulated_plastic_strain));
    return k_new;
}

//************************************************************************************
//************************************************************************************

double SmallStrainJ2Plasticity3D::GetPlasticPotential(const Properties& rMaterialProperties,
    const double accumulated_plastic_strain)
{
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[EXPONENTIAL_SATURATION_YIELD_STRESS] - rMaterialProperties[YIELD_STRESS];
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];

    // perfect plasticity energy value
    double plastic_free_energy = 0.;

    // linear hardening contribution
    if (hardening_modulus != 0.){
        plastic_free_energy += 0.5 *( hardening_modulus * std::pow(accumulated_plastic_strain, 2.0));
    }

    // exponential hardening contribution
    if (hardening_exponent != 0.) {
        plastic_free_energy += delta_k * (accumulated_plastic_strain
            + (1 / hardening_exponent) * std::exp( -hardening_exponent * accumulated_plastic_strain));
    }
    return plastic_free_energy;
}

//************************************************************************************
//************************************************************************************

double SmallStrainJ2Plasticity3D::GetAccumPlasticStrainRate(
        const double NormStressTrial,
        const Properties &rMaterialProperties,
        const double AccumulatedPlasticStrainOld
)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[EXPONENTIAL_SATURATION_YIELD_STRESS] - yield_stress;
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    const double tolerance = 1e-6 * yield_stress;
    const double mu = E / (2. * (1. + poisson_ratio));
    const double sqrt_two_thirds = std::sqrt(2. / 3.); // = 0.8164965809277260
    double dgamma = 0.0;
    double norm_yieldfunction = 1.0;
    double accumulated_plastic_strain = AccumulatedPlasticStrainOld;

    while (norm_yieldfunction > tolerance)
    {
        const double k_new = GetSaturationHardening(rMaterialProperties, AccumulatedPlasticStrainOld);
        const double kp_new = hardening_modulus + delta_k * (hardening_exponent * std::exp(-hardening_exponent * accumulated_plastic_strain));
        const double yieldfunction = - sqrt_two_thirds * k_new + NormStressTrial - 2. * mu * dgamma;
        const double derivative_yieldfunction = -2. * mu * (1. + kp_new / (3. * mu));
        dgamma -= yieldfunction / derivative_yieldfunction;
        accumulated_plastic_strain = AccumulatedPlasticStrainOld + sqrt_two_thirds * dgamma;
        norm_yieldfunction = std::abs(yieldfunction);
    }
    // TODO(@marandra): handle the case when no convergence is achieved:
    // limit the number of iterations and throw a warning/error

    return dgamma;
}

//************************************************************************************
//************************************************************************************

double SmallStrainJ2Plasticity3D::YieldFunction(
    const double NormDeviationStress,
    const Properties& rMaterialProperties,
    const double AccumulatedPlasticStrain
    )
{
    const double sqrt_two_thirds = std::sqrt(2. / 3.);
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[EXPONENTIAL_SATURATION_YIELD_STRESS] - yield_stress;
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    const double k_old = yield_stress + (hardening_modulus * AccumulatedPlasticStrain) +
        (delta_k) * (1. - std::exp(-hardening_exponent * AccumulatedPlasticStrain));

    return NormDeviationStress - k_old * sqrt_two_thirds;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::CalculateElasticMatrix(
    const Properties &rMaterialProperties, Matrix &rElasticMatrix)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    const double mu = E / (2. + 2. * poisson_ratio);

    if (rElasticMatrix.size1() != 6 || rElasticMatrix.size2() != 6)
        rElasticMatrix.resize(6, 6, false);
    rElasticMatrix.clear();

    rElasticMatrix(0, 0) = lambda + 2. * mu;
    rElasticMatrix(0, 1) = lambda;
    rElasticMatrix(0, 2) = lambda;
    rElasticMatrix(1, 0) = lambda;
    rElasticMatrix(1, 1) = lambda + 2. * mu;
    rElasticMatrix(1, 2) = lambda;
    rElasticMatrix(2, 0) = lambda;
    rElasticMatrix(2, 1) = lambda;
    rElasticMatrix(2, 2) = lambda + 2. * mu;
    rElasticMatrix(3, 3) = mu;
    rElasticMatrix(4, 4) = mu;
    rElasticMatrix(5, 5) = mu;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::CalculateTangentMatrix(
        const double DeltaGamma, const double NormStressTrial,
        const Vector &rYFNormalVector,
        const Properties &rMaterialProperties,
        const double AccumulatedPlasticStrain,
        Matrix &rTMatrix)
{
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[EXPONENTIAL_SATURATION_YIELD_STRESS] - rMaterialProperties[YIELD_STRESS];
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double mu = E / (2. + 2. * poisson_ratio);
    const double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));

    const double kp_new = hardening_modulus
            + delta_k * (hardening_exponent * std::exp(-hardening_exponent * AccumulatedPlasticStrain));

    const double theta_new = 1 - (2. * mu * DeltaGamma) / NormStressTrial;
    const double theta_new_b = 1. / (1. + kp_new / (3. * mu)) - (1. - theta_new);

    rTMatrix(0, 0) = volumetric_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (rYFNormalVector(0) * rYFNormalVector(0)));
    rTMatrix(0, 1) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (rYFNormalVector(0) * rYFNormalVector(1)));
    rTMatrix(0, 2) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (rYFNormalVector(0) * rYFNormalVector(2)));
    rTMatrix(0, 3) = -(2. * mu * theta_new_b * (rYFNormalVector(0) * rYFNormalVector(3)));
    rTMatrix(0, 4) = -(2. * mu * theta_new_b * (rYFNormalVector(0) * rYFNormalVector(4)));
    rTMatrix(0, 5) = -(2. * mu * theta_new_b * (rYFNormalVector(0) * rYFNormalVector(5)));

    rTMatrix(1, 0) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (rYFNormalVector(1) * rYFNormalVector(0)));
    rTMatrix(1, 1) = volumetric_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (rYFNormalVector(1) * rYFNormalVector(1)));
    rTMatrix(1, 2) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (rYFNormalVector(1) * rYFNormalVector(2)));
    rTMatrix(1, 3) = -(2. * mu * theta_new_b * (rYFNormalVector(1) * rYFNormalVector(3)));
    rTMatrix(1, 4) = -(2. * mu * theta_new_b * (rYFNormalVector(1) * rYFNormalVector(4)));
    rTMatrix(1, 5) = -(2. * mu * theta_new_b * (rYFNormalVector(1) * rYFNormalVector(5)));

    rTMatrix(2, 0) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (rYFNormalVector(2) * rYFNormalVector(0)));
    rTMatrix(2, 1) = volumetric_modulus + (2. * mu * theta_new * (-1. / 3.)) -
              (2. * mu * theta_new_b * (rYFNormalVector(2) * rYFNormalVector(1)));
    rTMatrix(2, 2) = volumetric_modulus + (2. * mu * theta_new * 2. / 3.) -
              (2. * mu * theta_new_b * (rYFNormalVector(2) * rYFNormalVector(2)));
    rTMatrix(2, 3) = -(2. * mu * theta_new_b * (rYFNormalVector(2) * rYFNormalVector(3)));
    rTMatrix(2, 4) = -(2. * mu * theta_new_b * (rYFNormalVector(2) * rYFNormalVector(4)));
    rTMatrix(2, 5) = -(2. * mu * theta_new_b * (rYFNormalVector(2) * rYFNormalVector(5)));

    rTMatrix(3, 0) = -(2. * mu * theta_new_b * (rYFNormalVector(3) * rYFNormalVector(0)));
    rTMatrix(3, 1) = -(2. * mu * theta_new_b * (rYFNormalVector(3) * rYFNormalVector(1)));
    rTMatrix(3, 2) = -(2. * mu * theta_new_b * (rYFNormalVector(3) * rYFNormalVector(2)));
    rTMatrix(3, 3) = mu * theta_new - (2. * mu * theta_new_b * (rYFNormalVector(3) * rYFNormalVector(3)));
    rTMatrix(3, 4) = -(2. *mu * theta_new_b * (rYFNormalVector(3) * rYFNormalVector(4)));
    rTMatrix(3, 5) = -(2. *mu * theta_new_b * (rYFNormalVector(3) * rYFNormalVector(5)));

    rTMatrix(4, 0) = -(2. * mu * theta_new_b * (rYFNormalVector(4) * rYFNormalVector(0)));
    rTMatrix(4, 1) = -(2. * mu * theta_new_b * (rYFNormalVector(4) * rYFNormalVector(1)));
    rTMatrix(4, 2) = -(2. * mu * theta_new_b * (rYFNormalVector(4) * rYFNormalVector(2)));
    rTMatrix(4, 3) = -(2. * mu * theta_new_b * (rYFNormalVector(4) * rYFNormalVector(3)));
    rTMatrix(4, 4) = mu * theta_new - (2. * mu * theta_new_b * (rYFNormalVector(4) * rYFNormalVector(4)));
    rTMatrix(4, 5) = -(2. * mu * theta_new_b * (rYFNormalVector(4) * rYFNormalVector(5)));

    rTMatrix(5, 0) = -(2. * mu * theta_new_b * (rYFNormalVector(5) * rYFNormalVector(0)));
    rTMatrix(5, 1) = -(2. * mu * theta_new_b * (rYFNormalVector(5) * rYFNormalVector(1)));
    rTMatrix(5, 2) = -(2. * mu * theta_new_b * (rYFNormalVector(5) * rYFNormalVector(2)));
    rTMatrix(5, 3) = -(2. * mu * theta_new_b * (rYFNormalVector(5) * rYFNormalVector(3)));
    rTMatrix(5, 4) = -(2. * mu * theta_new_b * (rYFNormalVector(5) * rYFNormalVector(4)));
    rTMatrix(5, 5) = mu * theta_new - (2. * mu * theta_new_b * (rYFNormalVector(5) * rYFNormalVector(5)));
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::GetLawFeatures(Features& rFeatures)
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

int SmallStrainJ2Plasticity3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(ISOTROPIC_HARDENING_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(EXPONENTIAL_SATURATION_YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(HARDENING_EXPONENT));

    return 0;
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mPlasticStrain", mPlasticStrain);
    rSerializer.save("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
}

//************************************************************************************
//************************************************************************************

void SmallStrainJ2Plasticity3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mPlasticStrain", mPlasticStrain);
    rSerializer.load("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
}

} /* namespace Kratos.*/
