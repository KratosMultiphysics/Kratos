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
#include "linear_j2_plasticity_3d.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearJ2Plasticity3D::LinearJ2Plasticity3D()
    : ConstitutiveLaw()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

LinearJ2Plasticity3D::LinearJ2Plasticity3D(const LinearJ2Plasticity3D &rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearJ2Plasticity3D::Clone() const
{
    LinearJ2Plasticity3D::Pointer p_clone(new LinearJ2Plasticity3D(*this));
    return p_clone;
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

LinearJ2Plasticity3D::~LinearJ2Plasticity3D()
{
}

//************************************************************************************
//************************************************************************************

bool LinearJ2Plasticity3D::Has(const Variable<bool>& rThisVariable)
{
    if(rThisVariable == INELASTIC_FLAG){
        return true;
    }
    return false;
}

//************************************************************************************
//************************************************************************************

bool LinearJ2Plasticity3D::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == PLASTIC_STRAIN){
        return true;
    }
    return false;
}



//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rThisVariable == PLASTIC_STRAIN){
        mAccumulatedPlasticStrain = rValue;
    }
}

//************************************************************************************
//************************************************************************************

bool& LinearJ2Plasticity3D::GetValue(
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

double& LinearJ2Plasticity3D::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == PLASTIC_STRAIN){
        rValue = mAccumulatedPlasticStrain;
    }

    return rValue;
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    mPlasticStrainOld = ZeroVector(this->GetStrainSize());
    mPlasticStrain = ZeroVector(this->GetStrainSize());
    mAccumulatedPlasticStrainOld = 0.0;
    mAccumulatedPlasticStrain = 0.0;
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    mPlasticStrainOld = mPlasticStrain;
    mAccumulatedPlasticStrainOld = mAccumulatedPlasticStrain;
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // The Properties of the material
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // The flags of the law
    Flags & r_constitutive_law_options=rValues.GetOptions();

    // The strain tensor
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if( r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // If we compute the tangent moduli or the stress
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ||
        r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR )) {
        Vector& r_stress_vector = rValues.GetStressVector();

        if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(r_strain_vector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
        }

        Matrix elastic_tensor;
        Matrix& tangent_tensor = rValues.GetConstitutiveMatrix();
        const double hardening_modulus = r_material_properties[ISOTROPIC_HARDENING_MODULUS];
        const double delta_k = r_material_properties[INFINITY_HARDENING_MODULUS];
        const double hardening_exponent = r_material_properties[HARDENING_EXPONENT];
        const double E = r_material_properties[YOUNG_MODULUS];
        const double poisson_ratio = r_material_properties[POISSON_RATIO];
        const double mu = E / (2. + 2. * poisson_ratio);
        const double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));
        const double sqrt_two_thirds = std::sqrt(2.0 / 3.0); // =0.8164965809277260
        double trial_yield_function;

        mPlasticStrain = mPlasticStrainOld;
        mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld;

        elastic_tensor.resize(6, 6, false);
        CalculateElasticMatrix(elastic_tensor, r_material_properties);
        Vector yield_tensionrial(6);
        noalias(yield_tensionrial) = prod(elastic_tensor, r_strain_vector - mPlasticStrainOld);

        // stress_trial_dev = sigma - 1/3 tr(sigma) * I
        Vector stress_trial_dev = yield_tensionrial;

        const double trace = 1.0 / 3.0 * (yield_tensionrial(0) + yield_tensionrial(1) + yield_tensionrial(2));
        stress_trial_dev(0) -= trace;
        stress_trial_dev(1) -= trace;
        stress_trial_dev(2) -= trace;
        const double norm_dev_stress = std::sqrt(stress_trial_dev(0) * stress_trial_dev(0) +
                                        stress_trial_dev(1) * stress_trial_dev(1) +
                                        stress_trial_dev(2) * stress_trial_dev(2) +
                                        2. * stress_trial_dev(3) * stress_trial_dev(3) +
                                        2. * stress_trial_dev(4) * stress_trial_dev(4) +
                                        2. * stress_trial_dev(5) * stress_trial_dev(5));
        trial_yield_function = this->YieldFunction(norm_dev_stress, r_material_properties);

        if (trial_yield_function <= 0.) {
            // ELASTIC
            mInelasticFlag = false;
            // We update the stress
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
                r_stress_vector = yield_tensionrial;
            }
            // We update the tangent tensor
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
                tangent_tensor = elastic_tensor;
            }
        } else {
            // INELASTIC
            mInelasticFlag = true;
            double dgamma = 0;
            Vector yield_function_normal_vector = stress_trial_dev / norm_dev_stress;
            if (delta_k != 0.0 && hardening_exponent != 0.0) {
                // Exponential softening
                dgamma = GetDeltaGamma(norm_dev_stress, r_material_properties);
            }
            else {
                // Linear softening
                dgamma = trial_yield_function /
                        (2. * mu * (1. + (hardening_modulus / (3. * mu))));
            }

            // We update the stress
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
                r_stress_vector(4) =
                    stress_trial_dev(4) - 2. * mu * dgamma * yield_function_normal_vector(4);
                r_stress_vector(5) =
                    stress_trial_dev(5) - 2. * mu * dgamma * yield_function_normal_vector(5);
            }

            mPlasticStrain(0) = mPlasticStrainOld(0) + dgamma * yield_function_normal_vector(0);
            mPlasticStrain(1) = mPlasticStrainOld(1) + dgamma * yield_function_normal_vector(1);
            mPlasticStrain(2) = mPlasticStrainOld(2) + dgamma * yield_function_normal_vector(2);
            mPlasticStrain(3) = mPlasticStrainOld(3) + dgamma * yield_function_normal_vector(3) * 2;
            mPlasticStrain(4) = mPlasticStrainOld(4) + dgamma * yield_function_normal_vector(4) * 2;
            mPlasticStrain(5) = mPlasticStrainOld(5) + dgamma * yield_function_normal_vector(5) * 2;
            mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld + sqrt_two_thirds * dgamma;

            // We update the tangent tensor
            if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
                CalculateTangentTensor(dgamma, norm_dev_stress, yield_function_normal_vector,
                                    r_material_properties, tangent_tensor);
            }
        }
    }
}

//************************************************************************************
//************************************************************************************

double& LinearJ2Plasticity3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == STRAIN_ENERGY){
        Vector& strain_vector = rParameterValues.GetStrainVector();
        if (rParameterValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(strain_vector) += rParameterValues.GetProcessInfo()[INITIAL_STRAIN];
        }
        const Properties& r_material_properties = rParameterValues.GetMaterialProperties();
        Matrix elastic_tensor(6, 6);
        CalculateElasticMatrix(elastic_tensor, r_material_properties);

        rValue = 0.5 * inner_prod(strain_vector - mPlasticStrain, prod(elastic_tensor, strain_vector - mPlasticStrain))
                 + GetPlasticPotential(r_material_properties);
    } else if(rThisVariable == PLASTIC_STRAIN){
        rValue = mAccumulatedPlasticStrain;
    }

    return(rValue);
}

//************************************************************************************
//************************************************************************************

Vector& LinearJ2Plasticity3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == STRAIN ||
        rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {

        const SizeType space_dimension = this->WorkingSpaceDimension();

        //1.-Compute total deformation gradient
        const Matrix& F = rParameterValues.GetDeformationGradientF();
        KRATOS_DEBUG_ERROR_IF(F.size1()!= space_dimension || F.size2() != space_dimension) << "expected size of F " << space_dimension << "x" << space_dimension << ", got " << F.size1() << "x" << F.size2() << std::endl;

        const Matrix C_tensor = prod(trans(F),F);
        ConstitutiveLawUtilities<6>::CalculateGreenLagrangianStrain(C_tensor, rValue);
    }

    return(rValue);
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

//************************************************************************************
//************************************************************************************

double LinearJ2Plasticity3D::GetSaturationHardening(const Properties& rMaterialProperties)
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

//************************************************************************************
//************************************************************************************

double LinearJ2Plasticity3D::GetPlasticPotential(const Properties& rMaterialProperties)
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

//************************************************************************************
//************************************************************************************

double LinearJ2Plasticity3D::GetDeltaGamma(
    const double NormStressTrial,
    const Properties& rMaterialProperties
    )
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
        const double kp_new = theta * hardening_modulus + delta_k * (hardening_exponent * std::exp(-hardening_exponent * mAccumulatedPlasticStrain));
        const double yieldfunction = - sqrt_two_thirds * k_new + NormStressTrial - 2. * mu * dgamma;
        const double derivative_yieldfunction = -2. * mu * (1. + kp_new / (3. * mu));
        dgamma -= yieldfunction / derivative_yieldfunction;
        mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld + sqrt_two_thirds * dgamma;
        norm_yieldfunction = std::abs(yieldfunction);
    }
    // TODO (marcelo): handle the case when no convergence is achieved.
    return dgamma;
}

//************************************************************************************
//************************************************************************************

double LinearJ2Plasticity3D::YieldFunction(
    const double NormDeviationStress,
    const Properties& rMaterialProperties
    )
{
    const double sqrt_two_thirds = std::sqrt(2.0 / 3.0);
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    const double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    const double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    const double k_old = yield_stress + (theta * hardening_modulus * mAccumulatedPlasticStrainOld) +
        (delta_k) * (1. - std::exp(-hardening_exponent * mAccumulatedPlasticStrainOld));

    return NormDeviationStress - k_old * sqrt_two_thirds;
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::CalculateElasticMatrix(
    Matrix &rElasticityTensor,
    const Properties &rMaterialProperties
    )
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    const double mu = E / (2. + 2. * poisson_ratio);

    if (rElasticityTensor.size1() != 6 || rElasticityTensor.size2() != 6)
        rElasticityTensor.resize(6, 6, false);
    rElasticityTensor.clear();

    rElasticityTensor(0, 0) = lambda + 2. * mu;
    rElasticityTensor(0, 1) = lambda;
    rElasticityTensor(0, 2) = lambda;
    rElasticityTensor(1, 0) = lambda;
    rElasticityTensor(1, 1) = lambda + 2. * mu;
    rElasticityTensor(1, 2) = lambda;
    rElasticityTensor(2, 0) = lambda;
    rElasticityTensor(2, 1) = lambda;
    rElasticityTensor(2, 2) = lambda + 2. * mu;
    rElasticityTensor(3, 3) = mu;
    rElasticityTensor(4, 4) = mu;
    rElasticityTensor(5, 5) = mu;
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::CalculateTangentTensor(
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

    const double kp_new = (theta * hardening_modulus) +  delta_k * (hardening_exponent * std::exp(-hardening_exponent * mAccumulatedPlasticStrain));

    const double theta_new = 1 - (2. * mu * DeltaGamma) / NormStressTrial;
    const double theta_new_b = 1. / (1. + kp_new / (3. * mu)) - (1. - theta_new);

    rElasticityTensor(0, 0) = volumetric_modulus + (2. *mu * theta_new * 2. / 3.) -
              (2. *mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(0)));
    rElasticityTensor(0, 1) = volumetric_modulus + (2. *mu * theta_new * (-1. / 3.)) -
              (2. *mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(1)));
    rElasticityTensor(0, 2) = volumetric_modulus + (2. *mu * theta_new * (-1. / 3.)) -
              (2. *mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(2)));
    rElasticityTensor(0, 3) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(3)));
    rElasticityTensor(0, 4) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(4)));
    rElasticityTensor(0, 5) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(0) * YieldFunctionNormalVector(5)));

    rElasticityTensor(1, 0) = volumetric_modulus + (2. *mu * theta_new * (-1. / 3.)) -
              (2. *mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(0)));
    rElasticityTensor(1, 1) = volumetric_modulus + (2. *mu * theta_new * 2. / 3.) -
              (2. *mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(1)));
    rElasticityTensor(1, 2) = volumetric_modulus + (2. *mu * theta_new * (-1. / 3.)) -
              (2. *mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(2)));
    rElasticityTensor(1, 3) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(3)));
    rElasticityTensor(1, 4) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(4)));
    rElasticityTensor(1, 5) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(1) * YieldFunctionNormalVector(5)));

    rElasticityTensor(2, 0) = volumetric_modulus + (2. *mu * theta_new * (-1. / 3.)) -
              (2. *mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(0)));
    rElasticityTensor(2, 1) = volumetric_modulus + (2. *mu * theta_new * (-1. / 3.)) -
              (2. *mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(1)));
    rElasticityTensor(2, 2) = volumetric_modulus + (2. *mu * theta_new * 2. / 3.) -
              (2. *mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(2)));
    rElasticityTensor(2, 3) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(3)));
    rElasticityTensor(2, 4) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(4)));
    rElasticityTensor(2, 5) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(2) * YieldFunctionNormalVector(5)));

    rElasticityTensor(3, 0) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(0)));
    rElasticityTensor(3, 1) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(1)));
    rElasticityTensor(3, 2) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(2)));
    rElasticityTensor(3, 3) = mu * theta_new - (2. * mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(3)));
    rElasticityTensor(3, 4) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(4)));
    rElasticityTensor(3, 5) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(3) * YieldFunctionNormalVector(5)));

    rElasticityTensor(4, 0) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(4) * YieldFunctionNormalVector(0)));
    rElasticityTensor(4, 1) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(4) * YieldFunctionNormalVector(1)));
    rElasticityTensor(4, 2) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(4) * YieldFunctionNormalVector(2)));
    rElasticityTensor(4, 3) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(4) * YieldFunctionNormalVector(3)));
    rElasticityTensor(4, 4) = mu * theta_new - (2. * mu * theta_new_b * (YieldFunctionNormalVector(4) * YieldFunctionNormalVector(4)));
    rElasticityTensor(4, 5) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(4) * YieldFunctionNormalVector(5)));

    rElasticityTensor(5, 0) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(5) * YieldFunctionNormalVector(0)));
    rElasticityTensor(5, 1) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(5) * YieldFunctionNormalVector(1)));
    rElasticityTensor(5, 2) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(5) * YieldFunctionNormalVector(2)));
    rElasticityTensor(5, 3) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(5) * YieldFunctionNormalVector(3)));
    rElasticityTensor(5, 4) = -(2. *mu * theta_new_b * (YieldFunctionNormalVector(5) * YieldFunctionNormalVector(4)));
    rElasticityTensor(5, 5) = mu * theta_new - (2. * mu * theta_new_b * (YieldFunctionNormalVector(5) * YieldFunctionNormalVector(5)));
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::GetLawFeatures(Features& rFeatures)
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

int LinearJ2Plasticity3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
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

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mInelasticFlag", mInelasticFlag);
    rSerializer.save("mPlasticStrain", mPlasticStrain);
    rSerializer.save("mPlasticStrainOld", mPlasticStrainOld);
    rSerializer.save("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
    rSerializer.save("mAccumulatedPlasticStrainOld", mAccumulatedPlasticStrainOld);
}

//************************************************************************************
//************************************************************************************

void LinearJ2Plasticity3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mInelasticFlag", mInelasticFlag);
    rSerializer.load("mPlasticStrain", mPlasticStrain);
    rSerializer.load("mPlasticStrainOld", mPlasticStrainOld);
    rSerializer.load("mAccumulatedPlasticStrain", mAccumulatedPlasticStrain);
    rSerializer.load("mAccumulatedPlasticStrainOld", mAccumulatedPlasticStrainOld);
}

} /* namespace Kratos.*/
