#include "linear_j2_plasticity_3d.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearJ2Plasticity3D::LinearJ2Plasticity3D()
    : ConstitutiveLaw()
{
}

//********************************COPY CONSTRUCTOR***********************************************
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

//********************************DESTRUCTOR***********************************************
//************************************************************************************

LinearJ2Plasticity3D::~LinearJ2Plasticity3D()
{
}

bool LinearJ2Plasticity3D::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        return true;
    }
    if(rThisVariable == INELASTIC_FLAG){
        return true;
    }
    return false;
}

double& LinearJ2Plasticity3D::GetValue(const Variable<double>& rThisVariable,
                                                        double& rValue)
{
    if(rThisVariable == STRAIN_ENERGY){
        rValue = mStrainEnergy;
        //TODO (marcelo); fix Bbar element to use CalculateValue to get StrainEnergy value
    }
    if(rThisVariable == INELASTIC_FLAG){
        rValue = mInelasticFlag;
    }
    return rValue;
}

void LinearJ2Plasticity3D::InitializeMaterial(const Properties& material_prop,
                                              const GeometryType& rElementGeometry,
                                              const Vector& rShapeFunctionsValues)
{
    mPlasticStrainOld = ZeroVector(this->GetStrainSize());
    mAccumulatedPlasticStrainOld = 0.0;
}

void LinearJ2Plasticity3D::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    mPlasticStrainOld = mPlasticStrain;
    mAccumulatedPlasticStrainOld = mAccumulatedPlasticStrain;
}

void LinearJ2Plasticity3D::CalculateMaterialResponsePK1(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearJ2Plasticity3D::CalculateMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearJ2Plasticity3D::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void LinearJ2Plasticity3D::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    Flags &Options = rValues.GetOptions();

    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    Vector& StrainVector = rValues.GetStrainVector();
    Vector& StressVector = rValues.GetStressVector();


    if (Options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    }

    if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(ConstitutiveMatrix, rMaterialProperties);
    }

    if (Options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        if (rValues.GetProcessInfo().Has(INITIAL_STRAIN)) {
            noalias(StrainVector) += rValues.GetProcessInfo()[INITIAL_STRAIN];
        }
        Matrix ElasticTensor;
        Matrix& TangentTensor = rValues.GetConstitutiveMatrix();
        double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
        double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
        double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
        double trial_yield_function;
        double E = rMaterialProperties[YOUNG_MODULUS];
        double poisson_ratio = rMaterialProperties[POISSON_RATIO];
        double mu = E / (2. + 2. * poisson_ratio);
        double volumetric_modulus = E / (3. * (1. - 2. * poisson_ratio));

        mPlasticStrain = mPlasticStrainOld;
        mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld;

        ElasticTensor.resize(6, 6);
        CalculateElasticMatrix(ElasticTensor, rMaterialProperties);
        Vector sigma_trial;
        sigma_trial.resize(6);
        sigma_trial = prod(ElasticTensor, StrainVector - mPlasticStrainOld);

        // StressTrialDev = sigma - 1/3 tr(sigma) * I
        Vector StressTrialDev;
        StressTrialDev.resize(6);
        StressTrialDev = sigma_trial;

        double trace =
            0.33333333333333333 * (sigma_trial(0) + sigma_trial(1) + sigma_trial(2));
        StressTrialDev(0) -= trace;
        StressTrialDev(1) -= trace;
        StressTrialDev(2) -= trace;
        double norm_dev_stress = std::sqrt(StressTrialDev(0) * StressTrialDev(0) +
                                           StressTrialDev(1) * StressTrialDev(1) +
                                           StressTrialDev(2) * StressTrialDev(2) +
                                           2. * StressTrialDev(3) * StressTrialDev(3) +
                                           2. * StressTrialDev(4) * StressTrialDev(4) +
                                           2. * StressTrialDev(5) * StressTrialDev(5));
        trial_yield_function = this->yieldFunction(norm_dev_stress, rMaterialProperties);
        double dgamma = 0;

        if (trial_yield_function <= 0.) {
            // ELASTIC
            mInelasticFlag = 0;
            StressVector = sigma_trial;
            TangentTensor = ElasticTensor;
        }
        else {
            // INELASTIC
            mInelasticFlag = 1;
            Vector YieldFunctionNormalVector = StressTrialDev / norm_dev_stress;
            if (delta_k != 0.0 && hardening_exponent != 0.0) {
                // Exponential softening
                dgamma = GetDeltaGamma(norm_dev_stress, rMaterialProperties);
            }
            else {
                // Linear softening
                dgamma = trial_yield_function /
                         (2. * mu * (1. + (hardening_modulus / (3. * mu))));
            }

            StressVector(0) =
                volumetric_modulus * (StrainVector(0) + StrainVector(1) + StrainVector(2)) +
                StressTrialDev(0) - 2. * mu * dgamma * YieldFunctionNormalVector(0);
            StressVector(1) =
                volumetric_modulus * (StrainVector(0) + StrainVector(1) + StrainVector(2)) +
                StressTrialDev(1) - 2. * mu * dgamma * YieldFunctionNormalVector(1);
            StressVector(2) =
                volumetric_modulus * (StrainVector(0) + StrainVector(1) + StrainVector(2)) +
                StressTrialDev(2) - 2. * mu * dgamma * YieldFunctionNormalVector(2);
            StressVector(3) =
                StressTrialDev(3) - 2. * mu * dgamma * YieldFunctionNormalVector(3);
            StressVector(4) =
                StressTrialDev(4) - 2. * mu * dgamma * YieldFunctionNormalVector(4);
            StressVector(5) =
                StressTrialDev(5) - 2. * mu * dgamma * YieldFunctionNormalVector(5);

            mPlasticStrain(0) = mPlasticStrainOld(0) + dgamma * YieldFunctionNormalVector(0);
            mPlasticStrain(1) = mPlasticStrainOld(1) + dgamma * YieldFunctionNormalVector(1);
            mPlasticStrain(2) = mPlasticStrainOld(2) + dgamma * YieldFunctionNormalVector(2);
            mPlasticStrain(3) = mPlasticStrainOld(3) + dgamma * YieldFunctionNormalVector(3) * 2;
            mPlasticStrain(4) = mPlasticStrainOld(4) + dgamma * YieldFunctionNormalVector(4) * 2;
            mPlasticStrain(5) = mPlasticStrainOld(5) + dgamma * YieldFunctionNormalVector(5) * 2;
            mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld + 0.8164965809277260 * dgamma;

            // Update derivative of the hardening-softening modulus

            CalculateTangentTensor(dgamma, norm_dev_stress, YieldFunctionNormalVector, rMaterialProperties, TangentTensor);
        }

    // Linear + exponential hardening
    mStrainEnergy = 0.5 * inner_prod(StrainVector - mPlasticStrain, prod(ElasticTensor, StrainVector - mPlasticStrain))
                    + GetPlasticPotential(rMaterialProperties);
    }
}

double& LinearJ2Plasticity3D::CalculateValue(Parameters& rParameterValues,
        const Variable<double>& rThisVariable, double& rValue) {

    //const Properties& MaterialProperties  = rParameterValues.GetMaterialProperties();
    //Vector& StrainVector                  = rParameterValues.GetStrainVector();
    //Vector& StressVector                  = rParameterValues.GetStressVector();
    //const double& E          = MaterialProperties[YOUNG_MODULUS];
    //const double& NU    = MaterialProperties[POISSON_RATIO];

    //if (rThisVariable == STRAIN_ENERGY)
    //{
    //    CalculateCauchyGreenStrain(rParameterValues, StrainVector);
    //    CalculatePK2Stress( StrainVector, StressVector, E, NU );

    //    rValue = 0.5 * inner_prod(StrainVector,StressVector); // Strain energy = 0.5*E:C:E
    //}

    if(rThisVariable == STRAIN_ENERGY){
        rValue = mStrainEnergy;
        KRATOS_WATCH("CALCULATE")
    }
    return( rValue );
}

void LinearJ2Plasticity3D::FinalizeMaterialResponsePK1(Parameters& rValues)
{
}

void LinearJ2Plasticity3D::FinalizeMaterialResponsePK2(Parameters& rValues)
{
}

void LinearJ2Plasticity3D::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
}

void LinearJ2Plasticity3D::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
}

double LinearJ2Plasticity3D::GetSaturationHardening(const Properties& rMaterialProperties)
{
    double yield_stress = rMaterialProperties[YIELD_STRESS];
    double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];
    //mAccumulatedPlasticStrain = mAccumulatedPlasticStrainOld;
        double k_new =
            yield_stress + (theta * hardening_modulus * mAccumulatedPlasticStrain) +
                delta_k * (1. - std::exp(-hardening_exponent * mAccumulatedPlasticStrain));
    return k_new;
}

double LinearJ2Plasticity3D::GetPlasticPotential(const Properties& rMaterialProperties)
{
    double theta = rMaterialProperties[REFERENCE_HARDENING_MODULUS];
    double hardening_modulus = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    double delta_k = rMaterialProperties[INFINITY_HARDENING_MODULUS];
    double hardening_exponent = rMaterialProperties[HARDENING_EXPONENT];

    double Wp_new = 0.5*(theta * hardening_modulus * std::pow(mAccumulatedPlasticStrain, 2.0)) +
                    delta_k * (mAccumulatedPlasticStrain -
                            (1/hardening_exponent) * (1- std::exp(-hardening_exponent * mAccumulatedPlasticStrain)));
    return Wp_new;
}

double LinearJ2Plasticity3D::GetDeltaGamma(double norm_s_trial,
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

double LinearJ2Plasticity3D::yieldFunction(const double norm_dev_stress,
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

void LinearJ2Plasticity3D::CalculateElasticMatrix(Matrix &D, const Properties &props)
{
    double E = props[YOUNG_MODULUS];
    double poisson_ratio = props[POISSON_RATIO];
    double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    double mu = E / (2. + 2. * poisson_ratio);

    D.clear();
    D.resize(6, 6, false);
    D = ZeroMatrix(6, 6);

    D(0, 0) = lambda + 2. * mu;
    D(0, 1) = lambda;
    D(0, 2) = lambda;
    D(1, 0) = lambda;
    D(1, 1) = lambda + 2. * mu;
    D(1, 2) = lambda;
    D(2, 0) = lambda;
    D(2, 1) = lambda;
    D(2, 2) = lambda + 2. * mu;
    D(3, 3) = mu;
    D(4, 4) = mu;
    D(5, 5) = mu;
}

void LinearJ2Plasticity3D::CalculateTangentTensor(
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
    D(0, 4) = -(2 * mu * theta_new_b * (N_new(0) * N_new(4)));
    D(0, 5) = -(2 * mu * theta_new_b * (N_new(0) * N_new(5)));

    D(1, 0) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
              (2 * mu * theta_new_b * (N_new(1) * N_new(0)));
    D(1, 1) = volumetric_modulus + (2 * mu * theta_new * 2. / 3.) -
              (2 * mu * theta_new_b * (N_new(1) * N_new(1)));
    D(1, 2) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
              (2 * mu * theta_new_b * (N_new(1) * N_new(2)));
    D(1, 3) = -(2 * mu * theta_new_b * (N_new(1) * N_new(3)));
    D(1, 4) = -(2 * mu * theta_new_b * (N_new(1) * N_new(4)));
    D(1, 5) = -(2 * mu * theta_new_b * (N_new(1) * N_new(5)));

    D(2, 0) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
              (2 * mu * theta_new_b * (N_new(2) * N_new(0)));
    D(2, 1) = volumetric_modulus + (2 * mu * theta_new * (-1. / 3.)) -
              (2 * mu * theta_new_b * (N_new(2) * N_new(1)));
    D(2, 2) = volumetric_modulus + (2 * mu * theta_new * 2. / 3.) -
              (2 * mu * theta_new_b * (N_new(2) * N_new(2)));
    D(2, 3) = -(2 * mu * theta_new_b * (N_new(2) * N_new(3)));
    D(2, 4) = -(2 * mu * theta_new_b * (N_new(2) * N_new(4)));
    D(2, 5) = -(2 * mu * theta_new_b * (N_new(2) * N_new(5)));

    D(3, 0) = -(2 * mu * theta_new_b * (N_new(3) * N_new(0)));
    D(3, 1) = -(2 * mu * theta_new_b * (N_new(3) * N_new(1)));
    D(3, 2) = -(2 * mu * theta_new_b * (N_new(3) * N_new(2)));
    D(3, 3) = mu * theta_new - (2 * mu * theta_new_b * (N_new(3) * N_new(3)));
    D(3, 4) = -(2 * mu * theta_new_b * (N_new(3) * N_new(4)));
    D(3, 5) = -(2 * mu * theta_new_b * (N_new(3) * N_new(5)));

    D(4, 0) = -(2 * mu * theta_new_b * (N_new(4) * N_new(0)));
    D(4, 1) = -(2 * mu * theta_new_b * (N_new(4) * N_new(1)));
    D(4, 2) = -(2 * mu * theta_new_b * (N_new(4) * N_new(2)));
    D(4, 3) = -(2 * mu * theta_new_b * (N_new(4) * N_new(3)));
    D(4, 4) = mu * theta_new - (2 * mu * theta_new_b * (N_new(4) * N_new(4)));
    D(4, 5) = -(2 * mu * theta_new_b * (N_new(4) * N_new(5)));

    D(5, 0) = -(2 * mu * theta_new_b * (N_new(5) * N_new(0)));
    D(5, 1) = -(2 * mu * theta_new_b * (N_new(5) * N_new(1)));
    D(5, 2) = -(2 * mu * theta_new_b * (N_new(5) * N_new(2)));
    D(5, 3) = -(2 * mu * theta_new_b * (N_new(5) * N_new(3)));
    D(5, 4) = -(2 * mu * theta_new_b * (N_new(5) * N_new(4)));
    D(5, 5) = mu * theta_new - (2 * mu * theta_new_b * (N_new(5) * N_new(5)));
}

void LinearJ2Plasticity3D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = 6;
    rFeatures.mSpaceDimension = 3;
}

int LinearJ2Plasticity3D::Check(const Properties& rMaterialProperties,
                                                 const GeometryType& rElementGeometry,
                                                 const ProcessInfo& rCurrentProcessInfo)
{
    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS] <= 0.0) {
        KRATOS_ERROR << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;
    }

    const double& nu = rMaterialProperties[POISSON_RATIO];
    const bool check = bool( (nu >0.499 && nu<0.501 ) || (nu < -0.999 && nu > -1.01 ) );
    if(POISSON_RATIO.Key() == 0 || check) {
        KRATOS_ERROR << "POISSON_RATIO has Key zero or invalid value " << std::endl;
    }

    if(YIELD_STRESS.Key() == 0 || rMaterialProperties[YIELD_STRESS] <= 0.0) {
        KRATOS_ERROR << "YIELD_STRESS has Key zero or invalid value " << std::endl;
    }

    if(REFERENCE_HARDENING_MODULUS.Key() == 0) {
        KRATOS_ERROR << "REFERENCE_HARDENING_MODULUS has Key zero" << std::endl;
    }

    if(ISOTROPIC_HARDENING_MODULUS.Key() == 0) {
        KRATOS_ERROR << "ISOTROPIC_HARDENING_MODULUS has Key zero" << std::endl;
    }

    if(INFINITY_HARDENING_MODULUS.Key() == 0) {
        KRATOS_ERROR << "INFINITY_HARDENING_MODULUS has Key zero" << std::endl;
    }

    if(HARDENING_EXPONENT.Key() == 0) {
        KRATOS_ERROR << "HARDENING_EXPONENT has Key zero" << std::endl;
    }

    return 0;
}
} /* namespace Kratos.*/
