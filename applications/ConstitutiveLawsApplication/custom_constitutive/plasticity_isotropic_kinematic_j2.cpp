// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Fernando Rastellini
//  Collaborators:   Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                   Marcelo Raschi
//

// System includes

// External includes

// Project includes
#include "plasticity_isotropic_kinematic_j2.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

PlasticityIsotropicKinematicJ2::PlasticityIsotropicKinematicJ2()
    : ConstitutiveLaw()
{
}

/********************************COPY CONSTRUCTOR***********************************/
/***********************************************************************************/

PlasticityIsotropicKinematicJ2::PlasticityIsotropicKinematicJ2(const PlasticityIsotropicKinematicJ2& rOther)
    : ConstitutiveLaw(rOther),
      mPlasticStrain(rOther.mPlasticStrain),
      mEquivalentPlasticStrain(rOther.mEquivalentPlasticStrain)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer PlasticityIsotropicKinematicJ2::Clone() const
{
    return Kratos::make_shared<PlasticityIsotropicKinematicJ2>(*this);
}

/********************************DESTRUCTOR*****************************************/
/***********************************************************************************/

PlasticityIsotropicKinematicJ2::~PlasticityIsotropicKinematicJ2()
{
}

/***********************************************************************************/
/***********************************************************************************/

bool PlasticityIsotropicKinematicJ2::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == EQUIVALENT_PLASTIC_STRAIN){
        return true;
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

double& PlasticityIsotropicKinematicJ2::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == EQUIVALENT_PLASTIC_STRAIN){
         rValue = mEquivalentPlasticStrain;
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rThisVariable == EQUIVALENT_PLASTIC_STRAIN){
        mEquivalentPlasticStrain = rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::InitializeMaterial(
    const Properties&   rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector&       rShapeFunctionsValues)
{
    // Initialization when material CL is created (or cloned)
    mPlasticStrain = ZeroVector(this->GetStrainSize());
    mEquivalentPlasticStrain = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::FinalizeSolutionStep(
    const Properties&   rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector&       rShapeFunctionsValues,
    const ProcessInfo&  rCurrentProcessInfo)
{
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Initializations (if any) at each time step

    // NOTHING DONE HERE
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues); //PK1 not yet implemented
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues); //PK2 not yet implemented
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues); //Kirchhoff not yet implemented
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    BoundedArrayType plastic_strain;
    double           equivalent_plastic_strain;

    switch ( VoigtSize )
    {
    case 3:
        KRATOS_ERROR << "The constitutive law is not adapted for the following VoigtSize: " << VoigtSize << std::endl;
        //this->CalculateResponse3(rValues, plastic_strain, equivalent_plastic_strain);
        break;
    case 4:
        KRATOS_ERROR << "The constitutive law is not adapted for the following VoigtSize: " << VoigtSize << std::endl;
        //this->CalculateResponse4(rValues, plastic_strain, equivalent_plastic_strain);
        break;
    case 6:
        this->CalculateResponse6(rValues, plastic_strain, equivalent_plastic_strain);
        break;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::CalculateResponse6(
    ConstitutiveLaw::Parameters& rValues,
    BoundedArrayType& rPlasticStrain,
    double& rEquivalentPlasticStrain
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // Values of the CL
    Flags&        r_cl_options = rValues.GetOptions();  // The flags of the law
    const Vector& r_strain_vector = rValues.GetStrainVector();
    Vector&       r_stress_vector = rValues.GetStressVector();
    MatrixType&   r_tangent_tensor = rValues.GetConstitutiveMatrix();

    // Material properties
    const double young = r_material_properties[YOUNG_MODULUS];
    const double poisson = r_material_properties[POISSON_RATIO];
    const double isotropic_hardening = r_material_properties[ISOTROPIC_HARDENING_MODULUS];
    const double bulk = young / (3. - 6.*poisson); //bulk_modulus
    const double shear = young / ( 2. + 2.*poisson); //shear_modulus
    const double sqrt_two_thirds = std::sqrt(2.0/3.0); // sqrt_two_thirds = 0.81649658

    // Initialization
    rPlasticStrain = mPlasticStrain;                      //initialize value with last converged step
    rEquivalentPlasticStrain = mEquivalentPlasticStrain;  //initialize value with last converged step

    MatrixType elastic_tensor(VoigtSize,VoigtSize);
    CalculateElasticMatrix6(r_material_properties, elastic_tensor);  //compute elastic matrix

    BoundedArrayType sigma_trial;
    noalias(sigma_trial) = prod(elastic_tensor, r_strain_vector - rPlasticStrain);  //elastic trial prediction

    // Calculate deviatoric stress vector   (  deviatoric_stress = sigma - 1/3 tr(sigma) * I  )
    BoundedArrayType deviatoric_stress = sigma_trial;
    const double pressure = (sigma_trial[0] + sigma_trial[1] + sigma_trial[2]) /3.0;
    deviatoric_stress[0] -= pressure;
    deviatoric_stress[1] -= pressure;
    deviatoric_stress[2] -= pressure;

    // Euclidean norm of deviatoric stress
    const double norm_dev_stress = std::sqrt(deviatoric_stress[0]*deviatoric_stress[0] +
                                             deviatoric_stress[1]*deviatoric_stress[1] +
                                             deviatoric_stress[2]*deviatoric_stress[2] +
                                          2.*deviatoric_stress[3]*deviatoric_stress[3] +
                                          2.*deviatoric_stress[4]*deviatoric_stress[4] +
                                          2.*deviatoric_stress[5]*deviatoric_stress[5]);

    //evaluate yield function
    double yield_function = this->YieldFunction(norm_dev_stress, r_material_properties);

    if (yield_function <= 0.) {

        if( r_cl_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
            // Update the stresses
            r_stress_vector = sigma_trial;
        }

        if( r_cl_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
            // Update the tangent tensor
            r_tangent_tensor = elastic_tensor;
        }

    } else {

        const BoundedArrayType yield_normal = deviatoric_stress / norm_dev_stress;

        // Linear hardening
        // calculate plastic multiplier (consistency condition)
        const double delta_gamma = yield_function /
                                 (2.* shear* (1. + (isotropic_hardening / (3.* shear))));  //Simo equation 3.3.7

        // Update plastic strains
        rPlasticStrain[0] += delta_gamma*yield_normal[0];
        rPlasticStrain[1] += delta_gamma*yield_normal[1];
        rPlasticStrain[2] += delta_gamma*yield_normal[2];
        rPlasticStrain[3] += delta_gamma*yield_normal[3]*2;
        rPlasticStrain[4] += delta_gamma*yield_normal[4]*2;
        rPlasticStrain[5] += delta_gamma*yield_normal[5]*2;
        rEquivalentPlasticStrain += sqrt_two_thirds*delta_gamma;

        if( r_cl_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
            // Update stresses
            const double minus_two_shear_delta_gamma = -2.* shear * delta_gamma;
            r_stress_vector[0] = bulk * (r_strain_vector[0] + r_strain_vector[1] + r_strain_vector[2])
                            + deviatoric_stress[0] + minus_two_shear_delta_gamma * yield_normal[0];
            r_stress_vector[1] = bulk * (r_strain_vector[0] + r_strain_vector[1] + r_strain_vector[2])
                            + deviatoric_stress[1] + minus_two_shear_delta_gamma * yield_normal[1];
            r_stress_vector[2] = bulk * (r_strain_vector[0] + r_strain_vector[1] + r_strain_vector[2])
                            + deviatoric_stress[2] + minus_two_shear_delta_gamma * yield_normal[2];
            r_stress_vector[3] = deviatoric_stress[3] + minus_two_shear_delta_gamma * yield_normal[3];
            r_stress_vector[4] = deviatoric_stress[4] + minus_two_shear_delta_gamma * yield_normal[4];
            r_stress_vector[5] = deviatoric_stress[5] + minus_two_shear_delta_gamma * yield_normal[5];
        }


        if( r_cl_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
            // Update the tangent tensor
            CalculateTangentTensor6(delta_gamma, norm_dev_stress, yield_normal,
                                    r_material_properties, r_tangent_tensor);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

double& PlasticityIsotropicKinematicJ2::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
         const Variable<double>& rThisVariable,
                         double& rValue)
{
    if(rThisVariable == EQUIVALENT_PLASTIC_STRAIN){
        rValue = mEquivalentPlasticStrain;
    }
    return(rValue);
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);  // not yet implemented
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);  // not yet implemented
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);  // not yet implemented
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    BoundedArrayType plastic_strain;
    double           equivalent_plastic_strain;

    switch ( VoigtSize )
    {
    case 3:
        KRATOS_ERROR << "The constitutive law is not adapted for the following VoigtSize: " << VoigtSize << std::endl;
        //this->CalculateResponse3(rValues, plastic_strain, equivalent_plastic_strain);
        break;
    case 4:
        KRATOS_ERROR << "The constitutive law is not adapted for the following VoigtSize: " << VoigtSize << std::endl;
        //this->CalculateResponse4(rValues, plastic_strain, equivalent_plastic_strain);
        break;
    case 6:
        this->CalculateResponse6(rValues, plastic_strain, equivalent_plastic_strain);
        break;
    }

    // Update member variables
    mPlasticStrain = plastic_strain;
    mEquivalentPlasticStrain = equivalent_plastic_strain;
}

/***********************************************************************************/
/***********************************************************************************/

double PlasticityIsotropicKinematicJ2::YieldFunction(
    const double NormDeviatoricStress,
    const Properties& rMaterialProperties
    )
{
    const double sqrt_two_thirds = std::sqrt(2.0 / 3.0);
    const double yield_stress = rMaterialProperties[YIELD_STRESS];
    const double isotropic_hardening = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double k_last = yield_stress + isotropic_hardening * mEquivalentPlasticStrain;

    return NormDeviatoricStress - sqrt_two_thirds * k_last;
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::CalculateElasticMatrix6(
    const Properties& rMaterialProperties,
    MatrixType&  rElasticityTensor
    )
{
    // Material properties
    const double young = rMaterialProperties[YOUNG_MODULUS];
    const double poisson = rMaterialProperties[POISSON_RATIO];
    //const double bulk = young / (3. - 6.*poisson); //bulk_modulus
    const double shear = young / ( 2. + 2.*poisson); //shear_modulus
    const double two_shear = 2. * shear;
    //const double mu_lame = shear;
    const double lambda_lame = shear * poisson / (0.5 - poisson);

    if (rElasticityTensor.size1() != 6 || rElasticityTensor.size2() != 6)
        rElasticityTensor.resize(6, 6, false);
    rElasticityTensor.clear(); //make zero

    // Elasticity tensor (Simo, equation 2.7.11)
    rElasticityTensor(0, 0) = lambda_lame + two_shear;
    rElasticityTensor(0, 1) = lambda_lame;
    rElasticityTensor(0, 2) = lambda_lame;
    rElasticityTensor(1, 0) = lambda_lame;
    rElasticityTensor(1, 1) = lambda_lame + two_shear;
    rElasticityTensor(1, 2) = lambda_lame;
    rElasticityTensor(2, 0) = lambda_lame;
    rElasticityTensor(2, 1) = lambda_lame;
    rElasticityTensor(2, 2) = lambda_lame + two_shear;
    rElasticityTensor(3, 3) = shear;
    rElasticityTensor(4, 4) = shear;
    rElasticityTensor(5, 5) = shear;
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::CalculateTangentTensor6(
    const double DeltaGamma,
    const double NormDeviatoricStress,
    const BoundedArrayType& rYieldNormal,
    const Properties& rMaterialProperties,
    MatrixType& rTangentTensor
    )
{
    // Material properties
    const double young = rMaterialProperties[YOUNG_MODULUS];
    const double poisson = rMaterialProperties[POISSON_RATIO];
    const double isotropic_hardening = rMaterialProperties[ISOTROPIC_HARDENING_MODULUS];
    const double bulk = young / (3. - 6.*poisson); //bulk_modulus
    const double shear = young / ( 2. + 2.*poisson); //shear_modulus

    const double k_prima = isotropic_hardening;
    const double theta = 1. - (2. * shear * DeltaGamma) / NormDeviatoricStress;
    const double theta_bar = 1. / (1. + k_prima / (3. * shear)) - (1. - theta);
    const double two_shear_theta = 2. * shear * theta;
    const double minus_two_shear_theta_bar = -2. * shear * theta_bar;

    // Calculate tangent tensor (Simo box 3.2)
    rTangentTensor(0, 0) = bulk + 2./3. * two_shear_theta
                         + minus_two_shear_theta_bar * rYieldNormal[0] * rYieldNormal[0];
    rTangentTensor(0, 1) = bulk - two_shear_theta / 3.
                         + minus_two_shear_theta_bar * rYieldNormal[0] * rYieldNormal[1];
    rTangentTensor(0, 2) = bulk - 1./3. * two_shear_theta
                         + minus_two_shear_theta_bar * rYieldNormal[0] * rYieldNormal[2];
    rTangentTensor(0, 3) = minus_two_shear_theta_bar * rYieldNormal[0] * rYieldNormal[3];
    rTangentTensor(0, 4) = minus_two_shear_theta_bar * rYieldNormal[0] * rYieldNormal[4];
    rTangentTensor(0, 5) = minus_two_shear_theta_bar * rYieldNormal[0] * rYieldNormal[5];

    rTangentTensor(1, 0) = bulk - two_shear_theta / 3.
                         + minus_two_shear_theta_bar * rYieldNormal[1] * rYieldNormal[0];
    rTangentTensor(1, 1) = bulk + 2./3. * two_shear_theta
                         + minus_two_shear_theta_bar * rYieldNormal[1] * rYieldNormal[1];
    rTangentTensor(1, 2) = bulk - two_shear_theta / 3.
                         + minus_two_shear_theta_bar * rYieldNormal[1] * rYieldNormal[2];
    rTangentTensor(1, 3) = minus_two_shear_theta_bar * rYieldNormal[1] * rYieldNormal[3];
    rTangentTensor(1, 4) = minus_two_shear_theta_bar * rYieldNormal[1] * rYieldNormal[4];
    rTangentTensor(1, 5) = minus_two_shear_theta_bar * rYieldNormal[1] * rYieldNormal[5];

    rTangentTensor(2, 0) = bulk - two_shear_theta / 3.
                         + minus_two_shear_theta_bar * rYieldNormal[2] * rYieldNormal[0];
    rTangentTensor(2, 1) = bulk - two_shear_theta / 3.
                         + minus_two_shear_theta_bar * rYieldNormal[2] * rYieldNormal[1];
    rTangentTensor(2, 2) = bulk + 2./3. * two_shear_theta
                         + minus_two_shear_theta_bar * rYieldNormal[2] * rYieldNormal[2];
    rTangentTensor(2, 3) = minus_two_shear_theta_bar * rYieldNormal[2] * rYieldNormal[3];
    rTangentTensor(2, 4) = minus_two_shear_theta_bar * rYieldNormal[2] * rYieldNormal[4];
    rTangentTensor(2, 5) = minus_two_shear_theta_bar * rYieldNormal[2] * rYieldNormal[5];

    rTangentTensor(3, 0) = minus_two_shear_theta_bar * rYieldNormal[3] * rYieldNormal[0];
    rTangentTensor(3, 1) = minus_two_shear_theta_bar * rYieldNormal[3] * rYieldNormal[1];
    rTangentTensor(3, 2) = minus_two_shear_theta_bar * rYieldNormal[3] * rYieldNormal[2];
    rTangentTensor(3, 3) = two_shear_theta * 0.5
                         + minus_two_shear_theta_bar * rYieldNormal[3] * rYieldNormal[3];
    rTangentTensor(3, 4) = minus_two_shear_theta_bar * rYieldNormal[3] * rYieldNormal[4];
    rTangentTensor(3, 5) = minus_two_shear_theta_bar * rYieldNormal[3] * rYieldNormal[5];

    rTangentTensor(4, 0) = minus_two_shear_theta_bar * rYieldNormal[4] * rYieldNormal[0];
    rTangentTensor(4, 1) = minus_two_shear_theta_bar * rYieldNormal[4] * rYieldNormal[1];
    rTangentTensor(4, 2) = minus_two_shear_theta_bar * rYieldNormal[4] * rYieldNormal[2];
    rTangentTensor(4, 3) = minus_two_shear_theta_bar * rYieldNormal[4] * rYieldNormal[3];
    rTangentTensor(4, 4) = two_shear_theta * 0.5
                         + minus_two_shear_theta_bar * rYieldNormal[4] * rYieldNormal[4];
    rTangentTensor(4, 5) = minus_two_shear_theta_bar * rYieldNormal[4] * rYieldNormal[5];

    rTangentTensor(5, 0) = minus_two_shear_theta_bar * rYieldNormal[5] * rYieldNormal[0];
    rTangentTensor(5, 1) = minus_two_shear_theta_bar * rYieldNormal[5] * rYieldNormal[1];
    rTangentTensor(5, 2) = minus_two_shear_theta_bar * rYieldNormal[5] * rYieldNormal[2];
    rTangentTensor(5, 3) = minus_two_shear_theta_bar * rYieldNormal[5] * rYieldNormal[3];
    rTangentTensor(5, 4) = minus_two_shear_theta_bar * rYieldNormal[5] * rYieldNormal[4];
    rTangentTensor(5, 5) = two_shear_theta * 0.5
                         + minus_two_shear_theta_bar * rYieldNormal[5] * rYieldNormal[5];
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = VoigtSize;
    rFeatures.mSpaceDimension = Dimension;
}

/***********************************************************************************/
/***********************************************************************************/

int PlasticityIsotropicKinematicJ2::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS));
    KRATOS_CHECK(rMaterialProperties.Has(ISOTROPIC_HARDENING_MODULUS));

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mPlasticStrain", mPlasticStrain);
    rSerializer.save("mEquivalentPlasticStrain", mEquivalentPlasticStrain);
}

/***********************************************************************************/
/***********************************************************************************/

void PlasticityIsotropicKinematicJ2::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mPlasticStrain", mPlasticStrain);
    rSerializer.load("mEquivalentPlasticStrain", mEquivalentPlasticStrain);
}

} /* namespace Kratos.*/
