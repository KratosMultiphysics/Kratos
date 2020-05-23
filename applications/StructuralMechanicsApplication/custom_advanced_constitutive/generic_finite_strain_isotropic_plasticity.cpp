// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//  Collaborator:    Alejandro Cornejo
//                   Lucia Barbu
//

// System includes

// Project includes
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_advanced_constitutive/generic_finite_strain_isotropic_plasticity.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_finite_strain_constitutive_law_integrator_plasticity.h"

// Hyperelastic behaviours
#include "custom_advanced_constitutive/hyper_elastic_isotropic_kirchhoff_3d.h"
#include "custom_advanced_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"

// Yield surfaces
#include "custom_advanced_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_advanced_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{
template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponsePK2(rValues);

    Vector& stress_vector                = rValues.GetStressVector();
    const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
    const double determinant_f           = rValues.GetDeterminantF();

    this->TransformStresses(stress_vector, deformation_gradient_f, determinant_f, ConstitutiveLaw::StressMeasure_PK2, ConstitutiveLaw::StressMeasure_PK1);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Auxiliar values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // We check the current step and NL iteration
    const ProcessInfo& r_current_process_info = rValues.GetProcessInfo();
    const bool first_computation = (r_current_process_info[NL_ITERATION_NUMBER] == 1 && r_current_process_info[STEP] == 1) ? true : false;

    if (first_computation) { // First computation always pure elastic for elemements not providing the strain
        if( r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
            this->CalculateValue(rValues, GREEN_LAGRANGE_STRAIN_VECTOR, r_strain_vector);
        }

        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Vector& r_stress_vector = rValues.GetStressVector();
            // We backup the deformation gradient
            const double& r_det_deformation_gradient_backup = rValues.GetDeterminantF();
            const Matrix& r_deformation_gradient_backup = rValues.GetDeformationGradientF();
            const Matrix& r_plastic_deformation_gradient = this->GetPlasticDeformationGradient();

            // We compute the elastic deformation gradient  Fe = plastic_indicator * inv(Fp)
            Matrix inverse_F_p ( Dimension, Dimension );
            double aux_det_Fp = 0;
            MathUtils<double>::InvertMatrix( r_plastic_deformation_gradient, inverse_F_p, aux_det_Fp);
            const Matrix elastic_deformation_gradient = prod(r_deformation_gradient_backup, inverse_F_p);

            rValues.SetDeterminantF(MathUtils<double>::DetMat(elastic_deformation_gradient));
            rValues.SetDeformationGradientF(elastic_deformation_gradient);
            Vector auxiliar_predictive_stress_vector;
            this->CalculateValue(rValues, PK2_STRESS_VECTOR, auxiliar_predictive_stress_vector);
            for (std::size_t i_voigt = 0; i_voigt < VoigtSize; ++i_voigt) {
                r_stress_vector[i_voigt] =  auxiliar_predictive_stress_vector[i_voigt];
            }

            // We revert the deformation gradient
            rValues.SetDeterminantF(r_det_deformation_gradient_backup);
            rValues.SetDeformationGradientF(r_deformation_gradient_backup);

            if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                this->CalculateValue(rValues, CONSTITUTIVE_MATRIX_PK2, r_constitutive_matrix);
            }
        }
    } else { // We check for plasticity
        // Integrate Stress plasticity
        Vector& integrated_stress_vector = rValues.GetStressVector();
        const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());

        this->CalculateValue(rValues, GREEN_LAGRANGE_STRAIN_VECTOR, r_strain_vector);

        // We compute the stress or the constitutive matrix
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {

            // We get some variables
            double threshold = this->GetThreshold();
            double plastic_dissipation = this->GetPlasticDissipation();
            Vector plastic_strain;
            Matrix plastic_deformation_gradient = this->GetPlasticDeformationGradient();

            // We backup the deformation gradient
            const double& r_det_deformation_gradient_backup = rValues.GetDeterminantF();
            const Matrix& r_deformation_gradient_backup = rValues.GetDeformationGradientF();

            // We compute the predicted stress vector
            BoundedArrayType predictive_stress_vector;
            if( r_constitutive_law_options.Is( ConstitutiveLaw::U_P_LAW ) ) {
                predictive_stress_vector = rValues.GetStressVector();
            } else {
                // We compute the elastic deformation gradient  Fe = plastic_indicator * inv(Fp)
                Matrix inverse_F_p ( Dimension, Dimension );
                double aux_det_Fp = 0;

                MathUtils<double>::InvertMatrix( plastic_deformation_gradient, inverse_F_p, aux_det_Fp);
                const Matrix elastic_deformation_gradient = prod(r_deformation_gradient_backup, inverse_F_p);

                rValues.SetDeterminantF(MathUtils<double>::DetMat(elastic_deformation_gradient));
                rValues.SetDeformationGradientF(elastic_deformation_gradient);
                Vector auxiliar_predictive_stress_vector;
                this->CalculateValue(rValues, PK2_STRESS_VECTOR, auxiliar_predictive_stress_vector);
                for (std::size_t i_voigt = 0; i_voigt < VoigtSize; ++i_voigt) {
                    predictive_stress_vector[i_voigt] =  auxiliar_predictive_stress_vector[i_voigt];
                }

                // We revert the deformation gradient
                rValues.SetDeterminantF(r_det_deformation_gradient_backup);
                rValues.SetDeformationGradientF(r_deformation_gradient_backup);

                // Check
                KRATOS_ERROR_IF(r_det_deformation_gradient_backup < std::numeric_limits<double>::epsilon()) << "Deformation gradient determinant (detF) < 0.0 : " << r_det_deformation_gradient_backup << std::endl;
            }

            // We compute the plastic strain
            const double& det_f = rValues.GetDeterminantF();
            const Matrix& deformation_gradient = rValues.GetDeformationGradientF();
            rValues.SetDeterminantF(MathUtils<double>::DetMat(plastic_deformation_gradient));
            rValues.SetDeformationGradientF(plastic_deformation_gradient);
            this->CalculateValue(rValues, GREEN_LAGRANGE_STRAIN_VECTOR, plastic_strain);
            rValues.SetDeterminantF(det_f);
            rValues.SetDeformationGradientF(deformation_gradient);

            // Check
            KRATOS_ERROR_IF(det_f < std::numeric_limits<double>::epsilon()) << "Deformation gradient determinant (detF) < 0.0 : " << det_f << std::endl;

            // Initialize Plastic Parameters
            double uniaxial_stress = 0.0, plastic_denominator = 0.0;
            BoundedArrayType yield_surface_derivative = ZeroVector(VoigtSize);     // DF/DS
            BoundedArrayType plastic_potential_derivative = ZeroVector(VoigtSize); // DG/DS
            const BoundedArrayType dummy_plastic_strain_increment = ZeroVector(VoigtSize);

            // Elastic Matrix
            this->CalculateValue(rValues, CONSTITUTIVE_MATRIX_PK2, r_constitutive_matrix);

            // Compute the plastic parameters
            const double plastic_indicator = TConstLawIntegratorType::CalculatePlasticParameters(
                predictive_stress_vector, r_strain_vector, uniaxial_stress,
                threshold, plastic_denominator, yield_surface_derivative, plastic_potential_derivative,
                plastic_dissipation, dummy_plastic_strain_increment,
                r_constitutive_matrix, rValues, characteristic_length,
                plastic_strain);

            if (plastic_indicator <= std::abs(1.0e-4 * threshold)) { // Elastic case
                noalias(integrated_stress_vector) = predictive_stress_vector;
            } else { // Plastic case
                // While loop backward euler
                /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
                TConstLawIntegratorType::IntegrateStressVector(
                    *this, GREEN_LAGRANGE_STRAIN_VECTOR, PK2_STRESS_VECTOR,
                    predictive_stress_vector, r_strain_vector, uniaxial_stress, threshold,
                    plastic_denominator, yield_surface_derivative, plastic_potential_derivative,
                    plastic_dissipation, plastic_deformation_gradient, mPreviousDeformationGradient,
                    r_constitutive_matrix, rValues, characteristic_length);
                noalias(integrated_stress_vector) = predictive_stress_vector;

                if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                    this->CalculateTangentTensor(rValues, ConstitutiveLaw::StressMeasure_PK2); // this modifies the ConstitutiveMatrix
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Auxiliar values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // We check the current step and NL iteration
    const ProcessInfo& r_current_process_info = rValues.GetProcessInfo();
    const bool first_computation = (r_current_process_info[NL_ITERATION_NUMBER] == 1 && r_current_process_info[STEP] == 1) ? true : false;

    if (first_computation) { // First computation always pure elastic for elemements not providing the strain
        if( r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
            this->CalculateValue(rValues, ALMANSI_STRAIN_VECTOR, r_strain_vector);
        }

        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Vector& r_stress_vector = rValues.GetStressVector();
            // We backup the deformation gradient
            const double& r_det_deformation_gradient_backup = rValues.GetDeterminantF();
            const Matrix& r_deformation_gradient_backup = rValues.GetDeformationGradientF();
            const Matrix& r_plastic_deformation_gradient = this->GetPlasticDeformationGradient();

            // We compute the elastic deformation gradient  Fe = plastic_indicator * inv(Fp)
            Matrix inverse_F_p ( Dimension, Dimension );
            double aux_det_Fp = 0;
            MathUtils<double>::InvertMatrix( r_plastic_deformation_gradient, inverse_F_p, aux_det_Fp);
            const Matrix elastic_deformation_gradient = prod(r_deformation_gradient_backup, inverse_F_p);

            rValues.SetDeterminantF(MathUtils<double>::DetMat(elastic_deformation_gradient));
            rValues.SetDeformationGradientF(elastic_deformation_gradient);
            Vector auxiliar_predictive_stress_vector;
            this->CalculateValue(rValues, KIRCHHOFF_STRESS_VECTOR, auxiliar_predictive_stress_vector);
            for (std::size_t i_voigt = 0; i_voigt < VoigtSize; ++i_voigt) {
                r_stress_vector[i_voigt] =  auxiliar_predictive_stress_vector[i_voigt];
            }

            // We revert the deformation gradient
            rValues.SetDeterminantF(r_det_deformation_gradient_backup);
            rValues.SetDeformationGradientF(r_deformation_gradient_backup);

            if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                this->CalculateValue(rValues, CONSTITUTIVE_MATRIX_KIRCHHOFF, r_constitutive_matrix);
            }
        }
    } else { // We check for plasticity
        // Integrate Stress plasticity
        Vector& integrated_stress_vector = rValues.GetStressVector();
        const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());

        if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            this->CalculateValue(rValues, ALMANSI_STRAIN_VECTOR, r_strain_vector);
        }

        // We compute the stress or the constitutive matrix
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {

            // We get some variables
            double threshold = this->GetThreshold();
            double plastic_dissipation = this->GetPlasticDissipation();
            Vector plastic_strain;
            Matrix plastic_deformation_gradient = this->GetPlasticDeformationGradient();

            // We backup the deformation gradient
            const double& r_det_deformation_gradient_backup = rValues.GetDeterminantF();
            const Matrix& r_deformation_gradient_backup = rValues.GetDeformationGradientF();

            // We compute the predicted stress vector
            BoundedArrayType predictive_stress_vector;
            if( r_constitutive_law_options.Is( ConstitutiveLaw::U_P_LAW ) ) {
                predictive_stress_vector = rValues.GetStressVector();
            } else {
                // We compute the elastic deformation gradient  Fe = plastic_indicator * inv(Fp)
                Matrix inverse_F_p ( Dimension, Dimension );
                double aux_det_Fp = 0;
                MathUtils<double>::InvertMatrix( plastic_deformation_gradient, inverse_F_p, aux_det_Fp);
                const Matrix elastic_deformation_gradient = prod(r_deformation_gradient_backup, inverse_F_p);

                rValues.SetDeterminantF(MathUtils<double>::DetMat(elastic_deformation_gradient));
                rValues.SetDeformationGradientF(elastic_deformation_gradient);
                Vector auxiliar_predictive_stress_vector;
                this->CalculateValue(rValues, KIRCHHOFF_STRESS_VECTOR, auxiliar_predictive_stress_vector);
                for (std::size_t i_voigt = 0; i_voigt < VoigtSize; ++i_voigt) {
                    predictive_stress_vector[i_voigt] =  auxiliar_predictive_stress_vector[i_voigt];
                }

                // We revert the deformation gradient
                rValues.SetDeterminantF(r_det_deformation_gradient_backup);
                rValues.SetDeformationGradientF(r_deformation_gradient_backup);
            }

            // We compute the plastic strain
            const double det_f = rValues.GetDeterminantF();
            const Matrix deformation_gradient = rValues.GetDeformationGradientF();
            rValues.SetDeterminantF(MathUtils<double>::DetMat(plastic_deformation_gradient));
            rValues.SetDeformationGradientF(plastic_deformation_gradient);
            this->CalculateValue(rValues, ALMANSI_STRAIN_VECTOR, plastic_strain);
            rValues.SetDeterminantF(det_f);
            rValues.SetDeformationGradientF(deformation_gradient);

            // Initialize Plastic Parameters
            double uniaxial_stress = 0.0, plastic_denominator = 0.0;
            BoundedArrayType yield_surface_derivative = ZeroVector(VoigtSize);     // DF/DS
            BoundedArrayType plastic_potential_derivative = ZeroVector(VoigtSize); // DG/DS
            const BoundedArrayType dummy_plastic_strain_increment = ZeroVector(VoigtSize);

            // Elastic Matrix
            this->CalculateValue(rValues, CONSTITUTIVE_MATRIX_KIRCHHOFF, r_constitutive_matrix);

            // Compute the plastic parameters
            const double plastic_indicator = TConstLawIntegratorType::CalculatePlasticParameters(
                predictive_stress_vector, r_strain_vector, uniaxial_stress,
                threshold, plastic_denominator, yield_surface_derivative, plastic_potential_derivative,
                plastic_dissipation, dummy_plastic_strain_increment,
                r_constitutive_matrix, rValues, characteristic_length,
                plastic_strain);

            if (plastic_indicator <= std::abs(1.0e-4 * threshold)) { // Elastic case
                noalias(integrated_stress_vector) = predictive_stress_vector;
            } else { // Plastic case
                // While loop backward euler
                /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
                TConstLawIntegratorType::IntegrateStressVector(
                    *this, ALMANSI_STRAIN_VECTOR, KIRCHHOFF_STRESS_VECTOR,
                    predictive_stress_vector, r_strain_vector, uniaxial_stress, threshold,
                    plastic_denominator, yield_surface_derivative, plastic_potential_derivative,
                    plastic_dissipation, plastic_deformation_gradient, mPreviousDeformationGradient,
                    r_constitutive_matrix, rValues, characteristic_length);
                noalias(integrated_stress_vector) = predictive_stress_vector;

                if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                    this->CalculateTangentTensor(rValues, ConstitutiveLaw::StressMeasure_Kirchhoff); // this modifies the ConstitutiveMatrix
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseKirchhoff(rValues);

    Vector& stress_vector       = rValues.GetStressVector();
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
    const double determinant_f  = rValues.GetDeterminantF();

    // Set to Cauchy Stress:
    stress_vector       /= determinant_f;
    constitutive_matrix /= determinant_f;
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::CalculateTangentTensor(
    ConstitutiveLaw::Parameters& rValues,
    const ConstitutiveLaw::StressMeasure& rStressMeasure
    )
{
    // Calculates the Tangent Constitutive Tensor by perturbation
    TangentOperatorCalculatorUtility::CalculateTangentTensorFiniteDeformation(rValues, this, rStressMeasure);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We construct the CL parameters
    ProcessInfo dummy_process_info;
    ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties, dummy_process_info);

    // We call the integrator
    double initial_threshold;
    TConstLawIntegratorType::GetInitialUniaxialThreshold(aux_param, initial_threshold);
    this->SetThreshold(initial_threshold);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Auxiliar values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // The characteristic length
    const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());

    this->CalculateValue(rValues, GREEN_LAGRANGE_STRAIN_VECTOR, r_strain_vector);

    // We get some variables
    double threshold = this->GetThreshold();
    double plastic_dissipation = this->GetPlasticDissipation();
    Vector plastic_strain;
    Matrix plastic_deformation_gradient = this->GetPlasticDeformationGradient();

    // We backup the deformation gradient
    const double& r_det_deformation_gradient_backup = rValues.GetDeterminantF();
    const Matrix& r_deformation_gradient_backup = rValues.GetDeformationGradientF();

    // We compute the predicted stress vector
    BoundedArrayType predictive_stress_vector;
    if( r_constitutive_law_options.Is( ConstitutiveLaw::U_P_LAW ) ) {
        predictive_stress_vector = rValues.GetStressVector();
    } else {
        // We compute the elastic deformation gradient  Fe = plastic_indicator * inv(Fp)
        Matrix inverse_F_p ( Dimension, Dimension );
        double aux_det_Fp = 0;

        MathUtils<double>::InvertMatrix( plastic_deformation_gradient, inverse_F_p, aux_det_Fp);
        const Matrix elastic_deformation_gradient = prod(r_deformation_gradient_backup, inverse_F_p);

        rValues.SetDeterminantF(MathUtils<double>::DetMat(elastic_deformation_gradient));
        rValues.SetDeformationGradientF(elastic_deformation_gradient);
        Vector auxiliar_predictive_stress_vector;
        this->CalculateValue(rValues, PK2_STRESS_VECTOR, auxiliar_predictive_stress_vector);
        for (std::size_t i_voigt = 0; i_voigt < VoigtSize; ++i_voigt) {
            predictive_stress_vector[i_voigt] =  auxiliar_predictive_stress_vector[i_voigt];
        }

        // We revert the deformation gradient
        rValues.SetDeterminantF(r_det_deformation_gradient_backup);
        rValues.SetDeformationGradientF(r_deformation_gradient_backup);

        // Check
        KRATOS_ERROR_IF(r_det_deformation_gradient_backup < std::numeric_limits<double>::epsilon()) << "Deformation gradient determinant (detF) < 0.0 : " << r_det_deformation_gradient_backup << std::endl;
    }

    // We compute the plastic strain
    const double& det_f = rValues.GetDeterminantF();
    const Matrix& deformation_gradient = rValues.GetDeformationGradientF();
    rValues.SetDeterminantF(MathUtils<double>::DetMat(plastic_deformation_gradient));
    rValues.SetDeformationGradientF(plastic_deformation_gradient);
    this->CalculateValue(rValues, GREEN_LAGRANGE_STRAIN_VECTOR, plastic_strain);
    rValues.SetDeterminantF(det_f);
    rValues.SetDeformationGradientF(deformation_gradient);

    // Check
    KRATOS_ERROR_IF(det_f < std::numeric_limits<double>::epsilon()) << "Deformation gradient determinant (detF) < 0.0 : " << det_f << std::endl;

    // Initialize Plastic Parameters
    double uniaxial_stress = 0.0, plastic_denominator = 0.0;
    BoundedArrayType yield_surface_derivative = ZeroVector(VoigtSize);     // DF/DS
    BoundedArrayType plastic_potential_derivative = ZeroVector(VoigtSize); // DG/DS
    const BoundedArrayType dummy_plastic_strain_increment = ZeroVector(VoigtSize);

    // Elastic Matrix
    this->CalculateValue(rValues, CONSTITUTIVE_MATRIX_PK2, r_constitutive_matrix);

    // Compute the plastic parameters
    const double plastic_indicator = TConstLawIntegratorType::CalculatePlasticParameters(
        predictive_stress_vector, r_strain_vector, uniaxial_stress,
        threshold, plastic_denominator, yield_surface_derivative, plastic_potential_derivative,
        plastic_dissipation, dummy_plastic_strain_increment,
        r_constitutive_matrix, rValues, characteristic_length,
        plastic_strain);

    if (plastic_indicator > std::abs(1.0e-4 * threshold)) { // Plastic case
        // While loop backward euler
        /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
        TConstLawIntegratorType::IntegrateStressVector(
            *this, GREEN_LAGRANGE_STRAIN_VECTOR, PK2_STRESS_VECTOR,
            predictive_stress_vector, r_strain_vector, uniaxial_stress, threshold,
            plastic_denominator, yield_surface_derivative, plastic_potential_derivative,
            plastic_dissipation, plastic_deformation_gradient, mPreviousDeformationGradient,
            r_constitutive_matrix, rValues, characteristic_length);
    }

    mPlasticDissipation = plastic_dissipation;
    mPlasticDeformationGradient = plastic_deformation_gradient;
    mPreviousDeformationGradient = deformation_gradient;
    mThreshold = threshold;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Auxiliar values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // The characteristic length
    const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());

    this->CalculateValue(rValues, ALMANSI_STRAIN_VECTOR, r_strain_vector);

    // We get some variables
    double threshold = this->GetThreshold();
    double plastic_dissipation = this->GetPlasticDissipation();
    Vector plastic_strain;
    Matrix plastic_deformation_gradient = this->GetPlasticDeformationGradient();

    // We backup the deformation gradient
    const double& r_det_deformation_gradient_backup = rValues.GetDeterminantF();
    const Matrix& r_deformation_gradient_backup = rValues.GetDeformationGradientF();

    // We compute the predicted stress vector
    BoundedArrayType predictive_stress_vector;
    if( r_constitutive_law_options.Is( ConstitutiveLaw::U_P_LAW ) ) {
        predictive_stress_vector = rValues.GetStressVector();
    } else {
        // We compute the elastic deformation gradient  Fe = plastic_indicator * inv(Fp)
        Matrix inverse_F_p ( Dimension, Dimension );
        double aux_det_Fp = 0;

        MathUtils<double>::InvertMatrix( plastic_deformation_gradient, inverse_F_p, aux_det_Fp);
        const Matrix elastic_deformation_gradient = prod(r_deformation_gradient_backup, inverse_F_p);

        rValues.SetDeterminantF(MathUtils<double>::DetMat(elastic_deformation_gradient));
        rValues.SetDeformationGradientF(elastic_deformation_gradient);
        Vector auxiliar_predictive_stress_vector;
        this->CalculateValue(rValues, KIRCHHOFF_STRESS_VECTOR, auxiliar_predictive_stress_vector);
        for (std::size_t i_voigt = 0; i_voigt < VoigtSize; ++i_voigt) {
            predictive_stress_vector[i_voigt] =  auxiliar_predictive_stress_vector[i_voigt];
        }

        // We revert the deformation gradient
        rValues.SetDeterminantF(r_det_deformation_gradient_backup);
        rValues.SetDeformationGradientF(r_deformation_gradient_backup);

        // Check
        KRATOS_ERROR_IF(r_det_deformation_gradient_backup < std::numeric_limits<double>::epsilon()) << "Deformation gradient determinant (detF) < 0.0 : " << r_det_deformation_gradient_backup << std::endl;
    }

    // We compute the plastic strain
    const double& det_f = rValues.GetDeterminantF();
    const Matrix& deformation_gradient = rValues.GetDeformationGradientF();
    rValues.SetDeterminantF(MathUtils<double>::DetMat(plastic_deformation_gradient));
    rValues.SetDeformationGradientF(plastic_deformation_gradient);
    this->CalculateValue(rValues, ALMANSI_STRAIN_VECTOR, plastic_strain);
    rValues.SetDeterminantF(det_f);
    rValues.SetDeformationGradientF(deformation_gradient);

    // Check
    KRATOS_ERROR_IF(det_f < std::numeric_limits<double>::epsilon()) << "Deformation gradient determinant (detF) < 0.0 : " << det_f << std::endl;

    // Initialize Plastic Parameters
    double uniaxial_stress = 0.0, plastic_denominator = 0.0;
    BoundedArrayType yield_surface_derivative = ZeroVector(VoigtSize);     // DF/DS
    BoundedArrayType plastic_potential_derivative = ZeroVector(VoigtSize); // DG/DS
    const BoundedArrayType dummy_plastic_strain_increment = ZeroVector(VoigtSize);

    // Elastic Matrix
    this->CalculateValue(rValues, CONSTITUTIVE_MATRIX_KIRCHHOFF, r_constitutive_matrix);

    // Compute the plastic parameters
    const double plastic_indicator = TConstLawIntegratorType::CalculatePlasticParameters(
        predictive_stress_vector, r_strain_vector, uniaxial_stress,
        threshold, plastic_denominator, yield_surface_derivative, plastic_potential_derivative,
        plastic_dissipation, dummy_plastic_strain_increment,
        r_constitutive_matrix, rValues, characteristic_length,
        plastic_strain);

    if (plastic_indicator > std::abs(1.0e-4 * threshold)) { // Plastic case
        // While loop backward euler
        /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
        TConstLawIntegratorType::IntegrateStressVector(
            *this, ALMANSI_STRAIN_VECTOR, KIRCHHOFF_STRESS_VECTOR,
            predictive_stress_vector, r_strain_vector, uniaxial_stress, threshold,
            plastic_denominator, yield_surface_derivative, plastic_potential_derivative,
            plastic_dissipation, plastic_deformation_gradient, mPreviousDeformationGradient,
            r_constitutive_matrix, rValues, characteristic_length);
    }

    mPlasticDissipation = plastic_dissipation;
    mPlasticDeformationGradient = plastic_deformation_gradient;
    mPreviousDeformationGradient = deformation_gradient;
    mThreshold = threshold;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    FinalizeMaterialResponseKirchhoff(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
bool GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
bool GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::Has(const Variable<Vector>& rThisVariable)
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
bool GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::Has(const Variable<Matrix>& rThisVariable)
{
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        mPlasticDissipation = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
void GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::SetValue(
    const Variable<Matrix>& rThisVariable,
    const Matrix& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == PLASTIC_DEFORMATION_GRADIENT) {
        mPlasticDeformationGradient = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
double& GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        rValue = mPlasticDissipation;
    } else {
        BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
Vector& GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        const Matrix C_plastic_tensor = prod(trans( mPlasticDeformationGradient), mPlasticDeformationGradient);
        ConstitutiveLawUtilities<VoigtSize>::CalculateHenckyStrain(C_plastic_tensor, rValue);
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
Matrix& GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == PLASTIC_DEFORMATION_GRADIENT) {
        rValue = mPlasticDeformationGradient;
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
double& GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == UNIAXIAL_STRESS) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        // Calculate the stress vector
        CalculateMaterialResponseCauchy(rParameterValues);
        const Vector& r_stress_vector = rParameterValues.GetStressVector();
        const Vector& r_strain_vector = rParameterValues.GetStrainVector();

        BoundedArrayType aux_stress_vector = r_stress_vector;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress( aux_stress_vector, r_strain_vector, rValue, rParameterValues);

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else if (rThisVariable == EQUIVALENT_PLASTIC_STRAIN) {
        // Calculate the stress vector
        Vector stress_vector;
        this->CalculateValue(rParameterValues, PK2_STRESS_VECTOR, stress_vector);

        // The plastic deformation gradient
        const Matrix& r_plastic_deformation_gradient = this->GetPlasticDeformationGradient();

        // We backup the deformation gradient
        const double& deformation_gradient_determinant_backup = rParameterValues.GetDeterminantF();
        const Matrix& r_deformation_gradient_backup = rParameterValues.GetDeformationGradientF();

        rParameterValues.SetDeterminantF(MathUtils<double>::DetMat(r_plastic_deformation_gradient));
        rParameterValues.SetDeformationGradientF(r_plastic_deformation_gradient);
        Vector plastic_strain;
        this->CalculateValue(rParameterValues, GREEN_LAGRANGE_STRAIN_VECTOR, plastic_strain);

        // Compute the equivalent plastic strain
        double uniaxial_stress;
        this->CalculateValue(rParameterValues, UNIAXIAL_STRESS, uniaxial_stress);
        TConstLawIntegratorType::CalculateEquivalentPlasticStrain(stress_vector, uniaxial_stress, plastic_strain, 0.0, rParameterValues, rValue);

        // We revert the deformation gradient
        rParameterValues.SetDeterminantF(deformation_gradient_determinant_backup);
        rParameterValues.SetDeformationGradientF(r_deformation_gradient_backup);

        return rValue;
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
Vector& GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
Matrix& GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == PLASTIC_STRAIN_TENSOR) {
        const Matrix C_plastic_tensor = prod(trans( mPlasticDeformationGradient), mPlasticDeformationGradient);
        Vector plastic_strain(VoigtSize);
        ConstitutiveLawUtilities<VoigtSize>::CalculateHenckyStrain(C_plastic_tensor, plastic_strain);
        rValue = MathUtils<double>::StrainVectorToTensor(plastic_strain);
    }else if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
        // The plastic deformation gradient
        const Matrix& r_plastic_deformation_gradient = this->GetPlasticDeformationGradient();

        // We backup the deformation gradient
        const Matrix& r_deformation_gradient_backup = rParameterValues.GetDeformationGradientF();

        // We compute the elastic deformation gradient  Fe = plastic_indicator * inv(Fp)
        Matrix inverse_F_p ( Dimension, Dimension );
        double aux_det_Fp = 0;
        MathUtils<double>::InvertMatrix( r_plastic_deformation_gradient, inverse_F_p, aux_det_Fp);
        const Matrix elastic_deformation_gradient = prod(r_deformation_gradient_backup, inverse_F_p);

        rParameterValues.SetDeformationGradientF(elastic_deformation_gradient);
        Vector elastic_stress;
        this->CalculateValue(rParameterValues, CAUCHY_STRESS_VECTOR, elastic_stress);

        rValue = MathUtils<double>::StressVectorToTensor(elastic_stress);

        // We revert the deformation gradient
        rParameterValues.SetDeformationGradientF(r_deformation_gradient_backup);

        return rValue;
    } else if (this->Has(rThisVariable)) {
        return this->GetValue(rThisVariable, rValue);
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TElasticBehaviourLaw, class TConstLawIntegratorType>
int GenericFiniteStrainIsotropicPlasticity<TElasticBehaviourLaw, TConstLawIntegratorType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    const int check_integrator = TConstLawIntegratorType::Check(rMaterialProperties);
    KRATOS_ERROR_IF_NOT(VoigtSize == this->GetStrainSize()) << "You are combining not compatible constitutive laws" << std::endl;
    if ((check_base + check_integrator) > 0) return 1;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

// Kirchhoff hyper elasticity
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>;

// Neo-Hookean hyper elasticity
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>;
} // namespace Kratos
