// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:   Alejandro Cornejo & Lucia Barbu
//  Collaborator:   Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class TangentOperatorCalculatorUtility
 * @ingroup StructuralMechanicsApplication
 * @brief An algorithm that derives numerically the constitutive tangent tensor at one GP
 * @details The procedure is defined in the PAPER "Caracterización de la delaminación en materiales
 compuestos mediante la teoría de mezclas serie/paralelo" X. Martinez, S. Oller y E. Barbero.
 * @authors Alejandro Cornejo & Lucia Barbu
 */
class TangentOperatorCalculatorUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TangentOperatorCalculatorUtility
    KRATOS_CLASS_POINTER_DEFINITION(TangentOperatorCalculatorUtility);

    /// Definition of size type
    typedef std::size_t SizeType;

    /// Definition of index type
    typedef std::size_t IndexType;

    /// Definition of the zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    // Definition of the perturbation coefficients
    static constexpr double PerturbationCoefficient1 = 1.0e-5;
    static constexpr double PerturbationCoefficient2 = 1.0e-10;

    // Definition of the perturbation threshold
    static constexpr double PerturbationThreshold = 1.0e-8;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    TangentOperatorCalculatorUtility()
    {
    }

    /// Destructor.
    virtual ~TangentOperatorCalculatorUtility() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Main method that computes the tangent tensor
     * @param rValues The properties of the CL
     * @param pConstitutiveLaw Pointer to the CL
     * @param rStressMeasure The stress measure of the law
     */
    static void CalculateTangentTensor(
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveLaw* pConstitutiveLaw,
        const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure_Cauchy,
        const bool ConsiderPertubationThreshold = true,
        const IndexType ApproximationOrder = 2
        )
    {
        // Ensure the proper flag
        Flags& cl_options = rValues.GetOptions();
        const bool use_element_provided_strain = cl_options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        if (use_element_provided_strain) {
            CalculateTangentTensorSmallDeformationProvidedStrain(rValues, pConstitutiveLaw, rStressMeasure, ConsiderPertubationThreshold, ApproximationOrder);
        } else {
            CalculateTangentTensorSmallDeformationNotProvidedStrain(rValues, pConstitutiveLaw, rStressMeasure, ConsiderPertubationThreshold, ApproximationOrder);
        }
    }

    /**
     * @brief Main method that computes the tangent tensor (provided strain)
     * @param rValues The properties of the CL
     * @param pConstitutiveLaw Pointer to the CL
     * @param rStressMeasure The stress measure of the law
     */
    static void CalculateTangentTensorSmallDeformationProvidedStrain(
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveLaw *pConstitutiveLaw,
        const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure_Cauchy,
        const bool ConsiderPertubationThreshold = true,
        const IndexType ApproximationOrder = 2
        )
    {
        // Converged values to be storaged
        const Vector unperturbed_strain_vector_gp = Vector(rValues.GetStrainVector());
        const Vector unperturbed_stress_vector_gp = Vector(rValues.GetStressVector());
        const auto &r_properties = rValues.GetMaterialProperties();
        const bool symmetrize_operator = (r_properties.Has(SYMMETRIZE_TANGENT_OPERATOR)) ? r_properties[SYMMETRIZE_TANGENT_OPERATOR] : false;

        // The number of components
        const SizeType num_components = unperturbed_stress_vector_gp.size();

        // The tangent tensor
        Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix();
        r_tangent_tensor.clear();
        Matrix auxiliary_tensor = ZeroMatrix(num_components,num_components);

        // Calculate the perturbation
        double pertubation = PerturbationThreshold;
        if (r_properties.Has(PERTURBATION_SIZE)) {
            pertubation = r_properties[PERTURBATION_SIZE];
            if (pertubation == -1.0)
                pertubation = std::sqrt(tolerance);
        } else {
            for (IndexType i_component = 0; i_component < num_components; ++i_component) {
                double component_perturbation;
                CalculatePerturbation(unperturbed_strain_vector_gp, i_component, component_perturbation);
                pertubation = std::max(component_perturbation, pertubation);
            }
            // We check that the perturbation has a threshold value of PerturbationThreshold
            if (ConsiderPertubationThreshold && pertubation < PerturbationThreshold) pertubation = PerturbationThreshold;
        }

        // Loop over components of the strain
        Vector& r_perturbed_strain = rValues.GetStrainVector();
        Vector& r_perturbed_integrated_stress = rValues.GetStressVector();
        if (ApproximationOrder == 1) {
            for (IndexType i_component = 0; i_component < num_components; ++i_component) {

                // Apply the perturbation
                PerturbateStrainVector(r_perturbed_strain, unperturbed_strain_vector_gp, pertubation, i_component);

                // We continue with the calculations
                IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

                // Compute tangent moduli
                const Vector delta_stress = r_perturbed_integrated_stress - unperturbed_stress_vector_gp;
                CalculateComponentsToTangentTensorFirstOrder(auxiliary_tensor, r_perturbed_strain-unperturbed_strain_vector_gp, delta_stress, i_component);

                // Reset the values to the initial ones
                noalias(r_perturbed_strain) = unperturbed_strain_vector_gp;
                noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;
            }
        } else if (ApproximationOrder == 2) {
            for (IndexType i_component = 0; i_component < num_components; ++i_component) {
                // Apply the perturbation (positive)
                PerturbateStrainVector(r_perturbed_strain, unperturbed_strain_vector_gp, pertubation, i_component);

                // We continue with the calculations
                IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

                // Compute stress (plus)
                const Vector strain_plus = r_perturbed_strain;
                const Vector stress_plus = r_perturbed_integrated_stress;

                // Reset the values to the initial ones
                noalias(r_perturbed_strain) = unperturbed_strain_vector_gp;
                noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;

                // Apply the perturbation (negative)
                PerturbateStrainVector(r_perturbed_strain, unperturbed_strain_vector_gp, - pertubation, i_component);

                // We continue with the calculations
                IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

                // Compute stress (minus)
                const Vector strain_minus = r_perturbed_strain;
                const Vector stress_minus = r_perturbed_integrated_stress;

                // Finally we compute the components
                CalculateComponentsToTangentTensorSecondOrder(auxiliary_tensor, strain_plus, strain_minus, stress_plus, stress_minus, i_component);

                // Reset the values to the initial ones
                noalias(r_perturbed_strain) = unperturbed_strain_vector_gp;
                noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;
            }
        } else if (ApproximationOrder == 4) { // Second order but with another approach. It is computed as first order but taking into account the second derivatives
            for (IndexType i_component = 0; i_component < num_components; ++i_component) {
                // Apply the perturbation (positive)
                PerturbateStrainVector(r_perturbed_strain, unperturbed_strain_vector_gp, pertubation, i_component);

                // We continue with the calculations
                IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

                // Compute stress (plus)
                const Vector strain_plus = r_perturbed_strain;
                const Vector stress_plus = r_perturbed_integrated_stress;

                // Reset the values to the initial ones
                noalias(r_perturbed_strain) = unperturbed_strain_vector_gp;
                noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;

                // Apply the perturbation twice
                PerturbateStrainVector(r_perturbed_strain, unperturbed_strain_vector_gp, 2.0*pertubation, i_component);

                // We continue with the calculations
                IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

                // Compute stress (minus)
                const Vector strain_2_plus = r_perturbed_strain;
                const Vector stress_2_plus = r_perturbed_integrated_stress;

                // Finally we compute the components
                const SizeType voigt_size = stress_plus.size();
                for (IndexType row = 0; row < voigt_size; ++row) {
                    auxiliary_tensor(row, i_component) = (stress_plus[row] - unperturbed_stress_vector_gp[row]) / pertubation - (stress_2_plus[row] - 2.0 * stress_plus[row] + unperturbed_stress_vector_gp[row]) / (2.0 * pertubation);
                }

                // Reset the values to the initial ones
                noalias(r_perturbed_strain) = unperturbed_strain_vector_gp;
                noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;
            }
        }
        if (symmetrize_operator)
            noalias(r_tangent_tensor) = 0.5*(auxiliary_tensor + trans(auxiliary_tensor));
        else
            noalias(r_tangent_tensor) = auxiliary_tensor;
    }

    /**
     * @brief Main method that computes the tangent tensor (not provided strain)
     * @param rValues The properties of the CL
     * @param pConstitutiveLaw Pointer to the CL
     * @param rStressMeasure The stress measure of the law
     */
    static void CalculateTangentTensorSmallDeformationNotProvidedStrain(
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveLaw *pConstitutiveLaw,
        const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure_Cauchy,
        const bool ConsiderPertubationThreshold = true,
        const IndexType ApproximationOrder = 2
        )
    {
        // Converged values to be storaged
        const Vector unperturbed_strain_vector_gp = Vector(rValues.GetStrainVector());
        const Vector unperturbed_stress_vector_gp = Vector(rValues.GetStressVector());

        const auto &r_properties = rValues.GetMaterialProperties();
        const bool symmetrize_operator = (r_properties.Has(SYMMETRIZE_TANGENT_OPERATOR)) ? r_properties[SYMMETRIZE_TANGENT_OPERATOR] : false;

        // The number of components
        const SizeType num_components = unperturbed_stress_vector_gp.size();

        // Converged values to be storaged (only used in case of elements that not provide the strain)
        const Matrix unperturbed_deformation_gradient_gp = Matrix(rValues.GetDeformationGradientF());
        const double det_unperturbed_deformation_gradient_gp = double(rValues.GetDeterminantF());

        // The tangent tensor
        Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix();
        r_tangent_tensor.clear();
        Matrix auxiliary_tensor = ZeroMatrix(num_components,num_components);

        const std::size_t size1 = unperturbed_deformation_gradient_gp.size1();
        const std::size_t size2 = unperturbed_deformation_gradient_gp.size2();

        KRATOS_ERROR_IF_NOT(ApproximationOrder == 1 || ApproximationOrder == 2) << "The approximation order for the perturbation is " << ApproximationOrder << ". Options are 1 and 2" << std::endl;

        // Calculate the perturbation
        double pertubation = PerturbationThreshold;
        if (r_properties.Has(PERTURBATION_SIZE)) {
            pertubation = r_properties[PERTURBATION_SIZE];
        } else {
            for (IndexType i_component = 0; i_component < num_components; ++i_component) {
                double component_perturbation;
                CalculatePerturbation(unperturbed_strain_vector_gp, i_component, component_perturbation);
                pertubation = std::max(component_perturbation, pertubation);
            }
            // We check that the perturbation has a threshold value of PerturbationThreshold
            if (ConsiderPertubationThreshold && pertubation < PerturbationThreshold) pertubation = PerturbationThreshold;
        }

        // Loop over components of the strain
        Vector& r_perturbed_strain = rValues.GetStrainVector();
        Vector& r_perturbed_integrated_stress = rValues.GetStressVector();
        Matrix& r_perturbed_deformation_gradient = const_cast<Matrix&>(rValues.GetDeformationGradientF());
        double& r_perturbed_det_deformation_gradient = const_cast<double&>(rValues.GetDeterminantF());

        if (ApproximationOrder == 1) {
            for (IndexType i_component = 0; i_component < size1; ++i_component) {
                for (IndexType j_component = i_component; j_component < size2; ++j_component) {
                    // Apply the perturbation
                    PerturbateDeformationGradient(r_perturbed_deformation_gradient, unperturbed_deformation_gradient_gp, pertubation, i_component, j_component);

                    // We continue with the calculations
                    IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

                    // Compute delta stress
                    const Vector delta_stress = r_perturbed_integrated_stress - unperturbed_stress_vector_gp;

                    // Finally we compute the components
                    const IndexType voigt_index = CalculateVoigtIndex(delta_stress.size(), i_component, j_component);
                    CalculateComponentsToTangentTensorFirstOrder(auxiliary_tensor, r_perturbed_strain-unperturbed_strain_vector_gp, delta_stress, voigt_index);

                    // Reset the values to the initial ones
                    noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;
                    noalias(r_perturbed_deformation_gradient) = unperturbed_deformation_gradient_gp;
                    r_perturbed_det_deformation_gradient = det_unperturbed_deformation_gradient_gp;
                }
            }
        } else if (ApproximationOrder == 2) {
            for (IndexType i_component = 0; i_component < size1; ++i_component) {
                for (IndexType j_component = i_component; j_component < size2; ++j_component) {
                    // Apply the perturbation (positive)
                    PerturbateDeformationGradient(r_perturbed_deformation_gradient, unperturbed_deformation_gradient_gp, pertubation, i_component, j_component);

                    // We continue with the calculations
                    IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

                    // Compute stress (plus)
                    const Vector strain_plus = r_perturbed_strain;
                    const Vector stress_plus = r_perturbed_integrated_stress;

                    // Reset the values to the initial ones
                    noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;
                    noalias(r_perturbed_deformation_gradient) = unperturbed_deformation_gradient_gp;
                    r_perturbed_det_deformation_gradient = det_unperturbed_deformation_gradient_gp;

                    // Apply the perturbation (negative)
                    PerturbateDeformationGradient(r_perturbed_deformation_gradient, unperturbed_deformation_gradient_gp, - pertubation, i_component, j_component);

                    // We continue with the calculations
                    IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

                    // Compute stress (minus)
                    const Vector strain_minus = r_perturbed_strain;
                    const Vector stress_minus = r_perturbed_integrated_stress;

                    // Finally we compute the components
                    const IndexType voigt_index = CalculateVoigtIndex(stress_plus.size(), i_component, j_component);
                    CalculateComponentsToTangentTensorSecondOrder(auxiliary_tensor, strain_plus, strain_minus, stress_plus, stress_minus, voigt_index);

                    // Reset the values to the initial ones
                    noalias(r_perturbed_integrated_stress) = unperturbed_stress_vector_gp;
                    noalias(r_perturbed_deformation_gradient) = unperturbed_deformation_gradient_gp;
                    r_perturbed_det_deformation_gradient = det_unperturbed_deformation_gradient_gp;
                }
            }
        }
        if (symmetrize_operator)
            noalias(r_tangent_tensor) = 0.5*(auxiliary_tensor + trans(auxiliary_tensor));
        else
            noalias(r_tangent_tensor) = auxiliary_tensor;
    }

    /**
     * @brief Main method that computes the tangent tensor (for finite deformation problems)
     * @param rValues The properties of the CL
     * @param pConstitutiveLaw Pointer to the CL
     * @param rStressMeasure The stress measure of the law
     */
    static void CalculateTangentTensorFiniteDeformation(
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveLaw* pConstitutiveLaw,
        const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure_PK2,
        const bool ConsiderPertubationThreshold = true,
        const IndexType ApproximationOrder = 2
        )
    {
        CalculateTangentTensorSmallDeformationNotProvidedStrain(rValues, pConstitutiveLaw, rStressMeasure, ConsiderPertubationThreshold, ApproximationOrder);
    }

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the pertubation
     * @param rStrainVector The vector of strains
     * @param Component Index of the component to compute
     * @param rPerturbation The resulting perturbation
     */
    static void CalculatePerturbation(
        const Vector& rStrainVector,
        const IndexType Component,
        double& rPerturbation
        )
    {
        double perturbation_1, perturbation_2;
        if (std::abs(rStrainVector[Component]) > tolerance) {
            perturbation_1 = PerturbationCoefficient1 * rStrainVector[Component];
        } else {
            double min_strain_component;
            GetMinAbsValue(rStrainVector, min_strain_component);
            perturbation_1 = PerturbationCoefficient1 * min_strain_component;
        }
        double max_strain_component;
        GetMaxAbsValue(rStrainVector, max_strain_component);
        perturbation_2 = PerturbationCoefficient2 * max_strain_component;
        rPerturbation = std::max(perturbation_1, perturbation_2);
    }

    /**
     * @brief This method computes the pertubation (for finite deformation problems)
     * @param rDeformationGradient The deformation gradient
     * @param ComponentI Index i of the component to compute
     * @param ComponentJ Index j of the component to compute
     * @param rPerturbation The resulting perturbation
     */
    static void CalculatePerturbationFiniteDeformation(
        const Matrix& rDeformationGradient,
        const IndexType ComponentI,
        const IndexType ComponentJ,
        double& rPerturbation
        )
    {
        double perturbation_1, perturbation_2;
        if (std::abs(rDeformationGradient(ComponentI, ComponentJ)) > tolerance) {
            perturbation_1 = PerturbationCoefficient1 * rDeformationGradient(ComponentI, ComponentJ);
        } else {
            double min_strain_component;
            GetMinAbsValue(rDeformationGradient, min_strain_component);
            perturbation_1 = PerturbationCoefficient1 * min_strain_component;
        }
        double max_strain_component;
        GetMaxAbsValue(rDeformationGradient, max_strain_component);
        perturbation_2 = PerturbationCoefficient2 * max_strain_component;
        rPerturbation = std::max(perturbation_1, perturbation_2);
    }

    /**
     * @brief This method perturbates the strain vector
     * @param rPerturbedStrainVector The strain vector to be perturbated
     * @param rStrainVectorGP It is the original strain vector
     * @param Perturbation The perturbation to be applied
     * @param Component Component of the vector to be perturbated
     */
    static void PerturbateStrainVector(
        Vector& rPerturbedStrainVector,
        const Vector& rStrainVectorGP,
        const double Perturbation,
        const IndexType Component
        )
    {
        noalias(rPerturbedStrainVector) = rStrainVectorGP;
        rPerturbedStrainVector[Component] += Perturbation;
    }

    /**
     * @brief This method perturbates the  deformation gradient
     * @param rPerturbedDeformationGradient The deformation gradient to be perturbated
     * @param rDeformationGradientGP It is the original deformation gradient
     * @param Perturbation The perturbation to be applied
     * @param ComponentI Component i of the matrix to be perturbated
     * @param ComponentJ Component j of the matrix to be perturbated
     */
    static void PerturbateDeformationGradient(
        Matrix& rPerturbedDeformationGradient,
        const Matrix& rDeformationGradientGP,
        const double Perturbation,
        const IndexType ComponentI,
        const IndexType ComponentJ
        )
    {
        Matrix aux_perturbation_matrix = IdentityMatrix(rDeformationGradientGP.size1());
        if (ComponentI == ComponentJ) {
            aux_perturbation_matrix(ComponentI, ComponentJ) += Perturbation;
        } else {
            aux_perturbation_matrix(ComponentI, ComponentJ) += 0.5 * Perturbation;
            aux_perturbation_matrix(ComponentJ, ComponentI) += 0.5 * Perturbation;
        }
        noalias(rPerturbedDeformationGradient) = prod(aux_perturbation_matrix, rDeformationGradientGP);
    }

    /**
     * @brief This method integrates the pertubated strain
     * @param rValues The properties of the CL
     * @param pConstitutiveLaw Pointer to the CL
     * @param rStressMeasure The stress measure of the law
     */
    static void IntegratePerturbedStrain(
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveLaw* pConstitutiveLaw,
        const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure_Cauchy
        )
    {
        Flags& cl_options = rValues.GetOptions();

        // In order to avoid recursivity...
        const bool flag_back_up_1 = cl_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_back_up_2 = cl_options.Is(ConstitutiveLaw::COMPUTE_STRESS);

        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        pConstitutiveLaw->CalculateMaterialResponse(rValues, rStressMeasure);

        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_back_up_1);
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_back_up_2);
    }

    /**
     * @brief This method computes the maximum absolut value (Vector)
     * @param rArrayValues The array containing the values
     * @param rMaxValue The value to be computed
     */
    static void GetMaxAbsValue(
        const Vector& rArrayValues,
        double& rMaxValue
        )
    {
        const SizeType dimension = rArrayValues.size();

        double aux = 0.0;
        for (IndexType i = 0; i < dimension; ++i) {
            if (std::abs(rArrayValues[i]) > aux) {
                aux = std::abs(rArrayValues[i]);
            }
        }

        rMaxValue = aux;
    }

    /**
     * @brief This method computes the maximum absolut value (Matrix)
     * @param rArrayValues The array containing the values
     * @param rMaxValue The value to be computed
     */
    static void GetMaxAbsValue(
        const Matrix& rMatrixValues,
        double& rMaxValue
        )
    {
        // Sizes of the matrix
        const SizeType size1 = rMatrixValues.size1();
        const SizeType size2 = rMatrixValues.size2();

        // The deformation gradient is an identity matrix for zero deformation
        const Matrix working_matrix = rMatrixValues - IdentityMatrix(size1);

        double aux = 0.0;
        for (IndexType i = 0; i < size1; ++i) {
            for (IndexType j = 0; j < size2; ++j) {
                if (std::abs(working_matrix(i, j)) > aux) {
                    aux = std::abs(working_matrix(i, j));
                }
            }
        }

        rMaxValue = aux;
    }

    /**
     * @brief This method computes the minimim absolut value (Vector)
     * @param rArrayValues The array containing the values
     * @param rMinValue The value to be computed
     */
    static void GetMinAbsValue(
        const Vector& rArrayValues,
        double& rMinValue
        )
    {
        const SizeType dimension = rArrayValues.size();

        double aux = std::numeric_limits<double>::max();
        for (IndexType i = 0; i < dimension; ++i) {
            if (std::abs(rArrayValues[i]) < aux) {
                aux = std::abs(rArrayValues[i]);
            }
        }

        rMinValue = aux;
    }

    /**
     * @brief This method computes the minimim absolut value (Matrix)
     * @param rArrayValues The array containing the values
     * @param rMinValue The value to be computed
     */
    static void GetMinAbsValue(
        const Matrix& rMatrixValues,
        double& rMinValue
        )
    {
       // Sizes of the matrix
        const SizeType size1 = rMatrixValues.size1();
        const SizeType size2 = rMatrixValues.size2();

        // The deformation gradient is an identity matrix for zero deformation
        const Matrix working_matrix = rMatrixValues - IdentityMatrix(size1);

        double aux = std::numeric_limits<double>::max();
        for (IndexType i = 0; i < size1; ++i) {
            for (IndexType j = 0; j < size2; ++j) {
                if (std::abs(working_matrix(i, j)) < aux) {
                    aux = std::abs(working_matrix(i, j));
                }
            }
        }

        rMinValue = aux;
    }

    /**
     * @brief This calculates the values to the tangent tensor
     * @param rTangentTensor The desired tangent tensor
     * @param rVectorStrain The perturbated strain considered
     * @param rDeltaStress The increment of stress
     * @param Component Index of the component to compute
     */
    static void CalculateComponentsToTangentTensorFirstOrder(
        Matrix& rTangentTensor,
        const Vector& rVectorStrain,
        const Vector& rDeltaStress,
        const IndexType Component
        )
    {
        const double perturbation = rVectorStrain[Component];
        const SizeType voigt_size = rDeltaStress.size();
        for (IndexType row = 0; row < voigt_size; ++row) {
            rTangentTensor(row, Component) = rDeltaStress[row] / perturbation;
        }
    }

    /**
     * @brief This calculates the values to the tangent tensor
     * @param rTangentTensor The desired tangent tensor
     * @param rVectorStrainPlus The positive perturbated strain considered
     * @param rVectorStrainMinus The negative perturbated strain considered
     * @param rStressPlus The stress with positive perturbation
     * @param rStressMinus The stress with negative perturbation
     * @param Component Index of the component to compute
     */
    static void CalculateComponentsToTangentTensorSecondOrder(
        Matrix& rTangentTensor,
        const Vector& rVectorStrainPlus,
        const Vector& rVectorStrainMinus,
        const Vector& rStressPlus,
        const Vector& rStressMinus,
        const IndexType Component
        )
    {
        const double perturbation = (rVectorStrainPlus[Component] - rVectorStrainMinus[Component]);
        const SizeType voigt_size = rStressPlus.size();
        for (IndexType row = 0; row < voigt_size; ++row) {
            rTangentTensor(row, Component) = (rStressPlus[row] - rStressMinus[row]) / perturbation;
        }
    }

    /**
     * @brief This calculates the proper Voigt index
     * @param rDeltaStress The increment of stress
     * @param Perturbation The pertubation considered
     * @param Component Index of the component to compute
     */
    static IndexType CalculateVoigtIndex(
        const SizeType VoigtSize,
        const IndexType ComponentI,
        const IndexType ComponentJ
        )
    {
        if (VoigtSize == 6) {
            switch(ComponentI) {
                case 0:
                    switch(ComponentJ) {
                        case 0:
                            return 0;
                        case 1:
                            return 3;
                        case 2:
                            return 5;
                        default:
                            return 0;
                    }
                case 1:
                    switch(ComponentJ) {
                        case 0:
                            return 3;
                        case 1:
                            return 1;
                        case 2:
                            return 4;
                        default:
                            return 0;
                    }
                case 2:
                    switch(ComponentJ) {
                        case 0:
                            return 5;
                        case 1:
                            return 4;
                        case 2:
                            return 2;
                        default:
                            return 0;
                    }
                default:
                    return 0;
            }
        } else {
            switch(ComponentI) {
                case 0:
                    switch(ComponentJ) {
                        case 0:
                            return 0;
                        case 1:
                            return 2;
                        default:
                            return 0;
                    }
                case 1:
                    switch(ComponentJ) {
                        case 0:
                            return 2;
                        case 1:
                            return 1;
                        default:
                            return 0;
                    }
                default:
                    return 0;
            }
        }
    }

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TangentOperatorCalculatorUtility &operator=(TangentOperatorCalculatorUtility const &rOther);
};
} // namespace Kratos.
