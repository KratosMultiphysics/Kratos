// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:   Alejandro Cornejo & Lucia Barbu
//  Collaborator:   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_TANGENT_OPERATOR_CALCULATOR_UTILITY_H_INCLUDED)
#define KRATOS_TANGENT_OPERATOR_CALCULATOR_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
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
        const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure_Cauchy
        )
    {
        // Ensure the proper flag
        Flags& cl_options = rValues.GetOptions();
        const bool use_element_provided_strain = cl_options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        if (use_element_provided_strain) {
            CalculateTangentTensorSmallDeformationProvidedStrain(rValues, pConstitutiveLaw, rStressMeasure);
        } else {
            CalculateTangentTensorSmallDeformationNotProvidedStrain(rValues, pConstitutiveLaw, rStressMeasure);
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
        const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure_Cauchy
        )
    {
        // Converged values to be storaged
        const Vector strain_vector_gp = rValues.GetStrainVector();
        const Vector stress_vector_gp = rValues.GetStressVector();

        // The number of components
        const SizeType num_components = stress_vector_gp.size();

        // The tangent tensor
        Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix();
        r_tangent_tensor.clear();
        Matrix auxiliar_tensor = ZeroMatrix(num_components,num_components);

        // Loop over components of the strain
        Vector& r_perturbed_strain = rValues.GetStrainVector();
        Vector& r_perturbed_integrated_stress = rValues.GetStressVector();
        for (IndexType i_component = 0; i_component < num_components; ++i_component) {
            // Calculate the perturbation
            double pertubation;
            CalculatePerturbation(r_perturbed_strain, i_component, pertubation);

            // We check that the perturbation has a threshold value of PerturbationThreshold
            //if (pertubation < PerturbationThreshold) pertubation = PerturbationThreshold;

            // Apply the perturbation
            PerturbateStrainVector(r_perturbed_strain, strain_vector_gp, pertubation, i_component);

            // We continue with the calculations
            IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

            const Vector& delta_stress = r_perturbed_integrated_stress - stress_vector_gp;
            CalculateComponentsToTangentTensor(auxiliar_tensor, delta_stress, pertubation, i_component);

            // Reset the values to the initial ones
            noalias(r_perturbed_strain) = strain_vector_gp;
            noalias(r_perturbed_integrated_stress) = stress_vector_gp;
        }
        noalias(r_tangent_tensor) = auxiliar_tensor;
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
        const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure_Cauchy
        )
    {
        // Converged values to be storaged
        const Vector stress_vector_gp = rValues.GetStressVector();

        // The number of components
        const SizeType num_components = stress_vector_gp.size();

        // Converged values to be storaged (only used in case of elements that not provide the strain)
        const Matrix unperturbed_deformation_gradient_gp = Matrix(rValues.GetDeformationGradientF());
        const double det_unperturbed_deformation_gradient_gp = double(rValues.GetDeterminantF());

        // The tangent tensor
        Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix();
        r_tangent_tensor.clear();
        Matrix auxiliar_tensor = ZeroMatrix(num_components,num_components);

        const std::size_t size1 = unperturbed_deformation_gradient_gp.size1();
        const std::size_t size2 = unperturbed_deformation_gradient_gp.size2();

        // Loop over components of the strain
        Vector& r_perturbed_integrated_stress = rValues.GetStressVector();
        Matrix& r_perturbed_deformation_gradient = const_cast<Matrix&>(rValues.GetDeformationGradientF());
        double& r_perturbed_det_deformation_gradient = const_cast<double&>(rValues.GetDeterminantF());
        for (IndexType i_component = 0; i_component < size1; ++i_component) {
            for (IndexType j_component = i_component; j_component < size2; ++j_component) {
                // Calculate the perturbation
                double pertubation;
                CalculatePerturbationFiniteDeformation(r_perturbed_deformation_gradient, i_component, j_component, pertubation);

                // We check that the perturbation has a threshold value of PerturbationThreshold
                //if (pertubation < PerturbationThreshold) pertubation = PerturbationThreshold;

                // Apply the perturbation
                PerturbateDeformationGradient(r_perturbed_deformation_gradient, unperturbed_deformation_gradient_gp, pertubation, i_component, j_component);

                // We continue with the calculations
                IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

                // Compute delta stress
                const Vector& delta_stress = r_perturbed_integrated_stress - stress_vector_gp;

                // Finally we compute the components
                const IndexType voigt_index = CalculateVoigtIndex(delta_stress, i_component, j_component);
                CalculateComponentsToTangentTensor(auxiliar_tensor, delta_stress, pertubation, voigt_index);

                // Reset the values to the initial ones
                noalias(r_perturbed_integrated_stress) = stress_vector_gp;
                noalias(r_perturbed_deformation_gradient) = unperturbed_deformation_gradient_gp;
                r_perturbed_det_deformation_gradient = det_unperturbed_deformation_gradient_gp;
            }
        }
        noalias(r_tangent_tensor) = auxiliar_tensor;
    }

    /**
     * @brief Main method that computes the tangent tensor (for finite deformation problems)
     * @param rValues The properties of the CL
     * @param pConstitutiveLaw Pointer to the CL
     * @param rStressMeasure The stress measure of the law
     */
    static void CalculateTangentTensorFiniteDeformation(
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveLaw *pConstitutiveLaw,
        const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure_PK2
        )
    {
        const Variable<Vector>& r_stress_variable = rStressMeasure == ConstitutiveLaw::StressMeasure_PK2 ? PK2_STRESS_VECTOR : KIRCHHOFF_STRESS_VECTOR;

        // Converged values to be storaged
        const Vector stress_vector_gp = rValues.GetStressVector();
        const Matrix unperturbed_deformation_gradient_gp = Matrix(rValues.GetDeformationGradientF());
        const double det_unperturbed_deformation_gradient_gp = double(rValues.GetDeterminantF());

        // The number of components
        const SizeType num_components = stress_vector_gp.size();

        const std::size_t size1 = unperturbed_deformation_gradient_gp.size1();
        const std::size_t size2 = unperturbed_deformation_gradient_gp.size2();

        double aux_det;
        Matrix inverse_perturbed_deformation_gradient(size1, size2);
        Matrix deformation_gradient_increment(size1, size2);

        // The tangent tensor
        Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix();
        r_tangent_tensor.clear();
        Matrix auxiliar_tensor = ZeroMatrix(num_components,num_components);

        // Loop over components of the strain
        Vector& r_perturbed_integrated_stress = rValues.GetStressVector();
        Matrix& r_perturbed_deformation_gradient = const_cast<Matrix&>(rValues.GetDeformationGradientF());
        double& r_perturbed_det_deformation_gradient = const_cast<double&>(rValues.GetDeterminantF());
        for (std::size_t i_component = 0; i_component < size1; ++i_component) {
            for (std::size_t j_component = i_component; j_component < size2; ++j_component) {
                // Calculate the perturbation
                double pertubation;
                CalculatePerturbationFiniteDeformation(r_perturbed_deformation_gradient, i_component, j_component, pertubation);

                // We check that the perturbation has a threshold value of PerturbationThreshold
                if (pertubation < PerturbationThreshold) pertubation = PerturbationThreshold;

                // Apply the perturbation
                PerturbateDeformationGradient(r_perturbed_deformation_gradient, unperturbed_deformation_gradient_gp, pertubation, i_component, j_component);

                // We continue with the calculations
                IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

                // We compute the new predictive stress vector // TODO: I need to check this, maybe I am overcomplicating everything
                MathUtils<double>::InvertMatrix(r_perturbed_deformation_gradient,inverse_perturbed_deformation_gradient, aux_det);
                deformation_gradient_increment = prod(inverse_perturbed_deformation_gradient, unperturbed_deformation_gradient_gp);
                rValues.SetDeterminantF(MathUtils<double>::DetMat(deformation_gradient_increment));
                rValues.SetDeformationGradientF(deformation_gradient_increment);
                Vector delta_stress;
                pConstitutiveLaw->CalculateValue(rValues, r_stress_variable, delta_stress);

                // Finally we compute the components
                const IndexType voigt_index = CalculateVoigtIndex(delta_stress, i_component, j_component);
                CalculateComponentsToTangentTensor(auxiliar_tensor, delta_stress, pertubation, voigt_index);

                // Reset the values to the initial ones
                noalias(r_perturbed_integrated_stress) = stress_vector_gp;
                noalias(r_perturbed_deformation_gradient) = unperturbed_deformation_gradient_gp;
                r_perturbed_det_deformation_gradient = det_unperturbed_deformation_gradient_gp;
            }
        }

        noalias(r_tangent_tensor) = auxiliar_tensor;
    }

    /**
     * @brief This method computes the pertubation
     * @param rStrainVector The vector of strains
     * @param Component Index of the component to compute
     * @param rPerturbation The resulting perturbation
     */
    static void CalculatePerturbation(
        const Vector& rStrainVector,
        const std::size_t Component,
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
        const std::size_t ComponentI,
        const std::size_t ComponentJ,
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
        const int Component
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
        const std::size_t ComponentI,
        const std::size_t ComponentJ
        )
    {
        noalias(rPerturbedDeformationGradient) = rDeformationGradientGP;
        rPerturbedDeformationGradient(ComponentI, ComponentJ) += Perturbation;
    }

    /**
     * @brief This method integrates the pertubated strain
     * @param rValues The properties of the CL
     * @param pConstitutiveLaw Pointer to the CL
     * @param rStressMeasure The stress measure of the law
     */
    static void IntegratePerturbedStrain(
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveLaw *pConstitutiveLaw,
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
        const std::size_t dimension = rArrayValues.size();
        std::vector<double> non_zero_values;

        for (std::size_t i = 1; i < dimension; ++i) {
            if (std::abs(rArrayValues[i]) > tolerance)
                non_zero_values.push_back(std::abs(rArrayValues[i]));
        }
        KRATOS_ERROR_IF(non_zero_values.size() == 0) << "The strain vector is full of 0's..." << std::endl;

        double aux = std::abs(non_zero_values[0]);
        for (std::size_t i = 1; i < non_zero_values.size(); ++i) {
            if (non_zero_values[i] > aux)
                aux = non_zero_values[i];
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
        const std::size_t size1 = rMatrixValues.size1();
        const std::size_t size2 = rMatrixValues.size2();
        std::vector<double> non_zero_values;

        for (std::size_t i = 0; i < size1; ++i) {
            for (std::size_t j = 0; j < size2; ++j) {
                if (std::abs(rMatrixValues(i, j)) > tolerance)
                    non_zero_values.push_back(std::abs(rMatrixValues(i, j)));
            }
        }
        KRATOS_ERROR_IF(non_zero_values.size() == 0) << "The deformation gradient is full of 0's..." << std::endl;

        double aux = std::abs(non_zero_values[0]);
        for (std::size_t i = 1; i < non_zero_values.size(); ++i) {
            if (non_zero_values[i] > aux)
                aux = non_zero_values[i];
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
        const std::size_t  dimension = rArrayValues.size();
        std::vector<double> non_zero_values;

        for (std::size_t i = 0; i < dimension; ++i) {
            if (std::abs(rArrayValues[i]) > tolerance)
                non_zero_values.push_back(std::abs(rArrayValues[i]));
        }
        KRATOS_ERROR_IF(non_zero_values.size() == 0) << "The strain vector is full of 0's..." << std::endl;

        double aux = std::abs(non_zero_values[0]);
        for (std::size_t i = 1; i < non_zero_values.size(); ++i) {
            if (non_zero_values[i] < aux)
                aux = non_zero_values[i];
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
        const std::size_t size1 = rMatrixValues.size1();
        const std::size_t size2 = rMatrixValues.size2();
        std::vector<double> non_zero_values;

        for (std::size_t i = 0; i < size1; ++i) {
            for (std::size_t j = 0; j < size2; ++j) {
                if (std::abs(rMatrixValues(i, j)) > tolerance)
                    non_zero_values.push_back(std::abs(rMatrixValues(i, j)));
            }
        }
        KRATOS_ERROR_IF(non_zero_values.size() == 0) << "The deformation gradient is full of 0's..." << std::endl;

        double aux = std::abs(non_zero_values[0]);
        for (std::size_t i = 1; i < non_zero_values.size(); ++i) {
            if (non_zero_values[i] < aux)
                aux = non_zero_values[i];
        }

        rMinValue = aux;
    }

    /**
     * @brief This calculates the values to the tangent tensor
     * @param rTangentTensor The desired tangent tensor
     * @param rDeltaStress The increment of stress
     * @param Perturbation The pertubation considered
     * @param Component Index of the component to compute
     */
    static void CalculateComponentsToTangentTensor(
        Matrix& rTangentTensor,
        const Vector& rDeltaStress,
        const double Perturbation,
        const IndexType Component
        )
    {
        const SizeType voigt_size = rDeltaStress.size();
        for (IndexType row = 0; row < voigt_size; ++row) {
            rTangentTensor(row, Component) = rDeltaStress[row] / Perturbation;
        }
    }

    /**
     * @brief This calculates the proper Voigt index
     * @param rDeltaStress The increment of stress
     * @param Perturbation The pertubation considered
     * @param Component Index of the component to compute
     */
    static IndexType CalculateVoigtIndex(
        const Vector& rDeltaStress,
        const IndexType ComponentI,
        const IndexType ComponentJ
        )
    {
        const SizeType voigt_size = rDeltaStress.size();
        if (voigt_size == 6) {
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
#endif // KRATOS_TANGENT_OPERATOR_CALCULATOR_PROCESS_H_INCLUDED  defined
