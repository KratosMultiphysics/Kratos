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
    static constexpr double PerturbationThreshold = 1.0e-6;

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
        ConstitutiveLaw *pConstitutiveLaw,
        const ConstitutiveLaw::StressMeasure& rStressMeasure = ConstitutiveLaw::StressMeasure_Cauchy
        )
    {
        // Converged values to be storaged
        const Vector unperturbed_strain_vector_gp = Vector(rValues.GetStrainVector());
        const Vector unperturbed_stress_vector_gp = Vector(rValues.GetStressVector());

        // Converged values to be storaged (only used in case of elements that not provide the strain)
        const Matrix& unperturbed_deformation_gradient_gp = rValues.GetDeformationGradientF();
        const double& det_unperturbed_deformation_gradient_gp = rValues.GetDeterminantF();
        Matrix equivalent_F;
        double det_equivalent_F;

        Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix();
        r_tangent_tensor.clear();

        // Ensure the proper flag
        Flags& cl_options = rValues.GetOptions();
        const bool use_element_provided_strain = cl_options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

        const SizeType num_components = unperturbed_strain_vector_gp.size();
        // Loop over components of the strain
        Vector& r_perturbed_strain = rValues.GetStrainVector();
        for (IndexType i_component = 0; i_component < num_components; ++i_component) {
            // Calculate the perturbation
            double pertubation;
            CalculatePerturbation(unperturbed_strain_vector_gp, i_component, pertubation);

            // We check that the perturbation has a threshold value of PerturbationThreshold
            if (pertubation < PerturbationThreshold) pertubation = PerturbationThreshold;

            // Apply the perturbation
            PerturbateStrainVector(r_perturbed_strain, unperturbed_strain_vector_gp, pertubation, i_component);

            // In case the element uses F as input instead of the strain vector
            if (!use_element_provided_strain) {
                noalias(equivalent_F) = ComputeEquivalentDeformationGradient(rValues);
                // Reset the values to the initial ones
                rValues.SetDeformationGradientF(equivalent_F);
                det_equivalent_F = MathUtils<double>::DetMat(equivalent_F);
                rValues.SetDeterminantF(det_equivalent_F);
            }

            // We continue with the calculations
            IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

            Vector& perturbed_integrated_stress = rValues.GetStressVector(); // now integrated
            const Vector delta_stress = perturbed_integrated_stress - unperturbed_stress_vector_gp;
            ComputeComponentsToTangentTensor(r_tangent_tensor, delta_stress, pertubation, i_component);

            // Reset the values to the initial ones
            noalias(r_perturbed_strain) = unperturbed_strain_vector_gp;
            noalias(perturbed_integrated_stress) = unperturbed_stress_vector_gp;

            // In case the element uses F as input instead of the strain vector
            if (!use_element_provided_strain){
                // Reset the values to the initial ones
                rValues.SetDeformationGradientF(unperturbed_deformation_gradient_gp);
                rValues.SetDeterminantF(det_unperturbed_deformation_gradient_gp);
            }
        }
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
        const Matrix& unperturbed_deformation_gradient_gp = rValues.GetDeformationGradientF();
        const double& det_unperturbed_deformation_gradient_gp = rValues.GetDeterminantF();

        // The size of the deformation gradient
        const SizeType size1 = unperturbed_deformation_gradient_gp.size1();
        const SizeType size2 = unperturbed_deformation_gradient_gp.size2();

        double aux_det;
        Matrix inverse_perturbed_deformation_gradient(size1, size2);
        Matrix perturbed_deformation_gradient(unperturbed_deformation_gradient_gp);
        Matrix deformation_gradient_increment(size1, size2);

        Matrix& tangent_tensor = rValues.GetConstitutiveMatrix();
        tangent_tensor.clear();

        // Loop over components of the strain
        for (IndexType i_component = 0; i_component < size1; ++i_component) {
            for (IndexType j_component = i_component; j_component < size2; ++j_component) { // Doing a symmetric perturbation

                // Calculate the perturbation
                double pertubation;
                CalculatePerturbationFiniteDeformation(unperturbed_deformation_gradient_gp, i_component, j_component, pertubation);

                // We check that the perturbation has a threshold value of PerturbationThreshold
                if (pertubation < PerturbationThreshold) pertubation = PerturbationThreshold;

                // Apply the perturbation
                PerturbateDeformationGradient(perturbed_deformation_gradient, unperturbed_deformation_gradient_gp, pertubation, i_component, j_component);

                // We continue with the calculations
                IntegratePerturbedStrain(rValues, pConstitutiveLaw, rStressMeasure);

                // We compute the new predictive stress vector
                MathUtils<double>::InvertMatrix(perturbed_deformation_gradient,inverse_perturbed_deformation_gradient, aux_det);
                deformation_gradient_increment = prod(inverse_perturbed_deformation_gradient, unperturbed_deformation_gradient_gp);
                rValues.SetDeterminantF(MathUtils<double>::DetMat(deformation_gradient_increment));
                rValues.SetDeformationGradientF(deformation_gradient_increment);
                Vector delta_stress;
                pConstitutiveLaw->CalculateValue(rValues, r_stress_variable, delta_stress);

                ComputeComponentsToTangentTensor(tangent_tensor, delta_stress, pertubation, i_component);

                // Reset the values to the initial ones
                rValues.SetDeformationGradientF(unperturbed_deformation_gradient_gp);
                rValues.SetDeterminantF(det_unperturbed_deformation_gradient_gp);
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
        const IndexType Component
        )
    {
        rPerturbedStrainVector = rStrainVectorGP;
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
        const SizeType ComponentI,
        const SizeType ComponentJ
        )
    {
        rPerturbedDeformationGradient = rDeformationGradientGP;
        rPerturbedDeformationGradient(ComponentI, ComponentJ) += Perturbation;
    }

    /**
     * @brief This method computes the equivalent deformation gradient for the elements which provide the deformation gradient as input
     * @param rValues The properties of the CL
     */
    static Matrix ComputeEquivalentDeformationGradient(ConstitutiveLaw::Parameters& rValues)
    {
        // We update the deformation gradient
        const Vector& r_strain_vector = rValues.GetStrainVector();
        const SizeType size = r_strain_vector.size();
        const SizeType F_size = (size == 6) ?  3 : 2;

        Matrix equivalent_F(F_size, F_size);

        for (IndexType i = 0; i < F_size; ++i) {
            equivalent_F(i, i) = 1.0 + r_strain_vector[i];
        }

        for (IndexType i = F_size; i < size; ++i) {
            const IndexType equivalent_i = (i == F_size) ? 0 : (i == 4) ? 1 : 0;
            const IndexType equivalent_j = (i == F_size) ? 1 : 2;
            equivalent_F(equivalent_i, equivalent_j) = 0.5 * r_strain_vector[i];
            equivalent_F(equivalent_j, equivalent_i) = 0.5 * r_strain_vector[i];
        }

        return equivalent_F;
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
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        pConstitutiveLaw->CalculateMaterialResponse(rValues, rStressMeasure);
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

        IndexType counter = 0;
        double aux = 0.0;
        for (IndexType i = 0; i < dimension; ++i) {
            if (std::abs(rArrayValues[i]) > aux) {
                aux = std::abs(rArrayValues[i]);
                ++counter;
            }
        }

        KRATOS_ERROR_IF(counter == 0) << "The strain vector is full of 0's..." << std::endl;

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

        IndexType counter = 0;
        double aux = 0.0;
        for (IndexType i = 0; i < size1; ++i) {
            for (IndexType j = 0; j < size2; ++j) {
                if (std::abs(working_matrix(i, j)) > aux) {
                    aux = std::abs(working_matrix(i, j));
                    ++counter;
                }
            }
        }

        KRATOS_ERROR_IF(counter == 0) << "The deformation gradient is an identity matrix" << std::endl;

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

        IndexType counter = 0;
        double aux = std::numeric_limits<double>::max();
        for (IndexType i = 0; i < dimension; ++i) {
            if (std::abs(rArrayValues[i]) < aux) {
                aux = std::abs(rArrayValues[i]);
                ++counter;
            }
        }

        KRATOS_ERROR_IF(counter == 0) << "The strain vector is full of 0's..." << std::endl;

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

        IndexType counter = 0;
        double aux = std::numeric_limits<double>::max();
        for (IndexType i = 0; i < size1; ++i) {
            for (IndexType j = 0; j < size2; ++j) {
                if (std::abs(working_matrix(i, j)) < aux) {
                    aux = std::abs(working_matrix(i, j));
                    ++counter;
                }
            }
        }

        KRATOS_ERROR_IF(counter == 0) << "The deformation gradient is an identity matrix" << std::endl;

        rMinValue = aux;
    }

    /**
     * @brief This assigns the values to the tangent tensor
     * @param rTangentTensor The desired tangent tensor
     * @param rDeltaStress The increment of stress
     * @param Perturbation The pertubation considered
     * @param Component Index of the component to compute
     */
    static void ComputeComponentsToTangentTensor(
        Matrix& rTangentTensor,
        const Vector& rDeltaStress,
        const double Perturbation,
        const IndexType Component
        )
    {
        const SizeType dimension = rDeltaStress.size();
        for (IndexType row = 0; row < dimension; ++row) {
            rTangentTensor(row, Component) = rDeltaStress[row] / Perturbation;
        }
    }

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
