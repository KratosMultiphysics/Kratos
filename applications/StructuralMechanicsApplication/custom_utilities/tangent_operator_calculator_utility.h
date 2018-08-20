// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:   Alejandro Cornejo & Lucia Barbu
//

#if !defined(KRATOS_TANGENT_OPERATOR_CALCULATOR_UTILITY_H_INCLUDED)
#define KRATOS_TANGENT_OPERATOR_CALCULATOR_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

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
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TangentOperatorCalculatorUtility
{
  public:
    /// Pointer definition of TangentOperatorCalculatorUtility
    KRATOS_CLASS_POINTER_DEFINITION(TangentOperatorCalculatorUtility);

    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// Constructor
    TangentOperatorCalculatorUtility()
    {
    }

    /// Destructor.
    virtual ~TangentOperatorCalculatorUtility() {}

    /**
     * @brief Main method that computes the tangent tensor
     * @param rValues The properties of the CL
     * @param pConstitutiveLaw Pointer to the CL
     */
    static void CalculateTangentTensor(
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveLaw *pConstitutiveLaw)
    {
        // Converged values to be storaged
        const Vector& strain_vector_gp = rValues.GetStrainVector();
        const Vector& stress_vector_gp = rValues.GetStressVector();

        Matrix& tangent_tensor = rValues.GetConstitutiveMatrix();
        tangent_tensor.clear();

        const std::size_t num_components = strain_vector_gp.size();
        // Loop over components of the strain
        for (std::size_t i_component = 0; i_component < num_components; ++i_component) {
            Vector &perturbed_strain = rValues.GetStrainVector();

            double pertubation;
            CalculatePerturbation(perturbed_strain, i_component, pertubation);
            PerturbateStrainVector(perturbed_strain, strain_vector_gp, pertubation, i_component);
            IntegratePerturbedStrain(rValues, pConstitutiveLaw);

            Vector& perturbed_integrated_stress = rValues.GetStressVector(); // now integrated
            const Vector& delta_stress = perturbed_integrated_stress - stress_vector_gp;
            AssignComponentsToTangentTensor(tangent_tensor, delta_stress, pertubation, i_component);

            // Reset the values to the initial ones
            noalias(perturbed_strain) = strain_vector_gp;
            noalias(perturbed_integrated_stress) = stress_vector_gp;
        }
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
        double& rPerturbation)
    {
        static constexpr double perturbation_coefficient_1 = 1.0e-5;
        static constexpr double perturbation_coefficient_2 = 1.0e-10;
        double perturbation_1, perturbation_2;
        if (std::abs(rStrainVector[Component]) > tolerance) {
            perturbation_1 = perturbation_coefficient_1 * rStrainVector[Component];
        } else {
            double min_strain_component;
            GetMinAbsValue(rStrainVector, min_strain_component);
            perturbation_1 = perturbation_coefficient_1 * min_strain_component;
        }
        double max_strain_component;
        GetMaxAbsValue(rStrainVector, max_strain_component);
        perturbation_2 = perturbation_coefficient_2 * max_strain_component;
        rPerturbation = std::max(perturbation_1, perturbation_2);
    }

    static void PerturbateStrainVector(
        Vector &rPerturbedStrainVector,
        const Vector& rStrainVectorGP,
        const double Perturbation,
        const int Component)
    {
        rPerturbedStrainVector = rStrainVectorGP;
        rPerturbedStrainVector[Component] += Perturbation;
    }

    /**
     * @brief This method integrates the pertubated strain
     * @param rValues The properties of the CL
     * @param pConstitutiveLaw Pointer to the CL
     */
    static void IntegratePerturbedStrain(
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveLaw *pConstitutiveLaw)
    {
        Flags& cl_options = rValues.GetOptions();
        // In order to avoid recursivity...
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        pConstitutiveLaw->CalculateMaterialResponseCauchy(rValues);
    }

    /**
     * @brief This method computes the maximum absolut value
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
     * @brief This method computes the minimim absolut value
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
     * @brief This assigns the values to the tangent tensor
     * @param rTangentTensor The desired tangent tensor
     * @param rDeltaStress The increment of stress
     * @param Perturbation The pertubation considered
     * @param Component Index of the component to compute
     */
    static void AssignComponentsToTangentTensor(
        Matrix& rTangentTensor,
        const Vector& rDeltaStress,
        const double Perturbation,
        const std::size_t Component
        )
    {
        const std::size_t  dimension = rDeltaStress.size();
        for (std::size_t row = 0; row < dimension; ++row) {
            rTangentTensor(row, Component) = rDeltaStress[row] / Perturbation;
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
