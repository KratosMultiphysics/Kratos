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

    static void CalculateTangentTensor(
        ConstitutiveLaw::Parameters &rValues,
        ConstitutiveLaw *pConstitutiveLaw)
    {
        // Converged values to be storaged
        const Vector strain_vector_gp = rValues.GetStrainVector();
        const Vector stress_vector_gp = rValues.GetStressVector();

        Matrix &tangent_tensor = rValues.GetConstitutiveMatrix();
        tangent_tensor.clear();

        const std::size_t num_components = strain_vector_gp.size();
        // Loop over components of the strain
        for (std::size_t Component = 0; Component < num_components; ++Component)
        {
            Vector &perturbed_strain = rValues.GetStrainVector();

            double Perturbation;
            CalculatePerturbation(perturbed_strain, Component, Perturbation);
            PerturbateStrainVector(perturbed_strain, strain_vector_gp, Perturbation, Component);
            IntegratePerturbedStrain(rValues, pConstitutiveLaw);

            Vector &perturbed_integrated_stress = rValues.GetStressVector(); // now integrated
            const Vector &delta_stress = perturbed_integrated_stress - stress_vector_gp;
            AssignComponentsToTangentTensor(tangent_tensor, delta_stress, Perturbation, Component);

            // Reset the values to the initial ones
            noalias(perturbed_strain) = strain_vector_gp;
            noalias(perturbed_integrated_stress) = stress_vector_gp;
        }
    }

    static void CalculatePerturbation(
        const Vector &rStrainVector,
        const int Component,
        double &rPerturbation)
    {
        double Pert1, Pert2;
        if (std::abs(rStrainVector[Component]) > tolerance)
        {
            Pert1 = 1.0e-5 * rStrainVector[Component];
        }
        else
        {
            double MinStrainComp;
            GetMinAbsValue(rStrainVector, MinStrainComp);
            Pert1 = 1.0e-5 * MinStrainComp;
        }
        double MaxStrainComp;
        GetMaxAbsValue(rStrainVector, MaxStrainComp);
        Pert2 = 1.0e-10 * MaxStrainComp;
        rPerturbation = std::max(Pert1, Pert2);
    }

    static void PerturbateStrainVector(
        Vector &rPerturbedStrainVector,
        const Vector &rStrainVectorGP,
        const double Perturbation,
        const int Component)
    {
        rPerturbedStrainVector = rStrainVectorGP;
        rPerturbedStrainVector[Component] += Perturbation;
    }

    static void IntegratePerturbedStrain(
        ConstitutiveLaw::Parameters &rValues,
        ConstitutiveLaw *pConstitutiveLaw)
    {
        Flags &cl_options = rValues.GetOptions();
        // In order to avoid recursivity...
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        pConstitutiveLaw->CalculateMaterialResponseCauchy(rValues);
    }

    static void GetMaxAbsValue(
        const Vector &rArrayValues,
        double &rMaxValue)
    {
        const int Dim = rArrayValues.size();
        std::vector<double> non_zero_values;

        for (std::size_t i = 1; i < Dim; ++i)
        {
            if (std::abs(rArrayValues[i]) > tolerance)
                non_zero_values.push_back(std::abs(rArrayValues[i]));
        }
        KRATOS_ERROR_IF(non_zero_values.size() == 0) << "The strain vector is full of 0's..." << std::endl;

        double aux = std::abs(non_zero_values[0]);
        for (std::size_t i = 1; i < non_zero_values.size(); ++i)
        {
            if (non_zero_values[i] > aux)
                aux = non_zero_values[i];
        }

        rMaxValue = aux;
    }

    static void GetMinAbsValue(
        const Vector &ArrayValues,
        double &rMinValue)
    {
        const int Dim = ArrayValues.size();
        std::vector<double> non_zero_values;

        for (std::size_t i = 0; i < Dim; ++i)
        {
            if (std::abs(ArrayValues[i]) > tolerance)
                non_zero_values.push_back(std::abs(ArrayValues[i]));
        }
        KRATOS_ERROR_IF(non_zero_values.size() == 0) << "The strain vector is full of 0's..." << std::endl;

        double aux = std::abs(non_zero_values[0]);
        for (std::size_t i = 1; i < non_zero_values.size(); ++i)
        {
            if (non_zero_values[i] < aux)
                aux = non_zero_values[i];
        }

        rMinValue = aux;
    }

    static void AssignComponentsToTangentTensor(
        Matrix &TangentTensor,
        const Vector &DeltaStress,
        const double Perturbation,
        const int Component)
    {
        const int Dim = DeltaStress.size();
        for (std::size_t row = 0; row < Dim; ++row)
        {
            TangentTensor(row, Component) = DeltaStress[row] / Perturbation;
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