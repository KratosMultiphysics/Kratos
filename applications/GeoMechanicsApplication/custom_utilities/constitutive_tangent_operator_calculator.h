//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Uxue Chasco
//
//
#pragma once

// System includes
#include <algorithm>
#include <cmath>
#include <limits>

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
 * @class ConstitutiveTangentOperatorCalculator
 * @brief An algorithm that derives numerically the constitutive tangent tensor at one GP
 */

class ConstitutiveTangentOperatorCalculator
{
public:
    ///@name Type Definitions
    ///@{

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

        ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Main method that computes the tangent tensor
     * @param rValues The properties of the Constitutive Law
     * @param pConstitutiveLaw Pointer to the Constitutive Law
     */
    static void CalculateTangentTensor(
        ConstitutiveLaw::Parameters& rParameters,
        const ConstitutiveLaw& rConstitutiveLaw)
    {
        Flags& r_options = rParameters.GetOptions();
        KRATOS_ERROR_IF_NOT(r_options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
            << "ConstitutiveTangentOperatorCalculator only supports element-provided strains."
            << std::endl;


        const Vector unperturbed_strain = rParameters.GetStrainVector();
        const Vector unperturbed_stress = rParameters.GetStressVector();
        const bool compute_stress = r_options.Is(ConstitutiveLaw::COMPUTE_STRESS);
        const bool compute_constitutive_tensor =
            r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        const SizeType num_components = unperturbed_strain.size();
        KRATOS_ERROR_IF(unperturbed_stress.size() != num_components)
            << "The strain and stress vectors must have the same size." << std::endl;

        Matrix& r_tangent_tensor = rParameters.GetConstitutiveMatrix();
        if (r_tangent_tensor.size1() != num_components ||
            r_tangent_tensor.size2() != num_components) {
            r_tangent_tensor.resize(num_components, num_components, false);
        }
        r_tangent_tensor.clear();

        const double perturbation = CalculatePerturbation(unperturbed_strain);

        auto& r_strain = rParameters.GetStrainVector();
        auto& r_stress = rParameters.GetStressVector();
        r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        for (IndexType component = 0; component < num_components; ++component) {
            noalias(r_strain) = unperturbed_strain;
            r_strain[component] += perturbation;
            noalias(r_stress) = unperturbed_stress;

            auto p_positive_perturbation_law = rConstitutiveLaw.Clone();
            p_positive_perturbation_law->CalculateMaterialResponseCauchy(rParameters);
            const Vector positive_stress = r_stress;

            noalias(r_strain) = unperturbed_strain;
            r_strain[component] -= perturbation;
            noalias(r_stress) = unperturbed_stress;

            auto p_negative_perturbation_law = rConstitutiveLaw.Clone();
            p_negative_perturbation_law->CalculateMaterialResponseCauchy(rParameters);
            const Vector negative_stress = r_stress;

            for (IndexType row = 0; row < num_components; ++row) {
                r_tangent_tensor(row, component) =
                    (positive_stress[row] - negative_stress[row]) / (2.0 * perturbation);
                // KRATOS_ERROR_IF_NOT(std::isfinite(r_tangent_tensor(row, component)))
                    // << "Non-finite numerical tangent at (" << row << ", " << component << ")\n"
                    // << "strain: " << unperturbed_strain << '\n'
                    // << "positive stress: " << positive_stress << '\n'
                    // << "negative stress: " << negative_stress << '\n'
                    // << "perturbation: " << perturbation << std::endl;
            }
        }

        // const double tangent_norm = norm_frobenius(r_tangent_tensor);
        // KRATOS_WARNING_IF("ConstitutiveTangentOperatorCalculator",
        //                   tangent_norm <= tolerance)
        //     << "Numerical tangent is zero or nearly zero.\n";

        // KRATOS_INFO("ConstitutiveTangentOperatorCalculator")
        //     << "Numerical tangent calculated.\n"
        //     << "norm: " << tangent_norm << '\n'
        //     << "strain: " << unperturbed_strain << '\n'
        //     << "stress: " << unperturbed_stress << '\n'
        //     << "perturbation: " << perturbation << '\n'
        //     << "tangent:\n"
        //     << r_tangent_tensor << std::endl;

        noalias(r_strain) = unperturbed_strain;
        noalias(r_stress) = unperturbed_stress;
        r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, compute_stress);
        r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR,
                      compute_constitutive_tensor);
    }

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
     * @param rPerturbation The resulting perturbation
     */
    static double CalculatePerturbation(const Vector& rStrainVector)
    {
        double maximum_absolute_strain = 0.0;
        for (const auto strain_component : rStrainVector) {
            maximum_absolute_strain =
                std::max(maximum_absolute_strain, std::abs(strain_component));
        }

        if (maximum_absolute_strain <= tolerance) {
            return PerturbationThreshold;
        }

        const double perturbation_1 =
            PerturbationCoefficient1 * maximum_absolute_strain;
        const double perturbation_2 =
            PerturbationCoefficient2 * maximum_absolute_strain;

        return std::max(
            PerturbationThreshold, std::max(perturbation_1, perturbation_2));
    }
};

} // namespace Kratos
