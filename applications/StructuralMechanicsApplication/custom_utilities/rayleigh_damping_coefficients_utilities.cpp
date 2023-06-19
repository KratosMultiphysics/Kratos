// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/rayleigh_damping_coefficients_utilities.h"

namespace Kratos
{
namespace RayleighDampingCoefficientsUtilities
{
Vector ComputeDampingCoefficients(Parameters ThisParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
        "echo_level"          : 0,
        "damping_ratio_0"     : 0.0,
        "damping_ratio_1"     : -1.0,
        "eigen_values_vector" : [0.0]
    })" );

    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Declaring coefficients
    Vector coefficients_vector = ZeroVector(2);

    // Getting echo level
    const int echo_level = ThisParameters["echo_level"].GetInt();

    // The eigen values vector
    const Vector& r_eigenvalue_vector = ThisParameters["eigen_values_vector"].GetVector();

    // We need at least 2 eigenvalues
    if (r_eigenvalue_vector.size() > 1) {
        KRATOS_INFO_IF("ComputeRayleighDampingCoefficientsProcess", echo_level > 0) <<  "The two first eigen values of the structure are:\n" << r_eigenvalue_vector[0] << " Hz\t and " << r_eigenvalue_vector[1] << " Hz" << std::endl;

        // The main damping ration
        array_1d<double, 2> damping_ratios;
        damping_ratios[0] = ThisParameters["damping_ratio_0"].GetDouble();
        damping_ratios[1] = ThisParameters["damping_ratio_0"].GetDouble();

        // We define the auxiliary alpha and beta
        double alpha = 0.0;
        double beta = 0.0;

        // Depending of the second damping ratio (if defined or not) we compute
        if (ThisParameters["damping_ratio_1"].GetDouble() > 0.0) {
            damping_ratios[1] = ThisParameters["damping_ratio_1"].GetDouble();
            BoundedMatrix<double, 2, 2> frequencies_matrix;
            frequencies_matrix(0, 0) = 0.5 / r_eigenvalue_vector[0];
            frequencies_matrix(0, 1) = 0.5 * r_eigenvalue_vector[0];
            frequencies_matrix(1, 0) = 0.5 / r_eigenvalue_vector[1];
            frequencies_matrix(1, 1) = 0.5 * r_eigenvalue_vector[1];

            // We compute the inverse
            double det;
            BoundedMatrix<double, 2, 2> inverse_frequencies_matrix;
            MathUtils<double>::InvertMatrix(frequencies_matrix, inverse_frequencies_matrix, det);

            // We can compute now
            const array_1d<double, 2> aux_solution = prod(inverse_frequencies_matrix, damping_ratios);
            alpha = aux_solution[0];
            beta = aux_solution[1];
        } else {
            beta = (2.0 * damping_ratios[0])/(r_eigenvalue_vector[0] + r_eigenvalue_vector[1]);
            alpha = r_eigenvalue_vector[0] * r_eigenvalue_vector[1] * beta;
        }

        // Assigning to vector
        coefficients_vector[0] = alpha;
        coefficients_vector[1] = beta;

        KRATOS_INFO_IF("RayleighDampingCoefficientsUtilities", echo_level > 0) <<  "With the damping ratios of:\t" << damping_ratios[0] << "\t and " << damping_ratios[1] << "\n the RAYLEIGH_ALPHA is: " << alpha << "\tand RAYLEIGH_BETA: " << beta << std::endl;
    }

    return coefficients_vector;

    KRATOS_CATCH("")
}
} // namespace RayleighDampingCoefficientsUtilities
} // namespace Kratos
