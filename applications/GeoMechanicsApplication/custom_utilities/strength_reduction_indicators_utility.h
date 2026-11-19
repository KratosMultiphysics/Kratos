// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#pragma once

#include "includes/kratos_export_api.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

#include <cstddef>
#include <fstream>
#include <string>

namespace Kratos
{

/// @brief Records the evolution of physical quantities as a function of the strength reduction
/// factor (SRF) of a c-phi reduction analysis, such that failure can be detected afterwards using
/// different criteria. Each recorded increment is appended as one tab-separated row to a text file.
/// The utility only reads the model state; it does not change the c-phi algorithm.
class KRATOS_API(GEO_MECHANICS_APPLICATION) StrengthReductionIndicatorsUtility
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(StrengthReductionIndicatorsUtility);

    StrengthReductionIndicatorsUtility(ModelPart& rModelPart, Parameters Settings);

    /// @brief max_i ||u_i|| over all nodes (DISPLACEMENT)
    [[nodiscard]] double CalculateMaxDisplacement() const;

    /// @brief sum_i ||u_i|| over all nodes (DISPLACEMENT)
    [[nodiscard]] double CalculateSumOfDisplacementNorms() const;

    /// @brief sum_i A_i ||u_i|| / sum_i A_i, where A_i is the tributary area (2D) or volume (3D)
    /// of node i. A_i is assembled from the integrals of the shape functions over the continuum
    /// elements. For elements where one of those integrals is not positive (e.g. the corner nodes
    /// of 6-noded triangles and 10-noded tetrahedra), the element measure is distributed equally
    /// among its nodes instead.
    [[nodiscard]] double CalculateAreaWeightedMeanDisplacement() const;

    /// @brief The incremental variants use du_i = DISPLACEMENT(buffer 0) - DISPLACEMENT(buffer 1),
    /// i.e. the displacement increment with respect to the last converged step. This is also valid
    /// for a step that did not converge.
    [[nodiscard]] double CalculateMaxIncrementalDisplacement() const;
    [[nodiscard]] double CalculateSumOfIncrementalDisplacementNorms() const;
    [[nodiscard]] double CalculateAreaWeightedMeanIncrementalDisplacement() const;

    /// @brief W = 1/2 sum_elements sum_GaussPoints (sigma : eps_e) w_g |J|, with
    /// eps_e = C^-1 sigma and C the isotropic linear elastic matrix (YOUNG_MODULUS, POISSON_RATIO).
    /// sigma is the CAUCHY_STRESS_VECTOR of the elements, which includes the initial stresses.
    [[nodiscard]] double CalculateTotalElasticEnergy() const;

    /// @brief Calculates all indicators and appends them as one row to the output file. The number
    /// of non-linear iterations and the residual norm are taken from the ProcessInfo
    /// (NL_ITERATION_NUMBER and RESIDUAL_NORM).
    void RecordIncrement(std::size_t Step, double StrengthReductionFactor, bool IsConverged);

    [[nodiscard]] const std::string& GetOutputFileName() const;

private:
    ModelPart&    mrModelPart;
    std::string   mOutputFileName;
    std::ofstream mOutputFile;

    [[nodiscard]] bool IsWritingRank() const;
};

} // namespace Kratos
