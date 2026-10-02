//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

#pragma once

// Project includes
#include "fluid_dynamics_application_variables.h"
#include "includes/condition.h"

namespace Kratos {

/**
 * @class SbmBoundaryIntegrationUtility
 * @brief Converts the compact modeler data into explicit 2D SBM quadrature data.
 * @details Both the classical shifted-boundary condition and the Gap-SBM
 *          condition use the same variational assembly. Their only difference
 *          is how the trace of the parent basis is evaluated:
 *
 *          - Classical SBM integrates on the surrogate face and obtains
 *            TraceShapeFunctions with the first-order Taylor shift.
 *          - Gap-SBM integrates on the reconstructed boundary and evaluates
 *            the extended parent basis directly at that point. In this case
 *            ShapeFunctions and TraceShapeFunctions are identical.
 *
 *          Keeping this conversion here prevents conditions from depending
 *          on magic matrix-column indices or duplicating geometric logic.
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) SbmBoundaryIntegrationUtility
{
public:
    using IndexType = Condition::IndexType;
    using SizeType = Condition::SizeType;
    using GeometryType = Condition::GeometryType;

    static constexpr SizeType NumberOfIntegrationPoints = 2;
    static constexpr SizeType ShiftedProjectionColumns = 3;
    static constexpr SizeType ReconstructedQuadratureColumns = 7;

    enum class IntegrationMode
    {
        ShiftedSurrogate,
        ReconstructedBoundary
    };

    struct IntegrationPointData
    {
        IntegrationMode Mode = IntegrationMode::ShiftedSurrogate;
        Point EvaluationPoint;
        array_1d<double, 2> Normal = ZeroVector(2);
        double Weight = 0.0;

        // N is used for the stress/pressure evaluated at the integration
        // point. TraceN represents the physical-boundary velocity trace.
        Vector ShapeFunctions;
        Vector TraceShapeFunctions;
        Matrix ShapeFunctionGradients;
    };

    /**
     * @brief Return the integration mode encoded by the projection matrix.
     * @note A 2x3 matrix stores projected points for classical SBM. A 2x7
     *       matrix stores reconstructed-boundary quadrature for Gap-SBM.
     */
    static IntegrationMode GetIntegrationMode(const Condition& rCondition);

    /**
     * @brief Validate the common SBM geometric data stored by the modeler.
     * @param RequireReconstructedBoundary Reject classical 2x3 data when true.
     */
    static void CheckConditionData(
        const Condition& rCondition,
        bool RequireReconstructedBoundary = false);

    /**
     * @brief Build all geometric and basis data required by an SBM condition.
     */
    static IntegrationPointData CalculateIntegrationPointData(
        const Condition& rCondition,
        IndexType IntegrationPointIndex);

private:
    static void CalculateParentShapeFunctions(
        const Condition& rCondition,
        const Point& rGlobalPoint,
        Vector& rN,
        Matrix& rDN_DX);
};

} // namespace Kratos
