//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_FLUID_LSS_SENSITIVITY_H)
#define KRATOS_FLUID_LSS_SENSITIVITY_H

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/kratos_parameters.h"
#include "utilities/parallel_utilities.h"
#include "utilities/openmp_utils.h"
#include "response_functions/adjoint_response_function.h"

// Application includes
#include "custom_utilities/fluid_adjoint_slip_utilities.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidLSSSensitivity
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(FluidLSSSensitivity);

    ///@}
    ///@name Life Cycle
    ///@{

    FluidLSSSensitivity()
    {
        KRATOS_TRY

        const int number_of_threads = ParallelUtilities::GetNumThreads();
        mTLS.resize(number_of_threads);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Operations
    ///@{

    void CalculateResidualSensitivity(
        Vector& rOutput,
        Element& rElement,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo)
    {
        const int thread_id = OpenMPUtils::ThisThread();
        TLS& r_tls = mTLS[thread_id];
        CalculateResidualSensitivity(rOutput, rElement, r_tls.mResiduals,
                r_tls.mAuxMatrix, r_tls.mResidualShapeDerivatives, r_tls.mRotatedResidualShapeDerivatives,
                r_tls.mDerivativeNodeIds, rFluidAdjointSlipUtilities, rProcessInfo);
    }

    void CalculateResidualSensitivity(
        Vector& rOutput,
        Condition& rCondition,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo)
    {
        const int thread_id = OpenMPUtils::ThisThread();
        TLS& r_tls = mTLS[thread_id];
        CalculateResidualSensitivity(rOutput, rCondition, r_tls.mResiduals,
                r_tls.mAuxMatrix, r_tls.mResidualShapeDerivatives, r_tls.mRotatedResidualShapeDerivatives,
                r_tls.mDerivativeNodeIds, rFluidAdjointSlipUtilities, rProcessInfo);
    }

    double CalculateResponseSensitivity(
        Element& rElement,
        AdjointResponseFunction& rResponse,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo)
    {
        const int thread_id = OpenMPUtils::ThisThread();
        TLS& r_tls = mTLS[thread_id];
        return CalculateResponseSensitivity(rElement, r_tls.mResidualShapeDerivatives, r_tls.mResponseShapeDerivatives, rResponse, rProcessInfo);
    }

    double CalculateResponseSensitivity(
        Condition& rCondition,
        AdjointResponseFunction& rResponse,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo)
    {
        const int thread_id = OpenMPUtils::ThisThread();
        TLS& r_tls = mTLS[thread_id];
        return CalculateResponseSensitivity(rCondition, r_tls.mResidualShapeDerivatives, r_tls.mResponseShapeDerivatives, rResponse, rProcessInfo);
    }

    ///@}
private:
    ///@name Private Classes
    ///@{

    struct TLS
    {
        Vector mResiduals;
        Matrix mAuxMatrix;
        Matrix mResidualShapeDerivatives;
        Matrix mRotatedResidualShapeDerivatives;
        Vector mResponseShapeDerivatives;
        std::vector<IndexType> mDerivativeNodeIds;
    };

    ///@}
    ///@name Private Members
    ///@{

    IndexType mShapeDerivativeNodeId;

    IndexType mShapeDerivativeDirection;

    IndexType mDimension;

    std::vector<TLS> mTLS;

    ///@}
    ///@name Private Static Operations
    ///@{

    template<class TEntityType>
    void CalculateResidualSensitivity(
        Vector& rOutput,
        TEntityType& rEntity,
        Vector& rResiduals,
        Matrix& rAuxMatrix,
        Matrix& rResidualShapeDerivatives,
        Matrix& rRotatedResidualShapeDerivatives,
        std::vector<IndexType>& rDerivativeNodeIds,
        const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
        const ProcessInfo& rProcessInfo) const
    {
        KRATOS_TRY

        // Calculates primal residual
        rEntity.CalculateLocalSystem(rAuxMatrix, rResiduals, rProcessInfo);
        rEntity.CalculateLocalVelocityContribution(rAuxMatrix, rResiduals, rProcessInfo);

        rEntity.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rResidualShapeDerivatives, rProcessInfo);
        rFluidAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedShapeVariableDerivatives(
            rRotatedResidualShapeDerivatives, rDerivativeNodeIds, rResiduals, rResidualShapeDerivatives, rEntity, rProcessInfo);

        if (rOutput.size() != rResiduals.size()) {
            rOutput.resize(rResiduals.size());
        }

        const auto& p_itr = std::find(rDerivativeNodeIds.begin(), rDerivativeNodeIds.end(), mShapeDerivativeNodeId);
        if (p_itr != rDerivativeNodeIds.end()) {
            noalias(rOutput) = row(rRotatedResidualShapeDerivatives, std::distance(rDerivativeNodeIds.begin(),  p_itr) * mDimension + mShapeDerivativeDirection);
        } else {
            noalias(rOutput) = ZeroVector(rResiduals.size());
        }

        KRATOS_CATCH("");
    }

    template<class TEntityType>
    double CalculateResponseSensitivity(
        TEntityType& rEntity,
        Matrix& rResidualShapeDerivatives,
        Vector& rResponseShapeDerivatives,
        AdjointResponseFunction& rResponseFunction,
        const ProcessInfo& rProcessInfo) const
    {
        KRATOS_TRY

        const auto& r_geometry = rEntity.GetGeometry();
        const auto& p_itr = std::find_if(r_geometry.begin(), r_geometry.end(), [&](const ModelPart::NodeType& rNode) -> bool {
            return rNode.Id() == mShapeDerivativeNodeId;
        });

        if (p_itr != r_geometry.end()) {
            rEntity.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rResidualShapeDerivatives, rProcessInfo);
            rResponseFunction.CalculatePartialSensitivity(rEntity, SHAPE_SENSITIVITY, rResidualShapeDerivatives, rResponseShapeDerivatives, rProcessInfo);
            return rResponseShapeDerivatives[std::distance(r_geometry.begin(), p_itr) * mDimension + mShapeDerivativeDirection];
        } else {
            return 0.0;
        }

        KRATOS_CATCH("");
    }


    ///@}
};

///@}

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_LSS_SENSITIVITY_H
