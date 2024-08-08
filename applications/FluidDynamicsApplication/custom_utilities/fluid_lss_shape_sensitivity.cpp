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
#include "custom_utilities/fluid_lss_sensitivity.h"
#include "custom_utilities/fluid_adjoint_slip_utilities.h"

// Include base h
#include "fluid_lss_shape_sensitivity.h"

namespace Kratos
{

FluidLSSShapeSensitivity::FluidLSSShapeSensitivity(
    Parameters Settings,
    const IndexType Dimension)
    : BaseType(),
      mDimension(Dimension)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
       "derivative_node_id"        : -1,
       "derivative_direction_index": -1
    })" );

    Settings.ValidateAndAssignDefaults(default_parameters);

    const int derivative_node_id = Settings["derivative_node_id"].GetInt();
    KRATOS_ERROR_IF(derivative_node_id <= 0) << "Derivative node id should be positive. [ derivative_node_id = " << derivative_node_id << " ].\n";
    mShapeDerivativeNodeId = derivative_node_id;

    const int derivative_direction_index = Settings["derivative_direction_index"].GetInt();
    KRATOS_ERROR_IF(derivative_direction_index < 0) << "Derivative direction should be positive. [ derivative_direction_index = " << derivative_direction_index << " ].\n";
    mShapeDerivativeDirection = derivative_direction_index;
    KRATOS_ERROR_IF(mShapeDerivativeDirection >= mDimension) << "Derivative direction should be less than the dimension. [ derivative_direction_index = " << mShapeDerivativeDirection << ", dimension = " << mDimension << " ].\n";

    const int number_of_threads = ParallelUtilities::GetNumThreads();
    mTLS.resize(number_of_threads);

    KRATOS_CATCH("");
}

const Variable<double>& FluidLSSShapeSensitivity::GetDerivativeVariable() const
{
    switch (mShapeDerivativeDirection) {
        case 0: return SHAPE_SENSITIVITY_X;
        case 1: return SHAPE_SENSITIVITY_Y;
        case 2: return SHAPE_SENSITIVITY_Z;
    }

    return Variable<double>::StaticObject();
}

void FluidLSSShapeSensitivity::CalculateResidualSensitivity(
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

void FluidLSSShapeSensitivity::CalculateResidualSensitivity(
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

double FluidLSSShapeSensitivity::CalculateResponseSensitivity(
    Element& rElement,
    AdjointResponseFunction& rResponse,
    const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
    const ProcessInfo& rProcessInfo)
{
    const int thread_id = OpenMPUtils::ThisThread();
    TLS& r_tls = mTLS[thread_id];
    return CalculateResponseSensitivity(rElement, r_tls.mResidualShapeDerivatives, r_tls.mResponseShapeDerivatives, rResponse, rProcessInfo);
}

double FluidLSSShapeSensitivity::CalculateResponseSensitivity(
    Condition& rCondition,
    AdjointResponseFunction& rResponse,
    const FluidAdjointSlipUtilities& rFluidAdjointSlipUtilities,
    const ProcessInfo& rProcessInfo)
{
    const int thread_id = OpenMPUtils::ThisThread();
    TLS& r_tls = mTLS[thread_id];
    return CalculateResponseSensitivity(rCondition, r_tls.mResidualShapeDerivatives, r_tls.mResponseShapeDerivatives, rResponse, rProcessInfo);
}


template<class TEntityType>
void FluidLSSShapeSensitivity::CalculateResidualSensitivity(
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
double FluidLSSShapeSensitivity::CalculateResponseSensitivity(
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

// template instantiations
template void FluidLSSShapeSensitivity::CalculateResidualSensitivity(Vector&, Condition&, Vector&, Matrix&, Matrix&, Matrix&, std::vector<IndexType>&, const FluidAdjointSlipUtilities&, const ProcessInfo&) const;
template void FluidLSSShapeSensitivity::CalculateResidualSensitivity(Vector&, Element&, Vector&, Matrix&, Matrix&, Matrix&, std::vector<IndexType>&, const FluidAdjointSlipUtilities&, const ProcessInfo&) const;

template double FluidLSSShapeSensitivity::CalculateResponseSensitivity(Condition&, Matrix&, Vector&, AdjointResponseFunction&, const ProcessInfo&) const;
template double FluidLSSShapeSensitivity::CalculateResponseSensitivity(Element&, Matrix&, Vector&, AdjointResponseFunction&, const ProcessInfo&) const;

} // namespace Kratos
