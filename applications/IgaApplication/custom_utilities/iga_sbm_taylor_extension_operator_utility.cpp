//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//

// System includes
#include <cmath>
#include <limits>

// External includes

// Project includes
#include "custom_utilities/iga_sbm_taylor_extension_operator_utility.h"
#include "custom_conditions/coupling_sbm_taylor_interface_6p_condition.h"
#include "iga_application_variables.h"
#include "includes/kratos_flags.h"

namespace Kratos
{

namespace {

using IndexType = IgaSbmTaylorExtensionOperatorUtility::IndexType;
using SizeType = IgaSbmTaylorExtensionOperatorUtility::SizeType;

// Shape-function-weighted position of rGeometry's own evaluation point.
array_1d<double, 3> PhysicalCenterOf(const Geometry<Node>& rGeometry)
{
    const Matrix& r_N = rGeometry.ShapeFunctionsValues();
    const SizeType n = rGeometry.PointsNumber();
    array_1d<double, 3> center = ZeroVector(3);
    for (IndexType i = 0; i < n; ++i) {
        center[0] += r_N(0, i) * rGeometry[i].X0();
        center[1] += r_N(0, i) * rGeometry[i].Y0();
        center[2] += r_N(0, i) * rGeometry[i].Z0();
    }
    return center;
}

// Single closest condition among rCandidates
Condition::Pointer FindClosestConditionByGeometry(
    const std::vector<Condition::Pointer>& rCandidates,
    const array_1d<double, 3>& rEvalPoint)
{
    Condition::Pointer best = nullptr;
    double best_dist = std::numeric_limits<double>::max();
    for (const auto& p_condition : rCandidates) {
        const double dist = norm_2(PhysicalCenterOf(p_condition->GetGeometry()) - rEvalPoint);
        if (dist < best_dist) {
            best_dist = dist;
            best = p_condition;
        }
    }
    KRATOS_ERROR_IF(best == nullptr)
        << "IgaSbmTaylorExtensionOperatorUtility: no candidate conditions to search." << std::endl;
    return best;
}

} 

void IgaSbmTaylorExtensionOperatorUtility::PrecomputeAndStoreDualCouplingTaylorData(
    const std::vector<Condition::Pointer>& rInterfaceConditions,
    const std::vector<Condition::Pointer>& rMasterSurrogateConditions,
    const std::vector<Condition::Pointer>& rSlaveSurrogateConditions)
{
    KRATOS_TRY

    for (const auto& p_interface_condition : rInterfaceConditions) {
        auto p_typed = dynamic_pointer_cast<CouplingSbmTaylorInterface6pCondition>(p_interface_condition);
        KRATOS_ERROR_IF(p_typed == nullptr)
            << "IgaSbmTaylorExtensionOperatorUtility: condition " << p_interface_condition->Id()
            << " is not a CouplingSbmTaylorInterface6pCondition." << std::endl;

        const array_1d<double, 3> eval_point_master =
            PhysicalCenterOf(p_interface_condition->GetGeometry().GetGeometryPart(0));
        const array_1d<double, 3> eval_point_slave =
            PhysicalCenterOf(p_interface_condition->GetGeometry().GetGeometryPart(1));

        Condition::Pointer p_master_source = FindClosestConditionByGeometry(rMasterSurrogateConditions, eval_point_master);
        Condition::Pointer p_slave_source = FindClosestConditionByGeometry(rSlaveSurrogateConditions, eval_point_slave);

        p_typed->SetShiftSources(p_master_source, p_slave_source);
    }

    KRATOS_CATCH("")
}

std::vector<Condition::Pointer> IgaSbmTaylorExtensionOperatorUtility::ReplaceWithInterfaceConditions(
    ModelPart& rCouplingModelPart,
    IndexType StartId)
{
    KRATOS_TRY

    std::vector<Condition::Pointer> old_conditions;
    old_conditions.reserve(rCouplingModelPart.NumberOfConditions());
    for (auto it = rCouplingModelPart.Conditions().begin(); it != rCouplingModelPart.Conditions().end(); ++it) {
        old_conditions.push_back(*(it.base()));
    }

    std::vector<Condition::Pointer> new_conditions;
    new_conditions.reserve(old_conditions.size());
    IndexType id = StartId;
    for (const auto& p_old_condition : old_conditions) {
        Condition::Pointer p_new_condition = rCouplingModelPart.CreateNewCondition(
            "CouplingSbmTaylorInterface6pCondition", id, p_old_condition->pGetGeometry(), p_old_condition->pGetProperties());
        new_conditions.push_back(p_new_condition);
        ++id;
    }

    for (const auto& p_old_condition : old_conditions) {
        p_old_condition->Set(TO_ERASE, true);
    }

    return new_conditions;

    KRATOS_CATCH("")
}

}  // namespace Kratos.
