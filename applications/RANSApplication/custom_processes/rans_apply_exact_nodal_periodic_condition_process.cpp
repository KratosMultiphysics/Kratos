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
#include <limits>
#include <string>
#include <functional>

// External includes

// Project includes
#include "fluid_dynamics_application_variables.h"
#include "includes/define.h"
#include "processes/reorder_and_optimize_modelpart_process.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "rans_apply_exact_nodal_periodic_condition_process.h"

namespace Kratos
{
RansApplyExactNodalPeriodicConditionProcess::RansApplyExactNodalPeriodicConditionProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.RecursivelyValidateAndAssignDefaults(GetDefaultParameters());

    mMasterModelPartName = rParameters["master_model_part_name"].GetString();
    mSlaveModelPartName = rParameters["slave_model_part_name"].GetString();
    mReorder = rParameters["reorder"].GetBool();
    mTolerance = rParameters["tolerance"].GetDouble();
    mEchoLevel = rParameters["echo_level"].GetInt();

    const auto& r_master_model_part = rModel.GetModelPart(mMasterModelPartName);
    const auto& r_slave_model_part = rModel.GetModelPart(mSlaveModelPartName);

    KRATOS_ERROR_IF(&r_master_model_part.GetRootModelPart() !=
                    &r_slave_model_part.GetRootModelPart())
        << "Master model part [ " << mMasterModelPartName << " ] and Slave model part [ "
        << mSlaveModelPartName << " ] has different root model parts. [ "
        << r_master_model_part.GetRootModelPart().Name()
        << " != " << r_slave_model_part.GetRootModelPart().Name() << " ]\n";

    const double eps = std::numeric_limits<double>::epsilon();

    // translation settings
    noalias(mTranslationDirection) =
        rParameters["translation_settings"]["direction"].GetVector();
    const double translation_direction_norm = norm_2(mTranslationDirection);
    if (translation_direction_norm > eps) {
        noalias(mTranslationDirection) = mTranslationDirection / translation_direction_norm;
    }
    mTranslationMagnitude = rParameters["translation_settings"]["magnitude"].GetDouble();

    // rotation settings
    noalias(mRotationCenter) = rParameters["rotation_settings"]["center"].GetVector();
    noalias(mRotationAxis) = rParameters["rotation_settings"]["axis"].GetVector();
    const double rotation_axis_norm = norm_2(mRotationAxis);
    if (rotation_axis_norm > eps) {
        noalias(mRotationAxis) = mRotationAxis / rotation_axis_norm;
    }
    mRotationAngle = rParameters["rotation_settings"]["angle"].GetDouble();

    KRATOS_ERROR_IF(mRotationAngle < eps && mTranslationMagnitude < eps)
        << "Either translation or rotation should be defined.\n";
    KRATOS_ERROR_IF(mRotationAngle > eps && rotation_axis_norm < eps)
        << "Rotation angle defined, but rotation axis is not defined.\n";
    KRATOS_ERROR_IF(mTranslationMagnitude > eps && translation_direction_norm < eps) << "Translation magnitude defined, but translation direction is not defined.\n";

    const int domain_size =
        r_master_model_part.GetRootModelPart().GetProcessInfo()[DOMAIN_SIZE];
    if (domain_size == 2) {
        KRATOS_WARNING_IF(this->Info(), (mRotationAxis[0] > eps || mRotationAxis[1] > eps))
            << "Using 2D rotation axis as [0.0, 0.0, 1.0].\n";
        mRotationAxis[0] = 0.0;
        mRotationAxis[0] = 0.0;
        mRotationAxis[1] = 1.0;
    }

    if (mRotationAngle > eps && mTranslationMagnitude > eps) {
        noalias(mRotationCenter) += mTranslationDirection * mTranslationMagnitude;
    }

    KRATOS_CATCH("");
}

void RansApplyExactNodalPeriodicConditionProcess::ExecuteInitialize()
{
    CreatePeriodicConditions();
    if (mReorder) {
        auto& r_model_part = mrModel.GetModelPart(mMasterModelPartName).GetRootModelPart();
        Parameters default_params(R"({})");
        ReorderAndOptimizeModelPartProcess reorder_process(r_model_part, default_params);
        reorder_process.Execute();
    }
}

std::string RansApplyExactNodalPeriodicConditionProcess::Info() const
{
    return std::string("RansApplyExactNodalPeriodicConditionProcess");
}

void RansApplyExactNodalPeriodicConditionProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansApplyExactNodalPeriodicConditionProcess::PrintData(std::ostream& rOStream) const
{
}

void RansApplyExactNodalPeriodicConditionProcess::CreatePeriodicConditions()
{
    KRATOS_TRY

    auto& r_base_model_part = mrModel.GetModelPart(mMasterModelPartName).GetRootModelPart();
    int condition_id = r_base_model_part.NumberOfConditions();
    auto p_properties = r_base_model_part.CreateNewProperties(
        r_base_model_part.NumberOfProperties() + 1);

    auto& r_master_model_part_nodes = mrModel.GetModelPart(mMasterModelPartName).Nodes();
    auto& r_slave_model_part_nodes = mrModel.GetModelPart(mSlaveModelPartName).Nodes();

    KRATOS_ERROR_IF(r_master_model_part_nodes.size() != r_slave_model_part_nodes.size())
        << this->Info() << " failed. Master and slave model part nodes mismatch.\n";

    const int number_of_nodes = r_master_model_part_nodes.size();
    const double eps = std::numeric_limits<double>::epsilon();
    const auto& translation = mTranslationDirection * mTranslationMagnitude;

    const std::function<array_1d<double, 3>(const array_1d<double, 3>&)> translation_and_rotation_operator =
        [&](const array_1d<double, 3>& rInitialPosition) -> const array_1d<double, 3> {
        return CalculateRotatedPosition(rInitialPosition + translation);
    };

    const std::function<array_1d<double, 3>(const array_1d<double, 3>&)> rotation_operator =
        [&](const array_1d<double, 3>& rInitialPosition) -> const array_1d<double, 3> {
        return CalculateRotatedPosition(rInitialPosition);
    };

    const std::function<array_1d<double, 3>(const array_1d<double, 3>&)> translation_operator =
        [&](const array_1d<double, 3>& rInitialPosition) -> const array_1d<double, 3> {
        return rInitialPosition + translation;
    };

    const auto& r_operator = (mRotationAngle > eps && mTranslationMagnitude > eps)
                                 ? translation_and_rotation_operator
                                 : (mRotationAngle > eps) ? rotation_operator : translation_operator;

    block_for_each(r_master_model_part_nodes, [&](ModelPart::NodeType& rMasterNode) {
        const auto& master_initial_position =
            rMasterNode.GetInitialPosition().Coordinates();
        const auto& master_final_position = r_operator(master_initial_position);
        int& r_master_patch_index = rMasterNode.FastGetSolutionStepValue(PATCH_INDEX);
        for (int j_node = 0; j_node < number_of_nodes; ++j_node) {
            const auto& r_slave_node = *(r_slave_model_part_nodes.begin() + j_node);
            if (norm_2(master_final_position -
                       r_slave_node.GetInitialPosition().Coordinates()) < mTolerance) {
                r_master_patch_index = r_slave_node.Id();
                break;
            }
        }

        KRATOS_ERROR_IF(r_master_patch_index == 0)
            << "Slave node is not found in " << mSlaveModelPartName << " near "
            << master_final_position << " with tolerance " << mTolerance
            << " for master node id=" << rMasterNode.Id() << " at "
            << rMasterNode.GetInitialPosition().Coordinates() << " in "
            << mMasterModelPartName << ".\n";

        if (static_cast<std::size_t>(r_master_patch_index) != rMasterNode.Id()) {
            std::vector<std::size_t> node_id_list = {
                rMasterNode.Id(), static_cast<std::size_t>(r_master_patch_index)};
#pragma omp critical
            {
                condition_id++;
                auto p_condition = r_base_model_part.CreateNewCondition(
                    "PeriodicCondition", condition_id, node_id_list, p_properties);
                p_condition->Set(PERIODIC);
                rMasterNode.Set(PERIODIC);
            }
        }
    });

    block_for_each(r_base_model_part.Conditions(), [&](ModelPart::ConditionType& rCondition) {
        if (rCondition.Is(PERIODIC)) {
            const auto& r_node_master = rCondition.GetGeometry()[0];
            auto& r_node_slave = rCondition.GetGeometry()[1];

            r_node_slave.SetLock();
            r_node_slave.FastGetSolutionStepValue(PATCH_INDEX) = r_node_master.Id();
            r_node_slave.Set(PERIODIC);
            r_node_slave.UnSetLock();
        }
    });

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Created periodic conditions between " << mMasterModelPartName << " and "
        << mSlaveModelPartName << " in " << r_base_model_part.Name() << ".\n";

    KRATOS_CATCH("");
}

array_1d<double, 3> RansApplyExactNodalPeriodicConditionProcess::CalculateRotatedPosition(
    const array_1d<double, 3>& rInitialPosition) const
{
    BoundedMatrix<double, 3, 3> rotation_matrix;
    const auto local_coord = rInitialPosition - mRotationCenter;
    CalculateRotationMatrix(rotation_matrix);
    const auto global_coord = mRotationCenter + prod(rotation_matrix, local_coord);

    return global_coord;
}

void RansApplyExactNodalPeriodicConditionProcess::CalculateRotationMatrix(
    BoundedMatrix<double, 3, 3>& rOutput) const
{
    const double ux = mRotationAxis[0];
    const double uy = mRotationAxis[1];
    const double uz = mRotationAxis[2];
    const double c2 = std::cos(mRotationAngle);
    const double s2 = std::sin(mRotationAngle);
    const double s1 = std::sin(mRotationAngle / 2.);

    rOutput(0, 0) = c2 + 2. * ux * ux * s1 * s1;
    rOutput(1, 1) = c2 + 2. * uy * uy * s1 * s1;
    rOutput(2, 2) = c2 + 2. * uz * uz * s1 * s1;
    rOutput(0, 1) = 2. * ux * uy * s1 * s1 - uz * s2;
    rOutput(1, 0) = 2. * ux * uy * s1 * s1 + uz * s2;
    rOutput(0, 2) = 2. * ux * uz * s1 * s1 + uy * s2;
    rOutput(2, 0) = 2. * ux * uz * s1 * s1 - uy * s2;
    rOutput(1, 2) = 2. * uy * uz * s1 * s1 - ux * s2;
    rOutput(2, 1) = 2. * uy * uz * s1 * s1 + ux * s2;
}

const Parameters RansApplyExactNodalPeriodicConditionProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "master_model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "slave_model_part_name"          : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "tolerance"                      : 1e-9,
            "translation_settings" :
            {
                "direction" : [0.0, 0.0, 0.0],
                "magnitude" : 0.0
            },
            "rotation_settings" :
            {
                "axis"   : [0.0, 0.0, 0.0],
                "center" : [0.0, 0.0, 0.0],
                "angle"  : 0.0
            },
            "echo_level"                     : 0,
            "reorder"                        : true
        })");

    return default_parameters;
}

} // namespace Kratos.
