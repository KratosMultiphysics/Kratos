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
#include <tuple>

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/define.h"
#include "includes/global_pointer_variables.h"
#include "processes/find_global_nodal_neighbours_for_entities_process.h"
#include "utilities/compute_neighbour_list_functor.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/pointer_communicator.h"

// Application includes
#include "custom_utilities/rans_variable_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_smooth_clip_scalar_variable_process.h"

namespace Kratos
{
RansSmoothClipScalarVariableProcess::RansSmoothClipScalarVariableProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mVariableName = rParameters["variable_name"].GetString();
    mModelPartName = rParameters["model_part_name"].GetString();
    mEchoLevel = rParameters["echo_level"].GetInt();
    mMinValue = rParameters["min_value"].GetDouble();
    mMaxValue = rParameters["max_value"].GetDouble();
    mNumberOfSweeps = rParameters["max_number_of_sweeps"].GetInt();
    mAlwaysFindNeighbourNodes = rParameters["always_find_neighbour_nodes"].GetBool();
    mInverseDistanceWeightingPowerParameter = rParameters["inverse_distance_weighting_power_parameter"].GetDouble();
    mIsNeighbourNodesInitialized = false;

    this->UpdateExecutionPointsList(rParameters["execution_points"].GetStringArray());

    KRATOS_CATCH("");
}

int RansSmoothClipScalarVariableProcess::Check()
{
    KRATOS_TRY

    const auto& r_scalar_variable = KratosComponents<Variable<double>>::Get(mVariableName);
    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    KRATOS_ERROR_IF(!r_model_part.HasNodalSolutionStepVariable(r_scalar_variable))
        << mVariableName << " is not found in nodal solution step variables list of "
        << mModelPartName << ".";

    return 0;

    KRATOS_CATCH("");
}

void RansSmoothClipScalarVariableProcess::FindNeighbourNodes()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    FindNodalNeighboursForEntitiesProcess<ModelPart::ElementsContainerType> find_neighbour_process(r_model_part.GetCommunicator().GetDataCommunicator(), r_model_part, NEIGHBOUR_NODES);
    find_neighbour_process.ClearNeighbours();
    find_neighbour_process.Execute();

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Neighbour node search is completed for " << r_model_part.FullName() << ".\n";

    KRATOS_CATCH("");
}

void RansSmoothClipScalarVariableProcess::ExecuteOperation()
{
    KRATOS_TRY

    using node_type = ModelPart::NodeType;

    if (!mIsNeighbourNodesInitialized || mAlwaysFindNeighbourNodes) {
        mIsNeighbourNodesInitialized = true;
        FindNeighbourNodes();
    }

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    auto& r_communicator = r_model_part.GetCommunicator();
    auto& r_data_communicator = r_communicator.GetDataCommunicator();
    auto& r_nodes = r_communicator.LocalMesh().Nodes();
    const auto& r_scalar_variable = KratosComponents<Variable<double>>::Get(mVariableName);

    // compute number of neighbours within range and the maximum count
    // first create the global pointer communicator
    auto constructor_functor =
        Kratos::ComputeNeighbourListFunctor<ModelPart::NodesContainerType, Variable<GlobalPointersVector<node_type>>>(
            r_nodes, NEIGHBOUR_NODES);
    GlobalPointerCommunicator<node_type> pointer_comm(r_data_communicator, constructor_functor);

    struct NodalDataHolder {
    public:
        NodalDataHolder()
            : mIsBeingClipped(false),
              mScalarValue(0.0),
              mCoordinates(0.0)
        {
        }

        NodalDataHolder(const bool IsBeingClipper,
                        const double ScalarValue,
                        const array_1d<double, 3>& rCoordinates)
            : mIsBeingClipped(IsBeingClipper),
              mScalarValue(ScalarValue),
              mCoordinates(rCoordinates)
        {
        }

        bool mIsBeingClipped;
        double mScalarValue;
        array_1d<double, 3> mCoordinates;

        void save(Serializer& rSerializer) const
        {
            rSerializer.save("mIsBeingClipped", mIsBeingClipped);
            rSerializer.save("mScalarValue", mScalarValue);
            rSerializer.save("mCoordinates", mCoordinates);
        }

        void load(Serializer& rSerializer)
        {
            rSerializer.load("mIsBeingClipped", mIsBeingClipped);
            rSerializer.load("mScalarValue", mScalarValue);
            rSerializer.load("mCoordinates", mCoordinates);
        }
    };

    // now create the nodal value getter
    auto nodal_data_proxy =
        pointer_comm.Apply([&r_scalar_variable](GlobalPointer<node_type> const& rGP) -> NodalDataHolder {
            return NodalDataHolder(rGP->Is(SELECTED),
                                   rGP->FastGetSolutionStepValue(r_scalar_variable),
                                   rGP->Coordinates());
        });

    auto accumulate_nodal_data = [&](node_type& rNode) {
        const auto& r_nodal_neighbours = rNode.GetValue(NEIGHBOUR_NODES).GetContainer();
        int& numer_of_smoothing_neighbour_nodes = rNode.GetValue(NUMBER_OF_SMOOTHING_NEIGHBOUR_NODES);

        numer_of_smoothing_neighbour_nodes = 0;
        for (const auto& r_gp : r_nodal_neighbours) {
            const double scalar_value = nodal_data_proxy.Get(r_gp).mScalarValue;
            if (scalar_value >= mMinValue and scalar_value <= mMaxValue) {
                numer_of_smoothing_neighbour_nodes += 1;
            }
        }
        rNode.Set(SELECTED, true);
    };

    double global_maximum_smoothing_neighbour_nodes = 0;
    unsigned int number_of_nodes_below_minimum, number_of_nodes_above_maximum;
    for (int i_smoothing_itr = 0; i_smoothing_itr < mNumberOfSweeps; ++i_smoothing_itr) {
        // make this -1 to identify which nodes require smoothed clipping
        // if > 0 needs smoothed mapping; =0 needs smoothed mapping, but no
        // valid neighbours found ; else no smoothed mapping is required.
        block_for_each(r_nodes, [](node_type& rNode) {
            rNode.SetValue(NUMBER_OF_SMOOTHING_NEIGHBOUR_NODES, -1);
            rNode.Set(SELECTED, false);
        });

        std::tie(number_of_nodes_below_minimum, number_of_nodes_above_maximum) =
            block_for_each<CombinedReduction<SumReduction<unsigned int>, SumReduction<unsigned int>>>(
                r_nodes, [&](node_type& rNode) -> std::tuple<unsigned int, unsigned int> {
                    double& r_value = rNode.FastGetSolutionStepValue(r_scalar_variable);

                    if (r_value < mMinValue) {
                        accumulate_nodal_data(rNode);
                        return std::make_tuple<unsigned int, unsigned int>(1, 0);
                    } else if (r_value > mMaxValue) {
                        accumulate_nodal_data(rNode);
                        return std::make_tuple<unsigned int, unsigned int>(0, 1);
                    }

                    return std::make_tuple<unsigned int, unsigned int>(0, 0);
                });

        // Stores followings
        // index - 0 : number_of_nodes_below_minimum
        // index - 1 : number_of_nodes_above_maximum
        std::vector<unsigned int> nodes_count = {number_of_nodes_below_minimum,
                                                number_of_nodes_above_maximum};
        const std::vector<unsigned int>& total_nodes_count =
            r_communicator.GetDataCommunicator().SumAll(nodes_count);

        const double maximum_smoothing_neighbour_nodes =
            block_for_each<MaxReduction<int>>(r_nodes, [](const node_type& rNode) -> int {
                return rNode.GetValue(NUMBER_OF_SMOOTHING_NEIGHBOUR_NODES);
            });

        global_maximum_smoothing_neighbour_nodes =
            r_data_communicator.MaxAll(maximum_smoothing_neighbour_nodes);

        if (global_maximum_smoothing_neighbour_nodes == -1) {
            // all the nodes have their scalar values within clipping range.
            // therefore existing
            break;
        } else if (global_maximum_smoothing_neighbour_nodes > 0) {
            // if found remaining values out of clipping range, then use inverse distance mapping
            block_for_each(r_nodes, [&](node_type& rNode) {
                if (rNode.GetValue(NUMBER_OF_SMOOTHING_NEIGHBOUR_NODES) == global_maximum_smoothing_neighbour_nodes) {
                    double& nodal_scalar_value =  rNode.FastGetSolutionStepValue(r_scalar_variable);

                    nodal_scalar_value = 0.0;
                    const auto& r_node_coordinate = rNode.Coordinates();
                    double total_inverse_distance_weighting = 0.0;
                    const auto& r_nodal_neighbours = rNode.GetValue(NEIGHBOUR_NODES).GetContainer();
                    for (const auto& r_gp : r_nodal_neighbours) {
                        const auto& nodal_data = nodal_data_proxy.Get(r_gp);
                        if (!nodal_data.mIsBeingClipped) {
                            const double nodal_distance = norm_2(nodal_data.mCoordinates - r_node_coordinate);
                            const double inverse_distance_weighting = std::pow(1 / nodal_distance, mInverseDistanceWeightingPowerParameter);

                            nodal_scalar_value += inverse_distance_weighting * nodal_data.mScalarValue;
                            total_inverse_distance_weighting += inverse_distance_weighting;
                        }
                    }
                    nodal_scalar_value /= total_inverse_distance_weighting;
                }
            });

        } else {
            block_for_each(r_nodes, [&](node_type& rNode) {
                double& nodal_scalar_value = rNode.FastGetSolutionStepValue(r_scalar_variable);
                if (nodal_scalar_value < mMinValue) {
                    rNode.FastGetSolutionStepValue(r_scalar_variable) = mMinValue;
                } else if (nodal_scalar_value > mMaxValue) {
                    rNode.FastGetSolutionStepValue(r_scalar_variable) = mMaxValue;
                }
            });
        }

        // update nodal proxy data
        nodal_data_proxy.Update();

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0 && (total_nodes_count[0] > 0 ||
                                                        total_nodes_count[1] > 0))
            << mVariableName << " is clipped between [ " << mMinValue << ", "
            << mMaxValue << " ] with " << global_maximum_smoothing_neighbour_nodes
            << " maximum number of neighbour nodes. [ " << total_nodes_count[0]
            << " nodes < " << mMinValue << " and " << total_nodes_count[1]
            << " nodes > " << mMaxValue << " out of "
            << r_model_part.GetCommunicator().GlobalNumberOfNodes() << " total nodes in "
            << mModelPartName << ", smoothing iteration = " << (i_smoothing_itr + 1)
            << "/" << mNumberOfSweeps << " ].\n";
    }

    r_communicator.SynchronizeVariable(r_scalar_variable);

    KRATOS_CATCH("");
}

std::string RansSmoothClipScalarVariableProcess::Info() const
{
    return std::string("RansSmoothClipScalarVariableProcess");
}

void RansSmoothClipScalarVariableProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansSmoothClipScalarVariableProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters RansSmoothClipScalarVariableProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"                           : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"                             : "PLEASE_SPECIFY_SCALAR_VARIABLE",
            "echo_level"                                : 0,
            "min_value"                                 : 1e-18,
            "max_value"                                 : 1e+30,
            "max_number_of_sweeps"                      : 10000,
            "inverse_distance_weighting_power_parameter": 2,
            "always_find_neighbour_nodes"               : false,
            "execution_points"                          : ["after_coupling_solve_step"]
        })");

    return default_parameters;
}

} // namespace Kratos.
