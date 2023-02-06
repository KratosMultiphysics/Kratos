//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <string>
#include <sstream>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes

// Include base h
#include "extract_boundary_nodes_process.h"

namespace Kratos {
const Parameters ExtractBoundaryNodesProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"             : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "boundary_sub_model_part_name": "Auto_Boundary",
            "echo_level"                  : 0
        })");

    return default_parameters;
}

ExtractBoundaryNodesProcess::ExtractBoundaryNodesProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    rParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    mModelPartName = rParameters["model_part_name"].GetString();
    mBoundaryNodesSubModelPartName = rParameters["boundary_sub_model_part_name"].GetString();
    mEchoLevel = rParameters["echo_level"].GetInt();
}

int ExtractBoundaryNodesProcess::Check()
{
    return 0;
}

void ExtractBoundaryNodesProcess::ExecuteInitialize()
{
    KRATOS_TRY

    auto& r_main_model_part = mrModel.GetModelPart(mModelPartName);
    auto& r_boundary_model_part = r_main_model_part.CreateSubModelPart(mBoundaryNodesSubModelPartName);

    using IndexType = std::size_t;
    using hashmap = std::unordered_map<std::vector<IndexType>, IndexType, KeyHasherRange<std::vector<IndexType>>, KeyComparorRange<std::vector<IndexType>>>;

    class HashMapReduction
    {
    public:
        typedef std::vector<std::vector<IndexType>> value_type;
        typedef hashmap return_type;

        hashmap mValue = hashmap(); // deliberately making the member value public, to allow one to change it as needed

        /// access to reduced value
        return_type GetValue() const
        {
            return mValue;
        }

        /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
        void LocalReduce(const value_type value){
            for (const auto& r_boundary: value) {
                mValue[r_boundary] += 1;
            }
        }

        /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
        void ThreadSafeReduce(const HashMapReduction& rOther)
        {
            const std::lock_guard<LockObject> scope_lock(ParallelUtilities::GetGlobalLock());
            for (const auto& r_value : rOther.GetValue()) {
                mValue[r_value.first] += r_value.second;
            }
        }
    };

    const IndexType domain_size = r_main_model_part.GetProcessInfo()[DOMAIN_SIZE];

    const hashmap& r_boundary_count_map = block_for_each<HashMapReduction>(r_main_model_part.Elements(), [&](const auto& rElement) {
        const auto& r_geometry = rElement.GetGeometry();

        KRATOS_ERROR_IF(domain_size != r_geometry.LocalSpaceDimension()) << "ExtractBoundaryNodesProcess: This function does only work"
            <<" for solid elements in 3D and surface elements in 2D!" << std::endl;

        const auto& boundaries = r_geometry.GenerateBoundariesEntities();

        std::vector<std::vector<IndexType>> boundary_nodes(boundaries.size());

        for(IndexType i = 0; i < boundaries.size(); ++i) {
            const auto& boundary = boundaries[i];
            std::vector<IndexType> boundary_node_ids(boundary.size());

            // Store node ids
            for(IndexType j = 0; j < boundary.size(); ++j) {
                boundary_node_ids[j] = boundary[j].Id();
            }

            //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
            std::sort(boundary_node_ids.begin(), boundary_node_ids.end());
            boundary_nodes.push_back(boundary_node_ids);
        }

        return boundary_nodes;
    });

    // now generate the boundary node edges / faces for the local rank
    std::vector<std::vector<IndexType>> boundaries;
    for (const auto& r_item : r_boundary_count_map) {
        if (r_item.second == 1) {
            boundaries.push_back(r_item.first);
        }
    }

    // now check the boundaries whether they are on the ghost nodes
    const auto& global_boundaries = block_for_each<AccumReduction<std::vector<IndexType>>>(boundaries, [&](const auto& rValue) {
        IndexType number_of_nodes_in_ghost_mesh = 0;
        for (const IndexType boundary_node_id : rValue) {
            for (const auto& r_node : r_main_model_part.GetCommunicator().GhostMesh().Nodes()) {
                if (r_node.Id() == boundary_node_id) {
                    number_of_nodes_in_ghost_mesh += 1;
                    break;
                }
            }
        }

        if (number_of_nodes_in_ghost_mesh == rValue.size()) {
            return std::vector<IndexType>();
        } else {
            return rValue;
        }
    });

    std::vector<IndexType> nodes_list;
    for (const auto& r_boundary : global_boundaries) {
        for (const auto& r_id: r_boundary) {
            nodes_list.push_back(r_id);
        }
    }

    r_boundary_model_part.AddNodes(nodes_list);

    KRATOS_CATCH("");
}

std::string ExtractBoundaryNodesProcess::Info() const
{
    std::stringstream msg;
    msg << "ExtractBoundaryNodesProcess [ ModelPartName = " << mModelPartName << ", BoundaryModelPartName = " << mBoundaryNodesSubModelPartName << " ]";
    return msg.str();
}

void ExtractBoundaryNodesProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

} // namespace Kratos.
