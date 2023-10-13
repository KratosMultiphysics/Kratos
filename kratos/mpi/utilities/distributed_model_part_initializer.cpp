//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//


// System includes


// External includes


// Project includes
#include "utilities/string_utilities.h"
#include "parallel_fill_communicator.h"
#include "model_part_communicator_utilities.h"
#include "distributed_model_part_initializer.h"


namespace Kratos {

namespace {

static constexpr char delim = ';';

void GetSubModelPartHierarchy(const ModelPart& rModelPart, std::string& rModelPartHierarchy)
{
    for (const auto& r_smp : rModelPart.SubModelParts()) {
        if (rModelPartHierarchy.size() > 0) { // this is not the first time sth is added
            rModelPartHierarchy.append(std::string(1, delim));
        }
        rModelPartHierarchy.append(r_smp.FullName());
        GetSubModelPartHierarchy(r_smp, rModelPartHierarchy);
    }
}

void RecursiveCreateModelParts(ModelPart& rModelPart, const std::string& rModelPartName)
{
    auto mp_names = StringUtilities::SplitStringByDelimiter(rModelPartName, '.');
    auto model_part_name = mp_names[0];

    ModelPart& model_part = rModelPart.HasSubModelPart(model_part_name) ? rModelPart.GetSubModelPart(model_part_name) : rModelPart.CreateSubModelPart(model_part_name);

    if (mp_names.size() > 1) {
        RecursiveCreateModelParts(model_part, std::string(rModelPartName).erase(0, model_part_name.size()+1));
    }
}

void CreateSubModelPartHierarchy(ModelPart& rModelPart, const std::string& rModelPartHierarchy)
{
    for (auto& smp_name : StringUtilities::SplitStringByDelimiter(rModelPartHierarchy, delim)) {
        smp_name.erase(0, rModelPart.Name().size()+1); // remove main-model-part name
        RecursiveCreateModelParts(rModelPart, smp_name);
    }
}

} // anonymous namespace

void DistributedModelPartInitializer::CopySubModelPartStructure()
{
    std::string model_part_hierarchy;
    int size_model_part_hierarchy;

    if (mrDataComm.Rank() == mSourceRank) {
        GetSubModelPartHierarchy(mrModelPart, model_part_hierarchy);
        size_model_part_hierarchy = model_part_hierarchy.size();
    }

    // broadcast size to allocate memory
    mrDataComm.Broadcast(size_model_part_hierarchy, mSourceRank);

    if (mrDataComm.Rank() != mSourceRank) {
        model_part_hierarchy.resize(size_model_part_hierarchy);
    }

    // broadcast ModelPart hierarchy
    mrDataComm.Broadcast(model_part_hierarchy, mSourceRank);

    if (mrDataComm.Rank() != mSourceRank) {
        CreateSubModelPartHierarchy(mrModelPart, model_part_hierarchy);
    }
}

void DistributedModelPartInitializer::Execute()
{
    ModelPartCommunicatorUtilities::SetMPICommunicator(mrModelPart, mrDataComm);

    CopySubModelPartStructure();

    // Compute communicaton plan and fill communicator meshes correctly
    ParallelFillCommunicator(mrModelPart, mrDataComm).Execute();
}

}  // namespace Kratos.
