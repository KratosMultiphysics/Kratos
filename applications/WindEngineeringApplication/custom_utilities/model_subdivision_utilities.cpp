//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Core includes
#include "utilities/parallel_utilities.h"

// Project includes
#include "model_subdivision_utilities.h"


namespace Kratos
{
namespace Wind
{


std::vector<Kratos::shared_ptr<ModelPart>> ModelSubdivisionUtilities::SortNodesBySlabs(
    ModelPart& rModelPart,
    const array_1d<double,3>& rBottomPoint,
    const array_1d<double,3>& rTopPoint,
    std::size_t numberOfSlabs,
    bool isRootSlabOpen)
{
    KRATOS_TRY

    // Construct slabs
    const ModelSubdivisionUtilities::SlabStack slabs(rBottomPoint, rTopPoint, numberOfSlabs, isRootSlabOpen);

    // Construct sub model parts
    std::vector<ModelPart*> sub_model_parts;
    sub_model_parts.reserve(slabs.size());
    for (std::size_t i_slab=0; i_slab<slabs.size(); ++i_slab) {
        sub_model_parts.push_back(&rModelPart.CreateSubModelPart(
            rModelPart.Name() + "_slab_" + std::to_string(i_slab)
        ));
    }

    // Assign nodes to sub model parts
    auto sort_function = [&slabs, &sub_model_parts](const Node& rNode) -> ModelPart*
    {
        KRATOS_TRY
        return slabs.IsInside(rNode) ?
            sub_model_parts[slabs.SlabIndex(rNode)]
            :
            nullptr;
        KRATOS_CATCH("");
    };

    ModelSubdivisionUtilities::SortNodes(
        rModelPart,
        sort_function,
        sub_model_parts
    );

    // Convert raw pointers to smart ones
    std::vector<Kratos::shared_ptr<ModelPart>> output;
    output.reserve(sub_model_parts.size());

    for (std::size_t i_model_part=0; i_model_part<sub_model_parts.size(); ++i_model_part) {
        output.push_back(rModelPart.SubModelParts()(sub_model_parts[i_model_part]->Name()));
    }

    return output;

    KRATOS_CATCH("");
}


void ModelSubdivisionUtilities::SortNodes(
    ModelPart& rModelPart,
    std::function<ModelPart*(const Node&)> sortFunction,
    std::vector<ModelPart*>& rSubModelParts)
{
    KRATOS_TRY

    ModelSubdivisionUtilities::ThreadSafeIndexSet index_sets(rSubModelParts);

    block_for_each(rModelPart.Nodes(),
        [&](const Node& rNode)
        {
            index_sets.Push(
                sortFunction(rNode),
                rNode.Id()
            );
        }
    );

    void (ModelPart::*AddNodes)(const std::vector<ModelPart::IndexType>&, ModelPart::IndexType) = &ModelPart::AddNodes;
    index_sets.Apply(AddNodes);

    KRATOS_CATCH("");
}


void ModelSubdivisionUtilities::SortElements(
    ModelPart& rModelPart,
    std::function<ModelPart*(const Element&)> sortFunction,
    std::vector<ModelPart*>& rSubModelParts)
{
    KRATOS_TRY

    ModelSubdivisionUtilities::ThreadSafeIndexSet index_sets(rSubModelParts);

    block_for_each(rModelPart.Elements(),
        [&](const Element& rElement)
        {
            index_sets.Push(
                sortFunction(rElement),
                rElement.Id()
            );
        }
    );

    void (ModelPart::*AddElements)(const std::vector<ModelPart::IndexType>&, ModelPart::IndexType) = &ModelPart::AddElements;
    index_sets.Apply(AddElements);

    KRATOS_CATCH("");
}


void ModelSubdivisionUtilities::SortConditions(
    ModelPart& rModelPart,
    std::function<ModelPart*(const Condition&)> sortFunction,
    std::vector<ModelPart*>& rSubModelParts)
{
    KRATOS_TRY

    ModelSubdivisionUtilities::ThreadSafeIndexSet index_sets(rSubModelParts);

    block_for_each(rModelPart.Conditions(),
        [&](const Condition& rCondition)
        {
            index_sets.Push(
                sortFunction(rCondition),
                rCondition.Id()
            );
        }
    );

    void (ModelPart::*AddConditions)(const std::vector<ModelPart::IndexType>&, ModelPart::IndexType) = &ModelPart::AddConditions;
    index_sets.Apply(AddConditions);

    KRATOS_CATCH("");
}


ModelSubdivisionUtilities::ThreadSafeIndexSet::ThreadSafeIndexSet(std::vector<ModelPart*>& rSubModelParts)
{
    for (ModelPart* p_model_part : rSubModelParts) {
        mIndexSets.emplace(std::make_pair(
            p_model_part,
            std::make_pair(std::vector<IndexType>(), make_unique<LockType>())
        ));
    }
}


ModelSubdivisionUtilities::Slab::Slab(
    const array_1d<double,3>& rBottomPoint,
    const array_1d<double,3>& rTopPoint,
    bool isOpen) :
        mBottomPlane {rBottomPoint, rBottomPoint - rTopPoint},
        mTopPlane {rTopPoint, rTopPoint - rBottomPoint},
        mIsOpen(isOpen)
{
    if (MathUtils<double>::Norm(this->Normal()) < 1e-15) {
        KRATOS_ERROR << "Degenerate slab!";
    }
}


ModelSubdivisionUtilities::SlabStack::SlabStack(
    const array_1d<double,3> rBottomPoint,
    const array_1d<double,3> rTopPoint,
    std::size_t numberOfSlabs,
    bool isOpen) :
        ModelSubdivisionUtilities::Slab(
            rBottomPoint,
            rTopPoint,
            isOpen)
{
    KRATOS_TRY

    if (numberOfSlabs == 0) {
        KRATOS_ERROR << "Number of slabs in a stack must be at least 1";
    }
    else if (1 < numberOfSlabs) {
        mInnerPlanes.reserve(numberOfSlabs - 1);

        // Compute the vector by which the inner planes are offset from each other
        const double normal_length = MathUtils<double>::Norm(this->Normal());
        const double slab_height = MathUtils<double>::Dot(rTopPoint - rBottomPoint, this->Normal()) / normal_length / numberOfSlabs;
        const array_1d<double,3> normal_segment = (slab_height / normal_length) * this->Normal();

        // Generate inner planes (boundaries are not included)
        for (std::size_t i_slab=0; i_slab<numberOfSlabs-1; ++i_slab) {
            mInnerPlanes.emplace_back( ModelSubdivisionUtilities::Slab::Plane {
                rBottomPoint + (i_slab + 1) * normal_segment,
                this->Normal()
            } );
        }
    }

    KRATOS_CATCH("");
}


std::size_t ModelSubdivisionUtilities::SlabStack::SlabIndex(const array_1d<double,3>& rPoint) const
{
    KRATOS_TRY

    // Perform a binary search on the inner planes
    auto it_top_plane = std::upper_bound(
        mInnerPlanes.begin(),
        mInnerPlanes.end(),
        rPoint,
        [](const array_1d<double,3>& r_point, const Plane& r_plane)
        {
            return !r_plane.IsOnPositiveSide(r_point, false);
        }
    );

    std::size_t slab_index = std::numeric_limits<std::size_t>::max();

    // The inner planes do not contain the boundaries, so those
    // cases have to be checked separately
    if (it_top_plane == mInnerPlanes.begin()) {
        if (this->IsBelow(rPoint)) {
            KRATOS_ERROR << rPoint << " is below the slab stack!";
        }
        else {
            slab_index = 0;
        }
    }
    else if (it_top_plane == mInnerPlanes.end()) {
        if (this->IsAbove(rPoint)) {
            KRATOS_ERROR << rPoint << " is above the slab stack!";
        }
        else {
            slab_index = mInnerPlanes.size();
        }
    }
    else {
        slab_index = std::distance(mInnerPlanes.begin(), it_top_plane);
    }

    return slab_index;

    KRATOS_CATCH("");
}


std::size_t ModelSubdivisionUtilities::SlabStack::size() const
{
    return mInnerPlanes.size() + 1;
}


} // namespace Wind
} // namespace Kratos