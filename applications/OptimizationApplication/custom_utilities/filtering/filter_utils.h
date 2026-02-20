//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>

// External includes
#include "nanoflann.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/model_part_utils.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TContainer>
class NanoFlannSingleContainerPositionAdapter {
public:
    ///@name Type definitions
    ///@{

    using ResultType = std::pair<unsigned int, double>;

    using ResultVectorType = std::vector<ResultType>;

    using PointerVectorType = std::vector<typename TContainer::value_type::Pointer>;

    ///@}
    ///@name Life cycle
    ///@{

    NanoFlannSingleContainerPositionAdapter(TContainer const * const pContainer)
        : mpContainer(pContainer)
    {}

    ///@}
    ///@name API member functions
    ///@{

    // Must return the number of data poins
    inline IndexType kdtree_get_point_count() const
    {
        return mpContainer->size();
    }

    // Must return the dim'th component of the idx'th point in the class:
    inline double kdtree_get_pt(
        const IndexType Index,
        const int Dimension) const
    {
        return OptimizationUtils::GetEntityPosition(*(mpContainer->begin() + Index))[Dimension];
    }

    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& bb) const
    {
        return false;
    }

    ///@}
    ///@name Public operations
    ///@{

    void GetResultingEntityPointersVector(
        PointerVectorType& rOutput,
        const ResultVectorType& rResultsVector) const
    {
        rOutput.resize(rResultsVector.size());
        std::transform(rResultsVector.begin(), rResultsVector.end(),
                        rOutput.begin(), [&](const auto& rResult) {
                            return *(mpContainer->ptr_begin() + rResult.first);
                        });
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    TContainer const * const mpContainer;

    ///@}
};

template<class TContainer>
class NanoFlannMultipleModelPartPositionAdapter {
public:
    ///@name Type definitions
    ///@{

    using ResultType = std::pair<unsigned int, double>;

    using ResultVectorType = std::vector<ResultType>;

    using PointerVectorType = std::vector<typename TContainer::value_type::Pointer>;

    ///@}
    ///@name Life cycle
    ///@{

    NanoFlannMultipleModelPartPositionAdapter(const std::vector<ModelPart*>& rModelParts)
    {
        IndexType total_size = 0;
        mModelPartOffsets.reserve(rModelParts.size());
        for (auto p_model_part : rModelParts) {
            mModelPartOffsets.push_back(total_size);
            const IndexType current_size = ModelPartUtils::GetContainer<TContainer>(*p_model_part).size();
            total_size += current_size;
        }

        mContainer.resize(total_size);

        IndexType offset = 0;
        for (auto p_model_part : rModelParts) {
            auto& r_container = ModelPartUtils::GetContainer<TContainer>(*p_model_part);
            std::copy(r_container.ptr_begin(), r_container.ptr_end(), mContainer.begin() + offset);
            offset += r_container.size();
        }
    }

    ///@}
    ///@name API member functions
    ///@{

    // Must return the number of data poins
    inline IndexType kdtree_get_point_count() const
    {
        return mContainer.size();
    }

    // Must return the dim'th component of the idx'th point in the class:
    inline double kdtree_get_pt(
        const IndexType Index,
        const int Dimension) const
    {
        return OptimizationUtils::GetEntityPosition(**(mContainer.begin() + Index))[Dimension];
    }

    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& bb) const
    {
        return false;
    }

    ///@}
    ///@name Public operations
    ///@{

    void GetResultingEntityPointersVector(
        PointerVectorType& rPointerVector,
        const ResultVectorType& rResultsVector) const
    {
        rPointerVector.resize(rResultsVector.size());
        std::transform(
            rResultsVector.begin(), rResultsVector.end(), rPointerVector.begin(),
            [&](const auto& rResult) { return mContainer[rResult.first]; });
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    PointerVectorType mContainer;

    std::vector<IndexType> mModelPartOffsets;

    ///@}
};

template<class TPointerVector>
struct NanoFlannKDTreeThreadLocalStorage
{
    // used to store the neighbour indices and squared distances
    typename std::vector<std::pair<unsigned int, double>> mNeighbourIndicesAndSquaredDistances;

    // used to store the list of neighbour entity pointers
    TPointerVector mNeighbourEntityPoints;

    // used to store the list of weights from filtering.
    std::vector<double> mListOfWeights;

    // used to store the list of weights in each component due to damping.
    std::vector<std::vector<double>> mListOfDampedWeights;
};

///@}

}