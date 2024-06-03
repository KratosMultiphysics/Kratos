// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef FLOOD_FILL_INCLUDE_H
#define FLOOD_FILL_INCLUDE_H

//// STL includes
#include <array>
#include <algorithm>
#include <stack>
#include <set>

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/utilities/mapping_utilities.h"

namespace queso {

///@name QuESo Classes
///@{

class BRepOperator;

/**
 * @class  FloodFill
 * @author Manuel Messmer
 * @brief Provides methods to robustly classify elements / cells as interior, exterior or trimmed.
*/
class FloodFill {
public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<IntersectionStatusType> StatusVectorType;
    typedef std::stack<IndexType> IndexStackType;
    typedef std::vector<bool> BoolVectorType;
    // GroupSetType holds the following information <Partition index, Element Indices, IsInsideCount>
    typedef std::tuple<IndexType, std::set<IndexType>, int > GroupSetType;
    typedef std::vector<GroupSetType> GroupSetVectorType;
    typedef std::pair<Vector3i, Vector3i> PartitionBoxType;
    typedef std::vector<PartitionBoxType> PartitionBoxVectorType;
    typedef std::pair<int, int> Partition1DBoxType;
    typedef std::vector<std::vector<std::set<IndexType>>> BoundaryIndicesVectorType;

    ///@}
    ///@name Life cycle
    ///@{

    /// @brief Constructor of FloodFill.
    /// @param pBrepOperator
    /// @param Parameters
    FloodFill(const BRepOperator* pBrepOperator, const Parameters& Parameters) :
        mpBrepOperator(pBrepOperator), mMapper(Parameters), mLowerBound(Parameters.LowerBoundXYZ()),
        mUpperBound(Parameters.UpperBoundXYZ()), mNumberOfElements( Parameters.NumberOfElements() )
    {
        // Obtain discretization of background mesh.
        mDelta[0] = std::abs(mUpperBound[0] - mLowerBound[0]) / (mNumberOfElements[0]);
        mDelta[1] = std::abs(mUpperBound[1] - mLowerBound[1]) / (mNumberOfElements[1]);
        mDelta[2] = std::abs(mUpperBound[2] - mLowerBound[2]) / (mNumberOfElements[2]);
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns a ptr to a vector that holds the states of each element. Vector is ordered according to index -> see: Mapper.
    /// @brief This function runs a flood fill repeatively and classifies each group based on the bounding elements that are trimmed. Each element that borders a trimmed
    ///        element is tested via local ray tracing and marked as inside or outside. The majority vote decides about the classification of each group.
    /// @return Unique<StatusVectorType>.
    Unique<StatusVectorType> ClassifyElements() const;

protected:

    ///@}
    ///@name Protected Operations
    ///@{

    /// @brief Only used for testing. Also provides the actual groups.
    /// @see ClassifyElements()
    /// @return std::pair<Unique<StatusVectorType>, Unique<GroupSetVectorType>>
    std::pair<Unique<StatusVectorType>, Unique<GroupSetVectorType>> ClassifyElementsForTest() const;

private:

    ///@}
    ///@name Private Operations
    ///@{

    /// @brief Interface for ClassifyElements() and ClassifyElementsForTest().
    /// @param rGroupsOutput.
    /// @return Unique<StatusVectorType>.
    Unique<StatusVectorType> ClassifyElements(GroupSetVectorType& rGroupsOutput) const;

    /// @brief Partitions the background grid in "n_threads" stripes along the directon with "n_element_max" and filles each domain (see: SinglePartitionFill()).
    ///        Afterwards the domains are merged (see: MergeGroups).
    /// @param rGroupSetVector Vector to group sets: see GroupSetVectorType.
    /// @param rPartitions Vector of partitions.
    /// @param rStates Global classification vector (enum: Inside / outside or trimmed).
    void PartitionedFill(GroupSetVectorType& rGroupSetVector, PartitionBoxVectorType& rPartitions, StatusVectorType& rStates) const;

    /// @brief Starts flood fill from 'Index' over given partition and adds found elements to 'rGroupSet'.
    ///        Marks trimmed elements in global vector: rStates.
    ///        Marks visited elements in rVisited.
    /// @param Index of start element.
    /// @param rGroupSet Current Group.
    /// @param rPartition Current Partition.
    /// @param rStates Global classification vector.
    /// @param rVisited Vector<bool> for all elements.
    void Fill(IndexType Index, GroupSetType& rGroupSet, const PartitionBoxType &rPartition,
        StatusVectorType& rStates, BoolVectorType& rVisited ) const;

    /// @brief Move from current 'Index' towards 'Direction' to next index.
    ///        If next is trimmed, it is marked in rStates. And Is_inside cout of 'rGroupSet' is increased/decreased.
    ///        If next index is not already visited and not trimmed, it is returned. Otherwise -1 is returned.
    /// @param Index Current index.
    /// @param Direction Move Direction: 0:+x, 1:-x, 2:+y, 3:-y, 4:+z, 5:-z
    /// @param rGroupSet Current group.
    /// @param rPartition Current partition.
    /// @param rStates Global classification vector.
    /// @param rVisited Vector<bool> for all elements.
    /// @return NextIndex.
    int Move(IndexType Index, IndexType Direction, GroupSetType& rGroupSet,
        const PartitionBoxType& rPartition, StatusVectorType& rStates, BoolVectorType& rVisited ) const;

    /// @brief Merge groups that emerged from PartitionedFill.
    /// @param rGroups Contains all groups from all partitiones.
    /// @param [out] rMergedGroups Output vector of groups.
    /// @param PartitionDir Direction along the partition was performed. Direction with "n_element_max".
    /// @param rPartitions Vector with all partitions.
    /// @param rStates Global classification vector.
    void MergeGroups(GroupSetVectorType& rGroups, GroupSetVectorType& rMergedGroups, IndexType PartitionDir,
        PartitionBoxVectorType& rPartitions, StatusVectorType& rStates) const;


    /// @brief Run flood fill to merge groups starting at 'GroupIndex'.
    /// @param GroupIndex Starting point.
    /// @param rGroupSetVector Vector that holds all groups.
    /// @param [out] rMergedGroups Output vector of groups.
    /// @param rBoundaryIndices Contains all indices of element that are on the boundary.
    /// @param PartitionDir Direction along the partition was performed. Direction with "n_element_max".
    /// @param rStates Global classification vector.
    /// @param rVisited Vector<bool> for all elements.
    void GroupFill(IndexType GroupIndex, GroupSetVectorType& rGroupSetVector, GroupSetVectorType& rMergedGroups,
        const BoundaryIndicesVectorType& rBoundaryIndices, IndexType PartitionDir, StatusVectorType& rStates, BoolVectorType& rVisited ) const;

    /// @brief Performs a local ray tracing of two adjacent elements. On that is part of a new group and its trimmed neighbour.
    ///        We up to 10 rays, from the center of the first element towards the intersected triangles of the cut element. Based on the orientation of each triangle
    ///        the ray indicates inside or outside. The majoriy decides about the classification of the tested element. Returns +1 if inside and -1 if outside.
    /// @param Index
    /// @param NextIndex
    /// @param rLowerOffset
    /// @param rUpperOffset
    /// @return int
    int GetIsInsideCount( IndexType Index, IndexType NextIndex, const PointType& rLowerOffset, const PointType& rUpperOffset ) const;

    /// @brief Returns the next index for a given direction.
    /// @param Direction
    /// @param Index Current index.
    /// @return int Next index.
    int GetNextIndex( IndexType Direction, IndexType Index ) const;

    /// @brief Returns the next index for a given direction. It also outputs 'offsets', which extend the size of the next element towards the direction
    ///        of the start element to catch trimming surfaces that are exactly on the boundary.
    /// @param Direction
    /// @param Index
    /// @param [out] rLowerBoundOffset
    /// @param [out] rUpperBoundOffset
    /// @return int NextIndex.
    int GetNextIndex( IndexType Direction, IndexType Index, PointType& rLowerBoundOffset, PointType& rUpperBoundOffset ) const;

    /// @brief Returns the next index for a given direction and partition. It also outputs 'offsets', which extend the size of the next element towards the direction
    ///        of the start element to catch trimming surfaces that are exactly on the boundary.
    /// @param Direction
    /// @param Index
    /// @param rPartition
    /// @param [out] rLowerBoundOffset
    /// @param [out] rUpperBoundOffset
    /// @return int NextIndex.
    int GetNextIndex(IndexType Direction, IndexType Index, const PartitionBoxType& rPartition,
        PointType& rLowerBoundOffset, PointType& rUpperBoundOffset) const;

    ///@}
    ///@name Members
    ///@{

    const BRepOperator* mpBrepOperator;
    Mapper mMapper;

    // The following parameters are global values, w.r.t. to the background mesh.
    const PointType mLowerBound;
    const PointType mUpperBound;
    const Vector3i mNumberOfElements;
    PointType mDelta;

}; // End class FloodFill
///@}
///@} QuESo Classes
} // End queso namespace

#endif // FLOOD_FILL_INCLUDE_H