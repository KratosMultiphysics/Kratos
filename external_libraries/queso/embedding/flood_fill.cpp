// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// External includes
#include <omp.h>
#include <map>
#include <cstdlib>
//// Project includes
#include "queso/embedding/flood_fill.h"
#include "queso/embedding/brep_operator.h"

namespace queso {

typedef FloodFill::StatusVectorType StatusVectorType;
typedef FloodFill::StatusVectorType StatusVectorType;
typedef FloodFill::GroupSetVectorType GroupSetVectorType;

Unique<StatusVectorType> FloodFill::ClassifyElements() const {
    GroupSetVectorType element_groups;
    return ClassifyElements(element_groups);
}

std::pair<Unique<StatusVectorType>, Unique<GroupSetVectorType>> FloodFill::ClassifyElementsForTest() const {
    GroupSetVectorType element_groups;
    auto p_states = ClassifyElements(element_groups);
    return std::make_pair( std::move(p_states), MakeUnique<GroupSetVectorType>(element_groups) );
}

Unique<StatusVectorType> FloodFill::ClassifyElements(GroupSetVectorType& rGroupsOutput) const {
    const IndexType total_num_elements = mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];

    // Set all elements as outside.
    StatusVectorType states(total_num_elements);
    std::fill(states.begin(), states.end(), IntersectionStatus::Outside );

    // Get max num elements. Partition will happen along side with max_elements.
    IndexType max_num_elements_per_dir = std::max<IndexType>( std::max<IndexType>(
            mNumberOfElements[0], mNumberOfElements[1] ), mNumberOfElements[2] );

    // Direction along we apply partition.
    IndexType partition_index = 0;
    std::array<IndexType, 2> directions;
    if( mNumberOfElements[0] == max_num_elements_per_dir ){
        partition_index = 0;
        directions[0] = 0; // +x
        directions[1] = 1; // -x
    } else if ( mNumberOfElements[1] == max_num_elements_per_dir ) {
        partition_index = 1;
        directions[0] = 2; // +y
        directions[1] = 3; // -y
    } else if ( mNumberOfElements[2] == max_num_elements_per_dir ) {
        partition_index = 2;
        directions[0] = 4; // +z
        directions[1] = 5; // -z
    }

    // Partition domain in "num_threads" partitions along the direction with 'max_num_elements_per_dir';
    std::vector<PartitionBoxType> partitions;
    IndexType num_threads = 1;
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }

    const IndexType partition_size = std::max<IndexType>( static_cast<IndexType>(std::ceil( static_cast<double>(mNumberOfElements[partition_index]) /
        static_cast<double>(num_threads) ) ), 1);

    Vector3i partition_lower_bound{0, 0, 0};
    Vector3i partition_upper_bound = mNumberOfElements;
    for( IndexType i = 0; i < num_threads; ++i ){
        if( partition_size*i < mNumberOfElements[partition_index] ){
            partition_lower_bound[partition_index] = partition_size*i;
            partition_upper_bound[partition_index] = std::min<IndexType>(partition_size*(i+1), mNumberOfElements[partition_index]);
            partitions.push_back( std::make_pair(partition_lower_bound,  partition_upper_bound) );
        }
    }

    // Start filling.
    GroupSetVectorType groups;
    PartitionedFill(groups, partitions, states);

    // Merge groups from all partitions.
    MergeGroups( groups, rGroupsOutput, partition_index, partitions, states );

    // Mark states
    for( auto& r_group : rGroupsOutput ){
        int inside_count = std::get<2>(r_group);
        auto state = (inside_count > 0) ? IntersectionStatus::Inside : IntersectionStatus::Outside;
        for( auto group_it = std::get<1>(r_group).begin(); group_it !=  std::get<1>(r_group).end(); ++group_it){
            states[*group_it] = state;
        }
    }

    return MakeUnique<StatusVectorType>(states);
}

void FloodFill::PartitionedFill(GroupSetVectorType& rGroupSetVector, PartitionBoxVectorType& rPartitions, StatusVectorType& rStates) const {
    // Loop through current partition
     const IndexType total_num_elements = mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];
    BoolVectorType visited(total_num_elements, false);
    #pragma omp parallel for firstprivate(visited) schedule(static, 1)
    for( int p_i = 0; p_i < static_cast<int>(rPartitions.size()); ++p_i){
        const auto& partition = rPartitions[p_i];
        for( IndexType i = partition.first[0]; i < partition.second[0]; ++i ){
            for( IndexType j = partition.first[1]; j < partition.second[1]; ++j ) {
                for( IndexType k = partition.first[2]; k < partition.second[2]; ++k ) {
                    const IndexType index = mMapper.GetVectorIndexFromMatrixIndices(i, j, k);
                    if( !visited[index] ) { // Unvisited
                        GroupSetType new_group; // Tuple: get<0> -> partition_index, get<1> -> index_set, get<2> -> is_inside_count.
                        std::get<0>(new_group) = p_i; // Partition index
                        Fill(index, new_group, partition, rStates, visited);
                        if( std::get<1>(new_group).size() > 0 ){
                            #pragma omp critical
                            rGroupSetVector.push_back(new_group);
                        }
                    }
                }
            }
        }
    }
}

void FloodFill::Fill(IndexType Index, GroupSetType& rGroupSet, const PartitionBoxType& rPartition,
        StatusVectorType& rStates, BoolVectorType& rVisited ) const {

    // Set Index as visited
    rVisited[Index] = true;
    const auto box = mMapper.GetBoundingBoxXYZFromIndex(Index);
    // Only start filling if current element is not trimmed.
    if( mpBrepOperator->IsTrimmed(box.first, box.second) ){
        rStates[Index] = IntersectionStatus::Trimmed;
    } else {
        // Tuple: get<0> -> partition_index, get<1> -> index_set, get<2> -> is_inside_count.
        std::get<1>(rGroupSet).insert(Index);

        // If box is not trimmed, run flood fill.
        IndexStackType index_stack;
        index_stack.push( Index );
        std::array<int, 6> new_indices;
        while( !index_stack.empty() ){
            /// 0:+x, 1:-x, 2:+y, 3:-y, 4:+z, 5:-z
            const IndexType current_index = index_stack.top();
            for( IndexType direction = 0; direction < 6; ++direction){
                new_indices[direction] = Move(current_index, direction, rGroupSet, rPartition, rStates, rVisited );
            }

            index_stack.pop();
            for( IndexType direction = 0; direction < 6; ++direction){
                if( new_indices[direction] > -1 ){
                    index_stack.push(new_indices[direction]);
                }
            }
        }
    }

}

int FloodFill::Move(IndexType Index, IndexType Direction, GroupSetType& rGroupSet,
        const PartitionBoxType& rPartition, StatusVectorType& rStates, BoolVectorType& rVisited ) const {

    const IndexType index = Index;
    const IndexType max_num_elements = mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];

    PointType lower_perturb{0.0, 0.0, 0.0};
    PointType upper_perturb{0.0, 0.0, 0.0};
    int next_index = GetNextIndex(Direction, index, rPartition, lower_perturb, upper_perturb);

    // Check if out-of-range
    if( next_index >= static_cast<int>(max_num_elements) || next_index < 0 ) {
        return -1;
    }

    // If next box is trimmed, add inside count.
    // Note that next box is slightly shifted towards original box to capture trims directly at the boundary.
    auto box_next = mMapper.GetBoundingBoxXYZFromIndex(next_index);
    if( mpBrepOperator->IsTrimmed( Math::Add(box_next.first, lower_perturb) , Math::Add(box_next.second, upper_perturb) )){
        // Tuple: get<0> -> partition_index, get<1> -> index_set, get<2> -> is_inside_count.
        std::get<2>(rGroupSet) += GetIsInsideCount(index, next_index, lower_perturb, upper_perturb);
        return -1;
    }

    // Already visited.
    if( rVisited[next_index] ) {
        return -1;
    }

    // Is trimmed.
    if( mpBrepOperator->IsTrimmed(box_next.first, box_next.second) ){
        rVisited[next_index] = true;
        rStates[next_index] = IntersectionStatus::Trimmed;
        return -1;
    }

    // Add next_index to current group.
    rVisited[next_index] = true;
    // Tuple: get<0> -> partition_index, get<1> -> index_set, get<2> -> is_inside_count.
    std::get<1>(rGroupSet).insert(next_index);

    return next_index;
}

void FloodFill::MergeGroups(GroupSetVectorType& rGroups, GroupSetVectorType& rMergedGroup,
        IndexType PartitionDir, PartitionBoxVectorType& rPartitions, StatusVectorType& rStates) const {

    const IndexType num_groups = rGroups.size();

    // Mapping of directions: 0:+x, 1:-x, 2:+y, 3:-y, 4:+z, 5:-z
    const std::array<IndexType, 2> walk_directions = {2*PartitionDir, (2*PartitionDir)+1};

    // Compute bounding box for each group (only along PartitionDir).
    std::vector<Partition1DBoxType> group_bounding_boxes(num_groups);
    std::fill(group_bounding_boxes.begin(), group_bounding_boxes.end(),
        std::make_pair( static_cast<int>(mNumberOfElements[PartitionDir]+1),  -1 ) );

    // We compute the indices of all element are located on the boundary of a group.
    BoundaryIndicesVectorType group_boundary_indices(num_groups);
    for( IndexType i = 0; i < num_groups; ++i){
        group_boundary_indices[i].resize(2);
    }

    #pragma omp parallel for
    for( int group_index = 0; group_index < static_cast<int>(num_groups);  ++group_index){
        auto& group_set = rGroups[group_index];
        // Tuple: get<0> -> partition_index, get<1> -> index_set, get<2> -> is_inside_count.
        const auto& index_set = std::get<1>(group_set);
        const auto partition_index = std::get<0>(group_set);

        auto& group_bounding_box = group_bounding_boxes[group_index];
        auto& boundary_indices = group_boundary_indices[group_index];
        PointType lower_offset, upper_offset;

        // Loop over elements in active group.
        for(auto& index : index_set ){
            const auto indices = mMapper.GetMatrixIndicesFromVectorIndex(index);

            //// Find boundary indices. We consider only boundary indices that also bound the bounding box of each element.
            // Upper bound
            if( static_cast<int>(indices[PartitionDir]) >= std::get<1>( group_bounding_box ) ){
                auto next_index = GetNextIndex(walk_directions[0], index, lower_offset, upper_offset );
                auto box_next = mMapper.GetBoundingBoxXYZFromIndex(next_index);
                if( !mpBrepOperator->IsTrimmed(Math::Add(box_next.first, lower_offset), Math::Add(box_next.second, upper_offset) ) ) {
                    if( static_cast<int>(indices[PartitionDir]) > std::get<1>( group_bounding_box ) ){
                        std::get<1>( group_bounding_box ) = indices[PartitionDir];
                        boundary_indices[1].clear();
                    }
                    boundary_indices[1].insert(index);
                }
            }
            // Lower bound
            if( static_cast<int>(indices[PartitionDir]) <= std::get<0>( group_bounding_box ) ){
                auto next_index = GetNextIndex(walk_directions[1], index, lower_offset, upper_offset );
                auto box_next = mMapper.GetBoundingBoxXYZFromIndex(next_index);
                if( !mpBrepOperator->IsTrimmed( Math::Add(box_next.first, lower_offset), Math::Add(box_next.second, upper_offset) ) ){
                    if( static_cast<int>(indices[PartitionDir]) < std::get<0>( group_bounding_box ) ){
                        std::get<0>( group_bounding_box ) = indices[PartitionDir];
                        boundary_indices[0].clear();
                    }
                    boundary_indices[0].insert(index);
                }
            }

            //// We add the GetIsInsideCount() to all elements at the partition boundaries.
            // Lower bound.
            if( indices[PartitionDir] == std::get<1>(rPartitions[partition_index])[PartitionDir]-1 ) {
                auto next_index = GetNextIndex(walk_directions[0], index, lower_offset, upper_offset );
                auto box_next = mMapper.GetBoundingBoxXYZFromIndex(next_index);
                if( mpBrepOperator->IsTrimmed( Math::Add(box_next.first, lower_offset), Math::Add(box_next.second, upper_offset) ) ){
                    std::get<2>(group_set) += GetIsInsideCount(index, next_index, lower_offset, upper_offset);
                }
            }
            // Upper bound.
            if( indices[PartitionDir] == std::get<0>(rPartitions[partition_index])[PartitionDir] ) {
                auto next_index = GetNextIndex(walk_directions[1], index, lower_offset, upper_offset );
                auto box_next = mMapper.GetBoundingBoxXYZFromIndex(next_index);
                if( mpBrepOperator->IsTrimmed( Math::Add( box_next.first, lower_offset), Math::Add( box_next.second, upper_offset) ) ){
                    std::get<2>(group_set) += GetIsInsideCount(index, next_index, lower_offset, upper_offset);
                }
            }
        }

    }

    // Run group fill. Must be single-thread.
    BoolVectorType visited(num_groups, false);
    for( IndexType group_index = 0; group_index < num_groups; ++group_index){
        if( !visited[group_index] ){
            GroupFill(group_index, rGroups, rMergedGroup, group_boundary_indices, PartitionDir, rStates, visited);
        }
    }

}

void FloodFill::GroupFill(IndexType GroupIndex, GroupSetVectorType& rGroupSetVector, GroupSetVectorType& rMergedGroups,
        const BoundaryIndicesVectorType& rBoundaryIndices, IndexType PartitionDir, StatusVectorType& rStates, BoolVectorType& rVisited ) const {

    // Mapping of directions: 0:+x, 1:-x, 2:+y, 3:-y, 4:+z, 5:-z
    const std::array<IndexType, 2> walk_directions = {2*PartitionDir, (2*PartitionDir)+1};

    // Set index as visited
    rVisited[GroupIndex] = true;

    // Add GroupIndex to stack. This is where group flood fill will start,
    IndexStackType index_stack;
    index_stack.push( GroupIndex );

    rMergedGroups.push_back( rGroupSetVector[GroupIndex] );

    // Map forward to backward direction.
    const std::array<IndexType, 6> map_direction = {1, 0, 3, 2, 5, 4};

    // Loop over stack and perform group flood fill.
    while( !index_stack.empty() ){
        const IndexType current_group_index = index_stack.top();
        index_stack.pop();
        const IndexType partition_index = std::get<0>(rGroupSetVector[current_group_index]);
        const auto& current_boundary_indices = rBoundaryIndices[current_group_index];

        for( IndexType i = 0; i < rGroupSetVector.size(); ++i ){
            auto other_partition_index = std::get<0>(rGroupSetVector[i]);
            if( !rVisited[i] && std::abs( static_cast<int>(partition_index-other_partition_index)) <=1 ){
                const auto& other_boundary_indices = rBoundaryIndices[i];
                bool are_neighbours = false;
                for( auto direction : walk_directions){
                    if( are_neighbours ){
                        break;
                    }
                    for( auto& iii : current_boundary_indices[map_direction[direction] % 2] ) { // -1{
                        auto next_index = GetNextIndex(direction, iii);
                        if( other_boundary_indices[direction % 2].find(next_index) != other_boundary_indices[direction % 2].end() ){
                            are_neighbours = true;
                            break;
                        }
                    }
                }

                if( are_neighbours ) {
                    // Add group to merged groups.
                    auto& merged_group = rMergedGroups[rMergedGroups.size()-1];
                    auto &new_set = std::get<1>(rGroupSetVector[i]);

                    std::get<1>(merged_group).insert( new_set.begin(), new_set.end() );
                    std::get<2>(merged_group) += std::get<2>(rGroupSetVector[i]);

                    rVisited[i] = true;
                    index_stack.push(i);
                }
            }
        }
    }
}

int FloodFill::GetIsInsideCount( IndexType Index, IndexType NextIndex, const PointType& rLowerOffset, const PointType& rUpperOffset ) const{
    const auto box_current = mMapper.GetBoundingBoxXYZFromIndex(Index);
    const auto box_next = mMapper.GetBoundingBoxXYZFromIndex(NextIndex);
    const PointType center_box = Math::AddAndMult(0.5, box_current.first, box_current.second);

    if( mpBrepOperator->OnBoundedSideOfClippedSection(center_box, Math::Add(box_next.first, rLowerOffset) , Math::Add(box_next.second, rUpperOffset) ) ) {
        return 1;
    } else {
        return -1;
    }
}

int FloodFill::GetNextIndex( IndexType Direction, IndexType Index ) const {
    PointType dummy_1, dummy_2;
    const PartitionBoxType partition_box = std::make_pair(Vector3i{0, 0, 0}, mNumberOfElements);
    return GetNextIndex(Direction, Index, partition_box, dummy_1, dummy_2);
}

int FloodFill::GetNextIndex( IndexType Direction, IndexType Index, PointType& rLowerBoundOffset, PointType& rUpperBoundOffset ) const {
    const PartitionBoxType partition_box = std::make_pair(Vector3i{0, 0, 0}, mNumberOfElements);
    return GetNextIndex(Direction, Index, partition_box, rLowerBoundOffset, rUpperBoundOffset);
}

int FloodFill::GetNextIndex(IndexType Direction, IndexType Index, const PartitionBoxType& rPartition, PointType& rLowerBoundOffset, PointType& rUpperBoundOffset) const {
    const auto indices = mMapper.GetMatrixIndicesFromVectorIndex(Index);
    Vector3i next_indices = indices;
    rLowerBoundOffset = {0.0, 0.0, 0.0};
    rUpperBoundOffset = {0.0, 0.0, 0.0};
    const double tolerance = 10*RelativeSnapTolerance(mDelta, SNAPTOL);
    switch(Direction){
        case 0:
            if( indices[0] < rPartition.second[0]-1 ){ next_indices[0] += 1; }
            else { return -1; }
            rLowerBoundOffset = {-tolerance, 0.0, 0.0};
            break;
        case 1:
            if( indices[0] > rPartition.first[0] ){ next_indices[0] -= 1; }
            else { return -1; }
            rUpperBoundOffset = {tolerance, 0.0, 0.0};
            break;
        case 2:
            if( indices[1] < rPartition.second[1]-1 ){ next_indices[1] += 1; }
            else { return -1; }
            rLowerBoundOffset = {0.0, -tolerance, 0.0};
            break;
        case 3:
            if( indices[1] > rPartition.first[1] ){ next_indices[1] -= 1; }
            else { return -1; }
            rUpperBoundOffset = {0.0, tolerance, 0.0};
            break;
        case 4:
            if( indices[2] < rPartition.second[2]-1 ){ next_indices[2] += 1; }
            else { return -1; }
            rLowerBoundOffset = {0.0, 0.0, -tolerance};
            break;
        case 5:
            if( indices[2] > rPartition.first[2] ){ next_indices[2] -= 1; }
            else { return -1; }
            rUpperBoundOffset = {0.0, 0.0, tolerance};
            break;
        default:
            QuESo_ERROR << " Direction is out-of-range.\n";
    }

    return mMapper.GetVectorIndexFromMatrixIndices(next_indices[0], next_indices[1], next_indices[2]);
}

} // End namespace queso