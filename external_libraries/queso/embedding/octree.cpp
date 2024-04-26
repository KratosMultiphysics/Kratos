// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// Project includes
#include "embedding/octree.h"
#include "embedding/trimmed_domain_base.h"

namespace queso {

///
/// Function Definition of Octree<TOperator>::Node
///

template<typename TOperator>
void Octree<TOperator>::Node::Refine(IndexType MinLevel, IndexType MaxLevel, const TOperator* pOperator){
    const PointType& r_lower_bound_xyz = mBoundsXYZ.first;
    const PointType& r_upper_bound_xyz = mBoundsXYZ.second;
    const PointType& r_lower_bound_uvw = mBoundsUVW.first;
    const PointType& r_upper_bound_uvw = mBoundsUVW.second;
    if( this->IsLeaf() ){
        const auto delta_xyz = (r_upper_bound_xyz - r_lower_bound_xyz) * 0.5;
        const PointType delta_x_xyz(delta_xyz[0], 0.0, 0.0);
        const PointType delta_y_xyz(0.0, delta_xyz[1], 0.0);
        const PointType delta_z_xyz(0.0, 0.0, delta_xyz[2]);
        const auto delta_uvw = (r_upper_bound_uvw - r_lower_bound_uvw) * 0.5;
        const PointType delta_x_uvw(delta_uvw[0], 0.0, 0.0);
        const PointType delta_y_uvw(0.0, delta_uvw[1], 0.0);
        const PointType delta_z_uvw(0.0, 0.0, delta_uvw[2]);

        //       d_________c
        //      /        /|                 y
        //     /        / |                ´|`
        //   h/_______g/  |b                |-->x
        //    |  a     |  /                /
        //    |        | /                Z
        //    |________|/
        //    e        f
        //
        if( (mLevel < MinLevel) || (mLevel < MaxLevel && mStatus == IntersectionStatus::Trimmed) ){
            // Corner a
            CreateNewNode(MinLevel, MaxLevel, 0, std::make_pair(r_lower_bound_xyz, r_lower_bound_xyz+delta_xyz),
                                                 std::make_pair(r_lower_bound_uvw, r_lower_bound_uvw+delta_uvw), pOperator);
            // Corner b (a+delta_x)
            CreateNewNode(MinLevel, MaxLevel, 1, std::make_pair(r_lower_bound_xyz+delta_x_xyz, r_lower_bound_xyz+delta_x_xyz+delta_xyz),
                                                 std::make_pair(r_lower_bound_uvw+delta_x_uvw, r_lower_bound_uvw+delta_x_uvw+delta_uvw), pOperator);
            // Corner c (g-delta_z)
            CreateNewNode(MinLevel, MaxLevel, 2, std::make_pair(r_upper_bound_xyz-delta_z_xyz-delta_xyz, r_upper_bound_xyz-delta_z_xyz),
                                                 std::make_pair(r_upper_bound_uvw-delta_z_uvw-delta_uvw, r_upper_bound_uvw-delta_z_uvw), pOperator);
            // Corner d (a+delta_y)
            CreateNewNode(MinLevel, MaxLevel, 3, std::make_pair(r_lower_bound_xyz+delta_y_xyz, r_lower_bound_xyz+delta_y_xyz+delta_xyz),
                                                 std::make_pair(r_lower_bound_uvw+delta_y_uvw, r_lower_bound_uvw+delta_y_uvw+delta_uvw), pOperator);
            // Corner e (a+delta_z)
            CreateNewNode(MinLevel, MaxLevel, 4, std::make_pair(r_lower_bound_xyz+delta_z_xyz, r_lower_bound_xyz+delta_z_xyz+delta_xyz),
                                                 std::make_pair(r_lower_bound_uvw+delta_z_uvw, r_lower_bound_uvw+delta_z_uvw+delta_uvw), pOperator);
            // Corner f (g-delta_y)
            CreateNewNode(MinLevel, MaxLevel, 5, std::make_pair(r_upper_bound_xyz-delta_y_xyz-delta_xyz, r_upper_bound_xyz-delta_y_xyz),
                                                 std::make_pair(r_upper_bound_uvw-delta_y_uvw-delta_uvw, r_upper_bound_uvw-delta_y_uvw), pOperator);
            // Corner g
            CreateNewNode(MinLevel, MaxLevel, 6, std::make_pair(r_upper_bound_xyz-delta_xyz, r_upper_bound_xyz),
                                                 std::make_pair(r_upper_bound_uvw-delta_uvw, r_upper_bound_uvw), pOperator);
            // Corner h (g-delta_x)
            CreateNewNode(MinLevel, MaxLevel, 7, std::make_pair(r_upper_bound_xyz-delta_x_xyz-delta_xyz, r_upper_bound_xyz-delta_x_xyz),
                                                 std::make_pair(r_upper_bound_uvw-delta_x_uvw-delta_uvw, r_upper_bound_uvw-delta_x_uvw), pOperator);

        }
    }
    else {
        for( IndexType i = 0; i < 8UL; ++i){
            if( mChildren[i] ){
                mChildren[i]->Refine(MinLevel, MaxLevel, pOperator);
            }
        }
    }

}

template<typename TOperator>
void Octree<TOperator>::Node::GetIntegrationPoints(IntegrationPointVectorType* pPoints, const Vector3i& rOrder, const TOperator* pOperator) const{
    if( this->IsLeaf() ){

        IntegrationPointVectorType integration_points_tmp{};
        // Note that QuadratureSingleElement::AssembleIPs clears integration_points_tmp.
        QuadratureSingleElement::AssembleIPs(integration_points_tmp, mBoundsUVW.first, mBoundsUVW.second, rOrder);
        if( mStatus == IntersectionStatus::Inside )
            pPoints->insert(pPoints->end(), integration_points_tmp.begin(), integration_points_tmp.end());
        else {
            for( auto& point : integration_points_tmp){
                const auto tmp_point = Mapping::PointFromParamToGlobal(point, mBoundsXYZ, mBoundsUVW);
                if( pOperator->IsInsideTrimmedDomain( tmp_point ) ){
                    pPoints->push_back(point);
                }
            }
        }
    } else {
        for( IndexType i = 0; i < 8UL; ++i){
            if( mChildren[i] ){ // If not nullptr
                mChildren[i]->GetIntegrationPoints(pPoints, rOrder, pOperator);
            }
        }
    }
}

template<typename TOperator>
void Octree<TOperator>::Node::NumberOfLeafs(IndexType& rValue) const {
    if( this->IsLeaf() ){
        ++rValue;
    } else {
        for( IndexType i = 0; i < 8UL; ++i){
            if( mChildren[i] ){ // If not nullptr
                mChildren[i]->NumberOfLeafs(rValue);
            }
        }
    }
}

template<typename TOperator>
void Octree<TOperator>::Node::NumberOfNodes(IndexType& rValue) const {
    ++rValue;
    if( !this->IsLeaf() ){
        for( IndexType i = 0; i < 8UL; ++i){
            if( mChildren[i] ){ // If not nullptr
                mChildren[i]->NumberOfNodes(rValue);
            }
        }
    }
}

template<typename TOperator>
void Octree<TOperator>::Node::CreateNewNode(IndexType MinLevel, IndexType MaxLevel, IndexType ChildIndex, const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW, const TOperator* pOperator){
    const auto status = pOperator->GetIntersectionState(rBoundsXYZ.first, rBoundsXYZ.second);
    if( status != IntersectionStatus::Outside ){
        mChildren[ChildIndex] = MakeUnique<Node>(rBoundsXYZ, rBoundsUVW, status, mLevel+1);
        ++mNumChildren;
        mChildren[ChildIndex]->Refine(MinLevel, MaxLevel, pOperator);
    }
}

template<typename TOperator>
bool Octree<TOperator>::Node::IsLeaf() const {
    return (mNumChildren == 0UL);
}

///
/// Function Definition of Octree<TOperator>
///

template<typename TOperator>
void Octree<TOperator>::Refine(IndexType MinLevel, IndexType MaxLevel){
    QuESo_ERROR_IF( MinLevel > MaxLevel ) << "MinLevel must be smaller/equal than MaxLevel. "
        << "Given MinLevel: " << MinLevel << ", MaxLevel: " << MaxLevel << ".\n";
    mMinLevel = MinLevel;
    mMaxLevel = MaxLevel;
    mpRoot->Refine(MinLevel, MaxLevel, mpOperator);
}


template<typename TOperator>
SizeType Octree<TOperator>::NumberOfLeafs() const{
    SizeType number_of_leafs = 0UL;
    mpRoot->NumberOfLeafs(number_of_leafs);
    return number_of_leafs;
}

template<typename TOperator>
SizeType Octree<TOperator>::NumberOfNodes() const{
    SizeType number_of_nodes = 0UL;
    mpRoot->NumberOfNodes(number_of_nodes);
    return number_of_nodes;
}

template<typename TOperator>
typename Octree<TOperator>::IntegrationPointVectorPtrType Octree<TOperator>::pGetIntegrationPoints(const Vector3i& rOrder) const{
    auto p_points = MakeUnique<IntegrationPointVectorType>();
    p_points->reserve(NumberOfLeafs());
    mpRoot->GetIntegrationPoints(p_points.get(), rOrder, mpOperator );
    return p_points;
}

template<typename TOperator>
void Octree<TOperator>::AddIntegrationPoints(IntegrationPointVectorType& rPoints, const Vector3i& rOrder) const{
    rPoints.reserve(NumberOfLeafs());
    mpRoot->GetIntegrationPoints(&rPoints, rOrder, mpOperator );
}

// Explicit instantiation Octree
template class Octree<TrimmedDomainBase>;

} // End namespace queso