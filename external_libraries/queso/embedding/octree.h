//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#ifndef OCTREE_INCLUDE_H
#define OCTREE_INCLUDE_H

//// STL includes
#include <array>

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/quadrature/single_element.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  Octree
 * @author Manuel Messmer
 * @tparam TOperator: Required member operations:
 *                   -GetIntersectionState(const PointType& rLowerBound, const PointType& rUpperBound)
 *                   -IsInsideTrimmedDomain(const PointType& rPoint)
*/
template<typename TOperator>
class Octree {
private:
    ///@name Type Definitions
    ///@{

    /**
     * @class  Octree::Node
     * @author Manuel Messmer
     * @brief Node to be used in octree. Each node represents an AABB.
    */
    class Node {
    public:
        ///@name Life Cycle
        ///@{

        /// @brief Contructor
        /// @param rBoundsXYZ Bounds of AABB in physical space.
        /// @param rBoundsUVW Bounds of AABB in parametric space..
        /// @param Status Options: 'Trimmed' or 'Inside'.
        /// @param Level Refinement Level (default - 0).
        Node( const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW, IntersectionStatusType Status, IndexType Level = 0UL) :
            mBoundsXYZ(rBoundsXYZ), mBoundsUVW(rBoundsUVW), mStatus(Status), mLevel(Level)
        {
            mChildren = {nullptr};
            mNumChildren = 0UL;
        }

        ///@}
        ///@name Operations
        ///@{

        /// @brief Recursive function (walks through octree) to split node into its 8 children.
        /// @param MinLevel Minimum refinement level (for all nodes).
        /// @param MaxLevel Maximum refinement level (for trimmed nodes).
        /// @param pOperator Operator to perform: GetIntersectonState().
        void Refine(IndexType MinLevel, IndexType MaxLevel, const TOperator* pOperator);

        /// @brief Recursive function (walks through octree) to get integration points.
        /// @tparam TElementType
        /// @details Distributes Gauss points (according to rOrder) and adds all points to pPoints that are
        ///          inside the domain.
        /// @param[out] pPoints IntegrationPointVectorType
        /// @param rOrder Order of Gauss quadrature.
        /// @param pOperator Operator to perfrom Inside/Outside test.
        template<typename TElementType>
        void GetIntegrationPoints(typename TElementType::IntegrationPointVectorType* pPoints, const Vector3i& rOrder, const TOperator* pOperator) const {
            if( this->IsLeaf() ){
                std::vector<typename TElementType::IntegrationPointType> integration_points_tmp{};
                // Note that QuadratureSingleElement::AssembleIPs clears integration_points_tmp.
                QuadratureSingleElement<TElementType>::AssembleIPs(integration_points_tmp, mBoundsUVW.first, mBoundsUVW.second, rOrder);
                if( mStatus == IntersectionStatus::Inside )
                    pPoints->insert(pPoints->end(), integration_points_tmp.begin(), integration_points_tmp.end());
                else {
                    for( auto& point : integration_points_tmp){
                        const auto tmp_point = Mapping::PointFromParamToGlobal(point.data(), mBoundsXYZ, mBoundsUVW);
                        if( pOperator->IsInsideTrimmedDomain( tmp_point ) ){
                            pPoints->push_back(point);
                        }
                    }
                }
            } else {
                for( IndexType i = 0; i < 8UL; ++i){
                    if( mChildren[i] ){ // If not nullptr
                        mChildren[i]->template GetIntegrationPoints<TElementType>(pPoints, rOrder, pOperator);
                    }
                }
            }
        }

        /// @brief Recursive function (walks through octree) to get total number of leaf nodes.
        /// @param[out] rValue // Return value
        void NumberOfLeafs(IndexType& rValue) const;

        /// @brief Recursive function (walks through octree) to get total number of nodes.
        /// @param[out] rValue // Return value
        void NumberOfNodes(IndexType& rValue) const;

        ///@}
    private:
        ///@name Private Operations
        ///@{

        /// @brief Create a new child node if aabb of this child node is NOT classified as Outside.
        /// @param MinLevel Minimum refinement level (for all nodes).
        /// @param MaxLevel Maximum refinement level (for trimmed nodes).
        /// @param ChildIndex Index of child node in mChildren[ChildIndex].
        /// @param rBoundsXYZ Bounds of AABB of new child in physical space.
        /// @param rBoundsUVW Bounds of AABB of new child in parametric space..
        /// @param pOperator Operator to perfrom Inside/Outside test.
        void CreateNewNode(IndexType MinLevel, IndexType MaxLevel, IndexType ChildIndex, const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW, const TOperator* pOperator);

        /// @brief Returns true if this node is a leaf node.
        /// @return bool
        bool IsLeaf() const;

        ///@}
        ///@name Private Member Variables
        ///@{
        BoundingBoxType mBoundsXYZ;
        BoundingBoxType mBoundsUVW;
        IntersectionStatus mStatus;
        std::array<Unique<Node>, 8> mChildren{};
        SizeType mLevel;
        SizeType mNumChildren;
        ///@}
    }; ///@} End Node Class
    ///@} End QuESo Classes

public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructor
    /// @param pOperator Operator to perform GetIntersectionState() and IsInsideTrimmedDomain().
    /// @param rBoundsXYZ Bounds of AABB of Root Node in physical space.
    /// @param rBoundsUVW Bounds of AABB of Root Node in parametric space.
    /// @param rParameters QuESo Parameters.
    Octree(const TOperator* pOperator, const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW)
        : mpOperator(pOperator) {

        // Create new root node.
        mpRoot = MakeUnique<Node>(rBoundsXYZ, rBoundsUVW, IntersectionStatus::Trimmed, 0UL);
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Refine octree up to given levels.
    /// @param MinLevel Minimum refinement level (for all nodes).
    /// @param MaxLevel Maximum refinement level (for trimmed nodes).
    void Refine(IndexType MinLevel, IndexType MaxLevel);

    /// @brief Returns total number of leaf nodes in octree.
    /// @return SizeType
    SizeType NumberOfLeafs() const;

    /// @brief Returns total number of nodes in octree.
    /// @return SizeType
    SizeType NumberOfNodes() const;

    /// @brief Returns ptr to integration points that are constructed on the leafs of the octree.
    ///        Standard Gauss quadrature rules (according to rOrder) are constructed on each leaf node.
    ///        But only points that are inside the domain are considered.
    ///        Also see: AddIntegrationPoints()
    /// @tparam TElementType
    /// @param rOrder Order of Gauss quadrature.
    /// @return IntegrationPointVectorPtrType
    template<typename TElementType>
    Unique<std::vector<typename TElementType::IntegrationPointType>> pGetIntegrationPoints(const Vector3i& rOrder) const {
        auto p_points = MakeUnique<std::vector<typename TElementType::IntegrationPointType>>();
        p_points->reserve(NumberOfLeafs());
        mpRoot->template GetIntegrationPoints<TElementType>(p_points.get(), rOrder, mpOperator );
        return p_points;
    }

    /// @brief Add integration points to rPoints. Points are constructed on the leafs of the octree.
    ///        Standard Gauss quadrature rules (according to rOrder) are constructed on each leaf node.
    ///        But only points that are inside the domain are considered.
    ///        Also see: pGetIntegrationPoints().
    /// @tparam TElementType
    /// @param rPoints Vector of integration points.
    /// @param rOrder Order of Gauss quadrature.
    template<typename TElementType>
    void AddIntegrationPoints(std::vector<typename TElementType::IntegrationPointType>& rPoints, const Vector3i& rOrder) const {
        rPoints.reserve(NumberOfLeafs());
        mpRoot->template GetIntegrationPoints<TElementType>(&rPoints, rOrder, mpOperator );
    }

    /// @brief Returns current refinement level of leaf nodes that are classified as inside.
    /// @return IndexType
    IndexType MinRefinementLevel(){
        return mMinLevel;
    }

    /// @brief Returns current refinement level of leaf nodes that are classified as trimmed.
    /// @return
    IndexType MaxRefinementLevel(){
        return mMaxLevel;
    }
    ///@}
private:

    ///@name Private Member Variables
    ///@{
    Unique<Node> mpRoot;
    const TOperator* mpOperator; // No ownership
    IndexType mMinLevel = 0;
    IndexType mMaxLevel = 0;
    ///@}

}; // End class Octree
///@} QuESo Classes
} // End queso namespace

#endif // OCTREE_INCLUDE_H