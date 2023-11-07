//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    klabra
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <type_traits>

// External includes

// Project includes
#include "search_structure.h"
#include "geometries/bounding_box.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
template<
std::size_t TDimension,
    class TPointType,
    class TPointerType,
    class TIteratorType,
    class TDistanceIteratorType,
    class TIteratorIteratorType = typename std::vector<TIteratorType>::iterator
    >
class TreeNode
{
public:

    ///@name Type Definitions

    /// Pointer definition of TreeNode
    KRATOS_CLASS_POINTER_DEFINITION(TreeNode);

    /// Define SizeType as std::size_t
    using SizeType = std::size_t;

    /// Define IndexType as std::size_t
    using IndexType = std::size_t;

    /// Define CoordinateType as double
    using CoordinateType = double;

    /// Define PointType as TPointType
    using PointType = TPointType;

    /// Define PointerType as TPointerType
    using PointerType = TPointerType;

    /// Define IteratorType as TIteratorType
    using IteratorType = TIteratorType;

    /// Define DistanceIteratorType as TDistanceIteratorType
    using DistanceIteratorType = TDistanceIteratorType;

    /// Define TreeNodeType as a TreeNode type with the given template arguments
    using TreeNodeType = TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType>;

    /// Define IteratorIteratorType as an iterator type for a vector of IteratorType
    using IteratorIteratorType = typename std::vector<IteratorType>::iterator;

    /// Define SearchStructureType as a SearchStructure type with the given template arguments
    using SearchStructureType = SearchStructure<IndexType, SizeType, CoordinateType, TIteratorType, IteratorIteratorType, TDimension>;

    virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const {}

    TreeNode() {}

    virtual ~TreeNode() {}

    virtual void SearchNearestPoint(PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance) {}

    virtual void SearchNearestPoint(PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance,
                                    SearchStructureType& Auxiliar) {}

    virtual void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                                DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults)
    {
        // must be implemented in derived classes.
        return;
    }

    virtual void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                                DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Auxiliar)
    {
        // must be implemented in derived classes.
        return;
    }

    virtual void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                                SizeType& NumberOfResults, SizeType const& MaxNumberOfResults)
    {
        // must be implemented in derived classes.
        return;
    }

    virtual void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                                SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Auxiliar)
    {
        // must be implemented in derived classes.
        return;
    }

    virtual void SearchInBox(PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& Results, SizeType& NumberOfResults,
                             SizeType const& MaxNumberOfResults ) // This corresponds with a AABB bounding-box
    {
        // must be implemented in derived classes.
        return;
    }

    // Note: Add OBB bounding-box and K-DOP bounding-box

    static IteratorType& NullIterator()
    {
        return msNull;
    }

    static PointerType& NullPointer()
    {
        return msNullPointer;
    }

    static TreeNode& NullLeaf()
    {
        return msNullLeaf;
    }

private:
    static IteratorType msNull;
    static PointerType msNullPointer;
    static TreeNode msNullLeaf;
};

template<std::size_t TDimension, class TPointType, class TPointerType, class TIteratorType, class TDistanceIteratorType, class TIteratorIteratorType>
typename TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType, TIteratorIteratorType>::IteratorType
TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType, TIteratorIteratorType>::msNull;

template<std::size_t TDimension, class TPointType, class TPointerType, class TIteratorType, class TDistanceIteratorType, class TIteratorIteratorType>
typename TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType, TIteratorIteratorType>::PointerType
TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType, TIteratorIteratorType>::msNullPointer;

template<std::size_t TDimension, class TPointType, class TPointerType, class TIteratorType, class TDistanceIteratorType, class TIteratorIteratorType>
TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType, TIteratorIteratorType>
TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType, TIteratorIteratorType>::msNullLeaf;

/**
 * @class Tree
 * @ingroup KratosCore
 * @brief A generic tree data structure for spatial partitioning.
 * @details This class implements a generic tree data structure for spatial partitioning.
 * @tparam TPartitionType The partitioning strategy type.
 * @author Carlos Labra
 */
template< class TPartitionType >
class Tree
{
public:
    ///@name Type Definitions
    ///@{

    // Helper template to extract ObjectType or default to void
    template <typename T, typename = void>
    struct GetObjectType {
        using type = void;
    };

    template <typename T>
    struct GetObjectType<T, std::void_t<typename T::ObjectType>> {
        using type = typename T::ObjectType;
    };

    /**
     * @class Partitions
     * @brief Class to represent partitions for the tree.
     */
    class Partitions
    {
    public:
        /**
         * @brief Constructor to initialize the number of partitions.
         * @param NumPartitions The number of partitions.
         */
        explicit Partitions( const std::size_t NumPartitions ) : mNumPartitions(NumPartitions) {}
        
        /// Destructor.
        ~Partitions() {};

        ///The number of partitions.
        std::size_t mNumPartitions;
    };

    /// Pointer definition of Tree
    KRATOS_CLASS_POINTER_DEFINITION(Tree);

    /// The partition type definition
    using PartitionType = TPartitionType;

    /// The leaf type definition
    using LeafType = typename PartitionType::LeafType;

    /// The point type definition
    using PointType = typename PartitionType::PointType;

    /// The iterator type definition
    using IteratorType = typename PartitionType::IteratorType;

    /// The object type
    using ObjectType = typename GetObjectType<PointType>::type;

    /// The distance iterator type definition
    using DistanceIteratorType = typename PartitionType::DistanceIteratorType;

    /// The pointer type definition
    using PointerType = typename PartitionType::PointerType;

    /// The distance function type definition
    using DistanceFunction = typename PartitionType::DistanceFunction;

    /// Dimension definition
    static constexpr std::size_t Dimension = PartitionType::Dimension;

    /// The node type definition
    using NodeType = TreeNode<Dimension,PointType,PointerType,IteratorType,DistanceIteratorType> ;

    /// The coordinate type definition
    using CoordinateType = typename NodeType::CoordinateType;

    /// The size type definition
    using SizeType = typename NodeType::SizeType;

    /// The index type definition
    using IndexType = typename NodeType::IndexType;

    /// The search structure type definition
    using SearchStructureType = typename PartitionType::SearchStructureType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Tree object
     * @param itPointsBegin Iterator to the first point
     * @param itPointsEnd Iterator to the last point
     * @param BucketSize Size of the bucket
     */
    Tree(
        IteratorType itPointsBegin,
        IteratorType itPointsEnd,
        SizeType BucketSize = 1
        ) : mBucketSize(BucketSize),
            mitPointsBegin(itPointsBegin),
            mitPointsEnd(itPointsEnd)
    {
        if(mitPointsBegin == mitPointsEnd)
            return;

        // The BB points
        auto& r_high_point = BoundingBoxHighPoint();
        auto& r_low_point = BoundingBoxLowPoint();

        for(SizeType i = 0 ; i < Dimension ; i++) {
            r_high_point[i] = (**mitPointsBegin)[i];
            r_low_point[i] = (**mitPointsBegin)[i];
        }

        for(IteratorType point_iterator = mitPointsBegin ; point_iterator != mitPointsEnd ; point_iterator++) {
            for(SizeType i = 0 ; i < Dimension ; i++) {
                if((**point_iterator)[i] > r_high_point[i]) {
                    r_high_point[i] = (**point_iterator)[i];
                } else if((**point_iterator)[i] < r_low_point[i]) {
                   r_low_point[i]  = (**point_iterator)[i];
                }
            }
        }

       mRoot = TPartitionType::Construct(mitPointsBegin, mitPointsEnd, r_high_point, r_low_point, mBucketSize);
    }

    /**
     * @brief Construct a new Tree object
     * @param itPointsBegin Iterator to the first point
     * @param itPointsEnd Iterator to the last point
     * @param Parts The partitions definition
     */
    Tree(
        IteratorType itPointsBegin,
        IteratorType itPointsEnd,
        Partitions Parts
        ) : mitPointsBegin(itPointsBegin), 
            mitPointsEnd(itPointsEnd)
    {
        if(mitPointsBegin == mitPointsEnd)
            return;

        // The BB points
        auto& r_high_point = BoundingBoxHighPoint();
        auto& r_low_point = BoundingBoxLowPoint();

        const SizeType NumPoints = SearchUtils::PointerDistance(mitPointsBegin,mitPointsEnd);
        mBucketSize = static_cast<std::size_t>( (double) NumPoints / (double) Parts.mNumPartitions ) + 1;

        r_high_point = **mitPointsBegin;
        r_low_point = **mitPointsBegin;
        for(IteratorType point_iterator = mitPointsBegin ; point_iterator != mitPointsEnd ; point_iterator++) {
            for(SizeType i = 0 ; i < Dimension ; i++) {
                if((**point_iterator)[i] > r_high_point[i]) {
                    r_high_point[i] = (**point_iterator)[i];
                } else if((**point_iterator)[i] < r_low_point[i]) {
                    r_low_point[i] = (**point_iterator)[i];
                }
            }
        }

        mRoot = TPartitionType::Construct(mitPointsBegin, mitPointsEnd, r_high_point, r_low_point, mBucketSize);
    }

    /// Destructor.
    virtual ~Tree()
    {
        delete mRoot;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Check if a point exists in the tree within a given tolerance.
     * @param ThisPoint The point to check.
     * @param Tolerance The tolerance for proximity check.
     * @return A pointer to the existing point or a null pointer.
     */
    PointerType ExistPoint(
        PointerType const& ThisPoint,
        CoordinateType const Tolerance = static_cast<CoordinateType>(10.0*DBL_EPSILON) 
        )
    {
        PointerType Result = *mitPointsBegin;
        CoordinateType ResultDistance = static_cast<CoordinateType>(DBL_MAX);
        // Searching the tree
        mRoot->SearchNearestPoint(*ThisPoint, Result, ResultDistance);
        if (ResultDistance<Tolerance*Tolerance)
            return Result;
        return NodeType::NullPointer();
    }

    /**
     * @brief Search for the nearest point to a given point.
     * @param ThisPoint The point for which to find the nearest point.
     * @return A pointer to the nearest point.
     */
    PointerType SearchNearestPoint(
        PointType const& ThisPoint,
        CoordinateType& rResultDistance
        )
    {
        PointerType Result = *mitPointsBegin;
        rResultDistance = static_cast<CoordinateType>(DBL_MAX); // DistanceFunction()(ThisPoint,**mitPointsBegin);

        // Searching the tree
        mRoot->SearchNearestPoint(ThisPoint, Result, rResultDistance);

        return Result;
    }

    /**
     * @brief Search for points within a given radius of a point.
     * @param ThisPoint The center point.
     * @param Radius The search radius.
     * @param Results Iterator to store the found points.
     * @param ResultsDistances Iterator to store the distances to found points.
     * @param MaxNumberOfResults Maximum number of results to return.
     * @return The number of points found within the specified radius.
     */
    PointerType SearchNearestPoint(PointType const& ThisPoint)
    {
        PointerType Result = *mitPointsBegin; // NULL ??
        CoordinateType rResultDistance = static_cast<CoordinateType>(DBL_MAX); // DistanceFunction()(ThisPoint,**mitPointsBegin);

        // Searching the tree
        mRoot->SearchNearestPoint(ThisPoint, Result, rResultDistance);

        return Result;
    }

    /**
     * @brief Search for points within a given radius of a point.
     * @param ThisPoint The center point.
     * @param Radius The search radius.
     * @param Results Iterator to store the found points.
     * @param ResultsDistances Iterator to store the distances to found points.
     * @param MaxNumberOfResults Maximum number of results to return.
     * @return The number of points found within the specified radius.
     */
    SizeType SearchInRadius(
        PointType const& ThisPoint,
        CoordinateType Radius,
        IteratorType Results,
        DistanceIteratorType ResultsDistances,
        SizeType MaxNumberOfResults
        )
    {
        // Using the square of radius for avoiding square root calculation during search
        const CoordinateType Radius2 = Radius * Radius;

        // Searching the tree
        SizeType NumberOfResults = 0;
        mRoot->SearchInRadius(ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults);

        return NumberOfResults;
    }

    /**
     * @brief Search for points within a given radius of a point.
     * @param ThisPoint The center point.
     * @param Radius The search radius.
     * @param Results Iterator to store the found points.
     * @param MaxNumberOfResults Maximum number of results to return.
     * @return The number of points found within the specified radius.
     */
    SizeType SearchInRadius(
        PointType const& ThisPoint,
        CoordinateType Radius,
        IteratorType Results,
        SizeType MaxNumberOfResults
        )
    {
        // Using the square of radius for avoiding square root calculation during search
        const CoordinateType Radius2 = Radius * Radius;

        // Searching the tree
        SizeType NumberOfResults = 0;
        mRoot->SearchInRadius(ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults);
        return NumberOfResults;
    }

    /**
     * @brief Search for points within a given axis-aligned box.
     * @param MinPointBox The minimum point of the bounding box.
     * @param MaxPointBox The maximum point of the bounding box.
     * @param Results Iterator to store the found points.
     * @param MaxNumberOfResults Maximum number of results to return.
     * @return The number of points found within the specified bounding box.
     */
    SizeType SearchInBox(
        PointType const& MinPointBox,
        PointType const& MaxPointBox,
        IteratorType Results,
        SizeType MaxNumberOfResults
        )
    {
        SizeType NumberOfResults = 0;
        mRoot->SearchInBox(MinPointBox,MaxPointBox,Results,NumberOfResults,MaxNumberOfResults);
        return NumberOfResults;
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get a reference to the low point of the bounding box.
     * @return A reference to the low point of the bounding box.
     */
    PointType& BoundingBoxLowPoint()
    {
        return mBoundingBox.GetMinPoint();
    }

    /**
     * @brief Get a reference to the high point of the bounding box.
     * @return A reference to the high point of the bounding box.
     */
    PointType& BoundingBoxHighPoint()
    {
        return mBoundingBox.GetMaxPoint();
    }

    /**
     * @brief Get the bounding box.
     * @details This function creates a bounding box using the low and high points and returns it.
     * @return The bounding box.
     */
    BoundingBox<PointType>& GetBoundingBox()
    {
        return mBoundingBox;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Tree";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Tree";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const
    {
        mRoot->PrintData(rOStream, "  ");
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    static LeafType msEmptyLeaf; /// The empty leaf

    ///@}
    ///@name Member Variables
    ///@{

    SizeType mBucketSize;                /// The bucket size

    BoundingBox<PointType> mBoundingBox; /// The bounding box of the tree

    IteratorType mitPointsBegin;         /// Iterator to the first point
    IteratorType mitPointsEnd;           /// Iterator to the last point

    NodeType* mRoot;                     /// The root node of the tree

    ///@}
    ///@name Private Operators:
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    Tree& operator=(Tree const& rOther);

    /// Copy constructor.
    Tree(Tree const& rOther);

    ///@}
}; // Class Tree

template< class TPartitionType >
typename Tree<TPartitionType>::LeafType Tree<TPartitionType>::msEmptyLeaf;

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TPartitionType>
inline std::istream& operator >> (std::istream& rIStream, Tree<TPartitionType>& rThis);

/// output stream function
template<class TPartitionType>
inline std::ostream& operator << (std::ostream& rOStream, const Tree<TPartitionType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.