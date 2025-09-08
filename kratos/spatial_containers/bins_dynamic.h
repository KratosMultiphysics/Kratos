//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Carlos Roig
//

#pragma once

// System includes
#include <array>
#include <cmath>
#include <algorithm>

// External includes

// Project includes
#include "tree.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

/**
 * @class BinsDynamic
 * @ingroup KratosCore
 * @brief A dynamic binning data structure template for organizing and querying points in multi-dimensional space.
 * @details The `BinsDynamic` class template provides a dynamic binning data structure for organizing and querying points in a multi-dimensional space. It is parameterized by the dimension of the space, the point type, the container type for storing points, and other optional template parameters for specifying point iterators, distance iterators, and distance functions.
 * This class inherits from `TreeNode` to leverage common functionalities for tree-based data structures.
 * @tparam TDimension            The dimensionality of the space.
 * @tparam TPointType            The type representing points in the space.
 * @tparam TContainerType        The container type for storing points.
 * @tparam TPointerType          The type of pointers or iterators to points in the container (default is inferred from `TContainerType`).
 * @tparam TIteratorType         The type of iterators for traversing the container (default is inferred from `TContainerType`).
 * @tparam TDistanceIteratorType The type of iterators for storing distance values (default is `std::vector<double>::iterator`).
 * @tparam TDistanceFunction     The type of the distance function used for querying (default is `Kratos::SearchUtils::SquaredDistanceFunction`).
 * @author Nelson Lafontaine
 * @author Carlos Roig
 */
template<
std::size_t TDimension,
    class TPointType,
    class TContainerType,
    class TPointerType = typename TContainerType::value_type,
    class TIteratorType = typename TContainerType::iterator,
    class TDistanceIteratorType = typename std::vector<double>::iterator,
    class TDistanceFunction = Kratos::SearchUtils::SquaredDistanceFunction<TDimension,TPointType>
    >
class BinsDynamic
    : public TreeNode<TDimension,TPointType, TPointerType, TIteratorType, TDistanceIteratorType, typename std::vector<TPointerType>::iterator >
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

    /// Pointer definition of BinsDynamic
    KRATOS_CLASS_POINTER_DEFINITION(BinsDynamic);

    enum { Dimension = TDimension };

    using PointType = TPointType;
    using BoundingBoxType = BoundingBox<PointType>;
    using ObjectType = typename GetObjectType<PointType>::type;
    using ContainerType = TContainerType;
    using IteratorType = TIteratorType;
    using DistanceIteratorType = TDistanceIteratorType;
    using PointerType = TPointerType;
    using DistanceFunction = TDistanceFunction;

    using TreeNodeType = TreeNode<Dimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType>;

    using CoordinateType = typename TreeNodeType::CoordinateType; // double
    using SizeType = typename TreeNodeType::SizeType;             // std::size_t
    using IndexType = typename TreeNodeType::IndexType;           // std::size_t

    using CoordinateArray = Tvector<CoordinateType, Dimension>;
    using SizeArray = Tvector<SizeType, Dimension>;
    using IndexArray = Tvector<IndexType, Dimension>;

    using IteratorIteratorType = typename TreeNodeType::IteratorIteratorType;
    using SearchStructureType = typename TreeNodeType::SearchStructureType;

    using LocalContainerType = std::vector<PointerType>;
    using LocalIterator = typename LocalContainerType::iterator;

    using CellType = Tvector<IndexType, Dimension>;
    using CellContainerType = std::vector<LocalContainerType>;

    using SearchNearestInRange = Kratos::SearchUtils::SearchNearestInRange<PointType, PointerType, LocalIterator, DistanceFunction, CoordinateType>;
    using SearchRadiusInRange = Kratos::SearchUtils::SearchRadiusInRange<PointType, LocalIterator, DistanceIteratorType, DistanceFunction, SizeType, CoordinateType, IteratorType>;
    using SearchBoxInRange = Kratos::SearchUtils::SearchBoxInRange<PointType, LocalIterator, SizeType, Dimension, IteratorType>;

    using CoordinateVectorType = std::vector<CoordinateType>;
    using IteratorVectorType = std::vector<IteratorType>;
    using DistanceIteratorVectorType = std::vector<DistanceIteratorType>;

    // Legacy typedef (to preserve compatibility in case someone was using these definitions)
    using PointVector = LocalContainerType;
    using PointIterator = LocalIterator;
    using LeafType = TreeNodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor for BinsDynamic.
     * @details This constructor initializes an empty BinsDynamic object with null iterators
     * and sets the number of cells to 0.
     */
    BinsDynamic()
        : mitPointsBegin(this->NullIterator()), mitPointsEnd(this->NullIterator()), mNumCells(0)
    {};

    /**
     * @brief Constructor for BinsDynamic with iterators and optional bucket size.
     * @param PointBegin   Iterator pointing to the beginning of the point data.
     * @param PointEnd     Iterator pointing to the end of the point data.
     * @param BucketSize   Optional parameter to specify the bucket size (default is 1).
     */
    BinsDynamic(IteratorType const& PointBegin, IteratorType const& PointEnd, SizeType BucketSize = 1)
        : mitPointsBegin(PointBegin), mitPointsEnd(PointEnd)
    {
        if (mitPointsBegin == mitPointsEnd)
            return;
        mNumCells = std::distance(mitPointsBegin, mitPointsEnd);
        CalculateBoundingBox();
        CalculateCellSize(mNumCells);
        AllocateCellsContainer();
        GenerateBins();
    }

    /**
     * @brief Constructor for BinsDynamic with iterators, min point, max point, and optional bucket size.
     * @param PointBegin   Iterator pointing to the beginning of the point data.
     * @param PointEnd     Iterator pointing to the end of the point data.
     * @param MinPoint     The minimum point of the bounding box.
     * @param MaxPoint     The maximum point of the bounding box.
     * @param BucketSize   Optional parameter to specify the bucket size (default is 1).
     */
    BinsDynamic(IteratorType const& PointBegin, IteratorType const& PointEnd, PointType const& MinPoint, PointType const& MaxPoint, SizeType BucketSize = 1)
        : mitPointsBegin(PointBegin), mitPointsEnd(PointEnd), mBoundingBox(MinPoint, MaxPoint)
    {
        if (mitPointsBegin == mitPointsEnd)
            return;

        mNumCells = std::distance(mitPointsBegin, mitPointsEnd);
        CalculateCellSize(mNumCells);
        AllocateCellsContainer();
        GenerateBins();
    }

    /**
     * @brief Constructor for BinsDynamic with min point, max point, and optional bucket size.
     * @param MinPoint     The minimum point of the bounding box.
     * @param MaxPoint     The maximum point of the bounding box.
     * @param BucketSize   Optional parameter to specify the bucket size (default is 1).
     */
    BinsDynamic(PointType const& MinPoint, PointType const& MaxPoint, SizeType BucketSize)
        : mBoundingBox(MinPoint, MaxPoint), mNumCells(0)
    {
        AssignCellSize(BucketSize);
        AllocateCellsContainer();
    }

    /**
     * @brief Constructor for BinsDynamic with iterators, box size, and optional bucket size.
     * @param PointBegin   Iterator pointing to the beginning of the point data.
     * @param PointEnd     Iterator pointing to the end of the point data.
     * @param BoxSize      The size of each cell in the bounding box.
     * @param BucketSize   Optional parameter to specify the bucket size (default is 1).
     */
    BinsDynamic(IteratorType const& PointBegin, IteratorType const& PointEnd, CoordinateType BoxSize, SizeType BucketSize = 1)
        : mitPointsBegin(PointBegin), mitPointsEnd(PointEnd)
    {
        if (mitPointsBegin == mitPointsEnd)
            return;
        mNumCells = std::distance(mitPointsBegin, mitPointsEnd);
        CalculateBoundingBox();
        AssignCellSize(BoxSize);
        AllocateCellsContainer();
        GenerateBins();
    }

    // Destructor
    ~BinsDynamic() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get an iterator pointing to the beginning of the point data.
     * @return An iterator pointing to the beginning of the point data.
     */
    IteratorType Begin()
    {
        return mitPointsBegin;
    }

    /**
     * @brief Get an iterator pointing to the end of the point data.
     * @return An iterator pointing to the end of the point data.
     */
    IteratorType End()
    {
        return mitPointsBegin;
    }

    /**
     * @brief Get the size of a cell in a specific dimension.
     * @param iDim The dimension for which to retrieve the cell size.
     * @return The size of a cell in the specified dimension.
     * @deprecated This function is deprecated. Use GetCellSize() instead.
     */
    KRATOS_DEPRECATED CoordinateType CellSize(SizeType const& iDim)
    {
        return mCellSize[iDim];
    }

    /**
     * @brief Get the number of cells in a specific dimension.
     * @param iDim The dimension for which to retrieve the number of cells.
     * @return The number of cells in the specified dimension.
     * @deprecated This function is deprecated. Use GetNumCells() instead.
     */
    KRATOS_DEPRECATED SizeType NumCell(SizeType const& iDim)
    {
        return mN[iDim];
    }

    /**
     * @brief Get the Cell Container object
     * @return CellContainerType& The Cell Container object
     */
    CellContainerType& GetCellContainer()
    {
        return mCells;
    }

    /**
     * @brief Get the Divisions object
     * @return SizeArray& Array containing the number of Cells in each dimension
     */
    SizeArray& GetDivisions()
    {
        return mN;
    }

    /**
     * @brief Get the Cell Size object
     * @return CoordinateArray& Array containing the size of the Cell in each dimension
     */
    CoordinateArray& GetCellSize()
    {
        return mCellSize;
    }

    /**
     * @brief Get a reference to the low point of the bounding box.
     * @return A reference to the low point of the bounding box.
     */
    PointType& GetMinPoint()
    {
        return mBoundingBox.GetMinPoint();
    }

    /**
     * @brief Get a reference to the high point of the bounding box.
     * @return A reference to the high point of the bounding box.
     */
    PointType& GetMaxPoint()
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
    ///@name Operations
    ///@{

    /**
    * @brief Calculate the bounding box of the tree node.
    * @details This function calculates the bounding box of the tree node based on the points contained within it. The bounding box is represented by the minimum and maximum points in each dimension.
    */
    void CalculateBoundingBox()
    {
        auto& r_min_point = GetMinPoint();
        auto& r_max_point = GetMaxPoint();
        noalias(r_min_point.Coordinates()) = (**mitPointsBegin).Coordinates();
        noalias(r_max_point.Coordinates()) = (**mitPointsBegin).Coordinates();

        for (IteratorType it_point = mitPointsBegin; it_point != mitPointsEnd; it_point++) {
            const auto& r_coordinates = (**it_point).Coordinates();
            for (SizeType i = 0; i < Dimension; i++) {
                if (r_coordinates[i] < r_min_point[i]) {
                    r_min_point[i] = r_coordinates[i];
                }
                if (r_coordinates[i] > r_max_point[i]) {
                    r_max_point[i] = r_coordinates[i];
                }
            }
        }
    }

    /**
     * @brief Calculates the cell size of the bins.
     * @details Calculates the cell size of the bins using an average approximation of the objects in the bins.
     * @param ApproximatedSize Approximate number of objects that will be stored in the bins
     */
    void CalculateCellSize(const std::size_t ApproximatedSize)
    {
        std::size_t average_number_of_cells = static_cast<std::size_t>(std::pow(static_cast<double>(ApproximatedSize), 1.00 / Dimension));

        std::array<double, 3> lengths;
        double average_length = 0.0;

        for (int i = 0; i < Dimension; i++) {
            lengths[i] = GetMaxPoint()[i] - GetMinPoint()[i];
            average_length += lengths[i];
        }
        average_length *= 1.0 / 3.0;

        if (average_length < std::numeric_limits<double>::epsilon()) {
            for(int i = 0; i < Dimension; i++) {
                mN[i] = 1;
            }
            return;
        }

        for (int i = 0; i < Dimension; i++) {
             mN[i] = static_cast<std::size_t>(lengths[i] / average_length * (double)average_number_of_cells) + 1;

            if (mN[i] > 1) {
                mCellSize[i] = lengths[i] / mN[i];
            } else {
                mCellSize[i] = average_length;
            }

            mInvCellSize[i] = 1.00 / mCellSize[i];
        }
    }

    /**
     * @brief Assign a uniform cell size for all dimensions.
     * @details This function sets the cell size for all dimensions to the same value.
     * @param BoxSize The size of each cell in all dimensions.
     */
    void AssignCellSize(CoordinateType BoxSize)
    {
        for (SizeType i = 0; i < Dimension; i++) {
            mCellSize[i] = BoxSize;
            mInvCellSize[i] = 1.0 / mCellSize[i];
            mN[i] = static_cast<SizeType>((GetMaxPoint()[i] - GetMinPoint()[i]) / mCellSize[i]) + 1;
        }
    }

    /**
     * @brief Allocate memory for the cell container.
     * @details This function allocates memory for the cell container based on the computed number of cells in each dimension.
     */
    void AllocateCellsContainer()
    {
        SizeType Size = 1;
        for (SizeType i = 0; i < Dimension; i++)
            Size *= mN[i];
        // Resize Global Container
        mCells.resize(Size);
    }

    /**
     * @brief Generate bins by assigning points to cells.
     * @details This function assigns points to cells based on their positions and populates the cell container with the points.
     */
    void GenerateBins()
    {
        for (IteratorType i_point = mitPointsBegin; i_point != mitPointsEnd; i_point++)
            mCells[CalculateIndex(**i_point)].push_back(*i_point);
    }

    /**
     * @brief Calculate the position of a coordinate in a specific dimension.
     * @param ThisCoord      The coordinate whose position is to be calculated.
     * @param ThisDimension  The dimension in which the coordinate is located.
     * @return The position index of the coordinate in the specified dimension.
     */
    IndexType CalculatePosition(CoordinateType const& ThisCoord, SizeType ThisDimension)
    {
        CoordinateType d_index = (ThisCoord - GetMinPoint()[ThisDimension]) * mInvCellSize[ThisDimension];
        IndexType index = static_cast<IndexType>((d_index < 0.00) ? 0.00 : d_index);
        return (index > mN[ThisDimension] - 1) ? mN[ThisDimension] - 1 : index;
    }

    /**
     * @brief Calculate the index of a point within the cell container.
     * @param ThisPoint The point whose index is to be calculated.
     * @return The index of the point within the cell container.
     */
    IndexType CalculateIndex(PointType const& ThisPoint)
    {
        IndexType Index = 0;
        for (SizeType iDim = Dimension - 1; iDim > 0; iDim--) {
            Index += CalculatePosition(ThisPoint[iDim], iDim);
            Index *= mN[iDim - 1];
        }
        Index += CalculatePosition(ThisPoint[0], 0);
        return Index;
    }

    /**
     * @brief Calculate the index of a cell within the cell container.
     * @param ThisIndex The index of the cell in each dimension.
     * @return The index of the cell within the cell container.
     */
    IndexType CalculateIndex(CellType const& ThisIndex)
    {
        IndexType Index = 0;
        for (SizeType iDim = Dimension - 1; iDim > 0; iDim--) {
            Index += ThisIndex[iDim];
            Index *= mN[iDim - 1];
        }
        Index += ThisIndex[0];
        return Index;
    }

    /**
     * @brief Calculate the cell index of a point.
     * @details This function calculates the cell index for a given point in each dimension.
     * @param ThisPoint The point whose cell index is to be calculated.
     * @return The cell index of the point in each dimension.
     */
    CellType CalculateCell(PointType const& ThisPoint)
    {
        CellType Cell;
        for (SizeType i = 0; i < Dimension; i++)
            Cell[i] = CalculatePosition(ThisPoint[i], i);
        return Cell;
    }

    /**
     * @brief Calculate the cell index of a point with an added radius.
     * @details This function calculates the cell index for a given point with an added radius in each dimension.
     * @param ThisPoint The point whose cell index is to be calculated.
     * @param Radius    The radius added to each coordinate before calculating the cell index.
     * @return The cell index of the point with the added radius in each dimension.
     */
    CellType CalculateCell(PointType const& ThisPoint, CoordinateType Radius)
    {
        CellType Cell;
        for (SizeType i = 0; i < Dimension; i++)
            Cell[i] = CalculatePosition(ThisPoint[i] + Radius, i);
        return Cell;
    }

    /**
     * @brief Add a point to the appropriate cell.
     * @details This function adds a point to the cell corresponding to its position.
     * @param ThisPoint The point to be added to the cell.
     */
    void AddPoint(PointerType const& ThisPoint)
    {
        mCells[CalculateIndex(*ThisPoint)].push_back(ThisPoint);
        mNumCells++;
    }

    //************************************************************************

    PointerType ExistPoint( PointerType const& ThisPoint, CoordinateType const Tolerance = static_cast<CoordinateType>(10.0*DBL_EPSILON) )
    {
        PointerType Nearest;
        CoordinateType Distance = static_cast<CoordinateType>(DBL_MAX);
        bool Found;
        SearchStructureType Box( CalculateCell(*ThisPoint,-Tolerance), CalculateCell(*ThisPoint,Tolerance), mN );
        SearchNearestInBox( *ThisPoint, Nearest, Distance, Box, Found );
        if(Found)
            return Nearest;
        return this->NullPointer();
    }

    //************************************************************************

    PointerType SearchNearestPoint( PointType const& ThisPoint )
    {
        if( mitPointsBegin == mitPointsEnd )
            return this->NullPointer();

        PointerType Result            = *mitPointsBegin;
        CoordinateType ResultDistance = static_cast<CoordinateType>(DBL_MAX);
        SearchStructureType Box( CalculateCell(ThisPoint), mN );
        SearchNearestPointLocal( ThisPoint, Result, ResultDistance, Box );
        return Result;
    }

    //************************************************************************

    PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType& ResultDistance )
    {
        if( mitPointsBegin == mitPointsEnd )
            return this->NullPointer();

        PointerType Result = *mitPointsBegin;
        ResultDistance     = static_cast<CoordinateType>(DBL_MAX);
        SearchStructureType Box( CalculateCell(ThisPoint), mN );
        SearchNearestPointLocal( ThisPoint, Result, ResultDistance, Box);
        return Result;
    }

    //************************************************************************

    // New Thread Safe!!!
    PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType& rResultDistance, SearchStructureType& Box )
    {
        PointerType Result = *mitPointsBegin; //static_cast<PointerType>(NULL);
        rResultDistance    = static_cast<CoordinateType>(DBL_MAX);
        Box.Set( CalculateCell(ThisPoint), mN );
        SearchNearestPointLocal( ThisPoint, Result, rResultDistance, Box);
        return Result;
    }

    //************************************************************************

    void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance ) override
    {
        SearchStructureType Box;
        Box.Set( CalculateCell(ThisPoint), mN );
        SearchNearestPointLocal(ThisPoint,rResult,rResultDistance,Box);
    }

    //************************************************************************

    void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box ) override
    {
        // This case is when BinStatic is a LeafType in Other Spacial Structure
        // Then, it is possible a better Result before this search
        Box.Set( CalculateCell(ThisPoint), mN );
        SearchNearestPointLocal( ThisPoint, rResult, rResultDistance, Box );
    }

    //************************************************************************

    void SearchNearestPoint( PointType* const& ThisPoints, SizeType const& NumberOfPoints, IteratorType &Results, std::vector<CoordinateType> ResultsDistances)
    {
        IndexPartition<SizeType>(NumberOfPoints).for_each(
            [&](SizeType iPoint)
            { Results[iPoint] = SearchNearestPoint(ThisPoints[iPoint],ResultsDistances[iPoint]); }
        );
    }

    //************************************************************************

    void SearchNearestPointLocal( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box )
    {
        if( mitPointsBegin == mitPointsEnd )
            return;

        bool Found = false;

        // set mBox
        Box.Set( CalculateCell(ThisPoint), mN );

        // initial search
        ++Box;
        SearchNearestInBox( ThisPoint, rResult, rResultDistance, Box, Found );
        // increase mBox and try again
        while(!Found) {
            ++Box;
            SearchNearestInBox( ThisPoint, rResult, rResultDistance, Box, Found );
        }

    }

    //************************************************************************

    SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results,
                             DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults )
    {
        CoordinateType Radius2 = Radius * Radius;
        SizeType NumberOfResults = 0;
        SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }

    //************************************************************************

    SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results,
                             DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
    {
        CoordinateType Radius2 = Radius * Radius;
        SizeType NumberOfResults = 0;
        Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }

    //************************************************************************

    void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults ) override
    {
        SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
    }

    //************************************************************************

    void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box ) override
    {
        Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
    }

    //************************************************************************

    void SearchInRadius( PointerType const& ThisPoints, SizeType const& NumberOfPoints, CoordinateVectorType const& Radius, IteratorVectorType Results,
                         DistanceIteratorVectorType ResultsDistances, std::vector<SizeType>& NumberOfResults, SizeType const& MaxNumberOfResults )
    {
        IndexPartition<SizeType>(NumberOfPoints).for_each(
            [&](SizeType iPoint)
            { NumberOfResults[iPoint] = SearchInRadius(ThisPoints[iPoint],Radius[iPoint],Results[iPoint],ResultsDistances[iPoint],MaxNumberOfResults); }
        );
    }

    // **** THREAD SAFE

    // Dimension = 1
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
            SearchRadiusInRange()(mCells[I].begin(),mCells[I].end(),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                SearchRadiusInRange()(mCells[I].begin(),mCells[I].end(),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 3
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                    SearchRadiusInRange()(mCells[I].begin(),mCells[I].end(),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
    }

    //************************************************************************

    SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, IteratorType Results, SizeType MaxNumberOfResults )
    {
        CoordinateType Radius2 = Radius * Radius;
        SizeType NumberOfResults = 0;
        SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }

    //************************************************************************

    SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, IteratorType Results,
                             SizeType MaxNumberOfResults, SearchStructureType& Box )
    {
        CoordinateType Radius2 = Radius * Radius;
        SizeType NumberOfResults = 0;
        Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }

    //************************************************************************

    void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults ) override
    {
        SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
    }

    //************************************************************************

    void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box ) override
    {
        Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
    }

    //************************************************************************

    // **** THREAD SAFE

    // Dimension = 1
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I++ )
            SearchRadiusInRange()(mCells[I].begin(),mCells[I].end(),ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I++ )
                SearchRadiusInRange()(mCells[I].begin(),mCells[I].end(),ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 3
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I++ )
                    SearchRadiusInRange()(mCells[I].begin(),mCells[I].end(),ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
    }

    //************************************************************************
    //************************************************************************

    // Dimension = 1
    void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance,
                             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, bool& Found )
    {
        Found = false;
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
            SearchNearestInRange()( mCells[I].begin(), mCells[I].end(), ThisPoint, ResultPoint, ResultDistance, Found );
    }

    // Dimension = 2
    void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance,
                             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, bool& Found )
    {
        Found = false;
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                SearchNearestInRange()( mCells[I].begin(), mCells[I].end(), ThisPoint, ResultPoint, ResultDistance, Found );
    }

    // Dimension = 3
    void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance,
                             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, bool& Found )
    {
        Found = false;
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                    SearchNearestInRange()( mCells[I].begin(), mCells[I].end(), ThisPoint, ResultPoint, ResultDistance, Found );
    }

    //************************************************************************
    //************************************************************************

    SizeType SearchInBox( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType Results,
                          SizeType MaxNumberOfResults )
    {
        SizeType NumberOfResults = 0;
        SearchStructureType Box( CalculateCell(SearchMinPoint), CalculateCell(SearchMaxPoint), mN );
        SearchInBoxLocal( SearchMinPoint, SearchMaxPoint, Results, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }

    //************************************************************************

    void SearchInBox(PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& Results, SizeType& NumberOfResults,
                     SizeType const& MaxNumberOfResults ) override
    {
        NumberOfResults = 0;
        SearchStructureType Box( CalculateCell(SearchMinPoint), CalculateCell(SearchMaxPoint), mN );
        SearchInBoxLocal( SearchMinPoint, SearchMaxPoint, Results, NumberOfResults, MaxNumberOfResults, Box );
    }

    //************************************************************************

    // Dimension = 1
    void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
                           SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                           SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
            SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,mCells[I].begin(),mCells[I].end(),ResultsPoint,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
                           SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                           SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,mCells[I].begin(),mCells[I].end(),ResultsPoint,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 3
    void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
                           SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                           SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                    SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,mCells[I].begin(),mCells[I].end(),ResultsPoint,NumberOfResults,MaxNumberOfResults);
    }

    //************************************************************************

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "BinsDynamic";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BinsDynamic";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const override
    {
        rOStream << Perfix << "Bin[" << SearchUtils::PointerDistance(mitPointsBegin, mitPointsEnd) << "] : " << std::endl;
        for(typename CellContainerType::const_iterator i_cell = mCells.begin() ; i_cell != mCells.end() ; i_cell++)
        {
            rOStream << Perfix << "[ " ;
            for(typename LocalContainerType::const_iterator i_point = i_cell->begin() ; i_point != i_cell->end() ; i_point++)
                rOStream << **i_point << "    ";
            rOStream << " ]" << std::endl;
        }
        rOStream << std::endl;
    }

    /// Print Size of Container
    void PrintSize( std::ostream& rout )
    {
        rout << " BinsSize: ";
        for(SizeType i = 0 ; i < Dimension ; i++)
            rout << "[" << mN[i] << "]";
        rout << std::endl;
    }

    /// Print Limits Points of the Container
    void PrintBox( std::ostream& rout )
    {
        rout << " BinsBox: Min [";
        GetMinPoint().Print(rout);
        rout <<       "];  Max [";
        GetMaxPoint().Print(rout);
        rout <<       "];  Size [";
        mCellSize.Print(rout);
        rout << "]" << std::endl;
    }

    /// Assignment operator.
    BinsDynamic& operator=(BinsDynamic const& rOther);

    /// Copy constructor.
    BinsDynamic(BinsDynamic const& rOther);

private:
    ///@name Member Variables
    ///@{

    IteratorType mitPointsBegin;  /// Iterator to the first point
    IteratorType mitPointsEnd;    /// Iterator to the last point

    BoundingBoxType mBoundingBox; /// The bounding box of the tree

    CoordinateArray mCellSize;    /// Array representing the size of each cell in each dimension.

    CoordinateArray mInvCellSize; /// Array representing the inverse of the cell size in each dimension.

    SizeArray mN;                 /// Array representing the number of cells in each dimension.

    SizeType mNumCells;           /// The total number of cells in the data structure.

    CellContainerType mCells;     /// Bins Access Vector ( vector<Iterator> )

    // Work Variables ( For non-copy of Search Variables )
    //BinBox SearchBox;

    ///@}
public:
    static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, PointType MaxPoint, PointType MinPoint, SizeType BucketSize)
    {
        const SizeType number_of_points = SearchUtils::PointerDistance(PointsBegin,PointsEnd);
        if (number_of_points == 0) {
            return NULL;
        } else {
            return new BinsDynamic( PointsBegin, PointsEnd, MinPoint, MaxPoint, BucketSize );
        }

    }

};

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template<
    std::size_t TDimension,
    class TPointType,
    class TContainerType,
    class TPointerType,
    class TIteratorType,
    class TDistanceIteratorType,
    class TDistanceFunction >
std::ostream & operator<<( std::ostream& rOStream,
                           BinsDynamic<TDimension,TPointType,TContainerType,TPointerType,TIteratorType,TDistanceIteratorType,TDistanceFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintSize(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}
///@}

}  // namespace Kratos.
