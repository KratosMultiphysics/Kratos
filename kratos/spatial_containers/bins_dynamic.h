//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Carlos Roig
//


#if !defined(KRATOS_BINS_DYNAMIC_CONTAINER_H_INCLUDE)
#define KRATOS_BINS_DYNAMIC_CONTAINER_H_INCLUDE

#include <array>
#include <cmath>
#include <algorithm>

#include "tree.h"

namespace Kratos
{

template<
std::size_t TDimension,
    class TPointType,
    class TContainerType,
    class TPointerType = typename TContainerType::value_type,
    class TIteratorType = typename TContainerType::iterator,
    class TDistanceIteratorType = typename std::vector<double>::iterator,
    class TDistanceFunction = Kratos::SearchUtils::SquaredDistanceFunction<TDimension,TPointType>
    >
class BinsDynamic : public TreeNode<TDimension,TPointType, TPointerType, TIteratorType, TDistanceIteratorType, typename std::vector<TPointerType>::iterator >
{

public:

    /// Pointer definition of BinsDynamic
    KRATOS_CLASS_POINTER_DEFINITION(BinsDynamic);

    typedef TreeNode<TDimension,TPointType,TPointerType,TIteratorType,TDistanceIteratorType> TreeNodeType;
    typedef TPointType                         PointType;
    typedef TContainerType                     ContainerType;
    typedef TIteratorType                      IteratorType;
    typedef TDistanceIteratorType              DistanceIteratorType;
    typedef TPointerType                       PointerType;
    typedef TDistanceFunction                  DistanceFunction;
    enum { Dimension = TDimension };

    typedef typename TreeNodeType::CoordinateType  CoordinateType;  // double
    typedef typename TreeNodeType::SizeType        SizeType;        // std::size_t
    typedef typename TreeNodeType::IndexType       IndexType;       // std::size_t

    typedef TreeNodeType LeafType;

    typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
    typedef typename TreeNodeType::SearchStructureType SearchStructureType;


    // Local Container ( PointPointer Container per Cell )
    // can be different to ContainerType
    // not always PointVector == ContainerType ( if ContainerType = C array )
    typedef std::vector<PointerType>       PointVector;
    typedef typename PointVector::iterator PointIterator;

    // Global Container
    typedef std::vector<PointVector>        CellContainerType;
    //typedef typename CellContainerType::iterator IteratorIteratorType;

    typedef Tvector<IndexType,TDimension>   CellType;

    typedef Kratos::SearchUtils::SearchNearestInRange<PointType,PointerType,PointIterator,DistanceFunction,CoordinateType> SearchNearestInRange;
    typedef Kratos::SearchUtils::SearchRadiusInRange<PointType,PointIterator,DistanceIteratorType,DistanceFunction,SizeType,CoordinateType,IteratorType> SearchRadiusInRange;
    typedef Kratos::SearchUtils::SearchBoxInRange<PointType,PointIterator,SizeType,TDimension,IteratorType> SearchBoxInRange;

    typedef std::vector<CoordinateType>         CoordinateVectorType;
    typedef std::vector<IteratorType>           IteratorVectorType;
    typedef std::vector<DistanceIteratorType>   DistanceIteratorVectorType;
    
    
public:

    //************************************************************************

    // constructor 1
    BinsDynamic() : mPointBegin(this->NullIterator()), mPointEnd(this->NullIterator()), mNumPoints(0)
    {};

    //************************************************************************

    BinsDynamic( IteratorType const& PointBegin, IteratorType const& PointEnd, SizeType BucketSize = 1 )
        : mPointBegin(PointBegin), mPointEnd(PointEnd)
    {
        if(mPointBegin==mPointEnd)
            return;
        mNumPoints = std::distance(mPointBegin,mPointEnd);
        CalculateBoundingBox();
        CalculateCellSize(mNumPoints);
        AllocateCellsContainer();
        GenerateBins();
    }

    //************************************************************************

    BinsDynamic( IteratorType const& PointBegin, IteratorType const& PointEnd, PointType const& MinPoint, PointType const& MaxPoint, SizeType BucketSize = 1 )
        : mPointBegin(PointBegin), mPointEnd(PointEnd)
    {
        if(mPointBegin==mPointEnd)
            return;

        mNumPoints = std::distance(mPointBegin,mPointEnd);
        for(SizeType i = 0 ; i < TDimension ; i++)
        {
            mMinPoint[i] = MinPoint[i];
            mMaxPoint[i] = MaxPoint[i];
        }
        CalculateCellSize(mNumPoints);
        AllocateCellsContainer();
        GenerateBins();
    }

    //************************************************************************

    BinsDynamic( PointType const& MinPoint, PointType const& MaxPoint, SizeType BucketSize )
        : mNumPoints(0)
    {
        for(SizeType i = 0 ; i < TDimension ; i++)
        {
            mMinPoint[i] = MinPoint[i];
            mMaxPoint[i] = MaxPoint[i];
        }
        AssingnCellSize(BucketSize);
        AllocateCellsContainer();
    }

    //************************************************************************

    BinsDynamic( IteratorType const& PointBegin, IteratorType const& PointEnd, CoordinateType BoxSize, SizeType BucketSize = 1 )
        : mPointBegin(PointBegin), mPointEnd(PointEnd)
    {
        if(mPointBegin==mPointEnd)
            return;
        mNumPoints = std::distance(mPointBegin,mPointEnd);
        CalculateBoundingBox();
        AssingnCellSize(BoxSize);
        AllocateCellsContainer();
        GenerateBins();
    }

    //************************************************************************

    // destructor
    virtual ~BinsDynamic() { }

    //************************************************************************

    IteratorType Begin()
    {
        return mPointBegin;
    }

    //************************************************************************

    IteratorType End()
    {
        return mPointBegin;
    }

    //************************************************************************

    CoordinateType CellSize( SizeType const& iDim )
    {
        return mCellSize[iDim];
    }

    //************************************************************************

    SizeType NumCell( SizeType const& iDim )
    {
        return mN[iDim];
    }

    //************************************************************************

    void CalculateBoundingBox()
    {
        for(SizeType i = 0 ; i < TDimension ; i++)
        {
            mMinPoint[i] = (**mPointBegin)[i];
            mMaxPoint[i] = (**mPointBegin)[i];
        }
        for(IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
            for(SizeType i = 0 ; i < TDimension ; i++)
            {
                if( (**Point)[i] < mMinPoint[i] ) mMinPoint[i] = (**Point)[i];
                if( (**Point)[i] > mMaxPoint[i] ) mMaxPoint[i] = (**Point)[i];
            }
    }

    //************************************************************************

    /** 
     * @brief Calculates the cell size of the bins.
     * 
     * Calculates the cell size of the bins using an average aproximation of the objects in the bins.
     * 
     * @param ApproximatedSize Aproximate number of objects that will be stored in the bins
     */
    void CalculateCellSize(std::size_t ApproximatedSize) 
    {
        std::size_t average_number_of_cells = static_cast<std::size_t>(std::pow(static_cast<double>(ApproximatedSize), 1.00 / Dimension));
        
        std::array<double, 3> lengths;
        double average_length = 0.00;
        
        for (int i = 0; i < Dimension; i++) {
            lengths[i] = mMaxPoint[i] - mMinPoint[i];
            average_length += lengths[i];
        }
        average_length *= 1.00 / 3.00;

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

    //************************************************************************

    void AssignCellSize( CoordinateType BoxSize )
    {
        for(SizeType i = 0 ; i < TDimension ; i++)
        {
            mCellSize[i] = BoxSize;
            mInvCellSize[i] = 1.00 / mCellSize[i];
            mN[i] = static_cast<SizeType>( (mMaxPoint[i]-mMinPoint[i]) / mCellSize[i]) + 1;
        }
    }

    //************************************************************************

    void AllocateCellsContainer()
    {
        SizeType Size = 1;
        for(SizeType i = 0 ; i < TDimension ; i++)
            Size *= mN[i];
        // Resize Global Container
        mPoints.resize(Size);
    }

    //************************************************************************

    void GenerateBins()
    {

        for(IteratorType i_point = mPointBegin ; i_point != mPointEnd ; i_point++)
            mPoints[CalculateIndex(**i_point)].push_back(*i_point);

    }

    //************************************************************************

    IndexType CalculatePosition( CoordinateType const& ThisCoord, SizeType ThisDimension )
    {
        CoordinateType d_index = (ThisCoord - mMinPoint[ThisDimension]) * mInvCellSize[ThisDimension];
        IndexType index = static_cast<IndexType>( (d_index < 0.00) ? 0.00 : d_index );
        return  (index > mN[ThisDimension]-1) ? mN[ThisDimension]-1 : index;
    }

    //************************************************************************

    IndexType CalculateIndex( PointType const& ThisPoint )
    {
        IndexType Index = 0;
        for(SizeType iDim = TDimension-1 ; iDim > 0 ; iDim--)
        {
            Index += CalculatePosition(ThisPoint[iDim],iDim);
            Index *= mN[iDim-1];
        }
        Index += CalculatePosition(ThisPoint[0],0);
        return Index;
    }

    //************************************************************************

    IndexType CalculateIndex( CellType const& ThisIndex )
    {
        IndexType Index = 0;
        for(SizeType iDim = TDimension-1 ; iDim > 0 ; iDim--)
        {
            Index += ThisIndex[iDim];
            Index *= mN[iDim-1];
        }
        Index += ThisIndex[0];
        return Index;
    }

    //************************************************************************

    CellType CalculateCell( PointType const& ThisPoint )
    {
        CellType Cell;
        for(SizeType i = 0 ; i < TDimension ; i++)
            Cell[i] = CalculatePosition(ThisPoint[i],i);
        return Cell;
    }

    CellType CalculateCell( PointType const& ThisPoint, CoordinateType Radius )
    {
        CellType Cell;
        for(SizeType i = 0 ; i < TDimension ; i++)
            Cell[i] = CalculatePosition(ThisPoint[i]+Radius,i);
        return Cell;
    }

    //************************************************************************

    void AddPoint( PointerType const& ThisPoint )
    {
        mPoints[CalculateIndex(*ThisPoint)].push_back(ThisPoint);
        mNumPoints++;
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
        if( mPointBegin == mPointEnd )
            return this->NullPointer();

        PointerType Result            = *mPointBegin;
        CoordinateType ResultDistance = static_cast<CoordinateType>(DBL_MAX);
        SearchStructureType Box( CalculateCell(ThisPoint), mN );
        SearchNearestPointLocal( ThisPoint, Result, ResultDistance, Box );
        return Result;
    }

    //************************************************************************

    PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType ResultDistance )
    {
        if( mPointBegin == mPointEnd )
            return this->NullPointer();

        PointerType Result = *mPointBegin;
        ResultDistance     = static_cast<CoordinateType>(DBL_MAX);
        SearchStructureType Box( CalculateCell(ThisPoint), mN );
        SearchNearestPointLocal( ThisPoint, Result, ResultDistance, Box);
        return Result;
    }

    //************************************************************************

    // New Thread Safe!!!
    PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType& rResultDistance, SearchStructureType& Box )
    {
        PointerType Result = *mPointBegin; //static_cast<PointerType>(NULL);
        rResultDistance    = static_cast<CoordinateType>(DBL_MAX);
        Box.Set( CalculateCell(ThisPoint), mN );
        SearchNearestPointLocal( ThisPoint, Result, rResultDistance, Box);
        return Result;
    }

    //************************************************************************

    void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance )
    {
        SearchStructureType Box;
        Box.Set( CalculateCell(ThisPoint), mN );
        SearchNearestPointLocal(ThisPoint,rResult,rResultDistance,Box);
    }

    //************************************************************************

    void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box )
    {
        // This case is when BinStatic is a LeafType in Other Spacial Structure
        // Then, it is possible a better Result before this search
        Box.Set( CalculateCell(ThisPoint), mN );
        SearchNearestPointLocal( ThisPoint, rResult, rResultDistance, Box );
    }
    
    //************************************************************************
        
    void SearchNearestPoint( PointType* const& ThisPoints, SizeType const& NumberOfPoints, IteratorType &Results, std::vector<CoordinateType> ResultsDistances)
    {
        #pragma omp parallel for
        for(int k=0; k< NumberOfPoints; k++)
            Results[k] = SearchNearestPoint(ThisPoints[k],ResultsDistances[k]);
    }

    //************************************************************************

    void SearchNearestPointLocal( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box )
    {
        if( mPointBegin == mPointEnd )
            return;

        bool Found = false;

        // set mBox
        Box.Set( CalculateCell(ThisPoint), mN );

        // initial search
        ++Box;
        SearchNearestInBox( ThisPoint, rResult, rResultDistance, Box, Found );
        // increase mBox and try again
        while(!Found)
        {
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
                         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
    {
        SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
    }

    //************************************************************************

    void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
    {
        Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
    }
    
    //************************************************************************

    void SearchInRadius( PointerType const& ThisPoints, SizeType const& NumberOfPoints, CoordinateVectorType const& Radius, IteratorVectorType Results,
                         DistanceIteratorVectorType ResultsDistances, std::vector<SizeType>& NumberOfResults, SizeType const& MaxNumberOfResults )
    {
        #pragma omp parallel for
        for(int k=0; k< NumberOfPoints; k++)
            NumberOfResults[k] = SearchInRadius(ThisPoints[k],Radius[k],Results[k],ResultsDistances[k],MaxNumberOfResults);
    }

    // **** THREAD SAFE

    // Dimension = 1
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
            SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 3
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                    SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
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
                         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
    {
        SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
    }

    //************************************************************************

    void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
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
            SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I++ )
                SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 3
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I++ )
                    SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
    }

    //************************************************************************
    //************************************************************************

    // Dimension = 1
    void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance,
                             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, bool& Found )
    {
        Found = false;
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
            SearchNearestInRange()( mPoints[I].begin(), mPoints[I].end(), ThisPoint, ResultPoint, ResultDistance, Found );
    }

    // Dimension = 2
    void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance,
                             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, bool& Found )
    {
        Found = false;
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                SearchNearestInRange()( mPoints[I].begin(), mPoints[I].end(), ThisPoint, ResultPoint, ResultDistance, Found );
    }

    // Dimension = 3
    void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance,
                             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, bool& Found )
    {
        Found = false;
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                    SearchNearestInRange()( mPoints[I].begin(), mPoints[I].end(), ThisPoint, ResultPoint, ResultDistance, Found );
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
                     SizeType const& MaxNumberOfResults )
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
            SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,mPoints[I].begin(),mPoints[I].end(),ResultsPoint,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
                           SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                           SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,mPoints[I].begin(),mPoints[I].end(),ResultsPoint,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 3
    void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
                           SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                           SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                    SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,mPoints[I].begin(),mPoints[I].end(),ResultsPoint,NumberOfResults,MaxNumberOfResults);
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
    virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const
    {
        rOStream << Perfix << "Bin[" << SearchUtils::PointerDistance(mPointBegin, mPointEnd) << "] : " << std::endl;
        for(typename CellContainerType::const_iterator i_cell = mPoints.begin() ; i_cell != mPoints.end() ; i_cell++)
        {
            rOStream << Perfix << "[ " ;
            for(typename PointVector::const_iterator i_point = i_cell->begin() ; i_point != i_cell->end() ; i_point++)
                rOStream << **i_point << "    ";
            rOStream << " ]" << std::endl;
        }
        rOStream << std::endl;
    }

    /// Print Size of Container
    void PrintSize( std::ostream& rout )
    {
        rout << " BinsSize: ";
        for(SizeType i = 0 ; i < TDimension ; i++)
            rout << "[" << mN[i] << "]";
        rout << std::endl;
    }

    /// Print Limits Points of the Container
    void PrintBox( std::ostream& rout )
    {
        rout << " BinsBox: Min [";
        mMinPoint.Print(rout);
        rout <<       "];  Max [";
        mMaxPoint.Print(rout);
        rout <<       "];  Size [";
        mCellSize.Print(rout);
        rout << "]" << std::endl;
    }

    /**
     * @brief Get the Cell Container object
     * 
     * @return CellContainerType& The Cell Container object
     */
    CellContainerType& GetCellContainer() {
        return mCells;
    }

    /**
     * @brief Get the Divisions object
     * 
     * @return SizeArray& Array containing the number of Cells in each dimension
     */
    SizeArray& GetDivisions() {
        return mN;
    }

    /**
     * @brief Get the Cell Size object
     * 
     * @return CoordinateArray& Array containing the size of the Cell in each dimension
     */
    CoordinateArray& GetCellSize() {
        return mCellSize;
    }

    /**
     * @brief Get the Min Point object
     * 
     * @return PointType& Min point of the bins
     */
    PointType& GetMinPoint() {
        return mMinPoint;
    }

    /**
     * @brief Get the Max Point object
     * 
     * @return PointType& Max point of the bins
     */
    PointType& GetMaxPoint() {
        return mMaxPoint;
    }

    /// Assignment operator.
    BinsDynamic& operator=(BinsDynamic const& rOther);

    /// Copy constructor.
    BinsDynamic(BinsDynamic const& rOther);

private:

    IteratorType     mPointBegin;
    IteratorType     mPointEnd;

    Tvector<CoordinateType,TDimension>  mMinPoint;
    Tvector<CoordinateType,TDimension>  mMaxPoint;
    Tvector<CoordinateType,TDimension>  mCellSize;
    Tvector<CoordinateType,TDimension>  mInvCellSize;
    Tvector<SizeType,TDimension>        mN;
    SizeType                            mNumPoints;

    // Bins Access Vector ( vector<Iterator> )
    CellContainerType mPoints;

    // Work Variables ( For non-copy of Search Variables )
    //BinBox SearchBox;

public:
    static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, PointType MaxPoint, PointType MinPoint, SizeType BucketSize)
    {

        SizeType number_of_points = SearchUtils::PointerDistance(PointsBegin,PointsEnd);
        if (number_of_points == 0)
            return NULL;
        else
        {
            return new BinsDynamic( PointsBegin, PointsEnd, MinPoint, MaxPoint, BucketSize );
        }

    }

};

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



}

#endif // KRATOS_BINS_DYNAMIC_CONTAINER_H_INCLUD
