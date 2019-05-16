//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    clabra
//



#if !defined(KRATOS_BINS_CONTAINER_H_INCLUDE)
#define KRATOS_BINS_CONTAINER_H_INCLUDE

#include "tree.h"


namespace Kratos
{



template<  std::size_t TDimension,
           class TPointType,
           class TContainerType,
           class TPointerType = typename TContainerType::value_type,
           class TIteratorType = typename TContainerType::iterator,
           class TDistanceIteratorType = typename std::vector<double>::iterator,
           class TDistanceFunction = Kratos::SearchUtils::SquaredDistanceFunction<TDimension,TPointType> >
class Bins : public TreeNode<TDimension,TPointType, TPointerType, TIteratorType, TDistanceIteratorType>
{


public:

    enum { Dimension = TDimension };

    typedef TPointType                                  PointType;
    typedef TContainerType                              ContainerType;
    typedef TIteratorType                               IteratorType;
    typedef TDistanceIteratorType                       DistanceIteratorType;
    typedef TPointerType                                PointerType;
    typedef TDistanceFunction                           DistanceFunction;

    typedef TreeNode<Dimension,PointType,PointerType,IteratorType,DistanceIteratorType> TreeNodeType;

    typedef typename TreeNodeType::SizeType             SizeType;
    typedef typename TreeNodeType::IndexType            IndexType;
    typedef typename TreeNodeType::CoordinateType       CoordinateType;

    typedef Tvector<CoordinateType,Dimension>           CoordinateArray;
    typedef Tvector<SizeType,Dimension>                 SizeArray;
    typedef Tvector<IndexType,Dimension>                IndexArray;

    typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
    typedef typename TreeNodeType::SearchStructureType  SearchStructureType;

    // Local Container ( PointPointer Container per Cell )
    // can be different to ContainerType
    // not always LocalContainerType == ContainerType ( if ContainerType = C array )
    typedef std::vector<PointerType>                    LocalContainerType;
    typedef typename LocalContainerType::iterator       LocalIterator;

    typedef Tvector<IndexType,TDimension>               CellType;

    typedef std::vector<IteratorType>                   IteratorVector;
    typedef typename IteratorVector::iterator           IteratorIterator;
    typedef typename IteratorVector::const_iterator     IteratorConstIterator;

    typedef Kratos::SearchUtils::SearchNearestInRange<PointType,PointerType,IteratorType,DistanceFunction,CoordinateType> SearchNearestInRange;
    typedef Kratos::SearchUtils::SearchRadiusInRange<PointType,IteratorType,DistanceIteratorType,DistanceFunction,SizeType,CoordinateType> SearchRadiusInRange;
    typedef Kratos::SearchUtils::SearchBoxInRange<PointType,IteratorType,SizeType,TDimension> SearchBoxInRange;

    // Legacy typedef ( to preserve compativility in case someone was using this definitions)
    typedef LocalContainerType                          PointVector;
    typedef LocalIterator                               PointIterator;
    typedef TreeNodeType                                LeafType;

    /// Pointer definition of Bins
    KRATOS_CLASS_POINTER_DEFINITION(Bins);


public:

    /**
     * @brief Default Constructor
     * 
     */
    Bins() : mPointBegin(this->NullIterator()), mPointEnd(this->NullIterator()) {};


    /**
     * @brief Constructs a new BinsStatic
     * 
     * Construct a new BinsStatic using a list of points and an automatically calculate cell size.
     * 
     * @param PointBegin Iterator to the first object of the bins
     * @param PointEnd Iterator to the last object of the bins
     * @param BucketSize Unused.
     */
    Bins( IteratorType const& PointBegin, IteratorType const& PointEnd, SizeType BucketSize = 1 )
        : mPointBegin(PointBegin), mPointEnd(PointEnd)
    {
        auto NumPoints = std::distance(mPointBegin, mPointEnd);

        if(mPointBegin==mPointEnd)
            return;

        CalculateBoundingBox();
        CalculateCellSize(NumPoints);
        AllocateCellsContainer();
        GenerateBins();
    }

    /**
     * @brief Constructs a new BinsStatic
     * 
     * Construct a new BinsStatic using a list of points and an automatically calculate cell size and a custom bounding box.
     * 
     * @param PointBegin Iterator to the first object of the bins
     * @param PointEnd Iterator to the last object of the bins
     * @param MinPoint Lower point of the custom bounding box
     * @param MaxPoint Upper point of the custom bounding box
     * @param BucketSize Unused
     */
    Bins( IteratorType const& PointBegin, IteratorType const& PointEnd, PointType const& MinPoint, PointType const& MaxPoint, SizeType BucketSize = 1 )
        : mPointBegin(PointBegin), mPointEnd(PointEnd)
    {
        auto NumPoints = std::distance(mPointBegin, mPointEnd);

        if(mPointBegin==mPointEnd)
            return;

        for(SizeType i = 0 ; i < TDimension ; i++)
        {
            mMinPoint[i] = MinPoint[i];
            mMaxPoint[i] = MaxPoint[i];
        }

        CalculateCellSize(NumPoints);
        AllocateCellsContainer();
        GenerateBins();
    }

    /**
     * @brief Constructs a new BinsStatic
     * 
     * Constructs a new BinsObjectDynamic using a list of objects and an user provided cell size.
     *
     * @param PointBegin Iterator to the first object of the bins
     * @param PointEnd Iterator to the last object of the bins
     * @param cellsize Size of the cells (equal for every dimension)
     * @param BucketSize Unsued
     */
    Bins( IteratorType const& PointBegin, IteratorType const& PointEnd, CoordinateType cellsize, SizeType BucketSize = 1 )
        : mPointBegin(PointBegin), mPointEnd(PointEnd)
    {
        if(mPointBegin==mPointEnd)
            return;

        CalculateBoundingBox();
        AssignCellSize(cellsize);
        AllocateCellsContainer();
        GenerateBins();
    }

    // Destructor
    ~Bins() override { }

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

    KRATOS_DEPRECATED CoordinateType CellSize( SizeType const& iDim )
    {
        return mCellSize[iDim];
    }

    //************************************************************************

    KRATOS_DEPRECATED SizeType NumCell( SizeType const& iDim )
    {
        return mN[iDim];
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

    //************************************************************************

private:

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
        mIndexCell.resize(Size+1);
        mIndexCellBegin = mIndexCell.begin();
        mIndexCellEnd   = mIndexCell.end();
    }

    //************************************************************************

    void GenerateBins( )
    {

        LocalContainerType TempPoint(mPointBegin,mPointEnd);

        // Reset index vector
        for( IteratorIterator Iter = mIndexCell.begin(); Iter != mIndexCell.end(); Iter++)
            *Iter = mPointBegin;

        // Update storage counter, storing ahead
        for( IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
            mIndexCell[ CalculateIndex(**Point) + 1 ]++;

        // Storage/reshufing pass 1

        // Update storage counter and store
        for( IteratorIterator Iter = mIndexCell.begin()+1 ; Iter != mIndexCell.end() ; Iter++)
            *Iter = *(Iter-1) + SearchUtils::PointerDistance(mPointBegin,*Iter);

        // Point pass 2
        // Store the points in lbin1

        // Update storage counter, storing in lbin1
        for( LocalIterator Point = TempPoint.begin() ; Point != TempPoint.end() ; Point++)
            *(mIndexCell[CalculateIndex(**Point)]++) = *Point;

        // Storage/reshuffing pass 2

        // Loop over bins, in reverse order
        for(IteratorIterator Iter = mIndexCell.end()-1; Iter != mIndexCell.begin(); Iter--)
            *Iter = *(Iter-1);
        mIndexCell[0] = mPointBegin;

    }

    //************************************************************************

    IndexType CalculatePosition( CoordinateType const& ThisCoord, SizeType ThisDimension )
    {
        CoordinateType d_index = (ThisCoord - mMinPoint[ThisDimension]) * mInvCellSize[ThisDimension];
        IndexType index = static_cast<SizeType>( (d_index < 0.00) ? 0.00 : d_index );
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

public:

    //************************************************************************
    //************************************************************************

    /**
     * @brief Return the closest point to ThisPoint in case it exists or a null pointer otherwise
     * 
     * @param ThisPoint Searched Point.
     * @param Tolerance Tolerance of the search.
     * @return PointerType a pointer to the nearest point in case it exists or nullptr otherwise.
     */
    PointerType ExistPoint( PointerType const& ThisPoint, CoordinateType const Tolerance = static_cast<CoordinateType>(10.0*DBL_EPSILON) )
    {
        PointerType Nearest;
        CoordinateType Distance = static_cast<CoordinateType>(DBL_MAX);
        bool Found;
        SearchStructureType Box( CalculateCell(*ThisPoint), mN, mIndexCellBegin );
        SearchNearestInBox( *ThisPoint, Nearest, Distance, Box, Found );
        if(Found)
            return Nearest;
        return this->NullPointer();
    }

    /**
     * @brief Return the nearest point to ThisPoint. This function can not return the same point.
     * 
     * @param ThisPoint Searched Point.
     * @return PointerType Pointer to the nearest element. Cannot return the same point as the one given as input.
     */
    PointerType SearchNearestPointInner( PointerType& ThisPoint )
    {
        PointerType Result            = *mPointBegin;                           //static_cast<PointerType>(NULL);
        CoordinateType ResultDistance = static_cast<CoordinateType>(DBL_MAX);
        SearchStructureType Box( CalculateCell(*ThisPoint), mN, mIndexCellBegin );
        SearchNearestPointLocalInner( ThisPoint, Result, ResultDistance, Box );
        return Result;
    }

    /**
     * @brief Return the nearest point to ThisPoint. This function can return the same point with distance 0.
     * 
     * @param ThisPoint Searched Point.
     * @return PointerType Pointer to the nearest element. ThisPoint in case it exists inside the bins.
     */
    PointerType SearchNearestPoint( PointType const& ThisPoint )
    {
        PointerType Result            = *mPointBegin;                           //static_cast<PointerType>(NULL);
        CoordinateType ResultDistance = static_cast<CoordinateType>(DBL_MAX);
        SearchStructureType Box( CalculateCell(ThisPoint), mN, mIndexCellBegin );
        SearchNearestPointLocal( ThisPoint, Result, ResultDistance, Box );
        return Result;
    }

    //************************************************************************

    PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType& rResultDistance )
    {
        PointerType Result = *mPointBegin;                           //static_cast<PointerType>(NULL);
        rResultDistance    = static_cast<CoordinateType>(DBL_MAX);
        SearchStructureType Box( CalculateCell(ThisPoint), mN, mIndexCellBegin );
        SearchNearestPointLocal( ThisPoint, Result, rResultDistance, Box);
        return Result;
    }

    //************************************************************************

    // New Thread Safe!!!
    PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType& rResultDistance, SearchStructureType& Box )
    {
        PointerType Result            = *mPointBegin;                           //static_cast<PointerType>(NULL);
        rResultDistance = static_cast<CoordinateType>(DBL_MAX);
        Box.Set( CalculateCell(ThisPoint), mN, mIndexCellBegin );
        SearchNearestPointLocal( ThisPoint, Result, rResultDistance, Box);
        return Result;
    }

    //************************************************************************
    //************************************************************************

    void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance ) override
    {
        SearchStructureType Box( CalculateCell(ThisPoint), mN, mIndexCellBegin );
        SearchNearestPointLocal(ThisPoint,rResult,rResultDistance,Box);
    }

    //************************************************************************

    void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box ) override
    {
        // This case is when BinStatic is a LeafType in Other Spacial Structure
        // Then, it is possible a better Result before this search
        Box.Set( CalculateCell(ThisPoint), mN, mIndexCellBegin );
        SearchNearestPointLocal( ThisPoint, rResult, rResultDistance, Box );
    }
    
    //************************************************************************
    
    void SearchNearestPoint( PointerType const& ThisPoints, SizeType const& NumberOfPoints, IteratorType &Results, std::vector<CoordinateType> ResultsDistances)
    {
        #pragma omp parallel for
        for(int k=0; k< NumberOfPoints; k++)
	        Results[k] = SearchNearestPoint((&(*ThisPoints))[k],ResultsDistances[k]);
    }

    //************************************************************************
    //************************************************************************

    // **** THREAD SAFE  -> The user pass the SearchStructure (BinBox)
    void SearchNearestPointLocal( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box )
    {
        if( mPointBegin == mPointEnd )
            return;

        bool Found;

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
    //************************************************************************

    // **** THREAD SAFE  -> The user pass the SearchStructure (BinBox)
    void SearchNearestPointLocalInner( PointerType& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box )
    {
        if( mPointBegin == mPointEnd )
            return;

        bool Found;

        // initial search
        ++Box;
        SearchNearestInBoxInner( ThisPoint, rResult, rResultDistance, Box, Found );
        // increase mBox and try again
        while(!Found)
        {
            ++Box;
            SearchNearestInBoxInner( ThisPoint, rResult, rResultDistance, Box, Found );
        }

    }

    //************************************************************************
    //************************************************************************

    SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results,
                             DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults )
    {
        CoordinateType Radius2 = Radius * Radius;
        SizeType NumberOfResults = 0;
        SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }

    //************************************************************************

    SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results,
                             DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
    {
        CoordinateType Radius2 = Radius * Radius;
        SizeType NumberOfResults = 0;
        Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }

    //************************************************************************

    void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults ) override
    {
        SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
    }

    //************************************************************************

    void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box ) override
    {
        Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
    }
    
    //************************************************************************
    
   void SearchInRadius( PointerType const& ThisPoints, SizeType const& NumberOfPoints, std::vector<CoordinateType> const& Radius, std::vector<IteratorType> Results,
                        std::vector<DistanceIteratorType> ResultsDistances, std::vector<SizeType>& NumberOfResults, SizeType const& MaxNumberOfResults )
    {
        #pragma omp parallel for
        for(int k=0; k< NumberOfPoints; k++)
	  NumberOfResults[k] = SearchInRadius((&(*ThisPoints))[k],Radius[k],Results[k],ResultsDistances[k],MaxNumberOfResults);
    }

    //************************************************************************

    // **** THREAD SAFE

    // Dimension = 1
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        SearchRadiusInRange()(*(Box.RowBegin),*(Box.RowEnd),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
            SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 3
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
            for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
                SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
    }

    //************************************************************************
    //************************************************************************

    SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, IteratorType Results, SizeType MaxNumberOfResults )
    {
        CoordinateType Radius2 = Radius * Radius;
        SizeType NumberOfResults = 0;
        SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }

    //************************************************************************

    SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, IteratorType Results,
                             SizeType MaxNumberOfResults, SearchStructureType& Box )
    {
        CoordinateType Radius2 = Radius * Radius;
        SizeType NumberOfResults = 0;
        Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }

    //************************************************************************

    void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults ) override
    {
        SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
    }

    //************************************************************************

    void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box ) override
    {
        Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
        SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
    }

    //************************************************************************

    // Dimension = 1
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        SearchRadiusInRange()(*(Box.RowBegin),*(Box.RowEnd),ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
            SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 3
    void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
                              SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                              SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
            for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
                SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
    }

    //************************************************************************
    //************************************************************************

    // Dimension = 1
    void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance,
                             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, bool& Found )
    {
        Found = false;
        SearchNearestInRange()( *(Box.RowBegin), *(Box.RowEnd), ThisPoint, ResultPoint, ResultDistance, Found );
    }

    // Dimension = 2
    void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance,
                             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, bool& Found )
    {
        Found = false;
        for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
            SearchNearestInRange()( Box.RowBegin[I], Box.RowEnd[I], ThisPoint, ResultPoint, ResultDistance, Found );
    }

    // Dimension = 3
    void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance,
                             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, bool& Found )
    {
        Found = false;
        for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
            for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
                SearchNearestInRange()( Box.RowBegin[I], Box.RowEnd[I], ThisPoint, ResultPoint, ResultDistance, Found );
    }

    //************************************************************************
    //************************************************************************

    // Dimension = 1
    void SearchNearestInBoxInner( PointerType& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance,
                                  SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, bool& Found )
    {
        Found = false;
        SearchNearestInnerInRange( *(Box.RowBegin), *(Box.RowEnd), ThisPoint, ResultPoint, ResultDistance, Found );
    }

    // Dimension = 2
    void SearchNearestInBoxInner( PointerType& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance,
                                  SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, bool& Found )
    {
        Found = false;
        for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
            SearchNearestInnerInRange( Box.RowBegin[I], Box.RowEnd[I], ThisPoint, ResultPoint, ResultDistance, Found );
    }

    // Dimension = 3
    void SearchNearestInBoxInner( PointerType& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance,
                                  SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, bool& Found )
    {
        Found = false;
        for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
            for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
                SearchNearestInnerInRange( Box.RowBegin[I], Box.RowEnd[I], ThisPoint, ResultPoint, ResultDistance, Found );
    }

    //************************************************************************
    //************************************************************************

    void SearchNearestInnerInRange( const IteratorType& RangeBegin, const IteratorType& RangeEnd, PointerType& ThisPoint,
                                    PointerType& Result, CoordinateType& Distance, bool& Found )
    {
        CoordinateType NewDistance;
        for(IteratorType Point = RangeBegin ; Point != RangeEnd ; Point++)
        {
            NewDistance = TDistanceFunction()(**Point,*ThisPoint);
            if( NewDistance < Distance && *Point != ThisPoint)
            {
                Result = *Point;
                Distance = NewDistance;
                Found = true;
            }
        }
    }


    //************************************************************************
    //************************************************************************

    SizeType SearchInBox( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType Results,
                          SizeType MaxNumberOfResults )
    {
        SizeType NumberOfResults = 0;
        SearchStructureType Box( CalculateCell(SearchMinPoint), CalculateCell(SearchMaxPoint), mN, mIndexCellBegin );
        SearchInBoxLocal( SearchMinPoint, SearchMaxPoint, Results, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }

    //************************************************************************

    void SearchInBox(PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& Results, SizeType& NumberOfResults,
                     SizeType const& MaxNumberOfResults ) override
    {
        NumberOfResults = 0;
        SearchStructureType Box( CalculateCell(SearchMinPoint), CalculateCell(SearchMaxPoint), mN, mIndexCellBegin );
        SearchInBoxLocal( SearchMinPoint, SearchMaxPoint, Results, NumberOfResults, MaxNumberOfResults, Box );
    }

    //************************************************************************

    // Dimension = 1
    void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
                           SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                           SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,*(Box.RowBegin),*(Box.RowEnd),ResultsPoint,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
                           SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                           SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
            SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,Box.RowBegin[I],Box.RowEnd[I],ResultsPoint,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 3
    void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
                           SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
                           SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
            for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
                SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,Box.RowBegin[I],Box.RowEnd[I],ResultsPoint,NumberOfResults,MaxNumberOfResults);
    }

    //************************************************************************
    //************************************************************************

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "BinsContainer";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BinsContainer";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const override
    {
        rOStream << Perfix << "Bin[" << SearchUtils::PointerDistance(mPointBegin, mPointEnd) << "] : " << std::endl;
        for(IteratorConstIterator i_cell = mIndexCell.begin() ; i_cell != mIndexCell.end()-1 ; i_cell++)
        {
            rOStream << Perfix << "[ " ;
            for(IteratorType i_point = *i_cell ; i_point != *(i_cell+1) ; i_point++)
                rOStream << **i_point << " ";
            rOStream << "]" << std::endl;
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

    /// Assignment operator.
    Bins& operator=(Bins const& rOther);

    /// Copy constructor.
    Bins(Bins const& rOther);

private:

    // Point Access Iterators (vector reordered!!)
    IteratorType    mPointBegin;
    IteratorType    mPointEnd;

    // Bin Parameters (Sizes,BoundingBox,...)
    PointType       mMinPoint;
    PointType       mMaxPoint;
    CoordinateArray mCellSize;
    CoordinateArray mInvCellSize;
    SizeArray       mN;

    // Bins Access Vector ( vector<Iterator> )
    IteratorVector           mIndexCell;
    IteratorIterator         mIndexCellBegin;
    IteratorIterator         mIndexCellEnd;

    // Work Variables ( For non-copy of Search Variables )
    //SearchStructureType mBox;

public:

    //TODO: check -- changed to avoid copy construction
//     static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, PointType MaxPoint, PointType MinPoint, SizeType BucketSize)
    static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, const PointType& MaxPoint, const PointType& MinPoint, SizeType BucketSize)
    {
        SizeType number_of_points = SearchUtils::PointerDistance(PointsBegin,PointsEnd);
        if (number_of_points == 0)
            return NULL;
        else
        {
            return new Bins( PointsBegin, PointsEnd, MinPoint, MaxPoint, BucketSize );
        }
    }


};

template< std::size_t TDimension, class TPointType, class TContainerType, class TPointerType,
          class TIteratorType, class TDistanceIteratorType, class TDistanceFunction >
std::ostream & operator<<( std::ostream& rOStream, Bins<TDimension,TPointType,TContainerType,TPointerType,TIteratorType,TDistanceIteratorType,TDistanceFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintSize(rOStream);
    rThis.PrintData(rOStream);
    return rOStream;
}



}

#endif // KRATOS_BINS_CONTAINER_H_INCLUDE
