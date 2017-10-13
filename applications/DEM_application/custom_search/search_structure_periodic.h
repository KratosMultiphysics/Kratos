// $Author: Guillermo Casas $

#if !defined(KRATOS_SEARCH_STRUCTURE_PERIODIC_H_INCLUDED )
#define  KRATOS_SEARCH_STRUCTURE_PERIODIC_H_INCLUDED
// System includes
#include <vector>
#include <cfloat>
// External includes

// Project includes
#include "includes/define.h"
#include "spatial_containers/search_structure.h"

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


/// TOOLS UTILS ///

template< class IndexType, class SizeType>
class SubBinAxisPeriodic
{
public:
    IndexType Min;
    IndexType Max;
    IndexType MaxSize;
    IndexType Block;

    SubBinAxisPeriodic() : Min(0), Max(0), MaxSize(0), Block(1) {}

    SubBinAxisPeriodic(IndexType const& Min_, IndexType const& Max_, IndexType const& MaxSize_, IndexType const& Block_)
    {
        Set(Min_,Max_,MaxSize_,Block_);
    }
    ~SubBinAxisPeriodic() {}
    void Set(IndexType const& iCell, IndexType const& MaxSize_, IndexType const& Block_)
    {
        Set(iCell,iCell,MaxSize_,Block_);
    }
    void Set(IndexType const& Min_, IndexType const& Max_, IndexType const& MaxSize_, IndexType const& Block_)
    {
        MaxSize = MaxSize_;
        Block = Block_;
        Min = Min_;
        Max = Max_;
//        Min = (Min_ >= 0) ?           Min_ : Min_ + MaxSize_;
//        Max = (Max_ < MaxSize_ - 1) ? Max_ : Max_ - MaxSize_ + 1;
//        Min = std::max( Min_, static_cast<IndexType>(0) );
//        Max = std::min( Max_, MaxSize-1 );
    }
    IndexType Begin()
    {
        return Min*Block;
    }
    IndexType End()
    {
        return Max*Block;
    }
    IndexType BeginIndex()
    {
        return Min;
    }
    IndexType EndIndex()
    {
        return Max;
    }
    SizeType Size()
    {
        if (Max >= Min){
            return static_cast<SizeType>(Max-Min);
        }
        else {
            return static_cast<SizeType>(MaxSize - 1 + (Max - Min));
        }
    }
    SubBinAxisPeriodic const& operator++()
    {
        if( Min > static_cast<IndexType>(0) ) Min--;
        if( Max < MaxSize-1 ) Max++;
        return *this;
    }
    SubBinAxisPeriodic const& operator--()
    {
        Min++;
        Max--;
        return *this;
    }
};


template<
class IndexType,
class SizeType,
class CoordinateType,
class IteratorType,
class IteratorIteratorType,
std::size_t Dimension >
class SearchStructurePeriodic
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SearchStructurePeriodic
    KRATOS_CLASS_POINTER_DEFINITION(SearchStructurePeriodic);

    typedef Tvector<IndexType,Dimension> IndexVector;
    typedef Tvector<SizeType,Dimension> SizeVector;

    typedef SearchStructurePeriodic<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,Dimension> ThisType;

public:
    // Bin
    //SubBinAxis<IndexType,SizeType> Axis[Dimension];
    SubBinAxisPeriodic<IndexType,SizeType> Axis[3];
    IteratorIteratorType RowBegin;
    IteratorIteratorType RowEnd;
    IteratorIteratorType DataBegin;
    // KDTree
    CoordinateType distance_to_partition;
    CoordinateType distance_to_partition2;
    CoordinateType residual_distance[Dimension];
    SizeType BucketCounter;

public:

    SearchStructurePeriodic() {}

    SearchStructurePeriodic( IndexVector const& Min_, IndexVector const& Max_, SizeVector const& MaxSize_, IteratorIteratorType const& IteratorBegin )
    {
        Set(Min_,Max_,MaxSize_,IteratorBegin);
    }

    SearchStructurePeriodic( IndexVector const& IndexCell, SizeVector const& MaxSize_, IteratorIteratorType const& IteratorBegin )
    {
        Set(IndexCell,IndexCell,MaxSize_,IteratorBegin);
    }

    SearchStructurePeriodic( IndexVector const& Min_, IndexVector const& Max_, SizeVector const& MaxSize_ )
    {
        Set(Min_,Max_,MaxSize_);
    }

    SearchStructurePeriodic( IndexVector const& IndexCell, SizeVector const& MaxSize_ )
    {
        Set(IndexCell,IndexCell,MaxSize_);
    }

    ~SearchStructurePeriodic() {}

    void Set( IndexVector const& IndexCell, SizeVector const& _MaxSize, IteratorIteratorType const& IteratorBegin )
    {
        Set( IndexCell, IndexCell, _MaxSize, IteratorBegin );
    }

    void Set( IndexVector const& Min_, IndexVector const& Max_, SizeVector const& MaxSize_, IteratorIteratorType const& IteratorBegin )
    {
        IndexType Block = 1;
        Axis[0].Set(Min_[0],Max_[0],MaxSize_[0],Block);

        for(SizeType i = 1; i < Dimension; i++)
        {
            Block *= MaxSize_[i-1];
            Axis[i].Set(Min_[i],Max_[i],MaxSize_[i],Block);
        }
        
        DataBegin = IteratorBegin;
        
        RowBegin = DataBegin + Axis[0].Min;
        RowEnd   = DataBegin + Axis[0].Max + 1;
    }

    void Set( IndexVector const& IndexCell, SizeVector const& MaxSize_ )
    {
        Set(IndexCell,IndexCell,MaxSize_);
    }

    void Set( IndexVector const& Min_, IndexVector const& Max_, SizeVector const& MaxSize_ )
    {
        IndexType Block = 1;
        Axis[0].Set(Min_[0],Max_[0],MaxSize_[0],Block);

        for(SizeType i = 1; i < Dimension; i++)
        {
            Block *= MaxSize_[i-1];
            Axis[i].Set(Min_[i],Max_[i],MaxSize_[i],Block);
        }
    }

    IndexType BeginRow(IndexType const& Idx)
    {
        return Idx + Axis[0].Min;
    }
    IndexType EndRow(IndexType const& Idx)
    {
        return Idx + Axis[0].Max+1;
    }

    SearchStructurePeriodic const& operator++()
    {
        for(SizeType i = 0; i < Dimension; i++)
            ++(Axis[i]);
        
        RowBegin = DataBegin + Axis[0].Min;
        RowEnd   = DataBegin + Axis[0].Max + 1;
        
        return *this;
    }

    SearchStructurePeriodic const& operator--()
    {
        for(SizeType i = 0; i < Dimension; i++)
            (Axis[i])--;
        
        RowBegin = DataBegin + Axis[0].Min;
        RowEnd   = DataBegin + Axis[0].Max + 1;
        
        return *this;
    }
};

}  // namespace Kratos.

#endif
