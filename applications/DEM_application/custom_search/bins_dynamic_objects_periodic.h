//  Author:    Guillermo Casas

#if !defined(KRATOS_BINS_DYNAMIC_OBJECTS_CONTAINER_PERIODIC_H_INCLUDED)
#define  KRATOS_BINS_DYNAMIC_OBJECTS_CONTAINER_PERIODIC_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "includes/define.h"
#include "spatial_containers/bins_dynamic_objects.h"
#include "spatial_containers/search_structure.h"
#include "search_structure_periodic.h"

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
template<class TConfigure>
class BinsObjectDynamicPeriodic: public BinsObjectDynamic<TConfigure> {
public:
///@name Type Definitions
///@{
enum { Dimension = TConfigure::Dimension };
typedef BinsObjectDynamic<TConfigure>                BaseClassBins;
typedef TConfigure                                   Configure;
typedef typename TConfigure::PointType               PointType;
typedef typename TConfigure::PointerType             PointerType;
typedef typename TConfigure::ContainerType           ContainerType;
typedef typename TConfigure::IteratorType            IteratorType;
typedef typename TConfigure::ResultContainerType     ResultContainerType;
typedef typename TConfigure::ResultIteratorType      ResultIteratorType;
typedef typename TConfigure::DistanceIteratorType    DistanceIteratorType;

typedef Cell<Configure> CellType;
typedef std::vector<CellType> CellContainerType;
typedef typename CellContainerType::iterator CellContainerIterator;

typedef TreeNode<Dimension, PointType, PointerType, IteratorType,  typename TConfigure::DistanceIteratorType> TreeNodeType;
typedef typename TreeNodeType::CoordinateType  CoordinateType;  // double
typedef typename TreeNodeType::SizeType        SizeType;        // std::size_t
typedef typename TreeNodeType::IndexType       IndexType;       // std::size_t


typedef Tvector<IndexType,Dimension>      IndexArray;
typedef Tvector<SizeType,Dimension>       SizeArray;
typedef Tvector<CoordinateType,Dimension> CoordinateArray;

///Contact Pair
typedef typename TConfigure::ContainerContactType  ContainerContactType;
typedef typename TConfigure::IteratorContactType IteratorContactType;

///typedef TreeNodeType LeafType;
typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
//typedef typename TreeNodeType::SearchStructureType  SearchStructureType;
typedef SearchStructurePeriodic<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,Dimension> SearchStructureType;


/// Pointer definition of BinsObjectDynamicPeriodic
KRATOS_CLASS_POINTER_DEFINITION(BinsObjectDynamicPeriodic);

///@}
///@name Life Cycle
///@{

/// Default constructor.
BinsObjectDynamicPeriodic(){}
/// Constructor de bins a bounding box

BinsObjectDynamicPeriodic (IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd)
  : BinsObjectDynamic<TConfigure>(ObjectsBegin, ObjectsEnd)
{}

BinsObjectDynamicPeriodic (IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd, const double domain_min[3], const double domain_max[3])
  : BinsObjectDynamic<TConfigure>(ObjectsBegin, ObjectsEnd)
{
    SetDomainLimits(domain_min, domain_max);
}

BinsObjectDynamicPeriodic (IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd, CoordinateType CellSize)
  : BinsObjectDynamic<TConfigure>(ObjectsBegin, ObjectsEnd, CellSize)
{}

BinsObjectDynamicPeriodic (const PointType& MinPoint, const PointType& MaxPoint, CoordinateType CellSize)
  : BinsObjectDynamic<TConfigure>(MinPoint, MaxPoint, CellSize)
{}

BinsObjectDynamicPeriodic (const PointType& MinPoint, const PointType& MaxPoint, SizeType NumPoints)
  : BinsObjectDynamic<TConfigure>(MinPoint, MaxPoint, NumPoints)
{}

/// Destructor.
virtual ~BinsObjectDynamicPeriodic(){}

IndexType CalculatePosition(CoordinateType const& ThisCoord, const SizeType& ThisDimension) override
{
CoordinateType d_index = ThisCoord;
if (ThisCoord < this->mDomainMin[ThisDimension]){
    CoordinateType domain_period = this->mDomainMax[ThisDimension] - this->mDomainMin[ThisDimension];
    d_index += domain_period;
}

else if (ThisCoord > this->mDomainMax[ThisDimension]){
    CoordinateType domain_period = this->mDomainMax[ThisDimension] - this->mDomainMin[ThisDimension];
    d_index -= domain_period;
}

d_index -= this->mMinPoint[ThisDimension];

IndexType index = static_cast<IndexType>(d_index * this->mInvCellSize[ThisDimension]);
return index;
}

SizeType SearchObjectsInRadiusExclusive(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results, const SizeType& MaxNumberOfResults) override
{
  PointType Low, High;
  SearchStructureType Box;
  SizeType NumberOfResults = 0;

  TConfigure::CalculateBoundingBox(ThisObject, Low, High, Radius);

  // Note that here we allow that Low > High:
  Box.Set(this->CalculateCell(Low), this->CalculateCell(High), this->mN );

  SearchInRadiusExclusivePeriodic(ThisObject, Radius, Results, NumberOfResults, MaxNumberOfResults, Box);

  return NumberOfResults;
}

SizeType SearchObjectsInRadiusExclusive(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results, DistanceIteratorType ResultDistances, const SizeType& MaxNumberOfResults) override
{
  return SearchObjectsInRadiusExclusive(ThisObject, Radius, Results, MaxNumberOfResults);
}

protected:

void CalculateBoundingBox() override
{
    PointType Low, High;
    TConfigure::CalculateBoundingBox(*this->mObjectsBegin, this->mMinPoint, this->mMaxPoint);

#ifdef _OPENMP
    SizeType number_of_threads = omp_get_max_threads();
#else
    SizeType number_of_threads = 1;
#endif

    std::vector<SizeType> node_partition;
    this->CreatePartition(number_of_threads, this->mObjectsSize, node_partition);

    std::vector<PointType> Max(number_of_threads);
    std::vector<PointType> Min(number_of_threads);

    for(SizeType k=0; k<number_of_threads; k++ )
    {
        Max[k] = this->mMaxPoint;
        Min[k] = this->mMinPoint;
    }

    IteratorType i_begin = this->mObjectsBegin;
    IteratorType i_end   = this->mObjectsEnd;

    for (IteratorType i_object = i_begin ; i_object != i_end ; i_object++ )
    {
        TConfigure::CalculateBoundingBox(*i_object, Low, High);
        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            this->mMaxPoint[i] = (this->mMaxPoint[i] < High[i]) ? High[i] : this->mMaxPoint[i];
            this->mMinPoint[i] = (this->mMinPoint[i] > Low[i])  ? Low[i]  : this->mMinPoint[i];
        }
    }

    PointType Epsilon = this->mMaxPoint - this->mMinPoint;

    for(SizeType i = 0 ; i < Dimension ; i++)
    {
        this->mMaxPoint[i] += Epsilon[i] * 0.01;
        this->mMinPoint[i] -= Epsilon[i] * 0.01;
    }

    for(SizeType i = 0 ; i < Dimension ; i++){
        this->mMaxPoint[i] = this->mDomainMax[i];
        this->mMinPoint[i] = this->mDomainMin[i];
    }
}

void SearchInRadiusExclusivePeriodic(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                  SearchStructureType& Box )
{

  PointType  MinCell, MaxCell;
  PointType  MinBox, MaxBox;

  for(SizeType i = 0; i < 3; i++)
  {
      MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * this->mCellSize[i] + this->mMinPoint[i];  //
      MaxBox[i] = MinBox[i] + this->mCellSize[i];
  }

  MinCell[2] = MinBox[2];
  MaxCell[2] = MaxBox[2];
  IndexType I_begin = Box.Axis[0].BeginIndex();
  IndexType II_begin = Box.Axis[1].BeginIndex();
  IndexType III_begin = Box.Axis[2].BeginIndex();
  IndexType I_end = Box.Axis[0].EndIndex();
  IndexType II_end = Box.Axis[1].EndIndex();
  IndexType III_end = Box.Axis[2].EndIndex();

  int III_box_size = int(Box.Axis[2].Size()) + 1;

  for(IndexType III = III_begin; III_box_size > 0; NextIndex(III, III_box_size, MinCell, MaxCell, Box.Axis, 2))
  {
      MinCell[1] = MinBox[1];
      MaxCell[1] = MaxBox[1];
      int II_box_size = int(Box.Axis[1].Size()) + 1;

      for(IndexType II = II_begin; II_box_size > 0; NextIndex(II, II_box_size, MinCell, MaxCell, Box.Axis, 1))
      {
          MinCell[0] = MinBox[0];
          MaxCell[0] = MaxBox[0];
          int I_box_size = int(Box.Axis[0].Size()) + 1;

          for(IndexType I = I_begin; I_box_size > 0; NextIndex(I, I_box_size, MinCell, MaxCell, Box.Axis, 0))
          {
              IndexType GlobalIndex = III * Box.Axis[2].Block + II * Box.Axis[1].Block + I * Box.Axis[0].Block;
              this->mCells[GlobalIndex].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);

//              if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
//              {
//                  this->mCells[GlobalIndex].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);
//              }
          }
      }
   }
}

void SearchInRadiusExclusivePeriodic(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                    SearchStructureType& Box )
{


    PointType  MinCell, MaxCell;
    PointType  MinBox, MaxBox;

    for(SizeType i = 0; i < 3; i++)
    {
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * this->mCellSize[i] + this->mMinPoint[i];  //
        MaxBox[i] = MinBox[i] + this->mCellSize[i];
    }

    IndexType I_begin = Box.Axis[0].BeginIndex();
    IndexType II_begin = Box.Axis[1].BeginIndex();
    IndexType III_begin = Box.Axis[2].BeginIndex();
    IndexType I_end = Box.Axis[0].EndIndex();
    IndexType II_end = Box.Axis[1].EndIndex();
    IndexType III_end = Box.Axis[2].EndIndex();

    int III_box_size = int(Box.Axis[2].Size()) + 1;

    for(IndexType III = III_begin; III_box_size > 0; NextIndex(III, III_box_size, MinCell, MaxCell, Box.Axis, 2))
    {
        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        int II_box_size = int(Box.Axis[1].Size()) + 1;

        for(IndexType II = II_begin; II_box_size > 0; NextIndex(II, II_box_size, MinCell, MaxCell, Box.Axis, 1))
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            int I_box_size = int(Box.Axis[0].Size()) + 1;

            for(IndexType I = I_begin; I_box_size > 0; NextIndex(I, I_box_size, MinCell, MaxCell, Box.Axis, 0))
            {
                IndexType GlobalIndex = III * Box.Axis[2].Block + II * Box.Axis[1].Block + I * Box.Axis[0].Block;
                this->mCells[GlobalIndex].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, ResultDistances, NumberOfResults, MaxNumberOfResults);
//                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
//                {
//                    this->mCells[GlobalIndex].SSearchObjectsInRadiusExclusiveearchObjectsInRadiusExclusive(ThisObject, Radius, Result, ResultDistances, NumberOfResults, MaxNumberOfResults);
//                }
            }
        }
    }
}

void FillObjectPeriodic(SearchStructureType& Box, const PointerType& i_object)
{
    PointType  MinCell, MaxCell;
    PointType  MinBox, MaxBox;

    for(SizeType i = 0; i < 3; i++)
    {
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * this->mCellSize[i] + this->mMinPoint[i];
        MaxBox[i] = MinBox[i] + this->mCellSize[i];
    }

    MinCell[2] = MinBox[2];
    MaxCell[2] = MaxBox[2];
    IndexType I_begin = Box.Axis[0].BeginIndex();
    IndexType II_begin = Box.Axis[1].BeginIndex();
    IndexType III_begin = Box.Axis[2].BeginIndex();

    int III_box_size = int(Box.Axis[2].Size()) + 1;

    for(IndexType III = III_begin; III_box_size > 0; NextIndex(III, III_box_size, MinCell, MaxCell, Box.Axis, 2))
    {
        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        int II_box_size = int(Box.Axis[1].Size()) + 1;

        for(IndexType II = II_begin; II_box_size > 0; NextIndex(II, II_box_size, MinCell, MaxCell, Box.Axis, 1))
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            int I_box_size = int(Box.Axis[0].Size()) + 1;

            for(IndexType I = I_begin; I_box_size > 0; NextIndex(I, I_box_size, MinCell, MaxCell, Box.Axis, 0))
            {
                IndexType GlobalIndex = III * Box.Axis[2].Block + II * Box.Axis[1].Block + I * Box.Axis[0].Block;
                this->mCells[GlobalIndex].Add(i_object);
//                if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
//                    this->mCells[GlobalIndex].Add(i_object);
            }
        }
    }
}


private:

double mDomainMin[3];
double mDomainMax[3];

void SetDomainLimits(const double domain_min[3], const double domain_max[3])
{
    mDomainMin[0] = domain_min[0];
    mDomainMin[1] = domain_min[1];
    mDomainMin[2] = domain_min[2];
    mDomainMax[0] = domain_max[0];
    mDomainMax[1] = domain_max[1];
    mDomainMax[2] = domain_max[2];

    for (SizeType i = 0 ; i < Dimension ; i++){
        this->mMaxPoint[i] = std::min(this->mMaxPoint[i], mDomainMax[i]);
        this->mMinPoint[i] = std::max(this->mMinPoint[i], mDomainMin[i]);
    }
}

inline void NextIndex(IndexType& Index, int& box_size_count, PointType& MinCell, PointType& MaxCell, const SubBinAxisPeriodic<IndexType,SizeType> Axis[3], const SizeType dimension)
{
    if (Index < this->mN[dimension] - 1){
        Index += 1;
        MinCell[dimension] += this->mCellSize[dimension];
        MaxCell[dimension] += this->mCellSize[dimension];
    }
    else {
        Index = 0;
        MinCell[dimension] += - MinCell[dimension];
        MaxCell[dimension] += - MaxCell[dimension] + this->mCellSize[dimension];
    }

    --box_size_count;
}


    ///@}
    ///@name Private Operators
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
    ///@}


public:
    /// Assignment operator.
    BinsObjectDynamicPeriodic<TConfigure> & operator=(const BinsObjectDynamicPeriodic<TConfigure> & rOther)
    {
        this-> mMinPoint            = rOther.mMinPoint;
        this-> mMaxPoint            = rOther.mMaxPoint;
        this->mObjectsBegin        = rOther.mObjectsBegin;
        this->mObjectsEnd          = rOther.mObjectsEnd;
        this->mObjectsSize         = rOther.mObjectsSize;
        this->mCellSize            = rOther.mCellSize;
        this->mInvCellSize         = rOther.mInvCellSize;
        this->mN                   = rOther.mN;
        this->mCells               = rOther.mCells;
        return *this;
    }

    /// Copy constructor.
    BinsObjectDynamicPeriodic(const BinsObjectDynamicPeriodic& rOther)
    {
        *this =  rOther;
    }

    /// Copy constructor.
    template<class T>
    BinsObjectDynamicPeriodic(const BinsObjectDynamicPeriodic<T>& rOther)
    {
        *this =  rOther;
    }

}; // Class BinsObjectDynamicPeriodic

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TConfigure>
inline std::istream& operator >> (std::istream& rIStream,
                                  BinsObjectDynamicPeriodic<TConfigure>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TConfigure>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BinsObjectDynamicPeriodic<TConfigure> & rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
