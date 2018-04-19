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
#include "DEM_application.h"
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

BinsObjectDynamicPeriodic (IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd, const array_1d<double, 3> domain_min, const array_1d<double, 3> domain_max)
{
    // : BinsObjectDynamic<TConfigure>(ObjectsBegin, ObjectsEnd)

    // we must repeat these operations here if we want to still derive from the regular BinsObjectDynamic (template) class
    this->mObjectsBegin = ObjectsBegin;
    this->mObjectsEnd = ObjectsEnd;
    this->mObjectsSize = SearchUtils::PointerDistance(this->mObjectsBegin, this->mObjectsEnd);
    SetDomainLimits(domain_min, domain_max);
    this->CalculateCellSize(this->mObjectsSize);
    this->AllocateContainer();
    GenerateBins();
}

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
  //
  Box.Set(this->CalculateCell(Low), this->CalculateCell(High), this->mN );

  SearchInRadiusExclusivePeriodic(ThisObject, Radius, Results, NumberOfResults, MaxNumberOfResults, Box);

  return NumberOfResults;
}

SizeType SearchObjectsInRadiusExclusive(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results, DistanceIteratorType ResultDistances, const SizeType& MaxNumberOfResults) override
{
  return SearchObjectsInRadiusExclusive(ThisObject, Radius, Results, MaxNumberOfResults);
}

protected:

void GenerateBins() override
{
    PointType Low, High;
    SearchStructureType Box;
    /// Fill container with objects
    int i = 0;
    for(IteratorType i_object = this->mObjectsBegin ; i_object != this->mObjectsEnd ; i_object++)
    {
        ++i;
        TConfigure::CalculateBoundingBox(*i_object, Low, High);
        Box.Set( this->CalculateCell(Low), this->CalculateCell(High), this->mN );
        FillObjectPeriodic(Box, *i_object);
    }
}

void SearchInRadiusExclusivePeriodic(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                  SearchStructureType& Box )
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

  int III_size = int(Box.Axis[2].Size()) + 1;

  for (IndexType III = III_begin; III_size > 0; NextIndex(III, III_size, MinCell, MaxCell, Box.Axis, 2))
  {
      MinCell[1] = MinBox[1];
      MaxCell[1] = MaxBox[1];
      int II_size = int(Box.Axis[1].Size()) + 1;

      for (IndexType II = II_begin; II_size > 0; NextIndex(II, II_size, MinCell, MaxCell, Box.Axis, 1))
      {
          MinCell[0] = MinBox[0];
          MaxCell[0] = MaxBox[0];
          int I_size = int(Box.Axis[0].Size()) + 1;

          for (IndexType I = I_begin; I_size > 0; NextIndex(I, I_size, MinCell, MaxCell, Box.Axis, 0))
          {
              IndexType GlobalIndex = III * Box.Axis[2].Block + II * Box.Axis[1].Block + I * Box.Axis[0].Block;
              //this->mCells[GlobalIndex].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);

              if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
              {
                  this->mCells[GlobalIndex].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);
              }
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
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * this->mCellSize[i] + this->mMinPoint[i];
        MaxBox[i] = MinBox[i] + this->mCellSize[i];
    }

    MinCell[2] = MinBox[2];
    MaxCell[2] = MaxBox[2];
    IndexType I_begin = Box.Axis[0].BeginIndex();
    IndexType II_begin = Box.Axis[1].BeginIndex();
    IndexType III_begin = Box.Axis[2].BeginIndex();

    int III_size = int(Box.Axis[2].Size()) + 1;

    for (IndexType III = III_begin; III_size > 0; NextIndex(III, III_size, MinCell, MaxCell, Box.Axis, 2))
    {
        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        int II_size = int(Box.Axis[1].Size()) + 1;

        for (IndexType II = II_begin; II_size > 0; NextIndex(II, II_size, MinCell, MaxCell, Box.Axis, 1))
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            int I_size = int(Box.Axis[0].Size()) + 1;

            for (IndexType I = I_begin; I_size > 0; NextIndex(I, I_size, MinCell, MaxCell, Box.Axis, 0))
            {
                IndexType GlobalIndex = III * Box.Axis[2].Block + II * Box.Axis[1].Block + I * Box.Axis[0].Block;
                //this->mCells[GlobalIndex].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, ResultDistances, NumberOfResults, MaxNumberOfResults);
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                {
                    this->mCells[GlobalIndex].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, ResultDistances, NumberOfResults, MaxNumberOfResults);
                }
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

    int III_size = int(Box.Axis[2].Size()) + 1;

    for (IndexType III = III_begin; III_size > 0; NextIndex(III, III_size, MinCell, MaxCell, Box.Axis, 2))
    {
        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        int II_size = int(Box.Axis[1].Size()) + 1;

        for (IndexType II = II_begin; II_size > 0; NextIndex(II, II_size, MinCell, MaxCell, Box.Axis, 1))
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            int I_size = int(Box.Axis[0].Size()) + 1;

            for (IndexType I = I_begin; I_size > 0; NextIndex(I, I_size, MinCell, MaxCell, Box.Axis, 0))
            {
                IndexType GlobalIndex = III * Box.Axis[2].Block + II * Box.Axis[1].Block + I * Box.Axis[0].Block;
                this->mCells[GlobalIndex].Add(i_object);
//                if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell)){
//                    this->mCells[GlobalIndex].Add(i_object);
//                }
            }
        }
    }
}


private:

array_1d<double, 3> mDomainMin;
array_1d<double, 3> mDomainMax;

void SetDomainLimits(const array_1d<double, 3> domain_min, const array_1d<double, 3> domain_max)
{
    mDomainMin[0] = domain_min[0];
    mDomainMin[1] = domain_min[1];
    mDomainMin[2] = domain_min[2];
    mDomainMax[0] = domain_max[0];
    mDomainMax[1] = domain_max[1];
    mDomainMax[2] = domain_max[2];

//    for (SizeType i = 0 ; i < Dimension ; i++){
//        this->mMaxPoint[i] = std::min(this->mMaxPoint[i], mDomainMax[i]);
//        this->mMinPoint[i] = std::max(this->mMinPoint[i], mDomainMin[i]);
//    }

    for (SizeType i = 0 ; i < Dimension ; i++){
        this->mMaxPoint[i] = mDomainMax[i];
        this->mMinPoint[i] = mDomainMin[i];
    }
}

inline void NextIndex(IndexType& Index, int& box_size_counter, PointType& MinCell, PointType& MaxCell, const SubBinAxisPeriodic<IndexType,SizeType> Axis[3], const SizeType dimension)
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

    --box_size_counter;
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
        this-> mMinPoint           = rOther.mMinPoint;
        this-> mMaxPoint           = rOther.mMaxPoint;
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
