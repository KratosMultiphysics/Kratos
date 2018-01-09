//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:       Nelson Lafontaine
//                      Carlos A. Roig

#if !defined(KRATOS_BINS_DYNAMIC_OBJECTS_CONTAINER_H_INCLUDED)
#define  KRATOS_BINS_DYNAMIC_OBJECTS_CONTAINER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <time.h>

// Project includes
#include "tree.h"
#include "cell.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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
class BinsObjectDynamic {
public:
  ///@name Type Definitions
  ///@{

  enum { Dimension = TConfigure::Dimension };

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
  typedef typename TreeNodeType::SearchStructureType  SearchStructureType;

  /// Pointer definition of BinsObjectDynamic
  KRATOS_CLASS_POINTER_DEFINITION(BinsObjectDynamic);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  BinsObjectDynamic() {}
  /// Constructor de bins a bounding box

  BinsObjectDynamic (IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd)
      : mObjectsBegin(ObjectsBegin), mObjectsEnd(ObjectsEnd) {

    mObjectsSize = SearchUtils::PointerDistance(mObjectsBegin,mObjectsEnd);
    CalculateBoundingBox();           // Calculate mMinPoint, mMaxPoint
    CalculateCellSize(mObjectsSize);  // Calculate number of Cells
    AllocateContainer();              // Allocate cell list
    GenerateBins();                   // Fill Cells with objects
  }

  BinsObjectDynamic (IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd, CoordinateType CellSize)
      : mObjectsBegin(ObjectsBegin), mObjectsEnd(ObjectsEnd) {

    mObjectsSize = SearchUtils::PointerDistance(mObjectsBegin,mObjectsEnd);
    CalculateBoundingBox();           // Calculate mMinPoint, mMaxPoint
    CalculateCellSize(CellSize);      // Calculate number of Cells
    AllocateContainer();              // Allocate cell list
    GenerateBins();                   // Fill Cells with objects
  }

  BinsObjectDynamic (const PointType& MinPoint, const PointType& MaxPoint, CoordinateType CellSize)
      : mObjectsSize(0), mObjectsBegin(0), mObjectsEnd(0) {

    for(SizeType i = 0; i < Dimension; i++) {
      mMinPoint[i] = MinPoint[i];
      mMaxPoint[i] = MaxPoint[i];
    }

    CalculateCellSize(CellSize);
    AllocateContainer();
  }

  BinsObjectDynamic (const PointType& MinPoint, const PointType& MaxPoint, SizeType NumPoints)
      : mObjectsSize(0), mObjectsBegin(0), mObjectsEnd(0) {

    for(SizeType i = 0; i < Dimension; i++) {
      mMinPoint[i] = MinPoint[i];
      mMaxPoint[i] = MaxPoint[i];
    }

    CalculateCellSize(NumPoints);
    AllocateContainer();
  }

  /// Destructor.
  virtual ~BinsObjectDynamic() {}

  /// Single search API

  /**
   * [SearchObjects description]
   * @param  ThisObject [description]
   * @param  Result     [description]
   * @return            [description]
   */
  SizeType SearchObjects(PointerType& ThisObject, ResultContainerType& Result) {
    PointType Low, High;
    SearchStructureType Box;

    TConfigure::CalculateBoundingBox(ThisObject, Low, High);
    Box.Set( CalculateCell(Low), CalculateCell(High), mN );
    SearchInBoxLocal(ThisObject, Result, Box );

    return Result.size();
  }

  /**
   * [SearchObjects description]
   * @param  ThisObject         [description]
   * @param  Result             [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  SizeType SearchObjects(PointerType& ThisObject, ResultIteratorType& Result, const SizeType& MaxNumberOfResults) {
    PointType Low, High;
    SearchStructureType Box;
    SizeType NumberOfResults = 0;

    TConfigure::CalculateBoundingBox(ThisObject, Low, High);
    Box.Set( CalculateCell(Low), CalculateCell(High), mN );
    SearchInBoxLocal(ThisObject, Result, NumberOfResults, MaxNumberOfResults, Box );

    return NumberOfResults;
  }

  /**
   * [SearchObjectsInCell description]
   * @param  ThisPoint          [description]
   * @param  Result             [description]
   * @return                    [description]
   */
  SizeType SearchObjectsInCell(const PointType& ThisPoint, ResultIteratorType Result) {
    /// Missing API for 'SearchObjectsInCell' without 'MaxNumberOfResults'
    KRATOS_ERROR << "Missing implementation of SearchObjectsInCell(PointerType, ResultIteratorType)" << std::endl;
  }

  /**
   * [SearchObjectsInCell description]
   * @param  ThisPoint          [description]
   * @param  Result             [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  SizeType SearchObjectsInCell(const PointType& ThisPoint, ResultIteratorType Result, const SizeType& MaxNumberOfResults) {
    IndexType icell = CalculateIndex(ThisPoint);

    if(mCells[icell].Size() < MaxNumberOfResults) {
      for(IteratorType i_object = mCells[icell].Begin() ; i_object != mCells[icell].End(); i_object++, Result++) {
        *Result = *i_object;
      }
      return mCells[icell].Size();
    } else {
      return -1;
    }
  }

  /**
   * [SearchObjectsExclusive description]
   * @param  ThisObject         [description]
   * @param  Result             [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  SizeType SearchObjectsExclusive(PointerType& ThisObject, ResultIteratorType& Result) {
    PointType Low, High;
    SearchStructureType Box;

    TConfigure::CalculateBoundingBox(ThisObject, Low, High);
    Box.Set( CalculateCell(Low), CalculateCell(High), mN );
    SearchObjectLocalExclusive(ThisObject, Result, Box );

    return Result.size();
  }

  /**
   * [SearchObjectsExclusive description]
   * @param  ThisObject         [description]
   * @param  Result             [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  SizeType SearchObjectsExclusive(PointerType& ThisObject, ResultIteratorType& Result, const SizeType& MaxNumberOfResults) {
    PointType Low, High;
    SearchStructureType Box;
    SizeType NumberOfResults = 0;

    TConfigure::CalculateBoundingBox(ThisObject, Low, High);
    Box.Set( CalculateCell(Low), CalculateCell(High), mN );
    SearchObjectLocalExclusive(ThisObject, Result, NumberOfResults, MaxNumberOfResults, Box );

    return NumberOfResults;
  }

  /**
   * [SearchObjectsInRadius description]
   * @param  ThisObject [description]
   * @param  Radius     [description]
   * @param  Results    [description]
   * @return            [description]
   */
  SizeType SearchObjectsInRadius(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results) {
    /// Missing API for 'SearchObjectsInRadius' without 'MaxNumberOfResults'
    KRATOS_ERROR << "Missing implementation of SearchObjectsInRadius(PointerType, const double, ResultIteratorType)" << std::endl;
  }

  /**
   * [SearchObjectsInRadius description]
   * @param  ThisObject         [description]
   * @param  Radius             [description]
   * @param  Results            [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  SizeType SearchObjectsInRadius(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results, const SizeType& MaxNumberOfResults) {
    PointType Low, High;
    SearchStructureType Box;
    SizeType NumberOfResults = 0;

    TConfigure::CalculateBoundingBox(ThisObject, Low, High, Radius);
    Box.Set( CalculateCell(Low), CalculateCell(High), mN );
    SearchInRadius(ThisObject, Radius, Results, NumberOfResults, MaxNumberOfResults, Box );

    return NumberOfResults;
  }

  /**
   * [SearchObjectsInRadius description]
   * @param  ThisObject      [description]
   * @param  Radius          [description]
   * @param  Results         [description]
   * @param  ResultDistances [description]
   * @return                 [description]
   */
  SizeType SearchObjectsInRadius(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results, DistanceIteratorType ResultDistances) {
    /// Missing API for 'SearchObjectsInRadius' without 'MaxNumberOfResults'
    KRATOS_ERROR << "Missing implementation of SearchObjectsInRadius(PointerType, const double, ResultIteratorType, DistanceIteratorType)" << std::endl;
  }

  /**
   * [SearchObjectsInRadius description]
   * @param  ThisObject         [description]
   * @param  Radius             [description]
   * @param  Results            [description]
   * @param  ResultDistances    [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  SizeType SearchObjectsInRadius(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results, DistanceIteratorType ResultDistances, const SizeType& MaxNumberOfResults) {
    PointType Low, High;
    SearchStructureType Box;
    SizeType NumberOfResults = 0;

    TConfigure::CalculateBoundingBox(ThisObject, Low, High, Radius);
    Box.Set( CalculateCell(Low), CalculateCell(High), mN );
    SearchInRadius(ThisObject, Radius, Results, ResultDistances, NumberOfResults, MaxNumberOfResults, Box );

    return NumberOfResults;
  }

  /**
   * [SearchObjectsInRadiusExclusive description]
   * @param  ThisObject         [description]
   * @param  Radius             [description]
   * @param  Results            [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  virtual SizeType SearchObjectsInRadiusExclusive(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results) {
    /// Missing API for 'SearchObjectsInRadiusExclusive' without 'MaxNumberOfResults'
    KRATOS_ERROR << "Missing implementation of SearchObjectsInRadiusExclusive(PointerType, const double, ResultIteratorType)" << std::endl;
  }

  /**
   * [SearchObjectsInRadiusExclusive description]
   * @param  ThisObject         [description]
   * @param  Radius             [description]
   * @param  Results            [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  virtual SizeType SearchObjectsInRadiusExclusive(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results, const SizeType& MaxNumberOfResults) {
    PointType Low, High;
    SearchStructureType Box;
    SizeType NumberOfResults = 0;

    TConfigure::CalculateBoundingBox(ThisObject, Low, High, Radius);
    Box.Set( CalculateCell(Low), CalculateCell(High), mN );
    SearchInRadiusExclusive(ThisObject, Radius, Results, NumberOfResults, MaxNumberOfResults, Box );

    return NumberOfResults;
  }

  /**
   * [SearchObjectsInRadiusExclusive description]
   * @param  ThisObject      [description]
   * @param  Radius          [description]
   * @param  Results         [description]
   * @param  ResultDistances [description]
   * @return                 [description]
   */
  virtual SizeType SearchObjectsInRadiusExclusive(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results, DistanceIteratorType ResultDistances) {
    /// Missing API for 'SearchObjectsInRadiusExclusive' without 'MaxNumberOfResults'
    KRATOS_ERROR << "Missing implementation of SearchObjectsInRadiusExclusive(PointerType, const double, ResultIteratorType, DistanceIteratorType)" << std::endl;
  }

  /**
   * [SearchObjectsInRadiusExclusive description]
   * @param  ThisObject         [description]
   * @param  Radius             [description]
   * @param  Results            [description]
   * @param  ResultDistances    [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  virtual SizeType SearchObjectsInRadiusExclusive(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results, DistanceIteratorType ResultDistances, const SizeType& MaxNumberOfResults) {
    PointType Low, High;
    SearchStructureType Box;
    SizeType NumberOfResults = 0;

    TConfigure::CalculateBoundingBox(ThisObject, Low, High, Radius);
    Box.Set( CalculateCell(Low), CalculateCell(High), mN );
    SearchInRadiusExclusive(ThisObject, Radius, Results, ResultDistances, NumberOfResults, MaxNumberOfResults, Box );

    return NumberOfResults;
  }

  /// Batch search API (needs to be extended with the missing functions)

  /**
   * [SearchObjectsInRadius description]
   * @param ThisObjects        [description]
   * @param NumberOfObjects    [description]
   * @param Radius             [description]
   * @param Results            [description]
   * @param NumberOfResults    [description]
   * @param MaxNumberOfResults [description]
   */
  void SearchObjectsInRadius(IteratorType const& ThisObjects, SizeType const& NumberOfObjects, std::vector<double>& Radius, std::vector<std::vector<PointerType> >& Results, std::vector<SizeType>& NumberOfResults, SizeType const& MaxNumberOfResults) {
    PointType Low, High;
    SearchStructureType Box;

    #pragma omp parallel for private(Low,High,Box)
    for(int i = 0; i < static_cast<int>(NumberOfObjects); i++) {
      ResultIteratorType ResultsPointer            = Results[i].begin();

      NumberOfResults[i] = 0;

      TConfigure::CalculateBoundingBox(ThisObjects[i], Low, High, Radius[i]);
      Box.Set( CalculateCell(Low), CalculateCell(High), mN );

      SearchInRadius(ThisObjects[i], Radius[i], ResultsPointer, NumberOfResults[i], MaxNumberOfResults, Box );
    }
  }

  /**
   * [SearchObjectsInRadius description]
   * @param ThisObjects        [description]
   * @param NumberOfObjects    [description]
   * @param Radius             [description]
   * @param Results            [description]
   * @param ResultsDistances   [description]
   * @param NumberOfResults    [description]
   * @param MaxNumberOfResults [description]
   */
  void SearchObjectsInRadius(IteratorType const& ThisObjects, SizeType const& NumberOfObjects, std::vector<double>& Radius, std::vector<std::vector<PointerType> >& Results, std::vector<std::vector<double> >& ResultsDistances, std::vector<SizeType>& NumberOfResults, SizeType const& MaxNumberOfResults) {
    PointType Low, High;
    SearchStructureType Box;

    #pragma omp parallel for private(Low,High,Box)
    for(int i = 0; i < static_cast<int>(NumberOfObjects); i++) {
      ResultIteratorType ResultsPointer            = Results[i].begin();
      DistanceIteratorType ResultsDistancesPointer = ResultsDistances[i].begin();

      NumberOfResults[i] = 0;

      TConfigure::CalculateBoundingBox(ThisObjects[i], Low, High, Radius[i]);
      Box.Set( CalculateCell(Low), CalculateCell(High), mN );

      SearchInRadius(ThisObjects[i], Radius[i], ResultsPointer, ResultsDistancesPointer, NumberOfResults[i], MaxNumberOfResults, Box );
    }
  }

  /**
   * [SearchObjectsInRadiusExclusive description]
   * @param ThisObjects        [description]
   * @param NumberOfObjects    [description]
   * @param Radius             [description]
   * @param Results            [description]
   * @param NumberOfResults    [description]
   * @param MaxNumberOfResults [description]
   */
virtual void SearchObjectsInRadiusExclusive(IteratorType const& ThisObjects, SizeType const& NumberOfObjects, std::vector<double>& Radius, std::vector<std::vector<PointerType> >& Results, std::vector<SizeType>& NumberOfResults, SizeType const& MaxNumberOfResults) {
    PointType Low, High;
    SearchStructureType Box;

    #pragma omp parallel for private(Low,High,Box)
    for(int i = 0; i < static_cast<int>(NumberOfObjects); i++) {
      ResultIteratorType ResultsPointer            = Results[i].begin();

      NumberOfResults[i] = 0;

      TConfigure::CalculateBoundingBox(ThisObjects[i], Low, High, Radius[i]);
      Box.Set( CalculateCell(Low), CalculateCell(High), mN );

      SearchInRadiusExclusive(ThisObjects[i], Radius[i], ResultsPointer, NumberOfResults[i], MaxNumberOfResults, Box );
    }
  }

  /**
   * [SearchObjectsInRadiusExclusive description]
   * @param ThisObjects        [description]
   * @param NumberOfObjects    [description]
   * @param Radius             [description]
   * @param Results            [description]
   * @param ResultsDistances   [description]
   * @param NumberOfResults    [description]
   * @param MaxNumberOfResults [description]
   */
virtual void SearchObjectsInRadiusExclusive(IteratorType const& ThisObjects, SizeType const& NumberOfObjects, std::vector<double>& Radius, std::vector<std::vector<PointerType> >& Results, std::vector<std::vector<double> >& ResultsDistances, std::vector<SizeType>& NumberOfResults, SizeType const& MaxNumberOfResults) {
    PointType Low, High;
    SearchStructureType Box;

    #pragma omp parallel for private(Low,High,Box)
    for(int i = 0; i < static_cast<int>(NumberOfObjects); i++) {
      ResultIteratorType ResultsPointer            = Results[i].begin();
      DistanceIteratorType ResultsDistancesPointer = ResultsDistances[i].begin();

      NumberOfResults[i] = 0;

      TConfigure::CalculateBoundingBox(ThisObjects[i], Low, High, Radius[i]);
      Box.Set( CalculateCell(Low), CalculateCell(High), mN );

      SearchInRadiusExclusive(ThisObjects[i], Radius[i], ResultsPointer, ResultsDistancesPointer, NumberOfResults[i], MaxNumberOfResults, Box );
    }
  }

  /// Contact search API

  /**
   * [SearchContact description]
   * NOTE[Charlie]: Why this function does not return the number of results like the others?
   * @param Result [description]
   */
  void SearchContact(ContainerContactType& Result) {
    for (CellContainerIterator icell = mCells.begin() ; icell!= mCells.end(); icell++) {
      icell->SearchContact(Result);
    }
  }

  /**
   * [SearchContact description]
   * @param  Result             [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  SizeType SearchContact(IteratorContactType& Result, const SizeType& MaxNumberOfResults ) {
    SizeType NumberOfResults = 0;
    for (CellContainerIterator icell = mCells.begin() ; icell!= mCells.end(); icell++) {
      icell->SearchContact(Result, NumberOfResults, MaxNumberOfResults);
    }
    return NumberOfResults;
  }

  /// Add/Remove

  /**
   * [AddObject description]
   * @param ThisObject [description]
   */
  virtual void AddObject(const PointerType& ThisObject) {
    PointType Low, High;
    SearchStructureType Box;

    TConfigure::CalculateBoundingBox(ThisObject, Low, High);
    Box.Set( CalculateCell(Low), CalculateCell(High), mN );
    FillObject(Box,ThisObject);

    mObjectsSize++;
  }

  /**
   * [RemoveObject description]
   * @param ThisObject [description]
   */
  void RemoveObject(const PointerType& ThisObject) {
    PointType Low, High;
    SearchStructureType Box;

    TConfigure::CalculateBoundingBox(ThisObject, Low, High);
    Box.Set( CalculateCell(Low), CalculateCell(High), mN );
    RemoveObjectLocal(Box,ThisObject);

    mObjectsSize--;
  }

  ///@}
  ///@name Operators
  ///@{


  ///@}
  ///@name Operations
  ///@{


  ///@}
  ///@name Access
  ///@{


  ///@}
  ///@name Inquiry
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  virtual std::string Info() const {
    return "BinsObjectDynamic" ;
  }

  /** Print information about this object.
   * Print information about this object.
   * @param rOStream [description]
   */
  virtual void PrintInfo(std::ostream& rOStream) const {
    rOStream << Info();
  }

  /** Print object's data.
   * Print object's data.
   * @param rOStream [description]
   * @param Perfix   [description]
   */
  virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const {
    rOStream << " BinsSize: ";
    for(SizeType i = 0 ; i < Dimension ; i++) {
        rOStream << "[" << mN[i] << "]";
    }
    rOStream << std::endl;
    rOStream << "  CellSize: ";
    for(SizeType i = 0 ; i < Dimension ; i++) {
        rOStream << "[" << mCellSize[i] << "]";
    }
    rOStream << std::endl;
    SizeType nn = 0;
    for(SizeType i = 0 ; i < mCells.size(); i++) {
        nn += mCells[i].Size();
    }
    rOStream << "NumPointers: " << nn << std::endl;
  }

  /** Print Size of Container
   * Print Size of Container
   * @param rout [description]
   */
  void PrintSize(std::ostream& rout) {
    rout << " BinsSize: ";
    for(SizeType i = 0 ; i < Dimension ; i++) {
      rout << "[" << mN[i] << "]";
    }
    rout << std::endl;
  }

  /** Print Limits Points of the Container
   * Print Limits Points of the Container
   * @param rout [description]
   */
  void PrintBox(std::ostream& rout) {
    rout << " BinsBox: Min [";
    mMinPoint.Print(rout);
    rout <<       "];  Max [";
    mMaxPoint.Print(rout);
    rout <<       "];  Size [";
    mCellSize.Print(rout);
    rout << "]" << std::endl;
  }

  /**
   * [GetCellContainer description]
   * @return [description]
   */
  CellContainerType& GetCellContainer() {
    return mCells;
  }

  /**
   * [GetDivisions description]
   * @return [description]
   */
  SizeArray& GetDivisions() {
    return mN;
  }

  /**
   * [GetCellSize description]
   * @return [description]
   */
  CoordinateArray& GetCellSize() {
    return mCellSize;
  }

  /**
   * [GetMinPoint description]
   * @return [description]
   */
  PointType& GetMinPoint() {
    return mMinPoint;
  }

  /**
   * [GetMaxPoint description]
   * @return [description]
   */
  PointType& GetMaxPoint() {
    return mMaxPoint;
  }

  /** Calculates the IndexArray (x[,y[,z]]) of the provided object.
   * Calculates the IndexArray (x[,y[,z]]) of the provided object.
   * The provided object must provide its coordinates through the [] operator.
   * @param  ThisObject Input Object
   * @return            Cell coordinates of 'ThisObject' in the bins
   */
  template<class GenericCoordType>
  IndexArray CalculateCell(const GenericCoordType& ThisObject) {
    IndexArray IndexCell;

    for(SizeType i = 0 ; i < Dimension ; i++) {
      IndexCell[i] = CalculatePosition(ThisObject[i],i);
    }

    return IndexCell;
  }

  /** Calculates the Index of the provided object.
   * Calculates the Index of the provided object.
   * The provided object must provide its coordinates through the [] operator.
   * @param  ThisObject Input Object
   * @return            Cell index of 'ThisObject' in the bins
   */
  template<class GenericCoordType>
  IndexType CalculateIndex(const GenericCoordType& ThisObject) {
    IndexType Index = 0;

    for(SizeType iDim = Dimension-1 ; iDim > 0 ; iDim--) {
      Index += CalculatePosition(ThisObject[iDim],iDim);
      Index *= mN[iDim-1];
    }

    Index += CalculatePosition(ThisObject[0],0);

    return Index;
  }

  /**
   * [CalculatePosition description]
   * @param  ThisCoord     [description]
   * @param  ThisDimension [description]
   * @return               [description]
   */
  virtual IndexType CalculatePosition(CoordinateType const& ThisCoord, const SizeType& ThisDimension) {
    CoordinateType d_index = (ThisCoord - mMinPoint[ThisDimension]) * mInvCellSize[ThisDimension];
    IndexType index = static_cast<IndexType>( (d_index < 0.00) ? 0.00 : d_index );

    return  (index > mN[ThisDimension]-1) ? mN[ThisDimension]-1 : index;
  }

protected:


    ///@}
    ///@name Friends
    ///@{


    ///@}

    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{

    /// It computes each object's boundinx box and uses it to find the max and min points
    virtual void CalculateBoundingBox()
    {
        PointType Low, High;
        TConfigure::CalculateBoundingBox(*mObjectsBegin,mMinPoint,mMaxPoint);

#ifdef _OPENMP
        SizeType number_of_threads = omp_get_max_threads();
#else
        SizeType number_of_threads = 1;
#endif

        std::vector<SizeType> node_partition;
        CreatePartition(number_of_threads, mObjectsSize, node_partition);

        std::vector<PointType> Max(number_of_threads);
        std::vector<PointType> Min(number_of_threads);

        for(SizeType k=0; k<number_of_threads; k++ )
        {
            Max[k] = mMaxPoint;
            Min[k] = mMinPoint;
        }

        IteratorType i_begin = mObjectsBegin;
        IteratorType i_end   = mObjectsEnd;

        for (IteratorType i_object = i_begin ; i_object != i_end ; i_object++ )
        {
            TConfigure::CalculateBoundingBox(*i_object, Low, High);
            for(SizeType i = 0 ; i < Dimension ; i++)
            {
                mMaxPoint[i] = (mMaxPoint[i] < High[i]) ? High[i] : mMaxPoint[i];
                mMinPoint[i] = (mMinPoint[i] > Low[i])  ? Low[i]  : mMinPoint[i];
            }
        }

        PointType Epsilon = mMaxPoint - mMinPoint;

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            mMaxPoint[i] += Epsilon[i] * 0.01;
            mMinPoint[i] -= Epsilon[i] * 0.01;
        }
    }


//************************************************************************
//************************************************************************

    virtual void CalculateCellSize()
    {

        CoordinateType delta[Dimension];
        CoordinateType alpha[Dimension];
        CoordinateType mult_delta = 1.00;
        SizeType index = 0;
        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            delta[i] = mMaxPoint[i] - mMinPoint[i];
            if ( delta[i] > delta[index] )
                index = i;
            delta[i] = (delta[i] == 0.00) ? 1.00 : delta[i];
        }

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            alpha[i] = delta[i] / delta[index];
            mult_delta *= alpha[i];
        }


        mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>(mObjectsSize/mult_delta), 1.00/Dimension) +1 );

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            if(i!=index)
            {
                mN[i] = static_cast<SizeType>(alpha[i] * mN[index]);
                mN[i] = ( mN[i] == 0 ) ? 1 : mN[i];
            }
        }

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            mCellSize[i] = delta[i] / mN[i];
            mInvCellSize[i] = 1.00 / mCellSize[i];
        }


    }

    virtual void CalculateCellSize(const CoordinateType& CellSize)
    {
        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            mCellSize[i] = CellSize;
            mInvCellSize[i] = 1.00 / mCellSize[i];
            mN[i] = static_cast<SizeType>( (mMaxPoint[i]-mMinPoint[i]) / mCellSize[i]) + 1;
        }
    }

    virtual void CalculateCellSize( const SizeType& NumPoints )
    {

        CoordinateType delta[Dimension];
        CoordinateType alpha[Dimension];
        CoordinateType mult_delta = 1.00;
        SizeType index = 0;
        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            delta[i] = mMaxPoint[i] - mMinPoint[i];
            if ( delta[i] > delta[index] )
                index = i;
            delta[i] = (delta[i] == 0.00) ? 1.00 : delta[i];
        }

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            alpha[i] = delta[i] / delta[index];
            mult_delta *= alpha[i];
        }

        mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>(NumPoints/mult_delta), 1.00/Dimension) +1 );

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            if(i!=index)
            {
                mN[i] = static_cast<SizeType>(alpha[i] * mN[index]);
                mN[i] = ( mN[i] == 0 ) ? 1 : mN[i];
            }
        }

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            mCellSize[i] = delta[i] / mN[i];
            mInvCellSize[i] = 1.00 / mCellSize[i];
        }

    }

//************************************************************************
//************************************************************************

    virtual void GenerateBins()
    {
        PointType Low, High;
        SearchStructureType Box;
        /// Fill container with objects
        for(IteratorType i_object = mObjectsBegin ; i_object != mObjectsEnd ; i_object++)
        {
            TConfigure::CalculateBoundingBox(*i_object, Low, High);
            Box.Set( CalculateCell(Low), CalculateCell(High), mN );
            FillObject(Box, *i_object);
        }
    }


//************************************************************************
//************************************************************************

// **** THREAD SAFE
// Dimension = 1
    void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;
        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0])
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                mCells[I].SearchObjects(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    mCells[I].SearchObjects(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
            }
        }
    }

    // Dimension = 3
    void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result,
                          SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    {
                        mCells[I].SearchObjects(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
                    }
                }
            }
        }
    }

//************************************************************************
//************************************************************************

    // **** THREAD SAFE

    // Dimension = 1
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Result,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;
        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
        {
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                mCells[I].SearchObjects(ThisObject, Result);
        }
    }

    // Dimension = 2
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Result,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];

        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    mCells[I].SearchObjects(ThisObject, Result);
            }
        }
    }

    // Dimension = 3
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Result,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];

        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    {
                        mCells[I].SearchObjects(ThisObject, Result);
                    }
                }
            }
        }
    }

//************************************************************************
//************************************************************************

// **** THREAD SAFE
// Dimension = 1
    void SearchObjectLocalExclusive(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                                SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;
        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0])
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                mCells[I].SearchObjectsExclusive(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchObjectLocalExclusive(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                                SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    mCells[I].SearchObjectsExclusive(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
            }
        }
    }

    // Dimension = 3
    void SearchObjectLocalExclusive(PointerType& ThisObject, ResultIteratorType& Result,
                                SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                                SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    {
                        mCells[I].SearchObjectsExclusive(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
                    }
                }
            }
        }
    }

//************************************************************************
//************************************************************************


    // **** THREAD SAFE

    // Dimension = 1
    void SearchObjectLocalExclusive(PointerType& ThisObject, ResultContainerType& Result,
                                SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;
        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                mCells[I].SearchObjectsExclusive(ThisObject, Result);

    }

    // Dimension = 2
    void SearchObjectLocalExclusive(PointerType& ThisObject, ResultContainerType& Result,
                                SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];

        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    mCells[I].SearchObjectsExclusive(ThisObject, Result);
            }
        }
    }

    // Dimension = 3
    void SearchObjectLocalExclusive(PointerType& ThisObject, ResultContainerType& Result,
                                SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];

        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
        {
            MinCell[2] = MinBox[2];
            MaxCell[2] = MaxBox[2];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    {
                        mCells[I].SearchObjectsExclusive(ThisObject, Result);
                    }
                }
            }
        }
    }

//************************************************************************
//************************************************************************


    // **** THREAD SAFE

    // Dimension = 1
    void SearchInRadius(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];

        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0])
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                mCells[I].SearchObjectsInRaius(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInRadius(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    mCells[I].SearchObjectsInRaius(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);
            }
        }
    }

    // Dimension = 3
    void SearchInRadius(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    {
                        mCells[I].SearchObjectsInRadius(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);
                    }
                }
            }
        }
    }

//************************************************************************
//************************************************************************


    // **** THREAD SAFE

    // Dimension = 1
    void SearchInRadius(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];

        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0])
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                mCells[I].SearchObjectsInRaius(ThisObject, Radius, Result, ResultDistances, NumberOfResults, MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInRadius(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    mCells[I].SearchObjectsInRaius(ThisObject, Radius, Result, ResultDistances, NumberOfResults, MaxNumberOfResults);
            }
        }
    }

    // Dimension = 3
    void SearchInRadius(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    {
                        mCells[I].SearchObjectsInRadius(ThisObject, Radius, Result, ResultDistances, NumberOfResults, MaxNumberOfResults);
                    }
                }
            }
        }
    }

//************************************************************************
//************************************************************************


    // **** THREAD SAFE

    // Dimension = 1
    virtual void SearchInRadiusExclusive(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];

        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0])
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                mCells[I].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);
    }

    // Dimension = 2
    virtual void SearchInRadiusExclusive(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    mCells[I].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);
            }
        }
    }

    // Dimension = 3
    virtual void SearchInRadiusExclusive(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    {
                        mCells[I].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);
                    }
                }
            }
        }
    }

//************************************************************************
//************************************************************************


    // **** THREAD SAFE

    // Dimension = 1
    virtual void SearchInRadiusExclusive(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];

        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0])
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                mCells[I].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, ResultDistances, NumberOfResults, MaxNumberOfResults);
    }

    // Dimension = 2
    virtual void SearchInRadiusExclusive(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    mCells[I].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, ResultDistances, NumberOfResults, MaxNumberOfResults);
            }
        }
    }

    // Dimension = 3
    virtual void SearchInRadiusExclusive(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    {
                        mCells[I].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, ResultDistances, NumberOfResults, MaxNumberOfResults);
                    }
                }
            }
        }
    }

//************************************************************************
//************************************************************************

    // Dimension = 1
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, const PointerType& i_object)
    {
        PointType  MinCell, MaxCell;

        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
        {
            if(TConfigure::IntersectionBox(i_object, MinCell, MaxCell))
                mCells[I].Add(i_object);
        }
    }


//************************************************************************
//************************************************************************

    // Dimension = 2
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, const PointerType& i_object)
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
                    mCells[I].Add(i_object);
            }
        }
    }


//************************************************************************
//************************************************************************

    // Dimension = 3
    virtual void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, const PointerType& i_object)
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2]+=mCellSize[2], MaxCell[2]+=mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
                        mCells[I].Add(i_object);
                }
            }
        }
    }

//************************************************************************
//************************************************************************

    // Dimension = 1
    void RemoveObjectLocal( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, const PointerType& i_object)
    {
        PointType  MinCell, MaxCell;

        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
        {
            if(TConfigure::IntersectionBox(i_object, MinCell, MaxCell))
                mCells[I].Remove(i_object);
        }
    }


//************************************************************************
//************************************************************************

    // Dimension = 2
    void RemoveObjectLocal( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, const PointerType& i_object)
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
                    mCells[I].Remove(i_object);
            }
        }
    }


//************************************************************************
//************************************************************************

    // Dimension = 3
    void RemoveObjectLocal( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, const PointerType& i_object)
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2]+=mCellSize[2], MaxCell[2]+=mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
                        mCells[I].Remove(i_object);
                }
            }
        }
    }


//************************************************************************
//************************************************************************

    void AllocateContainer()
    {
        SizeType Size = mN[0];
        for(SizeType i = 1 ; i < Dimension ; i++)
            Size *= mN[i];
        mCells.resize(Size);
    }

    inline void CreatePartition(SizeType number_of_threads, const SizeType number_of_rows, std::vector<SizeType>& partitions)
    {
        partitions.resize(number_of_threads+1);
        SizeType partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for(SizeType i = 1; i<number_of_threads; i++)
            partitions[i] = partitions[i-1] + partition_size ;
    }

    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    PointType    mMinPoint;
    PointType    mMaxPoint;

    SizeType     mObjectsSize;
    IteratorType mObjectsBegin;
    IteratorType mObjectsEnd;

    CoordinateArray  mCellSize;
    CoordinateArray  mInvCellSize;
    SizeArray        mN;

    CellContainerType mCells;  ///The bin

private:

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
    BinsObjectDynamic<TConfigure> & operator=(const BinsObjectDynamic<TConfigure> & rOther)
    {
        mMinPoint            = rOther.mMinPoint;
        mMaxPoint            = rOther.mMaxPoint;
        mObjectsBegin        = rOther.mObjectsBegin;
        mObjectsEnd          = rOther.mObjectsEnd;
        mObjectsSize         = rOther.mObjectsSize;
        mCellSize            = rOther.mCellSize;
        mInvCellSize         = rOther.mInvCellSize;
        mN                   = rOther.mN;
        mCells               = rOther.mCells;
        return *this;
    }

    /// Copy constructor.
    BinsObjectDynamic(const BinsObjectDynamic& rOther)
    {
        *this =  rOther;
    }

    /// Copy constructor.
    template<class T>
    BinsObjectDynamic(const BinsObjectDynamic<T>& rOther)
    {
        *this =  rOther;
    }

}; // Class BinsObjectDynamic

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TConfigure>
inline std::istream& operator >> (std::istream& rIStream,
                                  BinsObjectDynamic<TConfigure>& rThis)
{
    return rIStream;
}


/// output stream function
template<class TConfigure>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BinsObjectDynamic<TConfigure> & rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
