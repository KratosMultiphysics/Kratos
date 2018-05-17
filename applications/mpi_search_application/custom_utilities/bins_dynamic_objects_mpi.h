//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

#if !defined(KRATOS_BINS_DYNAMIC_OBJECTS_MPI_CONTAINER_H_INCLUDED)
#define  KRATOS_BINS_DYNAMIC_OBJECTS_MPI_CONTAINER_H_INCLUDED

// System includes
#include <unordered_set>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <time.h>

// OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// Mpi
#include "mpi.h"

// Project includes
#include "includes/define.h"
#include "includes/mpi_communicator.h"

#include "spatial_containers/tree.h"
#include "spatial_containers/cell.h"

#include "custom_configures/partition_configure.h"
#include "spatial_containers/bins_dynamic_objects.h"

// ACTIVATES AND DISABLES TIMER
#define CUSTOMTIMER 1

// Timer defines
#ifdef CUSTOMTIMER
  #define KRATOS_TIMER_START(t) Timer::Start(t);
  #define KRATOS_TIMER_STOP(t)  Timer::Stop(t);
#else
  #define KRATOS_TIMER_START(t)
  #define KRATOS_TIMER_STOP(t)
#endif

namespace Kratos {

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

/**
 * Mpi version of the Bins.
 */
template<class TObjectConfigure, class TPartConfigure = PartitionConfigure<TObjectConfigure>>
class BinsObjectDynamicMpi {

public:
  ///@name Type Definitions
  ///@{

  enum { Dimension = TPartConfigure::Dimension };

  /// Point
  typedef typename TObjectConfigure::PointType            PointType;

  /// Container
  typedef typename TObjectConfigure::ObjectType           ObjectType;
  typedef typename TObjectConfigure::ObjectContainerType  ObjectContainerType;
  typedef typename TObjectConfigure::PointerType          PointerType;
  typedef typename TObjectConfigure::ContainerType        ContainerType;
  typedef typename TObjectConfigure::IteratorType         IteratorType;
  typedef typename TObjectConfigure::DistanceIteratorType DistanceIteratorType;
  typedef typename TObjectConfigure::ResultContainerType  ResultContainerType;
  typedef typename TObjectConfigure::ResultIteratorType   ResultIteratorType;

  /// Search Structures
  typedef typename TPartConfigure::PointType              PartPointType;
  typedef typename TPartConfigure::ObjectType             PartObjectType;
  typedef typename TPartConfigure::PointerType            PartPtrType;

  /// Note: Please note that cell store TPartConfigure
  typedef Cell<TPartConfigure>                            CellType;
  typedef std::vector<CellType>                          CellContainerType;
  typedef std::vector<int>                               DomainContainerType;
  typedef typename CellContainerType::iterator           CellContainerIterator;

  typedef TreeNode<
    Dimension, PointType, PointerType, IteratorType,
    typename TPartConfigure::DistanceIteratorType>        TreeNodeType;
  typedef typename TreeNodeType::CoordinateType          CoordinateType;  // double
  typedef typename TreeNodeType::SizeType                SizeType;        // std::size_t
  typedef typename TreeNodeType::IndexType               IndexType;       // std::size_t

  typedef Tvector<IndexType,Dimension>                   IndexArray;
  typedef Tvector<SizeType,Dimension>                    SizeArray;
  typedef Tvector<CoordinateType,Dimension>              CoordinateArray;

  /// Contact Pair
  typedef typename TPartConfigure::ContainerContactType   ContainerContactType;
  typedef typename TPartConfigure::IteratorContactType    IteratorContactType;

  /// typedef TreeNodeType LeafType;
  typedef typename TreeNodeType::IteratorIteratorType    IteratorIteratorType;
  typedef typename TreeNodeType::SearchStructureType     SearchStructureType;

  /// Partition Bins configuration files
  typedef typename TPartConfigure::ResultContainerType    PartitionResultContainerType;

  /// Pointer definition of BinsObjectDynamic
  KRATOS_CLASS_POINTER_DEFINITION(BinsObjectDynamicMpi);

  ///@}
  ///@name Life Cycle
  ///@{

  /** Constructor using the objects.
   * Creates a bins by calculating the cell size based on the objects
   * from 'ObjectsBegin' to 'ObjectsEnd' and fills it with such objects.
   * @param ObjectsBegin Iterator pointing to the first object to add to the bins
   * @param ObjectsEnd   Iterator pointing to the last object to add to the bins
   */
  BinsObjectDynamicMpi(IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd) {
    if(TPartConfigure::Dimension != TObjectConfigure::Dimension) {
      KRATOS_ERROR << "'Partition' and 'Object' configuration must have the same dimension" << std::endl;
    }
    InitializeMpi();
    GenerateLocalBins(ObjectsBegin, ObjectsEnd);
    GeneratePartitionBins(ObjectsBegin, ObjectsEnd);
  }

  /** Constructor using the objects with a fixed cell size.
   * Creates a bins by using the cell size provided in 'CellSize' and filling
   * it with the objects from 'ObjectsBegin' to 'ObjectsEnd'.
   * @param ObjectsBegin Iterator pointing to the first object to add to the bins
   * @param ObjectsEnd   Iterator pointing to the last object to add to the bins
   * @param CellSize     Size of the cells (x,y,z) for the bins
   */
  BinsObjectDynamicMpi(IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd, CoordinateType CellSize) {
    if(TPartConfigure::Dimension != TObjectConfigure::Dimension) {
      KRATOS_ERROR << "'Partition' and 'Object' configuration must have the same dimension" << std::endl;
    }
    KRATOS_ERROR << "Not yet implemented" << std::endl;
  }

  /** Constructor using the bounding box with a fixed cell size.
   * Creates a bins using the bounding descrived by 'MinPoint' and 'MaxPoint'
   * with cells of size 'CellSize'.
   * @param ObjectsBegin Iterator pointing to the first object to add to the bins
   * @param ObjectsEnd   Iterator pointing to the last object to add to the bins
   * @param CellSize     Size of the cells (x,y,z) for the bins
   */
  BinsObjectDynamicMpi(const PointType& MinPoint, const PointType& MaxPoint, CoordinateType CellSize) {
    if(TPartConfigure::Dimension != TObjectConfigure::Dimension) {
      KRATOS_ERROR << "'Partition' and 'Object' configuration must have the same dimension" << std::endl;
    }
    KRATOS_ERROR << "Not yet implemented" << std::endl;
  }

  /** Constructor using the bounding box with a fixed cell size.
   * Creates a bins using the bounding descrived by 'MinPoint' and 'MaxPoint'
   * with cells of size calculated from the number of poitns 'NumPoints'.
   * @param ObjectsBegin Iterator pointing to the first object to add to the bins
   * @param ObjectsEnd   Iterator pointing to the last object to add to the bins
   * @param CellSize     Size of the cells (x,y,z) for the bins
   */
  BinsObjectDynamicMpi(const PointType& MinPoint, const PointType& MaxPoint, SizeType NumPoints) {
    if(TPartConfigure::Dimension != TObjectConfigure::Dimension) {
      KRATOS_ERROR << "'Partition' and 'Object' configuration must have the same dimension" << std::endl;
    }
    KRATOS_ERROR << "Not yet implemented" << std::endl;
  }

  /** Assignment operator.
   * Assignment operator.
   * @param rOther reference object
   */
  BinsObjectDynamicMpi<TObjectConfigure, TPartConfigure> & operator=(const BinsObjectDynamicMpi<TObjectConfigure, TPartConfigure> & rOther) {
    mpObjectBins    = rOther.mpObjectBins;

    return *this;
  }

  /** Copy constructor.
   * Copy constructor.
   * @param rOther reference object
   */
  BinsObjectDynamicMpi(const BinsObjectDynamicMpi& rOther) {
    *this = rOther;
  }

  /** Default constructor
   * Default constructor
   */
  virtual ~BinsObjectDynamicMpi() {
    // TODO: ObjectBins needs to be deleted!!!!!
    delete mpObjectBins;
  }


  /**
   * [SearchPartitionInRadius description]
   * @param Box           [description]
   * @param i_object      [description]
   * @param partitionList [description]
   * @param Radius        [description]
   */
  void SearchPartition(SearchStructureType & Box, const PointerType& i_object, std::unordered_set<std::size_t> & partitionList) {
    GetPartition(Box, i_object, partitionList);
  }

  /**
   * [SearchPartitionInRadius description]
   * @param Box           [description]
   * @param i_object      [description]
   * @param partitionList [description]
   * @param Radius        [description]
   */
  void SearchPartitionInRadius(SearchStructureType & Box, const PointerType& i_object, std::unordered_set<std::size_t> & partitionList, double Radius) {
    GetPartitionInRadius(Box, i_object, partitionList, Radius);
  }

  /** Searches objects in the cell.
   * Searches the objects in the same cell as 'ThisPoint'. the behaviour of this operation is highly variable
   * as a Point may be in cells of different size in different partitions.
   * This operation is not formally defined for distributed scenarios.
   * # IMPORTANT
   * At this point, this operation is LOCAL and does not provide interaction with the other processes
   * @param  ThisPoint          A Point in the space
   * @param  Result             Pointer to the list of results containing all objects up to 'MaxNumberOfResults'
   *                            that land in the same cell has 'Point'
   * @param  MaxNumberOfResults Maximum number of results.
   * @return                    The number of results found (0 to MaxNumberOfResults) or -1 in case of error.
   */
  SizeType SearchObjectsInCell(const PointType& ThisPoint, ResultIteratorType Result, const SizeType& MaxNumberOfResults) {
    KRATOS_ERROR << "Not yet implemented" << std::endl;
    return 0;
  }


  /** Searches objects in the same boundingbox
   * Searches up to 'MaxNumberOfResults' objects that intersect the boundingbox of 'ThisObject'
   * @param  ThisObject         Input object that determines the boundingbox to search from.
   * @param  Result             List of objects that intersect the boundingbox of 'ThisObject' up to 'MaxNumberOfResults'
   * @param  MaxNumberOfResults Maximum number of results.
   * @return                    The number of results found (0 to MaxNumberOfResults) or -1 in case of error.
   */
  SizeType SearchObjects(PointerType& ThisObject, ResultIteratorType& Result, const SizeType& MaxNumberOfResults, Communicator& Communicator) {
    KRATOS_ERROR << "Not yet implemented" << std::endl;
    return 0;
  }

  /**
   * [SearchObjectsInRadius description]
   * @param  ThisObject         [description]
   * @param  Radius             [description]
   * @param  Results            [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  void SearchObjectsInRadius(
      const IteratorType & ObjectsBegin,
      const IteratorType & ObjectsEnd,
      const SizeType & NumberOfObjects,
      std::vector<double> const & Radius,
      std::vector<std::vector<PointerType>> & Results,
      std::vector<std::vector<double>> & ResultsDistances,
      std::vector<SizeType>& NumberOfResults,
      SizeType const & MaxNumberOfResults,
      Communicator * MpiCommunicator = nullptr) {

    assert(NumberOfObjects > 0);
    assert(ObjectsEnd-ObjectsBegin == NumberOfObjects);

    std::vector<ObjectContainerType> RemoteResults(mpi_size);
    std::vector<ObjectContainerType> RemoteObjects(mpi_size);
    std::vector<ObjectContainerType> SearchResults(mpi_size);
    std::vector<ObjectContainerType> SendObjectToProcess(mpi_size);

    std::vector<std::vector<double>> SendRadiusToProcess(mpi_size, std::vector<double>(0));
    std::vector<std::vector<double>> RemoteRadius(mpi_size, std::vector<double>(0));
    std::vector<std::vector<double>> RemoteResultsSize(mpi_size, std::vector<double>(0));
    std::vector<std::vector<double>> RecvResultsPerObject(mpi_size, std::vector<double>(0));

    std::vector<bool> SentObjectsMap(NumberOfObjects*mpi_size, 0);

    // Local Search
    SearchInRadiusLocal(
      ObjectsBegin, ObjectsEnd, NumberOfObjects, Radius,
      Results, ResultsDistances, NumberOfResults,
      SendObjectToProcess, SendRadiusToProcess, SentObjectsMap,
      MaxNumberOfResults
    );

    // TODO: Transfer only the points
    TObjectConfigure::TransferObjects(MpiCommunicator, SendObjectToProcess, RemoteObjects);
    TObjectConfigure::TransferObjects(SendRadiusToProcess, RemoteRadius);

    // Remote search
    SearchInRadiusRemote(
      RemoteObjects, RemoteRadius,
      RemoteResults, RemoteResultsSize,
      MaxNumberOfResults,
      MpiCommunicator
    );

    // Here we transfer the objects that will be ghost (they have been found in a search with center "canditate to ghost for the other partition")
    TObjectConfigure::TransferObjects(MpiCommunicator, RemoteResults, SearchResults);
    TObjectConfigure::TransferObjects(RemoteResultsSize, RecvResultsPerObject);

    AssembleLocalAndRemoteResults(
      ObjectsBegin, ObjectsEnd, NumberOfObjects,
      Results, ResultsDistances, NumberOfResults,
      SentObjectsMap, SearchResults, RecvResultsPerObject,
      MpiCommunicator
    );

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
    return "BinsObjectDynamicMpi" ;
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const {
      rOStream << Info();
  }

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const {
    // rOStream << " BinsSize: ";
    // for(SizeType i = 0 ; i < Dimension ; i++) {
    //     rOStream << "[" << mN[i] << "]";
    // }
    // rOStream << std::endl;
    // rOStream << "  CellSize: ";
    // for(SizeType i = 0 ; i < Dimension ; i++) {
    //     rOStream << "[" << mCellSize[i] << "]";
    // }
    // rOStream << std::endl;
    // SizeType nn = 0;
    // for(SizeType i = 0 ; i < mCells.size(); i++) {
    //     nn += mCells[i].Size();
    // }
    // rOStream << "NumPointers: " << nn << std::endl;
  }

  /// Print Size of Container
  void PrintSize( std::ostream& rout ) {
    // rout << " BinsSize: ";
    // for(SizeType i = 0 ; i < Dimension ; i++) {
    //   rout << "[" << mN[i] << "]";
    // }
    // rout << std::endl;
  }

  /// Print Limits Points of the Container
  void PrintBox( std::ostream& rout ) {
    // rout << " BinsBox: Min [";
    // mMinPoint.Print(rout);
    // rout <<       "];  Max [";
    // mMaxPoint.Print(rout);
    // rout <<       "];  Size [";
    // mCellSize.Print(rout);
    // rout << "]" << std::endl;
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

  /** GetLocalCellContainer
   * GetLocalCellContainer
   * @return [description]
   */
  CellContainerType& GetLocalCellContainer() {
    return mpObjectBins->GetCellContainer();
  }

  /** GetLocalDivisions
   * GetLocalDivisions
   * @return [description]
   */
  SizeArray& GetLocalDivisions() {
    return mpObjectBins->GetDivisions();
  }

  /** GetLocalCellSize
   * GetLocalCellSize
   * @return [description]
   */
  CoordinateArray& GetLocalCellSize() {
    return mpObjectBins->GetCellSize();
  }

  /** GetLocalMinPoint
   * GetLocalMinPoint
   * @return [description]
   */
  PointType& GetLocalMinPoint() {
    return mpObjectBins->GetMinPoint();
  }

  /** GetLocalMaxPoint
   * GetLocalMaxPoint
   * @return [description]
   */
  PointType& GetLocalMaxPoint() {
    return mpObjectBins->GetMaxPoint();
  }

  /** Calculates the IndexArray (x[,y[,z]]) of the provided object.
   * Calculates the IndexArray (x[,y[,z]]) of the provided object.
   * The provided object must provide its coordinates through the [] operator.
   * @param  ThisObject Input Object
   * @return            Cell coordinates of 'ThisObject' in the bins
   */
  IndexArray CalculateCell(const PointType& ThisObject) {
    IndexArray IndexCell;

    for(SizeType i = 0 ; i < Dimension ; i++) {
      IndexCell[i] = CalculatePosition(ThisObject[i],i);
    }

    return IndexCell;
  }

  /*
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

  /** Calculates the position of a coordinate
   * Calculates the position of a coordinate inside the bins
   * @param  ThisCoord     [description]
   * @param  ThisDimension [description]
   * @return               [description]
   */
  virtual IndexType CalculatePosition(CoordinateType const& ThisCoord, const SizeType& ThisDimension) {
    CoordinateType d_index = (ThisCoord - mMinPoint[ThisDimension]) * mInvCellSize[ThisDimension];
    IndexType index = static_cast<IndexType>( (d_index < 0.00) ? 0.00 : d_index );

    return  (index > mN[ThisDimension]-1) ? mN[ThisDimension]-1 : index;
  }


private:

  /**
   * Default constructor its private and should never be called
   */
  BinsObjectDynamicMpi() {}

  double ReduceMaxRadius(IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd) {
    // Max Radius Ugly fix
    double local_max_radius = 0.0f;
    double max_radius = 0.0f;
    for (IteratorType ObjectItr = ObjectsBegin; ObjectItr != ObjectsEnd; ObjectItr++) {
      const double Radius = TObjectConfigure::GetObjectRadius(*ObjectItr, 0.0f);
      if(Radius > local_max_radius) local_max_radius = Radius;
    }

    MPI_Allreduce(&local_max_radius, &max_radius, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    return max_radius;
  }

  /**
   * Initializes mpi_rank and mpi_size for the given mpi communicator.
   * If no mpi communicator is provided MPI_COMM_WORLD is used as default.
   */
  void InitializeMpi(MPI_Comm mpiComm = MPI_COMM_WORLD) {
    MPI_Comm_rank(mpiComm, &mpi_rank);
    MPI_Comm_size(mpiComm, &mpi_size);
  }

  /**
   * Generate the bins for the local objects
   * @param ObjectsBegin Iterator to the first local object
   * @param ObjectsEnd   Iterator to the last local object
   */
  void GenerateLocalBins(IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd) {
    mpObjectBins = new BinsObjectDynamic<TObjectConfigure>(ObjectsBegin, ObjectsEnd);
  }

  /**
   * [CalculateCustomCellSize description]
   * @param magic [description]
   */
  void CalculateCustomCellSize(const std::vector<int> & magic) {
    double totalDiff[3];
    double max = 0;
    //int maxI = 0;
    //double size;

    for(SizeType i = 0; i < Dimension; i++) {
      totalDiff[i] = this->mMaxPoint[i]-this->mMinPoint[i];
      if(totalDiff[i] > max) {
        max = totalDiff[i];
        //maxI = i;
      }
    }

    //size = max / magic;

    for(SizeType i = 0; i < Dimension; i++) {
      this->mN[i] = magic[i];
      this->mCellSize[i] = totalDiff[i] / this->mN[i];
      this->mInvCellSize[i] = 1.00 / this->mCellSize[i];
    }
  }

  /**
   * Generate the bins for the partitions
   * @param ObjectsBegin Iterator to the first local object
   * @param ObjectsEnd   Iterator to the last local object
   * @param mpiComm      Communicator to use ( default MPI_COMM_WORLD )
   */
  void GeneratePartitionBins(IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd, MPI_Comm mpiComm = MPI_COMM_WORLD) {

    // Calculate the boundingBox Of this bins
    std::vector<double> localMinPoint(Dimension);
    std::vector<double> localMaxPoint(Dimension);

    std::vector<double> globalMinPoint(Dimension);
    std::vector<double> globalMaxPoint(Dimension);

    //auto & binsObjectMinPoint = mpObjectBins->GetMinPoint();
    //auto & binsObjectMaxPoint = mpObjectBins->GetMaxPoint();

    //double maxRadius = ReduceMaxRadius(ObjectsBegin, ObjectsEnd);

    for(int d = 0; d < Dimension; d++) {
      localMinPoint[d] = (*ObjectsBegin)->GetGeometry()[0][d];
      localMaxPoint[d] = (*ObjectsBegin)->GetGeometry()[0][d];
    }

    for(IteratorType objectItr = ObjectsBegin; objectItr != ObjectsEnd; objectItr++) {
      for(int d = 0; d < Dimension; d++) {
        localMinPoint[d] = localMinPoint[d] < (*objectItr)->GetGeometry()[0][d] ? localMinPoint[d] : (*objectItr)->GetGeometry()[0][d];
        localMaxPoint[d] = localMaxPoint[d] > (*objectItr)->GetGeometry()[0][d] ? localMaxPoint[d] : (*objectItr)->GetGeometry()[0][d];
      }
    }

    // for(std::size_t i = 0; i < Dimension; i++) {
    //   localMinPoint[i] = binsObjectMinPoint[i];
    //   localMaxPoint[i] = binsObjectMaxPoint[i];
    // }

    // Share the boundingbox of all bins
    MPI_Allreduce(&localMinPoint[0], &globalMinPoint[0], Dimension, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&localMaxPoint[0], &globalMaxPoint[0], Dimension, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    for(std::size_t i = 0; i < Dimension; i++) {
      this->mMinPoint[i] = globalMinPoint[i];
      this->mMaxPoint[i] = globalMaxPoint[i];
    }

    int findex = 1;
    std::vector<int> numberOfCells = {2, 2, 2};
    while(numberOfCells[0]*numberOfCells[1]*numberOfCells[2] != mpi_size) {
      numberOfCells[findex++] *= 2;
      findex %= 3;
    }

    this->CalculateCustomCellSize(numberOfCells);
    this->AllocateContainer();

    // Initialize the partition bins
    std::size_t numCells = 1;
    for(std::size_t i = 0; i < Dimension; i++) {
      numCells *= this->mN[i];
    }

    for(std::size_t i = 0; i < numCells; i++) {
      PartPtrType cellPartitions = PartPtrType(new PartObjectType());
      this->mCells[i].Add(cellPartitions);
    }

    // For every object in this bins, mark the local and global vectors of the partition bins in the
    // cell from empty (0) to occupied (1)
    PointType Low, High, Center;
    SearchStructureType Box;

    for(IteratorType objectItr = ObjectsBegin; objectItr != ObjectsEnd; objectItr++) {
      TPartConfigure::CalculateCenter(*objectItr, Center);
      Box.Set( this->CalculateCell(Center), this->CalculateCell(Center), this->mN );

      FillPartition(Box, *objectItr);
    }

    SynchronaizePartitionsInCell();
  }

  /** Synchronizes the partitions in every cell across all processes
   * Synchronizes the partitions in every cell across all processes. For every cell:
   * 1. Get the maximum number of partitions for each process in the cell
   * 2. Merge the local partitions with the remote partitions from every process
   */
  void SynchronaizePartitionsInCell() {
    for(std::size_t cellId = 0; cellId < this->mCells.size(); cellId++) {
      auto & cellPartitions = (*(this->mCells[cellId].GetObject(0)))();

      // Calculate the number of neighbours that will be sent to the process "p" from process "q"
      int sendSize = cellPartitions.size();
      int recvMaxSize = 0;
      std::vector<int> recvSize(mpi_size, 0);

      MPI_Allreduce(&sendSize, &recvMaxSize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allgather(&sendSize, 1, MPI_INT, &recvSize[0], 1, MPI_INT, MPI_COMM_WORLD);

      std::vector<int> sendBuffer(cellPartitions.begin(), cellPartitions.end());
      std::vector<int> recvBuffer(mpi_size*recvMaxSize);

      if(sendBuffer.size() > (unsigned)recvMaxSize) {
        KRATOS_ERROR << "Something very impossible happened: Trying to send " << sendBuffer.size() << " element from " << recvMaxSize << " maximum number of elements." << std::endl;
      }

      sendBuffer.resize(recvMaxSize);
      MPI_Allgather(&sendBuffer[0], recvMaxSize, MPI_INT, &recvBuffer[0], recvMaxSize, MPI_INT, MPI_COMM_WORLD);

      for(int i = 0; i < mpi_size; i++) {
        if(i != mpi_rank) {
          for(int j = 0; j < recvSize[i]; j++) {
            cellPartitions.insert(recvBuffer[i*recvMaxSize+j]);
          }
        }
      }
    }

    auto meanPartsPerCell = 0.0f;
    auto CellsWithOne = 0.0f;
    for(std::size_t i = 0; i < this->mCells.size(); i++) {
      meanPartsPerCell += (*(this->mCells[i].GetObject(0)))().size();
      if((*(this->mCells[i].GetObject(0)))().size() == 1) {
        CellsWithOne++;
      }
    }

    meanPartsPerCell /= this->mCells.size();
  }

  /** SearchInRadiusLocal
   * SearchInRadiusLocal
   * @param  ThisObject         [description]
   * @param  Radius             [description]
   * @param  Results            [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  void SearchInRadiusLocal(
      const IteratorType & ObjectsBegin,
      const IteratorType & ObjectsEnd,
      const SizeType & NumberOfObjects,
      std::vector<double> const & Radius,
      std::vector<std::vector<PointerType>> & Results,
      std::vector<std::vector<double>> & ResultsDistances,
      std::vector<SizeType> & NumberOfResults,
      std::vector<ObjectContainerType> & SendObjectToProcess,
      std::vector<std::vector<double>> & SendRadiusToProcess,
      std::vector<bool> & SentObjectsMap,
      SizeType const & MaxNumberOfResults) {

    // Local Search
    PointType Low, High;
    SearchStructureType Box;

    for(int i = 0; i < mpi_size; i++) {
      SendObjectToProcess[i].clear();
      SendRadiusToProcess[i].clear();
    }

    int TotalToSend = 0;
    int GlobalToSend = 0;

    int TotalNumObjects = NumberOfObjects;
    int GlobalNumObjects = 0;

    double maxRadius = ReduceMaxRadius(ObjectsBegin, ObjectsEnd);

    for(std::size_t i = 0; i < NumberOfObjects; i++) {
      auto ObjectItr = ObjectsBegin + i;
      auto extensionRadius = i < Radius.size() ? Radius[i] : 0.0f;
      auto ObjectRadius = TObjectConfigure::GetObjectRadius(*ObjectItr, extensionRadius);

      ResultIteratorType   ResultsPointer          = Results[i].begin();
      DistanceIteratorType ResultsDistancesPointer = ResultsDistances[i].begin();

      // Search the results
      NumberOfResults[i] = mpObjectBins->SearchObjectsInRadiusExclusive(
        *ObjectItr, ObjectRadius,
        ResultsPointer, ResultsDistancesPointer,
        MaxNumberOfResults
      );

      // Search the partitions
      TObjectConfigure::CalculateBoundingBox(*ObjectItr, Low, High);
      for(int i = 0; i < Dimension; i++) {
        Low[i] -= maxRadius;
        High[i] += maxRadius;
      }

      Box.Set( this->CalculateCell(Low), this->CalculateCell(High), this->mN);

      std::unordered_set<std::size_t> partitionList;
      this->GetPartitionInRadius(Box, *ObjectItr, partitionList, ObjectRadius);
      // this->GetPartitionInRadius(Box, *ObjectItr, partitionList, Radius[i]);

      // For each point with results < MaxResults and each process excluding ourself
      for(int j = 0; j < mpi_size; j++) {
        if(j != mpi_rank && partitionList.find(j) != partitionList.end() && (NumberOfResults[i] < MaxNumberOfResults)) {
          SendObjectToProcess[j].push_back(*ObjectItr);
          SendRadiusToProcess[j].push_back(ObjectRadius);

          SentObjectsMap[j*NumberOfObjects+i]=1;
          TotalToSend++;
        }
      }
    }

    MPI_Reduce(&TotalToSend, &GlobalToSend, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&TotalNumObjects, &GlobalNumObjects, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(mpi_rank == 0) {
      std::cout << "Sending " << GlobalToSend << " unic objects out of " << GlobalNumObjects << std::endl;
    }
  }

  /** SearchInRadiusRemote
   * SearchInRadiusRemote
   * @param  ThisObject         [description]
   * @param  Radius             [description]
   * @param  Results            [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  void SearchInRadiusRemote(
      std::vector<ObjectContainerType> & RemoteObjects,
      std::vector<std::vector<double>> & RemoteRadius,
      std::vector<ObjectContainerType> & RemoteResults,
      std::vector<std::vector<double>> & RemoteResultsSize,
      const SizeType & MaxNumberOfResults,
      Communicator * MpiCommunicator) {

    std::vector<int> RecvObjectSize(mpi_size, 0);

    int TotalToSend = 0;
    int GlobalToSend = 0;

    // Calculate remote points
    for(int i = 0; i < mpi_size; i++) {
      int destination = -1;

      if(i != mpi_rank) {
        if(MpiCommunicator) {
          Communicator::NeighbourIndicesContainerType communicator_ranks = MpiCommunicator->NeighbourIndices();
          int NumberOfRanks = MpiCommunicator->GetNumberOfColors();

          for(int j = 0; j < NumberOfRanks; j++) {
            if(i == communicator_ranks[j]) {
              destination = j;
            }
          }
        }

        if(destination == -1 && RemoteObjects[i].size() != 0) {
          KRATOS_WATCH(MpiCommunicator->NeighbourIndices());
          KRATOS_ERROR << "Process " << mpi_rank << " trying to process results from process " << i << " but can't find neighbour index" << std::endl;
        }

        if(destination != -1 ) {

          RecvObjectSize[i] = RemoteObjects[i].size();
          RemoteResults[i].reserve((MaxNumberOfResults+1)*RecvObjectSize[i]);
          RemoteResultsSize[i].resize(RemoteRadius[i].size());

          std::vector<double> TempResultsDistances(MaxNumberOfResults);

          for(int j = 0; j < RecvObjectSize[i]; j++) {

            DistanceIteratorType RemoteDistancesPointer = TempResultsDistances.begin(); //Useless
            ResultContainerType  TempResults(MaxNumberOfResults);
            ResultIteratorType   RemoteResultsPointer = TempResults.begin();

            *RemoteDistancesPointer = 0;

            // Remote Search
            RemoteResultsSize[i][j] = mpObjectBins->SearchObjectsInRadiusExclusive(
              RemoteObjects[i].GetContainer()[j], RemoteRadius[i][j],
              RemoteResultsPointer, RemoteDistancesPointer,
              MaxNumberOfResults
            );

            TotalToSend += RemoteResultsSize[i][j];

            for(ResultIteratorType result_it = TempResults.begin(); result_it != RemoteResultsPointer; ++result_it) {
              TObjectConfigure::UpdateLocalInterface(MpiCommunicator, *result_it, destination);
              RemoteResults[i].push_back(*result_it);
            }
          }
        }
      }
    }

    MPI_Reduce(&TotalToSend, &GlobalToSend, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(mpi_rank == 0) {
      std::cout << "Found total of " << GlobalToSend << " non-unique objects among all partitions" << std::endl;
    }
  }

  /** AssembleLocalAndRemoteResults
   * AssembleLocalAndRemoteResults
   * @param  ThisObject         [description]
   * @param  Radius             [description]
   * @param  Results            [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  void AssembleLocalAndRemoteResults(
      const IteratorType & ObjectsBegin,
      const IteratorType & ObjectsEnd,
      const SizeType & NumberOfObjects,
      std::vector<std::vector<PointerType>> & Results,
      std::vector<std::vector<double>> & ResultsDistances,
      std::vector<SizeType> & NumberOfResults,
      const std::vector<bool> & SentObjectsMap,
      std::vector<ObjectContainerType> & SearchResults,
      const std::vector<std::vector<double>> & RecvResultsPerObject,
      Communicator * MpiCommunicator) {

    for(int i = 0; i < mpi_size; i++) {
      if(i != mpi_rank) {
        int ParticleCounter = 0;
        int ResultCounter = 0;

        int source = -1;

        if(MpiCommunicator) {
          Communicator::NeighbourIndicesContainerType communicator_ranks = MpiCommunicator->NeighbourIndices();
          int NumberOfRanks = MpiCommunicator->GetNumberOfColors();
          for(int j = 0; j < NumberOfRanks; j++) {
            if(i == communicator_ranks[j]) {
              source = j;
            }
          }
        }

        for(size_t j = 0; j < NumberOfObjects; j++) {
          if(SentObjectsMap[i*NumberOfObjects+j]) {
            for(size_t k = 0; k < RecvResultsPerObject[i][ParticleCounter]; k++) {
              double dist;

              auto ObjectItr = ObjectsBegin + j;
              auto ObjectResult = SearchResults[i].GetContainer()[ResultCounter];

              Results[j][NumberOfResults[j]] = ObjectResult;
              TObjectConfigure::Distance((*ObjectItr), ObjectResult,dist);
              ResultsDistances[j][NumberOfResults[j]] = dist;
              NumberOfResults[j]++;

              TObjectConfigure::UpdateGhostInterface(MpiCommunicator, ObjectResult, source);

              ResultCounter++;
            }
            ParticleCounter++;
          }
        }
      }
    }
  }
  /**
   * [FillPartition description]
   * @param Box      [description]
   * @param i_object [description]
   */
  // void FillPartition(SearchStructureType Box, const PointerType& i_object) {
  //   PointType  MinCell, MaxCell;
  //
  //   MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * this->mCellSize[0] + this->mMinPoint[0];  //
  //   MaxCell[0] = MinCell[0] + this->mCellSize[0];
  //
  //   for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=this->mCellSize[0], MaxCell[0]+=this->mCellSize[0] ) {
  //     if(TObjectConfigure::IntersectionBox(i_object, MinCell, MaxCell)) {
  //       auto & local  = (*(this->mCells[I].GetObject(0)))();
  //       auto & global = (*(this->mCells[I].GetObject(1)))();
  //
  //       local[mpi_rank] = 1;
  //       global[mpi_rank] = 1;
  //     }
  //   }
  // }

  /** Adds the current partition to the parition cell intersecting i_object
   * Adds the current partition to the parition cell intersecting i_object
   * @param Box      i_object boundingbox
   * @param i_object Object
   */
  void FillPartition(SearchStructureType & Box, const PointerType& i_object) {
    PointType  MinCell, MaxCell;
    PointType  MinBox, MaxBox;

    for(SizeType i = 0; i < 3; i++) {
      MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * this->mCellSize[i] + this->mMinPoint[i];
      MaxBox[i] = MinBox[i] + this->mCellSize[i];
    }

    MinCell[2] = MinBox[2];
    MaxCell[2] = MaxBox[2];
    for(IndexType III = Box.Axis[2].Begin(); III <= Box.Axis[2].End(); III += Box.Axis[2].Block, MinCell[2]+=this->mCellSize[2], MaxCell[2]+=this->mCellSize[2]) {
      MinCell[1] = MinBox[1];
      MaxCell[1] = MaxBox[1];
      for(IndexType II = III + Box.Axis[1].Begin(); II <= III + Box.Axis[1].End(); II += Box.Axis[1].Block, MinCell[1]+=this->mCellSize[1], MaxCell[1]+=this->mCellSize[1]) {
        MinCell[0] = MinBox[0];
        MaxCell[0] = MaxBox[0];
        for(IndexType I = II + Box.Axis[0].Begin(); I <= II + Box.Axis[0].End(); I += Box.Axis[0].Block, MinCell[0]+=this->mCellSize[0], MaxCell[0]+=this->mCellSize[0]) {
          if(TPartConfigure::IntersectionBox(i_object, MinCell, MaxCell)) {
            auto & cellPartitions = (*(this->mCells[I].GetObject(0)))(); // std::vector<std::size_t>
            cellPartitions.insert(mpi_rank);
          }
        }
      }
    }
  }

  /** Returns all partitions in the cells intersecting i_object
   * Returns all partitions in the cells intersecting i_object
   * @param Box      i_object boundingbox
   * @param i_object Object
   */
  void GetPartition(SearchStructureType & Box, const PointerType& i_object, std::unordered_set<std::size_t> & partitionList) {
    PointType  MinCell, MaxCell;
    PointType  MinBox, MaxBox;

    for(SizeType i = 0; i < 3; i++) {
      MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * this->mCellSize[i] + this->mMinPoint[i];  //
      MaxBox[i] = MinBox[i] + this->mCellSize[i];
    }

    MinCell[2] = MinBox[2];
    MaxCell[2] = MaxBox[2];
    for(IndexType III = Box.Axis[2].Begin(); III <= Box.Axis[2].End(); III += Box.Axis[2].Block, MinCell[2]+=this->mCellSize[2], MaxCell[2]+=this->mCellSize[2]) {
      MinCell[1] = MinBox[1];
      MaxCell[1] = MaxBox[1];
      for(IndexType II = III + Box.Axis[1].Begin(); II <= III + Box.Axis[1].End(); II += Box.Axis[1].Block, MinCell[1]+=this->mCellSize[1], MaxCell[1]+=this->mCellSize[1]) {
        MinCell[0] = MinBox[0];
        MaxCell[0] = MaxBox[0];
        for(IndexType I = II + Box.Axis[0].Begin(); I <= II + Box.Axis[0].End(); I += Box.Axis[0].Block, MinCell[0]+=this->mCellSize[0], MaxCell[0]+=this->mCellSize[0]) {
          if(TPartConfigure::IntersectionBox(i_object, MinCell, MaxCell)) {
            auto & cellPartitions = (*(this->mCells[I].GetObject(0)))();
            partitionList.insert(cellPartitions.begin(), cellPartitions.end());
          }
        }
      }
    }
  }

  /** Returns all partitions in the cells intersecting i_object with and extended radius 'Radius'
   * Returns all partitions in the cells intersecting i_object with and extended radius 'Radius'
   * @param Box      i_object boundingbox
   * @param i_object Object
   * @param Radius   radius extension
   */
  void GetPartitionInRadius(SearchStructureType & Box, const PointerType& i_object, std::unordered_set<std::size_t> & partitionList, double Radius) {
    PointType  MinCell, MaxCell;
    PointType  MinBox, MaxBox;

    for(SizeType i = 0; i < 3; i++) {
      MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * this->mCellSize[i] + this->mMinPoint[i];  //
      MaxBox[i] = MinBox[i] + this->mCellSize[i];
    }

    MinCell[2] = MinBox[2];
    MaxCell[2] = MaxBox[2];
    for(IndexType III = Box.Axis[2].Begin(); III <= Box.Axis[2].End(); III += Box.Axis[2].Block, MinCell[2]+=this->mCellSize[2], MaxCell[2]+=this->mCellSize[2]) {
      MinCell[1] = MinBox[1];
      MaxCell[1] = MaxBox[1];
      for(IndexType II = III + Box.Axis[1].Begin(); II <= III + Box.Axis[1].End(); II += Box.Axis[1].Block, MinCell[1]+=this->mCellSize[1], MaxCell[1]+=this->mCellSize[1]) {
        MinCell[0] = MinBox[0];
        MaxCell[0] = MaxBox[0];
        for(IndexType I = II + Box.Axis[0].Begin(); I <= II + Box.Axis[0].End(); I += Box.Axis[0].Block, MinCell[0]+=this->mCellSize[0], MaxCell[0]+=this->mCellSize[0]) {
          if(TPartConfigure::IntersectionBox(i_object, MinCell, MaxCell, Radius) || true) {
            auto & cellPartitions = (*(this->mCells[I].GetObject(0)))();
            partitionList.insert(cellPartitions.begin(), cellPartitions.end());
          }
        }
      }
    }
  }

  /**
   * Creates the container
   */
  void AllocateContainer()
  {
      SizeType Size = mN[0];
      for(SizeType i = 1 ; i < Dimension ; i++)
          Size *= mN[i];
      mCells.resize(Size);
  }

  /**
   * Mpi related variables
   * @mpi_rank: id of the current process for the given mpi communicator
   * @mpi_size: number of processes for the given mpi communicator
   */
  int mpi_rank;
  int mpi_size;

  /**
   * Bins used in the BinsMpi.
   * @mPartitionBin: Stores where the partitions are located in the space
   * @mObjectBin: Stores where the objects of the partition/s associated to this process are located in the space
   */

  PointType    mMinPoint;
  PointType    mMaxPoint;

  SizeType     mObjectsSize;
  IteratorType mObjectsBegin;
  IteratorType mObjectsEnd;

  CoordinateArray  mCellSize;
  CoordinateArray  mInvCellSize;
  SizeArray        mN;

  CellContainerType mCells;  ///The bin

  BinsObjectDynamic<TObjectConfigure> * mpObjectBins;

};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TObjectConfigure, class TPartConfigure>
inline std::istream& operator >> (
    std::istream& rIStream,
    BinsObjectDynamicMpi<TObjectConfigure>& rThis) {

  return rIStream;
}


/// output stream function
template<class TObjectConfigure, class TPartConfigure>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const BinsObjectDynamicMpi<TObjectConfigure, TPartConfigure> & rThis) {

  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_BINS_DYNAMIC_OBJECTS_MPI_CONTAINER_H_INCLUDED defined
