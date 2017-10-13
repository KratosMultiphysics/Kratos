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

#if !defined(KRATOS_LLOYD_PARALLEL_PARTITIONER_H_INCLUDED)
#define  KRATOS_LLOYD_PARALLEL_PARTITIONER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

// Project includes
#include "mpi.h"
#include "spatial_containers/tree.h"
#include "spatial_containers/cell.h"

// Application includes
#include "custom_utilities/bins_dynamic_objects_mpi.h"
#include "processes/graph_coloring_process.h"

// Graph coloring
#include "processes/graph_coloring_process.h"

// TODO: This procedure seems unused. Maybe can be removed.
int compareFunction(const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}

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
class LloydParallelPartitioner
{
public:
    ///@name Type Definitions
    ///@{

    enum { Dimension = TConfigure::Dimension };

    // Point
    typedef TConfigure                                     Configure;
    typedef typename Configure::PointType                  PointType;
    // typedef typename TConfigure::ResultNumberIteratorType  ResultNumberIteratorType;

    // Container
    typedef typename Configure::PointerType               PointerType;
    typedef typename Configure::ContainerType             ContainerType;
    typedef typename ContainerType::iterator              IteratorType;
    typedef typename Configure::DistanceIteratorType      DistanceIteratorType;
    typedef typename Configure::ResultContainerType       ResultContainerType;
    typedef typename Configure::ElementsContainerType     ElementsContainerType;
    // typedef typename Configure::ResultPointerType         ResultPointerType;
    typedef typename Configure::ResultIteratorType        ResultIteratorType;
    typedef typename Configure::PointerContactType        PointerContactType;
    // typedef typename Configure::PointerTypeIterator       PointerTypeIterator;

    typedef WeakPointerVector<Element> ParticleWeakVector;

    // Search Structures
    typedef Cell<Configure> CellType;
    typedef std::vector<CellType> CellContainerType;
    typedef typename CellContainerType::iterator CellContainerIterator;

    typedef TreeNode<Dimension, PointType, PointerType, IteratorType, typename Configure::DistanceIteratorType> TreeNodeType;
    typedef typename TreeNodeType::CoordinateType         CoordinateType;  // double
    typedef typename TreeNodeType::SizeType               SizeType;        // std::size_t
    typedef typename TreeNodeType::IndexType              IndexType;       // std::size_t


    typedef Tvector<IndexType,Dimension>                  IndexArray;
    typedef Tvector<SizeType,Dimension>                   SizeArray;
    typedef Tvector<CoordinateType,Dimension>             CoordinateArray;

    ///Contact Pair
    typedef typename Configure::ContainerContactType      ContainerContactType;
    typedef typename Configure::IteratorContactType       IteratorContactType;

    ///typedef TreeNodeType LeafType;
    typedef typename TreeNodeType::IteratorIteratorType   IteratorIteratorType;
    typedef typename TreeNodeType::SearchStructureType    SearchStructureType;

    // Graph coloring process type
    typedef typename GraphColoringProcess::GraphType      GraphType;

    /// Pointer definition of BinsObjectDynamic
    KRATOS_CLASS_POINTER_DEFINITION(LloydParallelPartitioner);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LloydParallelPartitioner(IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd)
      : mNumberOfObjects(ObjectsEnd-ObjectsBegin), mObjectsBegin(ObjectsBegin), mObjectsEnd(ObjectsEnd)  {

      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

      // std::cout << "Begining partitioning" << std::endl;

      mpPartitionBins = new BinsObjectDynamicMpi<TConfigure>(mObjectsBegin, mObjectsEnd);
      mNumberOfCells = mpPartitionBins->GetCellContainer().size();

      if(mNumberOfCells < mpi_size) {
        KRATOS_ERROR << "Error: Number of cells in the bins must be at least equal to mpi_size. " << mNumberOfCells << std::endl;
      }

      if(mNumberOfCells % mpi_size) {
        // KRATOS_WARNING << "Warning: Number of cells is not multiple of mpi_size. Heavy imbalance may occur." << std::endl;
        std::cout << "Warning: Number of cells is not multiple of mpi_size. Heavy imbalance may occur. " << mNumberOfCells << std::endl;
      }

      if(mNumberOfCells < 10 * mpi_size) {
        // KRATOS_WARNING << "Warning: Number of cells is small. Partition Shape may be sub-optimal." << std::endl;
        std::cout << "Warning: Number of cells is small. Partition Shape may be sub-optimal. " << mNumberOfCells << std::endl;
      }
    }

    double ReduceMaxRadius(IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd) {
      // Max Radius Ugly fix
      double local_max_radius = 0.0f;
      double max_radius = 0.0f;
      for (IteratorType ObjectItr = ObjectsBegin; ObjectItr != ObjectsEnd; ObjectItr++) {
        const double Radius = TConfigure::GetObjectRadius(*ObjectItr, 0.0f);
        if(Radius > local_max_radius) local_max_radius = Radius;
      }

      MPI_Allreduce(&local_max_radius, &max_radius, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      return max_radius;
    }

    void SerialPartition() {
      std::vector<int> mpiSendObjectsPerCell(mNumberOfCells, 0);
      std::vector<int> mpiRecvObjectsPerCell(mNumberOfCells, 0);

      std::vector<int> CellPartition(mNumberOfCells, 0);
      std::vector<int> ObjectsPerPartition(mpi_size, 0);

      int mpiSendNumberOfObjects = mObjectsEnd - mObjectsBegin;
      int mpiRecvNumberOfObjects = 0;

      PointType ObjectCenter;
      PointType Low, High;

      SearchStructureType Box;

      // Calculate objects per cell
      for(std::size_t i = 0; i < (std::size_t)mNumberOfObjects; i++) {
        auto ObjectItr = mObjectsBegin + i;
        TConfigure::CalculateCenter(*ObjectItr, ObjectCenter);
        auto cellId = mpPartitionBins->CalculateIndex(ObjectCenter);
        mpiSendObjectsPerCell[cellId]++;
      }

      MPI_Allreduce(&mpiSendNumberOfObjects, &mpiRecvNumberOfObjects, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&mpiSendObjectsPerCell[0], &mpiRecvObjectsPerCell[0], mNumberOfCells, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      //int MeanObjectsPerPartition = mpiRecvNumberOfObjects / mpi_size;
      // std::cout << "mpiRecvNumberOfObjects: " << mpiRecvNumberOfObjects << " MeanObjectsPerPartition: " << MeanObjectsPerPartition << std::endl;

      // Assing each cell to the closest partition center
      // TODO: this is currently very unbalanced
      for(std::size_t cellId = 0; cellId < (std::size_t) mNumberOfCells; cellId++) {
        ObjectsPerPartition[cellId] += mpiRecvObjectsPerCell[cellId];
        CellPartition[cellId] = cellId;
      }

      std::cout << "Partititon " << mpi_rank << ": " << ObjectsPerPartition[mpi_rank] << std::endl;

      // Assign the partition to the objects based on their cell
      for(std::size_t i = 0; i < (std::size_t) mNumberOfObjects; i++) {
        auto ObjectItr = mObjectsBegin + i;

        TConfigure::CalculateCenter(*ObjectItr, ObjectCenter);

        auto cellId = mpPartitionBins->CalculateIndex(ObjectCenter);

        (*ObjectItr)->GetValue(PARTITION_INDEX) = CellPartition[cellId];
        for (unsigned int j = 0; j < (*ObjectItr)->GetGeometry().PointsNumber(); j++) {
          ModelPart::NodeType::Pointer NodePtr = (*ObjectItr)->GetGeometry().pGetPoint(j);
          NodePtr->FastGetSolutionStepValue(PARTITION_INDEX) = CellPartition[cellId];
        }
      }

      // std::cout << "Ending partitioning" << std::endl;
    }

    void VoronoiiPartition() {

      std::vector<int> mpiSendObjectsPerCell(mNumberOfCells, 0);
      std::vector<int> mpiRecvObjectsPerCell(mNumberOfCells, 0);

      std::vector<int> CellPartition(mNumberOfCells, 0);
      std::vector<double> CellDistances(mNumberOfCells, std::numeric_limits<double>::max());

      std::vector<double> mpiSendCellCenter(mNumberOfCells * Dimension, 0.0f);
      std::vector<double> mpiRecvCellCenter(mNumberOfCells * Dimension, 0.0f);

      std::vector<int> CellsPerPartition(mpi_size, 0);
      std::vector<int> ObjectsPerPartition(mpi_size, 0);
      std::vector<PointType> PartitionCenters(mpi_size);

      std::vector<double> mpiSendPartCenter(mpi_size * Dimension, 0.0f);
      std::vector<double> mpiRecvPartCenter(mpi_size * Dimension, 0.0f);

      std::vector<int> mpiSendPartNum(mpi_size, 0);
      std::vector<int> mpiRecvPartNum(mpi_size, 0);

      PointType ObjectCenter;

      // 1 - Calculate the centers of the cells based on the objects inside
      // TODO: Parallelize this (non-trivial)
      for(std::size_t i = 0; i < mNumberOfObjects; i++) {
        auto ObjectItr = mObjectsBegin + i;

        TConfigure::CalculateCenter(*ObjectItr, ObjectCenter);

        auto CellIndex = mpPartitionBins->CalculateIndex(ObjectCenter);

        mpiSendObjectsPerCell[CellIndex]++;
        for(int d = 0; d < Dimension; d++) {
          mpiSendCellCenter[CellIndex*Dimension+d] += ObjectCenter[d];
        }
      }

      // 1.1 - Communicate the number of objects per cell and the local sum of object coordinates
      MPI_Allreduce(&mpiSendObjectsPerCell[0], &mpiRecvObjectsPerCell[0], mNumberOfCells * Dimension, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&mpiSendCellCenter[0], &mpiRecvCellCenter[0], mNumberOfCells * Dimension, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      // 1.2 - Obtain the wheighted center of each cell with the data of all processes
      #pragma omp parallel for
      for(std::size_t cellId = 0; cellId < mNumberOfCells; cellId++) {
        for(int d = 0; d < Dimension; d++) {
          mpiRecvCellCenter[cellId*Dimension+d] /= mpiRecvObjectsPerCell[cellId];
        }
      }

      // 2 - Assign a random origin to each partition.
      auto minPoint = mpPartitionBins->GetMinPoint();
      auto maxPoint = mpPartitionBins->GetMaxPoint();
      auto boxSize  = maxPoint - minPoint;

      // Change this if we want real random
      // !!!!!MAKE SURE THIS IS THE SAME ON EVERY PARTITION OR IT WON'T WORK!!!!!
      std::srand(256);

      for(int i = 0; i < mpi_size; i++) {
        for(int d = 0; d < Dimension; d++) {
          PartitionCenters[i][d] = minPoint[d] + ((double)std::rand() / (double)RAND_MAX) * boxSize[d];
        }
      }

      // While not converged
      auto MaxIterations = 1e1;
      for(std::size_t iterations = 0; iterations < MaxIterations; iterations++ ) {

        // Assing each cell to the closest partition center
        for(std::size_t cellId = 0; cellId < mNumberOfCells; cellId++) {
          if(mpiRecvObjectsPerCell[cellId] != 0) {
            for(int i = 0; i < mpi_size; i++) {
              double cubeDistance = 0.0f;

              for(int d = 0; d < Dimension; d++) {
                // Manhattan distance shoudl prevent problems with the discretization of the space
                cubeDistance += std::abs(mpiRecvCellCenter[cellId*Dimension+d] - PartitionCenters[i][d]);
              }

              if(cubeDistance < CellDistances[cellId]) {
                CellDistances[cellId] = cubeDistance;
                CellPartition[cellId] = i;
              }
            }
          }
        }

        // At this point no synch should be needed

        // Update the center of the partitions
        for(int i = 0; i < mpi_size; i++) {
          CellsPerPartition[i] = 0;
          ObjectsPerPartition[i] = 0;
          for(int d = 0; d < Dimension; d++) {
            PartitionCenters[i][d] = 0.0f;
          }
        }

        // if(mpi_rank == 0) {
        //   std::cout << mNumberOfCells << std::endl;
        // }
        for(std::size_t cellId = 0; cellId < mNumberOfCells; cellId++) {
          if(mpiRecvObjectsPerCell[cellId] != 0) {
            CellsPerPartition[CellPartition[cellId]]++;
            ObjectsPerPartition[CellPartition[cellId]] += mpiRecvObjectsPerCell[cellId];
            for(int d = 0; d < Dimension; d++) {
              PartitionCenters[CellPartition[cellId]][d] += mpiRecvCellCenter[cellId*Dimension+d];
            }
          }
          // if(mpi_rank == 0) {
          //   for(int i = 0; i < mpi_size; i++) {
          //     std::cout << "Iteration: " << cellId << " Partition " << i << " has " << CellsPerPartition[i] << " Cells" << std::endl;
          //   }
          // }
        }

        for(std::size_t partId = 0; partId < mpi_size; partId++) {
          for(int d = 0; d < Dimension; d++) {
            PartitionCenters[partId][d] /= CellsPerPartition[partId]++;
          }
        }
      }

      if(mpi_rank == 0) {
        std::cout << mNumberOfCells << std::endl;
        for(int i = 0; i < mpi_size; i++) {
          std::cout << "Partition " << i << " has " << CellsPerPartition[i] << " Cells" << std::endl;
        }
      }

      // Assign the partition to the objects based on their cell
      for(std::size_t i = 0; i < mNumberOfObjects; i++) {
        auto ObjectItr = mObjectsBegin + i;

        TConfigure::CalculateCenter(*ObjectItr, ObjectCenter);

        auto CellIndex = mpPartitionBins->CalculateIndex(ObjectCenter);

        (*ObjectItr)->GetValue(PARTITION_INDEX) = CellPartition[CellIndex];
        for (unsigned int i = 0; i < (*ObjectItr)->GetGeometry().PointsNumber(); i++) {
          ModelPart::NodeType::Pointer NodePtr = (*ObjectItr)->GetGeometry().pGetPoint(i);
          NodePtr->FastGetSolutionStepValue(PARTITION_INDEX) = CellPartition[CellIndex];
        }
      }

      std::cout << "Ending partitioning" << std::endl;
    }

    void UpdateDomainGraph(IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd, GraphType & domainGraph) {
      PointType ObjectCenter;
      PointType Low, High;

      SearchStructureType Box;

      mObjectsBegin = ObjectsBegin;
      mObjectsEnd = ObjectsEnd;
      mNumberOfObjects = ObjectsEnd-ObjectsBegin;

      // Rebuild the bins
      free(mpPartitionBins);
      mpPartitionBins = new BinsObjectDynamicMpi<TConfigure>(mObjectsBegin, mObjectsEnd);

      // Assign the partition to the objects based on their cell
      double maxRadius = ReduceMaxRadius(mObjectsBegin, mObjectsEnd);

      for(std::size_t i = 0; i < (std::size_t) mNumberOfObjects; i++) {
        auto ObjectItr = mObjectsBegin + i;

        TConfigure::CalculateBoundingBox(*ObjectItr, Low, High);
        for(int i = 0; i < Dimension; i++) {
          Low[i] -= maxRadius;
          High[i] += maxRadius;
        }

        Box.Set( mpPartitionBins->CalculateCell(Low), mpPartitionBins->CalculateCell(High), mpPartitionBins->GetDivisions());

        std::unordered_set<std::size_t> partitionSet;
        auto ObjectRadius = TConfigure::GetObjectRadius(*ObjectItr, 0.0f);
        mpPartitionBins->SearchPartitionInRadius(Box, *ObjectItr, partitionSet, ObjectRadius);

        std::vector<std::size_t> partitionList(partitionSet.begin(), partitionSet.end());

        for(unsigned int i = 0; i < partitionList.size(); i++) {
          domainGraph(mpi_rank, mpi_rank) = 1;
          domainGraph(partitionList[i], mpi_rank) = 1;
          domainGraph(mpi_rank, partitionList[i]) = 1;
        }
      }
    }

    /// Destructor.
    virtual ~LloydParallelPartitioner() {
      delete mpPartitionBins;
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

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name MPI Variables
    ///@{

    int mpi_rank;
    int mpi_size;

    int mNumberOfObjects;
    int mNumberOfCells;

    IteratorType mObjectsBegin;
    IteratorType mObjectsEnd;

    BinsObjectDynamicMpi<TConfigure> * mpPartitionBins;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    inline void CreatePartition(SizeType number_of_threads, const SizeType number_of_rows, std::vector<SizeType>& partitions) {
      partitions.resize(number_of_threads+1);
      SizeType partition_size = number_of_rows / number_of_threads;
      partitions[0] = 0;
      partitions[number_of_threads] = number_of_rows;

      for(SizeType i = 1; i<number_of_threads; i++) {
        partitions[i] = partitions[i-1] + partition_size;
      }
    }


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
    LloydParallelPartitioner<TConfigure> & operator=(const LloydParallelPartitioner<TConfigure> & rOther) {
      mObjectsBegin = rOther.mObjectsBegin;
      mObjectsEnd = rOther.mObjectsEnd;

      mpPartitionBins = rOther.mpPartitionBins;

      return *this;
    }

    /// Copy constructor.
    LloydParallelPartitioner(const LloydParallelPartitioner& rOther) {
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
                                  LloydParallelPartitioner<TConfigure>& rThis)
{
    return rIStream;
}


/// output stream function
template<class TConfigure>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const LloydParallelPartitioner<TConfigure> & rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_LLOYD_PARALLEL_PARTITIONER_H_INCLUDED  defined
