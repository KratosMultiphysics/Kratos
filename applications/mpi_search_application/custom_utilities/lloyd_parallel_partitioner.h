/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Carlos Roig $
//   Date:                $Date: 13-08-2012 $
//   Revision:            $Revision: 1.1.1.1 $
//
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

int compareFunction(const void * a, const void * b)
{
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

    //Point
    typedef TConfigure                                      Configure;
    typedef typename Configure::PointType                  PointType;
    //typedef typename TConfigure::ResultNumberIteratorType   ResultNumberIteratorType;

    //Container
    typedef typename Configure::PointerType                PointerType;
    typedef typename Configure::ContainerType              ContainerType;
    typedef typename ContainerType::iterator               IteratorType;
    typedef typename Configure::DistanceIteratorType       DistanceIteratorType;
    typedef typename Configure::ResultContainerType        ResultContainerType;
    typedef typename Configure::ElementsContainerType       ElementsContainerType;
//     typedef typename Configure::ResultPointerType          ResultPointerType;
    typedef typename Configure::ResultIteratorType         ResultIteratorType;
    typedef typename Configure::PointerContactType         PointerContactType;
//     typedef typename Configure::PointerTypeIterator        PointerTypeIterator;

    typedef WeakPointerVector<Element> ParticleWeakVector;

    //Search Structures
    typedef Cell<Configure> CellType;
    typedef std::vector<CellType> CellContainerType;
    typedef typename CellContainerType::iterator CellContainerIterator;

    typedef TreeNode<Dimension, PointType, PointerType, IteratorType,  typename Configure::DistanceIteratorType> TreeNodeType;
    typedef typename TreeNodeType::CoordinateType  CoordinateType;  // double
    typedef typename TreeNodeType::SizeType        SizeType;        // std::size_t
    typedef typename TreeNodeType::IndexType       IndexType;       // std::size_t


    typedef Tvector<IndexType,Dimension>      IndexArray;
    typedef Tvector<SizeType,Dimension>       SizeArray;
    typedef Tvector<CoordinateType,Dimension> CoordinateArray;

    ///Contact Pair
    typedef typename Configure::ContainerContactType  ContainerContactType;
    typedef typename Configure::IteratorContactType IteratorContactType;

    ///typedef TreeNodeType LeafType;
    typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
    typedef typename TreeNodeType::SearchStructureType  SearchStructureType;

    PointerType CommunicationToken;

    /// Pointer definition of BinsObjectDynamic
    KRATOS_CLASS_POINTER_DEFINITION(LloydParallelPartitioner);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LloydParallelPartitioner()
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    }

    /// Destructor.
    virtual ~LloydParallelPartitioner() {
    }

    Vector MeanPoint;
    Vector Normal;
    Vector Plane;
    Vector Dot;

    //Lloyd based partitioning
    void LloydsBasedPartitioner(ModelPart& mModelPart, double MaxNodeRadius, int CalculateBoundary)
    {
        //Elements
        ContainerType& pElements = mModelPart.GetCommunicator().LocalMesh().ElementsArray();

        //Stuff
        double * SetCentroid = new double[mpi_size*Dimension];
        double * SendSetCentroid = new double[mpi_size*Dimension];

        //Boundary conditions
        MeanPoint.resize(mpi_size*Dimension,false);
        Normal.resize(mpi_size*Dimension,false);
        Plane.resize(mpi_size,false);
        Dot.resize(mpi_size,false);

        //Define algorthm iterations (maybe is a good idea pass this as a parameter)
        int NumIterations = 100;

        //Calcualate Partitions
        CalculatePartitions(pElements,NumIterations,SetCentroid,SendSetCentroid);

        //Calculate boundaries
        for(int i = 0; i < mpi_size; i++)
        {
            Dot[i] = 0;
            Plane[i] = 0;

            for(int j = 0; j < Dimension; j++)
            {
                MeanPoint[i*Dimension+j] = (SetCentroid[i*Dimension+j] + SetCentroid[mpi_rank*Dimension+j])/2;
                Normal[i*Dimension+j]    = (SetCentroid[i*Dimension+j] - SetCentroid[mpi_rank*Dimension+j]);
                Dot[i]                  += Normal[i*Dimension+j] * Normal[i*Dimension+j];
                Plane[i]                -= Normal[i*Dimension+j] * MeanPoint[i*Dimension+j];
            }

            Dot[i] = sqrt(Dot[i]);
        }

        delete [] SetCentroid;
        delete [] SendSetCentroid;
    }

    /**
     * Calcuate the partition for all elements and nodes. This function must correctly set
     * PARTICLE_INDEX variable for each element and node.
     * @param pElements: List of all elements to be calcualted.
     * @param iterations: Number of iterations.
     **/
    void CalculatePartitions(ContainerType& pElements, const int &NumIterations, double SetCentroid[], double SendSetCentroid[])
    {
        int * NumParticlesPerVirtualPartition = new int[mpi_size];
        int * SendNumParticlesPerVirtualPartition = new int[mpi_size];

        for(int iterations = 0; iterations < NumIterations; iterations++)
        {
            for(int i = 0; i < mpi_size; i++)
            {
                SendNumParticlesPerVirtualPartition[i] = 0;

                for(int j = 0; j < Dimension; j++)
                    SendSetCentroid[i*Dimension+j] = 0;
            }

            for (IteratorType particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it)
            {
                for (unsigned int i = 0; i < (*particle_pointer_it)->GetGeometry().PointsNumber(); i++)
                {
                    ModelPart::NodeType::Pointer i_nod = (*particle_pointer_it)->GetGeometry().pGetPoint(i);
                    int element_partition = i_nod->GetSolutionStepValue(PARTITION_INDEX);

                    for(int j = 0; j < Dimension; j++)
                    {
                        SendSetCentroid[element_partition*Dimension+j] += i_nod->Coordinate(j+1);
                    }

                    SendNumParticlesPerVirtualPartition[element_partition]++;
                }
            }

            MPI_Allreduce(SendSetCentroid,SetCentroid,mpi_size*Dimension,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            MPI_Allreduce(SendNumParticlesPerVirtualPartition,NumParticlesPerVirtualPartition,mpi_size,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

            //Set the centroids
            if(NumParticlesPerVirtualPartition[mpi_rank] == 0)
            {
                for(int j = 0; j < Dimension; j++)
                {
                    SendSetCentroid[mpi_rank*Dimension+j] = 0;
                    SetCentroid[mpi_rank*Dimension+j] = 0;
                }
            }
            else
            {
                for(int i = 0; i < mpi_size; i++)
                {
                    for(int j = 0; j < Dimension; j++)
                    {
                        SendSetCentroid[i*Dimension+j] /= NumParticlesPerVirtualPartition[i];
                        SetCentroid[i*Dimension+j] /= NumParticlesPerVirtualPartition[i];
                    }
                }
            }

            //Continue with the partitioning
            for (IteratorType particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it)
            {
                for (unsigned int i = 0; i < (*particle_pointer_it)->GetGeometry().PointsNumber(); i++)
                {
                    ModelPart::NodeType::Pointer i_nod = (*particle_pointer_it)->GetGeometry().pGetPoint(i);

                    int PartitionIndex = mpi_rank;

                    //Calcualte distance with my center first.
                    double distX = i_nod->X() - SetCentroid[mpi_rank*Dimension+0];
                    double distY = i_nod->Y() - SetCentroid[mpi_rank*Dimension+1];
                    double distZ = i_nod->Z() - SetCentroid[mpi_rank*Dimension+2];

                    distX *= distX;
                    distY *= distY;
                    distZ *= distZ;

                    double dist = (distX+distY+distZ);
                    double min_dist = dist;

                    for(int j = 0; j < mpi_size; j++)
                    {
                        distX = i_nod->X() - SetCentroid[j*Dimension+0];
                        distY = i_nod->Y() - SetCentroid[j*Dimension+1];
                        distZ = i_nod->Z() - SetCentroid[j*Dimension+2];

                        distX *= distX;
                        distY *= distY;
                        distZ *= distZ;

                        double TheyDist = (distX+distY+distZ);

                        if (TheyDist < min_dist)
                        {
                            min_dist = TheyDist;
                            PartitionIndex = j;
                        }
                    }

                    (*particle_pointer_it)->GetValue(PARTITION_INDEX) = PartitionIndex;
                    i_nod->FastGetSolutionStepValue(PARTITION_INDEX) = PartitionIndex;
                }
            }
        }

        delete [] NumParticlesPerVirtualPartition;
        delete [] SendNumParticlesPerVirtualPartition;

    }

    void CalculatePartitionInterface(ModelPart& mModelPart)
    {
        ContainerType pElements = mModelPart.GetCommunicator().LocalMesh().ElementsArray();
        ContainerType pElementsMarked;

        MPI_Barrier(MPI_COMM_WORLD);

        //Mark elements in boundary (Level-1)
        for (IteratorType particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it)
        {
            double Radius       = (*particle_pointer_it)->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
            int myRank          = (*particle_pointer_it)->GetGeometry()[0].FastGetSolutionStepValue(PARTITION_INDEX);

            (*particle_pointer_it)->GetGeometry()[0].FastGetSolutionStepValue(PARTITION_MASK) = 0;

            for(int j = 0; j < mpi_size; j++)
            {
                double NodeToCutPlaneDist = 0;

                NodeToCutPlaneDist += (*particle_pointer_it)->GetGeometry()[0].X() * Normal[j*Dimension+0];
                NodeToCutPlaneDist += (*particle_pointer_it)->GetGeometry()[0].Y() * Normal[j*Dimension+1];
                NodeToCutPlaneDist += (*particle_pointer_it)->GetGeometry()[0].Z() * Normal[j*Dimension+2];

                NodeToCutPlaneDist += Plane[j];
                NodeToCutPlaneDist /= Dot[j];

                NodeToCutPlaneDist  = fabs(NodeToCutPlaneDist);

                if (j != myRank)
                {
                    if( NodeToCutPlaneDist <= Radius*2)
                    {
                       (*particle_pointer_it)->GetGeometry()[0].FastGetSolutionStepValue(PARTITION_MASK) |= ((1 << j) | (1 << myRank));
                    }
                }
            }
        }
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

    PointType    mMinPoint;
    PointType    mMaxPoint;

    CoordinateArray  mCellSize;
    CoordinateArray  mInvCellSize;
    SizeArray        mN;

    CellContainerType mCells;  ///The bin

    ///@}
    ///@name MPI Variables
    ///@{

    int mpi_rank;
    int mpi_size;
    int mpi_total_elements;

    int ParticlesInDomain;
    int IdealParticlesInDomain;

    vector<int> mpi_connectivity;
    vector<vector<double> > mpi_MinPoints;
    vector<vector<double> > mpi_MaxPoints;

    vector<PointType> mMinBoundingBox;
    vector<PointType> mMaxBoundingBox;

    static double MyDimish;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

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
    LloydParallelPartitioner<TConfigure> & operator=(const LloydParallelPartitioner<TConfigure> & rOther)
    {
        mMinPoint            = rOther.mMinPoint;
        mMaxPoint            = rOther.mMaxPoint;
        mCellSize            = rOther.mCellSize;
        mInvCellSize         = rOther.mInvCellSize;
        mN                   = rOther.mN;
        mCells               = rOther.mCells;
        return *this;
    }

    /// Copy constructor.
    LloydParallelPartitioner(const LloydParallelPartitioner& rOther)
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
