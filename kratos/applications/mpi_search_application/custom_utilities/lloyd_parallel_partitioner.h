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
#include "bins_dynamic_objects_mpi.h"
    
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
    
    //Lloyd based partitioning
    void LloydsBasedParitioner(ModelPart& mModelPart, double MaxNodeRadius, int CalculateBoundry)
    {
        //Elements
        ContainerType& pElements = mModelPart.GetCommunicator().LocalMesh().ElementsArray();
        
        //Stuff
        double SetCentroid[mpi_size*Dimension], SendSetCentroid[mpi_size*Dimension];
        
        //Boundary conditions
        double MeanPoint[mpi_size*Dimension];
        double Normal[mpi_size*Dimension];
        double Plane[mpi_size];
        double Dot[mpi_size];
        
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
        
        //Move Elements
        mModelPart.GetCommunicator().TransferModelElements(mModelPart);
        
        //Calculate partitions interfaces
        CalculatePartitionInterface(mModelPart,Normal,Plane,Dot);
    }
    
    /**
     * Calcuate the partition for all elements and nodes. This function must correctly set
     * PARTICLE_INDEX variable for each element and node.
     * @param pElements: List of all elements to be calcualted.
     * @param iterations: Number of iterations.
     **/
    void CalculatePartitions(ContainerType& pElements, const int &NumIterations, double SetCentroid[], double SendSetCentroid[])
    {
        int NumParticlesPerVirtualPartition[mpi_size], SendNumParticlesPerVirtualPartition[mpi_size];
      
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
                    i_nod->GetSolutionStepValue(PARTITION_INDEX) = PartitionIndex;
                }
            }
        }
    }
        
    void CalculatePartitionInterface(ModelPart& mModelPart, double Normal[], double Plane[], double Dot[])
    {
        ContainerType pElements = mModelPart.GetCommunicator().LocalMesh().ElementsArray();
        ContainerType pElementsMarked;
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        //Mark elements in boundary (Level-1)
        for (IteratorType particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it)
        {         
            (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(INTERNAL_ENERGY) = 9999;
            
//             (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(OSS_SWITCH) = 0; 
            (*particle_pointer_it)->GetValue(OSS_SWITCH) = 0;
            
            for(int j = 0; j < mpi_size; j++)
            {
                double NodeToCutPlaneDist = 0;
                int plane = j;

                NodeToCutPlaneDist += (*particle_pointer_it)->GetGeometry()(0)->X() * Normal[plane*Dimension+0];
                NodeToCutPlaneDist += (*particle_pointer_it)->GetGeometry()(0)->Y() * Normal[plane*Dimension+1];
                NodeToCutPlaneDist += (*particle_pointer_it)->GetGeometry()(0)->Z() * Normal[plane*Dimension+2];
                NodeToCutPlaneDist += Plane[plane];
                NodeToCutPlaneDist /= Dot[plane];

                double Radius = (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                int myRank = (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX);
                
                if (plane == myRank)
                    NodeToCutPlaneDist += 1;
                
                if(fabs(NodeToCutPlaneDist) <= Radius)
                {
                    //Fill the marked element vector
                    pElementsMarked.push_back((*particle_pointer_it));
                    
                    (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(OSS_SWITCH) |= ((1 << plane) | (1 << myRank));
                    (*particle_pointer_it)->GetValue(OSS_SWITCH) |= ((1 << plane) | (1 << myRank));        
                }
            }
               
            (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(INTERNAL_ENERGY) = (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(OSS_SWITCH);
            (*particle_pointer_it)->GetValue(INTERNAL_ENERGY) = (*particle_pointer_it)->GetValue(OSS_SWITCH);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        std::cout << "Level2 Marking" << std::endl;
        
        //Mark elements in boundary (Level-2) (Only last iteration)
        std::vector<std::vector<PointerType> > SendMarkObjects(mpi_size, std::vector<PointerType>(0));
        std::vector<std::vector<PointerType> > RecvMarkObjects(mpi_size, std::vector<PointerType>(0));
      
        //Mark neighbours (Level-2)
        std::cout << pElements.end()-pElements.begin() << std::endl;
        BinsObjectDynamicMpi<Configure> particle_bin(pElements.begin(), pElements.end());
        
        int NumberOfMarkedElements = pElementsMarked.end() - pElementsMarked.begin();
        int MaximumNumberOfResults = 1000;
        
        std::vector<std::size_t>               NumberOfResults(NumberOfMarkedElements);
        std::vector<std::vector<PointerType> > Results(NumberOfMarkedElements, std::vector<PointerType>(MaximumNumberOfResults));
        std::vector<std::vector<double> >      ResultsDistances(NumberOfMarkedElements, std::vector<double>(MaximumNumberOfResults));
        std::vector<double>                    Radius(NumberOfMarkedElements);
        
        double new_extension = 1.01;
        for (IteratorType particle_pointer_it = pElementsMarked.begin(); particle_pointer_it != pElementsMarked.end(); ++particle_pointer_it)
        {    
              Radius[particle_pointer_it - pElementsMarked.begin()] = new_extension * ((1.0 +0.3) * (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS)); //if this is changed, then compobation before adding neighbours must change also.
        }
        
        particle_bin.SearchObjectsMpi(pElementsMarked,NumberOfMarkedElements,Radius,Results,ResultsDistances,NumberOfResults,MaximumNumberOfResults,mModelPart.GetCommunicator());
        
        for (IteratorType particle_pointer_it = pElementsMarked.begin(); particle_pointer_it != pElementsMarked.end(); ++particle_pointer_it)
        {    
            int p_index = particle_pointer_it-pElementsMarked.begin();
            
            for(unsigned int i = 0; i < NumberOfResults[p_index]; i++)
            {
                Results[p_index][i]->GetGeometry()(0)->GetSolutionStepValue(OSS_SWITCH) = (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(OSS_SWITCH);
                SendMarkObjects[Results[p_index][i]->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX)].push_back(Results[p_index][i]);
            }
        }

//         mModelPart.GetCommunicator().TransferModelElements(mModelPart,SendMarkObjects,RecvMarkObjects);
        
        int msgSendSize[mpi_size];
        int msgRecvSize[mpi_size];

        TConfigure::AsyncSendAndReceive(mModelPart.GetCommunicator(),SendMarkObjects,RecvMarkObjects,msgSendSize,msgRecvSize);
        
        for (IteratorType particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it)
        {
            for(int j = 0; j < mpi_size; j++)
            {
                if(j != mpi_rank)
                {
                    for(unsigned int k = 0; k < RecvMarkObjects[j].size(); k++)
                    {
                        if((*particle_pointer_it)->GetGeometry()(0)->Id() == RecvMarkObjects[j][k]->GetGeometry()(0)->Id())
                        {
                            (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(OSS_SWITCH) |= RecvMarkObjects[j][k]->GetGeometry()(0)->GetSolutionStepValue(OSS_SWITCH);
                            (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(INTERNAL_ENERGY) = (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(OSS_SWITCH);
                        }
                    }
                }
            }
        }
    }

    /// Destructor.
    virtual ~LloydParallelPartitioner() {
        char msg[12] = {'b','i','n','s','_','X','.','t','i','m','e','\0'};
        msg[5] = '0' + mpi_rank;
        Timer::SetOuputFile(msg);
        Timer::PrintTimingInformation();
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


