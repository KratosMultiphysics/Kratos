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


#if !defined(KRATOS_MORTON_PARALLEL_PARTITIONER_H_INCLUDED)
#define  KRATOS_MORTON_PARALLEL_PARTITIONER_H_INCLUDED



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

// #ifdef _OPENMP
// #include <omp.h>
// #endif

    
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
class MortonParallelPartitioner
{
public:
    ///@name Type Definitions
    ///@{

    enum { Dimension = TConfigure::Dimension };

    //Point
    typedef TConfigure                                      Configure;
    typedef typename TConfigure::PointType                  PointType;
    //typedef typename TConfigure::ResultNumberIteratorType   ResultNumberIteratorType;
    
    //Container
    typedef typename TConfigure::PointerType                PointerType;
    typedef typename TConfigure::ContainerType              ContainerType;
    typedef typename TConfigure::IteratorType               IteratorType;
    typedef typename TConfigure::DistanceIteratorType       DistanceIteratorType;
    typedef typename TConfigure::ResultContainerType        ResultContainerType;
    typedef typename TConfigure::ResultPointerType          ResultPointerType;
    typedef typename TConfigure::ResultIteratorType         ResultIteratorType;
    typedef typename TConfigure::PointerContactType         PointerContactType;
    typedef typename TConfigure::PointerTypeIterator        PointerTypeIterator;
    
    typedef WeakPointerVector<Element> ParticleWeakVector;

    //Search Structures
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
    
    PointerType CommunicationToken;

    /// Pointer definition of BinsObjectDynamic
    KRATOS_CLASS_POINTER_DEFINITION(MortonParallelPartitioner);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MortonParallelPartitioner() 
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    }
    
    void CommunicationPhase(ModelPart& mModelPart,
                            std::vector<std::vector<PointerType> > &SendObjects,
                            std::vector<std::vector<PointerType> > &RecvObjects) 
    {
        int msgSendSize[mpi_size];
        int msgRecvSize[mpi_size];
        
        //Let this loop ¡HERE! because of reasons.
        for(int i = 0; i < mpi_size; i++)
        {
            msgSendSize[i] = 0;
            msgRecvSize[i] = 0;
        }
      
        mModelPart.GetCommunicator().AsyncSendAndReceiveNodes(SendObjects,RecvObjects,msgSendSize,msgRecvSize);
//         TConfigure::AsyncSendAndRecive(SendObjects,RecvObjects,msgSendSize,msgRecvSize);
    }
    
    void Repartitionate(ModelPart& mModelPart, std::vector<std::vector<PointerType> > &RecvObjects)
    {
        ContainerType& ElementsLocal  = mModelPart.GetCommunicator().LocalMesh().ElementsArray();
        ContainerType& ElementsGlobal = mModelPart.ElementsArray();
        
        ModelPart::NodesContainerType& NodesLocal  = mModelPart.GetCommunicator().LocalMesh().Nodes();
        ModelPart::NodesContainerType& NodesGlobal = mModelPart.Nodes();
   
        ContainerType temp_particles_container_local;
        ContainerType temp_particles_container_global;
        
        ModelPart::NodesContainerType temp_nodes_container_local;
        ModelPart::NodesContainerType temp_nodes_container_global;
        
        temp_particles_container_local.reserve(ElementsLocal.size());
        temp_particles_container_global.reserve(ElementsGlobal.size());
        
        temp_nodes_container_local.reserve(NodesLocal.size());
        temp_nodes_container_global.reserve(NodesGlobal.size());

        temp_particles_container_local.swap(ElementsLocal);
        temp_particles_container_global.swap(ElementsGlobal);
        
        temp_nodes_container_local.swap(NodesLocal);
        temp_nodes_container_global.swap(NodesGlobal);

        //Keep Local elements and nodes from our domain
        for (IteratorType particle_pointer_it = temp_particles_container_local.begin();
             particle_pointer_it != temp_particles_container_local.end(); ++particle_pointer_it)
        {
            if((*particle_pointer_it)->GetValue(PARTITION_INDEX) == mpi_rank)
            {                 
                mModelPart.GetCommunicator().LocalMesh().Elements().push_back(*particle_pointer_it);
                for (unsigned int i = 0; i < (*particle_pointer_it)->GetGeometry().PointsNumber(); i++)
                {
                    ModelPart::NodeType::Pointer pNode = (*particle_pointer_it)->GetGeometry().pGetPoint(i);
                    mModelPart.GetCommunicator().LocalMesh().Nodes().push_back(pNode);
                }
            }
        }

        //Keep Global elements and nodes from our domain
        for (IteratorType particle_pointer_it = temp_particles_container_global.begin();
             particle_pointer_it != temp_particles_container_global.end(); ++particle_pointer_it)
        {
            if((*particle_pointer_it)->GetValue(PARTITION_INDEX) == mpi_rank)
            {                 
                mModelPart.Elements().push_back(*particle_pointer_it);
                for (unsigned int i = 0; i < (*particle_pointer_it)->GetGeometry().PointsNumber(); i++)
                {
                    ModelPart::NodeType::Pointer pNode = (*particle_pointer_it)->GetGeometry().pGetPoint(i);
                    mModelPart.Nodes().push_back(pNode);
                }
            }
        }

        //Add new elements and nodes
        for(int i = 0; i < mpi_size; i++)
        {
            std::vector<PointerType>& RecvProcessObjects = RecvObjects[i];
            
            for(size_t j = 0; j < RecvProcessObjects.size(); j++)
            {
                if(RecvProcessObjects[j])
                {  
                    mModelPart.GetCommunicator().LocalMesh().Elements().push_back(RecvProcessObjects[j]);
                    mModelPart.Elements().push_back(RecvProcessObjects[j]);
                    
                    for (unsigned int i = 0; i < RecvProcessObjects[j]->GetGeometry().PointsNumber(); i++)
                    {
                        ModelPart::NodeType::Pointer pNode = RecvProcessObjects[j]->GetGeometry().pGetPoint(i);
                        
                        pNode->GetSolutionStepValue(PARTITION_INDEX) = mpi_rank;
                        
                        mModelPart.Nodes().push_back(pNode);
                        mModelPart.GetCommunicator().LocalMesh().Nodes().push_back(pNode);
                    }
                }
            }
        }
        
        mModelPart.GetCommunicator().LocalMesh().Nodes().Unique();
        
        for (unsigned int i = 0; i < mModelPart.GetCommunicator().LocalMeshes().size(); i++)
            mModelPart.GetCommunicator().LocalMesh(i).Nodes().Unique();
            
        for (unsigned int i = 0; i < mModelPart.GetCommunicator().GhostMeshes().size(); i++)
            mModelPart.GetCommunicator().GhostMesh(i).Nodes().Unique();
    }
    
    //Paritcionament basat en l'algortime LLoyds per fer algo similar a una tesselacio de voronoi
    //el que vui consegir es el clustering de elements que no es fa amb morton pero amb la propietat que sigin
    //dominis no fixes.
    void LloydsBasedParitioner(ModelPart& mModelPart, double MaxNodeRadius, int CalculateBoundry)
    {
        std::vector<std::vector<PointerType> > SendObjects(mpi_size, std::vector<PointerType>(mModelPart.NumberOfElements()));
        std::vector<std::vector<PointerType> > RecvObjects(mpi_size, std::vector<PointerType>(mModelPart.NumberOfElements()));
    
        DefineKSets(mModelPart,SendObjects,RecvObjects,MaxNodeRadius,CalculateBoundry); //A l'algortime clasic definiriam N sets al azar. Aprofitem i ho fem amb els punts de cada domini
//         CommunicationPhase(mModelPart,SendObjects,RecvObjects); //Comuniquem els resultats
//         Repartitionate(mModelPart,RecvObjects); // I per ultim reparticionem
    }
    
    //Todo: Clean this mess
    void DefineKSets(ModelPart& mModelPart, 
                     std::vector<std::vector<PointerType> > &SendObjects, 
                     std::vector<std::vector<PointerType> > &RecvObjects,
                     double MaxNodeRadius,
                     int CalculateBoundry)
    {
        std::cout << "Entra al reparticionat" << std::endl;

        ContainerType& pElements = mModelPart.GetCommunicator().LocalMesh().ElementsArray();
      
        double SetCentroid[mpi_size*Dimension], SendSetCentroid[mpi_size*Dimension];
        int PartitionNumberOfElements[mpi_size], SendPartitionNumberOfElements[mpi_size];
        
        double MeanPoint[mpi_size*Dimension];
        double Normal[mpi_size*Dimension];
        double Plane[mpi_size];
        double Dot[mpi_size];
        
        double MidDistance[mpi_size];
        double FulDistance[mpi_size];
        
//         double thisDiferentialDistance = 0;
        
        //Mog: fix me
        static std::vector<double> LastDistance(0);
        
        if(LastDistance.size() == 0)
            for(int i = 0 ; i < mpi_size; i++)
                LastDistance.push_back(0);
        
//         double SelfWeight[mpi_size];
        double TheyWeight[mpi_size];
        
        int NewNumberOfObjects[mpi_size];
        
        int LocalParticles = mModelPart.GetCommunicator().LocalMesh().NumberOfElements();
        int TotalParticles = 0;
        
        MPI_Allreduce(&LocalParticles,&TotalParticles,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        
        SendPartitionNumberOfElements[mpi_rank] = LocalParticles;
        
        for(int i = 0; i < mpi_size; i++)
            NewNumberOfObjects[i] = 0;

        for(int i = 0; i < Dimension; i++)
            SendSetCentroid[mpi_rank*Dimension+i] = 0;
        
        for (IteratorType particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it)
        {
            for (unsigned int i = 0; i < (*particle_pointer_it)->GetGeometry().PointsNumber(); i++)
            {
                ModelPart::NodeType::Pointer i_nod = (*particle_pointer_it)->GetGeometry().pGetPoint(i);
                
                for(int j = 0; j < Dimension; j++)
                {
                    SendSetCentroid[mpi_rank*Dimension+j] += i_nod->Coordinate(j+1);
                    SetCentroid[mpi_rank*Dimension+j] = SendSetCentroid[mpi_rank*Dimension+j];
                }
            }
        }

        //Set the centroids
        if(LocalParticles == 0)
        {
            for(int j = 0; j < Dimension; j++)
                {
                    SendSetCentroid[mpi_rank*Dimension+j] = 0;
                    SetCentroid[mpi_rank*Dimension+j] = SendSetCentroid[mpi_rank*Dimension+j];
                }
        } 
        else
        {                   
            for(int i = 0; i < Dimension; i++)
            {
                SendSetCentroid[mpi_rank*Dimension+i] /= mModelPart.GetCommunicator().LocalMesh().NumberOfElements();
                SetCentroid[mpi_rank*Dimension+i] /= mModelPart.GetCommunicator().LocalMesh().NumberOfElements();
            }
        }
        
        std::cout << "Center of partition at: " << SetCentroid[mpi_rank*Dimension+0] << " " << SetCentroid[mpi_rank*Dimension+1] << " " << SetCentroid[mpi_rank*Dimension+2] << std::endl;
        
        MPI_Allgather(&SendPartitionNumberOfElements[mpi_rank],1,MPI_INT,PartitionNumberOfElements,1,MPI_INT,MPI_COMM_WORLD);
        MPI_Allgather(&SendSetCentroid[mpi_rank*Dimension],Dimension,MPI_DOUBLE,SetCentroid,Dimension,MPI_DOUBLE,MPI_COMM_WORLD);
        
//         std::cout << "Weights\t ";
        
        for(int i = 0; i < mpi_size; i++)
        {    
            FulDistance[i] = 0;
            
            if(i != mpi_rank)
            {
                for(int j = 0; j < Dimension; j++)
                {
                    double aux = SetCentroid[i*Dimension+j] - SetCentroid[mpi_rank*Dimension+j];
                    FulDistance[i] += (aux*aux); //aka D
                }
                
                MidDistance[i] = FulDistance[i]/2;
                
                double factor = 0;
//                 double reduction = 0.5f;

                if (PartitionNumberOfElements[i] > LocalParticles)
                    factor = (double)(PartitionNumberOfElements[i] - LocalParticles) / (double)PartitionNumberOfElements[i];
                else
                    factor = (double)(LocalParticles - PartitionNumberOfElements[i]) / (double)LocalParticles;
                
                static int cak = 5;
                if (LastDistance[i] == 0 || cak < 5)
                {
                    LastDistance[i] = MidDistance[i];
                }
                cak++;
                
                double thisDistance = LastDistance[i] * factor;   
                double newDistance  = fabs(LastDistance[i] - thisDistance);
                
//                 std::cout << "LD: " << LastDistance[i] << " TD: " << thisDistance << " ND: " <<  newDistance << " MD: " << MidDistance[i] << std::endl;

                // No se si aixo es calcula aixi
                if (PartitionNumberOfElements[i] > LocalParticles)
                {
                    TheyWeight[i] = 0.5f;//(MidDistance[i] / newDistance) * reduction;
//                     SelfWeight[i] = 0.5f;//(MidDistance[i] / (FulDistance[i]-newDistance)) * reduction;
                }
                else
                {
//                     SelfWeight[i] = 0.5f;//(MidDistance[i] / newDistance) * reduction;
                    TheyWeight[i] = 0.5f;//(MidDistance[i] / (FulDistance[i]-newDistance)) * reduction;
                }
                
                std::cout << TheyWeight[i] << " ";
                
                LastDistance[i] = newDistance;
            }
        }
        
//         std::cout << std::endl;

        //Calculate the perpendicular plane who cuts the new partition
        //TODO: Just for fun. This isn't working anymore with the new algorithm. ha. ha. ha.
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
                    if(j != mpi_rank)
                    {         
                        distX = i_nod->X() - SetCentroid[j*Dimension+0];
                        distY = i_nod->Y() - SetCentroid[j*Dimension+1];
                        distZ = i_nod->Z() - SetCentroid[j*Dimension+2];
                        
                        distX *= distX;
                        distY *= distY;
                        distZ *= distZ;
                        
                        double TheyDist = (distX+distY+distZ);// * TheyWeight[j];
                        double SelfDist = (distX+distY+distZ);// * SelfWeight[j];
                        
                        if (TheyDist < min_dist)
                        {
                            min_dist = TheyDist;
                            PartitionIndex = j;
                        }
                        
                        if (SelfDist < min_dist)
                        {
                            min_dist = SelfDist;
                            PartitionIndex = mpi_rank;
                        }
                    }
                }
                
                (*particle_pointer_it)->GetValue(PARTITION_INDEX) = PartitionIndex;
                (i_nod)->GetSolutionStepValue(PARTITION_INDEX) = PartitionIndex;
            }
        }
          
        for (IteratorType particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it)
        {
            for (unsigned int i = 0; i < (*particle_pointer_it)->GetGeometry().PointsNumber(); i++)
            {
                int PartitionIndex = (*particle_pointer_it)->GetValue(PARTITION_INDEX);
                if(PartitionIndex != mpi_rank)
                {
                    SendObjects[PartitionIndex][NewNumberOfObjects[PartitionIndex]] = *particle_pointer_it;
                    NewNumberOfObjects[PartitionIndex]++;
                }
            }
        }

        for(int i = 0; i < mpi_size; i++)
        {
            std::cout << "(" << mModelPart.GetCommunicator().LocalMesh().NumberOfElements() << ") " << "Partition: " << mpi_rank << " --> " << i << ": " << NewNumberOfObjects[i] << std::endl;
            SendObjects[i].resize(NewNumberOfObjects[i]);
        }
        
        //Comunication here!
        CommunicationPhase(mModelPart,SendObjects,RecvObjects);
        Repartitionate(mModelPart,RecvObjects);
        
        pElements = mModelPart.GetCommunicator().LocalMesh().ElementsArray();
         
        ContainerType pElementsMarked;
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        //Mark elements in boundary (Level-1)
        for (IteratorType particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it)
        {         
            (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(INTERNAL_ENERGY) = 9999;
            
            (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(OSS_SWITCH) = 0; 
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
        
        if(CalculateBoundry)
        {
            //This matrix are used to perform the level-2 marking reduction
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
            
            double new_extension = 2.0;
            for (IteratorType particle_pointer_it = pElementsMarked.begin(); particle_pointer_it != pElementsMarked.end(); ++particle_pointer_it)
            {    
                  Radius[particle_pointer_it - pElementsMarked.begin()] = new_extension * ((1.0 +0.3) * (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS)); //if this is changed, then compobation before adding neighbours must change also.
            }
            
            particle_bin.SearchObjectsMpi(mModelPart,pElementsMarked.begin(),NumberOfMarkedElements,Radius,Results,ResultsDistances,NumberOfResults,MaximumNumberOfResults,mModelPart.pGetCommunicator());
            
            for (IteratorType particle_pointer_it = pElementsMarked.begin(); particle_pointer_it != pElementsMarked.end(); ++particle_pointer_it)
            {    
                int p_index = particle_pointer_it-pElementsMarked.begin();
                
                for(int i = 0; i < NumberOfResults[p_index]; i++)
                {
                    Results[p_index][i]->GetGeometry()(0)->GetSolutionStepValue(OSS_SWITCH) = (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(OSS_SWITCH);
                    SendMarkObjects[Results[p_index][i]->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX)].push_back(Results[p_index][i]);
                }
            }

            CommunicationPhase(mModelPart,SendMarkObjects,RecvMarkObjects);
            
            for (IteratorType particle_pointer_it = pElements.begin(); particle_pointer_it != pElements.end(); ++particle_pointer_it)
            {
                for(int j = 0; j < mpi_size; j++)
                {
                    if(j != mpi_rank)
                    {
                        for(int k = 0; k < RecvMarkObjects[j].size(); k++)
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
        
        
    }

    /// Destructor.
    virtual ~MortonParallelPartitioner() {
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

//     IteratorType mObjectsBegin;
//     IteratorType mObjectsEnd;
//     SizeType     mObjectsSize;

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
    MortonParallelPartitioner<TConfigure> & operator=(const MortonParallelPartitioner<TConfigure> & rOther)
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
    MortonParallelPartitioner(const MortonParallelPartitioner& rOther)
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
                                  MortonParallelPartitioner<TConfigure>& rThis)
{
    return rIStream;
}


/// output stream function
template<class TConfigure>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MortonParallelPartitioner<TConfigure> & rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


