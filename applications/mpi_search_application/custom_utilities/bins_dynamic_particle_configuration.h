/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
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

#include "includes/serializer.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"

class ParticleSpatialConfigure
{
public:

    static const std::size_t Dim = 3;
    static const std::size_t Dimension = 3;

    typedef std::size_t  SizeType;
    typedef std::size_t  IndexType;

    typedef Kratos::Point<Dim>                  PointType;
    typedef Kratos::Tetrahedra3D4<PointType>    ObjectType;

    typedef PointType*                          PtrPointType;
    typedef ObjectType*                         PtrObjectType;

    typedef PtrPointType*                       PointVector;
    typedef PtrPointType*                       PointIterator;

    typedef PtrObjectType*                      ObjectVector;
    typedef PtrObjectType*                      ObjectIterator;

    typedef double*                             DistanceVector;
    typedef double*                             DistanceIterator;

    typedef Kratos::SearchUtils::SquaredDistanceFunction<Dim,PointType> DistanceFunction;
    
    static inline void MPI_TransferResults(std::vector<std::vector<PtrPointType> >& remoteResults,
                                           std::vector<std::vector<PtrPointType> >& SearchResults,
                                           int * NumberOfSendPoints,
                                           int * NumberOfRecvPoints,
                                           int * msgSendSize, 
                                           int * msgRecvSize
                                          ) 
    {
        int mpi_rank;
        int mpi_size;
      
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
      
        std::stringstream * serializer_buffer[mpi_size];
        std::string         messages[mpi_size];
        
        char * recvBuffers[mpi_size];
        
        for(int i = 0; i < mpi_size; i++)
        {
            if(mpi_rank != i && msgRecvSize[i])
            {
                Kratos::Serializer particleSerializer;
                particleSerializer.save("nodes",remoteResults[i]);
                  
                serializer_buffer[i] = (std::stringstream *)particleSerializer.pGetBuffer();
                messages[i] = std::string(serializer_buffer[i]->str());
                msgSendSize[i] = serializer_buffer[i]->str().size();
            }
        }
      
        MPI_Prepare_Communications(NumberOfSendPoints,NumberOfRecvPoints,msgSendSize,msgRecvSize);
        MPI_Async_SendAndRecive(messages,msgSendSize,msgRecvSize,recvBuffers);
        
        for (int i = 0; i < mpi_size; i++)
        { 
            if (i != mpi_rank)
            {
                Kratos::Serializer particleSerializer;
                serializer_buffer[0] = (std::stringstream *)particleSerializer.pGetBuffer();
                serializer_buffer[0]->write((char*)(recvBuffers[i]), msgRecvSize[i]);

                particleSerializer.load("nodes",SearchResults[i]);
            }
        }
    }
    
    static inline void MPI_TransferParticles(std::vector<std::vector<PtrPointType> >& SendPointToProcess,
                                             std::vector<std::vector<PtrPointType> >& SearchPetitions,
                                             int * NumberOfSendPoints,
                                             int * NumberOfRecvPoints,
                                             int * msgSendSize, 
                                             int * msgRecvSize
                                            )
    {
        int mpi_rank;
        int mpi_size;
      
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
      
        std::stringstream * serializer_buffer[mpi_size];
        std::string         messages[mpi_size];
        
        char * recvBuffers[mpi_size];

        //For each point and process fill the vector of send points
        for(int i = 0; i < mpi_size; i++)
        {
            if(i != mpi_rank && NumberOfSendPoints[i])
            {
                Kratos::Serializer particleSerializer;
                particleSerializer.save("nodes",SendPointToProcess[i]);
                
                serializer_buffer[i] = (std::stringstream *)particleSerializer.pGetBuffer();
                messages[i] = std::string(serializer_buffer[i]->str());
                msgSendSize[i] = serializer_buffer[i]->str().size();
            }
        }
        
        MPI_Prepare_Communications(NumberOfSendPoints,NumberOfRecvPoints,msgSendSize,msgRecvSize);
        MPI_Async_SendAndRecive(messages,msgSendSize,msgRecvSize,recvBuffers);
        
        for(int i = 0; i < mpi_size; i++) 
        {   
            if(i != mpi_rank && NumberOfRecvPoints[i])
            {                 
                Kratos::Serializer particleSerializer;
                serializer_buffer[0] = (std::stringstream *)particleSerializer.pGetBuffer();
                serializer_buffer[0]->write((char*)(recvBuffers[i]), msgRecvSize[i]);

                particleSerializer.load("nodes",SearchPetitions[i]);
            }
        }
    }
    
    static inline void MPI_Prepare_Communications(int * NumberOfSendElements,
                                                  int * NumberOfRecvElements,
                                                  int * msgSendSize,
                                                  int * msgRecvSize
                                                 )
    {
        int mpi_rank;
        int mpi_size;
      
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        
        //Message Size communication
        MPI_Alltoall(msgSendSize,1,MPI_INT,msgRecvSize,1,MPI_INT,MPI_COMM_WORLD);
        MPI_Alltoall(NumberOfSendElements,1,MPI_INT,NumberOfRecvElements,1,MPI_INT,MPI_COMM_WORLD);
    }
    
    static inline void MPI_Async_SendAndRecive(std::string messages[],
                                               int * msgSendSize,
                                               int * msgRecvSize,
                                               char * recvBuffers[]
                                              )
    {
        int mpi_rank;
        int mpi_size;
      
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        
        //Calculate number of communications
        int NumberOfCommunicationEvents = 0;
        int NumberOfCommunicationEventsIndex = 0;
        
        for(int j = 0; j < mpi_size; j++)
        {
            if(j != mpi_rank && msgRecvSize[j]) NumberOfCommunicationEvents++;
            if(j != mpi_rank && msgSendSize[j]) NumberOfCommunicationEvents++;
        }
        
        MPI_Request reqs[NumberOfCommunicationEvents];
        MPI_Status stats[NumberOfCommunicationEvents];
        
        //Set up all receive and send events
        for(int j = 0; j < mpi_size; j++)
        {
            if(j != mpi_rank && msgRecvSize[j])
            {
                recvBuffers[j] = (char *)malloc(sizeof(char) * msgRecvSize[j]);
                MPI_Irecv(recvBuffers[j],msgRecvSize[j],MPI_CHAR,j,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
            }

            if(j != mpi_rank && msgSendSize[j])
            {
                char * mpi_send_buffer = (char *)malloc(sizeof(char) * msgSendSize[j]);

                memcpy(mpi_send_buffer,messages[j].c_str(),msgSendSize[j]);
                MPI_Isend(mpi_send_buffer,msgSendSize[j],MPI_CHAR,j,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
            }
        }

        //wait untill all communications finish
        MPI_Waitall(NumberOfCommunicationEvents, reqs, stats);
    }
};

