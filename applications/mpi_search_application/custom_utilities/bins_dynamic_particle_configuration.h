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

    typedef Kratos::Point                  PointType;
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
    
    static inline void Save(std::vector<PtrPointType>& inputElements, std::string& buffer) 
    {
        Kratos::Serializer particleSerializer;
        std::stringstream * serializer_buffer;
        
        particleSerializer.save("nodes",inputElements);
          
        serializer_buffer = (std::stringstream *)particleSerializer.pGetBuffer();
        buffer = std::string(serializer_buffer->str());
    }

    static inline void Load(std::vector<PtrPointType>& outputElements, std::string& buffer) 
    {
        Kratos::Serializer particleSerializer;
        std::stringstream * serializer_buffer;
        
        serializer_buffer = (std::stringstream *)particleSerializer.pGetBuffer();
        serializer_buffer->write((char*)(buffer.c_str()), buffer.size());

        particleSerializer.load("nodes",outputElements);
    }
    
    template<class TObjectType>
    static inline void AsyncSendAndRecive(std::vector<TObjectType>& SendObjects,
                                          std::vector<TObjectType>& RecvObjects,
                                          int * msgSendSize,
                                          int * msgRecvSize)                                       
    {
        int mpi_rank;
        int mpi_size;
      
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  
        std::stringstream * serializer_buffer;
        std::string buffer[mpi_size];
        
        for(int i = 0; i < mpi_size; i++)
        {
            if(mpi_rank != i)
            {
                Kratos::Serializer particleSerializer;
                particleSerializer.save("nodes",SendObjects[i]);

                serializer_buffer = (std::stringstream *)particleSerializer.pGetBuffer();
                buffer[i] = std::string(serializer_buffer->str());
                msgSendSize[i] = buffer[i].size();
            }
        }
  
        MPI_Alltoall(msgSendSize,1,MPI_INT,msgRecvSize,1,MPI_INT,MPI_COMM_WORLD);
  
        int NumberOfCommunicationEvents = 0;
        int NumberOfCommunicationEventsIndex = 0;
        
        char * message[mpi_size];
        char * mpi_send_buffer[mpi_size];
        
        for(int j = 0; j < mpi_size; j++)
        {
            if(j != mpi_rank && msgRecvSize[j]) NumberOfCommunicationEvents++;
            if(j != mpi_rank && msgSendSize[j]) NumberOfCommunicationEvents++;
        }
        
        MPI_Request reqs[NumberOfCommunicationEvents];
        MPI_Status stats[NumberOfCommunicationEvents];

        //Set up all receive and send events
        for(int i = 0; i < mpi_size; i++)
        {
            if(i != mpi_rank && msgRecvSize[i])
            {
                message[i] = (char *)malloc(sizeof(char) * msgRecvSize[i]);

                MPI_Irecv(message[i],msgRecvSize[i],MPI_CHAR,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
            }

            if(i != mpi_rank && msgSendSize[i])
            {
                mpi_send_buffer[i] = (char *)malloc(sizeof(char) * msgSendSize[i]);

                memcpy(mpi_send_buffer[i],buffer[i].c_str(),msgSendSize[i]);
                MPI_Isend(mpi_send_buffer[i],msgSendSize[i],MPI_CHAR,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
            }
        }

        //wait untill all communications finish
        MPI_Waitall(NumberOfCommunicationEvents, reqs, stats);

        for(int i = 0; i < mpi_size; i++)
        { 
            if (i != mpi_rank && msgRecvSize[i])
            {
                Kratos::Serializer particleSerializer;
                std::stringstream * serializer_buffer;
                
                serializer_buffer = (std::stringstream *)particleSerializer.pGetBuffer();
                serializer_buffer->write(message[i], msgRecvSize[i]);
  
                particleSerializer.load("nodes",RecvObjects[i]);

            }
        }
        
        // Free buffers
        for(int i = 0; i < mpi_size; i++)
        {
            if(mpi_rank != i && msgRecvSize[i])
                free(message[i]);
            
            if(i != mpi_rank && msgSendSize[i])
                free(mpi_send_buffer[i]);
        }

    }
    
    
};

