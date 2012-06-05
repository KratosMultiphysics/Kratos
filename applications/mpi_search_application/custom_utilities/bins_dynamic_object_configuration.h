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


// System includes

// External includes
#include <boost/python.hpp>

// #include "spatial_containers/spatial_containers.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"

// includes from other applications
#include "../../DEM_application/custom_utilities/spheric_particle_hertzian.h"

namespace Kratos
{
template <class TParticle>
class ParticleConfigure
{

public:

    ///@name Type Definitions
    ///@{
    enum {Dimension = 3};
    typedef TParticle                                                   ParticleType;
    typedef Point< 3, double>                                           PointType;
    typedef typename ParticleType::DistanceIteratorType                 DistanceIteratorType;
    typedef typename ParticleType::Pointer                              PointerType;
    typedef typename std::vector<typename ParticleType::Pointer>        ContainerType;
    typedef typename std::vector<PointerType>::iterator                 IteratorType;
    typedef ContainerType                                               ResultContainerType;
    typedef IteratorType                                                ResultIteratorType;

    /// Contact Pairs
    typedef ContactPair<PointerType>                                    ContactPairType;
    typedef  std::vector<ContactPairType>                               ContainerContactType;
    typedef  typename ContainerContactType::iterator                    IteratorContactType;
    typedef  typename ContainerContactType::value_type                  PointerContactType;

    /// Pointer definition of SpatialContainersConfigure
    KRATOS_CLASS_POINTER_DEFINITION(ParticleConfigure);

    ///@}
    ///@name Life Cycle
    ///@{

    ParticleConfigure() {};
    virtual ~ParticleConfigure() {}

    //******************************************************************************************************************

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
    {
        rLowPoint = *(rObject->GetPointerToCenterNode());
        rHighPoint = *(rObject->GetPointerToCenterNode());
        double radius = rObject->GetRadius();
        for(std::size_t i = 0; i < 3; i++)
        {
            rLowPoint[i]  += -radius;
            rHighPoint[i] += radius;
        }
    }

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint, const double& Radius)
    {
        rLowPoint = *(rObject->GetPointerToCenterNode());
        rHighPoint = *(rObject->GetPointerToCenterNode());
        for(std::size_t i = 0; i < 3; i++)
        {
            rLowPoint[i]  += -Radius;
            rHighPoint[i] += Radius;
        }
    }

    //******************************************************************************************************************

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
    {
        array_1d<double, 3> rObj_2_to_rObj_1 = rObj_1->GetPosition() - rObj_2->GetPosition();
        double distance_2 = rObj_2_to_rObj_1[0] * rObj_2_to_rObj_1[0] + rObj_2_to_rObj_1[1] * rObj_2_to_rObj_1[1] + rObj_2_to_rObj_1[2] * rObj_2_to_rObj_1[2];
        //distance_2 is the inter-center distance squared (from the definition of distance in search-structure.h, with operator (,))
        double radius_1 = rObj_1->GetRadius();
        double radius_2 = rObj_2->GetRadius();
        double radius_sum = radius_1 + radius_2;
        bool intersect = (distance_2 - radius_sum * radius_sum) <= 0;
        return intersect;
    }

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2, double Radius)
    {
        array_1d<double, 3> rObj_2_to_rObj_1 = rObj_1->GetPosition() - rObj_2->GetPosition();
        double distance_2 = rObj_2_to_rObj_1[0] * rObj_2_to_rObj_1[0] + rObj_2_to_rObj_1[1] * rObj_2_to_rObj_1[1] + rObj_2_to_rObj_1[2] * rObj_2_to_rObj_1[2];
        //distance_2 is the inter-center distance squared (from the definition of distance in search-structure.h, with operator (,))
        double radius_1 = Radius;//Cambien el radi del objecte de cerca per el gran, aixi no tindria que petar res
        double radius_2 = rObj_2->GetRadius();
        double radius_sum = radius_1 + radius_2;
        bool intersect = (distance_2 - radius_sum * radius_sum) <= 0;
        return intersect;
    }

    //******************************************************************************************************************

    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {
//        double separation_from_particle_radius_ratio = 0.1;
        array_1d<double, 3> center_of_particle = rObject->GetPosition();
        double radius = rObject->GetRadius();
        bool intersect = (rLowPoint[0] - radius <= center_of_particle[0] && rLowPoint[1] - radius <= center_of_particle[1] && rLowPoint[2] - radius <= center_of_particle[2] &&
                          rHighPoint[0] + radius >= center_of_particle[0] && rHighPoint[1] + radius >= center_of_particle[1] && rHighPoint[2] + radius >= center_of_particle[2]);
        return  intersect;
    }

    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint, const double& Radius)
    {
//        double separation_from_particle_radius_ratio = 0.1;
        array_1d<double, 3> center_of_particle = rObject->GetPosition();
        double radius = Radius;//Cambien el radi del objecte de cerca per el gran, aixi no tindria que petar res
        bool intersect = (rLowPoint[0] - radius <= center_of_particle[0] && rLowPoint[1] - radius <= center_of_particle[1] && rLowPoint[2] - radius <= center_of_particle[2] &&
                          rHighPoint[0] + radius >= center_of_particle[0] && rHighPoint[1] + radius >= center_of_particle[1] && rHighPoint[2] + radius >= center_of_particle[2]);
        return  intersect;
    }

    //******************************************************************************************************************

    static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance)
    {
        array_1d<double, 3> center_of_particle1 = rObj_1->GetPosition();
        array_1d<double, 3> center_of_particle2 = rObj_2->GetPosition();

        distance = sqrt((center_of_particle1[0] - center_of_particle2[0]) * (center_of_particle1[0] - center_of_particle2[0]) +
                        (center_of_particle1[1] - center_of_particle2[1]) * (center_of_particle1[1] - center_of_particle2[1]) +
                        (center_of_particle1[2] - center_of_particle2[2]) * (center_of_particle1[2] - center_of_particle2[2]) );
    }
    
        
    static inline void MPI_TransferResults(std::vector<std::vector<PointerType> >& remoteResults,
                                           std::vector<std::vector<PointerType> >& SearchResults,
                                           int * NumberOfSendObjects,
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
      
        MPI_Prepare_Communications(NumberOfSendObjects,NumberOfRecvPoints,msgSendSize,msgRecvSize);
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
    
    static inline void MPI_TransferParticles(std::vector<std::vector<PointerType> >& SendObjectToProcess,
                                             std::vector<std::vector<PointerType> >& SearchPetitions,
                                             int * NumberOfSendObjects,
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
            if(i != mpi_rank && NumberOfSendObjects[i])
            {
                Kratos::Serializer particleSerializer;
                particleSerializer.save("nodes",SendObjectToProcess[i]);
                
                serializer_buffer[i] = (std::stringstream *)particleSerializer.pGetBuffer();
                messages[i] = std::string(serializer_buffer[i]->str());
                msgSendSize[i] = serializer_buffer[i]->str().size();
            }
        }
        
        MPI_Prepare_Communications(NumberOfSendObjects,NumberOfRecvPoints,msgSendSize,msgRecvSize);
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

    virtual std::string Info() const
    {
        return " Spatial Containers Configure for Particles";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

protected:

private:

    /// Assignment operator.
    ParticleConfigure& operator=(ParticleConfigure const& rOther);

    /// Copy constructor.
    ParticleConfigure(ParticleConfigure const& rOther);

    ///@}

}; // Class ParticleConfigure


template <class TParticle>
inline std::istream& operator >> (std::istream& rIStream, ParticleConfigure<TParticle> & rThis)
{
    return rIStream;
}

template <class TParticle>
inline std::ostream& operator << (std::ostream& rOStream, const ParticleConfigure<TParticle>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
}

