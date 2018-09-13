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

// Geometry includes
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"

// Spatial container includes
#include "utilities/spatial_containers_configure.h"

namespace Kratos
{
template <std::size_t TDimension>
class DiscreteParticleConfigure
{

public:

    ///@name Type Definitions
    ///@{
    enum { Dimension = TDimension,
                DIMENSION = TDimension,
                MAX_LEVEL = 16,
                MIN_LEVEL = 2
          };
      
    typedef  Point                        PointType;
    typedef  std::vector<double>::iterator                   DistanceIteratorType;
    typedef  ModelPart::ElementsContainerType::ContainerType ContainerType;
    typedef  ContainerType::value_type                       PointerType;
    typedef  ContainerType::iterator                         IteratorType;
    typedef  ModelPart::ElementsContainerType::ContainerType ResultContainerType;
    typedef  ResultContainerType::value_type                 ResultPointerType;
    typedef  ResultContainerType::iterator                   ResultIteratorType;
    typedef  ContactPair<PointerType>                        ContactPairType;
    typedef  std::vector<ContactPairType>                    ContainerContactType;
    typedef  ContainerContactType::iterator                  IteratorContactType;
    typedef  ContainerContactType::value_type                PointerContactType;
  
    typedef  std::vector<PointerType>::iterator              PointerTypeIterator;
    
    KRATOS_CLASS_POINTER_DEFINITION(DiscreteParticleConfigure);
    
    ///@}
    ///@name Life Cycle
    ///@{

    DiscreteParticleConfigure() {};
    virtual ~DiscreteParticleConfigure() {}

    //******************************************************************************************************************

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
    {
        rHighPoint = rLowPoint  = rObject->GetGeometry()[0];
        double radius = rObject->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);

        for(std::size_t i = 0; i < 3; i++)
        {
            rLowPoint[i]  += -radius;
            rHighPoint[i] += radius;
        }
    }

    static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint, const double& Radius)
    {
        rHighPoint = rLowPoint  = rObject->GetGeometry()[0];
        
        for(std::size_t i = 0; i < 3; i++)
        {
            rLowPoint[i]  += -Radius;
            rHighPoint[i] += Radius;
        }
    }

    //******************************************************************************************************************

    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
    {
        array_1d<double, 3> rObj_2_to_rObj_1 = rObj_1->GetGeometry()[0] - rObj_2->GetGeometry()[0];
        double distance_2 = inner_prod(rObj_2_to_rObj_1, rObj_2_to_rObj_1);
        //distance_2 is the inter-center distance squared (from the definition of distance in search-structure.h, with operator (,))
        const double& radius_1 = rObj_1->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
        const double& radius_2 = rObj_2->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
        double radius_sum      = radius_1 + radius_2;
        bool intersect         = (distance_2 - radius_sum * radius_sum) <= 0;
        
        return intersect;
    }


    static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2, double Radius)
    {
        array_1d<double, 3> rObj_2_to_rObj_1 = rObj_1->GetGeometry()[0] - rObj_2->GetGeometry()[0];
        double distance_2 = inner_prod(rObj_2_to_rObj_1, rObj_2_to_rObj_1);
        //distance_2 is the inter-center distance squared (from the definition of distance in search-structure.h, with operator (,))
        double radius_1 = Radius;//Cambien el radi del objecte de cerca per el gran, aixi no tindria que petar res
        const double& radius_2 = rObj_2->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
        double radius_sum = radius_1 + radius_2;
        bool intersect = (distance_2 - radius_sum * radius_sum) <= 0;
        
        return intersect;
    }

    //******************************************************************************************************************
    
    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint)
    {
//      double separation_from_particle_radius_ratio = 0.1;
        array_1d<double, 3> center_of_particle = rObject->GetGeometry()[0];
        const double& radius = rObject->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
        bool intersect = (rLowPoint[0] - radius <= center_of_particle[0] && rLowPoint[1] - radius <= center_of_particle[1] && rLowPoint[2] - radius <= center_of_particle[2] &&
            rHighPoint[0] + radius >= center_of_particle[0] && rHighPoint[1] + radius >= center_of_particle[1] && rHighPoint[2] + radius >= center_of_particle[2]);
        
        return  intersect;
    }

    static inline bool  IntersectionBox(const PointerType& rObject,  const PointType& rLowPoint, const PointType& rHighPoint, const double& Radius)
    {
//      double separation_from_particle_radius_ratio = 0.1;
        
        array_1d<double, 3> center_of_particle = rObject->GetGeometry()[0];
        double radius = Radius;//Cambien el radi del objecte de cerca per el gran, aixi no tindria que petar res
        bool intersect = (rLowPoint[0] - radius <= center_of_particle[0] && rLowPoint[1] - radius <= center_of_particle[1] && rLowPoint[2] - radius <= center_of_particle[2] &&
            rHighPoint[0] + radius >= center_of_particle[0] && rHighPoint[1] + radius >= center_of_particle[1] && rHighPoint[2] + radius >= center_of_particle[2]);
        
        return  intersect;
    }

    static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance)
    {
        array_1d<double, 3> center_of_particle1 = rObj_1->GetGeometry()[0];
        array_1d<double, 3> center_of_particle2 = rObj_2->GetGeometry()[0];

        distance = sqrt((center_of_particle1[0] - center_of_particle2[0]) * (center_of_particle1[0] - center_of_particle2[0]) +
                        (center_of_particle1[1] - center_of_particle2[1]) * (center_of_particle1[1] - center_of_particle2[1]) +
                        (center_of_particle1[2] - center_of_particle2[2]) * (center_of_particle1[2] - center_of_particle2[2]) );
    }
     
     //******************************************************************************************************************
    
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

    /// Turn back information as a string.
    virtual std::string Info() const {return " Spatial Containers Configure for Particles"; }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

protected:

private:

    /// Assignment operator.
    DiscreteParticleConfigure& operator=(DiscreteParticleConfigure const& rOther);

    /// Copy constructor.
    DiscreteParticleConfigure(DiscreteParticleConfigure const& rOther);

    ///@}

}; // Class ParticleConfigure
}

