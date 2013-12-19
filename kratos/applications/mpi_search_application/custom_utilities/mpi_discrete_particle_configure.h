//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine  $
//   Date:                $Date: 2012-05-24 $
//   Revision:            $Revision: 1.0 $
//


#if !defined(KRATOS_MPI_DISCRETE_PARTICLE__CONFIGURE_INCLUDED)
#define  KRATOS_MPI_DISCRETE_PARTICLE__CONFIGURE_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <cmath>
#include "utilities/spatial_containers_configure.h"
#include "includes/mpi_communicator.h"
#include "mpi.h"
#include "spatial_containers/spatial_search.h"

// External includes
#include "custom_utilities/discrete_particle_configure.h"

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

    
template <std::size_t TDimension>
class MpiDiscreteParticleConfigure : public DiscreteParticleConfigure<TDimension> 
{
public:

      /// Pointer definition of SpatialContainersConfigure
      KRATOS_CLASS_POINTER_DEFINITION(MpiDiscreteParticleConfigure);

      typedef Communicator::MeshType    MeshType;
      typedef ModelPart::ElementsContainerType  ElementsContainerType;

      ///@}
      ///@name Life Cycle
      ///@{

      MpiDiscreteParticleConfigure(){};
      virtual ~MpiDiscreteParticleConfigure(){}

      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      //******************************************************************************************************************

      static inline void ReduceIds(int& total_elements,const int& offset, int& first_element)
      {
          int mpi_rank;
          int mpi_size;
          
          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
          
          std::vector<int> reduceArray(mpi_size+1);
          
          MPI_Allgather(&total_elements,1,MPI_INT,&reduceArray[1],1,MPI_INT,MPI_COMM_WORLD);
          
          reduceArray[0] = 1 + offset;
          for(int i = 1; i <= mpi_size; i++)
              reduceArray[i] += reduceArray[i-1] + 1;
          
          first_element = reduceArray[mpi_rank];
      }

      template<class TObjectType>                            
      static inline void TransferObjects(Communicator& Communicator,
                                         TObjectType& SendObjects,
                                         TObjectType& RecvObjects)
      {
          Communicator.TransferObjects(SendObjects,RecvObjects);
      }
      
      template<class TObjectType>
      static inline void TransferObjects(MeshType& ghostMesh,
                                         std::vector<TObjectType>& SendObjects,
                                         std::vector<TObjectType>& RecvObjects,
                                         VariablesList* pVariablesList
                                        )                                       
      {
          int mpi_rank;
          int mpi_size;
      
          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
          
          int msgSendSize[mpi_size];
          int msgRecvSize[mpi_size];
        
          for(int i = 0; i < mpi_size; i++)
          {
              msgSendSize[i] = 0;
              msgRecvSize[i] = 0;
          }
          
          for(int i = 0; i < mpi_size; i++)
          {
              msgSendSize[i] = SendObjects[i].size();
          }
  
          MPI_Alltoall(msgSendSize,1,MPI_INT,msgRecvSize,1,MPI_INT,MPI_COMM_WORLD);
    
          int NumberOfCommunicationEvents = 0;
          int NumberOfCommunicationEventsIndex = 0;
          
          int * mpi_recv_buffer[mpi_size];
          int * mpi_send_buffer[mpi_size];
          
          double * mpi_recv_buffer_x[mpi_size];
          double * mpi_send_buffer_x[mpi_size];
          
          double * mpi_recv_buffer_y[mpi_size];
          double * mpi_send_buffer_y[mpi_size];
          
          double * mpi_recv_buffer_z[mpi_size];
          double * mpi_send_buffer_z[mpi_size];
          
          for(int j = 0; j < mpi_size; j++)
          {
              if(j != mpi_rank && msgRecvSize[j]) NumberOfCommunicationEvents+=4;
              if(j != mpi_rank && msgSendSize[j]) NumberOfCommunicationEvents+=4;
          }
          
          MPI_Request reqs[NumberOfCommunicationEvents];
          MPI_Status stats[NumberOfCommunicationEvents];

          //Set up all receive and send events
          for(int i = 0; i < mpi_size; i++)
          {
              if(i != mpi_rank && msgRecvSize[i])
              {
                  mpi_recv_buffer[i] = (int *)malloc(sizeof(int) * msgRecvSize[i]);
                  
                  mpi_recv_buffer_x[i] = (double *)malloc(sizeof(double) * msgRecvSize[i]);
                  mpi_recv_buffer_y[i] = (double *)malloc(sizeof(double) * msgRecvSize[i]);
                  mpi_recv_buffer_z[i] = (double *)malloc(sizeof(double) * msgRecvSize[i]);
                  
                  MPI_Irecv(mpi_recv_buffer[i],msgRecvSize[i],MPI_INT,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
                  
                  MPI_Irecv(mpi_recv_buffer_x[i],msgRecvSize[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
                  MPI_Irecv(mpi_recv_buffer_y[i],msgRecvSize[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
                  MPI_Irecv(mpi_recv_buffer_z[i],msgRecvSize[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
              }

              if(i != mpi_rank && msgSendSize[i])
              {
                  mpi_send_buffer[i] = (int *)malloc(sizeof(int) * msgSendSize[i]);
                  
                  mpi_send_buffer_x[i] = (double *)malloc(sizeof(double) * msgSendSize[i]);
                  mpi_send_buffer_y[i] = (double *)malloc(sizeof(double) * msgSendSize[i]);
                  mpi_send_buffer_z[i] = (double *)malloc(sizeof(double) * msgSendSize[i]);
                  
                  for (ElementsContainerType::iterator i_element = SendObjects[i].begin();
                      i_element != SendObjects[i].end(); ++i_element)
                  {
                      int index = i_element-SendObjects[i].begin();
                    
                      mpi_send_buffer[i][index] = i_element->Id();
                      
                      mpi_send_buffer_x[i][index] = i_element->GetGeometry()(0)->Coordinates()[0];
                      mpi_send_buffer_y[i][index] = i_element->GetGeometry()(0)->Coordinates()[1];
                      mpi_send_buffer_z[i][index] = i_element->GetGeometry()(0)->Coordinates()[2];
                  }

                  MPI_Isend(mpi_send_buffer[i],msgSendSize[i],MPI_INT,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
                  
                  MPI_Isend(mpi_send_buffer_x[i],msgSendSize[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
                  MPI_Isend(mpi_send_buffer_y[i],msgSendSize[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
                  MPI_Isend(mpi_send_buffer_z[i],msgSendSize[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
              }
          }
          
          //wait untill all communications finish
          MPI_Waitall(NumberOfCommunicationEvents, reqs, stats);

          MPI_Barrier(MPI_COMM_WORLD);
          
          for(int i = 0; i < mpi_size; i++)
          { 
              if (i != mpi_rank && msgRecvSize[i])
              {
//                   RecvObjects[i].resize(msgRecvSize[i]);
                  for(int j = 0; j < msgRecvSize[i]; j++)
                  {
                      Geometry< Node < 3 > >::PointsArrayType nodelist;
                      Node < 3 > ::Pointer pnew_node;
                      pnew_node = Node < 3 > ::Pointer(new Node<3>(mpi_recv_buffer[i][j], 
                                                                   mpi_recv_buffer_x[i][j], 
                                                                   mpi_recv_buffer_y[i][j], 
                                                                   mpi_recv_buffer_z[i][j]));
                      
                      // Ask Pooyan: Clone the whole variable list?
                      pnew_node->SetSolutionStepVariablesList(pVariablesList);
                   
                      nodelist.push_back(pnew_node);
                    
                      Element::Pointer p_ghost_element = Element::Pointer(new Element(mpi_recv_buffer[i][j], nodelist));
                      
                      RecvObjects[i].push_back(p_ghost_element);
                  }
              }

              MPI_Barrier(MPI_COMM_WORLD);
          }
          
          // Free buffers
          for(int i = 0; i < mpi_size; i++)
          {
              if(mpi_rank != i && msgRecvSize[i])
              {
                  free(mpi_recv_buffer[i]);
                  
                  free(mpi_recv_buffer_x[i]);
                  free(mpi_recv_buffer_y[i]);
                  free(mpi_recv_buffer_z[i]);
              }
              
              if(i != mpi_rank && msgSendSize[i])
              {
                  free(mpi_send_buffer[i]);
                  
                  free(mpi_send_buffer_x[i]);
                  free(mpi_send_buffer_y[i]);
                  free(mpi_send_buffer_z[i]);
              }
          }
      }
      
      template<class TObjectType>
      static inline void TransferObjects(std::vector<TObjectType>& SendObjects,
                                         std::vector<TObjectType>& RecvObjects)                                       
      {
          // TODO: REWRITE ALL THIS CODE
          
          int mpi_rank;
          int mpi_size;
      
          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  
          std::stringstream * serializer_buffer;
          std::vector<std::string> buffer(mpi_size);
          
          int msgSendSize[mpi_size];
          int msgRecvSize[mpi_size];
        
          for(int i = 0; i < mpi_size; i++)
          {
              msgSendSize[i] = 0;
              msgRecvSize[i] = 0;
          }
        
          for(int i = 0; i < mpi_size; i++)
          {
              if(mpi_rank != i)
              {
                  Kratos::Serializer particleSerializer;
                  particleSerializer.save("nodes",SendObjects[i]);
                  
                  serializer_buffer = (std::stringstream *)particleSerializer.pGetBuffer();
                  buffer[i] = std::string(serializer_buffer->str());
                  msgSendSize[i] = buffer[i].size()+1;
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
              if((i != mpi_rank) && msgRecvSize[i])
              {
                  message[i] = (char *)malloc(sizeof(char) * msgRecvSize[i]);

                  MPI_Irecv(message[i],msgRecvSize[i],MPI_CHAR,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
              }

              if((i != mpi_rank) && msgSendSize[i])
              {
                  mpi_send_buffer[i] = (char *)malloc(sizeof(char) * msgSendSize[i]);
                  memcpy(mpi_send_buffer[i],buffer[i].c_str(),msgSendSize[i]);

                  MPI_Isend(mpi_send_buffer[i],msgSendSize[i],MPI_CHAR,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
              }
          }
          
          //wait untill all communications finish
          MPI_Waitall(NumberOfCommunicationEvents, reqs, stats);

          MPI_Barrier(MPI_COMM_WORLD);

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

              MPI_Barrier(MPI_COMM_WORLD);
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
     
    //******************************************************************************************************************

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
    virtual std::string Info() const {return " Spatial Containers Configure for Particles"; }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{
      

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

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

    /// Assignment operator.
    MpiDiscreteParticleConfigure& operator=(MpiDiscreteParticleConfigure const& rOther);

    /// Copy constructor.
    MpiDiscreteParticleConfigure(MpiDiscreteParticleConfigure const& rOther);

    ///@}

    }; // Class ParticleConfigure

    ///@}

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// input stream function
    template <std::size_t TDimension>
    inline std::istream& operator >> (std::istream& rIStream, MpiDiscreteParticleConfigure<TDimension> & rThis){
        return rIStream;
        }

    /// output stream function
    template <std::size_t TDimension>
    inline std::ostream& operator << (std::ostream& rOStream, const MpiDiscreteParticleConfigure<TDimension>& rThis){
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
        }
        
    ///@}

}   // namespace Kratos.
#endif	/* MPI_DISCRETE_PARTICLE_CONFIGURE_H */
