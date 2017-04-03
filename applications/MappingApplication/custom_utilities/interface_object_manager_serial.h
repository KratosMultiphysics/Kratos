//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_INTERFACE_OBJECT_MANAGER_SERIAL_H_INCLUDED )
#define  KRATOS_INTERFACE_OBJECT_MANAGER_SERIAL_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "interface_object_manager_base.h"


namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

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
  class InterfaceObjectManagerSerial : public InterfaceObjectManagerBase
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of InterfaceObjectManagerSerial
      KRATOS_CLASS_POINTER_DEFINITION(InterfaceObjectManagerSerial);

      ///@}
      ///@name Life Cycle
      ///@{

      InterfaceObjectManagerSerial(ModelPart& rModelPart, int i_comm_rank, int i_comm_size,
                                   MapperUtilities::InterfaceObjectConstructionType i_interface_object_type,
                                   GeometryData::IntegrationMethod i_integration_method, const int i_echo_level,
                                   const double ApproximationTolerance) :
                                   InterfaceObjectManagerBase(
                                   rModelPart, i_comm_rank, i_comm_size, i_interface_object_type,
                                   i_integration_method, i_echo_level, ApproximationTolerance) { }

      /// Destructor.
      virtual ~InterfaceObjectManagerSerial() { }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      // **********************************************************************
      // Side we want to find neighbors for aka destination *******************
      // **********************************************************************
      void CheckResults() override {
          for (auto& interface_obj : m_interface_objects) {
              if (!interface_obj->NeighborFound()) {
                  std::cout << "MAPPER WARNING,\tPoint has not found a neighbor, [ "
                            << interface_obj->X() << " | " << interface_obj->Y() << " | "
                            << interface_obj->Z() << " ]" << std::endl;
              }
          }
      }

      void GetInterfaceObjectsSerialSearch(InterfaceObjectConfigure::ContainerType& candidate_send_objects) override {
          InitializeSizes();
          for (auto interface_obj : m_interface_objects) {
              if (!interface_obj->NeighborFound()) { // check if the interface object already found a neighbor
                  candidate_send_objects.push_back(interface_obj);
              }
          }
      }

      void PostProcessReceivedResults(const std::vector<double>& distances,
                                      const InterfaceObjectConfigure::ContainerType& candidate_send_objects) override {
          int i = 0;
          for (auto interface_obj : candidate_send_objects) {
              if (distances[i] > -0.5f) { // failed search has value "-1"
                  interface_obj->ProcessDistance(distances[i], m_comm_rank);
                  m_send_objects[m_comm_rank].push_back(interface_obj);
              }
              ++i;
          }
      }

      // **********************************************************************
      // Side where we search neighbors aka origin ****************************
      // **********************************************************************
      void StoreSearchResults(const std::vector<double>& distances,
                              const std::vector<InterfaceObject::Pointer> temp_closest_results,
                              const std::vector<std::vector<double>> temp_shape_functions,
                              const std::vector<array_1d<double,2>> temp_local_coordinates) override {
          for (std::size_t i = 0; i < distances.size(); ++i) {
              if (distances[i] > -0.5f) { // failed search has value "-1"
                  m_receive_objects[m_comm_rank].push_back(temp_closest_results[i]);
                  m_local_coordinates[m_comm_rank].push_back(temp_local_coordinates[i]);
                  m_shape_functions[m_comm_rank].push_back(temp_shape_functions[i]);
              }
          }
      }

      std::vector<InterfaceObject::Pointer>& GetDestinationInterfaceObjects() override {
          return m_send_objects[m_comm_rank];
      }

      std::vector<InterfaceObject::Pointer>& GetOriginInterfaceObjects() override {
          return m_receive_objects[m_comm_rank];
      }

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
      virtual std::string Info() const
      {
	      std::stringstream buffer;
        buffer << "InterfaceObjectManagerSerial" ;
        return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "InterfaceObjectManagerSerial";}

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

      void InitializeSizes() {
          int size = m_interface_objects.size();
          m_send_objects[m_comm_rank].reserve(size);
          m_receive_objects[m_comm_rank].reserve(size);
          m_shape_functions[m_comm_rank].reserve(size);
          m_local_coordinates[m_comm_rank].reserve(size);
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

      /// Assignment operator.
      InterfaceObjectManagerSerial& operator=(InterfaceObjectManagerSerial const& rOther);

    //   /// Copy constructor.
    //   InterfaceObjectManagerSerial(InterfaceObjectManagerSerial const& rOther){}


      ///@}

    }; // Class InterfaceObjectManagerSerial

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    InterfaceObjectManagerSerial& rThis)
    {
        return rIStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const InterfaceObjectManagerSerial& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_OBJECT_MANAGER_SERIAL_H_INCLUDED  defined
