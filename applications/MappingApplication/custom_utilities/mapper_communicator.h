//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher

#if !defined(KRATOS_MAPPER_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_MAPPER_COMMUNICATOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "interface_object_manager.h"
#include "interface_search_structure.h"


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
  class MapperCommunicator
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MapperCommunicator
      KRATOS_CLASS_POINTER_DEFINITION(MapperCommunicator);

      ///@}
      ///@name Life Cycle
      ///@{

      MapperCommunicator(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination) :
            m_model_part_origin(i_model_part_origin), m_model_part_destination(i_model_part_destination) {

      }

      /// Destructor.
      virtual ~MapperCommunicator() {

      }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      virtual void Initialize(InterfaceObjectManager::Pointer i_interface_object_manager_origin,
                              InterfaceObjectManager::Pointer i_interface_object_manager_destination,
                              bool bi_directional_search) {

          m_interface_object_manager_origin = i_interface_object_manager_origin;
          m_interface_object_manager_destination = i_interface_object_manager_destination;
          m_bidirectional_search = bi_directional_search;

          m_search_structure_map = InterfaceSearchStructure::Pointer ( new InterfaceSearchStructure(
             i_interface_object_manager_destination, i_interface_object_manager_origin,
             m_model_part_origin) );

          if (bi_directional_search) {
              m_search_structure_inverse_map = InterfaceSearchStructure::Pointer ( new InterfaceSearchStructure(
                 i_interface_object_manager_origin, i_interface_object_manager_destination,
                 m_model_part_destination) );
          }

          ComputeSearchStructure();
      }

      virtual void Initialize(InterfaceObjectManager::Pointer i_interface_object_manager_origin,
                              InterfaceObjectManager::Pointer i_interface_object_manager_destination,
                              bool bi_directional_search, double i_initial_search_radius,
                              int i_max_search_iterations) {

          m_interface_object_manager_origin = i_interface_object_manager_origin;
          m_interface_object_manager_destination = i_interface_object_manager_destination;
          m_bidirectional_search = bi_directional_search;

          m_search_structure_map = InterfaceSearchStructure::Pointer ( new InterfaceSearchStructure(
             i_interface_object_manager_destination, i_interface_object_manager_origin,
             i_initial_search_radius, i_max_search_iterations) );

          if (bi_directional_search) {
              m_search_structure_inverse_map = InterfaceSearchStructure::Pointer ( new InterfaceSearchStructure(
                 i_interface_object_manager_origin, i_interface_object_manager_destination,
                 i_initial_search_radius, i_max_search_iterations) );
          }

          ComputeSearchStructure();
      }

      virtual void ComputeSearchStructure() {
          m_search_structure_map->Clear();
          m_search_structure_map->Search();
          if (m_bidirectional_search)
              m_search_structure_inverse_map ->Clear();
              m_search_structure_inverse_map ->Search();
      }

      // Interface function for mapper developers; scalaer version
      virtual void TransferData(bool direction, const Variable<double>& origin_variable,
                                const Variable<double>& destination_variable,
                                const bool add_value, // TODO change order => variables first
                                bool use_other_search_structure = false) {
          TransferDataSerial(origin_variable, destination_variable,
                             direction, add_value, use_other_search_structure);
      }

      // Interface function for mapper developers; vector version
      virtual void TransferData(bool direction, const Variable< array_1d<double,3> >& origin_variable,
                                const Variable< array_1d<double,3> >& destination_variable,
                                const bool add_value,
                                bool use_other_search_structure = false) {
          TransferDataSerial(origin_variable, destination_variable,
                             direction, add_value, use_other_search_structure);
      }

      template <typename T>
      void TransferDataSerial(const Variable< T >& origin_variable,
                              const Variable< T >& destination_variable,
                              const bool direction, const bool add_value,
                              const bool use_other_search_structure) {
          if ((!m_bidirectional_search && !direction && !use_other_search_structure) ||
              (!m_bidirectional_search && direction && use_other_search_structure))
              KRATOS_ERROR << "MappingApplication; MapperCommunicator; \"TransferDataSerial\" Bi-Directional Search not defined!" << std::endl;

          // TODO here some implicit conversion is done, is that ok?
          std::vector<InterfaceObject::Pointer> point_list_1;
          std::vector<InterfaceObject::Pointer> point_list_2;

          SelectPointLists(point_list_1, point_list_2, direction, use_other_search_structure);

          ExchangeDataLocal(origin_variable, destination_variable,
                            point_list_1, point_list_2, add_value);
      }

      virtual void SelectPointLists(std::vector<InterfaceObject::Pointer>& point_list_1,
                            std::vector<InterfaceObject::Pointer>& point_list_2,
                            const bool direction, const bool use_other_search_structure) {
          // boolean "use_other_search_structure" is intended to be used for mortar mapping
          // to exchange data in both directions with the same match-partners
          if (direction) { // map; origin -> destination
              if (use_other_search_structure) { // use search structure of InverseMap
                  point_list_1 = m_interface_object_manager_origin->GetPointList();
                  point_list_2 = m_interface_object_manager_destination->GetPointListBins();
              } else { // use search structure of Map
                  point_list_1 = m_interface_object_manager_origin->GetPointListBins();
                  point_list_2 = m_interface_object_manager_destination->GetPointList();
              }
          } else { // inverse map; destination -> origin
              if (use_other_search_structure) { // use search structure of Map
                  point_list_1 = m_interface_object_manager_destination->GetPointList();
                  point_list_2 = m_interface_object_manager_origin->GetPointListBins();
              } else { // use search structure of InverseMap
                  point_list_1 = m_interface_object_manager_destination->GetPointListBins();
                  point_list_2 = m_interface_object_manager_origin->GetPointList();
              }
          }
          if (point_list_1.size() != point_list_2.size())
              KRATOS_ERROR << "MappingApplication; MapperCommunicator; \"SelectPointLists\" Size Mismatch!" << std::endl;
      }

      template <typename T>
      void ExchangeDataLocal(const Variable< T >& origin_variable,
                        const Variable< T >& destination_variable,
                        std::vector<InterfaceObject::Pointer>& point_list_1,
                        std::vector<InterfaceObject::Pointer>& point_list_2,
                        const bool add_value) {
          for (std::size_t i = 0; i < point_list_1.size(); ++i) {
              T data = point_list_1[i]->GetObjectValue(origin_variable);
              if (add_value) {
                  point_list_2[i]->AddObjectValue(destination_variable, data);
              } else {
                  point_list_2[i]->SetObjectValue(destination_variable, data);
              }
          }
      }

      virtual int MyPID() // Copy from "kratos/includes/communicator.h"
      {
          return 0;
      }

      virtual int TotalProcesses() // Copy from "kratos/includes/communicator.h"
      {
          return 1;
      }

      InterfaceObjectManager::Pointer GetInterfaceObjectManagerOrigin() {
          return m_interface_object_manager_origin;
      }

      InterfaceObjectManager::Pointer GetInterfaceObjectManagerDestination() {
          return m_interface_object_manager_destination;
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
        buffer << "MapperCommunicator" ;
        return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "MapperCommunicator";}

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
      ModelPart& m_model_part_origin;
      ModelPart& m_model_part_destination;

      InterfaceObjectManager::Pointer m_interface_object_manager_origin;
      InterfaceObjectManager::Pointer m_interface_object_manager_destination;

      InterfaceSearchStructure::Pointer m_search_structure_map;
      InterfaceSearchStructure::Pointer m_search_structure_inverse_map;

      bool m_bidirectional_search; // TODO from Mapper or from *.json?

      // int m_omp_threshold_num_nodes = 1000; // TODO constexpr???

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
      MapperCommunicator& operator=(MapperCommunicator const& rOther);

    //   /// Copy constructor.
    //   MapperCommunicator(MapperCommunicator const& rOther){}


      ///@}

    }; // Class MapperCommunicator

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    MapperCommunicator& rThis)
{
    return rIStream;
}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MapperCommunicator& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_COMMUNICATOR_H_INCLUDED  defined
