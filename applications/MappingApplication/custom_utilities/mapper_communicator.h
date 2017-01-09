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
#include "interface_object_manager_serial.h"
#include "interface_search_structure.h"
#include "mapper_utilities.h"


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

      MapperCommunicator(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination,
                         Parameters& i_json_parameters) :
            m_model_part_origin(i_model_part_origin),
            m_model_part_destination(i_model_part_destination) {

          m_initial_search_radius = i_json_parameters["search_radius"].GetDouble();
          m_max_search_iterations = i_json_parameters["search_iterations"].GetInt();
          m_echo_level = i_json_parameters["echo_level"].GetInt();
      }

      /// Destructor.
      virtual ~MapperCommunicator() { }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      virtual void InitializeOrigin(MapperUtilities::InterfaceObjectConstructionType i_interface_object_type_origin,
                      GeometryData::IntegrationMethod i_integration_method_origin = GeometryData::NumberOfIntegrationMethods) {

          m_interface_object_manager_origin = InterfaceObjectManagerBase::Pointer(
              new InterfaceObjectManagerSerial(m_model_part_origin, MyPID(),
                                               TotalProcesses(),
                                               i_interface_object_type_origin,
                                               i_integration_method_origin,
                                               m_echo_level) );

          // Save for updating the interface
          m_interface_object_type_origin = i_interface_object_type_origin;
          m_integration_method_origin = i_integration_method_origin;
      }

      virtual void InitializeDestination(MapperUtilities::InterfaceObjectConstructionType i_interface_object_type_destination,
                      GeometryData::IntegrationMethod i_integration_method_destination = GeometryData::NumberOfIntegrationMethods) {

          m_interface_object_manager_destination = InterfaceObjectManagerBase::Pointer(
              new InterfaceObjectManagerSerial(m_model_part_destination, MyPID(),
                                               TotalProcesses(),
                                               i_interface_object_type_destination,
                                               i_integration_method_destination,
                                               m_echo_level) );

          // Save for updating the interface
          m_interface_object_type_destination = i_interface_object_type_destination;
          m_integration_method_destination = i_integration_method_destination;
      }

      void Initialize() {
          InitializeSearchStructure();
          InvokeSearch(m_initial_search_radius, m_max_search_iterations);
      }

      void UpdateInterface(Kratos::Flags& options, double i_initial_search_radius) {
          if (options.Is(MapperFlags::REMESHED)) { // recompute the managers and the search structure
              InitializeOrigin(m_interface_object_type_origin, m_integration_method_origin);
              InitializeDestination(m_interface_object_type_destination, m_integration_method_destination);
              InitializeSearchStructure();
              // TODO update number of nodes / conditions in mapper.h ...?
          } else { // clear the managers
              m_interface_object_manager_origin->Clear();
              m_interface_object_manager_destination->Clear();
          }

          if (i_initial_search_radius < 0.0f) {
              i_initial_search_radius = MapperUtilities::ComputeSearchRadius(m_model_part_origin,
                                                                             m_model_part_destination);
          }
          m_initial_search_radius = i_initial_search_radius; // update the search radius

          InvokeSearch(m_initial_search_radius, m_max_search_iterations);
      }



      // Interface function for mapper developers; scalar version
      virtual void TransferNodalData(const Variable<double>& origin_variable,
                                     const Variable<double>& destination_variable,
                                     Kratos::Flags& options,
                                     double factor = 1.0f) {
          TransferDataSerial(origin_variable, destination_variable,
                             options, factor);
      }

      // Interface function for mapper developers; vector version
      virtual void TransferNodalData(const Variable< array_1d<double,3> >& origin_variable,
                                     const Variable< array_1d<double,3> >& destination_variable,
                                     Kratos::Flags& options,
                                     double factor = 1.0f) {
          TransferDataSerial(origin_variable, destination_variable,
                             options, factor);
      }

      // Interface function for mapper developers; scalar version
      virtual void TransferInterpolatedData(const Variable<double>& origin_variable,
                                            const Variable<double>& destination_variable,
                                            Kratos::Flags& options,
                                            double factor = 1.0f) {
          options.Set(MapperFlags::INTERPOLATE_VALUES);
          TransferDataSerial(origin_variable, destination_variable,
                             options, factor);
      }

      // Interface function for mapper developers; vector version
      virtual void TransferInterpolatedData(const Variable< array_1d<double,3> >& origin_variable,
                                            const Variable< array_1d<double,3> >& destination_variable,
                                            Kratos::Flags& options,
                                            double factor = 1.0f) {
          options.Set(MapperFlags::INTERPOLATE_VALUES);
          TransferDataSerial(origin_variable, destination_variable,
                             options, factor);
      }

      // Interface function for mapper developer
      virtual void TransferShapeFunctions(Kratos::Flags& options) {


      }

      virtual int MyPID() // Copy from "kratos/includes/communicator.h"
      {
          return 0;
      }

      virtual int TotalProcesses() // Copy from "kratos/includes/communicator.h"
      {
          return 1;
      }

      InterfaceObjectManagerBase::Pointer GetInterfaceObjectManagerOrigin() {
          return m_interface_object_manager_origin;
      }

      InterfaceObjectManagerBase::Pointer GetInterfaceObjectManagerDestination() {
          return m_interface_object_manager_destination;
      }

      void PrintTime(const std::string& mapper_name,
                     const std::string& function_name,
                     const double& elapsed_time) {
          if (m_echo_level == 1 && MyPID() == 0) {
              std::cout  << "MAPPER TIMER: \"" << mapper_name << "\", \"" << function_name
                        << "\" took " <<  elapsed_time << " seconds" << std::endl;
          }
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
      virtual std::string Info() const {
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

      InterfaceObjectManagerBase::Pointer m_interface_object_manager_origin;
      InterfaceObjectManagerBase::Pointer m_interface_object_manager_destination;

      MapperUtilities::InterfaceObjectConstructionType m_interface_object_type_origin;
      GeometryData::IntegrationMethod m_integration_method_origin;
      MapperUtilities::InterfaceObjectConstructionType m_interface_object_type_destination;
      GeometryData::IntegrationMethod m_integration_method_destination;

      InterfaceSearchStructure::Pointer m_search_structure;

      double m_initial_search_radius;
      int m_max_search_iterations;

      int m_echo_level = 0;

      ///@}
      ///@name Protected Operators
      ///@{


      ///@}
      ///@name Protected Operations
      ///@{

      template <typename T>
      void TransferDataSerial(const Variable< T >& origin_variable,
                              const Variable< T >& destination_variable,
                              Kratos::Flags& options,
                              double factor) {
          if (options.Is(MapperFlags::SWAP_SIGN)) {
              factor *= (-1);
          }

          ExchangeDataLocal(origin_variable, destination_variable, options, factor);
      }

      template <typename T>
      void ExchangeDataLocal(const Variable< T >& origin_variable,
                             const Variable< T >& destination_variable,
                             Kratos::Flags& options,
                             const double factor) {

          std::vector< T > values;
          if (options.IsNot(MapperFlags::INVERSE_DIRECTION)) {
              m_interface_object_manager_origin->FillBufferWithValues(values, origin_variable, options);
              m_interface_object_manager_destination->ProcessValues(values, destination_variable, options, factor);
          } else {
              KRATOS_ERROR << "this direction is not yet implemented" << std::endl;
          }
      }

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

      virtual void InitializeSearchStructure() {
          m_search_structure = InterfaceSearchStructure::Pointer ( new InterfaceSearchStructure(
            m_interface_object_manager_destination, m_interface_object_manager_origin, m_echo_level) );
      }

      virtual void InvokeSearch(const double i_initial_search_radius,
                                const int i_max_search_iterations) {
          m_search_structure->Search(i_initial_search_radius, i_max_search_iterations);
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
