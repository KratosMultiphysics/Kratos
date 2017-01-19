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

#if !defined(KRATOS_MAPPER_FACTORY_H_INCLUDED )
#define  KRATOS_MAPPER_FACTORY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "containers/flags.h"
#include "includes/kratos_parameters.h"
#include "mapper_utilities.h"

#include "custom_utilities/nearest_neighbor_mapper.h"
#include "custom_utilities/nearest_element_mapper.h"
// #include "custom_utilities/approximate_mortar_mapper.h"


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
  class MapperFactory
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MapperFactory
      KRATOS_CLASS_POINTER_DEFINITION(MapperFactory);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      MapperFactory(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination,
                    Parameters& i_json_parameters) :
                    m_model_part_origin(i_model_part_origin),
                    m_model_part_destination(i_model_part_destination),
                    m_json_parameters(i_json_parameters) {

          CheckAndValidateJson();
          ReadAndCheckInterfaceModelParts();
          ConstructMapper();
      }

      /// Destructor.
      virtual ~MapperFactory() { }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void UpdateInterface(Kratos::Flags& options, double search_radius) {
          double start_time = MapperUtilities::GetCurrentTime();
          m_p_mapper->UpdateInterface(options, search_radius);
          double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

          m_p_mapper->GetMapperCommunicator()->PrintTime(m_mapper_type,
                                                         "UpdateInterface",
                                                         elapsed_time);
      }


      /* This function maps a variable from Origin to Destination */
      void Map(const Variable<double>& origin_variable,
               const Variable<double>& destination_variable,
               Kratos::Flags& options) {
          double start_time = MapperUtilities::GetCurrentTime();
          m_p_mapper->Map(origin_variable, destination_variable, options);
          double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

          m_p_mapper->GetMapperCommunicator()->PrintTime(m_mapper_type,
                                                         "Map",
                                                         elapsed_time);
      }

      /* This function maps a variable from Origin to Destination */
      void Map(const Variable< array_1d<double,3> >& origin_variable,
               const Variable< array_1d<double,3> >& destination_variable,
               Kratos::Flags& options) {
          double start_time = MapperUtilities::GetCurrentTime();
          m_p_mapper->Map(origin_variable, destination_variable, options);
          double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

          m_p_mapper->GetMapperCommunicator()->PrintTime(m_mapper_type,
                                                         "Map",
                                                         elapsed_time);
      }


      /* This function maps from Destination to Origin */
      void InverseMap(const Variable<double>& origin_variable,
                      const Variable<double>& destination_variable,
                      Kratos::Flags& options) {
          double start_time = MapperUtilities::GetCurrentTime();
          m_p_mapper->InverseMap(origin_variable, destination_variable, options);
          double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

          m_p_mapper->GetMapperCommunicator()->PrintTime(m_mapper_type,
                                                         "InverseMap",
                                                         elapsed_time);
      }

      /* This function maps from Destination to Origin */
      void InverseMap(const Variable< array_1d<double,3> >& origin_variable,
                      const Variable< array_1d<double,3> >& destination_variable,
                      Kratos::Flags& options) {
          double start_time = MapperUtilities::GetCurrentTime();
          m_p_mapper->InverseMap(origin_variable, destination_variable, options);
          double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

          m_p_mapper->GetMapperCommunicator()->PrintTime(m_mapper_type,
                                                         "InverseMap",
                                                         elapsed_time);
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
         buffer << "MapperFactory" ;
         return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "MapperFactory";}

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

      Mapper::Pointer m_p_mapper;
      std::string m_mapper_type;

      ModelPart& m_model_part_origin;
      ModelPart& m_model_part_destination;

      ModelPart* m_p_interface_model_part_origin;
      ModelPart* m_p_interface_model_part_destination;

      Parameters& m_json_parameters;
      Parameters m_default_settings = Parameters( R"(
      {
             "mapper_type"                           : "",
             "interface_submodel_part_origin"        : "",
             "interface_submodel_part_destination"   : "",
             "search_radius"                         : -1.0,
             "search_iterations"                     : 5,
             "non_conforming_interface"              : false,
             "echo_level"                            : 0
       }  )" );

       bool m_compute_search_radius;

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      void CheckAndValidateJson() {
          // Check if the three basic parameters are present
          if (!m_json_parameters.Has("mapper_type")) {
              KRATOS_ERROR << "MappingApplication; MapperFactory; no "
                           << "\"mapper_type\" defined in json" << std::endl;
          }

          if (!m_json_parameters.Has("interface_submodel_part_origin")) {
              KRATOS_ERROR << "MappingApplication; MapperFactory; no "
                           << "\"interface_submodel_part_origin\" "
                           << "defined in json" << std::endl;
          }

          if (!m_json_parameters.Has("interface_submodel_part_destination")) {
              KRATOS_ERROR << "MappingApplication; MapperFactory; no "
                           << "\"interface_submodel_part_destination\" "
                           << "defined in json" << std::endl;
          }

          // Check if there is a valid input for the search parameters
          m_compute_search_radius = true;
          if (m_json_parameters.Has("search_radius")) {
              m_compute_search_radius = false;

              if (m_json_parameters["search_radius"].GetDouble() < 0.0f) {
                  KRATOS_ERROR << "MappingApplication; MapperFactory; Invalid "
                               << "Search Radius specified" << std::endl;
              }
          }

          if (m_json_parameters.Has("search_iterations")) {
              if (m_json_parameters["search_iterations"].GetInt() < 1) {
                  KRATOS_ERROR << "MappingApplication; MapperFactory; Number of "
                               << "specified Search Iterations too small" << std::endl;
              }
          }

          m_json_parameters.RecursivelyValidateAndAssignDefaults(m_default_settings);
      }

      void ReadAndCheckInterfaceModelParts() {
        // TODO discuss pointer stuff with someone
        std::string name_interface_submodel_part = m_json_parameters["interface_submodel_part_origin"].GetString();
        m_p_interface_model_part_origin = &m_model_part_origin.GetSubModelPart(name_interface_submodel_part);

        name_interface_submodel_part = m_json_parameters["interface_submodel_part_destination"].GetString();
        m_p_interface_model_part_destination = &m_model_part_destination.GetSubModelPart(name_interface_submodel_part);

        if (MapperUtilities::ComputeNumberOfNodes(*m_p_interface_model_part_origin) < 1 &&
            MapperUtilities::ComputeNumberOfConditions(*m_p_interface_model_part_origin) < 1)
            KRATOS_ERROR << "MappingApplication; MapperFactory; Neither nodes nor "
                         << "conditions found in the origin model part" << std::endl;

        if (MapperUtilities::ComputeNumberOfNodes(*m_p_interface_model_part_destination) < 1 &&
            MapperUtilities::ComputeNumberOfConditions(*m_p_interface_model_part_destination) < 1)
            KRATOS_ERROR << "MappingApplication; MapperFactory; Neither nodes nor "
                         << "conditions found in the destination model part" << std::endl;

        int domain_size_origin = m_p_interface_model_part_origin->GetProcessInfo()[DOMAIN_SIZE];
        int domain_size_destination = m_p_interface_model_part_destination->GetProcessInfo()[DOMAIN_SIZE];

        if (domain_size_origin != domain_size_destination) {
            KRATOS_ERROR << "MappingApplication; MapperFactory; Domain sizes of the "
                         << "two model parts are not compatible" << std::endl;
        }

        // Compute the search radius in case it was not specified, can only be done after the modelparts are read
        if (m_compute_search_radius) {
          double search_radius = MapperUtilities::ComputeSearchRadius(*m_p_interface_model_part_origin,
                                                                      *m_p_interface_model_part_destination);
          m_json_parameters["search_radius"].SetDouble(search_radius);
        }
      }

      void ConstructMapper() {
          m_mapper_type = m_json_parameters["mapper_type"].GetString();

          double start_time = MapperUtilities::GetCurrentTime();

          if (m_mapper_type == "NearestNeighbor") {
              m_p_mapper = Mapper::Pointer(new NearestNeighborMapper(*m_p_interface_model_part_origin,
                                                                     *m_p_interface_model_part_destination,
                                                                     m_json_parameters));
          } else if (m_mapper_type == "NearestElement") {
              m_p_mapper = Mapper::Pointer(new NearestElementMapper(*m_p_interface_model_part_origin,
                                                                    *m_p_interface_model_part_destination,
                                                                    m_json_parameters));

          } /*else if (m_mapper_type == "Barycentric") {
              m_p_mapper = Mapper::Pointer(new BarycentricMapper(*m_p_interface_model_part_origin,
                                                                 *m_p_interface_model_part_destination,
                                                                 m_json_parameters));

          } *//*else if (m_mapper_type == "RBF") {
              m_p_mapper = Mapper::Pointer(new RBFMapper(*m_p_interface_model_part_origin,
                                                         *m_p_interface_model_part_destination,
                                                         m_json_parameters));

          } *//*else if (m_mapper_type == "ApproximateMortar") {
              m_p_mapper = Mapper::Pointer(new ApproximateMortarMapper(*m_p_interface_model_part_origin,
                                                                       *m_p_interface_model_part_destination,
                                                                       m_json_parameters));

          } *//*else if (m_mapper_type == "Mortar") {
              m_p_mapper = Mapper::Pointer(new MortarMapper(*m_p_interface_model_part_origin,
                                                            *m_p_interface_model_part_destination,
                                                            m_json_parameters));

          } *//*else if (m_mapper_type == "IGA") {
              m_p_mapper = Mapper::Pointer(new IGAMapper(*m_p_interface_model_part_origin,
                                                         *m_p_interface_model_part_destination,
                                                         m_json_parameters));

          } */else {
              KRATOS_ERROR << "MappingApplication; MapperFactory; "
                           << "\"ConstructMapper\" Selected Mapper "
                           << "not implemented" << std::endl;
          }

          double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

          m_p_mapper->GetMapperCommunicator()->PrintTime(m_mapper_type,
                                                         "Mapper Construction",
                                                         elapsed_time);

      }

      ///@}
      ///@name Private  Access
      ///@{


      ///@}
      ///@name Private Inquiry
      ///@{s


      ///@}
      ///@name Un accessible methods
      ///@{

      /// Assignment operator.
      MapperFactory& operator=(MapperFactory const& rOther);

    //   /// Copy constructor.
    //   MapperFactory(MapperFactory const& rOther){}


      ///@}

    }; // Class MapperFactory

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    MapperFactory& rThis)
    {
        return rIStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MapperFactory& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_FACTORY_H_INCLUDED  defined
