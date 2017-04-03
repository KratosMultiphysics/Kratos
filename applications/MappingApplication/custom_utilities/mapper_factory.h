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
      MapperFactory(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination,
                    Parameters& rJsonParameters) :
                    mrModelPartOrigin(rModelPartOrigin),
                    mrModelPartDestination(rModelPartDestination),
                    mrJsonParameters(rJsonParameters) {

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

      void UpdateInterface(Kratos::Flags& rOptions, double SearchRadius) {
          double start_time = MapperUtilities::GetCurrentTime();
          mpMapper->UpdateInterface(rOptions, SearchRadius);
          double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

          mpMapper->GetMapperCommunicator()->PrintTime(mMapperType,
                                                         "UpdateInterface",
                                                         elapsed_time);
      }


      /* This function maps a variable from Origin to Destination */
      void Map(const Variable<double>& rOriginVariable,
               const Variable<double>& rDestinationVariable,
               Kratos::Flags& rOptions) {
          double start_time = MapperUtilities::GetCurrentTime();
          mpMapper->Map(rOriginVariable, rDestinationVariable, rOptions);
          double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

          mpMapper->GetMapperCommunicator()->PrintTime(mMapperType,
                                                         "Map",
                                                         elapsed_time);
      }

      /* This function maps a variable from Origin to Destination */
      void Map(const Variable< array_1d<double,3> >& rOriginVariable,
               const Variable< array_1d<double,3> >& rDestinationVariable,
               Kratos::Flags& rOptions) {
          double start_time = MapperUtilities::GetCurrentTime();
          mpMapper->Map(rOriginVariable, rDestinationVariable, rOptions);
          double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

          mpMapper->GetMapperCommunicator()->PrintTime(mMapperType,
                                                         "Map",
                                                         elapsed_time);
      }


      /* This function maps from Destination to Origin */
      void InverseMap(const Variable<double>& rOriginVariable,
                      const Variable<double>& rDestinationVariable,
                      Kratos::Flags& rOptions) {
          double start_time = MapperUtilities::GetCurrentTime();
          mpMapper->InverseMap(rOriginVariable, rDestinationVariable, rOptions);
          double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

          mpMapper->GetMapperCommunicator()->PrintTime(mMapperType,
                                                         "InverseMap",
                                                         elapsed_time);
      }

      /* This function maps from Destination to Origin */
      void InverseMap(const Variable< array_1d<double,3> >& rOriginVariable,
                      const Variable< array_1d<double,3> >& rDestinationVariable,
                      Kratos::Flags& rOptions) {
          double start_time = MapperUtilities::GetCurrentTime();
          mpMapper->InverseMap(rOriginVariable, rDestinationVariable, rOptions);
          double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

          mpMapper->GetMapperCommunicator()->PrintTime(mMapperType,
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

      Mapper::Pointer mpMapper;
      std::string mMapperType;

      ModelPart& mrModelPartOrigin;
      ModelPart& mrModelPartDestination;

      ModelPart* mpInterfaceModelPartOrigin;
      ModelPart* mpInterfaceModelPartDestination;

      Parameters& mrJsonParameters;
      Parameters mDefaultParameters = Parameters( R"(
      {
             "mapper_type"                           : "",
             "interface_submodel_part_origin"        : "",
             "interface_submodel_part_destination"   : "",
             "search_radius"                         : -1.0,
             "search_iterations"                     : 3,
             "approximation_tolerance"               : -1.0,
             "echo_level"                            : 0
       }  )" );

       bool mComputeSearchRadius;

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      void CheckAndValidateJson() {
          // Check if the three basic parameters are present
          if (!mrJsonParameters.Has("mapper_type")) {
              KRATOS_ERROR << "No \"mapper_type\" defined in json" << std::endl;
          }

          if (!mrJsonParameters.Has("interface_submodel_part_origin")) {
              KRATOS_ERROR << "No \"interface_submodel_part_origin\" "
                           << "defined in json" << std::endl;
          }

          if (!mrJsonParameters.Has("interface_submodel_part_destination")) {
              KRATOS_ERROR << "No \"interface_submodel_part_destination\" "
                           << "defined in json" << std::endl;
          }

          // Check if there is a valid input for the search parameters
          mComputeSearchRadius = true;
          if (mrJsonParameters.Has("search_radius")) {
              mComputeSearchRadius = false;

              if (mrJsonParameters["search_radius"].GetDouble() < 0.0f) {
                  KRATOS_ERROR << "Invalid Search Radius specified" << std::endl;
              }
          }

          if (mrJsonParameters.Has("search_iterations")) {
              if (mrJsonParameters["search_iterations"].GetInt() < 1) {
                  KRATOS_ERROR << "Number of specified Search Iterations too small" << std::endl;
              }
          }

          if (mrJsonParameters.Has("approximation_tolerance")) {
              if (mrJsonParameters["approximation_tolerance"].GetDouble() < 0.0f) {
                  KRATOS_ERROR << "Invalid Tolerance for Approximations specified" << std::endl;
              }
          }

          mrJsonParameters.RecursivelyValidateAndAssignDefaults(mDefaultParameters);
      }

      void ReadAndCheckInterfaceModelParts() {
        // TODO discuss pointer stuff with someone
        std::string name_interface_submodel_part = mrJsonParameters["interface_submodel_part_origin"].GetString();
        mpInterfaceModelPartOrigin = &mrModelPartOrigin.GetSubModelPart(name_interface_submodel_part);

        name_interface_submodel_part = mrJsonParameters["interface_submodel_part_destination"].GetString();
        mpInterfaceModelPartDestination = &mrModelPartDestination.GetSubModelPart(name_interface_submodel_part);

        if (MapperUtilities::ComputeNumberOfNodes(*mpInterfaceModelPartOrigin) < 1 &&
            MapperUtilities::ComputeNumberOfConditions(*mpInterfaceModelPartOrigin) < 1)
            KRATOS_ERROR << "Neither nodes nor conditions found "
                         << "in the origin model part" << std::endl;

        if (MapperUtilities::ComputeNumberOfNodes(*mpInterfaceModelPartDestination) < 1 &&
            MapperUtilities::ComputeNumberOfConditions(*mpInterfaceModelPartDestination) < 1)
            KRATOS_ERROR << "Neither nodes nor conditions found "
                         << "in the destination model part" << std::endl;

        // Compute the search radius in case it was not specified, can only be done after the modelparts are read
        if (mComputeSearchRadius) {
          double search_radius = MapperUtilities::ComputeSearchRadius(*mpInterfaceModelPartOrigin,
                                                                      *mpInterfaceModelPartDestination);
          mrJsonParameters["search_radius"].SetDouble(search_radius);
        }
      }

      void ConstructMapper() {
          mMapperType = mrJsonParameters["mapper_type"].GetString();

          double start_time = MapperUtilities::GetCurrentTime();

          if (mMapperType == "NearestNeighbor") {
              mpMapper = Mapper::Pointer(new NearestNeighborMapper(*mpInterfaceModelPartOrigin,
                                                                   *mpInterfaceModelPartDestination,
                                                                    mrJsonParameters));
          } else if (mMapperType == "NearestElement") {
              mpMapper = Mapper::Pointer(new NearestElementMapper(*mpInterfaceModelPartOrigin,
                                                                  *mpInterfaceModelPartDestination,
                                                                   mrJsonParameters));

          } /*else if (mMapperType == "Barycentric") {
              mpMapper = Mapper::Pointer(new BarycentricMapper(*mpInterfaceModelPartOrigin,
                                                                 *mpInterfaceModelPartDestination,
                                                                 mrJsonParameters));

          } *//*else if (mMapperType == "RBF") {
              mpMapper = Mapper::Pointer(new RBFMapper(*mpInterfaceModelPartOrigin,
                                                         *mpInterfaceModelPartDestination,
                                                         mrJsonParameters));

          } *//*else if (mMapperType == "Mortar") {
              mpMapper = Mapper::Pointer(new MortarMapper(*mpInterfaceModelPartOrigin,
                                                            *mpInterfaceModelPartDestination,
                                                            mrJsonParameters));

          } *//*else if (mMapperType == "IGA") {
              mpMapper = Mapper::Pointer(new IGAMapper(*mpInterfaceModelPartOrigin,
                                                         *mpInterfaceModelPartDestination,
                                                         mrJsonParameters));

          } */else {
              KRATOS_ERROR << "Selected Mapper not implemented" << std::endl;
          }

          double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

          mpMapper->GetMapperCommunicator()->PrintTime(mMapperType,
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