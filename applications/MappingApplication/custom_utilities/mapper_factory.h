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

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

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

/// Python Interface of the MappingApplication
/** This class constructs the mappers and exposes them to Python
* Some checks are performed to see if the Input (ModelParts and JSON-Parameters) are valid
* Also the additional timing information is implemented here (echo_level = 1)
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
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
        mrJsonParameters(rJsonParameters)
    {
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

    void UpdateInterface(Kratos::Flags& rOptions, double SearchRadius)
    {
        double start_time = MapperUtilities::GetCurrentTime();
        mpMapper->UpdateInterface(rOptions, SearchRadius);
        double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

        mpMapper->pGetMapperCommunicator()->PrintTime(mMapperType,
                "UpdateInterface",
                elapsed_time);
    }


    /* This function maps a variable from Origin to Destination */
    void Map(const Variable<double>& rOriginVariable,
             const Variable<double>& rDestinationVariable,
             Kratos::Flags& rOptions)
    {
        double start_time = MapperUtilities::GetCurrentTime();
        mpMapper->Map(rOriginVariable, rDestinationVariable, rOptions);
        double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

        mpMapper->pGetMapperCommunicator()->PrintTime(mMapperType,
                "Map",
                elapsed_time);
    }

    /* This function maps a variable from Origin to Destination */
    void Map(const Variable< array_1d<double, 3> >& rOriginVariable,
             const Variable< array_1d<double, 3> >& rDestinationVariable,
             Kratos::Flags& rOptions)
    {
        double start_time = MapperUtilities::GetCurrentTime();
        mpMapper->Map(rOriginVariable, rDestinationVariable, rOptions);
        double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

        mpMapper->pGetMapperCommunicator()->PrintTime(mMapperType,
                "Map",
                elapsed_time);
    }


    /* This function maps from Destination to Origin */
    void InverseMap(const Variable<double>& rOriginVariable,
                    const Variable<double>& rDestinationVariable,
                    Kratos::Flags& rOptions)
    {
        double start_time = MapperUtilities::GetCurrentTime();
        mpMapper->InverseMap(rOriginVariable, rDestinationVariable, rOptions);
        double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

        mpMapper->pGetMapperCommunicator()->PrintTime(mMapperType,
                "InverseMap",
                elapsed_time);
    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                    const Variable< array_1d<double, 3> >& rDestinationVariable,
                    Kratos::Flags& rOptions)
    {
        double start_time = MapperUtilities::GetCurrentTime();
        mpMapper->InverseMap(rOriginVariable, rDestinationVariable, rOptions);
        double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

        mpMapper->pGetMapperCommunicator()->PrintTime(mMapperType,
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
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperFactory";
    }

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

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ReadAndCheckInterfaceModelParts()
    {
        int echo_level = 0;
        // read the echo_level temporarily, bcs the mrJsonParameters have not yet been validated and defaults assigned
        if (mrJsonParameters.Has("echo_level"))
        {
            echo_level = std::max(echo_level, mrJsonParameters["echo_level"].GetInt());
        }

        int comm_rank_origin = mrModelPartOrigin.GetCommunicator().MyPID();
        int comm_rank_destination = mrModelPartDestination.GetCommunicator().MyPID();

        if (mrJsonParameters.Has("interface_submodel_part_origin"))
        {
            std::string name_interface_submodel_part = mrJsonParameters["interface_submodel_part_origin"].GetString();
            mpInterfaceModelPartOrigin = &mrModelPartOrigin.GetSubModelPart(name_interface_submodel_part);

            if (echo_level >= 3 && comm_rank_origin == 0)
            {
                std::cout << "Mapper: SubModelPart used for Origin-ModelPart" << std::endl;
            }
        }
        else
        {
            mpInterfaceModelPartOrigin = &mrModelPartOrigin;

            if (echo_level >= 3 && comm_rank_origin == 0)
            {
                std::cout << "Mapper: Main ModelPart used for Origin-ModelPart" << std::endl;
            }
        }

        if (mrJsonParameters.Has("interface_submodel_part_destination"))
        {
            std::string name_interface_submodel_part = mrJsonParameters["interface_submodel_part_destination"].GetString();
            mpInterfaceModelPartDestination = &mrModelPartDestination.GetSubModelPart(name_interface_submodel_part);

            if (echo_level >= 3 && comm_rank_destination == 0)
            {
                std::cout << "Mapper: SubModelPart used for Destination-ModelPart" << std::endl;
            }
        }
        else
        {
            mpInterfaceModelPartDestination = &mrModelPartDestination;

            if (echo_level >= 3 && comm_rank_destination == 0)
            {
                std::cout << "Mapper: Main ModelPart used for Destination-ModelPart" << std::endl;
            }
        }

        const int num_nodes_origin = MapperUtilities::ComputeNumberOfNodes(*mpInterfaceModelPartOrigin);
        const int num_conditions_origin = MapperUtilities::ComputeNumberOfConditions(*mpInterfaceModelPartOrigin);
        const int num_elements_origin = MapperUtilities::ComputeNumberOfElements(*mpInterfaceModelPartOrigin);

        const int num_nodes_destination = MapperUtilities::ComputeNumberOfNodes(*mpInterfaceModelPartDestination);
        const int num_conditions_destination = MapperUtilities::ComputeNumberOfConditions(*mpInterfaceModelPartDestination);
        const int num_elements_destination = MapperUtilities::ComputeNumberOfElements(*mpInterfaceModelPartDestination);

        // Check if the ModelPart contains entities
        if (num_nodes_origin + num_conditions_origin + num_elements_origin < 1)
        {
            KRATOS_ERROR << "Neither Nodes nor Conditions nor Elements found "
                         << "in the Origin ModelPart" << std::endl;
        }

        if (num_conditions_origin > 0 && num_elements_origin > 0)
        {
            KRATOS_ERROR << "Origin ModelPart contains both Conditions and Elements "
                         << "which is not permitted" << std::endl;
        }

        if (num_nodes_destination + num_conditions_destination + num_elements_destination < 1)
        {
            KRATOS_ERROR << "Neither Nodes nor Conditions nor Elements found "
                         << "in the Destination ModelPart" << std::endl;
        }

        if (num_conditions_destination > 0 && num_elements_destination > 0)
        {
            KRATOS_ERROR << "Destination ModelPart contains both Conditions and Elements "
                         << "which is not permitted" << std::endl;
        }
    }

    void ConstructMapper()
    {
        double start_time = MapperUtilities::GetCurrentTime();

        if (!mrJsonParameters.Has("mapper_type"))
        {
            KRATOS_ERROR << "No \"mapper_type\" defined in json" << std::endl;
        }

        mMapperType = mrJsonParameters["mapper_type"].GetString();

        if (mMapperType == "NearestNeighbor")
        {
            if (mrJsonParameters.Has("approximation_tolerance"))
            {
                KRATOS_ERROR << "Invalid Parameter \"approximation_tolerance\" "
                             << "specified for Nearest Neighbor Mapper" << std::endl;
            }

            mpMapper = Mapper::Pointer(new NearestNeighborMapper(*mpInterfaceModelPartOrigin,
                                       *mpInterfaceModelPartDestination,
                                       mrJsonParameters));
        }
        else if (mMapperType == "NearestElement")
        {
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

          } */else
        {
            KRATOS_ERROR << "Selected Mapper \"" << mMapperType << "\" not implemented" << std::endl;
        }

        double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

        mpMapper->pGetMapperCommunicator()->PrintTime(mMapperType,
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