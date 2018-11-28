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
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

#include "custom_mappers/mapper.h"

namespace Kratos
{
///@addtogroup MappingApplication
///@{

///@name Kratos Classes
///@{

/// Python Interface of the MappingApplication
/** This class constructs the mappers and exposes them to Python
* Some checks are performed to see if the Input (ModelParts and JSON-Parameters) are valid
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

    /// Destructor.
    virtual ~MapperFactory() = default;

    ///@}
    ///@name Operations
    ///@{


    // template<class TSparseSpace, class TDenseSpace>
    // static typename Mapper<TSparseSpace, TDenseSpace>::Pointer CreateMapper(
    //     ModelPart& rModelPartOrigin,
    //     ModelPart& rModelPartDestination,
    //     Parameters MapperSettings);

    // template<class TSparseSpace, class TDenseSpace>
    // static void Register(const std::string& rMapperName,
    //               typename Mapper<TSparseSpace, TDenseSpace>::Pointer pMapperPrototype);

    template<class TSparseSpace, class TDenseSpace>
    static typename Mapper<TSparseSpace, TDenseSpace>::Pointer CreateMapper(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters MapperSettings)
    {
        ModelPart& r_interface_model_part_origin = ReadInterfaceModelPart(rModelPartOrigin, MapperSettings, "origin");
        ModelPart& r_interface_model_part_destination = ReadInterfaceModelPart(rModelPartDestination, MapperSettings, "destination");

        const std::string mapper_name = MapperSettings["mapper_type"].GetString();

        const auto& mapper_list = GetRegisteredMappersList<TSparseSpace, TDenseSpace>();

        if (mapper_list.find(mapper_name) != mapper_list.end()) {
            // Removing Parameters that are not needed by the Mapper
            MapperSettings.RemoveValue("mapper_type");
            MapperSettings.RemoveValue("interface_submodel_part_origin");
            MapperSettings.RemoveValue("interface_submodel_part_destination");

            // TODO check why this works, Clone currently returns a unique ptr!!!
            return mapper_list.at(mapper_name)->Clone(r_interface_model_part_origin,
                                                      r_interface_model_part_destination,
                                                      MapperSettings);
        }
        else {
            std::stringstream err_msg;
            err_msg << "The requested Mapper \"" << mapper_name <<"\" is not not available!\n"
                    << "The following Mappers are available:" << std::endl;

            for (auto const& registered_mapper : mapper_list)
                err_msg << "\t" << registered_mapper.first << "\n";

            KRATOS_ERROR << err_msg.str() << std::endl;
        }
    }

    template<class TSparseSpace, class TDenseSpace>
    static void Register(const std::string& rMapperName,
                  typename Mapper<TSparseSpace, TDenseSpace>::Pointer pMapperPrototype)
    {
        GetRegisteredMappersList<TSparseSpace, TDenseSpace>().insert(
            std::make_pair(rMapperName, pMapperPrototype));
    }


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

private:
    ///@name Private Operations
    ///@{

    /// Default constructor.
    MapperFactory() {}

    static ModelPart& ReadInterfaceModelPart(ModelPart& rModelPart,
                                             Parameters InterfaceParameters,
                                             const std::string& InterfaceSide);

    // template<class TSparseSpace, class TDenseSpace>
    // static std::unordered_map<std::string, typename Mapper<TSparseSpace,
    //     TDenseSpace>::Pointer>& GetRegisteredMappersList();

    template<class TSparseSpace, class TDenseSpace>
    static std::unordered_map<std::string, typename Mapper<TSparseSpace,
        TDenseSpace>::Pointer>& GetRegisteredMappersList()
    {
        static std::unordered_map<std::string, typename Mapper<TSparseSpace, TDenseSpace>::Pointer> registered_mappers;

        return registered_mappers;
    }

    ///@}

}; // Class MapperFactory

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/* typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> UblasDenseSpaceType;

// template<>
// inline MPI_Datatype MapperUtilitiesMPI::GetMPIDatatype<int>(const int& rValue)
// {
//     return MPI_INT ;
// }

template<>
std::unordered_map<std::string, typename Mapper<SparseSpaceType,
    UblasDenseSpaceType>::Pointer>& MapperFactory::GetRegisteredMappersList<SparseSpaceType, UblasDenseSpaceType>()
{
    static std::unordered_map<std::string, typename Mapper<SparseSpaceType, UblasDenseSpaceType>::Pointer> registered_mappers;

    return registered_mappers;
} */

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