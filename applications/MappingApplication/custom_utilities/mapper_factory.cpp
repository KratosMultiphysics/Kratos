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


// System includes


// External includes 


// Project includes
#include "mapper_factory.h"


namespace Kratos
{
    Mapper::Pointer MapperFactory::CreateMapper(ModelPart& rModelPartOrigin, 
                                                ModelPart& rModelPartDestination,
                                                Parameters JsonParameters)
    {
        ModelPart& r_interface_model_part_origin = ReadInterfaceModelPart(rModelPartOrigin, JsonParameters, "origin");
        ModelPart& r_interface_model_part_destination = ReadInterfaceModelPart(rModelPartDestination, JsonParameters, "destination");

        const std::string mapper_name = JsonParameters["mapper_type"].GetString();

        const auto& mapper_list = MapperFactory::GetRegisteredMappersList();

        if (mapper_list.find(mapper_name) != mapper_list.end())
        {
            return mapper_list.at(mapper_name)->Clone(r_interface_model_part_origin,
                                                      r_interface_model_part_destination,
                                                      JsonParameters);
        }
        else
        {
            std::stringstream err_msg;
            err_msg << "The requested Mapper \"" << mapper_name <<"\" is not registered! The following mappers are registered:" << std::endl;
            for (auto const& registered_mapper : mapper_list)
            {
                err_msg << registered_mapper.first << ", ";
            }
            KRATOS_ERROR << err_msg.str() << std::endl;  
        }
    }

    void MapperFactory::Register(const std::string MapperName,
                                 Mapper::Pointer pMapperPrototype)
    {
        GetRegisteredMappersList().insert(make_pair(MapperName, pMapperPrototype));
    }

    ModelPart& MapperFactory::ReadInterfaceModelPart(ModelPart& rModelPart,
                                                     Parameters InterfaceParameters,
                                                     const std::string& InterfaceSide)
    {
        int echo_level = 0;
        // read the echo_level temporarily, bcs the mJsonParameters have not yet been validated and defaults assigned
        if (InterfaceParameters.Has("echo_level"))
        {
            echo_level = std::max(echo_level, InterfaceParameters["echo_level"].GetInt());
        }

        int comm_rank = rModelPart.GetCommunicator().MyPID();

        std::string key_sub_model_part = "interface_submodel_part_";
        key_sub_model_part.append(InterfaceSide);

        if (InterfaceParameters.Has(key_sub_model_part))
        {
            const std::string name_interface_submodel_part = InterfaceParameters[key_sub_model_part].GetString();
            
            if (echo_level >= 3 && comm_rank == 0)
            {
                std::cout << "Mapper: SubModelPart used for " << InterfaceSide << "-ModelPart" << std::endl;
            }

            return rModelPart.GetSubModelPart(name_interface_submodel_part);
        }
        else
        {
            if (echo_level >= 3 && comm_rank == 0)
            {
                std::cout << "Mapper: Main ModelPart used for " << InterfaceSide << "-ModelPart" << std::endl;
            }

            return rModelPart;
        }
    }

    std::unordered_map<std::string, Mapper::Pointer>& MapperFactory::GetRegisteredMappersList()
    {
        static std::unordered_map<std::string, Mapper::Pointer> registered_mappers;

        return registered_mappers;
    }
  
}  // namespace Kratos.

