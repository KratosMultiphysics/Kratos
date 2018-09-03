//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Pooyan Dadvand
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include <iostream>
#include <string>
#include <sstream>


namespace Kratos
{
    void Model::AddModelPart( ModelPart* pModelPart) //TODO: DEPRECATED. to be removed. this is TEMPORARY
    {
        KRATOS_TRY

        //TODO: flat map should disappear in the future!!
        auto search = mflat_map.find(pModelPart->Name());
        if( search == mflat_map.end())
        {
            mflat_map[pModelPart->Name()] = pModelPart;

            //walk the submodelparts
            for(auto& part : pModelPart->SubModelParts())
                AddModelPartRawPointer(&part);
        }
        else
        {
            if(&(*search->second) != &*(pModelPart))
                KRATOS_ERROR << "trying to add to the Model two DISTINCT model parts with the same name. This should be possible (and it will be in the future) if they belong to two different root model_parts, but it is currently disallowed";

        }

        //add the root model part to the list
        ModelPart& root_model_part = pModelPart->GetRootModelPart();
        mRootModelPartMap[root_model_part.Name()] = &root_model_part;

        KRATOS_CATCH("")
    }

    void Model::AddModelPartRawPointer( ModelPart* pModelPart) //TODO: DEPRECATED. to be removed. this is TEMPORARY
    {
        KRATOS_TRY

        //TODO: flat map should disappear in the future!!
        auto search = mflat_map.find(pModelPart->Name());
        if( search == mflat_map.end())
        {
            mflat_map[pModelPart->Name()] = pModelPart;

            //walk the submodelparts
            for(auto& part : pModelPart->SubModelParts())
                AddModelPartRawPointer(&part);
        }
        else
        {
            if(&(*search->second) != &*(pModelPart))
                KRATOS_ERROR << "trying to add to the Model two DISTINCT model parts with the same name. This should be possible (and it will be in the future) if they belong to two different root model_parts, but it is currently disallowed";

        }

        //add the root model part to the list
        ModelPart& root_model_part = pModelPart->GetRootModelPart();
        mRootModelPartMap[root_model_part.Name()] = &root_model_part;

        KRATOS_CATCH("")
    }

    ModelPart& Model::GetModelPart(const std::string& rFullModelPartName)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF( rFullModelPartName.empty() ) << "Attempting to find a "
            << "ModelPart with empty name (\"\")!" << std::endl;

        auto search = mflat_map.find(rFullModelPartName);
        if(search != mflat_map.end()) {
            // TODO enable the warning
            // KRATOS_WARNING_IF("Model", (search->second)->IsSubModelPart())
            //     << "Getting a SubModelPart from the Model without "
            //     << "specifying the RootModelPart is deprecated and will be removed\n"
            //     << "Please use e.g \"RootModelPart.SubModelPart.SubSubModelPart\" "
            //     << "as input for this function" << std::endl;
            return *(search->second);
        }
        else //look for it in the "root_map" which is where it is suppossed to be finally
        {
            std::vector< std::string > subparts_list;
            GetSubPartsList(rFullModelPartName, subparts_list);

            //token 0 is the root
            auto search = mRootModelPartMap.find(subparts_list[0]);
            if(search != mRootModelPartMap.end())
            {
                ModelPart* mp = search->second;
                for(unsigned int i=1; i<subparts_list.size(); i++)
                {
                    KRATOS_ERROR_IF_NOT(mp->HasSubModelPart(subparts_list[i]))
                        << "The ModelPart named : \"" << subparts_list[i]
                        << "\" was not found as SubModelPart of : \""
                        << subparts_list[i-1] << "\". The total input string was \""
                        << rFullModelPartName << "\"" << std::endl;
                    mp = &(mp->GetSubModelPart(subparts_list[i]));
                }
                return *mp;
            }
            else
            {
                KRATOS_ERROR << "The ModelPart named : \"" << subparts_list[0]
                    << "\" was not found as root-ModelPart. The total input string was \""
                    << rFullModelPartName << "\"" << std::endl;
            }

        }

        KRATOS_CATCH("")
    }

    bool Model::HasModelPart(const std::string& rFullModelPartName)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF( rFullModelPartName.empty() ) << "Attempting to find a "
            << "ModelPart with empty name (\"\")!" << std::endl;

        std::vector< std::string > subparts_list;
        GetSubPartsList(rFullModelPartName, subparts_list);

        //token 0 is the root
        auto search = mRootModelPartMap.find(subparts_list[0]);
        if(search != mRootModelPartMap.end())
        {
            ModelPart* mp = search->second;
            for(unsigned int i=1; i<subparts_list.size(); i++)
            {
                if(!mp->HasSubModelPart(subparts_list[i]))
                    return false;
                mp = &(mp->GetSubModelPart(subparts_list[i]));
            }
            return true;
        }
        else
        {
            return false;
        }

        KRATOS_CATCH("")
    }

    std::string Model::Info() const
    {
        std::stringstream ss;
        for(auto it = mRootModelPartMap.begin(); it!=mRootModelPartMap.end(); it++)
        {
            ss<< it->second->Info() << std::endl << std::endl;
        }
        return ss.str();
    }

    /// Print information about this object.

    void Model::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.

    void Model::PrintData(std::ostream& rOStream) const
    {
    }

    void Model::GetSubPartsList(const std::string& rFullModelPartName,
                                std::vector<std::string>& rSubPartsList)
    {
        std::istringstream iss(rFullModelPartName);
        std::string token;
        rSubPartsList.clear();
        while (std::getline(iss, token, '.'))
        {
            rSubPartsList.push_back(token);
        }
    }

}  // namespace Kratos.


