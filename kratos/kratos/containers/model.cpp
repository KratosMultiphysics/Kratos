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
    void Model::AddModelPart( ModelPart::Pointer pmodel_part) //TODO: DEPRECATED. to be removed. this is TEMPORARY
    {
        KRATOS_TRY
        
        //TODO: flat map should disappear in the future!!
        auto search = mflat_map.find(pmodel_part->Name());
        if( search == mflat_map.end())
        {
            mflat_map[pmodel_part->Name()] = pmodel_part.get();
        }
        else
        {
            if(&(*search->second) != &*(pmodel_part.get()))
                KRATOS_ERROR << "trying to add to the Model two DISTINCT model parts with the same name. This should be possible (and it will be in the future) if they belong to two different root model_parts, but it is currently disallowed";
                
        }
        
        //add the root model part to the list
        ModelPart& root_model_part = pmodel_part->GetRootModelPart();
        mroot_map[root_model_part.Name()] = &root_model_part;
        
        KRATOS_CATCH("")
    }
    
    ModelPart& Model::GetModelPart(std::string name)
    {
        KRATOS_TRY
        
        auto search = mflat_map.find(name);
        if(search != mflat_map.end()) {
            //TODO: issue a warning message telling that this behaviour is DEPRECATED and that it will be removed
            return *(search->second);
        }
        else //look for it in the "root_map" which is where it is suppossed to be finally
        {
            std::vector< std::string > subparts_list;
            std::istringstream iss(name);
            std::string token;
            while (std::getline(iss, token, '.'))
            {
                subparts_list.push_back(token);
            }
            
            //token 0 is the root
            auto search = mroot_map.find(subparts_list[0]);
            if(search != mroot_map.end()) 
            {
                ModelPart* mp = search->second;
                for(unsigned int i=1; i<subparts_list.size(); i++)
                {
                    if(! mp->HasSubModelPart(subparts_list[i]))
                        KRATOS_ERROR << "The model part named : \"" << subparts_list[i] << "\" was not found as submodelpart of : \"" << subparts_list[i-1] << "\" . total input string was \"" << name << "\"" ;
                    mp = &(mp->GetSubModelPart(subparts_list[i]));
                }
                return *mp;
            }
            else
                KRATOS_ERROR << "The model part named : \"" << subparts_list[0] << "\" was not found as root model part.  total input string was \"" << name << "\"" ;
            
        }
            
        KRATOS_CATCH("")
    }

  std::string Model::Info() const
    {
        std::stringstream ss;
        for(auto it = mroot_map.begin(); it!=mroot_map.end(); it++)
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

  
}  // namespace Kratos.


