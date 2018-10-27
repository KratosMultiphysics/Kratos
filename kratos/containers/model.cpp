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
    void Model::Reset()
    {
        mRootModelPartMap.clear();
        //mListOfVariablesLists.clear(); //this has to be done AFTER clearing the RootModelParts
    }
    
    ModelPart& Model::CreateModelPart( const std::string ModelPartName, ModelPart::IndexType NewBufferSize ) 
    {
        KRATOS_TRY
                
        auto search = mRootModelPartMap.find(ModelPartName);
        if( search == mRootModelPartMap.end()) {
            auto pvar_list = Kratos::make_unique<VariablesList>();

             ModelPart* pmodel_part = new ModelPart(ModelPartName, NewBufferSize, pvar_list.get(), *this );
            mRootModelPartMap[ModelPartName] = std::unique_ptr<ModelPart>(pmodel_part); //note that i create it separately since Model is friend of ModelPart but unique_ptr is not

            GetListOfVariableLists().insert(std::move(pvar_list));
            return *(mRootModelPartMap[ModelPartName].get());
        } else {
            KRATOS_ERROR << "trying to create a root modelpart with name " << ModelPartName << " however a ModelPart with the same name already exists";
        }

        KRATOS_CATCH("")
    }

    void Model::DeleteModelPart( const std::string ModelPartName  ) 
    {
        KRATOS_TRY

        if(this->HasModelPart(ModelPartName))
            mRootModelPartMap.erase(ModelPartName); //NOTE: the corresponding variable list should NOT be removed
        else
            KRATOS_WARNING("Info") << "attempting to delete inexisting modelpart : " << ModelPartName << std::endl;

        
        KRATOS_CATCH("")
    }
    
    void Model::RenameModelPart( const std::string OldName, const std::string NewName )
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(this->HasModelPart(OldName)) << "The Old Name is not in model (as a root model part). Required old name was : " << OldName << std::endl;

        KRATOS_ERROR_IF(this->HasModelPart(NewName)) << "The New Name is already existing in model. Proposed name was : " << NewName << std::endl;
        
        mRootModelPartMap[OldName]->Name() = NewName; //change the name of the existing modelpart
        
        CreateModelPart(NewName);
        
        mRootModelPartMap[NewName].swap(mRootModelPartMap[OldName]);
        
        mRootModelPartMap.erase(OldName);

        KRATOS_CATCH("")
    }
    
    
    ModelPart& Model::GetModelPart(const std::string& rFullModelPartName)
    {
        KRATOS_TRY
        
        KRATOS_ERROR_IF( rFullModelPartName.empty() ) << "Attempting to find a "
            << "ModelPart with empty name (\"\")!" << std::endl;
            
        std::vector< std::string > subparts_list = GetSubPartsList(rFullModelPartName);
        
        
        if(subparts_list.size() == 1) //it is a root model part
        {
            auto search = mRootModelPartMap.find(subparts_list[0]);
            if(search != mRootModelPartMap.end())
            {
                return *(search->second);
            }
            else //let's also search it as a flat name - a feature that SHOULD BE DEPRECATED
            {
                for(auto it = mRootModelPartMap.begin(); it!=mRootModelPartMap.end(); it++)
                {
                     ModelPart* pmodel_part = RecursiveSearchByName(subparts_list[0], (it->second.get()));
                     if(pmodel_part != nullptr) //give back the first one that was found
                         return *pmodel_part;
                }
                
                //if we are here we did not find it
                KRATOS_ERROR << "The ModelPart named : \"" << subparts_list[0]
                     << "\" was not found either as root-ModelPart or as a flat name. The total input string was \""
                     << rFullModelPartName << "\"" << std::endl;
            }
        }
        else //it is a submodelpart with the full name provided
        {
            auto search = mRootModelPartMap.find(subparts_list[0]);
            if(search != mRootModelPartMap.end())
            {
                ModelPart* p_model_part = (search->second).get();
                for(unsigned int i=1; i<subparts_list.size(); ++i)
                {
                    KRATOS_ERROR_IF_NOT(p_model_part->HasSubModelPart(subparts_list[i]))
                        << "The ModelPart named : \"" << subparts_list[i]
                        << "\" was not found as SubModelPart of : \""
                        << subparts_list[i-1] << "\". The total input string was \""
                        << rFullModelPartName << "\"" << std::endl;
                    p_model_part = &p_model_part->GetSubModelPart(subparts_list[i]);
                }
                return *p_model_part;
            }
            else 
            {
                KRATOS_ERROR << "root model part " << rFullModelPartName << " not found" << std::endl;
            }
            
        }

        KRATOS_CATCH("")
    }

    bool Model::HasModelPart(const std::string& rFullModelPartName)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF( rFullModelPartName.empty() ) << "Attempting to find a "
            << "ModelPart with empty name (\"\")!" << std::endl;

        std::vector< std::string > subparts_list =  GetSubPartsList(rFullModelPartName);

        //token 0 is the root
        auto search = mRootModelPartMap.find(subparts_list[0]);
        if(search != mRootModelPartMap.end())
        {
            ModelPart* mp = (search->second).get();
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
             ss<< *((it->second).get()) << std::endl << std::endl;
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

    std::vector<std::string> Model::GetSubPartsList(const std::string& rFullModelPartName)
    {
        std::vector<std::string> rSubPartsList;
        std::istringstream iss(rFullModelPartName);
        std::string token;
        rSubPartsList.clear();
        while (std::getline(iss, token, '.'))
        {
            rSubPartsList.push_back(token);
        }
        return rSubPartsList;
    }
    
    ModelPart* Model::RecursiveSearchByName(const std::string& ModelPartName, ModelPart* pModelPart)
    {

        for(auto& part : pModelPart->SubModelParts())
        {
            if(part.Name() == ModelPartName)
                return &part;
            else
            {
                ModelPart* pmodel_part = RecursiveSearchByName(ModelPartName, &part);
                if(pmodel_part != nullptr)
                    return pmodel_part;
            }
        }
        return nullptr;
    }

    void Model::save(Serializer& rSerializer) const
    {
        //we construct auxiliary arrays to avoid having to serialize sets and maps of unique_ptrs
        std::vector<VariablesList* > aux_var_lists;
        std::vector<std::string> aux_names;
        //std::vector<ModelPart* > aux_model_part_pointers;
        aux_var_lists.reserve(GetListOfVariableLists().size());
        aux_names.reserve(mRootModelPartMap.size());
        //aux_model_part_pointers.reserve(mRootModelPartMap.size());

        for(auto it = mRootModelPartMap.begin(); it!=mRootModelPartMap.end(); ++it)
        {
            aux_names.push_back(it->first);
            //aux_model_part_pointers.push_back((it->second).get());
        }

        for(auto it = GetListOfVariableLists().begin(); it!=GetListOfVariableLists().end(); ++it)
            aux_var_lists.push_back(it->get());

        rSerializer.save("ListOfVariablesLists", aux_var_lists);
        rSerializer.save("ModelPartNames", aux_names);

        for(auto it = mRootModelPartMap.begin(); it!=mRootModelPartMap.end(); ++it)
        {
            rSerializer.save(it->first, (it->second).get());
        }
        //rSerializer.save("ModelPartPointers", aux_model_part_pointers);
    }

    void Model::load(Serializer& rSerializer)
    {
        //we construct auxiliary arrays to avoid having to serialize sets and maps of unique_ptrs
        std::vector<VariablesList* > aux_var_lists;
        std::vector<std::string> aux_names;
        // std::vector<ModelPart* > aux_model_part_pointers;

        rSerializer.load("ListOfVariablesLists", aux_var_lists);
        rSerializer.load("ModelPartNames", aux_names);
        // rSerializer.load("ModelPartPointers", aux_model_part_pointers);

        for(IndexType i=0; i<aux_var_lists.size(); ++i) {
            auto p_aux_list = std::unique_ptr<VariablesList>(aux_var_lists[i]);
            GetListOfVariableLists().insert(std::move(p_aux_list)); //NOTE: the ordering may be changed since the pointers are changed, however it should not matter
        }

        for(IndexType i=0; i<aux_names.size(); ++i) {
            //NOTE: CreateModelPart CANNOT be used here
            auto dummy_list = Kratos::make_unique<VariablesList>();
            ModelPart* pmodel_part = new ModelPart(aux_names[i], 1, dummy_list.get(), *this );
            rSerializer.load("MP", pmodel_part);
            mRootModelPartMap.insert(std::make_pair(aux_names[i],std::unique_ptr<ModelPart>(pmodel_part)));
        }


    }


}  // namespace Kratos.


