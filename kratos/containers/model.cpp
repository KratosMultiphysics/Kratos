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
#include "containers/model.h"

namespace Kratos {
namespace { // empty namespace for internal functions
/**
 * @brief This method gets the names of all parent-modelparts given a submodelpart-name
 * @param rModelPart The SubModelPart for which the parents-modelpart-names are to be extracted
 * @param rModelPartNames The names of the ModelParts
 * @TODO remove this function when the flat-map is removed (it will no longer be needed)
 */
void GetNameWithAscendants(const ModelPart& rModelPart, std::vector<std::string>& rModelPartNames)
{
    rModelPartNames.insert(rModelPartNames.begin(), rModelPart.Name()); // "push_front"
    if (rModelPart.IsSubModelPart()) {
        GetNameWithAscendants(*rModelPart.GetParentModelPart(), rModelPartNames);
    }
}
} // empty namespace for internal functions

void Model::Reset()
{
    mRootModelPartMap.clear();
}

ModelPart& Model::CreateModelPart( const std::string ModelPartName, ModelPart::IndexType NewBufferSize )
{
    KRATOS_TRY

    auto search = mRootModelPartMap.find(ModelPartName);
    if( search == mRootModelPartMap.end()) {
        auto pvar_list = Kratos::make_intrusive<VariablesList>();

            ModelPart* pmodel_part = new ModelPart(ModelPartName, NewBufferSize, pvar_list, *this );
        mRootModelPartMap[ModelPartName] = std::unique_ptr<ModelPart>(pmodel_part); //note that i create it separately since Model is friend of ModelPart but unique_ptr is not

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

    std::vector< std::string > subparts_list = SplitSubModelPartHierarchy(rFullModelPartName);


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
                if (pmodel_part != nullptr) { //give back the first one that was found
                    // Get the names of the parent-modelparts to print them in the warning
                    std::vector<std::string> model_part_names;
                    GetNameWithAscendants(*pmodel_part, model_part_names);

                    std::stringstream msg;
                    msg << model_part_names[0];
                    for (std::size_t i=1; i<model_part_names.size(); ++i) {
                        msg << "." << model_part_names[1];
                    }

                    KRATOS_WARNING("Model") << "DEPRECATION_WARNING: The ModelPart \"" << subparts_list[0] << "\"\nis retrieved from the Model by using the flat-map!\nThis will be removed end of November 2019\nPlease prepend the Parent-ModelPart-Names like this:\n\"" << msg.str() << "\"" << std::endl;

                    return *pmodel_part;
                }
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

const ModelPart& Model::GetModelPart(const std::string& rFullModelPartName) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF( rFullModelPartName.empty() ) << "Attempting to find a "
        << "ModelPart with empty name (\"\")!" << std::endl;

    std::vector< std::string > subparts_list = SplitSubModelPartHierarchy(rFullModelPartName);


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
                if (pmodel_part != nullptr) { //give back the first one that was found
                    // Get the names of the parent-modelparts to print them in the warning
                    std::vector<std::string> model_part_names;
                    GetNameWithAscendants(*pmodel_part, model_part_names);

                    std::stringstream msg;
                    msg << model_part_names[0];
                    for (std::size_t i=1; i<model_part_names.size(); ++i) {
                        msg << "." << model_part_names[1];
                    }

                    KRATOS_WARNING("Model") << "DEPREATION_WARNING: The ModelPart \"" << subparts_list[0] << "\"\nis retrieved from the Model by using the flat-map!\nThis will be removed end of November 2019\nPlease prepend the Parent-ModelPart-Names like this:\n\"" << msg.str() << "\"" << std::endl;

                    return *pmodel_part;
                }
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

bool Model::HasModelPart(const std::string& rFullModelPartName) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF( rFullModelPartName.empty() ) << "Attempting to find a "
        << "ModelPart with empty name (\"\")!" << std::endl;

    std::vector< std::string > subparts_list = SplitSubModelPartHierarchy(rFullModelPartName);

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

std::vector<std::string> Model::GetModelPartNames() const
{
    std::vector<std::string> model_parts_names;

    // We fill the vector
    for (auto& mps : mRootModelPartMap) {
        const std::string& r_root_mp_name = mps.first;
        model_parts_names.push_back(r_root_mp_name);

        // First level of submodelparts
        auto& p_root_mp = mps.second;
        if (p_root_mp->NumberOfSubModelParts() > 0) {
            const std::vector<std::string>& sub_model_part_names = p_root_mp->GetSubModelPartNames();
            for (auto& r_sub_name : sub_model_part_names) {
                model_parts_names.push_back(r_root_mp_name + "." + r_sub_name);
            }

            // Second level of submodelparts
            for (auto& r_sub_mp : p_root_mp->SubModelParts()) {
                if (r_sub_mp.NumberOfSubModelParts() > 0) {
                    const std::string& r_sub_name = r_sub_mp.Name();
                    const std::vector<std::string>& sub_sub_model_part_names = r_sub_mp.GetSubModelPartNames();
                    for (auto& r_sub_sub_name : sub_sub_model_part_names) {
                        model_parts_names.push_back(r_root_mp_name + "." + r_sub_name + "." + r_sub_sub_name);
                    }
                }
            }
        }
    }

    return model_parts_names;
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

std::vector<std::string> Model::SplitSubModelPartHierarchy(const std::string& rFullModelPartName) const
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

ModelPart* Model::RecursiveSearchByName(const std::string& ModelPartName, ModelPart* pModelPart) const
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
    std::vector<std::string> aux_names;
    aux_names.reserve(mRootModelPartMap.size());

    for(auto it = mRootModelPartMap.begin(); it!=mRootModelPartMap.end(); ++it)
    {
        aux_names.push_back(it->first);
    }

    rSerializer.save("ModelPartNames", aux_names);

    for(auto it = mRootModelPartMap.begin(); it!=mRootModelPartMap.end(); ++it)
    {
        rSerializer.save(it->first, (it->second).get());
    }
}

void Model::load(Serializer& rSerializer)
{
    //we construct auxiliary arrays to avoid having to serialize sets and maps of unique_ptrs
    std::vector<std::string> aux_names;

    rSerializer.load("ModelPartNames", aux_names);

    for(IndexType i=0; i<aux_names.size(); ++i) {
        //NOTE: CreateModelPart CANNOT be used here
        auto dummy_list = Kratos::make_intrusive<VariablesList>();
        ModelPart* pmodel_part = new ModelPart(aux_names[i], 1, dummy_list, *this );
        rSerializer.load(aux_names[i], pmodel_part);
        mRootModelPartMap.insert(std::make_pair(aux_names[i],std::unique_ptr<ModelPart>(pmodel_part)));
    }
}

}  // namespace Kratos.
