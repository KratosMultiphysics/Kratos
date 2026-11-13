//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   
//

// System includes
#include <algorithm>

// External includes 

// Project includes
#include "includes/kratos_components.h"
#include "containers/variables_list.h"


namespace Kratos
{
    void VariablesList::save(Serializer& rSerializer) const
    {
        std::size_t size = mVariables.size();
        rSerializer.save("Size", size);
        for (std::size_t i = 0; i < size; i++)
        {
            rSerializer.save("VariableName", mVariables[i]->Name());
        }

        std::size_t dof_size = mDofVariables.size();
        rSerializer.save("DofSize", dof_size);
        for (std::size_t i = 0; i < dof_size; i++)
        {
            rSerializer.save("DofVariableName", mDofVariables[i]->Name());
            if(mDofReactions[i] == nullptr){
                rSerializer.save("HasReaction", false);
            }
            else{
                rSerializer.save("HasReaction", true);
                rSerializer.save("DofReactionName", mDofReactions[i]->Name());
            }
        }
    }

    void VariablesList::load(Serializer& rSerializer)
    {
        // The archive is read into a temporary list, so that a rejected archive leaves this list unchanged
        VariablesList loaded_list;
        std::size_t size;
        rSerializer.load("Size", size);
        std::string name;
        for (std::size_t i = 0; i < size; i++)
        {
            rSerializer.load("VariableName", name);
            loaded_list.Add(*KratosComponents<VariableData>::pGet(name));
        }
        rSerializer.load("DofSize", size);
        for (std::size_t i = 0; i < size; i++)
        {
            rSerializer.load("DofVariableName", name);
            bool has_reaction;
            rSerializer.load("HasReaction", has_reaction);
            
            if(has_reaction){
                std::string reaction_name;
                rSerializer.load("DofReactionName", reaction_name);
                loaded_list.AddDof(KratosComponents<VariableData>::pGet(name),
                                   KratosComponents<VariableData>::pGet(reaction_name));

            }
            else{
                loaded_list.AddDof(KratosComponents<VariableData>::pGet(name), nullptr);
            }

        }

        // An unreferenced list may be extended, so the archive must start with its variables. A shared list holds data
        // blocks with this layout and must match exactly (also an empty one referenced only by sub model parts).
        // Compared by key, as an inline constexpr variable has one instance per shared library.
        const auto is_same = [](const VariableData* pA, const VariableData* pB) { return *pA == *pB; };
        const auto& r_loaded_dofs = loaded_list.mDofVariables;
        const auto& r_loaded_reactions = loaded_list.mDofReactions;
        const bool is_unreferenced = use_count() <= 1;
        const std::size_t current_size = mVariables.size();
        const bool is_matching_variables = is_unreferenced
            ? current_size <= loaded_list.mVariables.size()
              && std::equal(mVariables.begin(), mVariables.end(), loaded_list.mVariables.begin(), is_same)
            : std::equal(mVariables.begin(), mVariables.end(), loaded_list.mVariables.begin(),
                         loaded_list.mVariables.end(), is_same);
        // Each Dof stores its index in this list, so the common dofs must keep their order. Extra loaded dofs are
        // appended (restart adds dofs after loading); extra current dofs (added later) keep the loaded indices valid.
        const std::size_t common_dofs = std::min(mDofVariables.size(), r_loaded_dofs.size());
        const bool is_matching = is_matching_variables
            && std::equal(mDofVariables.begin(), mDofVariables.begin() + common_dofs, r_loaded_dofs.begin(), is_same);
        KRATOS_ERROR_IF_NOT(is_matching)
            << "The archive does not match this VariablesList, which may be shared by data containers: it must hold "
            << "the same variables and start with the same dofs, in the same order [ current variables = "
            << this->size() << ", loaded variables = " << loaded_list.size() << ", current dofs = "
            << mDofVariables.size() << ", loaded dofs = " << r_loaded_dofs.size() << " ].\n";
        for (std::size_t i = 0; i < common_dofs; ++i) {
            KRATOS_ERROR_IF(mDofReactions[i] != nullptr && r_loaded_reactions[i] != nullptr
                            && !is_same(mDofReactions[i], r_loaded_reactions[i]))
                << "The archive gives the dof " << mDofVariables[i]->Name() << " the reaction "
                << r_loaded_reactions[i]->Name() << ", but this VariablesList has "
                << mDofReactions[i]->Name() << ".\n";
        }

        for (std::size_t i = current_size; i < loaded_list.mVariables.size(); ++i) {
            Add(*loaded_list.mVariables[i]);
        }
        for (std::size_t i = 0; i < r_loaded_dofs.size(); ++i) {
            if (i >= mDofVariables.size()) {
                AddDof(r_loaded_dofs[i], r_loaded_reactions[i]);
            } else if (mDofReactions[i] == nullptr) {
                mDofReactions[i] = r_loaded_reactions[i];
            }
        }
    }

    
}  // namespace Kratos.


