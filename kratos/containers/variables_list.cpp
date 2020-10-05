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

// External includes 

// Project includes
#include "includes/kratos_components.h"
#include "containers/variables_list.h"


namespace Kratos
{
    void VariablesList::ApplyVisitor(VariablesListVisitorBase& TheVisitor){
        for (auto p_variable : mVariables)
        {
            p_variable->AcceptVisitor(TheVisitor);
        }
    }

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
        std::size_t size;
        rSerializer.load("Size", size);
        std::string name;
        for (std::size_t i = 0; i < size; i++)
        {
            rSerializer.load("VariableName", name);
            Add(*KratosComponents<VariableData>::pGet(name));
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
                AddDof(KratosComponents<VariableData>::pGet(name), KratosComponents<VariableData>::pGet(reaction_name));

            }
            else{
                AddDof(KratosComponents<VariableData>::pGet(name), nullptr);
            }

        }

    }
    void VariablesList::AddVariable(const std::string& rVarName)
    {
        //ugly but needed: here we need to get the component type by checking the Kratos components
        if(KratosComponents<Variable<bool>>::Has(rVarName))
            Add(KratosComponents<Variable<bool>>::Get(rVarName));
        else if(KratosComponents<Variable<int>>::Has(rVarName))
            Add(KratosComponents<Variable<int>>::Get(rVarName));
        else if(KratosComponents<Variable<unsigned int>>::Has(rVarName))
            Add(KratosComponents<Variable<unsigned int>>::Get(rVarName));
        else if(KratosComponents<Variable<double>>::Has(rVarName))
            Add(KratosComponents<Variable<double>>::Get(rVarName));
        else if(KratosComponents<Variable<array_1d<double,3>>>::Has(rVarName))
            Add(KratosComponents<Variable<array_1d<double,3>>>::Get(rVarName));
        else if(KratosComponents<Variable<Vector>>::Has(rVarName))
            Add(KratosComponents<Variable<Vector>>::Get(rVarName));
        else if(KratosComponents<Variable<Matrix>>::Has(rVarName))
            Add(KratosComponents<Variable<Matrix>>::Get(rVarName));
        else
        {
            KRATOS_ERROR << "variable " << rVarName << " not found as recognized variable type" << std::endl;
        }
    }
    
}  // namespace Kratos.


