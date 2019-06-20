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
//                   Jordi Cotela
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "containers/generic_variables_list.h"
#include "includes/kratos_components.h"

namespace Kratos
{
    void GenericVariablesList::AddVariable(const std::string& rVarName)
    {
        //ugly but needed: here we need to get the component type by checking the Kratos components
        if(KratosComponents<Variable<bool>>::Has(rVarName))
            mVariables.push_back(KratosComponents<Variable<bool>>::Get(rVarName));
        else if(KratosComponents<Variable<int>>::Has(rVarName))
            mVariables.push_back(KratosComponents<Variable<int>>::Get(rVarName));
        else if(KratosComponents<Variable<unsigned int>>::Has(rVarName))
            mVariables.push_back(KratosComponents<Variable<unsigned int>>::Get(rVarName));
        else if(KratosComponents<Variable<double>>::Has(rVarName))
            mVariables.push_back(KratosComponents<Variable<double>>::Get(rVarName));
        else if(KratosComponents<Variable<array_1d<double,3>>>::Has(rVarName))
            mVariables.push_back(KratosComponents<Variable<array_1d<double,3>>>::Get(rVarName));
        else if(KratosComponents<Variable<Vector>>::Has(rVarName))
            mVariables.push_back(KratosComponents<Variable<Vector>>::Get(rVarName));
        else if(KratosComponents<Variable<Matrix>>::Has(rVarName))
            mVariables.push_back(KratosComponents<Variable<Matrix>>::Get(rVarName));
        else if(KratosComponents<variable_component_type>::Has(rVarName))
            mVariables.push_back(KratosComponents<variable_component_type>::Get(rVarName));
        else
        {
            KRATOS_ERROR << "variable " << rVarName << " not found as recognized variable type" << std::endl;
        }
        

    }


}  // namespace Kratos.


