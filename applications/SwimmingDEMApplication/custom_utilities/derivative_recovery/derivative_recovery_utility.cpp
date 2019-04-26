//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Guillermo Casas (gcasas@cimne.upc.edu)
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "utilities/openmp_utils.h"
#include "processes/find_nodal_h_process.h"

// Application includes
#include "../../swimming_DEM_application.h"
#include "derivative_recovery_utility.h"

namespace Kratos
{

DerivativeRecoveryUtility::DerivativeRecoveryUtility(
    ModelPart& rModelPart,
    Parameters rParameters,
    RecoveryVariablesContainer& rVariablesContainer)
    : mStoreFullGradient(false), mrModelPart(rModelPart), mrVariablesContainer(rVariablesContainer)
{
    this->CheckDefaultsAndSettings(rParameters);
}

DerivativeRecoveryUtility::DerivativeRecoveryUtility(
    Model &rModel,
    Parameters rParameters,
    RecoveryVariablesContainer& rVariablesContainer)
    : mStoreFullGradient(false), mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString())), mrVariablesContainer(rVariablesContainer)
{
    this->CheckDefaultsAndSettings(rParameters);
}

void DerivativeRecoveryUtility::Recover()
{

    KRATOS_TRY;
    // for (auto pair : mrVariablesContainer.GetVariables("gradient")){
    //     this->CalculateGradient(pair.first, pair.second);
    // }

    // for (auto pair : mrVariablesContainer.GetVariables("divergence")){
    //     this->CalculateDivergence(pair.first, pair.second);
    // }

    // for (auto pair : mrVariablesContainer.GetVariables("rotational")){
    //     this->CalculateRotational(pair.first, pair.second);
    // }

    for (auto variable_name_pair : mrVariablesContainer.GetVariables("material_derivative")){
        auto variable = variable_name_pair.first;
        auto derivative_variable = variable_name_pair.second;
        this->CalculateMaterialDerivative(variable, derivative_variable);
    }

    // for (auto pair : mrVariablesContainer.GetVariables("laplacian")){
    //     this->CalculateLaplacian(pair.first, pair.second);
    // }

    KRATOS_CATCH("");
}

void DerivativeRecoveryUtility::CalculateMaterialDerivative(const std::string VariableName, const std::string DerivativeVariableName)
{
    if (KratosComponents<DoubleVarType>::Has(VariableName)){
        const DoubleVarType& variable = KratosComponents<DoubleVarType>::Get(VariableName);
        if (KratosComponents<DoubleVarType>::Has(DerivativeVariableName)){
            const DoubleVarType& derivative_variable = KratosComponents<DoubleVarType>::Get(DerivativeVariableName);
            this->CalculateMaterialDerivative(variable, derivative_variable);
        }

        else {
            KRATOS_THROW_ERROR(std::invalid_argument, "DerivativeRecoveryUtility: The derivative variable is not registered in Kratos with a type matching that of the variable to derivate: ", DerivativeVariableName);
        }
    }

    else if (KratosComponents<ArrayVarType>::Has(VariableName)){
        const ArrayVarType& variable = KratosComponents<ArrayVarType>::Get(VariableName);
        if (KratosComponents<ArrayVarType>::Has(DerivativeVariableName)){
            const ArrayVarType& derivative_variable = KratosComponents<ArrayVarType>::Get(DerivativeVariableName);
            this->CalculateMaterialDerivative(variable, derivative_variable);
        }

        else {
            KRATOS_THROW_ERROR(std::invalid_argument, "DerivativeRecoveryUtility: The derivative variable is not registered in Kratos with a type matching that of the variable to derivate: ", DerivativeVariableName);
        }
    }

    else {
        KRATOS_THROW_ERROR(std::invalid_argument, VariableName, "is not registered as any type of compatible variable: DOUBLE or ARRAY_1D");
    }
}

// template <>
// void DerivativeRecoveryUtility::AddPartialTimeDerivative(const ScalarVariableType& rVariable, const ScalarVariableType& rTimeDerivativeVariable)
// {
//     this->AddPartialTimeDerivative(rVariable, rTimeDerivativeVariable);
// }

// template <>
// void CalculateMaterialDerivative(const ScalarVariableType& rVariable, const ScalarVariableType& rTimeDerivativeVariable)
// {
//     this->CalculateMaterialDerivative(rVariable, rTimeDerivativeVariable);
// }

// template <>
// void CalculateMaterialDerivative(const VectorVariableType& rVariable, const VectorVariableType& rTimeDerivativeVariable)
// {
//     this->CalculateMaterialDerivative(rVariable, rTimeDerivativeVariable);
// }

}  // namespace Kratos.
