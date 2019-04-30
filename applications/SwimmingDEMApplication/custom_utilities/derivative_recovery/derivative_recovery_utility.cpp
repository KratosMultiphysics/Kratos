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
#include "json/json_fwd.hpp"

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
    Parameters rParameters)
    : mrModelPart(rModelPart)
{
    this->CheckDefaultVariablesAreInSettings(rParameters);
    this->CheckDefaultsAndSettings(rParameters["settings"]);
    this->ReadAllVariablePairs(rParameters["variables_for_recovery"]);
    mStoreFullGradient = rParameters["store_full_gradient_option"].GetBool();
}

DerivativeRecoveryUtility::DerivativeRecoveryUtility(
    Model& rModel,
    Parameters rParameters)
    : mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    this->CheckDefaultVariablesAreInSettings(rParameters);
    this->CheckDefaultsAndSettings(rParameters["settings"]);
    this->ReadAllVariablePairs(rParameters["variables_for_recovery"]);
    mStoreFullGradient = rParameters["store_full_gradient_option"].GetBool();
}

void DerivativeRecoveryUtility::ReadAllVariablePairs(Parameters rVariablesForRecovery)
{
    this->ReadVariablePairs("gradient", rVariablesForRecovery);
    this->ReadVariablePairs("material_derivative", rVariablesForRecovery);
}

void DerivativeRecoveryUtility::ReadVariablePairs(std::string OperatorName, Parameters rVariablesForRecovery)
{
    const Parameters derivative_pairs = rVariablesForRecovery[OperatorName];
    for (auto pair_it = derivative_pairs.begin();
        pair_it != derivative_pairs.end(); ++pair_it){
        auto key = pair_it.name();
        mVariablesContainer.AddRecoveryPair(OperatorName, key, derivative_pairs[key].GetString());
    }
}

void DerivativeRecoveryUtility::CheckDefaultVariablesAreInSettings(Parameters rParameters)
{
    Parameters default_parameters(
    R"({
        "model_part_name" : "FluidModelPart",
        "recoverer_name" : "DerivativeRecoveryUtility",
        "store_full_gradient_option" : true,
        "settings" : {},
        "variables_for_recovery" : {}
    })"
    );

    rParameters.ValidateAndAssignDefaults(default_parameters);
    Parameters default_variables_to_recover(

    R"({
        "gradient" : {},
        "divergence" : {},
        "laplacian" : {},
        "material_derivative" : {},
        "rotational" : {}
    })"
    );

    rParameters["variables_for_recovery"].ValidateAndAssignDefaults(default_variables_to_recover);
}

void DerivativeRecoveryUtility::Recover()
{
    KRATOS_TRY;

    for (auto pair : mVariablesContainer.GetVariables("gradient")){
        this->CalculateGradient(pair.first, pair.second);
    }

    for (auto pair : mVariablesContainer.GetVariables("divergence")){
        this->CalculateDivergence(pair.first, pair.second);
    }

    for (auto pair : mVariablesContainer.GetVariables("rotational")){
        this->CalculateRotational(pair.first, pair.second);
    }

    for (auto pair : mVariablesContainer.GetVariables("material_derivative")){
        this->CalculateMaterialDerivative(pair.first, pair.second);
    }

    for (auto pair : mVariablesContainer.GetVariables("laplacian")){
        this->CalculateLaplacian(pair.first, pair.second);
    }

    KRATOS_CATCH("");
}

void DerivativeRecoveryUtility::CalculateGradient(const std::string VariableName, const std::string DerivativeVariableName)
{
    bool was_not_able_to_recover = true;
    was_not_able_to_recover = !this->CalculateGradientIfPossible<DoubleVarType, ArrayVarType>(VariableName, DerivativeVariableName);
    was_not_able_to_recover = !this->CalculateGradientIfPossible<ComponentVarType, ArrayVarType>( VariableName, DerivativeVariableName);

    if (was_not_able_to_recover) {
        KRATOS_THROW_ERROR(std::invalid_argument, VariableName, "is not registered as any type of compatible variable: DOUBLE or SCALAR COMPONENT");
    }
}

void DerivativeRecoveryUtility::CalculateDivergence(const std::string VariableName, const std::string DerivativeVariableName)
{
    bool was_not_able_to_recover = true;
    was_not_able_to_recover = !this->CalculateDivergenceIfPossible<ArrayVarType, DoubleVarType>(VariableName, DerivativeVariableName);

    if (was_not_able_to_recover) {
        KRATOS_THROW_ERROR(std::invalid_argument, VariableName, "is not registered as any type of compatible variable: DOUBLE or SCALAR COMPONENT");
    }
}

void DerivativeRecoveryUtility::CalculateRotational(const std::string VariableName, const std::string DerivativeVariableName)
{
    bool was_not_able_to_recover = true;
    was_not_able_to_recover = !this->CalculateRotationalIfPossible<ArrayVarType, ArrayVarType>(VariableName, DerivativeVariableName);

    if (was_not_able_to_recover) {
        KRATOS_THROW_ERROR(std::invalid_argument, VariableName, "is not registered as any type of compatible variable: ARRAY_1D");
    }
}

void DerivativeRecoveryUtility::CalculateMaterialDerivative(const std::string VariableName, const std::string DerivativeVariableName)
{
    bool was_not_able_to_recover = true;
    was_not_able_to_recover = !this->CalculateMaterialDerivativeIfPossible<DoubleVarType, DoubleVarType>(VariableName, DerivativeVariableName);
    was_not_able_to_recover = !this->CalculateMaterialDerivativeIfPossible<ComponentVarType, DoubleVarType>(VariableName, DerivativeVariableName);
    was_not_able_to_recover = !this->CalculateMaterialDerivativeIfPossible<ArrayVarType, ArrayVarType>(VariableName, DerivativeVariableName);

    if (was_not_able_to_recover) {
        KRATOS_THROW_ERROR(std::invalid_argument, VariableName, "is not registered as any type of compatible variable: DOUBLE, SCALAR COMPONENT or ARRAY_1D");
    }
}

void DerivativeRecoveryUtility::CalculateLaplacian(const std::string VariableName, const std::string DerivativeVariableName)
{
    bool was_not_able_to_recover = true;
    was_not_able_to_recover = !this->CalculateLaplacianIfPossible<DoubleVarType, DoubleVarType>(VariableName, DerivativeVariableName);
    was_not_able_to_recover = !this->CalculateLaplacianIfPossible<ComponentVarType, DoubleVarType>(VariableName, DerivativeVariableName);
    was_not_able_to_recover = !this->CalculateLaplacianIfPossible<ArrayVarType, ArrayVarType>(VariableName, DerivativeVariableName);

    if (was_not_able_to_recover) {
        KRATOS_THROW_ERROR(std::invalid_argument, VariableName, "is not registered as any type of compatible variable: DOUBLE, SCALAR COMPONENT or ARRAY_1D");
    }
}

// template <>
// void DerivativeRecoveryUtility::AddPartialTimeDerivative(const DoubleVarType& rVariable, const DoubleVarType& rTimeDerivativeVariable)
// {
//     this->AddPartialTimeDerivative(rVariable, rTimeDerivativeVariable);
// }

// template <>
// void CalculateMaterialDerivative(const DoubleVarType& rVariable, const DoubleVarType& rTimeDerivativeVariable)
// {
//     this->CalculateMaterialDerivative(rVariable, rTimeDerivativeVariable);
// }

// template <>
// void CalculateMaterialDerivative(const ArrayVarType& rVariable, const ArrayVarType& rTimeDerivativeVariable)
// {
//     this->CalculateMaterialDerivative(rVariable, rTimeDerivativeVariable);
// }

}  // namespace Kratos.
