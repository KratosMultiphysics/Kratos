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

const std::set<std::string> RecoveryVariablesContainer::smOperators {
    "gradient",
    "divergence",
    "rotational",
    "material_derivative",
    "laplacian"};

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
    this->ReadVariablePairs("divergence", rVariablesForRecovery);
    this->ReadVariablePairs("rotational", rVariablesForRecovery);
    this->ReadVariablePairs("material_derivative", rVariablesForRecovery);
    this->ReadVariablePairs("laplacian", rVariablesForRecovery);
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

    for (const auto& var : mVariablesContainer.GetVariables()){
        this->SetCurrentVariableRecovery(var.first);

        // This grants access to info. about all the derivatives during the calculation
        // of any particular one.

        this->InitializeRecovery();

        this->CalculateGradient();
        this->CalculateDivergence();
        this->CalculateRotational();
        this->CalculateMaterialDerivative();
        this->CalculateLaplacian();

        this->FinalizeRecovery();
    }

    KRATOS_CATCH("");
}

void DerivativeRecoveryUtility::CalculateGradient()
{
    if (this->MustRecover("gradient")){
        const auto& var_name = this->GetCurrentVariableName();
        const auto& deriv_var_name = this->GetDerivativeVariableName("gradient");
        bool was_able_to_recover = this->CalculateGradientIfPossible<DoubleVarType, ArrayVarType>(var_name, deriv_var_name);
        if (! was_able_to_recover){
            was_able_to_recover = this->CalculateGradientIfPossible<ComponentVarType, ArrayVarType>(var_name, deriv_var_name);
        }
        if (! was_able_to_recover){
            was_able_to_recover = this->CalculateGradientIfPossible<ArrayVarType, TensorVarType>(var_name, deriv_var_name);
        }
        if (! was_able_to_recover) {
            KRATOS_THROW_ERROR(std::invalid_argument, var_name, " is not registered as any type of compatible variable: DOUBLE or SCALAR COMPONENT");
        }
    }
}

void DerivativeRecoveryUtility::CalculateDivergence()
{
    if (this->MustRecover("divergence")){
        const auto& var_name = this->GetCurrentVariableName();
        const auto& deriv_var_name = this->GetDerivativeVariableName("divergence");
        bool was_able_to_recover = this->CalculateDivergenceIfPossible<ArrayVarType, DoubleVarType>(var_name, deriv_var_name);

        if (! was_able_to_recover) {
            KRATOS_THROW_ERROR(std::invalid_argument, var_name, " is not registered as any type of compatible variable: DOUBLE or SCALAR COMPONENT");
        }
    }
}

void DerivativeRecoveryUtility::CalculateRotational()
{
    if (this->MustRecover("rotational")){
        const auto& var_name = this->GetCurrentVariableName();
        const auto& deriv_var_name = this->GetDerivativeVariableName("rotational");
        bool was_able_to_recover = this->CalculateRotationalIfPossible<ArrayVarType, ArrayVarType>(var_name, deriv_var_name);

        if (! was_able_to_recover) {
            KRATOS_THROW_ERROR(std::invalid_argument, var_name, " is not registered as any type of compatible variable: DOUBLE or SCALAR COMPONENT");
        }
    }
}

void DerivativeRecoveryUtility::CalculateMaterialDerivative()
{
    if (this->MustRecover("material_derivative")){
        const auto& var_name = this->GetCurrentVariableName();
        const auto& deriv_var_name = this->GetDerivativeVariableName("material_derivative");
        bool was_able_to_recover = this->CalculateMaterialDerivativeIfPossible<DoubleVarType, DoubleVarType>(var_name, deriv_var_name);
        if (! was_able_to_recover){
            was_able_to_recover = this->CalculateMaterialDerivativeIfPossible<ComponentVarType, DoubleVarType>(var_name, deriv_var_name);
        }
        if (! was_able_to_recover){
            was_able_to_recover = this->CalculateMaterialDerivativeIfPossible<ArrayVarType, ArrayVarType>(var_name, deriv_var_name);
        }
        if (! was_able_to_recover) {
            KRATOS_THROW_ERROR(std::invalid_argument, var_name, " is not registered as any type of compatible variable: DOUBLE or SCALAR COMPONENT");
        }
    }
}

void DerivativeRecoveryUtility::CalculateLaplacian()
{
    if (this->MustRecover("laplacian")){
        const auto& var_name = this->GetCurrentVariableName();
        const auto& deriv_var_name = this->GetDerivativeVariableName("laplacian");
        bool was_able_to_recover = this->CalculateLaplacianIfPossible<DoubleVarType, DoubleVarType>(var_name, deriv_var_name);
        if (! was_able_to_recover){
            was_able_to_recover = this->CalculateLaplacianIfPossible<ComponentVarType, DoubleVarType>(var_name, deriv_var_name);
        }
        if (! was_able_to_recover){
            was_able_to_recover = this->CalculateLaplacianIfPossible<ArrayVarType, ArrayVarType>(var_name, deriv_var_name);
        }
        if (! was_able_to_recover) {
            KRATOS_THROW_ERROR(std::invalid_argument, var_name, " is not registered as any type of compatible variable: DOUBLE or SCALAR COMPONENT");
        }
    }
}

}  // namespace Kratos.
