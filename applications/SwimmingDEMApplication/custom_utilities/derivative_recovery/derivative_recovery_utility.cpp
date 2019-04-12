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

    for (auto pair : mrVariablesContainer.GetVariables("gradient")){
        this->CalculateGradient(pair.first, pair.second);
    }

    for (auto pair : mrVariablesContainer.GetVariables("divergence")){
        this->CalculateDivergence(pair.first, pair.second);
    }

    for (auto pair : mrVariablesContainer.GetVariables("rotational")){
        this->CalculateRotational(pair.first, pair.second);
    }

    for (auto pair : mrVariablesContainer.GetVariables("material derivative")){
        this->CalculateMaterialDerivative(pair.first, pair.second);
    }

    for (auto pair : mrVariablesContainer.GetVariables("Laplacian")){
        this->CalculateLaplacian(pair.first, pair.second);
    }

    KRATOS_CATCH("");
}

}  // namespace Kratos.
