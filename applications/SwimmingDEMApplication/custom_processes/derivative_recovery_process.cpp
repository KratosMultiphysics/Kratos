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
#include "../swimming_DEM_application.h"
#include "derivative_recovery_process.h"

namespace Kratos
{

DerivativeRecoveryProcess::DerivativeRecoveryProcess(
    ModelPart& rModelPart,
    Parameters Param)
    : Process(), mrModelPart(rModelPart), mMaterialDerivativeContainer(MATERIAL_ACCELERATION)
{
    this->CheckDefaultsAndProcessSettings(Param);
    mStoreFullGradient = Param['store_full_gradient_option'].GetBool();
}

DerivativeRecoveryProcess::DerivativeRecoveryProcess(
    Model &rModel,
    Parameters Param)
    : Process(), mrModelPart(rModel.GetModelPart(Param["model_part_name"].GetString())), mMaterialDerivativeContainer(MATERIAL_ACCELERATION)
{
    this->CheckDefaultsAndProcessSettings(Param);
    mStoreFullGradient = Param['store_full_gradient_option'].GetBool();
}

void DerivativeRecoveryProcess::CheckDefaultsAndProcessSettings(Parameters Param)
{
    Parameters default_parameters( R"(
    {
        "model_part_name"                        : "",
        "distance_factor"                        : 2.0,
        "distance_threshold"                     : 0.001,
        "continuous_distance"                    : true,
        "check_at_each_time_step"                : true,
        "avoid_almost_empty_elements"            : true,
        "deactivate_full_negative_elements"      : true,
        "recover_original_distance_at_each_step" : false
    }  )" );

    Param.ValidateAndAssignDefaults(default_parameters);
}

void DerivativeRecoveryProcess::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();
}

void DerivativeRecoveryProcess::ExecuteInitialize() {

    KRATOS_TRY;


    KRATOS_CATCH("");
}

void DerivativeRecoveryProcess::ExecuteBeforeSolutionLoop() {
    this->ExecuteInitializeSolutionStep();
    this->ExecuteFinalizeSolutionStep();
}

void DerivativeRecoveryProcess::ExecuteInitializeSolutionStep() {

}

void DerivativeRecoveryProcess::ExecuteFinalizeSolutionStep() {
    this->CalculateVectorMaterialDerivative();
}

/* Protected functions ****************************************************/

/* Private functions ****************************************************/

};  // namespace Kratos.
