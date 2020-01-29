//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes

// Project includes
#include "factories/linear_solver_factory.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "processes/variational_distance_calculation_process.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"

// Include base h
#include "rans_wall_distance_calculation_process.h"

namespace Kratos
{
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
RansWallDistanceCalculationProcess<TSparseSpace, TDenseSpace, TLinearSolver>::RansWallDistanceCalculationProcess(
    Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"          : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "max_iterations"           : 10,
            "echo_level"               : 0,
            "wall_flag_variable_name"  : "STRUCTURE",
            "wall_flag_variable_value" : true,
            "linear_solver_settings" : {
                "solver_type"     : "amgcl"
            }
        })");

    mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mpLinearSolver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(
        mrParameters["linear_solver_settings"]);

    mMaxIterations = mrParameters["max_iterations"].GetInt();
    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mWallFlagVariableName = mrParameters["wall_flag_variable_name"].GetString();
    mWallFlagVariableValue = mrParameters["wall_flag_variable_value"].GetBool();

    KRATOS_CATCH("");
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
int RansWallDistanceCalculationProcess<TSparseSpace, TDenseSpace, TLinearSolver>::Check()
{
    KRATOS_TRY


    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, DISTANCE);

    return 0.0;

    KRATOS_CATCH("");
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void RansWallDistanceCalculationProcess<TSparseSpace, TDenseSpace, TLinearSolver>::ExecuteInitialize()
{
    CalculateWallDistances();
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
std::string RansWallDistanceCalculationProcess<TSparseSpace, TDenseSpace, TLinearSolver>::Info() const
{
    return std::string("RansWallDistanceCalculationProcess");
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void RansWallDistanceCalculationProcess<TSparseSpace, TDenseSpace, TLinearSolver>::PrintInfo(
    std::ostream& rOStream) const
{
    rOStream << this->Info();
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void RansWallDistanceCalculationProcess<TSparseSpace, TDenseSpace, TLinearSolver>::PrintData(
    std::ostream& rOStream) const
{
}

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
void RansWallDistanceCalculationProcess<TSparseSpace, TDenseSpace, TLinearSolver>::CalculateWallDistances()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const Flags& r_wall_flag = KratosComponents<Flags>::Get(mWallFlagVariableName);

    VariableUtils variable_utilities;

    variable_utilities.SetScalarVar(DISTANCE, 1.0, r_model_part.Nodes());
    variable_utilities.SetScalarVarForFlag(DISTANCE, 0.0, r_model_part.Nodes(),
                                           r_wall_flag, mWallFlagVariableValue);

    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];

    if (domain_size == 2)
    {
        VariationalDistanceCalculationProcess<2, TSparseSpace, TDenseSpace, TLinearSolver> distance_calculation_process(
            r_model_part, mpLinearSolver, mMaxIterations);
        distance_calculation_process.Execute();
    }
    else if (domain_size == 3)
    {
        VariationalDistanceCalculationProcess<3, TSparseSpace, TDenseSpace, TLinearSolver> distance_calculation_process(
            r_model_part, mpLinearSolver, mMaxIterations);
        distance_calculation_process.Execute();
    }
    else
    {
        KRATOS_ERROR << "Unknown domain size = " << domain_size;
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Wall distances calculated in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

// template instantiations

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;

template class RansWallDistanceCalculationProcess<SparseSpaceType, LocalSpaceType, LinearSolver<SparseSpaceType, LocalSpaceType>>;

} // namespace Kratos.
