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
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_nut_y_plus_wall_function_update_process.h"

namespace Kratos
{
RansNutYPlusWallFunctionUpdateProcess::RansNutYPlusWallFunctionUpdateProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "linear_region_settings":{
                "c1"              : 0.1,
                "beta_zero"       : 0.0708
            },
            "log_region_settings":{
                "von_karman"      : 0.41
            },
            "region_type" : "assumes_both_regions",
            "echo_level"  : 0,
            "min_value"   : 1e-18
        })");

    rParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mEchoLevel = rParameters["echo_level"].GetInt();
    mModelPartName = rParameters["model_part_name"].GetString();
    mVonKarman = rParameters["von_karman"].GetDouble();
    mBetaZero = rParameters["c1"].GetDouble();
    mC1 = rParameters["beta_zero"].GetDouble();
    mMinValue = rParameters["min_value"].GetDouble();

    const std::string& r_region_type = rParameters["region_type"].GetString();

    if (r_region_type == "assumes_both_regions")
    {
        mRegionType = RegionTypes::ASSUMES_BOTH_REGIONS;
        KRATOS_INFO(this->Info()) << "Wall function assumes either linear or "
                                     "log region is present in "
                                  << mModelPartName << ".\n";
    }
    else if (r_region_type == "assumes_log_region_only")
    {
        mRegionType = RegionTypes::ASSUMES_LOG_REGION_ONLY;
        KRATOS_INFO(this->Info())
            << "Wall function assumes only log region is present in "
            << mModelPartName << ".\n";
    }
    else
    {
        KRATOS_ERROR << "Unknown region type. Only supported region types are "
                        "\"assumes_both_regions\" and "
                        "\"assumes_log_region_only\". [ region_type = "
                     << r_region_type << " ].\n";
    }

    KRATOS_CATCH("");
}

RansNutYPlusWallFunctionUpdateProcess::RansNutYPlusWallFunctionUpdateProcess(
    Model& rModel,
    const std::string& rModelPartName,
    const double VonKarman,
    const double MinValue,
    const int EchoLevel)
    : mrModel(rModel),
      mModelPartName(rModelPartName),
      mVonKarman(VonKarman),
      mMinValue(MinValue),
      mEchoLevel(EchoLevel)
{
    mRegionType = RegionTypes::ASSUMES_LOG_REGION_ONLY;
    KRATOS_INFO(this->Info())
        << "Wall function assumes only log region is present in "
        << mModelPartName << ".\n";
    mBetaZero = 0.0;
    mC1 = 0.0;
}

RansNutYPlusWallFunctionUpdateProcess::RansNutYPlusWallFunctionUpdateProcess(
    Model& rModel,
    const std::string& rModelPartName,
    const double VonKarman,
    const double BetaZero,
    const double C1,
    const double MinValue,
    const int EchoLevel)
    : mrModel(rModel),
      mModelPartName(rModelPartName),
      mVonKarman(VonKarman),
      mBetaZero(BetaZero),
      mC1(C1),
      mMinValue(MinValue),
      mEchoLevel(EchoLevel)
{
    mRegionType = RegionTypes::ASSUMES_BOTH_REGIONS;
    KRATOS_INFO(this->Info())
        << "Wall function assumes either linear or log region is present in "
        << mModelPartName << ".\n";
}

int RansNutYPlusWallFunctionUpdateProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, KINEMATIC_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, VISCOSITY);

    return 0;

    KRATOS_CATCH("");
}

void RansNutYPlusWallFunctionUpdateProcess::ExecuteInitialize()
{
    CalculateConditionNeighbourCount();
}

void RansNutYPlusWallFunctionUpdateProcess::CalculateConditionNeighbourCount()
{
    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    VariableUtils().SetNonHistoricalVariableToZero(
        NUMBER_OF_NEIGHBOUR_CONDITIONS, r_model_part.Nodes());

    const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        ConditionType& r_cond = *(r_model_part.ConditionsBegin() + i_cond);
        ConditionGeometryType& r_geometry = r_cond.GetGeometry();
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
        {
            NodeType& r_node = r_geometry[i_node];
            r_node.SetLock();
            r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS) += 1;
            r_node.UnSetLock();
        }
    }

    r_model_part.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_CONDITIONS);

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Calculated number of neighbour conditions in " << mModelPartName << ".\n";
}

void RansNutYPlusWallFunctionUpdateProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    if (!mIsInitialized)
    {
        this->Execute();
        mIsInitialized = true;
    }

    KRATOS_CATCH("");
}

void RansNutYPlusWallFunctionUpdateProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    VariableUtils().SetHistoricalVariableToZero(TURBULENT_VISCOSITY,
                                                r_model_part.Nodes());

    const double y_plus_limit = r_model_part.GetProcessInfo()[RANS_Y_PLUS_LIMIT];

    const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        ConditionType& r_cond = *(r_model_part.ConditionsBegin() + i_cond);
        ConditionType::GeometryType& r_geometry = r_cond.GetGeometry();

        double y_plus = r_cond.GetValue(RANS_Y_PLUS);
        if (mRegionType == RegionTypes::ASSUMES_LOG_REGION_ONLY)
        {
            y_plus = std::max(y_plus, y_plus_limit);
        }

        if (y_plus >= y_plus_limit)
        {
            for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
            {
                NodeType& r_node = r_geometry[i_node];
                const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
                r_node.SetLock();
                r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) +=
                    mVonKarman * y_plus * nu;
                r_node.UnSetLock();
            }
        }
        else
        {
            for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
            {
                NodeType& r_node = r_geometry[i_node];
                const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
                r_node.SetLock();
                r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) +=
                    mC1 * mBetaZero * nu * std::pow(y_plus, 4) / 6.0;
                r_node.UnSetLock();
            }
        }
    }
    r_model_part.GetCommunicator().AssembleCurrentData(TURBULENT_VISCOSITY);

    const int number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
        double& r_nut = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        const double number_of_neighbour_conditions =
            static_cast<double>(r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS));
        r_nut = RansCalculationUtilities::SoftMax(
            r_nut / number_of_neighbour_conditions, mMinValue);
        r_node.FastGetSolutionStepValue(VISCOSITY) =
            r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) + r_nut;
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Calculated wall function based nu_t for " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

std::string RansNutYPlusWallFunctionUpdateProcess::Info() const
{
    return std::string("RansNutYPlusWallFunctionUpdateProcess");
}

void RansNutYPlusWallFunctionUpdateProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutYPlusWallFunctionUpdateProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
