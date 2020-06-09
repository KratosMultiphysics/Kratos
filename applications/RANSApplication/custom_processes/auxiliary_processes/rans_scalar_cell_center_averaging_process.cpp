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
#include "includes/define.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_scalar_cell_center_averaging_process.h"

namespace Kratos
{
RansScalarCellCenterAveragingProcess::RansScalarCellCenterAveragingProcess(Model& rModel,
                                                                           Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "echo_level"           : 0,
            "model_part_name"      : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "input_variable_name"  : "PLEASE_SPECIFY_AVERAGING_VARIABLE",
            "output_variable_name" : "PLEASE_SPECIFY_OUTPUT_VARIABLE"
        })");

    mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mModelPartName = mrParameters["model_part_name"].GetString();
    mAveragingInputVariableName = mrParameters["input_variable_name"].GetString();
    mAveragingOutputVariableName = mrParameters["output_variable_name"].GetString();
    mEchoLevel = mrParameters["echo_level"].GetInt();

    KRATOS_CATCH("");
}

int RansScalarCellCenterAveragingProcess::Check()
{
    KRATOS_TRY


    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const Variable<double>& input_variable =
        KratosComponents<Variable<double>>::Get(mAveragingInputVariableName);
    const Variable<double>& output_variable =
        KratosComponents<Variable<double>>::Get(mAveragingOutputVariableName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, input_variable);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, output_variable);

    return 0.0;

    KRATOS_CATCH("");
}

void RansScalarCellCenterAveragingProcess::ExecuteInitialize()
{
    FindNodalConditionNeighbourCount();
    CalculateCellCenterAverage();
}

void RansScalarCellCenterAveragingProcess::Execute()
{
    CalculateCellCenterAverage();
}

std::string RansScalarCellCenterAveragingProcess::Info() const
{
    return std::string("RansScalarCellCenterAveragingProcess");
}

void RansScalarCellCenterAveragingProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansScalarCellCenterAveragingProcess::PrintData(std::ostream& rOStream) const
{
}

void RansScalarCellCenterAveragingProcess::FindNodalConditionNeighbourCount()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    VariableUtils().SetNonHistoricalVariableToZero(
        NUMBER_OF_NEIGHBOUR_CONDITIONS, r_model_part.Nodes());

    const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
    for (int i_condition = 0; i_condition < number_of_conditions; ++i_condition)
    {
        Condition& r_condition = *(r_model_part.ConditionsBegin() + i_condition);
        Condition::GeometryType& r_geometry = r_condition.GetGeometry();
        const int number_of_nodes = r_geometry.PointsNumber();

        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = r_geometry[i_node];
            const int current_value = r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS);
            r_node.SetLock();
            r_node.SetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS, current_value + 1);
            r_node.UnSetLock();
        }
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Calculated nodal neighbour conditions count for nodes in "
        << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

void RansScalarCellCenterAveragingProcess::CalculateCellCenterAverage()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const Variable<double>& r_input_variable =
        KratosComponents<Variable<double>>::Get(mAveragingInputVariableName);
    const Variable<double>& r_output_variable =
        KratosComponents<Variable<double>>::Get(mAveragingOutputVariableName);

    VariableUtils variable_utilities;

    variable_utilities.SetHistoricalVariableToZero(r_output_variable,
                                                   r_model_part.Nodes());
    variable_utilities.SetNonHistoricalVariableToZero(r_output_variable,
                                                      r_model_part.Nodes());

    const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
    for (int i_condition = 0; i_condition < number_of_conditions; ++i_condition)
    {
        Condition& r_condition = *(r_model_part.ConditionsBegin() + i_condition);
        Condition::GeometryType& r_condition_geometry = r_condition.GetGeometry();
        const int number_of_nodes = r_condition_geometry.PointsNumber();

        KRATOS_ERROR_IF(!r_condition.Has(NEIGHBOUR_ELEMENTS))
            << "Parent element not found for condition id=" << r_condition.Id()
            << ". Please run \"FindConditionParentProcess\" for "
            << mModelPartName << ".\n";

        const Element::WeakPointer& p_parent_element = r_condition.GetValue(NEIGHBOUR_ELEMENTS)(0);
        const Element::GeometryType& r_parent_element_geometry =
            p_parent_element->GetGeometry();

        Vector gauss_weights;
        Matrix shape_functions;
        Element::GeometryType::ShapeFunctionsGradientsType shape_function_derivatives;

        RansCalculationUtilities::CalculateGeometryData(
            r_parent_element_geometry, GeometryData::GI_GAUSS_1, gauss_weights,
            shape_functions, shape_function_derivatives);

        const Vector& gauss_shape_functions = row(shape_functions, 0);
        const double value = RansCalculationUtilities::EvaluateInPoint(
            r_parent_element_geometry, r_input_variable, gauss_shape_functions);

        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = r_condition_geometry[i_node];

            const double current_output_value = r_node.GetValue(r_output_variable);
            r_node.SetLock();
            r_node.SetValue(r_output_variable, current_output_value + value);
            r_node.UnSetLock();
        }
    }

    const int number_of_nodes = r_model_part.NumberOfNodes();
    unsigned int number_of_calculated_nodes = 0;
#pragma omp parallel for reduction(+ : number_of_calculated_nodes)
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
        const double number_of_neighbouring_conditions =
            static_cast<double>(r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS));
        r_node.FastGetSolutionStepValue(r_output_variable) =
            r_node.GetValue(r_output_variable) / number_of_neighbouring_conditions;
        number_of_calculated_nodes++;
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Cell centered averaging is stored at " << r_output_variable.Name()
        << " using " << r_input_variable.Name() << " for "
        << number_of_calculated_nodes << " nodes in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

} // namespace Kratos.
