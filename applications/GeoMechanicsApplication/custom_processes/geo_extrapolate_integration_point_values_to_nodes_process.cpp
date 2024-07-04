//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Jonathan Nuttall
//                   Wijtze Pieter Kikstra
//                   Richard Faasse

#include "custom_processes/geo_extrapolate_integration_point_values_to_nodes_process.h"
#include "containers/model.h"
#include "custom_utilities/linear_nodal_extrapolator.h"
#include "utilities/atomic_utilities.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
GeoExtrapolateIntegrationPointValuesToNodesProcess::GeoExtrapolateIntegrationPointValuesToNodesProcess(
    Model& rModel, Parameters ThisParameters)
    : GeoExtrapolateIntegrationPointValuesToNodesProcess(
          rModel.GetModelPart(ThisParameters["model_part_name"].GetString()), ThisParameters)
{
}

GeoExtrapolateIntegrationPointValuesToNodesProcess::GeoExtrapolateIntegrationPointValuesToNodesProcess(
    ModelPart& rMainModelPart, Parameters ThisParameters)
    : mrModelPart(rMainModelPart), mpExtrapolator(std::make_unique<LinearNodalExtrapolator>())
{
    const Parameters default_parameters =
        GeoExtrapolateIntegrationPointValuesToNodesProcess::GetDefaultParameters();
    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    FillVariableLists(ThisParameters);
}

GeoExtrapolateIntegrationPointValuesToNodesProcess::~GeoExtrapolateIntegrationPointValuesToNodesProcess() = default;

const Parameters GeoExtrapolateIntegrationPointValuesToNodesProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "model_part_name"            : "",
        "echo_level"                 : 0,
        "list_of_variables"          : []
    })");
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::FillVariableLists(const Parameters& rParameters)
{
    for (const std::string& r_variable_name : rParameters["list_of_variables"].GetStringArray()) {
        KRATOS_ERROR_IF_NOT(TryAddVariableToList(r_variable_name, mDoubleVariables) ||
                            TryAddVariableToList(r_variable_name, mArrayVariables) ||
                            TryAddVariableToList(r_variable_name, mVectorVariables) ||
                            TryAddVariableToList(r_variable_name, mMatrixVariables))
            << "Only double, array, vector and matrix variables are allowed in the "
               "variables list."
            << std::endl;
    }
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::Execute()
{
    ExecuteBeforeSolutionLoop();
    ExecuteFinalizeSolutionStep();
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::ExecuteBeforeSolutionLoop()
{
    InitializeAverageVariablesForElements();
    InitializeVectorAndMatrixZeros();
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::InitializeAverageVariablesForElements() const
{
    VariableUtils().SetNonHistoricalVariable(mrAverageVariable, 0.0, mrModelPart.Nodes());

    block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
        if (rElement.IsActive()) {
            for (auto& node : rElement.GetGeometry()) {
                auto& node_var_to_update = node.GetValue(mrAverageVariable);
                AtomicAdd(node_var_to_update, 1.0);
            }
        }
    });
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::InitializeVectorAndMatrixZeros()
{
    if (mrModelPart.Elements().empty()) return;

    auto first_element = mrModelPart.Elements().begin();

    const auto integration_points_number =
        first_element->GetGeometry().IntegrationPointsNumber(first_element->GetIntegrationMethod());

    const ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
    InitializeZerosOfVectorVariables(*first_element, integration_points_number, r_process_info);
    InitializeZerosOfMatrixVariables(*first_element, integration_points_number, r_process_info);
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::InitializeZerosOfVectorVariables(
    Element& rFirstElement, SizeType NumberOfIntegrationPoints, const ProcessInfo& rProcessInfo)
{
    for (const auto p_var : mVectorVariables) {
        std::vector<Vector> values_on_integration_points(NumberOfIntegrationPoints);
        rFirstElement.CalculateOnIntegrationPoints(*p_var, values_on_integration_points, rProcessInfo);
        mZeroValuesOfVectorVariables.try_emplace(p_var, ZeroVector(values_on_integration_points[0].size()));
    }
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::InitializeZerosOfMatrixVariables(
    Element& rFirstElement, SizeType NumberOfIntegrationPoints, const ProcessInfo& rProcessInfo)
{
    for (const auto p_var : mMatrixVariables) {
        std::vector<Matrix> values_on_integration_points(NumberOfIntegrationPoints);
        rFirstElement.CalculateOnIntegrationPoints(*p_var, values_on_integration_points, rProcessInfo);
        mZeroValuesOfMatrixVariables.try_emplace(
            p_var, ZeroMatrix(values_on_integration_points[0].size1(),
                              values_on_integration_points[0].size2()));
    }
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::ExecuteFinalizeSolutionStep()
{
    InitializeVariables();
    CacheExtrapolationMatricesForElements();

    block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
        if (rElement.IsActive()) {
            AddIntegrationPointContributionsForAllVariables(rElement, GetCachedExtrapolationMatrixFor(rElement));
        }
    });
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::CacheExtrapolationMatricesForElements()
{
    // This is specifically a single-thread range-based for-loop and no block_for_each.
    // Running this in parallel would lead to issues, since multiple threads might try to
    // write to the same cache at the same time.
    for (const auto& rElement : mrModelPart.Elements()) {
        if (!ExtrapolationMatrixIsCachedFor(rElement)) {
            CacheExtrapolationMatrixFor(
                rElement, mpExtrapolator->CalculateElementExtrapolationMatrix(
                              rElement.GetGeometry(), rElement.GetIntegrationMethod()));
        }
    }
}

bool GeoExtrapolateIntegrationPointValuesToNodesProcess::ExtrapolationMatrixIsCachedFor(const Element& rElement) const
{
    return mExtrapolationMatrixMap.count(typeid(rElement).hash_code()) > 0;
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::CacheExtrapolationMatrixFor(const Element& rElement,
                                                                                     const Matrix& rExtrapolationMatrix)
{
    mExtrapolationMatrixMap[typeid(rElement).hash_code()] = rExtrapolationMatrix;
}

const Matrix& GeoExtrapolateIntegrationPointValuesToNodesProcess::GetCachedExtrapolationMatrixFor(const Element& rElement) const
{
    return mExtrapolationMatrixMap.at(typeid(rElement).hash_code());
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::AddIntegrationPointContributionsForAllVariables(
    Element& rElem, const Matrix& rExtrapolationMatrix) const
{
    const SizeType integration_points_number =
        rElem.GetGeometry().IntegrationPointsNumber(rElem.GetIntegrationMethod());

    for (const auto p_var : mDoubleVariables) {
        AddIntegrationContributionsToNodes(rElem, *p_var, rExtrapolationMatrix,
                                           integration_points_number, AtomicAdd<double>);
    }
    for (const auto p_var : mArrayVariables) {
        AddIntegrationContributionsToNodes(rElem, *p_var, rExtrapolationMatrix,
                                           integration_points_number, AtomicAdd<double, 3>);
    }
    for (const auto p_var : mVectorVariables) {
        AddIntegrationContributionsToNodes(rElem, *p_var, rExtrapolationMatrix,
                                           integration_points_number, AtomicAddVector<Vector, Vector>);
    }
    for (const auto p_var : mMatrixVariables) {
        AddIntegrationContributionsToNodes(rElem, *p_var, rExtrapolationMatrix,
                                           integration_points_number, AtomicAddMatrix<Matrix, Matrix>);
    }
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::InitializeVariables()
{
    const array_1d<double, 3> zero_array = ZeroVector(3);

    block_for_each(mrModelPart.Nodes(), [this, zero_array](Node& rNode) {
        for (const auto p_var : mDoubleVariables) {
            rNode.FastGetSolutionStepValue(*p_var) = 0.0;
        }
        for (const auto p_var : mArrayVariables) {
            rNode.FastGetSolutionStepValue(*p_var) = zero_array;
        }
        for (const auto p_var : mVectorVariables) {
            rNode.FastGetSolutionStepValue(*p_var) = mZeroValuesOfVectorVariables[p_var];
        }
        for (const auto p_var : mMatrixVariables) {
            rNode.FastGetSolutionStepValue(*p_var) = mZeroValuesOfMatrixVariables[p_var];
        }
    });
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::ExecuteFinalize()
{
    block_for_each(mrModelPart.Nodes(),
                   [this](Node& rNode) { rNode.GetData().Erase(mrAverageVariable); });
}

std::string GeoExtrapolateIntegrationPointValuesToNodesProcess::Info() const
{
    return "GeoExtrapolateIntegrationPointValuesToNodesProcess";
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

} // namespace Kratos.
