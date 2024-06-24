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
#include "custom_utilities/nodal_extrapolator.h"
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
    ModelPart& rMainModelPart, Parameters rParameters)
    : mrModelPart(rMainModelPart), mpExtrapolator(std::make_unique<NodalExtrapolator>())
{
    const Parameters default_parameters =
        GeoExtrapolateIntegrationPointValuesToNodesProcess::GetDefaultParameters();
    rParameters.ValidateAndAssignDefaults(default_parameters);

    FillVariableLists(rParameters);
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
    InitializeVectorAndMatrixSizesOfVariables();
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

    mrModelPart.GetCommunicator().AssembleNonHistoricalData(mrAverageVariable);
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::InitializeVectorAndMatrixSizesOfVariables()
{
    if (mrModelPart.Elements().empty()) return;

    auto first_element = mrModelPart.Elements().begin();

    const auto integration_points_number =
        first_element->GetGeometry().IntegrationPointsNumber(first_element->GetIntegrationMethod());

    const ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
    InitializeDefaultsOfVectorVariables(*first_element, integration_points_number, r_process_info);
    InitializeDefaultsOfMatrixVariables(*first_element, integration_points_number, r_process_info);
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::InitializeDefaultsOfVectorVariables(
    Element& rFirstElement, SizeType NumberOfIntegrationPoints, const ProcessInfo& rProcessInfo)
{
    for (const auto p_var : mVectorVariables) {
        std::vector<Vector> values_on_integration_points(NumberOfIntegrationPoints);
        rFirstElement.CalculateOnIntegrationPoints(*p_var, values_on_integration_points, rProcessInfo);
        mDefaultValuesOfVectorVariables.insert({p_var, ZeroVector(values_on_integration_points[0].size())});
    }
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::InitializeDefaultsOfMatrixVariables(
    Element& rFirstElement, SizeType NumberOfIntegrationPoints, const ProcessInfo& rProcessInfo)
{
    for (const auto p_var : mMatrixVariables) {
        std::vector<Matrix> values_on_integration_points(NumberOfIntegrationPoints);
        rFirstElement.CalculateOnIntegrationPoints(*p_var, values_on_integration_points, rProcessInfo);
        mDefaultValuesOfMatrixVariables.try_emplace(
            p_var, ZeroMatrix(values_on_integration_points[0].size1(),
                              values_on_integration_points[0].size2()));
    }
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::ExecuteFinalizeSolutionStep()
{
    InitializeVariables();

    block_for_each(mrModelPart.Elements(), [this](Element& rElem) {
        if (rElem.IsActive()) {
            auto& r_this_geometry = rElem.GetGeometry();
            const GeometryData::IntegrationMethod this_integration_method = rElem.GetIntegrationMethod();
            const Matrix extrapolation_matrix = GetExtrapolationMatrix(rElem);

            const SizeType integration_points_number =
                r_this_geometry.IntegrationPoints(this_integration_method).size();
            AddIntegrationContributionsForAllVariableLists(rElem, integration_points_number, extrapolation_matrix);
        }
    });

    AssembleNodalDataForAllVariableLists();
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::AssembleNodalDataForAllVariableLists()
{
    AssembleNodalData(mDoubleVariables);
    AssembleNodalData(mArrayVariables);
    AssembleNodalData(mVectorVariables);
    AssembleNodalData(mMatrixVariables);
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::AddIntegrationContributionsForAllVariableLists(
    Element& rElem, SizeType NumberOfIntegrationPoints, const Matrix& rExtrapolationMatrix)
{
    for (const auto p_var : mDoubleVariables) {
        AddIntegrationContributionsToNodes(rElem, *p_var, rExtrapolationMatrix,
                                           NumberOfIntegrationPoints, AtomicAdd<double>);
    }
    for (const auto p_var : mArrayVariables) {
        AddIntegrationContributionsToNodes(rElem, *p_var, rExtrapolationMatrix,
                                           NumberOfIntegrationPoints, AtomicAdd<double, 3>);
    }
    for (const auto p_var : mVectorVariables) {
        AddIntegrationContributionsToNodes(rElem, *p_var, rExtrapolationMatrix,
                                           NumberOfIntegrationPoints, AtomicAddVector<Vector, Vector>);
    }
    for (const auto p_var : mMatrixVariables) {
        AddIntegrationContributionsToNodes(rElem, *p_var, rExtrapolationMatrix,
                                           NumberOfIntegrationPoints, AtomicAddMatrix<Matrix, Matrix>);
    }
}

Matrix GeoExtrapolateIntegrationPointValuesToNodesProcess::GetExtrapolationMatrix(const Element& rElement) const
{
    if (!ExtrapolationMatrixIsCachedFor(rElement)) {
        CacheExtrapolationMatrixFor(rElement, mpExtrapolator->CalculateElementExtrapolationMatrix(
                                                  rElement.GetGeometry(), rElement.GetIntegrationMethod()));
    }

    return GetCachedExtrapolationMatrixFor(rElement);
}

bool GeoExtrapolateIntegrationPointValuesToNodesProcess::ExtrapolationMatrixIsCachedFor(const Element& rElement) const
{
    return mExtrapolationMatrixMap.count(typeid(rElement).hash_code()) > 0;
}

Matrix GeoExtrapolateIntegrationPointValuesToNodesProcess::GetCachedExtrapolationMatrixFor(const Element& rElement) const
{
    return mExtrapolationMatrixMap.at(typeid(rElement).hash_code());
}

void GeoExtrapolateIntegrationPointValuesToNodesProcess::CacheExtrapolationMatrixFor(const Element& rElement,
                                                                                     const Matrix& rExtrapolationMatrix) const
{
    mExtrapolationMatrixMap[typeid(rElement).hash_code()] = rExtrapolationMatrix;
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
            rNode.FastGetSolutionStepValue(*p_var) = mDefaultValuesOfVectorVariables[p_var];
        }
        for (const auto p_var : mMatrixVariables) {
            rNode.FastGetSolutionStepValue(*p_var) = mDefaultValuesOfMatrixVariables[p_var];
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
