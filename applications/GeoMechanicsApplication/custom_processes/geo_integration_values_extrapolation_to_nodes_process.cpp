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
//

// Project includes
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/triangle_2d_3.h"

// Include the point locator
#include "custom_processes/geo_integration_values_extrapolation_to_nodes_process.h"
#include "custom_utilities/nodal_extrapolator.h"
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
GeoIntegrationValuesExtrapolationToNodesProcess::GeoIntegrationValuesExtrapolationToNodesProcess(Model& rModel, Parameters ThisParameters)
    : GeoIntegrationValuesExtrapolationToNodesProcess(
          rModel.GetModelPart(ThisParameters["model_part_name"].GetString()), ThisParameters)
{
}

GeoIntegrationValuesExtrapolationToNodesProcess::GeoIntegrationValuesExtrapolationToNodesProcess(
    ModelPart& rMainModelPart, Parameters ThisParameters)
    : mrModelPart(rMainModelPart),
      mrAverageVariable(KratosComponents<Variable<double>>::Get(ThisParameters["average_variable"].GetString()))
{
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    GetVariableLists(ThisParameters);
}

GeoIntegrationValuesExtrapolationToNodesProcess::~GeoIntegrationValuesExtrapolationToNodesProcess() = default;

const Parameters GeoIntegrationValuesExtrapolationToNodesProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "model_part_name"            : "",
        "echo_level"                 : 0,
        "average_variable"           : "NODAL_AREA",
        "list_of_variables"          : []
    })");
}

void GeoIntegrationValuesExtrapolationToNodesProcess::GetVariableLists(const Parameters& rParameters)
{
    for (const std::string& r_variable_name : rParameters["list_of_variables"].GetStringArray()) {
        KRATOS_ERROR_IF_NOT(TryAddVariableToList(r_variable_name, mDoubleVariable) ||
                            TryAddVariableToList(r_variable_name, mArrayVariable) ||
                            TryAddVariableToList(r_variable_name, mVectorVariable) ||
                            TryAddVariableToList(r_variable_name, mMatrixVariable))
            << "Only double, array, vector and matrix variables are allowed in the "
               "variables list."
            << std::endl;
    }
}

void GeoIntegrationValuesExtrapolationToNodesProcess::Execute()
{
    ExecuteBeforeSolutionLoop();
    ExecuteFinalizeSolutionStep();
}

void GeoIntegrationValuesExtrapolationToNodesProcess::ExecuteBeforeSolutionLoop()
{
    InitializeAverageVariablesForElements();
    InitializeVectorAndMatrixSizesOfVariables();
}

void GeoIntegrationValuesExtrapolationToNodesProcess::InitializeVectorAndMatrixSizesOfVariables()
{
    if (ModelPartContainsAtLeastOneElement()) {
        auto r_first_element = mrModelPart.Elements().begin();

        const GeometryData::IntegrationMethod this_integration_method =
            r_first_element->GetIntegrationMethod();
        const GeometryType::IntegrationPointsArrayType& integration_points =
            r_first_element->GetGeometry().IntegrationPoints(this_integration_method);
        const SizeType integration_points_number = integration_points.size();

        const ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        InitializeSizesOfVectorVariables(*r_first_element, integration_points_number, r_process_info);
        InitializeSizesOfMatrixVariables(*r_first_element, integration_points_number, r_process_info);
    }
}

bool GeoIntegrationValuesExtrapolationToNodesProcess::ModelPartContainsAtLeastOneElement() const
{
    return !mrModelPart.Elements().empty();
}

void GeoIntegrationValuesExtrapolationToNodesProcess::InitializeSizesOfVectorVariables(
    Element& r_first_element, SizeType integration_points_number, const ProcessInfo& r_process_info)
{
    for (const auto p_var : mVectorVariable) {
        std::vector<Vector> values_on_integration_points(integration_points_number);
        r_first_element.CalculateOnIntegrationPoints(*p_var, values_on_integration_points, r_process_info);
        mSizesOfVectorVariables.insert({p_var, values_on_integration_points[0].size()});
    }
}

void GeoIntegrationValuesExtrapolationToNodesProcess::InitializeSizesOfMatrixVariables(
    Element& r_first_element, SizeType integration_points_number, const ProcessInfo& r_process_info)
{
    for (const auto p_var : mMatrixVariable) {
        std::vector<Matrix> values_on_integration_points(integration_points_number);
        r_first_element.CalculateOnIntegrationPoints(*p_var, values_on_integration_points, r_process_info);
        std::pair<SizeType, SizeType> matrix_dimensions(values_on_integration_points[0].size1(),
                                                        values_on_integration_points[0].size2());
        mSizesOfMatrixVariables.try_emplace(p_var, matrix_dimensions);
    }
}

void GeoIntegrationValuesExtrapolationToNodesProcess::InitializeAverageVariablesForElements() const
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

void GeoIntegrationValuesExtrapolationToNodesProcess::ExecuteFinalizeSolutionStep()
{
    InitializeVariables();

    block_for_each(mrModelPart.Elements(), [this](Element& rElem) {
        if (rElem.IsActive()) {
            auto& r_this_geometry = rElem.GetGeometry();
            const GeometryData::IntegrationMethod this_integration_method = rElem.GetIntegrationMethod();
            const Matrix extrapolation_matrix =
                GetExtrapolationMatrix(rElem, r_this_geometry, this_integration_method);

            const SizeType integration_points_number =
                r_this_geometry.IntegrationPoints(this_integration_method).size();
            AddIntegrationContributionsForAllVariableLists(rElem, integration_points_number, extrapolation_matrix);
        }
    });

    AssembleNodalDataForAllVariableLists();
}

void GeoIntegrationValuesExtrapolationToNodesProcess::AssembleNodalDataForAllVariableLists()
{
    AssembleNodalData(mDoubleVariable);
    AssembleNodalData(mArrayVariable);
    AssembleNodalData(mVectorVariable);
    AssembleNodalData(mMatrixVariable);
}

void GeoIntegrationValuesExtrapolationToNodesProcess::AddIntegrationContributionsForAllVariableLists(
    Element& rElem, const SizeType integration_points_number, const Matrix& extrapolation_matrix)
{
    for (const auto p_var : mDoubleVariable) {
        AddIntegrationContributionsToNodes(rElem, *p_var, extrapolation_matrix, integration_points_number);
    }
    for (const auto p_var : mArrayVariable) {
        AddIntegrationContributionsToNodes(rElem, *p_var, extrapolation_matrix, integration_points_number);
    }
    for (const auto p_var : mVectorVariable) {
        AddIntegrationContributionsToNodes(rElem, *p_var, extrapolation_matrix, integration_points_number);
    }
    for (const auto p_var : mMatrixVariable) {
        AddIntegrationContributionsToNodes(rElem, *p_var, extrapolation_matrix, integration_points_number);
    }
}

Matrix GeoIntegrationValuesExtrapolationToNodesProcess::GetExtrapolationMatrix(
    const Element&                         rElem,
    GeometricalObject::GeometryType&       r_this_geometry,
    const GeometryData::IntegrationMethod& this_integration_method)
{
    if (!ExtrapolationMatrixIsCachedFor(rElem)) {
        CacheExtrapolationMatrixFor(rElem, mpExtrapolator->CalculateElementExtrapolationMatrix(
                                               r_this_geometry, this_integration_method));
    }

    return GetCachedExtrapolationMatrixFor(rElem);
}

bool GeoIntegrationValuesExtrapolationToNodesProcess::ExtrapolationMatrixIsCachedFor(const Element& rElem) const
{
    return mExtrapolationMatrixMap.count(typeid(rElem).hash_code()) > 0;
}

Matrix GeoIntegrationValuesExtrapolationToNodesProcess::GetCachedExtrapolationMatrixFor(const Element& rElem)
{
    return mExtrapolationMatrixMap[typeid(rElem).hash_code()];
}

void GeoIntegrationValuesExtrapolationToNodesProcess::CacheExtrapolationMatrixFor(const Element& rElem,
                                                                                  const Matrix& rExtrapolationMatrix)
{
    mExtrapolationMatrixMap[typeid(rElem).hash_code()] = rExtrapolationMatrix;
}

void GeoIntegrationValuesExtrapolationToNodesProcess::InitializeVariables()
{
    const array_1d<double, 3> zero_array = ZeroVector(3);

    block_for_each(mrModelPart.Nodes(), [this, zero_array](Node& rNode) {
        for (const auto p_var : mDoubleVariable) {
            rNode.FastGetSolutionStepValue(*p_var) = 0.0;
        }
        for (const auto p_var : mArrayVariable) {
            rNode.FastGetSolutionStepValue(*p_var) = zero_array;
        }
        for (const auto p_var : mVectorVariable) {
            const Vector zero_vector               = ZeroVector(mSizesOfVectorVariables[p_var]);
            rNode.FastGetSolutionStepValue(*p_var) = zero_vector;
        }
        for (const auto p_var : mMatrixVariable) {
            const Matrix zero_matrix =
                ZeroMatrix(mSizesOfMatrixVariables[p_var].first, mSizesOfMatrixVariables[p_var].second);
            rNode.FastGetSolutionStepValue(*p_var) = zero_matrix;
        }
    });
}

void GeoIntegrationValuesExtrapolationToNodesProcess::ExecuteFinalize()
{
    block_for_each(mrModelPart.Nodes(),
                   [this](Node& rNode) { rNode.GetData().Erase(mrAverageVariable); });
}

} // namespace Kratos.
