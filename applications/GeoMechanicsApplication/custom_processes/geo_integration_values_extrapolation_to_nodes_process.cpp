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

#include "custom_processes/geo_integration_values_extrapolation_to_nodes_process.h"
#include "containers/model.h"
#include "custom_utilities/nodal_extrapolator.h"
#include "utilities/atomic_utilities.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
GeoIntegrationValuesExtrapolationToNodesProcess::GeoIntegrationValuesExtrapolationToNodesProcess(Model& rModel, Parameters rParameters)
    : GeoIntegrationValuesExtrapolationToNodesProcess(
          rModel.GetModelPart(rParameters["model_part_name"].GetString()), rParameters)
{
}

GeoIntegrationValuesExtrapolationToNodesProcess::GeoIntegrationValuesExtrapolationToNodesProcess(
    ModelPart& rMainModelPart, Parameters rParameters)
    : mrModelPart(rMainModelPart),
      mrAverageVariable(KratosComponents<Variable<double>>::Get(rParameters["average_variable"].GetString())),
      mpExtrapolator(std::make_unique<NodalExtrapolator>())
{
    const Parameters default_parameters = GetDefaultParameters();
    rParameters.ValidateAndAssignDefaults(default_parameters);

    GetVariableLists(rParameters);
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
        KRATOS_ERROR_IF_NOT(TryAddVariableToList(r_variable_name, mDoubleVariables) ||
                            TryAddVariableToList(r_variable_name, mArrayVariables) ||
                            TryAddVariableToList(r_variable_name, mVectorVariables) ||
                            TryAddVariableToList(r_variable_name, mMatrixVariables))
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
    Element& rFirstElement, SizeType NumberOfIntegrationPoints, const ProcessInfo& rProcessInfo)
{
    for (const auto p_var : mVectorVariables) {
        std::vector<Vector> values_on_integration_points(NumberOfIntegrationPoints);
        rFirstElement.CalculateOnIntegrationPoints(*p_var, values_on_integration_points, rProcessInfo);
        mSizesOfVectorVariables.insert({p_var, values_on_integration_points[0].size()});
    }
}

void GeoIntegrationValuesExtrapolationToNodesProcess::InitializeSizesOfMatrixVariables(
    Element& rFirstElement, SizeType NumberOfIntegrationPoints, const ProcessInfo& rProcessInfo)
{
    for (const auto p_var : mMatrixVariables) {
        std::vector<Matrix> values_on_integration_points(NumberOfIntegrationPoints);
        rFirstElement.CalculateOnIntegrationPoints(*p_var, values_on_integration_points, rProcessInfo);
        std::pair<SizeType, SizeType> matrix_dimensions(values_on_integration_points[0].size1(),
                                                        values_on_integration_points[0].size2());
        mSizesOfMatrixVariables.try_emplace(p_var, matrix_dimensions);
    }
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
    AssembleNodalData(mDoubleVariables);
    AssembleNodalData(mArrayVariables);
    AssembleNodalData(mVectorVariables);
    AssembleNodalData(mMatrixVariables);
}

void GeoIntegrationValuesExtrapolationToNodesProcess::AddIntegrationContributionsForAllVariableLists(
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

Matrix GeoIntegrationValuesExtrapolationToNodesProcess::GetExtrapolationMatrix(
    const Element& rElement, GeometricalObject::GeometryType& rGeometry, const GeometryData::IntegrationMethod& rIntegrationMethod)
{
    if (!ExtrapolationMatrixIsCachedFor(rElement)) {
        CacheExtrapolationMatrixFor(
            rElement, mpExtrapolator->CalculateElementExtrapolationMatrix(rGeometry, rIntegrationMethod));
    }

    return GetCachedExtrapolationMatrixFor(rElement);
}

bool GeoIntegrationValuesExtrapolationToNodesProcess::ExtrapolationMatrixIsCachedFor(const Element& rElement) const
{
    return mExtrapolationMatrixMap.count(typeid(rElement).hash_code()) > 0;
}

Matrix GeoIntegrationValuesExtrapolationToNodesProcess::GetCachedExtrapolationMatrixFor(const Element& rElement)
{
    return mExtrapolationMatrixMap[typeid(rElement).hash_code()];
}

void GeoIntegrationValuesExtrapolationToNodesProcess::CacheExtrapolationMatrixFor(const Element& rElement,
                                                                                  const Matrix& rExtrapolationMatrix)
{
    mExtrapolationMatrixMap[typeid(rElement).hash_code()] = rExtrapolationMatrix;
}

void GeoIntegrationValuesExtrapolationToNodesProcess::InitializeVariables()
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
            rNode.FastGetSolutionStepValue(*p_var) = ZeroVector(mSizesOfVectorVariables[p_var]);
        }
        for (const auto p_var : mMatrixVariables) {
            rNode.FastGetSolutionStepValue(*p_var) =
                ZeroMatrix(mSizesOfMatrixVariables[p_var].first, mSizesOfMatrixVariables[p_var].second);
        }
    });
}

void GeoIntegrationValuesExtrapolationToNodesProcess::ExecuteFinalize()
{
    block_for_each(mrModelPart.Nodes(),
                   [this](Node& rNode) { rNode.GetData().Erase(mrAverageVariable); });
}

} // namespace Kratos.
