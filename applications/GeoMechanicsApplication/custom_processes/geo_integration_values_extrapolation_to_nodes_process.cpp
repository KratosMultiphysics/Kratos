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
      mrAverageVariable(KratosComponents<Variable<double>>::Get(ThisParameters["average_variable"].GetString())),
      mpExtrapolator(std::make_unique<NodalExtrapolator>())
{
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mExtrapolateNonHistorical = ThisParameters["extrapolate_non_historical"].GetBool();
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
        "list_of_variables"          : [],
        "extrapolate_non_historical" : true
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

    for (const auto p_var : mDoubleVariable) {
        KRATOS_INFO("GeoIntegrationValuesExtrapolationToNodesProcess")
            << "Double variable: " << p_var->Name() << std::endl;
    }
    for (const auto p_var : mArrayVariable) {
        KRATOS_INFO("GeoIntegrationValuesExtrapolationToNodesProcess")
            << "Array variable: " << p_var->Name() << std::endl;
    }

    for (const auto p_var : mVectorVariable) {
        KRATOS_INFO("GeoIntegrationValuesExtrapolationToNodesProcess")
            << "Vector variable: " << p_var->Name() << std::endl;
    }
    for (const auto p_var : mMatrixVariable) {
        KRATOS_INFO("GeoIntegrationValuesExtrapolationToNodesProcess")
            << "Matrix variable: " << p_var->Name() << std::endl;
    }
}

void GeoIntegrationValuesExtrapolationToNodesProcess::Execute()
{
    ExecuteBeforeSolutionLoop();
    ExecuteFinalizeSolutionStep();
}

void GeoIntegrationValuesExtrapolationToNodesProcess::ExecuteBeforeSolutionLoop()
{
    // We initialize the average variable
    VariableUtils().SetNonHistoricalVariable(mrAverageVariable, 0.0, mrModelPart.Nodes());

    // We initialize the map of coincident and maps of sizes
    InitializeMaps();
}

void GeoIntegrationValuesExtrapolationToNodesProcess::InitializeMaps()
{
    // The list of elements
    auto& r_elements_array = mrModelPart.Elements();
    auto  it_elem_begin    = r_elements_array.begin();

    FillAverageVariableForElements();

    mrModelPart.GetCommunicator().AssembleNonHistoricalData(mrAverageVariable);

    // The process info
    const ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();

    // First we check if the model part contains at least one element
    if (r_elements_array.size() != 0) {
        // The first iterator of elements
        const auto& r_this_geometry_begin = it_elem_begin->GetGeometry();

        // Auxiliar values
        const GeometryData::IntegrationMethod this_integration_method = it_elem_begin->GetIntegrationMethod();
        const GeometryType::IntegrationPointsArrayType& integration_points =
            r_this_geometry_begin.IntegrationPoints(this_integration_method);
        const SizeType integration_points_number = integration_points.size();

        // We init the vector sizes
        for (const auto p_var : mVectorVariable) {
            std::vector<Vector> aux_result(integration_points_number);
            it_elem_begin->CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
            mSizeVectors.insert({p_var, aux_result[0].size()});
        }

        // We init the matrix sizes
        for (const auto p_var : mMatrixVariable) {
            std::vector<Matrix> aux_result(integration_points_number);
            it_elem_begin->CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
            std::pair<SizeType, SizeType> aux_pair(aux_result[0].size1(), aux_result[0].size2());
            mSizeMatrixes.insert({p_var, aux_pair});
        }
    }
}

void GeoIntegrationValuesExtrapolationToNodesProcess::FillAverageVariableForElements() const
{
    block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
        if (rElement.IsActive()) {
            for (auto& node : rElement.GetGeometry()) {
                auto& node_var_to_update = node.GetValue(mrAverageVariable);
                AtomicAdd(node_var_to_update, 1.0);
            }
        }
    });
}

void GeoIntegrationValuesExtrapolationToNodesProcess::ExecuteFinalizeSolutionStep()
{
    // We initialize the values
    InitializeVariables();

    // The process info
    const ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();

    // The list of elements
    auto& r_elements_array = mrModelPart.Elements();

    // Auxiliar values
    block_for_each(r_elements_array, [this, &r_process_info](Element& rElem) {
        // Only active elements
        if (rElem.IsActive()) {
            auto& r_this_geometry = rElem.GetGeometry();
            const GeometryData::IntegrationMethod this_integration_method = rElem.GetIntegrationMethod();
            auto integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const SizeType integration_points_number = integration_points.size();

            const Matrix extrapolation_matrix =
                GetExtrapolationMatrix(rElem, r_this_geometry, this_integration_method);

            for (const auto p_var : mDoubleVariable) {
                AddIntegrationContributionsToNodes(rElem, *p_var, r_process_info,
                                                   extrapolation_matrix, integration_points_number);
            }
            for (const auto p_var : mArrayVariable) {
                AddIntegrationContributionsToNodes(rElem, *p_var, r_process_info,
                                                   extrapolation_matrix, integration_points_number);
            }
            for (const auto p_var : mVectorVariable) {
                AddIntegrationContributionsToNodes(rElem, *p_var, r_process_info,
                                                   extrapolation_matrix, integration_points_number);
            }
            for (const auto p_var : mMatrixVariable) {
                AddIntegrationContributionsToNodes(rElem, *p_var, r_process_info,
                                                   extrapolation_matrix, integration_points_number);
            }
        }
    });

    // Assemble nodal data
    if (mExtrapolateNonHistorical) {
        for (const auto p_var : mDoubleVariable) {
            mrModelPart.GetCommunicator().AssembleNonHistoricalData(*p_var);
        }
        for (const auto p_var : mArrayVariable) {
            mrModelPart.GetCommunicator().AssembleNonHistoricalData(*p_var);
        }
        for (const auto p_var : mVectorVariable) {
            mrModelPart.GetCommunicator().AssembleNonHistoricalData(*p_var);
        }
        for (const auto p_var : mMatrixVariable) {
            mrModelPart.GetCommunicator().AssembleNonHistoricalData(*p_var);
        }
    } else {
        for (const auto p_var : mDoubleVariable) {
            mrModelPart.GetCommunicator().AssembleCurrentData(*p_var);
        }
        for (const auto p_var : mArrayVariable) {
            mrModelPart.GetCommunicator().AssembleCurrentData(*p_var);
        }
        for (const auto p_var : mVectorVariable) {
            mrModelPart.GetCommunicator().AssembleCurrentData(*p_var);
        }
        for (const auto p_var : mMatrixVariable) {
            mrModelPart.GetCommunicator().AssembleCurrentData(*p_var);
        }
    }
}

Matrix GeoIntegrationValuesExtrapolationToNodesProcess::GetExtrapolationMatrix(
    const Element&                         rElem,
    GeometricalObject::GeometryType&       r_this_geometry,
    const GeometryData::IntegrationMethod& this_integration_method)
{
    Matrix result;
    // check if the element type hash is in the extrapolation matrix map
    if (mExtrapolationMatrixMap.count(typeid(rElem).hash_code()) > 0) {
        result = mExtrapolationMatrixMap[typeid(rElem).hash_code()];
    } else {
        // calculate the extrapolation matrix
        result = mpExtrapolator->CalculateElementExtrapolationMatrix(r_this_geometry, this_integration_method);
        mExtrapolationMatrixMap[typeid(rElem).hash_code()] = result;
    }
    return result;
}

void GeoIntegrationValuesExtrapolationToNodesProcess::InitializeVariables()
{
    // Initializing some auxiliar values
    array_1d<double, 3> zero_array = ZeroVector(3);

    // The list of nodes
    auto& r_nodes_array = mrModelPart.Nodes();

    // Initialize values
    block_for_each(r_nodes_array, [this, zero_array](Node& rNode) {
        if (mExtrapolateNonHistorical) {
            // We initialize the doubles values
            for (const auto p_var : mDoubleVariable) {
                rNode.SetValue(*p_var, 0.0);
            }
            // We initialize the arrays values
            for (const auto p_var : mArrayVariable) {
                rNode.SetValue(*p_var, zero_array);
            }
            // We initialize the vectors values
            for (const auto p_var : mVectorVariable) {
                const Vector zero_vector = ZeroVector(mSizeVectors[p_var]);
                rNode.SetValue(*p_var, zero_vector);
            }
            // We initialize the matrix values
            for (const auto p_var : mMatrixVariable) {
                const Matrix zero_matrix =
                    ZeroMatrix(mSizeMatrixes[p_var].first, mSizeMatrixes[p_var].second);
                rNode.SetValue(*p_var, zero_matrix);
            }
        } else {
            // We initialize the doubles values
            for (const auto p_var : mDoubleVariable) {
                rNode.FastGetSolutionStepValue(*p_var) = 0.0;
            }
            // We initialize the arrays values
            for (const auto p_var : mArrayVariable) {
                rNode.FastGetSolutionStepValue(*p_var) = zero_array;
            }
            // We initialize the vectors values
            for (const auto p_var : mVectorVariable) {
                const Vector zero_vector               = ZeroVector(mSizeVectors[p_var]);
                rNode.FastGetSolutionStepValue(*p_var) = zero_vector;
            }
            // We initialize the matrix values
            for (const auto p_var : mMatrixVariable) {
                const Matrix zero_matrix =
                    ZeroMatrix(mSizeMatrixes[p_var].first, mSizeMatrixes[p_var].second);
                rNode.FastGetSolutionStepValue(*p_var) = zero_matrix;
            }
        }
    });
}

void GeoIntegrationValuesExtrapolationToNodesProcess::ExecuteFinalize()
{
    // The list of nodes
    auto& r_nodes_array = mrModelPart.Nodes();

    // Remove average variable
    block_for_each(r_nodes_array, [this](Node& rNode) {
        auto& data = rNode.GetData();
        data.Erase(mrAverageVariable);

        // We erase the doubles values
        for (const auto p_var : mDoubleVariable) {
            if (mExtrapolateNonHistorical) data.Erase(*p_var);
        }

        // We erase the arrays values
        for (const auto p_var : mArrayVariable) {
            if (mExtrapolateNonHistorical) data.Erase(*p_var);
        }

        // We erase the vectors values
        for (const auto p_var : mVectorVariable) {
            if (mExtrapolateNonHistorical) data.Erase(*p_var);
        }

        // We erase the matrix values
        for (const auto p_var : mMatrixVariable) {
            if (mExtrapolateNonHistorical) data.Erase(*p_var);
        }
    });
}

} // namespace Kratos.
