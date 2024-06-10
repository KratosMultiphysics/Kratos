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

    mExtrapolateNonHistorical = ThisParameters["extrapolate_non_historical"].GetBool();
    GetVariableLists(ThisParameters);
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
    // We initialize the average variable
    VariableUtils().SetNonHistoricalVariable(mrAverageVariable, 0.0, mrModelPart.Nodes());

    // We initialize the map of coincident and maps of sizes
    InitializeMaps();
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
    block_for_each(r_elements_array, TLSType(), [this, &r_process_info](Element& rElem, TLSType& rTls) {
        // Only active elements
        if (rElem.IsActive()) {
            auto& r_this_geometry = rElem.GetGeometry();
            const GeometryData::IntegrationMethod this_integration_method = rElem.GetIntegrationMethod();
            auto integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const SizeType integration_points_number = integration_points.size();
            const SizeType number_of_nodes           = r_this_geometry.size();
            Matrix         extrapolation_matrix;
            // check if the element type hash is in the extrapolation matrix map

            if (mExtrapolationMatrixMap.count(typeid(rElem).hash_code()) > 0) {
                extrapolation_matrix = mExtrapolationMatrixMap[typeid(rElem).hash_code()];
            } else {
                // calculate the extrapolation matrix
                extrapolation_matrix = CalculateElementExtrapolationMatrix(
                    rElem, r_this_geometry, integration_points_number, integration_points,
                    this_integration_method, number_of_nodes, rTls);
                mExtrapolationMatrixMap[typeid(rElem).hash_code()] = extrapolation_matrix;
            }

            // We add the doubles values
            for (const auto p_var : mDoubleVariable) {
                std::vector<double> aux_result(integration_points_number);
                rElem.CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
                for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                    double aux_sol = 0.;
                    // double aux_sol = inner_prod( row(extrapolation_matrix, i_node), aux_result);
                    for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                        aux_sol += extrapolation_matrix(i_node, i_gauss_point) * aux_result[i_gauss_point];
                    }
                    aux_sol /= r_this_geometry[i_node].GetValue(mrAverageVariable);
                    auto& aux_value = mExtrapolateNonHistorical
                                          ? r_this_geometry[i_node].GetValue(*p_var)
                                          : r_this_geometry[i_node].FastGetSolutionStepValue(*p_var);
                    AtomicAdd(aux_value, aux_sol);
                }
            }

            // We add the arrays values
            for (const auto p_var : mArrayVariable) {
                std::vector<array_1d<double, 3>> aux_result(integration_points_number);
                rElem.CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        double nodal_area = r_this_geometry[i_node].GetValue(mrAverageVariable);
                        array_1d<double, 3>& aux_value =
                            mExtrapolateNonHistorical
                                ? r_this_geometry[i_node].GetValue(*p_var)
                                : r_this_geometry[i_node].FastGetSolutionStepValue(*p_var);
                        const array_1d<double, 3>& aux_sol = extrapolation_matrix(i_node, i_gauss_point) *
                                                             aux_result[i_gauss_point] / nodal_area;
                        AtomicAdd(aux_value, aux_sol);
                    }
                }
            }

            // We add the vectors values
            for (const auto p_var : mVectorVariable) {
                std::vector<Vector> aux_result(integration_points_number);
                rElem.CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        double  nodal_area    = r_this_geometry[i_node].GetValue(mrAverageVariable);
                        Vector& aux_value     = mExtrapolateNonHistorical
                                                    ? r_this_geometry[i_node].GetValue(*p_var)
                                                    : r_this_geometry[i_node].FastGetSolutionStepValue(*p_var);
                        const Vector& aux_sol = extrapolation_matrix(i_node, i_gauss_point) *
                                                aux_result[i_gauss_point] / nodal_area;
                        for (IndexType i_comp = 0; i_comp < aux_sol.size(); ++i_comp) {
                            AtomicAdd(aux_value[i_comp], aux_sol[i_comp]);
                        }
                    }
                }
            }

            // We add the matrix values
            for (const auto p_var : mMatrixVariable) {
                std::vector<Matrix> aux_result(integration_points_number);
                rElem.CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        double  nodal_area    = r_this_geometry[i_node].GetValue(mrAverageVariable);
                        Matrix& aux_value     = mExtrapolateNonHistorical
                                                    ? r_this_geometry[i_node].GetValue(*p_var)
                                                    : r_this_geometry[i_node].FastGetSolutionStepValue(*p_var);
                        const Matrix& aux_sol = extrapolation_matrix(i_node, i_gauss_point) *
                                                aux_result[i_gauss_point] / nodal_area;
                        for (IndexType i_comp = 0; i_comp < aux_sol.size1(); ++i_comp) {
                            for (IndexType j_comp = 0; j_comp < aux_sol.size2(); ++j_comp) {
                                AtomicAdd(aux_value(i_comp, j_comp), aux_sol(i_comp, j_comp));
                            }
                        }
                    }
                }
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

Matrix GeoIntegrationValuesExtrapolationToNodesProcess::CalculateElementExtrapolationMatrix(
    [[maybe_unused]] Element&                 rElem,
    GeometryType&                             r_this_geometry,
    SizeType                                  integration_points_number,
    GeometryType::IntegrationPointsArrayType& integration_points,
    GeometryData::IntegrationMethod           this_integration_method,
    SizeType                                  number_of_nodes,
    TLSType&                                  rTls) const
{
    KRATOS_ERROR_IF(r_this_geometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Triangle &&
                    r_this_geometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Quadrilateral);

    KRATOS_ERROR_IF(r_this_geometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Triangle &&
                    (number_of_nodes != 3 && number_of_nodes != 6));

    KRATOS_ERROR_IF(r_this_geometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Quadrilateral &&
                    (number_of_nodes != 4 && number_of_nodes != 8));

    // Sofar this works for 3, 4, 6 and 8 node planar elements
    // for 2 and 3 node line elements the extension is straightforward.
    // for volume elements ( hexa, tetra, wedge ) the midside node interpolation step is more elaborate
    GeometryType*                 p_low_order_geometry = &r_this_geometry;
    std::unique_ptr<GeometryType> p_new_low_order_geometry;
    switch (r_this_geometry.PointsNumber()) {
    case 6:
        p_new_low_order_geometry = std::make_unique<Triangle2D3<Node>>(
            r_this_geometry(0), r_this_geometry(1), r_this_geometry(2));
        p_low_order_geometry = p_new_low_order_geometry.get();
        break;
    case 8:
        p_new_low_order_geometry = std::make_unique<Quadrilateral2D4<Node>>(
            r_this_geometry(0), r_this_geometry(1), r_this_geometry(2), r_this_geometry(3));
        p_low_order_geometry = p_new_low_order_geometry.get();
        break;
    default:
        // deliberately empty
        break;
    }

    // calculate extrapolation matrix towards corner nodes
    SizeType number_of_low_order_nodes = p_low_order_geometry->PointsNumber();
    rTls.N.resize(number_of_low_order_nodes);
    Matrix quasi_mass_mat = ZeroMatrix(number_of_low_order_nodes, number_of_low_order_nodes);
    Matrix node_coefficient(number_of_low_order_nodes, integration_points_number);
    rTls.vector_J = r_this_geometry.DeterminantOfJacobian(rTls.vector_J, this_integration_method);
    for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
        // local_coordinates --> isoparametric coordinates or for triangles area coordinates
        const array_1d<double, 3>& r_local_coordinates = integration_points[i_gauss_point].Coordinates();
        // shape function for this i.p.
        p_low_order_geometry->ShapeFunctionsValues(rTls.N, r_local_coordinates);
        quasi_mass_mat += outer_prod(rTls.N, rTls.N) * rTls.vector_J[i_gauss_point] *
                          integration_points[i_gauss_point].Weight();
        column(node_coefficient, i_gauss_point) =
            rTls.N * rTls.vector_J[i_gauss_point] * integration_points[i_gauss_point].Weight();
    }

    double MetricDet;
    Matrix quasi_mass_mat_inverse;
    MathUtils<double>::InvertMatrix(quasi_mass_mat, quasi_mass_mat_inverse, MetricDet, -1.);

    auto extrapolation_matrix = Matrix{prod(quasi_mass_mat_inverse, node_coefficient)};

    // add rows for midside nodes if needed ( mean value of corner nodes )
    if (r_this_geometry.PointsNumber() > p_low_order_geometry->PointsNumber()) {
        extrapolation_matrix.resize(r_this_geometry.PointsNumber(), extrapolation_matrix.size2());
        std::size_t n_filled = p_low_order_geometry->PointsNumber();
        for (std::size_t i_row = 0; i_row < r_this_geometry.PointsNumber() - n_filled; ++i_row) {
            row(extrapolation_matrix, n_filled + i_row) =
                0.5 * (row(extrapolation_matrix, i_row) + row(extrapolation_matrix, (i_row + 1) % n_filled));
        }
    }

    return extrapolation_matrix;
}

void GeoIntegrationValuesExtrapolationToNodesProcess::ExecuteFinalize()
{
    // The list of nodes
    auto& r_nodes_array = mrModelPart.Nodes();

    // Remove average variable
    block_for_each(r_nodes_array, [&](Node& rNode) {
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

void GeoIntegrationValuesExtrapolationToNodesProcess::InitializeMaps()
{
    // The list of elements
    auto& r_elements_array = mrModelPart.Elements();
    auto  it_elem_begin    = r_elements_array.begin();

    // Some definitions
    struct TLSType {
        Vector vector_J;
        Vector N;
    };

    // Fill the average value
    block_for_each(r_elements_array, [this](Element& rElement) {
        // Only active elements
        if (rElement.IsActive()) {
            // The geometry of the element
            auto& r_this_geometry = rElement.GetGeometry();

            const SizeType number_of_nodes = r_this_geometry.size();
            for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                auto& node_var_to_update = r_this_geometry[i_node].GetValue(mrAverageVariable);
                AtomicAdd(node_var_to_update, 1.0);
            }
        }
    });

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

void GeoIntegrationValuesExtrapolationToNodesProcess::InitializeVariables()
{
    // Initializing some auxiliar values
    array_1d<double, 3> zero_array = ZeroVector(3);

    // The list of nodes
    auto& r_nodes_array = mrModelPart.Nodes();

    // Initialize values
    block_for_each(r_nodes_array, [&](Node& rNode) {
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

} // namespace Kratos.
