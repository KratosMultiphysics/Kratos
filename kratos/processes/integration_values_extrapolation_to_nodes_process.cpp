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
//

// System includes

// External includes

// Project includes

// Include the point locator
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"
#include "processes/integration_values_extrapolation_to_nodes_process.h"

namespace Kratos
{
IntegrationValuesExtrapolationToNodesProcess::IntegrationValuesExtrapolationToNodesProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : IntegrationValuesExtrapolationToNodesProcess(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()), ThisParameters)
{
}

/***********************************************************************************/
/***********************************************************************************/

IntegrationValuesExtrapolationToNodesProcess::IntegrationValuesExtrapolationToNodesProcess(
    ModelPart& rMainModelPart,
    Parameters ThisParameters
    ) : mrModelPart(rMainModelPart)
{
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = ThisParameters["echo_level"].GetInt();
    mExtrapolateNonHistorical = ThisParameters["extrapolate_non_historical"].GetBool();
    mAreaAverage = ThisParameters["area_average"].GetBool();

    // The average variable
    mpAverageVariable = &(KratosComponents<Variable<double>>::Get(ThisParameters["average_variable"].GetString()));

    // We get the list of variables
    for (const std::string& r_variable_name : ThisParameters["list_of_variables"].GetStringArray()){
        if (KratosComponents<Variable<double>>::Has(r_variable_name)){
            mDoubleVariable.push_back(&(KratosComponents< Variable<double>>::Get(r_variable_name)));
        } else if (KratosComponents<Variable<array_1d<double, 3> >>::Has(r_variable_name)) {
            mArrayVariable.push_back(&(KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name)));
        } else if (KratosComponents<Variable<Vector>>::Has(r_variable_name)) {
            mVectorVariable.push_back(&(KratosComponents<Variable<Vector>>::Get(r_variable_name)));
        } else if (KratosComponents<Variable<Matrix>>::Has(r_variable_name)) {
            mMatrixVariable.push_back(&(KratosComponents<Variable<Matrix>>::Get(r_variable_name)));
        } else {
            KRATOS_ERROR << "Only double, array, vector and matrix variables are allowed in the variables list." ;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void IntegrationValuesExtrapolationToNodesProcess::Execute()
{
    // We execute all the necessary steps
    ExecuteBeforeSolutionLoop();
    ExecuteFinalizeSolutionStep();
//     ExecuteFinalize();
}

/***********************************************************************************/
/***********************************************************************************/

void IntegrationValuesExtrapolationToNodesProcess::ExecuteBeforeSolutionLoop()
{
    // We initialize the average variable
    VariableUtils().SetNonHistoricalVariable(*mpAverageVariable, 0.0, mrModelPart.Nodes());

    // We initialize the map of coincident and maps of sizes
    InitializeMaps();
}

/***********************************************************************************/
/***********************************************************************************/

void IntegrationValuesExtrapolationToNodesProcess::ExecuteFinalizeSolutionStep()
{
    // We initialize the values
    InitializeVariables();

    // The process info
    const ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();

    // The list of elements
    auto& r_elements_array = mrModelPart.Elements();

    // Auxiliar value
    struct TLSType
    {
        Vector vector_J, N;
    };

    block_for_each(r_elements_array, TLSType(), [&](Element& rElem, TLSType& rTls){
        // Only active elements
        if (rElem.IsActive()) {
            auto& r_this_geometry = rElem.GetGeometry();

            // Auxiliar values
            const GeometryData::IntegrationMethod this_integration_method = rElem.GetIntegrationMethod();
            const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const SizeType integration_points_number = integration_points.size();
            const SizeType number_of_nodes = r_this_geometry.size();

            // Prepare matrix of coeffients
            Matrix node_coefficient(number_of_nodes, integration_points_number);

            // Point elements only have one node
            if (r_this_geometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Point) {
                // Definition of node coefficient
                rTls.vector_J = r_this_geometry.DeterminantOfJacobian(rTls.vector_J , this_integration_method );
                if (rTls.N.size() != number_of_nodes )
                    rTls.N.resize(number_of_nodes);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    const double area_coeff = mAreaAverage ? integration_points[i_gauss_point].Weight() * rTls.vector_J[i_gauss_point] : 1.0;
                    const array_1d<double, 3>& r_local_coordinates = integration_points[i_gauss_point].Coordinates();
                    r_this_geometry.ShapeFunctionsValues( rTls.N, r_local_coordinates );
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        const double average_variable_value = r_this_geometry[i_node].GetValue(*mpAverageVariable);
                        const double coeff_coincident_node = std::abs(average_variable_value) > std::numeric_limits<double>::epsilon() ? area_coeff/average_variable_value : area_coeff;
                        node_coefficient(i_node, i_gauss_point) = coeff_coincident_node * std::abs(rTls.N[i_node]);
                    }
                }
            } else { // Point geometry, only one node
                const double gp_coefficient = 1.0/static_cast<double>(integration_points_number);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    node_coefficient(0, i_gauss_point) = gp_coefficient;
                }
            }

            // We add the doubles values
            for ( const auto p_var : mDoubleVariable) {
                std::vector<double> aux_result(integration_points_number);
                rElem.CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        double& aux_value = (mExtrapolateNonHistorical) ? r_this_geometry[i_node].GetValue(*p_var) : r_this_geometry[i_node].FastGetSolutionStepValue(*p_var);
                        AtomicAdd(aux_value, node_coefficient(i_node, i_gauss_point) * aux_result[i_gauss_point]);
                    }
                }
            }

            // We add the arrays values
            for ( const auto p_var : mArrayVariable) {
                std::vector<array_1d<double, 3>> aux_result(integration_points_number);
                rElem.CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        array_1d<double, 3>& aux_value = (mExtrapolateNonHistorical) ? r_this_geometry[i_node].GetValue(*p_var) : r_this_geometry[i_node].FastGetSolutionStepValue(*p_var);
                        const array_1d<double, 3>& aux_sol = node_coefficient(i_node, i_gauss_point) * aux_result[i_gauss_point];
                        AtomicAdd(aux_value, aux_sol);
                    }
                }
            }

            // We add the vectors values
            for ( const auto p_var : mVectorVariable) {
                std::vector<Vector> aux_result(integration_points_number);
                rElem.CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        Vector& aux_value = (mExtrapolateNonHistorical) ? r_this_geometry[i_node].GetValue(*p_var) : r_this_geometry[i_node].FastGetSolutionStepValue(*p_var);
                        const Vector& aux_sol = node_coefficient(i_node, i_gauss_point) * aux_result[i_gauss_point];
                        for (IndexType i_comp = 0; i_comp < aux_sol.size(); ++i_comp) {
                            AtomicAdd(aux_value[i_comp], aux_sol[i_comp]);
                        }
                    }
                }
            }

            // We add the matrix values
            for ( const auto p_var : mMatrixVariable) {
                std::vector<Matrix> aux_result(integration_points_number);
                rElem.CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        Matrix& aux_value = (mExtrapolateNonHistorical) ? r_this_geometry[i_node].GetValue(*p_var) : r_this_geometry[i_node].FastGetSolutionStepValue(*p_var);
                        const Matrix& aux_sol = node_coefficient(i_node, i_gauss_point) * aux_result[i_gauss_point];
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
    if (mExtrapolateNonHistorical)
    {
        for (const auto p_var : mDoubleVariable) {
            mrModelPart.GetCommunicator().AssembleNonHistoricalData(*p_var);
        }
        for (const auto p_var : mArrayVariable)  {
            mrModelPart.GetCommunicator().AssembleNonHistoricalData(*p_var);
        }
        for (const auto p_var : mVectorVariable) {
            mrModelPart.GetCommunicator().AssembleNonHistoricalData(*p_var);
        }
        for (const auto p_var : mMatrixVariable) {
            mrModelPart.GetCommunicator().AssembleNonHistoricalData(*p_var);
        }
    }
    else
    {
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

/***********************************************************************************/
/***********************************************************************************/

void IntegrationValuesExtrapolationToNodesProcess::ExecuteFinalize()
{
    // The list of nodes
    auto& r_nodes_array = mrModelPart.Nodes();

    // Remove average variable
    block_for_each(r_nodes_array, [&](Node& rNode){
        auto& data = rNode.GetData();
        data.Erase(*mpAverageVariable);

        // We erase the doubles values
        for ( const auto p_var : mDoubleVariable) {
            if (mExtrapolateNonHistorical) data.Erase(*p_var);
        }

        // We erase the arrays values
        for ( const auto p_var : mArrayVariable) {
            if (mExtrapolateNonHistorical) data.Erase(*p_var);
        }

        // We erase the vectors values
        for ( const auto p_var : mVectorVariable) {
            if (mExtrapolateNonHistorical) data.Erase(*p_var);
        }

        // We erase the matrix values
        for ( const auto p_var : mMatrixVariable) {
            if (mExtrapolateNonHistorical) data.Erase(*p_var);
        }
    });
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters IntegrationValuesExtrapolationToNodesProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"            : "",
        "echo_level"                 : 0,
        "area_average"               : true,
        "average_variable"           : "NODAL_AREA",
        "list_of_variables"          : [],
        "extrapolate_non_historical" : true
    })" );
    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

void IntegrationValuesExtrapolationToNodesProcess::InitializeMaps()
{
    // The list of elements
    auto& r_elements_array = mrModelPart.Elements();
    auto it_elem_begin = r_elements_array.begin();

    // Some definitions
    struct TLSType
    {
        Vector vector_J, N;
    };

    // Fill the average value
    block_for_each(r_elements_array, TLSType(), [&](Element& rElem, TLSType& rTls){
        // Only active elements
        if (rElem.IsActive()) {
            // The geometry of the element
            auto& r_this_geometry = rElem.GetGeometry();

            // Point elements only have one node
            if (r_this_geometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Point) {
                // Auxiliar values
                const GeometryData::IntegrationMethod this_integration_method = rElem.GetIntegrationMethod();
                const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const SizeType integration_points_number = integration_points.size();
                const SizeType number_of_nodes = r_this_geometry.size();

                // The jacobian of the geometry
                rTls.vector_J = r_this_geometry.DeterminantOfJacobian(rTls.vector_J , this_integration_method );
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    const array_1d<double, 3>& r_local_coordinates = integration_points[i_gauss_point].Coordinates();
                    if (rTls.N.size() != number_of_nodes )
                        rTls.N.resize(number_of_nodes);
                    r_this_geometry.ShapeFunctionsValues( rTls.N, r_local_coordinates );
                    const double area_coeff = mAreaAverage ? integration_points[i_gauss_point].Weight() * rTls.vector_J[i_gauss_point] : 1.0;
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        AtomicAdd(r_this_geometry[i_node].GetValue(*mpAverageVariable), std::abs(rTls.N[i_node]) * area_coeff);
                    }
                }
            }
        }
    });

    mrModelPart.GetCommunicator().AssembleNonHistoricalData(*mpAverageVariable);

    // The process info
    const ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();

    // First we check if the model part contains at least one element
    if (r_elements_array.size() != 0) {
        // The first iterator of elements
        auto& r_this_geometry_begin = it_elem_begin->GetGeometry();

        // Auxiliar values
        const GeometryData::IntegrationMethod this_integration_method = it_elem_begin->GetIntegrationMethod();
        const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry_begin.IntegrationPoints(this_integration_method);
        const SizeType integration_points_number = integration_points.size();

        // We init the vector sizes
        for ( const auto p_var : mVectorVariable) {
            std::vector<Vector> aux_result(integration_points_number);
            it_elem_begin->CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
            mSizeVectors.insert({p_var, aux_result[0].size()});
        }

        // We init the matrix sizes
        for ( const auto p_var : mMatrixVariable) {
            std::vector<Matrix> aux_result(integration_points_number);
            it_elem_begin->CalculateOnIntegrationPoints(*p_var, aux_result, r_process_info);
            std::pair<SizeType, SizeType> aux_pair(aux_result[0].size1(), aux_result[0].size2());
            mSizeMatrixes.insert({p_var, aux_pair});
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void IntegrationValuesExtrapolationToNodesProcess::InitializeVariables()
{
    // Initializing some auxiliar values
    array_1d<double, 3> zero_array = ZeroVector(3);

    // The list of nodes
    auto& r_nodes_array = mrModelPart.Nodes();

    // Initialize values
    block_for_each(r_nodes_array, [&](Node& rNode){
        if (mExtrapolateNonHistorical)
        {
            // We initialize the doubles values
            for ( const auto p_var : mDoubleVariable) {
                rNode.SetValue(*p_var, 0.0);
            }
            // We initialize the arrays values
            for ( const auto p_var : mArrayVariable) {
                rNode.SetValue(*p_var, zero_array);
            }
            // We initialize the vectors values
            for ( const auto p_var : mVectorVariable) {
                const Vector zero_vector = ZeroVector(mSizeVectors[p_var]);
                rNode.SetValue(*p_var, zero_vector);
            }
            // We initialize the matrix values
            for ( const auto p_var : mMatrixVariable) {
                const Matrix zero_matrix = ZeroMatrix(mSizeMatrixes[p_var].first, mSizeMatrixes[p_var].second);
                rNode.SetValue(*p_var, zero_matrix);
            }
        }
        else
        {
            // We initialize the doubles values
            for ( const auto p_var : mDoubleVariable) {
                rNode.FastGetSolutionStepValue(*p_var) = 0.0;
            }
            // We initialize the arrays values
            for ( const auto p_var : mArrayVariable) {
                rNode.FastGetSolutionStepValue(*p_var) = zero_array;
            }
            // We initialize the vectors values
            for ( const auto p_var : mVectorVariable) {
                const Vector zero_vector = ZeroVector(mSizeVectors[p_var]);
                rNode.FastGetSolutionStepValue(*p_var) = zero_vector;
            }
            // We initialize the matrix values
            for ( const auto p_var : mMatrixVariable) {
                const Matrix zero_matrix = ZeroMatrix(mSizeMatrixes[p_var].first, mSizeMatrixes[p_var].second);
                rNode.FastGetSolutionStepValue(*p_var) = zero_matrix;
            }
        }
    });
}

}  // namespace Kratos.
