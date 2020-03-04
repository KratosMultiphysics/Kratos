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
    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"            : "",
        "echo_level"                 : 0,
        "area_average"               : true,
        "average_variable"           : "NODAL_AREA",
        "list_of_variables"          : [],
        "extrapolate_non_historical" : true
    })");

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = ThisParameters["echo_level"].GetInt();
    mExtrapolateNonHistorical = ThisParameters["extrapolate_non_historical"].GetBool();
    mAreaAverage = ThisParameters["area_average"].GetBool();

    // The average variable
    mpAverageVariable = &(KratosComponents<Variable<double>>::Get(ThisParameters["average_variable"].GetString()));

    // We get the list of variables
    const SizeType n_variables = ThisParameters["list_of_variables"].size();

    for (IndexType p_var = 0; p_var < n_variables; ++p_var){
        const std::string& r_variable_name = ThisParameters["list_of_variables"].GetArrayItem(p_var).GetString();

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
    // We execute all the necesary steps
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
    const auto it_elem_begin = r_elements_array.begin();

    // Auxiliar value
    Vector vector_J, N;

    #pragma omp parallel for private(vector_J, N)
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = it_elem_begin + i;

        // Only active elements. Detect if the element is active or not. If the user did not make any choice the element
        // NOTE: Is active by default
        const bool element_is_active = it_elem->IsDefined(ACTIVE) ? it_elem->Is(ACTIVE) : true;
        if (element_is_active) {
            auto& r_this_geometry = it_elem->GetGeometry();

            // Auxiliar values
            const GeometryData::IntegrationMethod this_integration_method = it_elem->GetIntegrationMethod();
            const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const SizeType integration_points_number = integration_points.size();
            const SizeType number_of_nodes = r_this_geometry.size();

            // Definition of node coefficient
            vector_J = r_this_geometry.DeterminantOfJacobian(vector_J , this_integration_method );
            Matrix node_coefficient(number_of_nodes, integration_points_number);
            if (N.size() != number_of_nodes )
                N.resize(number_of_nodes);
            for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                const double area_coeff = mAreaAverage ? integration_points[i_gauss_point].Weight() * vector_J[i_gauss_point] : 1.0;
                const array_1d<double, 3>& r_local_coordinates = integration_points[i_gauss_point].Coordinates();
                r_this_geometry.ShapeFunctionsValues( N, r_local_coordinates );
                for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                    const double average_variable_value = r_this_geometry[i_node].GetValue(*mpAverageVariable);
                    const double coeff_coincident_node = std::abs(average_variable_value) > std::numeric_limits<double>::epsilon() ? area_coeff/average_variable_value : area_coeff;
                    node_coefficient(i_node, i_gauss_point) = coeff_coincident_node * std::abs(N[i_node]);
                }
            }

            // We add the doubles values
            for ( const auto p_var : mDoubleVariable) {
                std::vector<double> aux_result(integration_points_number);
                it_elem->GetValueOnIntegrationPoints(*p_var, aux_result, r_process_info);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        double& aux_value = (mExtrapolateNonHistorical) ? r_this_geometry[i_node].GetValue(*p_var) : r_this_geometry[i_node].FastGetSolutionStepValue(*p_var);
                        #pragma omp atomic
                        aux_value += node_coefficient(i_node, i_gauss_point) * aux_result[i_gauss_point];
                    }
                }
            }

            // We add the arrays values
            for ( const auto p_var : mArrayVariable) {
                std::vector<array_1d<double, 3>> aux_result(integration_points_number);
                it_elem->GetValueOnIntegrationPoints(*p_var, aux_result, r_process_info);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        array_1d<double, 3>& aux_value = (mExtrapolateNonHistorical) ? r_this_geometry[i_node].GetValue(*p_var) : r_this_geometry[i_node].FastGetSolutionStepValue(*p_var);
                        const array_1d<double, 3>& aux_sol = node_coefficient(i_node, i_gauss_point) * aux_result[i_gauss_point];
                        for (IndexType i_comp = 0; i_comp < 3; ++i_comp) {
                            #pragma omp atomic
                            aux_value[i_comp] += aux_sol[i_comp];
                        }
                    }
                }
            }

            // We add the vectors values
            for ( const auto p_var : mVectorVariable) {
                std::vector<Vector> aux_result(integration_points_number);
                it_elem->GetValueOnIntegrationPoints(*p_var, aux_result, r_process_info);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        Vector& aux_value = (mExtrapolateNonHistorical) ? r_this_geometry[i_node].GetValue(*p_var) : r_this_geometry[i_node].FastGetSolutionStepValue(*p_var);
                        const Vector& aux_sol = node_coefficient(i_node, i_gauss_point) * aux_result[i_gauss_point];
                        for (IndexType i_comp = 0; i_comp < aux_sol.size(); ++i_comp) {
                            #pragma omp atomic
                            aux_value[i_comp] += aux_sol[i_comp];
                        }
                    }
                }
            }

            // We add the matrix values
            for ( const auto p_var : mMatrixVariable) {
                std::vector<Matrix> aux_result(integration_points_number);
                it_elem->GetValueOnIntegrationPoints(*p_var, aux_result, r_process_info);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        Matrix& aux_value = (mExtrapolateNonHistorical) ? r_this_geometry[i_node].GetValue(*p_var) : r_this_geometry[i_node].FastGetSolutionStepValue(*p_var);
                        const Matrix& aux_sol = node_coefficient(i_node, i_gauss_point) * aux_result[i_gauss_point];
                        for (IndexType i_comp = 0; i_comp < aux_sol.size1(); ++i_comp) {
                            for (IndexType j_comp = 0; j_comp < aux_sol.size2(); ++j_comp) {
                                #pragma omp atomic
                                aux_value(i_comp, j_comp) += aux_sol(i_comp, j_comp);
                            }
                        }
                    }
                }
            }
        }
    }

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
    const auto it_node_begin = r_nodes_array.begin();

    // Remove average variable
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;
        auto& data = it_node->Data();
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
    }
}

/***********************************************************************************/
/***********************************************************************************/

void IntegrationValuesExtrapolationToNodesProcess::InitializeMaps()
{
    // The list of elements
    auto& r_elements_array = mrModelPart.Elements();
    auto it_elem_begin = r_elements_array.begin();

    // Some definitions
    Vector vector_J, N;

    // Fill the average value
    #pragma omp parallel for private(vector_J, N)
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = it_elem_begin + i;
        // Only active elements. Detect if the element is active or not. If the user did not make any choice the element
        // NOTE: Is active by default
        const bool element_is_active = it_elem->IsDefined(ACTIVE) ? it_elem->Is(ACTIVE) : true;
        if (element_is_active) {
            // The geometry of the element
            auto& r_this_geometry = it_elem->GetGeometry();

            // Auxiliar values
            const GeometryData::IntegrationMethod this_integration_method = it_elem->GetIntegrationMethod();
            const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const SizeType integration_points_number = integration_points.size();
            const SizeType number_of_nodes = r_this_geometry.size();

            // The jacobian of the geometry
            vector_J = r_this_geometry.DeterminantOfJacobian(vector_J , this_integration_method );
            for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                const array_1d<double, 3>& r_local_coordinates = integration_points[i_gauss_point].Coordinates();
                if (N.size() != number_of_nodes )
                    N.resize(number_of_nodes);
                r_this_geometry.ShapeFunctionsValues( N, r_local_coordinates );
                const double area_coeff = mAreaAverage ? integration_points[i_gauss_point].Weight() * vector_J[i_gauss_point] : 1.0;
                for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                    #pragma omp atomic
                    r_this_geometry[i_node].GetValue(*mpAverageVariable) += std::abs(N[i_node]) * area_coeff;
                }
            }
        }
    }

    mrModelPart.GetCommunicator().AssembleNonHistoricalData(*mpAverageVariable);

    // The process info
    const ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();

    // First we check if the model part constains at least one element
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
            it_elem_begin->GetValueOnIntegrationPoints(*p_var, aux_result, r_process_info);
            mSizeVectors.insert({p_var, aux_result[0].size()});
        }

        // We init the matrix sizes
        for ( const auto p_var : mMatrixVariable) {
            std::vector<Matrix> aux_result(integration_points_number);
            it_elem_begin->GetValueOnIntegrationPoints(*p_var, aux_result, r_process_info);
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
    const auto it_node_begin = r_nodes_array.begin();

    // Initialize values
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;

        // We initialize the doubles values
        for ( const auto p_var : mDoubleVariable) {
            if (mExtrapolateNonHistorical) it_node->SetValue(*p_var, 0.0);
            else it_node->FastGetSolutionStepValue(*p_var) = 0.0;
        }

        // We initialize the arrays values
        for ( const auto p_var : mArrayVariable) {
            if (mExtrapolateNonHistorical) it_node->SetValue(*p_var, zero_array);
            else it_node->FastGetSolutionStepValue(*p_var) = zero_array;
        }

        // We initialize the vectors values
        for ( const auto p_var : mVectorVariable) {
            const Vector zero_vector = ZeroVector(mSizeVectors[p_var]);
            if (mExtrapolateNonHistorical) it_node->SetValue(*p_var, zero_vector);
            else it_node->FastGetSolutionStepValue(*p_var) = zero_vector;
        }

        // We initialize the matrix values
        for ( const auto p_var : mMatrixVariable) {
            const Matrix zero_matrix = ZeroMatrix(mSizeMatrixes[p_var].first, mSizeMatrixes[p_var].second);
            if (mExtrapolateNonHistorical) it_node->SetValue(*p_var, zero_matrix);
            else it_node->FastGetSolutionStepValue(*p_var) = zero_matrix;
        }
    }
}

}  // namespace Kratos.
