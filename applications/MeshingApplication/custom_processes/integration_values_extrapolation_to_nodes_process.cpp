// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes

// Include the point locator
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"
#include "custom_processes/integration_values_extrapolation_to_nodes_process.h"

namespace Kratos
{
IntegrationValuesExtrapolationToNodesProcess::IntegrationValuesExtrapolationToNodesProcess(
        ModelPart& rMainModelPart,
        Parameters ThisParameters
        ) : mrThisModelPart(rMainModelPart)
{
    Parameters default_parameters = Parameters(R"(
    {
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
    mAverageVariable = KratosComponents< Variable<double> >::Get(ThisParameters["average_variable"].GetString());

    // We get the list of variables
    const SizeType n_variables = ThisParameters["list_of_variables"].size();

    for (IndexType i_var = 0; i_var < n_variables; ++i_var){
        const std::string& variable_name = ThisParameters["list_of_variables"].GetArrayItem(i_var).GetString();

        if(KratosComponents<Variable<double>>::Has(variable_name)){
            Variable<double> variable = KratosComponents< Variable<double> >::Get(variable_name);
            mDoubleVariable.push_back(variable);
        } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(variable_name)) {
            Variable<array_1d< double, 3>> variable = KratosComponents< Variable<array_1d< double, 3>> >::Get(variable_name);
            mArrayVariable.push_back(variable);
        } else if (KratosComponents< Variable< Vector > >::Has(variable_name)) {
            Variable<Vector> variable = KratosComponents< Variable<Vector> >::Get(variable_name);
            mVectorVariable.push_back(variable);
        } else if (KratosComponents< Variable< Matrix > >::Has(variable_name)) {
            Variable<Matrix> variable = KratosComponents< Variable<Matrix> >::Get(variable_name);
            mMatrixVariable.push_back(variable);
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
    VariableUtils().SetNonHistoricalVariable(mAverageVariable, 0.0, mrThisModelPart.Nodes());

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
    const ProcessInfo& process_info = mrThisModelPart.GetProcessInfo();

    // The list of elements
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i) {
        auto it_elem = elements_array.begin() + i;

        // Only active elements. Detect if the element is active or not. If the user did not make any choice the element
        // NOTE: Is active by default
        const bool element_is_active = ((it_elem)->IsDefined(ACTIVE)) ? (it_elem)->Is(ACTIVE) : true;
        if (element_is_active) {
            auto& r_this_geometry = it_elem->GetGeometry();

            // Auxiliar values
            const GeometryData::IntegrationMethod this_integration_method = it_elem->GetIntegrationMethod();
            const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const SizeType integration_points_number = integration_points.size();
            const SizeType number_of_nodes = r_this_geometry.size();

            // Definition of node coefficient
            Matrix node_coefficient(number_of_nodes, integration_points_number);
            if (integration_points_number == 1) { // In case of just one GP the extrapolation it is just one
                const array_1d<double, 3>& local_coordinates = integration_points[0].Coordinates();
                Vector N( number_of_nodes );
                r_this_geometry.ShapeFunctionsValues( N, local_coordinates );
                for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                    node_coefficient(i_node, 0) = N[i_node];
                }
            } else { // Otherwise we need to build a matrix to invert or approximate the inverse
                Matrix shape_function_matrix(integration_points_number, number_of_nodes);
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    const array_1d<double, 3>& local_coordinates = integration_points[i_gauss_point].Coordinates();
                    Vector N( number_of_nodes );
                    r_this_geometry.ShapeFunctionsValues( N, local_coordinates );
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        shape_function_matrix(i_gauss_point, i_node) = N[i_node];
                    }
                }
                if (integration_points_number == number_of_nodes) {
                    double det;
                    MathUtils<double>::InvertMatrix(shape_function_matrix, node_coefficient, det);
                } else { // TODO: Try to use the QR utility
                    KRATOS_WARNING_IF("IntegrationValuesExtrapolationToNodesProcess", mEchoLevel > 0) << "Number of integration points higher than the number of nodes in element: " << it_elem->Id() << ". This is costly and could lose accuracy" << std::endl;
                    double det;
                    MathUtils<double>::GeneralizedInvertMatrix(shape_function_matrix, node_coefficient, det);
                }
            }

            // We add the doubles values
            for ( const auto& i_var : mDoubleVariable) {
                std::vector<double> aux_result(integration_points_number);
                it_elem->GetValueOnIntegrationPoints(i_var, aux_result, process_info);
                Vector vector_J;
                vector_J = r_this_geometry.DeterminantOfJacobian(vector_J , this_integration_method );
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    const double area_coeff = mAreaAverage ? integration_points[i_gauss_point].Weight() * vector_J[i_gauss_point] : 1.0;
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        double& aux_value = (mExtrapolateNonHistorical) ? r_this_geometry[i_node].GetValue(i_var) : r_this_geometry[i_node].FastGetSolutionStepValue(i_var);
                        #pragma omp atomic
                        aux_value += area_coeff * node_coefficient(i_node, i_gauss_point) * aux_result[i_gauss_point];
                    }
                }
            }

            // We add the arrays values
            for ( const auto& i_var : mArrayVariable) {
                std::vector<array_1d<double, 3>> aux_result(integration_points_number);
                it_elem->GetValueOnIntegrationPoints(i_var, aux_result, process_info);
                Vector vector_J;
                vector_J = r_this_geometry.DeterminantOfJacobian(vector_J , this_integration_method );
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    const double area_coeff = mAreaAverage ? integration_points[i_gauss_point].Weight() * vector_J[i_gauss_point] : 1.0;
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        array_1d<double, 3>& aux_value = (mExtrapolateNonHistorical) ? r_this_geometry[i_node].GetValue(i_var) : r_this_geometry[i_node].FastGetSolutionStepValue(i_var);
                        const array_1d<double, 3>& aux_sol = area_coeff * node_coefficient(i_node, i_gauss_point) * aux_result[i_gauss_point];
                        for (IndexType i_comp = 0; i_comp < 3; ++i_comp) {
                            #pragma omp atomic
                            aux_value[i_comp] += aux_sol[i_comp];
                        }
                    }
                }
            }

            // We add the vectors values
            for ( const auto& i_var : mVectorVariable) {
                std::vector<Vector> aux_result(integration_points_number);
                it_elem->GetValueOnIntegrationPoints(i_var, aux_result, process_info);
                Vector vector_J;
                vector_J = r_this_geometry.DeterminantOfJacobian(vector_J , this_integration_method );
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    const double area_coeff = mAreaAverage ? integration_points[i_gauss_point].Weight() * vector_J[i_gauss_point] : 1.0;
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        Vector& aux_value = (mExtrapolateNonHistorical) ? r_this_geometry[i_node].GetValue(i_var) : r_this_geometry[i_node].FastGetSolutionStepValue(i_var);
                        const Vector& aux_sol = area_coeff * node_coefficient(i_node, i_gauss_point) * aux_result[i_gauss_point];
                        for (IndexType i_comp = 0; i_comp < aux_sol.size(); ++i_comp) {
                            #pragma omp atomic
                            aux_value[i_comp] += aux_sol[i_comp];
                        }
                    }
                }
            }

            // We add the matrix values
            for ( const auto& i_var : mMatrixVariable) {
                std::vector<Matrix> aux_result(integration_points_number);
                it_elem->GetValueOnIntegrationPoints(i_var, aux_result, process_info);
                Vector vector_J;
                vector_J = r_this_geometry.DeterminantOfJacobian(vector_J , this_integration_method );
                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                    const double area_coeff = mAreaAverage ? integration_points[i_gauss_point].Weight() * vector_J[i_gauss_point] : 1.0;
                    for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                        Matrix& aux_value = (mExtrapolateNonHistorical) ? r_this_geometry[i_node].GetValue(i_var) : r_this_geometry[i_node].FastGetSolutionStepValue(i_var);
                        const Matrix& aux_sol = area_coeff * node_coefficient(i_node, i_gauss_point) * aux_result[i_gauss_point];
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

    // The list of nodes
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();

    // We ponderate the values
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;

        const double average_variable_value = it_node->GetValue(mAverageVariable);
        const double coeff_coincident_node = average_variable_value > std::numeric_limits<double>::epsilon() ? 1.0/it_node->GetValue(mAverageVariable) : 1.0;

        // We initialize the doubles values
        for ( const auto& i_var : mDoubleVariable) {
            if (mExtrapolateNonHistorical) it_node->GetValue(i_var) *= coeff_coincident_node;
            else it_node->FastGetSolutionStepValue(i_var) *= coeff_coincident_node;
        }

        // We initialize the arrays values
        for ( const auto& i_var : mArrayVariable) {
            if (mExtrapolateNonHistorical) it_node->GetValue(i_var) *= coeff_coincident_node;
            else it_node->FastGetSolutionStepValue(i_var) *= coeff_coincident_node;
        }

        // We initialize the vectors values
        for ( const auto& i_var : mVectorVariable) {
            if (mExtrapolateNonHistorical) it_node->GetValue(i_var) *= coeff_coincident_node;
            else it_node->FastGetSolutionStepValue(i_var) *= coeff_coincident_node;
        }

        // We initialize the matrix values
        for ( const auto& i_var : mMatrixVariable) {
            if (mExtrapolateNonHistorical) it_node->GetValue(i_var) *= coeff_coincident_node;
            else it_node->FastGetSolutionStepValue(i_var) *= coeff_coincident_node;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void IntegrationValuesExtrapolationToNodesProcess::ExecuteFinalize()
{
    // The list of nodes
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();

    // Remove average variable
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        auto& data = it_node->Data();
        data.Erase(mAverageVariable);

        // We erase the doubles values
        for ( const auto& i_var : mDoubleVariable) {
            if (mExtrapolateNonHistorical) data.Erase(i_var);
        }

        // We erase the arrays values
        for ( const auto& i_var : mArrayVariable) {
            if (mExtrapolateNonHistorical) data.Erase(i_var);
        }

        // We erase the vectors values
        for ( const auto& i_var : mVectorVariable) {
            if (mExtrapolateNonHistorical) data.Erase(i_var);
        }

        // We erase the matrix values
        for ( const auto& i_var : mMatrixVariable) {
            if (mExtrapolateNonHistorical) data.Erase(i_var);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void IntegrationValuesExtrapolationToNodesProcess::InitializeMaps()
{
    // The list of elements
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    // Fill the average value
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i) {
        auto it_elem = elements_array.begin() + i;
        // Only active elements. Detect if the element is active or not. If the user did not make any choice the element
        // NOTE: Is active by default
        const bool element_is_active = ((it_elem)->IsDefined(ACTIVE)) ? (it_elem)->Is(ACTIVE) : true;
        if (element_is_active) {
            // The geometry of the element
            auto& r_this_geometry = it_elem->GetGeometry();

            // Auxiliar values
            const GeometryData::IntegrationMethod this_integration_method = it_elem->GetIntegrationMethod();
            const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const SizeType integration_points_number = integration_points.size();
            const SizeType number_of_nodes = r_this_geometry.size();

            // The jacobian of the geometry
            Vector vector_J;
            vector_J = r_this_geometry.DeterminantOfJacobian(vector_J , this_integration_method );
            for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point) {
                const array_1d<double, 3>& local_coordinates = integration_points[i_gauss_point].Coordinates();
                Vector N( number_of_nodes );
                r_this_geometry.ShapeFunctionsValues( N, local_coordinates );
                const double area_coeff = mAreaAverage ? integration_points[i_gauss_point].Weight() * vector_J[i_gauss_point] : 1.0;
                for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                    #pragma omp atomic
                    r_this_geometry[i_node].GetValue(mAverageVariable) += N[i_node] * area_coeff;
                }
            }
        }
    }

    // The process info
    const ProcessInfo& process_info = mrThisModelPart.GetProcessInfo();

    // The first iterator of elements
    auto it_elem_begin = elements_array.begin();
    auto& r_this_geometry_begin = it_elem_begin->GetGeometry();

    // Auxiliar values
    const GeometryData::IntegrationMethod this_integration_method = it_elem_begin->GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry_begin.IntegrationPoints(this_integration_method);
    const SizeType integration_points_number = integration_points.size();

    // We init the vector sizes
    for ( const auto& i_var : mVectorVariable) {
        std::vector<Vector> aux_result(integration_points_number);
        it_elem_begin->GetValueOnIntegrationPoints(i_var, aux_result, process_info);
        mSizeVectors.insert({i_var, aux_result[0].size()});
    }

    // We init the matrix sizes
    for ( const auto& i_var : mMatrixVariable) {
        std::vector<Matrix> aux_result(integration_points_number);
        it_elem_begin->GetValueOnIntegrationPoints(i_var, aux_result, process_info);
        std::pair<SizeType, SizeType> aux_pair(aux_result[0].size1(), aux_result[0].size2());
        mSizeMatrixes.insert({i_var, aux_pair});
    }
}

/***********************************************************************************/
/***********************************************************************************/

void IntegrationValuesExtrapolationToNodesProcess::InitializeVariables()
{
    // Initializing some auxiliar values
    array_1d<double, 3> zero_array(3, 0.0);

    // The list of nodes
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();

    // Initialize values
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;

        // We initialize the doubles values
        for ( const auto& i_var : mDoubleVariable) {
            if (mExtrapolateNonHistorical) it_node->SetValue(i_var, 0.0);
            else it_node->FastGetSolutionStepValue(i_var) = 0.0;
        }

        // We initialize the arrays values
        for ( const auto& i_var : mArrayVariable) {
            if (mExtrapolateNonHistorical) it_node->SetValue(i_var, zero_array);
            else it_node->FastGetSolutionStepValue(i_var) = zero_array;
        }

        // We initialize the vectors values
        for ( const auto& i_var : mVectorVariable) {
            const Vector zero_vector = ZeroVector(mSizeVectors[i_var]);
            if (mExtrapolateNonHistorical) it_node->SetValue(i_var, zero_vector);
            else it_node->FastGetSolutionStepValue(i_var) = zero_vector;
        }

        // We initialize the matrix values
        for ( const auto& i_var : mMatrixVariable) {
            const Matrix zero_matrix = ZeroMatrix(mSizeMatrixes[i_var].first, mSizeMatrixes[i_var].second);
            if (mExtrapolateNonHistorical) it_node->SetValue(i_var, zero_matrix);
            else it_node->FastGetSolutionStepValue(i_var) = zero_matrix;
        }
    }
}

}  // namespace Kratos.
