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
        "echo_level"                 : 1,
        "list_of_variables"          : [],
        "extrapolate_non_historical" : true
    })");
    ThisParameters.ValidateAndAssignDefaults(default_parameters);
    
    mEchoLevel = ThisParameters["echo_level"].GetInt();
    mExtrapolateNonHistorical = ThisParameters["extrapolate_non_historical"].GetBool();

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

    // We initialize the map of coincident and maps of sizes
    InitializeMaps();
}

/***********************************************************************************/
/***********************************************************************************/

void IntegrationValuesExtrapolationToNodesProcess::Execute()
{
    // We initialize the values
    InitializeVariables();

    // The list of elements
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i) {

    }
}

/***********************************************************************************/
/***********************************************************************************/

void IntegrationValuesExtrapolationToNodesProcess::InitializeMaps()
{
    // The list of nodes
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();

    // Initialize map
    for(IndexType i = 0; i < nodes_array.size(); ++i) {
        auto it_node = nodes_array.begin() + i;
        mCoincidentMap.insert({it_node->Id(), 0});
    }

    // The list of elements
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    // Fill the map
    for(IndexType i = 0; i < elements_array.size(); ++i) {
        auto it_elem = elements_array.begin() + i;
        for (auto& node : it_elem->GetGeometry()) {
            mCoincidentMap[node.Id()] += 1;
        }
    }

    // The process info
    const ProcessInfo& process_info = mrThisModelPart.GetProcessInfo();

    // The first iterator of elements
    auto it_elem_begin = elements_array.begin();
    auto& r_geo_begin = it_elem_begin->GetGeometry();

    // Auxiliar values
    const GeometryData::IntegrationMethod this_integration_method = it_elem_begin->GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geo_begin.IntegrationPoints(this_integration_method);
    const SizeType integration_points_number = integration_points.size();

    // We init the vector sizes
    for ( const auto& i_var : mVectorVariable) {
        std::vector<Vector> aux_result;
        it_elem_begin->GetValueOnIntegrationPoints(i_var, aux_result, process_info);
        mSizeVectors.insert({i_var, aux_result[0].size()});
    }

    // We init the matrix sizes
    for ( const auto& i_var : mMatrixVariable) {
        std::vector<Matrix> aux_result;
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
