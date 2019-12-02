//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Saransh Saxena
//

// System includes

// External includes

// Project includes
#include "custom_processes/simple_error_calculator_process.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
template<SizeType TDim>
SimpleErrorCalculatorProcess<TDim>::SimpleErrorCalculatorProcess(ModelPart& rThisModelPart, Parameters ThisParameters):mrThisModelPart(rThisModelPart)
{
    //WIP
    Parameters default_parameters = Parameters(R"(
    {
        "minimal_size"                        : 0.01,
        "maximal_size"                        : 10.0,
        "refinement_strategy"                 : "Simple_Error_Calculator",
        "reference_variable_name"             : "ERROR_RATIO",
        "echo_level"                          : 0
    })"
    );

    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    mMinSize = ThisParameters["minimal_size"].GetDouble();
    mMaxSize = ThisParameters["maximal_size"].GetDouble();
    mEchoLevel = ThisParameters["echo_level"].GetInt();
    mReferenceVariable = ThisParameters["reference_variable_name"].GetString();
}

template<SizeType TDim>
void SimpleErrorCalculatorProcess<TDim>::Execute()
{
    // Initialize the metric
    // a) Check for Metric Scalar in Meshing Application
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has("METRIC_SCALAR")) << "Import Meshing Application" <<std::endl;
    KRATOS_ERROR_IF_NOT(KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has("NODAL_TEMP_GRADIENT")) << "Error Creating Temperature Gradient" <<std::endl;
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has("NODAL_AREA")) << "ERROR:: NODAL_AREA Variable doesn't exist" <<std::endl;
    const double& scalar_variable = KratosComponents<Variable<double>>::Get("METRIC_SCALAR");

    // b) Retrive Nodes and Elements from the Model Part
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    KRATOS_DEBUG_ERROR_IF(nodes_array.size() == 0) <<  "ERROR:: Empty list of nodes" << std::endl;
    ElementsArrayType& elements_array = mrThisModelPart.Elements();
    KRATOS_DEBUG_ERROR_IF(elements_array.size() == 0) <<  "ERROR:: Empty list of elements" << std::endl;

    // c) Initialize Metric Scalar for all nodes to 0
    const auto it_node_begin = nodes_array.begin();
    if (it_node_begin->Has(scalar_variable) == false) {
        // Iterate over the nodes
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i)
            (it_node_begin + i)->SetValue(scalar_variable, 0.0);
    }

    // d) Call the Error Estimation Calculator
    ErrorEstimatorImplementation();

    // e) Call the MetricScalar Function 
    CalculateMetricScalar();
}

template<SizeType TDim>
void SimpleErrorCalculatorProcess<TDim>::ErrorEstimatorImplementation()
{
    // a) Obtain Nodes and Elements from Model Part
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    // b) Loop over Elements and calculate RHS
    const int number_nodes = static_cast<int>(nodes_array.size());
    const int number_elements = static_cast<int>(elements_array.size());
    KRATOS_DEBUG_ERROR_IF(number_nodes == 0) <<  "ERROR:: Empty list of nodes" << std::endl;
    KRATOS_DEBUG_ERROR_IF(number_elements == 0) <<  "ERROR:: Empty list of elements" << std::endl;

    const auto it_elem_begin = elements_array.begin();
 
    // Loop over the elements
    for (int i_elem = 0; i_elem < number_elements; ++i_elem) {
        auto it_elem = it_elem_begin + i_elem;
        auto r_geometry = it_elem->GetGeometry();
        auto n_nodes = r_geometry.size();

        Vector GaussWeights;
        Vector DetJ;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;

        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_1);
        const unsigned int NumGPoints = integration_points.size();

        r_geometry.ShapeFunctionsIntegrationPointsGradients(ShapeDerivatives, DetJ, GeometryData::GI_GAUSS_1);
        
        if (ShapeFunctions.size1() != NumGPoints || ShapeFunctions.size2() != n_nodes) {
            ShapeFunctions.resize(NumGPoints, n_nodes, false);
        }
        ShapeFunctions = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

        if (GaussWeights.size() != NumGPoints) {
            GaussWeights.resize(NumGPoints, false);
        }

        for (unsigned int g = 0; g < NumGPoints; g++)
            GaussWeights[g] = DetJ[g] * integration_points[g].Weight();
        
        for (unsigned int g = 0; g < NumGPoints; g++) {
            const auto& rDN_DX = ShapeDerivatives[g];
            const auto& Ncontainer = ShapeFunctions[g];

            std::vector<double,TDim> GaussPointTGrad;
            for (unsigned int k = 0; k < TDim; k++) {
                GaussPointTGrad[k] = 0.0;
            }

            for (unsigned int j = 0; j < TDim; j++) {
                for (unsigned int i_node = 0; i_node < n_nodes; i_node++) {
                    GaussPointTGrad[j] += rDN_DX(i_node,j)*r_geometry[i_node].FastGetSolutionStepValue(TEMPERATURE); 
                }
                GaussPointTGrad[j] *= GaussWeights[g];
            }

            for (int i_node = 0; i_node < n_nodes; i_node++) {
                for (unsigned int i_dim = 0; i_dim < TDim; i_dim++) {
                    r_geometry[i_node].FastGetSolutionStepValue(NODAL_TEMP_GRADIENT)[i_dim] += Ncontainer[i_node]*GaussPointTGrad[i_dim]*rgeometry[i_node].FastGetSolutionStepValue(NODAL_AREA);
                }
            } 
        }                

}
