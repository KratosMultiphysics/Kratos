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
#include "convection_diffusion_application_variables.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
    
template<SizeType TDim>
SimpleErrorCalculatorProcess<TDim>::SimpleErrorCalculatorProcess(ModelPart& rThisModelPart, Parameters ThisParameters):mrThisModelPart(rThisModelPart)
{
    //WIP
	mrThisModelPart.AddNodalSolutionStepVariable(NODAL_TEMP_GRADIENT);
	mrThisModelPart.AddNodalSolutionStepVariable(NODAL_ERROR_PROJ);
            
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
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has("NODAL_AREA")) << "ERROR:: NODAL_AREA Variable doesn't exist" <<std::endl;
    //KRATOS_ERROR_IF_NOT(KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has("NODAL_TEMP_GRADIENT")) << "ERROR:: NODAL_TEMP_GRADIENT Does not Exist" <<std::endl;
    //KRATOS_ERROR_IF_NOT(KratosComponents<Variable<array_1d<double, 3>>>::Has("NODAL_TEMP_GRADIENT")) << "ERROR:: NODAL_TEMP_GRADIENT Does not Exist" <<std::endl;
    
    const Variable<double>& scalar_variable = KratosComponents<Variable<double>>::Get("METRIC_SCALAR");
    const Variable<array_1d<double, 3>>& TDim_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get("NODAL_TEMP_GRADIENT");

    // b) Retrive Nodes and Elements from the Model Part
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    ElementsArrayType& elements_array = mrThisModelPart.Elements();
    KRATOS_DEBUG_ERROR_IF(nodes_array.size() == 0) <<  "ERROR:: Empty list of nodes" << std::endl;
    KRATOS_DEBUG_ERROR_IF(elements_array.size() == 0) <<  "ERROR:: Empty list of elements" << std::endl;

    // c) Initialize Metric Scalar and Nodal Temperature Gradient for all nodes to 0
    const auto it_node_begin = nodes_array.begin();
    if (it_node_begin->Has(scalar_variable) == false) {
        // Iterate over the nodes
        //#pragma omp parallel for
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i)
            (it_node_begin + i)->SetValue(scalar_variable, 0.0);
    }

    if (it_node_begin->Has(TDim_variable) == false) {
        // Iterate over the nodes
        //#pragma omp parallel for
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            (it_node_begin + i)->SetValue(TDim_variable, ZeroVector(3));
        }
    }

    Vector nodal_area = ZeroVector(nodes_array.size());
    CalculateNodalArea(nodal_area);

    // d) Call the Nodal Temperature Gradient Calculator Function
    CalculateNodalTempGradient(nodal_area);

    // e) Call the Nodal Error Function
    CalculateNodalError(nodal_area);

    // e) Call the Metric Scalar Function 
    //CalculateMetricScalar();
}

template<SizeType TDim>
void SimpleErrorCalculatorProcess<TDim>::CalculateNodalTempGradient(Vector& nodal_area)
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
    const auto it_node_begin = nodes_array.begin();
    Matrix nodal_grad = ZeroMatrix(number_nodes, 3);
    //std::vector<array_1d<double,3>> nodal_grad(number_nodes);

    //for (unsigned int i_grad = 0; i_grad < number_nodes; i_grad++)
    //{
    //    /* code */
    //}
    
    const Variable<array_1d<double, 3>>& nodal_temp_grad = KratosComponents<Variable<array_1d<double, 3>>>::Get("NODAL_TEMP_GRADIENT");
 
    // Loop over the elements
    for (unsigned int i_elem = 0; i_elem < number_elements; i_elem++) {
        auto it_elem = it_elem_begin + i_elem;
        auto r_geometry = it_elem->GetGeometry();
        const auto n_nodes = r_geometry.size();

        Vector GaussWeights;
        Vector DetJ;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        unsigned int NumGPoints = 0;
        CalculateGeomData(r_geometry, ShapeFunctions, ShapeDerivatives, DetJ, GaussWeights, NumGPoints);
        
        for (unsigned int g = 0; g < NumGPoints; g++) {
            const auto& rDN_DX = ShapeDerivatives[g];
            const Vector& Ncontainer = row(ShapeFunctions, g);

            Vector GaussPointTGrad(3);
            for (unsigned int k = 0; k < 3; k++) {
                GaussPointTGrad[k] = 0.0;
            }

            for (unsigned int j = 0; j < TDim; j++) {
                for (unsigned int i_node = 0; i_node < n_nodes; i_node++) {
                    GaussPointTGrad[j] += rDN_DX(i_node,j)*r_geometry[i_node].FastGetSolutionStepValue(TEMPERATURE); 
                }
                GaussPointTGrad[j] *= GaussWeights[g];
            }

            for (int i_node = 0; i_node < n_nodes; i_node++) {
                for (int j = 0; j < TDim; j++) {
                    nodal_grad(r_geometry[i_node].Id()-1, j) += Ncontainer[i_node]*GaussPointTGrad[j]/nodal_area[r_geometry[i_node].Id()-1];
                }
            } 
        }                
    }

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); i++) {
        auto it_node = it_node_begin + i;
		const Vector& nodal_grad_row = row(nodal_grad, i);
        it_node->SetValue(nodal_temp_grad, nodal_grad_row);
        //array_1d<double>& nodal_temp_grad = (it_node_begin + i)->FastGetSolutionStepValue(NODAL_TEMP_GRADIENT);
        //nodal_temp_grad = nodal_grad[i];
        //(it_node_begin + i)->SetValue(NODAL_TEMP_GRADIENT_X, nodal_grad(i,0));
        //(it_node_begin + i)->SetValue(NODAL_TEMP_GRADIENT_Y, nodal_grad(i,1));
        //(it_node_begin + i)->SetValue(NODAL_TEMP_GRADIENT_Z, nodal_grad(i,2));
    }

}

template<SizeType TDim>
void SimpleErrorCalculatorProcess<TDim>::CalculateNodalError(Vector& nodal_area)
{
    // a) Obtain Nodes and Elements from Model Part
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    // b) Check Nodes and Elements
    const int number_nodes = static_cast<int>(nodes_array.size());
    const int number_elements = static_cast<int>(elements_array.size());
    KRATOS_DEBUG_ERROR_IF(number_nodes == 0) <<  "ERROR:: Empty list of nodes" << std::endl;
    KRATOS_DEBUG_ERROR_IF(number_elements == 0) <<  "ERROR:: Empty list of elements" << std::endl;

    const auto it_elem_begin = elements_array.begin();
    const auto it_node_begin = nodes_array.begin();
    double global_gw = 0.0;
    double global_del_sigma = 0.0;
    Vector element_del_sigma = ZeroVector(number_elements);
    Vector element_error_norm_squared = ZeroVector(number_elements);
    Vector element_gw = ZeroVector(number_elements);
    Vector nodal_error_proj = ZeroVector(number_nodes);
    Matrix element_error = ZeroMatrix(number_elements, 3);

    // c) Loop over the elements to calculate Element Error
    for (int i_elem = 0; i_elem < number_elements; i_elem++) {
        auto it_elem = it_elem_begin + i_elem;
        auto r_geometry = it_elem->GetGeometry();
        const auto n_nodes = r_geometry.size();
        unsigned int elem_index = i_elem;

        Vector GaussWeights;
        Vector DetJ;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        unsigned int NumGPoints = 0;
        CalculateGeomData(r_geometry, ShapeFunctions, ShapeDerivatives, DetJ, GaussWeights, NumGPoints);
        
        for (unsigned int g = 0; g < NumGPoints; g++) {
            element_gw[elem_index] += GaussWeights[g];
            global_gw += GaussWeights[g];
            const auto& rDN_DX = ShapeDerivatives[g];
            const Vector& Ncontainer = row(ShapeFunctions, g);

            Vector GaussPointTGrad(3);
            Vector NodalTempGrad(3);
            for (unsigned int k = 0; k < 3; k++) {
                GaussPointTGrad[k] = 0.0;
                NodalTempGrad[k] = 0.0;
            }

            for (unsigned int j = 0; j < 3; j++) {
                for (unsigned int i_node = 0; i_node < n_nodes; i_node++) {
                    GaussPointTGrad[j] += rDN_DX(i_node,j)*r_geometry[i_node].FastGetSolutionStepValue(TEMPERATURE); 
                }
                GaussPointTGrad[j] *= GaussWeights[g];
            }

            for (unsigned int i_node = 0; i_node < n_nodes; i_node++) {
                if (TDim == 2) {
                    //array_1d<double, 3> tempGrad = r_geometry[i_node].FastGetSolutionStepValue(NODAL_TEMP_GRADIENT);
                    //NodalTempGrad[0] += Ncontainer[i_node]*tempGrad[0];
                    //NodalTempGrad[1] += Ncontainer[i_node]*tempGrad[1];
                    NodalTempGrad[0] += Ncontainer[i_node]*r_geometry[i_node].FastGetSolutionStepValue(NODAL_TEMP_GRADIENT_X);
                    NodalTempGrad[1] += Ncontainer[i_node]*r_geometry[i_node].FastGetSolutionStepValue(NODAL_TEMP_GRADIENT_Y);
                }
                else if (TDim == 3) {
                    //array_1d<double, 3> tempGrad = r_geometry[i_node].FastGetSolutionStepValue(NODAL_TEMP_GRADIENT);
                    //NodalTempGrad[0] += Ncontainer[i_node]*tempGrad[0];
                    //NodalTempGrad[1] += Ncontainer[i_node]*tempGrad[1];
                    //NodalTempGrad[2] += Ncontainer[i_node]*tempGrad[2];
                    NodalTempGrad[0] += Ncontainer[i_node]*r_geometry[i_node].FastGetSolutionStepValue(NODAL_TEMP_GRADIENT_X);
                    NodalTempGrad[1] += Ncontainer[i_node]*r_geometry[i_node].FastGetSolutionStepValue(NODAL_TEMP_GRADIENT_Y);
                    NodalTempGrad[2] += Ncontainer[i_node]*r_geometry[i_node].FastGetSolutionStepValue(NODAL_TEMP_GRADIENT_Z);
                }
            }

            for (unsigned int j = 0; j < TDim; j++) {
                NodalTempGrad[j] *= GaussWeights[g];
                element_error(elem_index, j) += (NodalTempGrad[j] - GaussPointTGrad[j]);               
            }
        }

        for (unsigned int j = 0; j < TDim; j++) {
            element_error_norm_squared[elem_index] += element_error(elem_index, j)*element_error(elem_index, j);
        }
        element_del_sigma[elem_index] = element_error_norm_squared[elem_index]/element_gw[elem_index];
        global_del_sigma += element_error_norm_squared[elem_index];
    }
    global_del_sigma = global_del_sigma/global_gw;

    // d) Loop over the elements to calculate Nodal Error Projections
    for (int i_elem = 0; i_elem < number_elements; i_elem++) {
        auto it_elem = it_elem_begin + i_elem;
        auto r_geometry = it_elem->GetGeometry();
        const auto n_nodes = r_geometry.size();
        const auto elem_index = i_elem;

        Vector GaussWeights;
        Vector DetJ;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
		unsigned int NumGPoints = 0;
        CalculateGeomData(r_geometry, ShapeFunctions, ShapeDerivatives, DetJ, GaussWeights, NumGPoints);

        const Vector& Ncontainer = row(ShapeFunctions, 0);

        for (unsigned int i_node = 0; i_node < n_nodes; i_node++) {
            r_geometry[i_node].FastGetSolutionStepValue(NODAL_ERROR_PROJ) += Ncontainer[i_node]*element_del_sigma[elem_index]/(global_del_sigma*nodal_area[r_geometry[i_node].Id()-1]);
        }        
    }
}

template<SizeType TDim>
void SimpleErrorCalculatorProcess<TDim>::CalculateNodalArea(Vector& nodal_area)
{
    // a) Obtain Nodes and Elements from Model Part
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    // b) Loop over Elements and calculate nodal area
    const int number_nodes = static_cast<int>(nodes_array.size());
    const int number_elements = static_cast<int>(elements_array.size());

    const auto it_elem_begin = elements_array.begin();
    for (unsigned int i_elem = 0; i_elem < number_elements; i_elem++) {
        auto it_elem = it_elem_begin + i_elem;
        auto r_geometry = it_elem->GetGeometry();
        const auto n_nodes = r_geometry.size();

        for (unsigned int i_node = 0; i_node < n_nodes; i_node++) {
            nodal_area[r_geometry[i_node].Id()-1] += + r_geometry.Area()/r_geometry.size();
        }
    }

}

template<SizeType TDim>
void SimpleErrorCalculatorProcess<TDim>::CalculateGeomData(GeometryType& r_geom, Matrix& ShapeFunctions, ShapeFunctionDerivativesArrayType& ShapeDerivatives,Vector& DetJ, Vector& GaussWeights, unsigned int& NumGPoints)
{
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(GeometryData::GI_GAUSS_1);
    NumGPoints = integration_points.size();
    const auto n_nodes = r_geom.size();

    r_geom.ShapeFunctionsIntegrationPointsGradients(ShapeDerivatives, DetJ, GeometryData::GI_GAUSS_1);
        
    if (ShapeFunctions.size1() != NumGPoints || ShapeFunctions.size2() != n_nodes) {
        ShapeFunctions.resize(NumGPoints, n_nodes, false);
    }
    ShapeFunctions = r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

    if (GaussWeights.size() != NumGPoints) {
        GaussWeights.resize(NumGPoints, false);
    }

    for (unsigned int g = 0; g < NumGPoints; g++) {
        GaussWeights[g] = DetJ[g] * integration_points[g].Weight();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class SimpleErrorCalculatorProcess<2>;
template class SimpleErrorCalculatorProcess<3>;

}; // namespace Kratos