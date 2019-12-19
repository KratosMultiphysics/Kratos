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
#include <iostream>
#include <fstream>

// External includes

// Project includes
#include "utilities/variable_utils.h"
#include "includes/element.h"
#include "processes/find_nodal_neighbours_process.h"

// Application includes
#include "convection_diffusion_application_variables.h"
#include "custom_utilities/meshing_utilities.h"
#include "meshing_application_variables.h"

// Include base h
#include "custom_processes/metrics_temperature_gradient_process.h"

namespace Kratos
{

template <std::size_t TDim>
MetricsTemperatureGradientProcess<TDim>::MetricsTemperatureGradientProcess(ModelPart &rThisModelPart, Parameters ThisParameters) : mrThisModelPart(rThisModelPart)
{
    Parameters default_parameters = Parameters(R"(
    {
        "minimal_size"                        : 0.01,
        "maximal_size"                        : 10.0,
        "nodal_averaging"                     : true,
        "reference_variable_name"             : "ERROR_RATIO",
        "historical_results"                  : true,
        "echo_level"                          : 0
    })");

    ThisParameters.ValidateAndAssignDefaults(default_parameters);
    mMinSize = ThisParameters["minimal_size"].GetDouble();
    mMaxSize = ThisParameters["maximal_size"].GetDouble();
    mEchoLevel = ThisParameters["echo_level"].GetInt();
    mNodalAveragingH = ThisParameters["nodal_averaging"].GetBool();
    mHistoricalResults = ThisParameters["historical_results"].GetBool();
}

template <std::size_t TDim>
void MetricsTemperatureGradientProcess<TDim>::Execute()
{
    KRATOS_TRY

    this->CreateMap();
    // Initialize the metric
    // a) Check for Metric Scalar in Meshing Application
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has("METRIC_SCALAR")) << "Import Meshing Application" << std::endl;
    const Variable<double> &scalar_variable = KratosComponents<Variable<double>>::Get("METRIC_SCALAR");

    // c) Initialize Metric Scalar and Nodal Temperature Gradient for all nodes to 0
    VariableUtils().SetNonHistoricalVariableToZero(scalar_variable, mrThisModelPart.Nodes());
    if (mHistoricalResults)
    {
        VariableUtils()
            .SetHistoricalVariableToZero(NODAL_TEMP_GRADIENT, mrThisModelPart.Nodes());
        VariableUtils()
            .SetHistoricalVariableToZero(NODAL_ERROR_PROJ, mrThisModelPart.Nodes());
    }
    else
    {
        VariableUtils()
            .SetNonHistoricalVariableToZero(NODAL_TEMP_GRADIENT, mrThisModelPart.Nodes());
        VariableUtils()
            .SetNonHistoricalVariableToZero(NODAL_ERROR_PROJ, mrThisModelPart.Nodes());
    }

    Vector nodal_area(mrThisModelPart.NumberOfNodes(), 0.0);
    this->CalculateNodalArea(nodal_area);

    // d) Call the Nodal Temperature Gradient Calculator Function
    this->CalculateNodalTempGradient(nodal_area);

    // // e) Call the Nodal Error Function
    this->CalculateNodalError(nodal_area);

    // e) Call the Metric Scalar Function
    this->CalculateMetric();

    KRATOS_CATCH("");
}

template <std::size_t TDim>
void MetricsTemperatureGradientProcess<TDim>::CalculateNodalTempGradient(Vector &nodal_area)
{
    //b) Loop over Elements and calculate RHS
    const int number_nodes = mrThisModelPart.NumberOfNodes();
    const int number_elements = mrThisModelPart.NumberOfElements();
    Matrix nodal_grad = ZeroMatrix(number_nodes, 3);

    // Loop over the elements
    for (unsigned int i_elem = 0; i_elem < number_elements; i_elem++)
    {
        ModelPart::ElementType &r_element = *(mrThisModelPart.ElementsBegin() + i_elem);
        ModelPart::ElementType::GeometryType &r_geometry = r_element.GetGeometry();
        const auto n_nodes = r_geometry.PointsNumber();

        Vector GaussWeights;
        Vector DetJ;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        unsigned int NumGPoints = 0;
        CalculateGeomData(r_geometry, ShapeFunctions, ShapeDerivatives, DetJ, GaussWeights, NumGPoints);
        for (unsigned int g = 0; g < NumGPoints; g++)
        {
            const Matrix &rDN_DX = ShapeDerivatives[g];
            const Vector &Ncontainer = row(ShapeFunctions, g);
            Vector GaussPointTGrad(3, 0.0);
            for (unsigned int j = 0; j < TDim; j++)
            {
                for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
                {
                    //const int n_id = r_geometry[i_node].Id();
                    auto nodal_temperature = r_geometry[i_node].FastGetSolutionStepValue(TEMPERATURE);
                    GaussPointTGrad[j] += rDN_DX(i_node, j) * nodal_temperature;
                }
                GaussPointTGrad[j] *= GaussWeights[g];
            }
            for (int i_node = 0; i_node < n_nodes; i_node++)
            {
                const int n_id = r_geometry[i_node].Id();
                auto &rNodalArea = nodal_area[mNodeMap.find(n_id)->second];
                if (rNodalArea == 1.0)
                {
                    nodal_grad(mNodeMap.find(n_id)->second, 0) += 0.0;
                    nodal_grad(mNodeMap.find(n_id)->second, 1) += 0.0;
                    nodal_grad(mNodeMap.find(n_id)->second, 2) += 0.0;
                }
                else
                {
                    nodal_grad(mNodeMap.find(n_id)->second, 0) += Ncontainer[i_node] * GaussPointTGrad[0] / rNodalArea;
                    nodal_grad(mNodeMap.find(n_id)->second, 1) += Ncontainer[i_node] * GaussPointTGrad[1] / rNodalArea;
                    nodal_grad(mNodeMap.find(n_id)->second, 2) += Ncontainer[i_node] * GaussPointTGrad[2] / rNodalArea;
                }
            }
        }
    }

    //#pragma omp parallel for
    for (int i = 0; i < number_nodes; i++)
    {
        ModelPart::NodeType &r_node = *(mrThisModelPart.NodesBegin() + i);
        const Vector &nodal_grad_row = row(nodal_grad, mNodeMap.find(r_node.Id())->second);

        if (mHistoricalResults)
        {
            r_node.FastGetSolutionStepValue(NODAL_TEMP_GRADIENT) = nodal_grad_row;
        }
        else
        {
            r_node.SetValue(NODAL_TEMP_GRADIENT, nodal_grad_row);
        }
    }
}

template <std::size_t TDim>
void MetricsTemperatureGradientProcess<TDim>::CalculateNodalError(Vector &nodal_area)
{
    // a) Obtain Nodes and Elements from Model Part
    const int number_nodes = mrThisModelPart.NumberOfNodes();
    const int number_elements = mrThisModelPart.NumberOfElements();
    
    // KRATOS_INFO(this->Info()) << "Calculating nodal error for " << number_nodes << " nodes in " << mrThisModelPart.Name() << "\n";

    double global_gw = 0.0;
    double global_del_sigma = 0.0;

    // Debug Options
    // double max_Nodal_error = 0.0;
    // double max_N_tempGrad = 0.0;
    // double max_GP_tempGrad = 0.0;
    // double max_elem_del_sigma = 0.0;
    double max_mElementError = 0.0;

    Vector element_del_sigma = ZeroVector(number_elements);
    Vector element_error_norm_squared = ZeroVector(number_elements);
    Vector element_gw = ZeroVector(number_elements);
    Vector nodal_error_proj = ZeroVector(number_nodes);
    Matrix element_error = ZeroMatrix(number_elements, 3);

    // c) Loop over the elements to calculate Element Error
    for (int i_elem = 0; i_elem < number_elements; i_elem++)
    {
        ModelPart::ElementType &r_element = *(mrThisModelPart.ElementsBegin() + i_elem);
        ModelPart::ElementType::GeometryType &r_geometry = r_element.GetGeometry();
        const auto n_nodes = r_geometry.PointsNumber();

        Vector GaussWeights;
        Vector DetJ;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        unsigned int NumGPoints = 0;
        CalculateGeomData(r_geometry, ShapeFunctions, ShapeDerivatives, DetJ, GaussWeights, NumGPoints);

        for (unsigned int g = 0; g < NumGPoints; g++)
        {
            element_gw[i_elem] += GaussWeights[g];
            global_gw += GaussWeights[g];
            const Matrix &rDN_DX = ShapeDerivatives[g];
            const Vector &Ncontainer = row(ShapeFunctions, g);

            Vector GaussPointTGrad(3, 0.0);
            Vector NodalTempGrad(3, 0.0);

            for (unsigned int j = 0; j < TDim; j++)
            {
                for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
                {
                    GaussPointTGrad[j] += rDN_DX(i_node, j) * r_geometry[i_node].FastGetSolutionStepValue(TEMPERATURE);
                }
            }

            for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
            {
                if (mHistoricalResults)
                {
                    NodalTempGrad[0] += Ncontainer[i_node] * r_geometry[i_node].FastGetSolutionStepValue(NODAL_TEMP_GRADIENT_X);
                    NodalTempGrad[1] += Ncontainer[i_node] * r_geometry[i_node].FastGetSolutionStepValue(NODAL_TEMP_GRADIENT_Y);
                    if (TDim == 3)
                    {
                        NodalTempGrad[2] += Ncontainer[i_node] * r_geometry[i_node].FastGetSolutionStepValue(NODAL_TEMP_GRADIENT_Z);
                    }
                }
                else
                {
                    NodalTempGrad[0] += Ncontainer[i_node] * r_geometry[i_node].GetValue(NODAL_TEMP_GRADIENT_X);
                    NodalTempGrad[1] += Ncontainer[i_node] * r_geometry[i_node].GetValue(NODAL_TEMP_GRADIENT_Y);
                    if (TDim == 3)
                    {
                        NodalTempGrad[2] += Ncontainer[i_node] * r_geometry[i_node].GetValue(NODAL_TEMP_GRADIENT_Z);
                    }
                }
            }

            for (unsigned int j = 0; j < TDim; j++)
            {
                NodalTempGrad[j] *= GaussWeights[g];
                GaussPointTGrad[j] *= GaussWeights[g];                
                element_error(i_elem, j) += (NodalTempGrad[j] - GaussPointTGrad[j]);
            }
        }

        for (unsigned int j = 0; j < TDim; j++)
        {
            element_error_norm_squared[i_elem] += element_error(i_elem, j) * element_error(i_elem, j);
        }
        
        element_del_sigma[i_elem] = std::pow(element_error_norm_squared[i_elem] / element_gw[i_elem], 0.5);
        global_del_sigma += element_error_norm_squared[i_elem];
        
    }
    global_del_sigma = std::pow(global_del_sigma / global_gw, 0.5);

    if (mElementError.size() != number_elements)
    {
        mElementError.resize(number_elements, false);
        mElementError.clear();
    }

    for (int i_elem = 0; i_elem < number_elements; i_elem++)
    {
        mElementError[i_elem] = element_del_sigma[i_elem]/global_del_sigma;

        if (mElementError[i_elem] > max_mElementError)
        {
            max_mElementError = mElementError[i_elem];
        }
        
    }
    
    KRATOS_WATCH(max_mElementError);

    // d) Loop over the elements to calculate Nodal Error Projections
    // for (int i_elem = 0; i_elem < number_elements; i_elem++)
    // {
    //     ModelPart::ElementType &r_element = *(mrThisModelPart.ElementsBegin() + i_elem);
    //     ModelPart::ElementType::GeometryType &r_geometry = r_element.GetGeometry();
    //     const auto n_nodes = r_geometry.PointsNumber();

    //     Vector GaussWeights;
    //     Vector DetJ;
    //     Matrix ShapeFunctions;
    //     ShapeFunctionDerivativesArrayType ShapeDerivatives;
    //     unsigned int NumGPoints = 0;
    //     CalculateGeomData(r_geometry, ShapeFunctions, ShapeDerivatives, DetJ, GaussWeights, NumGPoints);

    //     const Vector &Ncontainer = row(ShapeFunctions, 0);

    //     for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
    //     {
    //         const int n_id = r_geometry[i_node].Id();
    //         auto& rNodalArea = nodal_area[mNodeMap.find(n_id)->second];
    //         if ( rNodalArea != 0.0 || rNodalArea != 1.0)
    //         {
    //             if (mHistoricalResults)
    //             {
    //                 r_geometry[i_node].FastGetSolutionStepValue(NODAL_ERROR_PROJ) += Ncontainer[i_node] * element_del_sigma[i_elem] / (global_del_sigma * rNodalArea);
    //             }
    //             else
    //             {
    //                 r_geometry[i_node].GetValue(NODAL_ERROR_PROJ) += Ncontainer[i_node] * element_del_sigma[i_elem] / (global_del_sigma * rNodalArea);
    //             }
    //         }
    //     }
    // }
}

template <std::size_t TDim>
void MetricsTemperatureGradientProcess<TDim>::CalculateMetric()
{
    const double tolerance = std::numeric_limits<double>::epsilon();

    ElementsArrayType& elements_array = mrThisModelPart.Elements();
    const int number_elements = mrThisModelPart.NumberOfElements();

    // Debug Options
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < number_elements; i_elem++)
    {
        auto it_elem = elements_array.begin() + i_elem;
        
        // Compute the current element size h
        MeshingUtilities::ComputeElementSize(it_elem);

        const double element_error = mElementError[i_elem];
        const double coeff = std::abs(element_error) < tolerance ? 1.0 : 1.0/element_error;
        double new_element_size = coeff*it_elem->GetValue(ELEMENT_H);

        // Check if element sizes are in specified limits. If not, set them to the limit case
        if(new_element_size < mMinSize)
            new_element_size = mMinSize;

        if(new_element_size > mMaxSize)
            new_element_size = mMaxSize;

        it_elem->SetValue(ELEMENT_H, new_element_size);
    }

    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    const int number_nodes = mrThisModelPart.NumberOfNodes();

    {
        FindNodalNeighboursProcess find_neighbours(mrThisModelPart);
        if (nodes_array.begin()->Has(NEIGHBOUR_ELEMENTS)) 
        {
            find_neighbours.ClearNeighbours();
        }
        find_neighbours.Execute();
    }

    #pragma omp parallel for
    for (int i_node = 0; i_node < number_nodes; i_node++)
    {
        auto it_node = nodes_array.begin() + i_node;
        /**************************************************************************
        ** Determine nodal element size h:
        ** if mAverageNodalH == true : the nodal element size is averaged from the element size of neighboring elements
        ** if mAverageNodalH == false: the nodal element size is the minimum element size from neighboring elements
        */
        double h_min = 0.0;
        auto& neigh_elements = it_node->GetValue(NEIGHBOUR_ELEMENTS);
        for(auto i_neighbour_elements = neigh_elements.begin(); i_neighbour_elements != neigh_elements.end(); i_neighbour_elements++){
            const double element_h = i_neighbour_elements->GetValue(ELEMENT_H);
            if(mNodalAveragingH == false) {
                if(h_min == 0.0 || h_min > element_h)
                    h_min = element_h;
            } else {
                h_min += element_h;
            }
        }

        // Average Nodal H
        if(mNodalAveragingH) h_min = h_min/static_cast<double>(neigh_elements.size());

        // Setting value
        it_node->SetValue(METRIC_SCALAR, h_min);

        //KRATOS_INFO_IF("MetricErrorProcess", mEchoLevel > 2) << "Node " << it_node->Id() << " has metric: "<< h_min << std::endl;
        KRATOS_INFO("MetricTemperatureGradientProcess") << "Node " << it_node->Id() << " has metric: "<< h_min << std::endl;
    }
    
}

template <std::size_t TDim>
void MetricsTemperatureGradientProcess<TDim>::CreateMap()
{
    int node_index = 0;
    int elem_index = 0;

    std::fstream OutFile("Node_Map.txt", std::ios::out);
    OutFile << "node_index"
            << "\t\t"
            << "r_node.Id()"
            << "\n";

    for (const ModelPart::NodeType &r_node : mrThisModelPart.Nodes())
    {
        OutFile << node_index << "\t\t" << r_node.Id() << "\n";
        mNodeMap.insert(std::make_pair<int, int>(r_node.Id(), node_index++));
    }

    for (const ModelPart::ElementType &r_elem : mrThisModelPart.Elements())
    {
        mElemMap.insert(std::make_pair<int, int>(r_elem.Id(), elem_index++));
    }
    OutFile.close();
    KRATOS_INFO(this->Info()) << "Creating Node Map for " << mrThisModelPart.NumberOfNodes() << " nodes in " << mrThisModelPart.Name() << "\n";
}

template <std::size_t TDim>
void MetricsTemperatureGradientProcess<TDim>::CalculateNodalArea(Vector &rNodalArea) const
{
    std::fstream OutFile("Nodal_Area.txt", std::ios::out);
    OutFile << "r_node.Id()"
            << "\t\t"
            << "Nodal_Area"
            << "\n";
    int count_zeros = 0;
    int count_zero_node_elem = 0;

    const int number_nodes = mrThisModelPart.NumberOfNodes();
    const int number_elements = mrThisModelPart.NumberOfElements();
    KRATOS_INFO(this->Info()) << "Calculating nodal area of " << number_nodes << " nodes in " << mrThisModelPart.Name() << "\n";

    if (rNodalArea.size() != number_nodes)
    {
        rNodalArea.resize(number_nodes);
        rNodalArea.clear();
    }

    // b) Loop over Elements and calculate nodal area
    for (int i_elem = 0; i_elem < number_elements; ++i_elem)
    {
        const ModelPart::ElementType &r_element = *(mrThisModelPart.ElementsBegin() + i_elem);
        const ModelPart::ElementType::GeometryType &r_geometry = r_element.GetGeometry();
        const int n_nodes = r_geometry.PointsNumber();

        for (int i_node = 0; i_node < n_nodes; i_node++)
        {
            const int n_id = r_geometry[i_node].Id();
            if (n_nodes > 0)
            {
                rNodalArea[mNodeMap.find(n_id)->second] += r_geometry.Area() / static_cast<double>(n_nodes);
            }
            else
            {
                count_zero_node_elem++;
            }
        }
    }

    for (int i_node = 0; i_node < number_nodes; i_node++)
    {
        const ModelPart::NodeType &r_node = *(mrThisModelPart.NodesBegin() + i_node);
        if (rNodalArea[i_node] == 0.0)
        {
            rNodalArea[i_node] = 1.0;
            count_zeros++;
        }
        OutFile << r_node.Id() << "\t\t" << rNodalArea[i_node] << "\n";
    }
    OutFile << "Number of 0s"
            << "\t\t" << count_zeros << "\n";
    OutFile << "Number of 0 node elements"
            << "\t\t" << count_zero_node_elem << "\n";
    OutFile.close();
}

template <std::size_t TDim>
void MetricsTemperatureGradientProcess<TDim>::CalculateGeomData(GeometryType &r_geom, Matrix &ShapeFunctions, ShapeFunctionDerivativesArrayType &ShapeDerivatives, Vector &DetJ, Vector &GaussWeights, unsigned int &NumGPoints)
{
    const GeometryType::IntegrationPointsArrayType &integration_points = r_geom.IntegrationPoints(GeometryData::GI_GAUSS_1);
    NumGPoints = integration_points.size();
    const auto n_nodes = r_geom.size();

    r_geom.ShapeFunctionsIntegrationPointsGradients(ShapeDerivatives, DetJ, GeometryData::GI_GAUSS_1);

    if (ShapeFunctions.size1() != NumGPoints || ShapeFunctions.size2() != n_nodes)
    {
        ShapeFunctions.resize(NumGPoints, n_nodes, false);
    }
    ShapeFunctions = r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

    if (GaussWeights.size() != NumGPoints)
    {
        GaussWeights.resize(NumGPoints, false);
    }

    for (unsigned int g = 0; g < NumGPoints; g++)
    {
        GaussWeights[g] = DetJ[g] * integration_points[g].Weight();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class MetricsTemperatureGradientProcess<2>;
template class MetricsTemperatureGradientProcess<3>;

}; // namespace Kratos