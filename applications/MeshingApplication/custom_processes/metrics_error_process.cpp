// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// Project includes
#include "custom_processes/metrics_error_process.h"

namespace Kratos
{
template<SizeType TDim>
MetricErrorProcess<TDim>::MetricErrorProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters
        ):mrThisModelPart(rThisModelPart)
{               
    Parameters default_parameters = Parameters(R"(
    {
        "minimal_size"                        : 0.01,
        "maximal_size"                        : 1.0,
        "error_strategy_parameters":
        {
            "target_error"                        : 0.01,
            "set_number_of_elements"              : false,
            "number_of_elements"                  : 1000,
            "average_nodal_h"                     : false
        },
        "echo_level"                          : 0
    })"
    );

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mMinSize = ThisParameters["minimal_size"].GetDouble();
    mMaxSize = ThisParameters["maximal_size"].GetDouble();

    mSetElementNumber = ThisParameters["error_strategy_parameters"]["set_number_of_elements"].GetBool();
    mElementNumber = ThisParameters["error_strategy_parameters"]["number_of_elements"].GetInt();
    mTargetError = ThisParameters["error_strategy_parameters"]["error"].GetDouble();
    mAverageNodalH = ThisParameters["error_strategy_parameters"]["average_nodal_h"].GetBool();

    mEchoLevel = ThisParameters["echo_level"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void MetricErrorProcess<TDim>::Execute()
{
    /******************************************************************************
    --1-- Calculate new element size --1--
    ******************************************************************************/
    CalculateElementSize();

    /******************************************************************************
    --2-- Calculate metric (for each node) --2--
    ******************************************************************************/
    CalculateMetric();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void MetricErrorProcess<TDim>::CalculateElementSize()
{
    // Get error values
    const ProcessInfo& r_process_info = mrThisModelPart.GetProcessInfo();
    const double energy_norm_overall = r_process_info[ENERGY_NORM_OVERALL];
    const double error_overall = r_process_info[ERROR_OVERALL];

    // Some initializations
    const double tolerance = std::numeric_limits<double>::epsilon();

     // Loop over all elements:
    ElementsArrayType& elements_array = mrThisModelPart.Elements();
    const int num_elem = static_cast<int>(elements_array.size());

    // Compute new element size
    #pragma omp parallel for
    for(int i_elem = 0; i_elem < num_elem; ++i_elem){
        auto it_elem = elements_array.begin() + i_elem;

        //Compute the current element size h
        ComputeElementSize(it_elem);

        // Compute new element size
        const double element_error = it_elem->GetValue(ELEMENT_ERROR);
        const double coeff = std::abs(element_error) < tolerance ? 1.0 : 1.0/element_error;
        double new_element_size = coeff * it_elem->GetValue(ELEMENT_H);

        // If a target number for elements is given: use this, else: use current element number
        // if(mSetElementNumber == true && mElementNumber<mrThisModelPart.Elements().size())
        if(mSetElementNumber == true)
            new_element_size *= std::sqrt((std::pow(energy_norm_overall, 2) + std::pow(error_overall, 2))/mElementNumber) * mTargetError;
        else
            new_element_size *= std::sqrt((energy_norm_overall*energy_norm_overall+error_overall*error_overall)/mrThisModelPart.Elements().size())*mTargetError;


        // Check if element sizes are in specified limits. If not, set them to the limit case
        if(new_element_size < mMinSize)
            new_element_size = mMinSize;
        if(new_element_size > mMaxSize)
            new_element_size = mMaxSize;

        it_elem->SetValue(ELEMENT_H, new_element_size);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void MetricErrorProcess<TDim>::CalculateMetric()
{
    // Getting metric variable
    const Variable<Vector>& metric_variable = KratosComponents<Variable<Vector>>::Get("MMG_METRIC");

    // Iteration over all nodes
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());

    #pragma omp parallel for
    for(int i_node = 0; i_node < num_nodes; ++i_node) {
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
            if(mAverageNodalH == false) {
                if(h_min == 0.0 || h_min > element_h)
                    h_min = element_h;
            } else {
                h_min += element_h;
            }
        }

        if(mAverageNodalH == true)
            h_min = h_min/static_cast<double>(neigh_elements.size());

        // Set metric
        BoundedMatrix<double, TDim, TDim> metric_matrix = ZeroMatrix(TDim, TDim);
        for(IndexType i = 0;i < TDim; ++i)
            metric_matrix(i,i) = 1.0/std::pow(h_min, 2);

        // Transform metric matrix to a vector
        Vector metric(3 * (TDim - 1));

        metric[0] = metric_matrix(0, 0);
        metric[1] = metric_matrix(0, 1);

        if (TDim == 2) {
            metric[2] = metric_matrix(1, 1);
        } else  {
            metric[2] = metric_matrix(0, 2);
            metric[3] = metric_matrix(1, 1);
            metric[4] = metric_matrix(1, 2);
            metric[5] = metric_matrix(2, 2);
        }

        // Setting value
        it_node->SetValue(metric_variable, metric);

        KRATOS_INFO_IF("MetricErrorProcess", mEchoLevel > 2) << "Node " << it_node->Id() << " has metric: "<< metric << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void MetricErrorProcess<TDim>::ComputeElementSize(ElementItType itElement)
{
    auto& this_geometry = itElement->GetGeometry();

    // Here we compute the element size. This process is designed for triangles and tetrahedra, so we only specify for this geometries. Otherwise we take the length (and we throw a warning)
    if (this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3){ // Triangular elements
        itElement->SetValue(ELEMENT_H, 2.0 * this_geometry.Circumradius());
    } else if(this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4){ // Tetrahedral elements
        itElement->SetValue(ELEMENT_H,std::pow(12.0 * this_geometry.Volume()/std::sqrt(2.0), 1.0/3.0));
    } else { // In any othe case just considers the length of the element
        KRATOS_WARNING("MetricErrorProcess") << "This process is designed for tetrahedra (3D) and triangles (2D). Error expected" << std::endl;
        itElement->SetValue(ELEMENT_H, this_geometry.Length());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class MetricErrorProcess<2>;
template class MetricErrorProcess<3>;

};// namespace Kratos.
