// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix
//                   Anna Rehr
//

// System includes

// External includes

// Project includes
#include "meshing_application_variables.h"
#include "processes/find_nodal_neighbours_process.h"
#include "custom_processes/metrics_error_process.h"
#include "custom_utilities/meshing_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
template<SizeType TDim>
MetricErrorProcess<TDim>::MetricErrorProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters
        ):mrThisModelPart(rThisModelPart)
{
    ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mMinSize = ThisParameters["minimal_size"].GetDouble();
    mMaxSize = ThisParameters["maximal_size"].GetDouble();

    mSetElementNumber = ThisParameters["error_strategy_parameters"]["set_target_number_of_elements"].GetBool();
    mElementNumber = ThisParameters["error_strategy_parameters"]["target_number_of_elements"].GetInt();
    mTargetError = ThisParameters["error_strategy_parameters"]["target_error"].GetDouble();
    mAverageNodalH = ThisParameters["error_strategy_parameters"]["perform_nodal_h_averaging"].GetBool();

    mEchoLevel = ThisParameters["echo_level"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void MetricErrorProcess<TDim>::Execute()
{
    /******************************************************************************
    --1-- Initialize metric --1--
    ******************************************************************************/
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();
    KRATOS_DEBUG_ERROR_IF(r_nodes_array.size() == 0) <<  "ERROR:: Empty list of nodes" << std::endl;
    if (!r_nodes_array.begin()->Has(METRIC_SCALAR)) {
        VariableUtils().SetNonHistoricalVariableToZero(METRIC_SCALAR, r_nodes_array);
    }

    /******************************************************************************
    --2-- Calculate new element size --2--
    ******************************************************************************/
    CalculateElementSize();

    /******************************************************************************
    --3-- Calculate metric (for each node) --3--
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
    ElementsArrayType& r_elements_array = mrThisModelPart.Elements();

    // Auxiliar member variables
    const auto number_of_elements = mrThisModelPart.NumberOfElements();

    // Compute new element size
    block_for_each(r_elements_array,
        [this,&tolerance,&energy_norm_overall,&error_overall,&number_of_elements](Element& rElement) {

        //Compute the current element size h
        MeshingUtilities::ComputeElementSize(rElement);

        // Compute new element size
        const double element_error = rElement.GetValue(ELEMENT_ERROR);
        const double coeff = std::abs(element_error) < tolerance ? 1.0 : 1.0/element_error;
        double new_element_size = coeff * rElement.GetValue(ELEMENT_H);

        // If a target number for elements is given: use this, else: use current element number
        // if(mSetElementNumber && mElementNumber < number_of_elements)
        if(mSetElementNumber)
            new_element_size *= std::sqrt((std::pow(energy_norm_overall, 2) + std::pow(error_overall, 2))/mElementNumber) * mTargetError;
        else
            new_element_size *= std::sqrt((energy_norm_overall*energy_norm_overall+error_overall*error_overall)/number_of_elements)*mTargetError;

        // Check if element sizes are in specified limits. If not, set them to the limit case
        if(new_element_size < mMinSize)
            new_element_size = mMinSize;
        if(new_element_size > mMaxSize)
            new_element_size = mMaxSize;

        rElement.SetValue(ELEMENT_H, new_element_size);
    });
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void MetricErrorProcess<TDim>::CalculateMetric()
{
    // Array of nodes
    NodesArrayType& r_nodes_array = mrThisModelPart.Nodes();

    // We do a find of neighbours
    {
        FindNodalNeighboursProcess find_neighbours(mrThisModelPart);
        if (r_nodes_array.begin()->Has(NEIGHBOUR_ELEMENTS)) find_neighbours.ClearNeighbours();
        find_neighbours.Execute();
    }

    // Iteration over all nodes
    KRATOS_DEBUG_ERROR_IF(r_nodes_array.size() == 0) <<  "ERROR:: Empty list of nodes" << std::endl;

    // Auxiliar variables
    const auto average_nodal_h = mAverageNodalH;
    const auto echo_level = mEchoLevel;

    block_for_each(r_nodes_array,
        [&average_nodal_h,&echo_level](NodeType& rNode) {
        /**************************************************************************
        ** Determine nodal element size h:
        ** if average_nodal_h == true : the nodal element size is averaged from the element size of neighboring elements
        ** if average_nodal_h == false: the nodal element size is the minimum element size from neighboring elements
        */
        double h_min = 0.0;
        auto& neigh_elements = rNode.GetValue(NEIGHBOUR_ELEMENTS);
        for(auto i_neighbour_elements = neigh_elements.begin(); i_neighbour_elements != neigh_elements.end(); i_neighbour_elements++){
            const double element_h = i_neighbour_elements->GetValue(ELEMENT_H);
            if(average_nodal_h == false) {
                if(h_min == 0.0 || h_min > element_h)
                    h_min = element_h;
            } else {
                h_min += element_h;
            }
        }

        // Average Nodal H
        if(average_nodal_h) h_min = h_min/static_cast<double>(neigh_elements.size());

        // Setting value
        rNode.SetValue(METRIC_SCALAR, h_min);

        KRATOS_INFO_IF("MetricErrorProcess", echo_level > 2) << "Node " << rNode.Id() << " has metric: "<< h_min << std::endl;
    });
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
const Parameters MetricErrorProcess<TDim>::GetDefaultParameters() const
{
    /**
     * We configure using the following parameters:
     * minimal_size: The minimal size to consider on the remeshing
     * maximal_size: The maximal size to consider on the remeshing
     * target_error: The target error
     * set_target_number_of_elements: If the number of elements will be forced or not
     * target_number_of_elements: The estimated/desired number of elements
     * perform_nodal_h_averaging: If the nodal size to consider will be averaged over the mesh
     * echo_level: The verbosity
     */
    const Parameters default_parameters = Parameters(R"(
    {
        "minimal_size"                        : 0.01,
        "maximal_size"                        : 1.0,
        "error_strategy_parameters":
        {
            "target_error"                        : 0.01,
            "set_target_number_of_elements"       : false,
            "target_number_of_elements"           : 1000,
            "perform_nodal_h_averaging"           : false
        },
        "echo_level"                          : 0
    })"
    );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class MetricErrorProcess<2>;
template class MetricErrorProcess<3>;

};// namespace Kratos.
