// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Riccardo Tosi
//

// System includes

// External includes

// Project includes
#include "processes/find_nodal_neighbours_process.h"
#include "custom_processes/metrics_divergencefree_process.h"
#include "custom_utilities/meshing_utilities.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
template<SizeType TDim>
MetricDivergenceFreeProcess<TDim>::MetricDivergenceFreeProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters
        ):mrThisModelPart(rThisModelPart)
{
    /**
     * We configure using the following parameters:
     * minimal_size: The minimal size to consider on the remeshing
     * maximal_size: The maximal size to consider on the remeshing
     * target_divergencefree: The target divergence
     * set_target_number_of_elements: If the number of elements will be forced or not
     * target_number_of_elements: The estimated/desired number of elements
     * perform_nodal_h_averaging: If the nodal size to consider will be averaged over the mesh
     * echo_level: The verbosity
     */
    Parameters default_parameters = Parameters(R"(
    {
        "minimal_size"                        : 0.01,
        "maximal_size"                        : 10.0,
        "refinement_strategy"                 : "maximum_strategy",
        "maximum_strategy":
        {
            "target_refinement_coefficient"       : 0.01,
            "refinement_coefficient"              : 2,
            "reference_variable_name"             : "DIVERGENCE",
            "set_target_number_of_elements"       : false,
            "target_number_of_elements"           : 1000,
            "perform_nodal_h_averaging"           : false
        },
        "echo_level"                          : 0
    })"
    );

    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    mMinSize = ThisParameters["minimal_size"].GetDouble();
    mMaxSize = ThisParameters["maximal_size"].GetDouble();
    mRefinementStrategy = ThisParameters["refinement_strategy"].GetString();

    // maximum_strategy
    mReferenceVariable = ThisParameters["maximum_strategy"]["reference_variable_name"].GetString();
    mSetElementNumber = ThisParameters["maximum_strategy"]["set_target_number_of_elements"].GetBool();
    mTargetRefinementCoefficient = ThisParameters["maximum_strategy"]["target_refinement_coefficient"].GetDouble();
    mRefinementCoefficient = ThisParameters["maximum_strategy"]["refinement_coefficient"].GetDouble();
    mElementNumber = ThisParameters["maximum_strategy"]["target_number_of_elements"].GetInt();
    mAverageNodalH = ThisParameters["maximum_strategy"]["perform_nodal_h_averaging"].GetBool();
    mDivergenceFreeMaxValue = -1;

    mEchoLevel = ThisParameters["echo_level"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void MetricDivergenceFreeProcess<TDim>::Execute()
{
    /******************************************************************************
    --1-- Initialize metric --1--
    ******************************************************************************/
    // Tensor variable definition
    const Variable<TensorArrayType>& tensor_variable = KratosComponents<Variable<TensorArrayType>>::Get("METRIC_TENSOR_"+std::to_string(TDim)+"D");

    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    KRATOS_DEBUG_ERROR_IF(nodes_array.size() == 0) <<  "ERROR:: Empty list of nodes" << std::endl;
    ElementsArrayType& elements_array = mrThisModelPart.Elements();
    KRATOS_DEBUG_ERROR_IF(elements_array.size() == 0) <<  "ERROR:: Empty list of elements" << std::endl;

    if (nodes_array.begin()->Has(tensor_variable) == false) {
        const TensorArrayType zero_array(3 * (TDim - 1), 0.0);

        // Iteration over the nodes
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i)
            (nodes_array.begin() + i)->SetValue(tensor_variable, zero_array);
    }

    // /******************************************************************************
    // --2-- Initialize what is needed for each refinement strategy --2--
    // ******************************************************************************/

    // Selection refinement strategy
    if (mRefinementStrategy == "maximum_strategy") {
        // Iteration over all elements to find maximum value of divergence
        const int number_elements = static_cast<int>(elements_array.size());
        mDivergenceFreeMaxValue = -1;
        const auto& r_reference_var = KratosComponents<Variable<double>>::Get(mReferenceVariable);
        for(int i_elem = 0; i_elem < number_elements; ++i_elem) {
            auto it_elem = elements_array.begin() + i_elem;
            const double divergencefree_elem_value = abs(it_elem->GetValue(r_reference_var));
            if (divergencefree_elem_value > mDivergenceFreeMaxValue) {
                mDivergenceFreeMaxValue = divergencefree_elem_value;
            KRATOS_WATCH(mDivergenceFreeMaxValue);
            }
        }
    }

    // /******************************************************************************
    // --3-- Calculate metric --3--
    // ******************************************************************************/
    CalculateMetric();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void MetricDivergenceFreeProcess<TDim>::CalculateMetric()
{
    // Array of nodes
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();

    // Array of elements
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    // Tensor variable definition
    const Variable<TensorArrayType>& tensor_variable = KratosComponents<Variable<TensorArrayType>>::Get("METRIC_TENSOR_"+std::to_string(TDim)+"D");

    // We do a find of neighbours
    {
        FindNodalNeighboursProcess find_neighbours(mrThisModelPart);
        if (nodes_array.begin()->Has(NEIGHBOUR_ELEMENTS)) find_neighbours.ClearNeighbours();
        find_neighbours.Execute();
    }

    // Iteration over all nodes
    const int num_nodes = static_cast<int>(nodes_array.size());
    KRATOS_DEBUG_ERROR_IF(num_nodes == 0) <<  "ERROR:: Empty list of nodes" << std::endl;

    // Selection refinement strategy
    if (mRefinementStrategy == "maximum_strategy") {

        // Reference variable
        const auto& r_reference_var = KratosComponents<Variable<double>>::Get(mReferenceVariable);

        // Divergence over full domain
        double divergencedomain_value = 0;

        const int number_elements = static_cast<int>(elements_array.size());
        for(int i_elem = 0; i_elem < number_elements; ++i_elem) {
            auto it_elem = elements_array.begin() + i_elem;
            divergencedomain_value += pow(it_elem->GetValue(r_reference_var),2);
        }
        divergencedomain_value = sqrt(divergencedomain_value);

        #pragma omp parallel for
        for(int i_node = 0; i_node < num_nodes; ++i_node) {
            auto it_node = nodes_array.begin() + i_node;

            // MinSize by default
            double element_size = mMinSize;

            // Initialize Divergence
            double divergencefree_elem_value = 0;

            // Estimate element size
            const double nodal_h = it_node->GetValue(NODAL_H);

            // Get divergence value on the node: cycle the neigh elements and take maximum DIVERGENCE value
            auto& neigh_elements = it_node->GetValue(NEIGHBOUR_ELEMENTS);
            for(auto i_neighbour_elements = neigh_elements.begin(); i_neighbour_elements != neigh_elements.end(); i_neighbour_elements++) {
                if (i_neighbour_elements->GetValue(r_reference_var) > divergencefree_elem_value) {
                    divergencefree_elem_value = i_neighbour_elements->GetValue(r_reference_var);
                }
            }



            // Set element size
            // if (divergencefree_elem_value >= mTargetRefinementCoefficient*mDivergenceFreeMaxValue) {
            //     element_size = nodal_h/mRefinementCoefficient;
            // }
            // else {
            //     element_size = nodal_h;
            // }
            element_size = mRefinementCoefficient*nodal_h*(divergencedomain_value/sqrt(number_elements)/divergencefree_elem_value);


            // Set metric
            BoundedMatrix<double, TDim, TDim> metric_matrix = ZeroMatrix(TDim, TDim);
            for(IndexType i = 0;i < TDim; ++i)
                metric_matrix(i,i) = 1.0/std::pow(element_size, 2);

            // Transform metric matrix to a vector
            const TensorArrayType metric = MathUtils<double>::StressTensorToVector<MatrixType, TensorArrayType>(metric_matrix);

            // Setting value
            it_node->SetValue(tensor_variable, metric);

            KRATOS_INFO_IF("MetricDivergenceFreeProcess", mEchoLevel > 2) << "Node " << it_node->Id() << " has metric: "<< metric << std::endl;

        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class MetricDivergenceFreeProcess<2>;
template class MetricDivergenceFreeProcess<3>;

};// namespace Kratos.
