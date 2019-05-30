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
     * echo_level: The verbosity
     */
    Parameters default_parameters = Parameters(R"(
    {
        "minimal_size"                        : 0.01,
        "maximal_size"                        : 10.0,
        "refinement_strategy"                 : "maximum_strategy",
        "mean_distribution_strategy":
        {
            "target_refinement_coefficient"       : 0.9,
            "refinement_bound"                    : 2.0,
            "reference_variable_name"             : "DIVERGENCE",
            "reference_norm_name"                 : "DIVERGENCE_H1SEMINORM",
            "reference_volume_name"               : "AUX_VOLUME"
        },

        "maximum_strategy":
        {
            "target_refinement_coefficient"       : 0.1,
            "refinement_coefficient"              : 2.0,
            "reference_variable_name"             : "DIVERGENCE"
        },
        "echo_level"                          : 0
    })"
    );

    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    mMinSize = ThisParameters["minimal_size"].GetDouble();
    mMaxSize = ThisParameters["maximal_size"].GetDouble();
    mRefinementStrategy = ThisParameters["refinement_strategy"].GetString();
    mEchoLevel = ThisParameters["echo_level"].GetInt();

    // Mean strategy
    mMeanStrategyReferenceVariable = ThisParameters["mean_distribution_strategy"]["reference_variable_name"].GetString();
    mMeanStrategyReferenceNorm = ThisParameters["mean_distribution_strategy"]["reference_norm_name"].GetString();
    mMeanStrategyReferenceVolume = ThisParameters["mean_distribution_strategy"]["reference_volume_name"].GetString();
    mMeanStrategyTargetRefinementCoefficient = ThisParameters["mean_distribution_strategy"]["target_refinement_coefficient"].GetDouble();
    mMeanStrategyRefinementBound = ThisParameters["mean_distribution_strategy"]["refinement_bound"].GetDouble();
    mMeanStrategyDivergenceFreeInterpolationValue = -1;

    // Maximum strategy
    mMaxStrategyReferenceVariable = ThisParameters["maximum_strategy"]["reference_variable_name"].GetString();
    mMaxStrategyTargetRefinementCoefficient = ThisParameters["maximum_strategy"]["target_refinement_coefficient"].GetDouble();
    mMaxStrategyRefinementCoefficient = ThisParameters["maximum_strategy"]["refinement_coefficient"].GetDouble();
    mDivergenceFreeMaxValue = -1;

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
    InitializeRefinementStrategy();

    // // Selection refinement strategy
    // if (mRefinementStrategy == "maximum_strategy") {
    //     // Iteration over all elements to find maximum value of divergence
    //     const int number_elements = static_cast<int>(elements_array.size());
    //     mDivergenceFreeMaxValue = -1;
    //     const auto& r_reference_var = KratosComponents<Variable<double>>::Get(mMaxStrategyReferenceVariable);
    //     for(int i_elem = 0; i_elem < number_elements; ++i_elem) {
    //         auto it_elem = elements_array.begin() + i_elem;
    //         const double divergencefree_elem_value = abs(it_elem->GetValue(r_reference_var));
    //         if (divergencefree_elem_value > mDivergenceFreeMaxValue) {
    //             mDivergenceFreeMaxValue = divergencefree_elem_value;
    //         KRATOS_WATCH(mDivergenceFreeMaxValue);
    //         }
    //     }
    // }

    // /******************************************************************************
    // --3-- Calculate metric --3--
    // ******************************************************************************/
    CalculateMetric();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void MetricDivergenceFreeProcess<TDim>::InitializeRefinementStrategy()
{
    // Prepare arrays
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    // Maximum refinement strategy
    if (mRefinementStrategy == "maximum_strategy") {
        // Iteration over all elements to find maximum value of reference variable
        const int number_elements = static_cast<int>(elements_array.size());
        mDivergenceFreeMaxValue = -1;
        const auto& r_reference_var = KratosComponents<Variable<double>>::Get(mMaxStrategyReferenceVariable);
        for(int i_elem = 0; i_elem < number_elements; ++i_elem) {
            auto it_elem = elements_array.begin() + i_elem;
            const double divergencefree_elem_value = abs(it_elem->GetValue(r_reference_var));
            if (divergencefree_elem_value > mDivergenceFreeMaxValue) {
                mDivergenceFreeMaxValue = divergencefree_elem_value;
            KRATOS_WATCH(mDivergenceFreeMaxValue);
            }
        }
    }

    // Mean strategy
    else if (mRefinementStrategy == "mean_distribution_strategy") {
        // Find of neighbours
        {
            FindNodalNeighboursProcess find_neighbours(mrThisModelPart);
            if (nodes_array.begin()->Has(NEIGHBOUR_ELEMENTS)) find_neighbours.ClearNeighbours();
            find_neighbours.Execute();
        }

        // Iteration over nodes
        const int number_nodes = static_cast<int>(nodes_array.size());

        mMeanStrategyDivergenceFreeInterpolationValue = -1;
        const auto& r_reference_var  = KratosComponents<Variable<double>>::Get(mMeanStrategyReferenceVariable);
        const auto& r_reference_vol  = KratosComponents<Variable<double>>::Get(mMeanStrategyReferenceVolume);
        for(int i_node = 0; i_node < number_nodes; ++i_node) {
            auto it_node = nodes_array.begin() + i_node;

            // Initialize Divergence
            double divergencefree_interp_value = 0;
            double local_volume = 0;

            // Get divergence value on the node: cycle the neigh elements and take sum of DIVERGENCE value
            auto& neigh_elements = it_node->GetValue(NEIGHBOUR_ELEMENTS);
            for(auto i_neighbour_elements = neigh_elements.begin(); i_neighbour_elements != neigh_elements.end(); i_neighbour_elements++) {
                divergencefree_interp_value += i_neighbour_elements->GetValue(r_reference_var)*i_neighbour_elements->GetValue(r_reference_vol);
                local_volume += i_neighbour_elements->GetValue(r_reference_vol);
            }
            divergencefree_interp_value /= local_volume;
            if (divergencefree_interp_value > mMeanStrategyDivergenceFreeInterpolationValue) {
                mMeanStrategyDivergenceFreeInterpolationValue = divergencefree_interp_value;
            KRATOS_WATCH(mMeanStrategyDivergenceFreeInterpolationValue);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void MetricDivergenceFreeProcess<TDim>::CalculateMetric()
{
    // Preapre arrays
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    // Tensor variable definition
    const Variable<TensorArrayType>& tensor_variable = KratosComponents<Variable<TensorArrayType>>::Get("METRIC_TENSOR_"+std::to_string(TDim)+"D");

    // Find of neighbours
    {
        FindNodalNeighboursProcess find_neighbours(mrThisModelPart);
        if (nodes_array.begin()->Has(NEIGHBOUR_ELEMENTS)) find_neighbours.ClearNeighbours();
        find_neighbours.Execute();
    }

    // Iteration over all nodes
    const int number_nodes = static_cast<int>(nodes_array.size());
    KRATOS_DEBUG_ERROR_IF(number_nodes == 0) <<  "ERROR:: Empty list of nodes" << std::endl;

    // Maximum refinement strategy
    if (mRefinementStrategy == "maximum_strategy") {

        // Reference variable
        const auto& r_reference_var = KratosComponents<Variable<double>>::Get(mMaxStrategyReferenceVariable);

        #pragma omp parallel for
        for(int i_node = 0; i_node < number_nodes; ++i_node) {
            auto it_node = nodes_array.begin() + i_node;

            // Initialize Divergence
            double divergencefree_elem_value = 0;

            // Estimate element size
            const double nodal_h = it_node->GetValue(NODAL_H);

            // Get divergence value on the node: cycle the neigh elements and take MAXIMUM value
            auto& neigh_elements = it_node->GetValue(NEIGHBOUR_ELEMENTS);
            for(auto i_neighbour_elements = neigh_elements.begin(); i_neighbour_elements != neigh_elements.end(); i_neighbour_elements++) {
                if (i_neighbour_elements->GetValue(r_reference_var) > divergencefree_elem_value) {
                    divergencefree_elem_value = i_neighbour_elements->GetValue(r_reference_var);
                }
            }

            // Set element size
            double element_size;
            if (divergencefree_elem_value >= mMaxStrategyTargetRefinementCoefficient*mDivergenceFreeMaxValue) {
                element_size = nodal_h/mMaxStrategyRefinementCoefficient;
            }
            else {
                element_size = nodal_h;
            }

            // Check with max and min size
            if (element_size < mMinSize) element_size = mMinSize;
            if (element_size > mMaxSize) element_size = mMaxSize;

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

    // Mean refinement strategy
    else if (mRefinementStrategy == "mean_distribution_strategy") {

        // Reference variable
        const auto& r_reference_var  = KratosComponents<Variable<double>>::Get(mMeanStrategyReferenceVariable);
        const auto& r_reference_norm = KratosComponents<Variable<double>>::Get(mMeanStrategyReferenceNorm);
        const auto& r_reference_vol  = KratosComponents<Variable<double>>::Get(mMeanStrategyReferenceVolume);

        // Divergence over full domain
        double divergencedomain_value = 0;

        const int number_elements = static_cast<int>(elements_array.size());
        for(int i_elem = 0; i_elem < number_elements; ++i_elem) {
            auto it_elem = elements_array.begin() + i_elem;
            divergencedomain_value += std::pow(it_elem->GetValue(r_reference_var),2);
        }
        divergencedomain_value = std::sqrt(divergencedomain_value);

        #pragma omp parallel for
        for(int i_node = 0; i_node < number_nodes; ++i_node) {
            auto it_node = nodes_array.begin() + i_node;

            // Initialize Divergence
            double divergencefree_interp_value = 0;
            double local_volume = 0;

            // Get divergence interpolation value on the node
            auto& neigh_elements = it_node->GetValue(NEIGHBOUR_ELEMENTS);
            for(auto i_neighbour_elements = neigh_elements.begin(); i_neighbour_elements != neigh_elements.end(); i_neighbour_elements++) {
                divergencefree_interp_value += i_neighbour_elements->GetValue(r_reference_var)*i_neighbour_elements->GetValue(r_reference_vol);
                local_volume += i_neighbour_elements->GetValue(r_reference_vol);
            }
            divergencefree_interp_value /= local_volume;

            // Estimate element size
            const double nodal_h = it_node->GetValue(NODAL_H);

            // Set element size
            double element_size;
            double factor = mMeanStrategyTargetRefinementCoefficient*divergencedomain_value/sqrt(number_elements)/divergencefree_interp_value;
            if (factor > mMeanStrategyRefinementBound) factor = mMeanStrategyRefinementBound;
            if (factor < 1.0/mMeanStrategyRefinementBound) factor = 1.0/mMeanStrategyRefinementBound;
            element_size = factor*nodal_h;

            // Check with max and min size
            if (element_size < mMinSize) element_size = mMinSize;
            if (element_size > mMaxSize) element_size = mMaxSize;

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
