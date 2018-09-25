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
#include "custom_processes/nodal_values_interpolation_process.h"

namespace Kratos
{
template<SizeType TDim>
NodalValuesInterpolationProcess<TDim>::NodalValuesInterpolationProcess(
        ModelPart& rOriginMainModelPart,
        ModelPart& rDestinationMainModelPart,
        Parameters ThisParameters
        ):mrOriginMainModelPart(rOriginMainModelPart),
          mrDestinationMainModelPart(rDestinationMainModelPart),
          mThisParameters(ThisParameters)
{
    Parameters default_parameters = Parameters(R"(
    {
        "echo_level"                 : 1,
        "framework"                  : "Eulerian",
        "max_number_of_searchs"      : 1000,
        "interpolate_non_historical" : true,
        "extrapolate_contour_values" : true,
        "step_data_size"             : 0,
        "buffer_size"                : 0
    })");
    mThisParameters.ValidateAndAssignDefaults(default_parameters);

    KRATOS_INFO_IF("NodalValuesInterpolationProcess", mThisParameters["echo_level"].GetInt() > 0) << "Step data size: " << mThisParameters["step_data_size"].GetInt() << " Buffer size: " << mThisParameters["buffer_size"].GetInt() << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void NodalValuesInterpolationProcess<TDim>::Execute()
{
    // We create the locator
    BinBasedFastPointLocator<TDim> point_locator = BinBasedFastPointLocator<TDim>(mrOriginMainModelPart);
    point_locator.UpdateSearchDatabase();

    // Iterate in the nodes
    NodesArrayType& nodes_array = mrDestinationMainModelPart.Nodes();
    const SizeType num_nodes = nodes_array.end() - nodes_array.begin();

    if (mThisParameters["interpolate_non_historical"].GetBool())
        GetListNonHistoricalVariables();

    /* Nodes */
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(num_nodes); ++i) {
        auto it_node = nodes_array.begin() + i;

        Vector shape_functions;
        Element::Pointer p_element;

        const array_1d<double, 3>& coordinates = it_node->Coordinates();
        const bool is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, mThisParameters["max_number_of_searchs"].GetInt(), 5.0e-2);

        if (is_found == false) {
            if (mThisParameters["echo_level"].GetInt() > 0 || ConvertFramework(mThisParameters["framework"].GetString()) == FrameworkEulerLagrange::LAGRANGIAN) { // NOTE: In the case we are in a Lagrangian framework this is serious and should print a message
                KRATOS_WARNING("NodalValuesInterpolationProcess") << "WARNING: Node "<< it_node->Id() << " not found (interpolation not posible)" << "\n\t X:"<< it_node->X() << "\t Y:"<< it_node->Y() << "\t Z:"<< it_node->Z() << std::endl;
                KRATOS_WARNING_IF("NodalValuesInterpolationProcess", ConvertFramework(mThisParameters["framework"].GetString()) == FrameworkEulerLagrange::LAGRANGIAN) << "WARNING: YOU ARE IN A LAGRANGIAN FRAMEWORK THIS IS DANGEROUS" << std::endl;
            }
        } else {
            if (mThisParameters["interpolate_non_historical"].GetBool())
                CalculateData(*(it_node.base()), p_element, shape_functions);
            for(IndexType i_step = 0; i_step < mThisParameters["buffer_size"].GetInt(); ++i_step)
                CalculateStepData(*(it_node.base()), p_element, shape_functions, i_step);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void NodalValuesInterpolationProcess<TDim>::GetListNonHistoricalVariables()
{
    // We iterate over the model parts (in order to have the most extended possible list of variables)
    for (auto& submodel : mrOriginMainModelPart.SubModelParts()) {
        auto it_node = submodel.Nodes().begin();

        const auto& double_components = KratosComponents<Variable<double>>::GetComponents();

        for (auto& comp : double_components) {
            if (it_node->Has(*(comp.second))) {
                mListDoublesVariables.insert(*(comp.second));
            }
        }
        const auto& array_components = KratosComponents<Variable<array_1d<double, 3>>>::GetComponents();

        for (auto& comp : array_components) {
            if (it_node->Has(*(comp.second))) {
                mListArraysVariables.insert(*(comp.second));
            }
        }
        const auto& vector_components = KratosComponents<Variable<Vector>>::GetComponents();

        for (auto& comp : vector_components) {
            if (it_node->Has(*(comp.second))) {
                mListVectorVariables.insert(*(comp.second));
            }
        }
        const auto& matrix_components = KratosComponents<Variable<Matrix>>::GetComponents();

        for (auto& comp : matrix_components) {
            if (it_node->Has(*(comp.second))) {
                mListMatrixVariables.insert(*(comp.second));
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void NodalValuesInterpolationProcess<TDim>::CalculateData(
    NodeType::Pointer pNode,
    const Element::Pointer& pElement,
    const Vector& rShapeFunctions
    )
{
    // The nodal data (non-historical)
    DataValueContainer& data = pNode->Data();

    // The nodal data (non-historical) of each node of the original mesh
    GeometryType& geom = pElement->GetGeometry();
    const SizeType number_of_nodes = geom.size();
    std::vector<DataValueContainer> node_data(number_of_nodes);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        node_data[i] = geom[i].Data();
    }

    // Now we interpolate the values of each node
    double aux_coeff;
    for (auto& var : mListDoublesVariables) {
        aux_coeff = 0.0;
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            if (node_data[i].Has(var)) aux_coeff += rShapeFunctions[i];
        }
        if (aux_coeff > 0.0) {
            aux_coeff = 1.0/aux_coeff;
            double aux_value = 0.0;
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                aux_value += rShapeFunctions[i] * node_data[i].GetValue(var);
            }
            data.SetValue(var, aux_coeff * aux_value);
        }
    }
    for (auto& var : mListArraysVariables) {
            aux_coeff = 0.0;
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                if (node_data[i].Has(var)) aux_coeff += rShapeFunctions[i];
            }
        if (aux_coeff > 0.0)  {
            if (aux_coeff > 0.0) aux_coeff = 1.0/aux_coeff;
            array_1d<double, 3> aux_value(3, 0.0);
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                aux_value += rShapeFunctions[i] * node_data[i].GetValue(var);
            }
            data.SetValue<array_1d<double, 3>>(var, aux_coeff * aux_value);
        }
    }
    for (auto& var : mListVectorVariables) {
            aux_coeff = 0.0;
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                if (node_data[i].Has(var)) aux_coeff += rShapeFunctions[i];
            }
        if (aux_coeff > 0.0)  {
            if (aux_coeff > 0.0) aux_coeff = 1.0/aux_coeff;
            Vector aux_value = ZeroVector(node_data[0].GetValue(var).size());
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                aux_value += rShapeFunctions[i] * node_data[i].GetValue(var);
            }
            data.SetValue<Vector>(var, aux_coeff * aux_value);
        }
    }
    for (auto& var : mListMatrixVariables) {
            aux_coeff = 0.0;
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                if (node_data[i].Has(var)) aux_coeff += rShapeFunctions[i];
            }
        if (aux_coeff > 0.0)  {
            if (aux_coeff > 0.0) aux_coeff = 1.0/aux_coeff;
            Matrix aux_value = ZeroMatrix(node_data[0].GetValue(var).size1(), node_data[0].GetValue(var).size2());
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                aux_value += rShapeFunctions[i] * node_data[i].GetValue(var);
            }
            data.SetValue<Matrix>(var, aux_coeff * aux_value);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void NodalValuesInterpolationProcess<TDim>::CalculateStepData(
    NodeType::Pointer pNode,
    const Element::Pointer& pElement,
    const Vector& rShapeFunctions,
    const IndexType Step
    )
{
    // The nodal data (historical)
    double* step_data = pNode->SolutionStepData().Data(Step);
    for (IndexType j = 0; j < mThisParameters["step_data_size"].GetInt(); ++j)
        step_data[j] = 0;

    // The nodal data (historical) of each node of the original mesh
    GeometryType& geom = pElement->GetGeometry();
    const SizeType number_of_nodes = geom.size();
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const double* nodal_data = geom[i].SolutionStepData().Data(Step);
        // Now we interpolate the values of each node
        for (IndexType j = 0; j < mThisParameters["step_data_size"].GetInt(); ++j) {
            step_data[j] += rShapeFunctions[i] * nodal_data[j];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class NodalValuesInterpolationProcess<2>;
template class NodalValuesInterpolationProcess<3>;

}  // namespace Kratos.
