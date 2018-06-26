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
          mrDestinationMainModelPart(rDestinationMainModelPart)
{
    Parameters DefaultParameters = Parameters(R"(
    {
    "echo_level"                 : 1,
    "framework"                  : "Eulerian",
    "max_number_of_searchs"      : 1000,
    "interpolate_non_historical" : true,
    "step_data_size"             : 0,
    "buffer_size"                : 0
    })");
    ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
    
    mEchoLevel = ThisParameters["echo_level"].GetInt();
    mFramework = ConvertFramework(ThisParameters["framework"].GetString());
    mMaxNumberOfResults = ThisParameters["max_number_of_searchs"].GetInt();
    mInterpolateNonHistorical = ThisParameters["interpolate_non_historical"].GetBool();
    mStepDataSize = ThisParameters["step_data_size"].GetInt();
    mBufferSize   = ThisParameters["buffer_size"].GetInt();

    KRATOS_INFO_IF("NodalValuesInterpolationProcess", mEchoLevel > 0) << "Step data size: " << mStepDataSize << " Buffer size: " << mBufferSize << std::endl;
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
    
    if (mInterpolateNonHistorical)
        GetListNonHistoricalVariables();

    /* Nodes */
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(num_nodes); ++i) {
        auto it_node = nodes_array.begin() + i;
        
        Vector shape_functions;
        Element::Pointer p_element;
        
        const array_1d<double, 3>& coordinates = it_node->Coordinates();
        const bool is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, mMaxNumberOfResults, 5.0e-2);
        
        if (is_found == false) {
            if (mEchoLevel > 0 || mFramework == FrameworkEulerLagrange::LAGRANGIAN) { // NOTE: In the case we are in a Lagrangian framework this is serious and should print a message
                KRATOS_WARNING("NodalValuesInterpolationProcess") << "WARNING: Node "<< it_node->Id() << " not found (interpolation not posible)" << "\n\t X:"<< it_node->X() << "\t Y:"<< it_node->Y() << "\t Z:"<< it_node->Z() << std::endl;
                KRATOS_WARNING_IF("NodalValuesInterpolationProcess", mFramework == FrameworkEulerLagrange::LAGRANGIAN) << "WARNING: YOU ARE IN A LAGRANGIAN FRAMEWORK THIS IS DANGEROUS" << std::endl;
            }
        } else {
            if (mInterpolateNonHistorical)
                CalculateData(*(it_node.base()), p_element, shape_functions);
            for(IndexType i_step = 0; i_step < mBufferSize; ++i_step)
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

template<>
void NodalValuesInterpolationProcess<2>::CalculateData(
    NodeType::Pointer pNode,
    const Element::Pointer& pElement,
    const Vector& rShapeFunctions
    )
{
    // The nodal data (non-historical)
    auto& data = pNode->Data();

    // The nodal data (non-historical) of each node of the original mesh
    const auto& node_data_0 = pElement->GetGeometry()[0].Data();
    const auto& node_data_1 = pElement->GetGeometry()[1].Data();
    const auto& node_data_2 = pElement->GetGeometry()[2].Data();

    // Now we interpolate the values of each node
    double aux_coeff;
    for (auto& var : mListDoublesVariables) {
        aux_coeff = 0.0;
        if (node_data_0.Has(var)) aux_coeff += rShapeFunctions[0];
        if (node_data_1.Has(var)) aux_coeff += rShapeFunctions[1];
        if (node_data_2.Has(var)) aux_coeff += rShapeFunctions[2];
        if (aux_coeff > 0.0) {
            aux_coeff = 1.0/aux_coeff;
            data.SetValue(var,
      aux_coeff * (rShapeFunctions[0] * (node_data_0.GetValue(var))
                 + rShapeFunctions[1] * (node_data_1.GetValue(var))
                 + rShapeFunctions[2] * (node_data_2.GetValue(var))));
        }
    }
    for (auto& var : mListArraysVariables) {
            aux_coeff = 0.0;
            if (node_data_0.Has(var)) aux_coeff += rShapeFunctions[0];
            if (node_data_1.Has(var)) aux_coeff += rShapeFunctions[1];
            if (node_data_2.Has(var)) aux_coeff += rShapeFunctions[2];
        if (aux_coeff > 0.0)  {
            if (aux_coeff > 0.0) aux_coeff = 1.0/aux_coeff;
            data.SetValue<array_1d<double, 3>>(var,
      aux_coeff * (rShapeFunctions[0] * (node_data_0.GetValue(var))
                 + rShapeFunctions[1] * (node_data_1.GetValue(var))
                 + rShapeFunctions[2] * (node_data_2.GetValue(var))));
        }
    }
    for (auto& var : mListVectorVariables) {
            aux_coeff = 0.0;
            if (node_data_0.Has(var)) aux_coeff += rShapeFunctions[0];
            if (node_data_1.Has(var)) aux_coeff += rShapeFunctions[1];
            if (node_data_2.Has(var)) aux_coeff += rShapeFunctions[2];
        if (aux_coeff > 0.0)  {
            if (aux_coeff > 0.0) aux_coeff = 1.0/aux_coeff;
            data.SetValue<Vector>(var,
      aux_coeff * (rShapeFunctions[0] * (node_data_0.GetValue(var))
                 + rShapeFunctions[1] * (node_data_1.GetValue(var))
                 + rShapeFunctions[2] * (node_data_2.GetValue(var))));
        }
    }
    for (auto& var : mListMatrixVariables) {
            aux_coeff = 0.0;
            if (node_data_0.Has(var)) aux_coeff += rShapeFunctions[0];
            if (node_data_1.Has(var)) aux_coeff += rShapeFunctions[1];
            if (node_data_2.Has(var)) aux_coeff += rShapeFunctions[2];
        if (aux_coeff > 0.0)  {
            if (aux_coeff > 0.0) aux_coeff = 1.0/aux_coeff;
            data.SetValue<Matrix>(var,
      aux_coeff * (rShapeFunctions[0] * (node_data_0.GetValue(var))
                 + rShapeFunctions[1] * (node_data_1.GetValue(var))
                 + rShapeFunctions[2] * (node_data_2.GetValue(var))));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void NodalValuesInterpolationProcess<3>::CalculateData(
    NodeType::Pointer pNode,
    const Element::Pointer& pElement,
    const Vector& rShapeFunctions
    )
{
    // The nodal data (non-historical)
    auto& data = pNode->Data();

    // The nodal data (non-historical) of each node of the original mesh
    const auto& node_data_0 = pElement->GetGeometry()[0].Data();
    const auto& node_data_1 = pElement->GetGeometry()[1].Data();
    const auto& node_data_2 = pElement->GetGeometry()[2].Data();
    const auto& node_data_3 = pElement->GetGeometry()[3].Data();

    // Now we interpolate the values of each node
    double aux_coeff;
    for (auto& var : mListDoublesVariables) {
        aux_coeff = 0.0;
        if (node_data_0.Has(var)) aux_coeff += rShapeFunctions[0];
        if (node_data_1.Has(var)) aux_coeff += rShapeFunctions[1];
        if (node_data_2.Has(var)) aux_coeff += rShapeFunctions[2];
        if (node_data_3.Has(var)) aux_coeff += rShapeFunctions[3];
        if (aux_coeff > 0.0)  {
            if (aux_coeff > 0.0) aux_coeff = 1.0/aux_coeff;
            data.SetValue(var,
      aux_coeff * (rShapeFunctions[0] * (node_data_0.GetValue(var))
                 + rShapeFunctions[1] * (node_data_1.GetValue(var))
                 + rShapeFunctions[2] * (node_data_2.GetValue(var))
                 + rShapeFunctions[3] * (node_data_3.GetValue(var))));
        }
    }
    for (auto& var : mListArraysVariables) {
        aux_coeff = 0.0;
        if (node_data_0.Has(var)) aux_coeff += rShapeFunctions[0];
        if (node_data_1.Has(var)) aux_coeff += rShapeFunctions[1];
        if (node_data_2.Has(var)) aux_coeff += rShapeFunctions[2];
        if (node_data_3.Has(var)) aux_coeff += rShapeFunctions[3];
        if (aux_coeff > 0.0)  {
            if (aux_coeff > 0.0) aux_coeff = 1.0/aux_coeff;
            data.SetValue<array_1d<double, 3>>(var,
      aux_coeff * (rShapeFunctions[0] * (node_data_0.GetValue(var))
                 + rShapeFunctions[1] * (node_data_1.GetValue(var))
                 + rShapeFunctions[2] * (node_data_2.GetValue(var))
                 + rShapeFunctions[3] * (node_data_3.GetValue(var))));
        }
    }
    for (auto& var : mListVectorVariables) {
        aux_coeff = 0.0;
        if (node_data_0.Has(var)) aux_coeff += rShapeFunctions[0];
        if (node_data_1.Has(var)) aux_coeff += rShapeFunctions[1];
        if (node_data_2.Has(var)) aux_coeff += rShapeFunctions[2];
        if (node_data_3.Has(var)) aux_coeff += rShapeFunctions[3];
        if (aux_coeff > 0.0)  {
            data.SetValue<Vector>(var,
      aux_coeff * (rShapeFunctions[0] * (node_data_0.GetValue(var))
                 + rShapeFunctions[1] * (node_data_1.GetValue(var))
                 + rShapeFunctions[2] * (node_data_2.GetValue(var))
                 + rShapeFunctions[3] * (node_data_3.GetValue(var))));
        }
    }
    for (auto& var : mListMatrixVariables) {
        aux_coeff = 0.0;
        if (node_data_0.Has(var)) aux_coeff += rShapeFunctions[0];
        if (node_data_1.Has(var)) aux_coeff += rShapeFunctions[1];
        if (node_data_2.Has(var)) aux_coeff += rShapeFunctions[2];
        if (node_data_3.Has(var)) aux_coeff += rShapeFunctions[3];
        if (aux_coeff > 0.0)  {
            if (aux_coeff > 0.0) aux_coeff = 1.0/aux_coeff;
            data.SetValue<Matrix>(var,
      aux_coeff * (rShapeFunctions[0] * (node_data_0.GetValue(var))
                 + rShapeFunctions[1] * (node_data_1.GetValue(var))
                 + rShapeFunctions[2] * (node_data_2.GetValue(var))
                 + rShapeFunctions[3] * (node_data_3.GetValue(var))));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void NodalValuesInterpolationProcess<2>::CalculateStepData(
    NodeType::Pointer pNode,
    const Element::Pointer& pElement,
    const Vector& rShapeFunctions,
    const IndexType Step
    )
{
    // The nodal data (historical)
    double* step_data = pNode->SolutionStepData().Data(Step);

    // The nodal data (historical) of each node of the original mesh
    // NOTE: This just works with triangle (you are going to have problems with anything else)
    const double* node_data_0 = pElement->GetGeometry()[0].SolutionStepData().Data(Step);
    const double* node_data_1 = pElement->GetGeometry()[1].SolutionStepData().Data(Step);
    const double* node_data_2 = pElement->GetGeometry()[2].SolutionStepData().Data(Step);

    // Now we interpolate the values of each node
    for (IndexType j = 0; j < mStepDataSize; ++j) {
        step_data[j] = rShapeFunctions[0] * node_data_0[j]
                     + rShapeFunctions[1] * node_data_1[j]
                     + rShapeFunctions[2] * node_data_2[j];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void NodalValuesInterpolationProcess<3>::CalculateStepData(
    NodeType::Pointer pNode,
    const Element::Pointer& pElement,
    const Vector& rShapeFunctions,
    const IndexType Step
    )
{
    // The nodal data (historical)
    double* step_data = pNode->SolutionStepData().Data(Step);
    
    // The nodal data (historical) of each node of the original mesh
    // NOTE: This just works with tetrahedron (you are going to have problems with anything else)
    const double* node_data_0 = pElement->GetGeometry()[0].SolutionStepData().Data(Step);
    const double* node_data_1 = pElement->GetGeometry()[1].SolutionStepData().Data(Step);
    const double* node_data_2 = pElement->GetGeometry()[2].SolutionStepData().Data(Step);
    const double* node_data_3 = pElement->GetGeometry()[3].SolutionStepData().Data(Step);
    
    // Now we interpolate the values of each node
    for (IndexType j = 0; j < mStepDataSize; ++j) {
        step_data[j] = rShapeFunctions[0] * node_data_0[j]
                     + rShapeFunctions[1] * node_data_1[j]
                     + rShapeFunctions[2] * node_data_2[j]
                     + rShapeFunctions[3] * node_data_3[j];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class NodalValuesInterpolationProcess<2>;
template class NodalValuesInterpolationProcess<3>;

}  // namespace Kratos.
