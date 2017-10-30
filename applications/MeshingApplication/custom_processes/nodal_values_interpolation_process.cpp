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

// System includes

// External includes

// Project includes

// Include the point locator
#include "custom_processes/nodal_values_interpolation_process.h"

namespace Kratos
{
template<unsigned int TDim>
NodalValuesInterpolationProcess<TDim>::NodalValuesInterpolationProcess(
        ModelPart& rOriginMainModelPart,
        ModelPart& rDestinationMainModelPart,
        Parameters ThisParameters
        ):mrOriginMainModelPart(rOriginMainModelPart),
          mrDestinationMainModelPart(rDestinationMainModelPart)
{
    Parameters DefaultParameters = Parameters(R"(
    {
    "echo_level"            : 1, 
    "framework"             : "Eulerian", 
    "max_number_of_searchs" : 1000, 
    "step_data_size"        : 0, 
    "buffer_size"           : 0
    })");
    ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
    
    mEchoLevel = ThisParameters["echo_level"].GetInt();
    mFramework = ConvertFramework(ThisParameters["framework"].GetString());
    mMaxNumberOfResults = ThisParameters["max_number_of_searchs"].GetInt();
    mStepDataSize = ThisParameters["step_data_size"].GetInt();
    mBufferSize   = ThisParameters["buffer_size"].GetInt();

    if (mEchoLevel > 0) std::cout << "Step data size: " << mStepDataSize << " Buffer size: " << mBufferSize << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void NodalValuesInterpolationProcess<TDim>::Execute()
{
    // We create the locator
    BinBasedFastPointLocator<TDim> point_locator = BinBasedFastPointLocator<TDim>(mrOriginMainModelPart);
    point_locator.UpdateSearchDatabase();
    
    // Iterate in the nodes
    NodesArrayType& nodes_array = mrDestinationMainModelPart.Nodes();
    auto num_nodes = nodes_array.end() - nodes_array.begin();
    
    /* Nodes */
//         #pragma omp parallel for 
    for(unsigned int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        
        Vector shape_functions;
        Element::Pointer p_element;
        
        const array_1d<double, 3>& coordinates = it_node->Coordinates();
        const bool is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, mMaxNumberOfResults, 5.0e-2);
        
        if (is_found == false)
        {
            if (mEchoLevel > 0 || mFramework == Lagrangian) // NOTE: In the case we are in a Lagrangian framework this is serious and should print a message
            {
                std::cout << "WARNING: Node "<< it_node->Id() << " not found (interpolation not posible)" << std::endl;
                std::cout << "\t X:"<< it_node->X() << "\t Y:"<< it_node->Y() << "\t Z:"<< it_node->Z() << std::endl;
                
                if (mFramework == Lagrangian)
                {
                    std::cout << "WARNING: YOU ARE IN A LAGRANGIAN FRAMEWORK THIS IS DANGEROUS" << std::endl;
                }
            }
        }
        else
        {
            for(unsigned int i_step = 0; i_step < mBufferSize; ++i_step)
            {
                CalculateStepData(*(it_node.base()), p_element, shape_functions, i_step);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void NodalValuesInterpolationProcess<2>::CalculateStepData(
    NodeType::Pointer pNode,
    const Element::Pointer& pElement,
    const Vector& ShapeFunctions,
    const unsigned int Step
    )
{
    double* step_data = pNode->SolutionStepData().Data(Step);
    
    const double* node_data_0 = pElement->GetGeometry()[0].SolutionStepData().Data(Step);
    const double* node_data_1 = pElement->GetGeometry()[1].SolutionStepData().Data(Step);
    const double* node_data_2 = pElement->GetGeometry()[2].SolutionStepData().Data(Step);
    
    for (unsigned int j = 0; j < mStepDataSize; ++j)
    {
        step_data[j] = ShapeFunctions[0] * node_data_0[j]
                     + ShapeFunctions[1] * node_data_1[j]
                     + ShapeFunctions[2] * node_data_2[j];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void NodalValuesInterpolationProcess<3>::CalculateStepData(
    NodeType::Pointer pNode,
    const Element::Pointer& pElement,
    const Vector& ShapeFunctions,
    const unsigned int Step
    )
{
    double* step_data = pNode->SolutionStepData().Data(Step);
    
    // NOTE: This just works with tetrahedron (you are going to have problems with anything else)
    const double* node_data_0 = pElement->GetGeometry()[0].SolutionStepData().Data(Step);
    const double* node_data_1 = pElement->GetGeometry()[1].SolutionStepData().Data(Step);
    const double* node_data_2 = pElement->GetGeometry()[2].SolutionStepData().Data(Step);
    const double* node_data_3 = pElement->GetGeometry()[3].SolutionStepData().Data(Step);
    
    for (unsigned int j = 0; j < mStepDataSize; ++j)
    {
        step_data[j] = ShapeFunctions[0] * node_data_0[j]
                     + ShapeFunctions[1] * node_data_1[j]
                     + ShapeFunctions[2] * node_data_2[j]
                     + ShapeFunctions[3] * node_data_3[j];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
FrameworkEulerLagrange NodalValuesInterpolationProcess<TDim>::ConvertFramework(const std::string& str)
{
    if(str == "Lagrangian") 
    {
        return Lagrangian;
    }
    else if(str == "Eulerian") 
    {
        return Eulerian;
    }
    else
    {
        return Eulerian;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class NodalValuesInterpolationProcess<2>;
template class NodalValuesInterpolationProcess<3>;

}  // namespace Kratos.
