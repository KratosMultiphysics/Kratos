//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi 
//                   Vicente Mataix Ferrandiz
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/variable_utils.h"
#include "utilities/geometry_utilities.h"
#include "processes/compute_nodal_gradient_process.h"

namespace Kratos
{
template<bool THistorical> 
void ComputeNodalGradientProcess<THistorical>::Execute()
{
    KRATOS_TRY;
    
    // Set to zero
    ClearGradient();
    
    const auto& it_element_begin = mrModelPart.ElementsBegin();
    const auto& r_first_element_geometry = it_element_begin->GetGeometry();
    const std::size_t dimension = r_first_element_geometry.WorkingSpaceDimension();
    const std::size_t local_space_dimension = r_first_element_geometry.LocalSpaceDimension();
    const std::size_t number_of_nodes = r_first_element_geometry.PointsNumber();
    
    // The integration points
    const auto& integration_method = r_first_element_geometry.GetDefaultIntegrationMethod();
    const auto& integration_points = r_first_element_geometry.IntegrationPoints(integration_method);
    const std::size_t number_of_integration_points = integration_points.size();
    
    Matrix DN_DX = ZeroMatrix(number_of_nodes, dimension);
    Vector N = ZeroVector(number_of_nodes);
    Matrix J0 = ZeroMatrix(dimension, local_space_dimension);
    
    #pragma omp parallel for firstprivate(DN_DX,  N, J0)
    for(int i=0; i<static_cast<int>(mrModelPart.Elements().size()); ++i) {
        auto it_elem = it_element_begin + i;
        auto& r_geometry = it_elem->GetGeometry();
    
        Vector values(number_of_nodes);
        if (mrOriginVariableDoubleList.size() > 0) {
            for(std::size_t i=0; i<number_of_nodes; ++i)
                values[i] = r_geometry[i].FastGetSolutionStepValue(mrOriginVariableDoubleList[0]);
        } else {
            for(std::size_t i=0; i<number_of_nodes; ++i)
                values[i] = r_geometry[i].FastGetSolutionStepValue(mrOriginVariableComponentsList[0]);
        }
        
        // The containers of the shape functions and the local gradients
        const auto& rNcontainer = r_geometry.ShapeFunctionsValues(integration_method);
        const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(integration_method);
        
        for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
            // Getting the shape functions
            noalias(N) = row(rNcontainer, point_number);
            
            // Getting the jacobians and local gradients
            GeometryUtils::JacobianOnInitialConfiguration(r_geometry, integration_points[point_number], J0);
            double detJ0;
            Matrix InvJ0;
            MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);
            const Matrix& rDN_De = rDN_DeContainer[point_number];
            GeometryUtils::ShapeFunctionsGradients(rDN_De, InvJ0, DN_DX);
            
            const Vector grad = prod(trans(DN_DX), values);
            const double gauss_point_volume = integration_points[point_number].Weight() * detJ0;
            
            for(std::size_t i=0; i<number_of_nodes; ++i) {
                for(std::size_t k=0; k<dimension; ++k) {
                    double& val = GetGradient(r_geometry, i,k);
                    
                    #pragma omp atomic
                    val += N[i] * gauss_point_volume*grad[k];
                }
                
                double& vol = r_geometry[i].GetValue(mrAreaVariable);
                
                #pragma omp atomic
                vol += N[i] * gauss_point_volume;
            }
        }
    }
    
    PonderateGradient();

    KRATOS_CATCH("")
}
    
/***********************************************************************************/
/***********************************************************************************/
    
template<>
ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    Variable<double>& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    // We push the list of double variables
    mrOriginVariableDoubleList.push_back(rOriginVariable);
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rGradientVariable, mrModelPart.Nodes());
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rAreaVariable )) << "Missing variable " <<  rAreaVariable << std::endl;
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    Variable<double>& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    // We push the list of double variables
    mrOriginVariableDoubleList.push_back(rOriginVariable);
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rGradientVariable )) << "Missing variable " << rGradientVariable << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rAreaVariable )) << "Missing variable " <<  rAreaVariable << std::endl;
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    ComponentType& rOriginVariable,
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    // We push the components list
    mrOriginVariableComponentsList.push_back(rOriginVariable);
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rGradientVariable, mrModelPart.Nodes());
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rAreaVariable )) << "Missing variable " <<  rAreaVariable << std::endl;
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    ComponentType& rOriginVariable,
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    // We push the components list
    mrOriginVariableComponentsList.push_back(rOriginVariable);

    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rGradientVariable )) << "Missing variable " << rGradientVariable << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rAreaVariable )) << "Missing variable " <<  rAreaVariable << std::endl;
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>::ClearGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->SetValue(mrAreaVariable, 0.0);
        it_node->FastGetSolutionStepValue(mrGradientVariable).clear();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>::ClearGradient()
{
    const array_1d<double, 3> aux_zero_vector = ZeroVector(3);
    
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->SetValue(mrAreaVariable, 0.0);
        it_node->SetValue(mrGradientVariable, aux_zero_vector);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
double& ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>::GetGradient(
    Element::GeometryType& rThisGeometry,
    unsigned int i, 
    unsigned int k
    )
{
    double& val = rThisGeometry[i].FastGetSolutionStepValue(mrGradientVariable)[k];
    
    return val;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
double& ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>::GetGradient(
    Element::GeometryType& rThisGeometry,
    unsigned int i, 
    unsigned int k
    )
{
    double& val = rThisGeometry[i].GetValue(mrGradientVariable)[k];
    
    return val;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrGradientVariable) /= it_node->GetValue(mrAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i)
    {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->GetValue(mrGradientVariable) /= it_node->GetValue(mrAreaVariable);
    }
}
/***********************************************************************************/
/***********************************************************************************/

template class ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>;
template class ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>;

} /* namespace Kratos.*/
