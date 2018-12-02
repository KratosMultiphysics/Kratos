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
#include <string>
#include <iostream>
#include <algorithm>

/* External includes */

/* Project includes */
#include "utilities/variable_utils.h"
#include "utilities/geometry_utilities.h"
#include "processes/compute_nodal_gradient_process.h"

namespace Kratos
{
template<class TVarType, HistoricalValues THist> 
void ComputeNodalGradientProcess<TVarType, THist>::Execute()
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
        for(std::size_t i=0; i<number_of_nodes; ++i)
            values[i] = r_geometry[i].FastGetSolutionStepValue(mrOriginVariable);
        
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
ComputeNodalGradientProcess<Variable<double>, Historical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    Variable<double>& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rGradientVariable, mrModelPart.Nodes());
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rAreaVariable )) << "Missing variable " <<  rAreaVariable << std::endl;
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<Variable<double>, NonHistorical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    Variable<double>& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rGradientVariable )) << "Missing variable " << rGradientVariable << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rAreaVariable )) << "Missing variable " <<  rAreaVariable << std::endl;
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<component_type, Historical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    component_type& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rGradientVariable, mrModelPart.Nodes());
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rAreaVariable )) << "Missing variable " <<  rAreaVariable << std::endl;
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<component_type, NonHistorical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    component_type& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rGradientVariable )) << "Missing variable " << rGradientVariable << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rAreaVariable )) << "Missing variable " <<  rAreaVariable << std::endl;
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ComputeNodalGradientProcess<Variable<double>, Historical>::ClearGradient()
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

template<>
void ComputeNodalGradientProcess<component_type, Historical>::ClearGradient()
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
void ComputeNodalGradientProcess<Variable<double>, NonHistorical>::ClearGradient()
{
    const array_1d<double, 3> AuxZeroVector = ZeroVector(3);
    
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->SetValue(mrAreaVariable, 0.0);
        it_node->SetValue(mrGradientVariable, AuxZeroVector);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<component_type, NonHistorical>::ClearGradient()
{
    const array_1d<double, 3> AuxZeroVector = ZeroVector(3);
    
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->SetValue(mrAreaVariable, 0.0);
        it_node->SetValue(mrGradientVariable, AuxZeroVector);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
double& ComputeNodalGradientProcess<Variable<double>, Historical>::GetGradient(
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
double& ComputeNodalGradientProcess<component_type, Historical>::GetGradient(
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
double& ComputeNodalGradientProcess<Variable<double>, NonHistorical>::GetGradient(
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
double& ComputeNodalGradientProcess<component_type, NonHistorical>::GetGradient(
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
void ComputeNodalGradientProcess<Variable<double>, Historical>::PonderateGradient()
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
void ComputeNodalGradientProcess<component_type, Historical>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i)
    {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrGradientVariable) /= it_node->GetValue(mrAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<Variable<double>, NonHistorical>::PonderateGradient()
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

template <>
void ComputeNodalGradientProcess<component_type, NonHistorical>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->GetValue(mrGradientVariable) /= it_node->GetValue(mrAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class ComputeNodalGradientProcess<Variable<double>, Historical>;
template class ComputeNodalGradientProcess<Variable<double>, NonHistorical>;
template class ComputeNodalGradientProcess<component_type, Historical>;
template class ComputeNodalGradientProcess<component_type, NonHistorical>;

} /* namespace Kratos.*/
