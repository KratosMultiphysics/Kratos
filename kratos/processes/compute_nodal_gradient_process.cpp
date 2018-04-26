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
template< int TDim, class TVarType, HistoricalValues THist> 
void ComputeNodalGradientProcess<TDim, TVarType, THist>::Execute()
{
    KRATOS_TRY;
    
    // Set to zero
    ClearGradient();
    
    BoundedMatrix<double,TDim+1,  TDim> DN_DX;
    array_1d<double,TDim+1> N;
    double Volume;
    
    #pragma omp parallel for private(DN_DX,  N,  Volume)
    for(int i=0; i<static_cast<int>(mrModelPart.Elements().size()); ++i) {
        auto it_elem = mrModelPart.ElementsBegin()+i;
        Element::GeometryType& geom = it_elem->GetGeometry();
        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
        
        array_1d<double, TDim+1> values;
        for(unsigned int i=0; i<TDim+1; ++i)
            values[i] = geom[i].FastGetSolutionStepValue(mrOriginVariable);
        
        const array_1d<double,TDim> grad = prod(trans(DN_DX), values);
        
        for(unsigned int i=0; i<TDim+1; ++i) {
            for(unsigned int k=0; k<TDim; ++k) {
                double& val = GetGradient(geom, i,k);
                
                #pragma omp atomic
                val += N[i]*Volume*grad[k];
            }
            
            double& vol = geom[i].FastGetSolutionStepValue(mrAreaVariable);
            
            #pragma omp atomic
            vol += N[i]*Volume;
        }
    }
    
    PonderateGradient();

    KRATOS_CATCH("")
}
    
/***********************************************************************************/
/***********************************************************************************/
    
template<>
ComputeNodalGradientProcess<2, Variable<double>, Historical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    Variable<double>& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rGradientVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rAreaVariable, mrModelPart.Nodes());
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    Variable<double>& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rGradientVariable )) << "Missing variable " << rGradientVariable;
    VariableUtils().CheckVariableExists(rAreaVariable, mrModelPart.Nodes());
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<3, Variable<double>, Historical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    Variable<double>& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rGradientVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rAreaVariable, mrModelPart.Nodes());
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    Variable<double>& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rGradientVariable )) << "Missing variable " << rGradientVariable;
    VariableUtils().CheckVariableExists(rAreaVariable, mrModelPart.Nodes());
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<2, component_type, Historical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    component_type& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rGradientVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rAreaVariable, mrModelPart.Nodes());
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<2, component_type, NonHistorical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    component_type& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rGradientVariable )) << "Missing variable " << rGradientVariable;
    VariableUtils().CheckVariableExists(rAreaVariable, mrModelPart.Nodes());
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<3, component_type, Historical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    component_type& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rGradientVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rAreaVariable, mrModelPart.Nodes());
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<3, component_type, NonHistorical>::ComputeNodalGradientProcess(
    ModelPart& rModelPart, 
    component_type& rOriginVariable, 
    Variable<array_1d<double,3> >& rGradientVariable, 
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY
    
    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    KRATOS_ERROR_IF_NOT(rModelPart.Nodes().begin()->Has( rGradientVariable )) << "Missing variable " << rGradientVariable;
    VariableUtils().CheckVariableExists(rAreaVariable, mrModelPart.Nodes());
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ComputeNodalGradientProcess<2, Variable<double>, Historical>::ClearGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
        it_node->FastGetSolutionStepValue(mrGradientVariable).clear();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ComputeNodalGradientProcess<3, Variable<double>, Historical>::ClearGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
        it_node->FastGetSolutionStepValue(mrGradientVariable).clear();
    }
}

template<>
void ComputeNodalGradientProcess<2, component_type, Historical>::ClearGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
        it_node->FastGetSolutionStepValue(mrGradientVariable).clear();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ComputeNodalGradientProcess<3, component_type, Historical>::ClearGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
        it_node->FastGetSolutionStepValue(mrGradientVariable).clear();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>::ClearGradient()
{
    const array_1d<double, 3> AuxZeroVector = ZeroVector(3);
    
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
        it_node->SetValue(mrGradientVariable, AuxZeroVector);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>::ClearGradient()
{
    const array_1d<double, 3> AuxZeroVector = ZeroVector(3);
    
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
        it_node->SetValue(mrGradientVariable, AuxZeroVector);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<2, component_type, NonHistorical>::ClearGradient()
{
    const array_1d<double, 3> AuxZeroVector = ZeroVector(3);
    
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
        it_node->SetValue(mrGradientVariable, AuxZeroVector);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<3, component_type, NonHistorical>::ClearGradient()
{
    const array_1d<double, 3> AuxZeroVector = ZeroVector(3);
    
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
        it_node->SetValue(mrGradientVariable, AuxZeroVector);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
double& ComputeNodalGradientProcess<2, Variable<double>, Historical>::GetGradient(
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
double& ComputeNodalGradientProcess<3, Variable<double>, Historical>::GetGradient(
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
double& ComputeNodalGradientProcess<2, component_type, Historical>::GetGradient(
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
double& ComputeNodalGradientProcess<3, component_type, Historical>::GetGradient(
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
double& ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>::GetGradient(
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
double& ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>::GetGradient(
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
double& ComputeNodalGradientProcess<2, component_type, NonHistorical>::GetGradient(
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
double& ComputeNodalGradientProcess<3, component_type, NonHistorical>::GetGradient(
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
void ComputeNodalGradientProcess<2, Variable<double>, Historical>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrGradientVariable) /= it_node->FastGetSolutionStepValue(mrAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<3, Variable<double>, Historical>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrGradientVariable) /= it_node->FastGetSolutionStepValue(mrAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<2, component_type, Historical>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i)
    {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrGradientVariable) /= it_node->FastGetSolutionStepValue(mrAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<3, component_type, Historical>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i)
    {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrGradientVariable) /= it_node->FastGetSolutionStepValue(mrAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i)
    {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->GetValue(mrGradientVariable) /= it_node->FastGetSolutionStepValue(mrAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->GetValue(mrGradientVariable) /= it_node->FastGetSolutionStepValue(mrAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<2, component_type, NonHistorical>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->GetValue(mrGradientVariable) /= it_node->FastGetSolutionStepValue(mrAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<3, component_type, NonHistorical>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->GetValue(mrGradientVariable) /= it_node->FastGetSolutionStepValue(mrAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class ComputeNodalGradientProcess<2, Variable<double>, Historical>;
template class ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>;
template class ComputeNodalGradientProcess<3, Variable<double>, Historical>;
template class ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>;
template class ComputeNodalGradientProcess<2, component_type, Historical>;
template class ComputeNodalGradientProcess<2, component_type, NonHistorical>;
template class ComputeNodalGradientProcess<3, component_type, Historical>;
template class ComputeNodalGradientProcess<3, component_type, NonHistorical>;

} /* namespace Kratos.*/
