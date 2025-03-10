//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Bastian Devresse, https://github.com/bdevresse
//
//
//  This file includes functionalities from OpenSubdiv in case we cannot use OpenSubdiv.
//  This includes subdivision and ...
//
//


// Project includes
#include "subdivision_surfaces/catmull_clark.h"
#include "expression/literal_flat_expression.h"

// Include base h
#include "subdivision_utilities.h"

// other includes 
#include <string.h>

namespace Kratos {


template<class TContainerType>
ContainerExpression<TContainerType> SDSUtils::ProjectForward(
    const ContainerExpression<TContainerType>& rInputExpression)
{
    KRATOS_TRY

    ContainerExpression<TContainerType> output_container(*rInputExpression.pGetModelPart());

    return output_container;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> SDSUtils::ProjectBackward(
    const ContainerExpression<TContainerType>& rInputExpression)
{
    KRATOS_TRY

    ContainerExpression<TContainerType> output_container(*rInputExpression.pGetModelPart());

    return output_container;

    KRATOS_CATCH("");
}

std::vector<double> SDSUtils::CalculateMappingRelation(
    ModelPart& rControlPolygon,
    const ModelPart& rControlledMesh,
    const std::string SubdivScheme,
    const bool FixFreeEdges)
{
    KRATOS_TRY

    // const IndexType number_of_entities = rControlledMesh.NumberOfNodes();
    // SizeType num_control_pts = rControlPolygon.NumberOfNodes();
    // std::vector<IndexType> r_shape(num_control_pts);
    // r_shape = {0,1,2,3,4,5,6,7};
    // auto p_flat_data_expression = LiteralFlatExpression<double>::Create(number_of_entities, r_shape);
    // output_container.SetExpression(p_flat_data_expression);
    // auto& r_output_expression = *p_flat_data_expression;

    // const SizeType num_control_pts = rControlPolygon.NumberOfNodes();
    // const SizeType num_fe_nodes = rControlledMesh.NumberOfNodes();
    std::vector<double> output_data;
    if(SubdivScheme == "catmull_clark") {
        CatmullClarkSDS subdiv_surface = CatmullClarkSDS(rControlPolygon, rControlledMesh);
        subdiv_surface.CreateMappingMatrix(output_data, rControlPolygon, rControlledMesh, FixFreeEdges);
    }

    // KRATOS_INFO("Finished CreateMappingMatrix, now starting data allocation for output_container") << std::endl;
    // KRATOS_INFO("r_output_expression.size()") << r_output_expression.size() << std::endl;
    // r_output_expression.Info();

    // SizeType num_fe_nodes = rControlledMesh.NumberOfNodes();
    // for(IndexType fe_idx = 0; fe_idx < num_fe_nodes; ++fe_idx) {
    //     for(IndexType cp_idx = 0; cp_idx < num_control_pts; ++cp_idx) {
    //         KRATOS_INFO("r_output_expression.GetItemShape()") << r_output_expression.GetItemShape() << std::endl;
    //         KRATOS_INFO("fe_idx ") << fe_idx << " , cp_idx : " << cp_idx << " , fe_idx * num_control_pts + cp_idx : " << fe_idx * num_control_pts + cp_idx << std::endl;
    //         KRATOS_INFO("output_data[fe_idx * num_control_pts + cp_idx] ") << output_data[fe_idx * num_control_pts + cp_idx] << std::endl;
    //         r_output_expression.SetData(fe_idx, cp_idx, output_data[fe_idx * num_control_pts + cp_idx]);
    //         KRATOS_INFO("r_output_expression.GetItemComponentCount()") << r_output_expression.GetItemComponentCount() << std::endl;
    //     }
    // }

    return output_data;

    KRATOS_CATCH("");
}



// template instantiations
#define KRATOS_INSTANTIATE_CATMULL_CLARK_SDS_UTIL_METHODS(CONTAINER_TYPE)                                                               \
    template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> SDSUtils::ProjectForward(             \
        const ContainerExpression<CONTAINER_TYPE>&);                                                                    \
    template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> SDSUtils::ProjectBackward(            \
        const ContainerExpression<CONTAINER_TYPE>&);                                                                    
    // template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> SDSUtils::CalculateMappingRelation(   
    //     const ContainerExpression<CONTAINER_TYPE>&, ModelPart&,                                                                   
    //     const ModelPart&, const std::string, const bool);

KRATOS_INSTANTIATE_CATMULL_CLARK_SDS_UTIL_METHODS(ModelPart::NodesContainerType)
KRATOS_INSTANTIATE_CATMULL_CLARK_SDS_UTIL_METHODS(ModelPart::ConditionsContainerType)
KRATOS_INSTANTIATE_CATMULL_CLARK_SDS_UTIL_METHODS(ModelPart::ElementsContainerType)

#undef KRATOS_INSTANTIATE_SIGMOIDAL_PROJECTION_UTIL_METHODS

} // namespace kratos
