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

template<class TContainerType>
ContainerExpression<TContainerType> SDSUtils::CalculateMappingRelation(
    const ContainerExpression<TContainerType>& rInputExpression,
    ModelPart& rControlPolygon,
    const ModelPart& rControlledMesh,
    const std::string SubdivScheme,
    const bool FixFreeEdges)
{
    KRATOS_TRY

    ContainerExpression<TContainerType> output_container(*rInputExpression.pGetModelPart());
    std::vector<double> output_data;
    if(SubdivScheme == "catmull_clark") {
        CatmullClarkSDS subdiv_surface = CatmullClarkSDS(rControlPolygon, rControlledMesh);
        subdiv_surface.CreateMappingMatrix(output_data, rControlPolygon, rControlledMesh, FixFreeEdges);
    }
    // CreateMappingMatrix(output_data, rControlPolygon, rControlledMesh, FixFreeEdges);

    return output_container;

    KRATOS_CATCH("");
}



// template instantiations
#define KRATOS_INSTANTIATE_CATMULL_CLARK_SDS_UTIL_METHODS(CONTAINER_TYPE)                                                               \
    template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> SDSUtils::ProjectForward(             \
        const ContainerExpression<CONTAINER_TYPE>&);                                                                    \
    template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> SDSUtils::ProjectBackward(            \
        const ContainerExpression<CONTAINER_TYPE>&);                                                                    \
    template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> SDSUtils::CalculateMappingRelation(   \
        const ContainerExpression<CONTAINER_TYPE>&, ModelPart&,                                                                   \
        const ModelPart&, const std::string, const bool);

KRATOS_INSTANTIATE_CATMULL_CLARK_SDS_UTIL_METHODS(ModelPart::NodesContainerType)
KRATOS_INSTANTIATE_CATMULL_CLARK_SDS_UTIL_METHODS(ModelPart::ConditionsContainerType)
KRATOS_INSTANTIATE_CATMULL_CLARK_SDS_UTIL_METHODS(ModelPart::ElementsContainerType)

#undef KRATOS_INSTANTIATE_SIGMOIDAL_PROJECTION_UTIL_METHODS

} // namespace kratos
