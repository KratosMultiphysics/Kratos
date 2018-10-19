//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "custom_utilities/interface_preprocessor.h"
#include "custom_mappers/nearest_element_mapper.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(InterfacePreprocessor_NodeBasedLocalSystems, KratosMappingApplicationSerialTestSuite)
{
    Node<3>::Pointer p_point1(new Node<3>(1, 0.00, 0.00, 0.00));
    Node<3>::Pointer p_point2(new Node<3>(2, 0.00, 10.00, 0.00));
    Node<3>::Pointer p_point3(new Node<3>(3, 10.00, 10.00, 0.00));
    Node<3>::Pointer p_point4(new Node<3>(4, 10.00, 0.00, 0.00));

    Quadrilateral2D4<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4);

    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions" : 3,
        "element_name"        : "Element2D3N",
        "create_skin_sub_model_part": false
    }  )");

    StructuredMeshGeneratorProcess(geometry, model_part, mesher_parameters).Execute();

    InterfacePreprocessor::MapperLocalSystemPointerVectorPointer p_mapper_local_systems(
        Kratos::make_shared<InterfacePreprocessor::MapperLocalSystemPointerVector>());

    const Kratos::unique_ptr<MapperLocalSystem> p_ref_local_system(Kratos::make_unique<NearestElementLocalSystem>());

    KRATOS_CHECK(p_ref_local_system->UseNodesAsBasis()); // this is true in the NearestElementLocalSystem

    InterfacePreprocessor interface_preprocess(model_part, p_mapper_local_systems);

    interface_preprocess.CreateMapperLocalSystems(p_ref_local_system);

    KRATOS_CHECK_EQUAL(model_part.NumberOfNodes(), p_mapper_local_systems->size());

    for (const auto& rp_local_sys : (*p_mapper_local_systems))
        KRATOS_CHECK_EQUAL(typeid(*p_ref_local_system), typeid(*rp_local_sys));
}

KRATOS_TEST_CASE_IN_SUITE(InterfacePreprocessor_UpdateInterface, KratosMappingApplicationSerialTestSuite)
{
    /*
    This test checks if the CreateMapperLocalSystems can be called multiple times
    This is needed in case the Interface is updated => remeshed
    */

    Node<3>::Pointer p_point1(new Node<3>(1, 0.00, 0.00, 0.00));
    Node<3>::Pointer p_point2(new Node<3>(2, 0.00, 10.00, 0.00));
    Node<3>::Pointer p_point3(new Node<3>(3, 10.00, 10.00, 0.00));
    Node<3>::Pointer p_point4(new Node<3>(4, 10.00, 0.00, 0.00));

    Quadrilateral2D4<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4);

    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Generated");

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions" : 3,
        "element_name"        : "Element2D3N",
        "create_skin_sub_model_part": false
    }  )");

    StructuredMeshGeneratorProcess(geometry, model_part, mesher_parameters).Execute();

    InterfacePreprocessor::MapperLocalSystemPointerVectorPointer p_mapper_local_systems(
        Kratos::make_shared<InterfacePreprocessor::MapperLocalSystemPointerVector>());

    const Kratos::unique_ptr<MapperLocalSystem> p_ref_local_system(Kratos::make_unique<NearestElementLocalSystem>());

    KRATOS_CHECK(p_ref_local_system->UseNodesAsBasis()); // this is true in the NearestElementLocalSystem

    InterfacePreprocessor interface_preprocess(model_part, p_mapper_local_systems);

    interface_preprocess.CreateMapperLocalSystems(p_ref_local_system);

    KRATOS_CHECK_EQUAL(model_part.NumberOfNodes(), p_mapper_local_systems->size());

    for (const auto& rp_local_sys : (*p_mapper_local_systems))
        KRATOS_CHECK_EQUAL(typeid(*p_ref_local_system), typeid(*rp_local_sys));

    model_part.GetMesh().Clear();

    Parameters mesher_parameters_2(R"(
    {
        "number_of_divisions" : 5,
        "element_name"        : "Element2D3N",
        "create_skin_sub_model_part": false
    }  )");

    StructuredMeshGeneratorProcess(geometry, model_part, mesher_parameters_2).Execute();

    interface_preprocess.CreateMapperLocalSystems(p_ref_local_system);

    KRATOS_CHECK_EQUAL(model_part.NumberOfNodes(), p_mapper_local_systems->size());

    for (const auto& rp_local_sys : (*p_mapper_local_systems))
        KRATOS_CHECK_EQUAL(typeid(*p_ref_local_system), typeid(*rp_local_sys));
}

}  // namespace Testing
}  // namespace Kratos