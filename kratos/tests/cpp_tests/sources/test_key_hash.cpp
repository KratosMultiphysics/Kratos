//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "includes/key_hash.h"
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"

namespace Kratos {
    namespace Testing {

    /**
     *  Here the VariableHasher is test
     */
    KRATOS_TEST_CASE_IN_SUITE(VariableHasher, KratosCoreFastSuite)
    {
        VariableHasher<Variable<double>> variable_hasher;

        KRATOS_CHECK_EQUAL(variable_hasher(TEMPERATURE), 9132515808512);
    }

    /**
     *  Here the pVariableHasher is test
     */
    KRATOS_TEST_CASE_IN_SUITE(pVariableHasher, KratosCoreFastSuite)
    {
        pVariableHasher<Variable<double>> variable_hasher;

        KRATOS_CHECK_EQUAL(variable_hasher(&TEMPERATURE), 9132515808512);
    }

    /**
     *  Here the IndexedObjectHasher is test
     */
    KRATOS_TEST_CASE_IN_SUITE(IndexedObjectHasher, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("test");
        auto p_node = r_model_part.CreateNewNode(1, 1., 0, 0);
        auto& r_node = *p_node;
        IndexedObjectHasher<Node<3>> indexed_object_hasher;

        KRATOS_CHECK_EQUAL(indexed_object_hasher(r_node), 1);
    }

    /**
     *  Here the IndexedObjectPointerHasher is test
     */
    KRATOS_TEST_CASE_IN_SUITE(IndexedObjectPointerHasher, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("test");
        auto p_node = r_model_part.CreateNewNode(1, 1., 0, 0);
        IndexedObjectPointerHasher<Node<3>::Pointer> indexed_object_hasher;

        KRATOS_CHECK_EQUAL(indexed_object_hasher(p_node), 1);
    }

    /**
     *  Here the VectorIntegerHasher is test
     */
    KRATOS_TEST_CASE_IN_SUITE(VectorIndexHasher, KratosCoreFastSuite)
    {
        VectorIndexHasher<std::vector<std::size_t>> vector_integer_hasher;

        std::vector<std::size_t> vector_to_generate_hash_1({6, 1});
        std::vector<std::size_t> vector_to_generate_hash_2({4, 3, 6, 9});
        KRATOS_CHECK_EQUAL(vector_integer_hasher(vector_to_generate_hash_1), 175247769174);
        KRATOS_CHECK_EQUAL(vector_integer_hasher(vector_to_generate_hash_2), 706246304684648);
    }

    /**
     *  Here the DofPointerHasher is test
     */
    KRATOS_TEST_CASE_IN_SUITE(DofPointerHasher, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("test");
        auto p_node = r_model_part.CreateNewNode(1, 1., 0, 0);
        p_node->AddDof(DISPLACEMENT_X, REACTION_X);
        DofPointerHasher dof_pointer_hasher;

        KRATOS_CHECK_EQUAL(dof_pointer_hasher(p_node->pGetDof(DISPLACEMENT_X)), 9560262452381);
    }

    /**
     *  Here the PairHasher is test
     */
    KRATOS_TEST_CASE_IN_SUITE(PairHasher, KratosCoreFastSuite)
    {
        PairHasher<std::size_t, std::size_t> pair_hasher;
        std::pair<std::size_t, std::size_t> pair_test({1,2});

        KRATOS_CHECK_EQUAL(pair_hasher(pair_test), 175247769363);
    }

}  // namespace Testing.
}  // namespace Kratos.
