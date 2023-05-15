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
        IndexedObjectHasher<Node> indexed_object_hasher;

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
        IndexedObjectPointerHasher<Node::Pointer> indexed_object_hasher;

        KRATOS_CHECK_EQUAL(indexed_object_hasher(p_node), 1);
    }

    /**
     *  Here the VectorIntegerHasher is test
     */
    KRATOS_TEST_CASE_IN_SUITE(VectorIndexHasher, KratosCoreFastSuite)
    {
        VectorIndexHasher<std::vector<std::size_t>> vector_integer_hasher;

        std::vector<std::size_t> vector_to_generate_hash_1({6, 1});
        std::vector<std::size_t> vector_to_generate_hash_1_duplicated({6, 1});
        std::vector<std::size_t> vector_to_generate_hash_2({4, 3, 6, 9});
        std::vector<std::size_t> vector_to_generate_hash_2_duplicated({4, 3, 6, 9});
        KRATOS_CHECK_EQUAL(vector_integer_hasher(vector_to_generate_hash_1), vector_integer_hasher(vector_to_generate_hash_1_duplicated));
        KRATOS_CHECK_EQUAL(vector_integer_hasher(vector_to_generate_hash_2), vector_integer_hasher(vector_to_generate_hash_2_duplicated));
    }

    /**
     *  Here the DofPointerHasher is test
     */
    KRATOS_TEST_CASE_IN_SUITE(DofPointerHasher, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("test");
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(REACTION);
        auto p_node = r_model_part.CreateNewNode(1, 1., 0, 0);
        p_node->AddDof(DISPLACEMENT_X, REACTION_X);
        DofPointerHasher dof_pointer_hasher;

        KRATOS_CHECK_EQUAL(dof_pointer_hasher(p_node->pGetDof(DISPLACEMENT_X)), dof_pointer_hasher(r_model_part.GetNode(1).pGetDof(DISPLACEMENT_X)));
    }

    /**
     *  Here the PairHasher is test
     */
    KRATOS_TEST_CASE_IN_SUITE(PairHasher, KratosCoreFastSuite)
    {
        PairHasher<std::size_t, std::size_t> pair_hasher;
        std::pair<std::size_t, std::size_t> pair_test({1,2});
        std::pair<std::size_t, std::size_t> pair_test_duplicated({1,2});

        KRATOS_CHECK_EQUAL(pair_hasher(pair_test), pair_hasher(pair_test_duplicated));
    }

    /**
     *  Here the VariableHasher is test
     */
    KRATOS_TEST_CASE_IN_SUITE(VariableComparator, KratosCoreFastSuite)
    {
        VariableComparator<Variable<double>> variable_comparor;

        KRATOS_CHECK_IS_FALSE(variable_comparor(TEMPERATURE, PRESSURE));
    }

    /**
     *  Here the pVariableComparator is test
     */
    KRATOS_TEST_CASE_IN_SUITE(pVariableComparator, KratosCoreFastSuite)
    {
        pVariableComparator<Variable<double>> variable_comparor;

        KRATOS_CHECK_IS_FALSE(variable_comparor(&TEMPERATURE, &PRESSURE));
    }

    /**
     *  Here the IndexedObjectComparator is test
     */
    KRATOS_TEST_CASE_IN_SUITE(IndexedObjectComparator, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("test");
        auto p_node1 = r_model_part.CreateNewNode(1, 1., 0, 0);
        auto p_node2 = r_model_part.CreateNewNode(2, 2., 0, 0);
        auto& r_node1 = *p_node1;
        auto& r_node2 = *p_node2;
        IndexedObjectComparator<Node> indexed_object_comparor;

        KRATOS_CHECK(indexed_object_comparor(r_node1, r_model_part.GetNode(1)));
        KRATOS_CHECK_IS_FALSE(indexed_object_comparor(r_node1, r_node2));
    }

    /**
     *  Here the IndexedObjectPointerComparator is test
     */
    KRATOS_TEST_CASE_IN_SUITE(IndexedObjectPointerComparator, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("test");
        auto p_node1 = r_model_part.CreateNewNode(1, 1., 0, 0);
        auto p_node2 = r_model_part.CreateNewNode(2, 2., 0, 0);
        IndexedObjectPointerComparator<Node::Pointer> indexed_object_comparor;

        KRATOS_CHECK(indexed_object_comparor(p_node1, r_model_part.pGetNode(1)));
        KRATOS_CHECK_IS_FALSE(indexed_object_comparor(p_node1, p_node2));
    }

    /**
     *  Here the SharedPointerComparator is test
     */
    KRATOS_TEST_CASE_IN_SUITE(SharedPointerComparator, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("test");
        auto p_node1 = r_model_part.CreateNewNode(1, 1., 0, 0);
        auto p_node2 = r_model_part.CreateNewNode(2, 2., 0, 0);
        SharedPointerComparator<Node::Pointer> shared_pointer_comparor;

        KRATOS_CHECK(shared_pointer_comparor(p_node1, r_model_part.pGetNode(1)));
        KRATOS_CHECK_IS_FALSE(shared_pointer_comparor(p_node1, p_node2));
    }

    /**
     *  Here the VectorIndexComparor is test
     */
    KRATOS_TEST_CASE_IN_SUITE(VectorIndexComparor, KratosCoreFastSuite)
    {
        VectorIndexComparor<std::vector<std::size_t>> vector_integer_comparor;

        std::vector<std::size_t> vector_1({6, 1});
        std::vector<std::size_t> vector_1_duplicated({6, 1});
        std::vector<std::size_t> vector_2({4, 3, 6, 9});

        KRATOS_CHECK(vector_integer_comparor(vector_1, vector_1_duplicated));
        KRATOS_CHECK_IS_FALSE(vector_integer_comparor(vector_1, vector_2));
    }

    /**
     *  Here the DofPointerComparor is test
     */
    KRATOS_TEST_CASE_IN_SUITE(DofPointerComparor, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("test");
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(REACTION);
        auto p_node = r_model_part.CreateNewNode(1, 1., 0, 0);
        p_node->AddDof(DISPLACEMENT_X, REACTION_X);
        p_node->AddDof(DISPLACEMENT_Y, REACTION_Y);
        DofPointerComparor dof_pointer_comparor;

        KRATOS_CHECK_IS_FALSE(dof_pointer_comparor(p_node->pGetDof(DISPLACEMENT_X), p_node->pGetDof(DISPLACEMENT_Y)));
    }

    /**
     *  Here the PairComparor is test
     */
    KRATOS_TEST_CASE_IN_SUITE(PairComparor, KratosCoreFastSuite)
    {
        PairComparor<std::size_t, std::size_t> pair_comparor;
        std::pair<std::size_t, std::size_t> pair_test({1,2});
        std::pair<std::size_t, std::size_t> pair_test_duplicated({1,2});
        std::pair<std::size_t, std::size_t> pair_test_inverted({2,1});

        KRATOS_CHECK(pair_comparor(pair_test, pair_test_duplicated));
        KRATOS_CHECK_IS_FALSE(pair_comparor(pair_test, pair_test_inverted));
    }

}  // namespace Testing.
}  // namespace Kratos.
