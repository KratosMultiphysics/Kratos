//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Philipp Bucher
//                   Vicente Mataix Ferrandiz
//

// System includes
#include <numeric>

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "tests/test_utilities/cpp_tests_utilities.h"

namespace Kratos::Testing {

    void GenerateGenericModelPart(ModelPart& rModelPart)
    {
        Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

        CppTestsUtilities::Create2DGeometry(rModelPart, "Element2D3N");

        rModelPart.CreateNewCondition("LineCondition2D2N", 1, {{1,2}}, p_elem_prop);
        rModelPart.CreateNewCondition("LineCondition2D2N", 2, {{1,4}}, p_elem_prop);
        rModelPart.CreateNewCondition("LineCondition2D2N", 3, {{2,5}}, p_elem_prop);
        rModelPart.CreateNewCondition("LineCondition2D2N", 4, {{5,6}}, p_elem_prop);

        std::vector<Node::Pointer> condition_nodes_3 (2);
        condition_nodes_3[0] = rModelPart.pGetNode(5);
        condition_nodes_3[1] = rModelPart.pGetNode(6);

        rModelPart.CreateNewGeometry("Line2D2", 1, {{1,2}});
        rModelPart.CreateNewGeometry("Line2D2", 2, rModelPart.pGetCondition(1)->pGetGeometry());
        rModelPart.CreateNewGeometry("Line2D2", 3, rModelPart.pGetCondition(2)->pGetGeometry());
        rModelPart.CreateNewGeometry("Line2D2", 4, PointerVector<Node>{condition_nodes_3});

        std::vector<Node::Pointer> element_nodes_0 (3);
        element_nodes_0[0] = rModelPart.pGetNode(1);
        element_nodes_0[1] = rModelPart.pGetNode(2);
        element_nodes_0[2] = rModelPart.pGetNode(3);

        rModelPart.CreateNewGeometry("Triangle2D3", PointerVector<Node>{element_nodes_0});
        rModelPart.CreateNewGeometry("Triangle2D3", rModelPart.pGetElement(1)->pGetGeometry());
        rModelPart.CreateNewGeometry("Triangle2D3", "Geometry_7", {{2,5,3}});
        rModelPart.CreateNewGeometry("Triangle2D3", "Geometry_8", {{5,6,3}});
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartDataValueContainer, KratosCoreFastSuite)
    {
        Model model;
        ModelPart& r_model_part = model.CreateModelPart("Main");
        r_model_part.SetValue(DENSITY,1.0);
        KRATOS_EXPECT_TRUE(r_model_part.Has(DENSITY));
        KRATOS_EXPECT_FALSE(r_model_part.Has(TEMPERATURE));
        KRATOS_EXPECT_DOUBLE_EQ(r_model_part.GetValue(DENSITY),1.0);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartFlag, KratosCoreFastSuite)
    {
        Model model;
        ModelPart& r_model_part = model.CreateModelPart("Main");

        r_model_part.Set(ACTIVE);
        KRATOS_EXPECT_TRUE(r_model_part.Is(ACTIVE));
        KRATOS_EXPECT_FALSE(r_model_part.Is(BOUNDARY));
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartSubModelPartsIterator, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        r_model_part.CreateSubModelPart("Inlet1");
        r_model_part.CreateSubModelPart("Inlet2");
        r_model_part.CreateSubModelPart("Outlet");
        r_model_part.CreateSubModelPart("AnotherOutlet");

        KRATOS_EXPECT_EQ(r_model_part.NumberOfSubModelParts(), 4);

        std::size_t id = 1;
        for(auto i_SubModelPart = r_model_part.SubModelPartsBegin() ; i_SubModelPart != r_model_part.SubModelPartsEnd() ; i_SubModelPart++){
            i_SubModelPart->CreateNewNode(id++, 0.00,0.00,0.00);
        }

        KRATOS_EXPECT_EQ(r_model_part.NumberOfNodes(), 4);
        KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("Inlet1").NumberOfNodes(), 1);
        KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("Outlet").NumberOfNodes(), 1);

        const auto& r_const_ref = r_model_part;
        auto r_smp_names = r_const_ref.GetSubModelPartNames();
        std::sort(r_smp_names.begin(), r_smp_names.end());
        const std::vector<std::string> r_smp_ref_names {"AnotherOutlet", "Inlet1", "Inlet2", "Outlet"};

        for (std::size_t i=0; i<r_smp_names.size(); ++i) {
            KRATOS_EXPECT_EQ(r_smp_names[i], r_smp_ref_names[i]);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartAddNodalSolutionStepVariable, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        r_model_part.AddNodalSolutionStepVariable(VELOCITY);

        r_model_part.CreateNewNode(123, 0.00,0.00,0.00);

        KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_model_part.AddNodalSolutionStepVariable(PRESSURE),
            "Error: Attempting to add the variable \"PRESSURE\" to the model part with name \"Main\" which is not empty");

        r_model_part.AddNodalSolutionStepVariable(VELOCITY); // Adding the same Variable twice is fine bcs it wont do anything
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartHasNodalSolutionStepVariable, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");


        r_model_part.AddNodalSolutionStepVariable(VELOCITY);

        KRATOS_EXPECT_TRUE(r_model_part.HasNodalSolutionStepVariable(VELOCITY));
        KRATOS_EXPECT_FALSE(r_model_part.HasNodalSolutionStepVariable(PRESSURE));
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartFullName, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        ModelPart& r_sub_model_part = r_model_part.CreateSubModelPart("SubModelPart");
        ModelPart& r_sub_sub_model_part = r_sub_model_part.CreateSubModelPart("SubSubModelPart");

        KRATOS_EXPECT_EQ(r_model_part.FullName(), "Main");
        KRATOS_EXPECT_EQ(r_sub_model_part.FullName(), "Main.SubModelPart");
        KRATOS_EXPECT_EQ(r_sub_sub_model_part.FullName(), "Main.SubModelPart.SubSubModelPart");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartEmptyName, KratosCoreFastSuite)
    {
        Model current_model;

        // Constructor with name
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(current_model.CreateModelPart(""),
            "Error: Please don't use empty names (\"\") when creating a ModelPart");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartBaseCreation, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Check results
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfNodes() == 6);
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfElements() == 4);
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfConditions() == 4);
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfGeometries() == 8);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartAddGeometries, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Generate submodelpart (I)
        auto& r_sub_model_1 = r_model_part.CreateSubModelPart("Sub1");
        auto& r_sub_sub_model_1 = r_sub_model_1.CreateSubModelPart("SubSub1");
        KRATOS_EXPECT_TRUE(r_sub_model_1.NumberOfGeometries() == 0);
        KRATOS_EXPECT_TRUE(r_sub_sub_model_1.NumberOfGeometries() == 0);

        // Copy one
        r_sub_model_1.AddGeometries(r_model_part.GeometriesBegin(), r_model_part.GeometriesEnd());

        // Check results
        KRATOS_EXPECT_TRUE(r_sub_model_1.NumberOfGeometries() == 8);
        KRATOS_EXPECT_TRUE(r_sub_sub_model_1.NumberOfGeometries() == 0);

        // Copy one
        r_sub_sub_model_1.AddGeometries(r_model_part.GeometriesBegin(), r_model_part.GeometriesEnd());

        // Check results
        KRATOS_EXPECT_TRUE(r_sub_model_1.NumberOfGeometries() == 8);
        KRATOS_EXPECT_TRUE(r_sub_sub_model_1.NumberOfGeometries() == 8);

        // Generate submodelpart (II)
        auto& r_sub_model_2 = r_model_part.CreateSubModelPart("Sub2");
        auto& r_sub_sub_model_2 = r_sub_model_2.CreateSubModelPart("SubSub2");
        KRATOS_EXPECT_TRUE(r_sub_model_2.NumberOfGeometries() == 0);
        KRATOS_EXPECT_TRUE(r_sub_sub_model_2.NumberOfGeometries() == 0);

        // Copy one
        const std::vector<std::size_t> indexes = {1,2,3,4};
        r_sub_model_2.AddGeometries(indexes);

        // Check results
        KRATOS_EXPECT_TRUE(r_sub_model_2.NumberOfGeometries() == 4);
        KRATOS_EXPECT_TRUE(r_sub_sub_model_2.NumberOfGeometries() == 0);

        // Copy one
        r_sub_sub_model_2.AddGeometries(indexes);

        // Check results
        KRATOS_EXPECT_TRUE(r_sub_model_2.NumberOfGeometries() == 4);
        KRATOS_EXPECT_TRUE(r_sub_sub_model_2.NumberOfGeometries() == 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartTable, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        Table<double>::Pointer p_table = Kratos::make_shared<Table<double>>();

        r_model_part.AddTable(0, p_table);

        Table<double>& r_table_model_part = r_model_part.GetTable(0);

        for (std::size_t i = 0; i < 6; ++i)
            r_table_model_part.PushBack(static_cast<double>(i), 2.0 * static_cast<double>(i));

        double nearest = (r_table_model_part.GetNearestRow(2.1))[0];
        KRATOS_EXPECT_DOUBLE_EQ(nearest, 4.0);
        KRATOS_EXPECT_DOUBLE_EQ(r_table_model_part.GetValue(2.1), 4.2);
        KRATOS_EXPECT_DOUBLE_EQ(r_table_model_part(2.1), 4.2);
        KRATOS_EXPECT_DOUBLE_EQ(r_table_model_part.GetDerivative(2.1), 2.0);

        auto& r_data = r_table_model_part.Data();
        KRATOS_EXPECT_EQ(r_data.size(), 6);

        ModelPart& r_sub_model_part = r_model_part.CreateSubModelPart("Sub");

        r_sub_model_part.AddTable(1, p_table);

        Table<double>& r_table_sub_model_part = r_sub_model_part.GetTable(1);

        nearest = (r_table_sub_model_part.GetNearestRow(2.1))[0];
        KRATOS_EXPECT_DOUBLE_EQ(nearest, 4.0);
        KRATOS_EXPECT_DOUBLE_EQ(r_table_sub_model_part.GetValue(2.1), 4.2);
        KRATOS_EXPECT_DOUBLE_EQ(r_table_sub_model_part(2.1), 4.2);
        KRATOS_EXPECT_DOUBLE_EQ(r_table_sub_model_part.GetDerivative(2.1), 2.0);

        r_data = r_table_sub_model_part.Data();
        KRATOS_EXPECT_EQ(r_data.size(), 6);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveElements, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Set flags
        r_model_part.pGetElement(1)->Set(TO_ERASE, true);
        r_model_part.pGetElement(2)->Set(TO_ERASE, true);

        // Call method
        r_model_part.RemoveElements(TO_ERASE);

        // Check results
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfNodes() == 6);
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfElements() == 2);
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfConditions() == 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartAddElementsWithNodes, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        ModelPart& r_sub_model_part = r_model_part.CreateSubModelPart("SubMain");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Call method
        auto aux_util = AuxiliarModelPartUtilities(r_sub_model_part);
        std::vector<std::size_t> list_elements {{1, 2}};
        aux_util.AddElementsWithNodes(list_elements);

        // Check results
        KRATOS_EXPECT_TRUE(r_sub_model_part.NumberOfNodes() == 4);
        KRATOS_EXPECT_TRUE(r_sub_model_part.NumberOfElements() == 2);
        KRATOS_EXPECT_TRUE(r_sub_model_part.NumberOfConditions() == 0);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartAddConditionsWithNodes, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        ModelPart& r_sub_model_part = r_model_part.CreateSubModelPart("SubMain");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Call method
        auto aux_util = AuxiliarModelPartUtilities(r_sub_model_part);
        std::vector<std::size_t> list_conditions {{1, 2}};
        aux_util.AddConditionsWithNodes(list_conditions);

        // Check results
        KRATOS_EXPECT_TRUE(r_sub_model_part.NumberOfNodes() == 3);
        KRATOS_EXPECT_TRUE(r_sub_model_part.NumberOfElements() == 0);
        KRATOS_EXPECT_TRUE(r_sub_model_part.NumberOfConditions() == 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveElementsAndBelongings, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Set flags
        r_model_part.pGetElement(1)->Set(TO_ERASE, true);
        r_model_part.pGetElement(2)->Set(TO_ERASE, true);

        // Call method
        auto aux_util = AuxiliarModelPartUtilities(r_model_part);
        aux_util.RemoveElementsAndBelongings(TO_ERASE);

        // Check results
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfNodes() == 4);
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfElements() == 2);
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfConditions() == 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveConditions, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Set flags
        r_model_part.pGetCondition(1)->Set(TO_ERASE, true);
        r_model_part.pGetCondition(2)->Set(TO_ERASE, true);

        // Call method
        r_model_part.RemoveConditions(TO_ERASE);

        // Check results
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfNodes() == 6);
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfConditions() == 2);
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfElements() == 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveConditionsAndBelongings, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Set flags
        r_model_part.pGetCondition(1)->Set(TO_ERASE, true);
        r_model_part.pGetCondition(2)->Set(TO_ERASE, true);

        // Call method
        auto aux_util = AuxiliarModelPartUtilities(r_model_part);
        aux_util.RemoveConditionsAndBelongings(TO_ERASE);

        // Check results
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfNodes() == 3);
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfConditions() == 2);
        KRATOS_EXPECT_TRUE(r_model_part.NumberOfElements() == 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartEnsureModelPartOwnsProperties, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        Properties::Pointer p_elem_prop = Kratos::make_shared<Properties>(1);
        for (auto& r_cond : r_model_part.Conditions()) {
            r_cond.SetProperties(p_elem_prop);
        }
        for (auto& r_elem : r_model_part.Elements()) {
            r_elem.SetProperties(p_elem_prop);
        }

        KRATOS_EXPECT_EQ(r_model_part.NumberOfProperties(), 1);

        // Call method
        auto aux_util = AuxiliarModelPartUtilities(r_model_part);
        aux_util.EnsureModelPartOwnsProperties(false);

        // Check results
        KRATOS_EXPECT_EQ(r_model_part.NumberOfProperties(), 2);

        aux_util.EnsureModelPartOwnsProperties(true);

        // Check results
        KRATOS_EXPECT_EQ(r_model_part.NumberOfProperties(), 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartGetSubModelPart, KratosCoreFastSuite)
    {
        Model model;

        auto& model_part = model.CreateModelPart("Main");
        const auto& r_const_model_part = model_part;

        // Checking SubModelPart
        model_part.CreateSubModelPart("Inlet1");

        ModelPart& smp = model_part.GetSubModelPart("Inlet1");
        KRATOS_EXPECT_EQ("Inlet1", smp.Name());

        const ModelPart& c_smp = r_const_model_part.GetSubModelPart("Inlet1");
        KRATOS_EXPECT_EQ("Inlet1", c_smp.Name());

        KRATOS_EXPECT_EXCEPTION_IS_THROWN(model_part.GetSubModelPart("Random"),
            "Error: There is no sub model part with name \"Random\" in model part \"Main\"\nThe following sub model parts are available:");

        KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_const_model_part.GetSubModelPart("Random"),
            "Error: There is no sub model part with name \"Random\" in model part \"Main\"\nThe following sub model parts are available:");

        // Checking SubSubModelPart
        auto& ssmp = smp.CreateSubModelPart("sub_inlet");

        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1"));
        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1.sub_inlet"));
        KRATOS_EXPECT_TRUE(r_const_model_part.HasSubModelPart("Inlet1"));
        KRATOS_EXPECT_TRUE(r_const_model_part.HasSubModelPart("Inlet1.sub_inlet"));

        KRATOS_EXPECT_EQ("sub_inlet", model_part.GetSubModelPart("Inlet1.sub_inlet").Name());
        KRATOS_EXPECT_EQ("sub_inlet", r_const_model_part.GetSubModelPart("Inlet1.sub_inlet").Name());

        KRATOS_EXPECT_EXCEPTION_IS_THROWN(model_part.GetSubModelPart("Inlet1.random_sub_inlet"),
            "Error: There is no sub model part with name \"random_sub_inlet\" in model part \"Main.Inlet1\"\nThe following sub model parts are available:\n\t\"sub_inlet\"");

        KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_const_model_part.GetSubModelPart("Inlet1.random_sub_inlet"),
            "Error: There is no sub model part with name \"random_sub_inlet\" in model part \"Main.Inlet1\"\nThe following sub model parts are available:\n\t\"sub_inlet\"");

        // Checking SubSubSubModelPart
        ssmp.CreateSubModelPart("tiny_inlet");

        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1.sub_inlet.tiny_inlet"));
        KRATOS_EXPECT_TRUE(r_const_model_part.HasSubModelPart("Inlet1.sub_inlet.tiny_inlet"));

        KRATOS_EXPECT_EQ("tiny_inlet", model_part.GetSubModelPart("Inlet1.sub_inlet.tiny_inlet").Name());
        KRATOS_EXPECT_EQ("tiny_inlet", r_const_model_part.GetSubModelPart("Inlet1.sub_inlet.tiny_inlet").Name());

        KRATOS_EXPECT_EXCEPTION_IS_THROWN(model_part.GetSubModelPart("Inlet1.sub_inlet.big_inlet"),
            "Error: There is no sub model part with name \"big_inlet\" in model part \"Main.Inlet1.sub_inlet\"\nThe following sub model parts are available:\n\t\"tiny_inlet\"");

        KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_const_model_part.GetSubModelPart("Inlet1.sub_inlet.big_inlet"),
            "Error: There is no sub model part with name \"big_inlet\" in model part \"Main.Inlet1.sub_inlet\"\nThe following sub model parts are available:\n\t\"tiny_inlet\"");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartHasSubModelPart, KratosCoreFastSuite)
    {
        Model model;

        auto& model_part = model.CreateModelPart("Main");

        model_part.CreateSubModelPart("Inlet1");

        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1"));

        ModelPart& smp = model_part.GetSubModelPart("Inlet1");
        auto& ssmp = smp.CreateSubModelPart("SubInlet");
        ssmp.CreateSubModelPart("SuperSubInlet");

        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1.SubInlet"));
        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1.SubInlet.SuperSubInlet"));
        KRATOS_EXPECT_TRUE(smp.HasSubModelPart("SubInlet"));
        KRATOS_EXPECT_TRUE(smp.HasSubModelPart("SubInlet.SuperSubInlet"));
        KRATOS_EXPECT_TRUE(ssmp.HasSubModelPart("SuperSubInlet"));

        KRATOS_EXPECT_FALSE(model_part.HasSubModelPart("SubInlet"));
        KRATOS_EXPECT_FALSE(model_part.HasSubModelPart("SuperSubInlet"));
        KRATOS_EXPECT_FALSE(smp.HasSubModelPart("SuperSubInlet"));

        KRATOS_EXPECT_FALSE(model_part.HasSubModelPart("Random"));
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartCreateSubModelPart, KratosCoreFastSuite)
    {
        Model model;

        auto& model_part = model.CreateModelPart("Main");

        // Checking SubModelPart
        auto& smp = model_part.CreateSubModelPart("Inlet1");
        KRATOS_EXPECT_EQ("Inlet1", smp.Name());

        // Checking SubSubModelPart
        // here Inlet1 exists already
        auto& ssmp = model_part.CreateSubModelPart("Inlet1.sub_inlet");

        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1"));
        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1.sub_inlet"));
        KRATOS_EXPECT_TRUE(smp.HasSubModelPart("sub_inlet"));

        // here InletCustom does NOT exists yet
        auto& ssmp_2 = model_part.CreateSubModelPart("InletCustom.ccc_sub_inlet");

        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("InletCustom"));
        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("InletCustom.ccc_sub_inlet"));
        KRATOS_EXPECT_TRUE(model_part.GetSubModelPart("InletCustom").HasSubModelPart("ccc_sub_inlet"));
        KRATOS_EXPECT_EQ(ssmp_2.Name(), "ccc_sub_inlet");

        // Checking SubSubSubModelPart
        auto& sssmp = model_part.CreateSubModelPart("Inlet1.sub_inlet.aabbcc");

        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1"));
        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1.sub_inlet"));
        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1.sub_inlet.aabbcc"));
        KRATOS_EXPECT_TRUE(ssmp.HasSubModelPart("aabbcc"));
        KRATOS_EXPECT_EQ(sssmp.Name(), "aabbcc");

        // here nothing exists yet
        auto& sssmp_2 = model_part.CreateSubModelPart("Fancy.xxx_sub_inlet.uztr");

        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Fancy"));
        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Fancy.xxx_sub_inlet"));
        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Fancy.xxx_sub_inlet.uztr"));
        KRATOS_EXPECT_TRUE(model_part.GetSubModelPart("Fancy").HasSubModelPart("xxx_sub_inlet"));
        KRATOS_EXPECT_TRUE(model_part.GetSubModelPart("Fancy").HasSubModelPart("xxx_sub_inlet.uztr"));
        KRATOS_EXPECT_EQ(sssmp_2.Name(), "uztr");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartRemoveSubModelPart, KratosCoreFastSuite)
    {
        Model model;

        auto& model_part = model.CreateModelPart("Main");

        // Checking SubModelPart
        auto& smp = model_part.CreateSubModelPart("Inlet1");
        model_part.CreateSubModelPart("Inlet1.sub_inlet.sub_sub_inlet");
        auto& ssmp_2 = model_part.CreateSubModelPart("test1.InletCustom.ccc_sub_inlet");

        ssmp_2.RemoveSubModelPart("ccc_sub_inlet");
        KRATOS_EXPECT_FALSE(ssmp_2.HasSubModelPart("ccc_sub_inlet"));

        model_part.RemoveSubModelPart("Inlet1.sub_inlet.sub_sub_inlet");
        KRATOS_EXPECT_FALSE(smp.HasSubModelPart("sub_inlet.sub_sub_inlet"));
        KRATOS_EXPECT_TRUE(smp.HasSubModelPart("sub_inlet"));
        KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1.sub_inlet"));

        model_part.RemoveSubModelPart("test1");
        KRATOS_EXPECT_FALSE(model_part.HasSubModelPart("test1"));
        KRATOS_EXPECT_FALSE(model_part.HasSubModelPart("test1.InletCustom"));
        KRATOS_EXPECT_FALSE(model_part.HasSubModelPart("test1.InletCustom"));

        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            model_part.RemoveSubModelPart("Inlet1.sub_inlet.sub_sub_inlet.test"),
            "Error: There is no sub model part with name \"sub_sub_inlet\" in model part \"Main.Inlet1.sub_inlet\"");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartSubRangeAddition, KratosCoreFastSuite)
    {
        Model model;
        auto& r_model_part = model.CreateModelPart("test");
        Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);

        for (IndexType i = 0; i < 16; ++i) {
            r_model_part.CreateNewNode(i + 1, 0.0, 0.0, 0.0);
        }

        KRATOS_EXPECT_EQ(r_model_part.Nodes().size(), 16);
        r_model_part.Nodes().shrink_to_fit();
        KRATOS_EXPECT_EQ(r_model_part.Nodes().capacity(), 16);

        auto& r_sub_model_part_1 = r_model_part.CreateSubModelPart("sub1");

        r_sub_model_part_1.AddNodes(r_model_part.Nodes().begin(), r_model_part.Nodes().end());
        // above should not do anything to the r_model_part, hence capacity should be the same.
        KRATOS_EXPECT_EQ(r_model_part.Nodes().capacity(), 16);

        KRATOS_EXPECT_EQ(r_sub_model_part_1.Nodes().size(), 16);
        r_sub_model_part_1.Nodes().shrink_to_fit();
        KRATOS_EXPECT_EQ(r_sub_model_part_1.Nodes().capacity(), 16);

        auto& r_sub_model_part_2 = r_sub_model_part_1.CreateSubModelPart("sub2");
        r_sub_model_part_2.AddNodes(r_model_part.Nodes().begin() + 8, r_model_part.Nodes().end());
        // above should not do anything to the r_model_part, hence capacity should be the same.
        KRATOS_EXPECT_EQ(r_model_part.Nodes().size(), 16);
        KRATOS_EXPECT_EQ(r_model_part.Nodes().capacity(), 16);
        KRATOS_EXPECT_EQ(r_sub_model_part_1.Nodes().size(), 16);
        // the capacity of the r_sub_model_part_1 may be changed because, the r_sub_model_part_2
        // used the iterators of the r_model_part. Then r_sub_model_part_1 is not a subset of r_model_part
        // because even though they are pointing to the same memory locations, the intrusive_ptrs memory locations
        // are not a subset. Therefore not checking for the capacity.
        KRATOS_EXPECT_EQ(r_sub_model_part_2.Nodes().size(), 8);
        r_sub_model_part_2.Nodes().shrink_to_fit();
        KRATOS_EXPECT_EQ(r_sub_model_part_2.Nodes().capacity(), 8);

        // now we add using a sub range to sub2
        auto& r_sub_model_part_3 = r_sub_model_part_1.CreateSubModelPart("sub3");
        r_sub_model_part_3.AddNodes(r_model_part.Nodes().begin() + 4, r_model_part.Nodes().end());
        // above should not do anything to the r_model_part, hence capacity should be the same.
        KRATOS_EXPECT_EQ(r_model_part.Nodes().size(), 16);
        KRATOS_EXPECT_EQ(r_model_part.Nodes().capacity(), 16);
        KRATOS_EXPECT_EQ(r_sub_model_part_1.Nodes().size(), 16);

        // again here the capacity changes because not r_sub_model_part_1 contains intrusive_ptrs
        // which are not a sub set of the r_model_part
        auto r_sub_model_part_1_capacity = r_sub_model_part_1.Nodes().capacity();
        KRATOS_EXPECT_EQ(r_sub_model_part_2.Nodes().size(), 8);
        KRATOS_EXPECT_EQ(r_sub_model_part_2.Nodes().capacity(), 8);
        KRATOS_EXPECT_EQ(r_sub_model_part_3.Nodes().size(), 12);
        r_sub_model_part_3.Nodes().shrink_to_fit();
        KRATOS_EXPECT_EQ(r_sub_model_part_3.Nodes().capacity(), 12);

        auto& r_sub_model_part_4 = r_sub_model_part_3.CreateSubModelPart("sub4");
        r_sub_model_part_4.AddNodes(r_sub_model_part_3.Nodes().begin() + 4, r_sub_model_part_3.Nodes().end() - 1);

        // now there shouldn't be any change in the size or the capacity, because
        // now we added r_sub_model_part_3 items to the r_sub_model_part_4. Then it will add them to the
        // r_sub_model_part_4, and when it checks IsSubSet for added nodes with r_sub_model_part_3, then it will
        // be a subset, hence all the parent model part additions are ignored.
        KRATOS_EXPECT_EQ(r_model_part.Nodes().size(), 16);
        KRATOS_EXPECT_EQ(r_model_part.Nodes().capacity(), 16);
        KRATOS_EXPECT_EQ(r_sub_model_part_1.Nodes().size(), 16);
        KRATOS_EXPECT_EQ(r_sub_model_part_1.Nodes().capacity(), r_sub_model_part_1_capacity);
        KRATOS_EXPECT_EQ(r_sub_model_part_2.Nodes().size(), 8);
        KRATOS_EXPECT_EQ(r_sub_model_part_2.Nodes().capacity(), 8);
        KRATOS_EXPECT_EQ(r_sub_model_part_3.Nodes().size(), 12);
        KRATOS_EXPECT_EQ(r_sub_model_part_3.Nodes().capacity(), 12);
        KRATOS_EXPECT_EQ(r_sub_model_part_4.Nodes().size(), 7);
        r_sub_model_part_4.Nodes().shrink_to_fit();
        KRATOS_EXPECT_EQ(r_sub_model_part_4.Nodes().capacity(), 7);
    }
}  // namespace Kratos::Testing.
