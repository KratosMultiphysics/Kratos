//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Philipp Bucher
//                   Vicente Mataix Ferrandiz
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/cpp_tests_utilities.h"

namespace Kratos {
  namespace Testing {

    typedef Node NodeType;

    void GenerateGenericModelPart(ModelPart& rModelPart)
    {
        Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

        CppTestsUtilities::Create2DGeometry(rModelPart, "Element2D3N");

        rModelPart.CreateNewCondition("LineCondition2D2N", 1, {{1,2}}, p_elem_prop);
        rModelPart.CreateNewCondition("LineCondition2D2N", 2, {{1,4}}, p_elem_prop);
        rModelPart.CreateNewCondition("LineCondition2D2N", 3, {{2,5}}, p_elem_prop);
        rModelPart.CreateNewCondition("LineCondition2D2N", 4, {{5,6}}, p_elem_prop);

        std::vector<NodeType::Pointer> condition_nodes_3 (2);
        condition_nodes_3[0] = rModelPart.pGetNode(5);
        condition_nodes_3[1] = rModelPart.pGetNode(6);

        rModelPart.CreateNewGeometry("Line2D2", 1, {{1,2}});
        rModelPart.CreateNewGeometry("Line2D2", 2, rModelPart.pGetCondition(1)->pGetGeometry());
        rModelPart.CreateNewGeometry("Line2D2", 3, rModelPart.pGetCondition(2)->pGetGeometry());
        rModelPart.CreateNewGeometry("Line2D2", 4, PointerVector<NodeType>{condition_nodes_3});

        std::vector<NodeType::Pointer> element_nodes_0 (3);
        element_nodes_0[0] = rModelPart.pGetNode(1);
        element_nodes_0[1] = rModelPart.pGetNode(2);
        element_nodes_0[2] = rModelPart.pGetNode(3);

        rModelPart.CreateNewGeometry("Triangle2D3", PointerVector<NodeType>{element_nodes_0});
        rModelPart.CreateNewGeometry("Triangle2D3", rModelPart.pGetElement(1)->pGetGeometry());
        rModelPart.CreateNewGeometry("Triangle2D3", "Geometry_7", {{2,5,3}});
        rModelPart.CreateNewGeometry("Triangle2D3", "Geometry_8", {{5,6,3}});
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartDataValueContainer, KratosCoreFastSuite)
    {
        Model model;
        ModelPart& r_model_part = model.CreateModelPart("Main");
        r_model_part.SetValue(DENSITY,1.0);
        KRATOS_CHECK(r_model_part.Has(DENSITY));
        KRATOS_CHECK_IS_FALSE(r_model_part.Has(TEMPERATURE));
        KRATOS_CHECK_DOUBLE_EQUAL(r_model_part.GetValue(DENSITY),1.0);
    }
    
    KRATOS_TEST_CASE_IN_SUITE(ModelPartFlag, KratosCoreFastSuite)
    {
        Model model;
        ModelPart& r_model_part = model.CreateModelPart("Main");
        
        r_model_part.Set(ACTIVE);
        KRATOS_CHECK(r_model_part.Is(ACTIVE));
        KRATOS_CHECK_IS_FALSE(r_model_part.Is(BOUNDARY));
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartSubModelPartsIterator, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        r_model_part.CreateSubModelPart("Inlet1");
        r_model_part.CreateSubModelPart("Inlet2");
        r_model_part.CreateSubModelPart("Outlet");
        r_model_part.CreateSubModelPart("AnotherOutlet");

        KRATOS_CHECK_EQUAL(r_model_part.NumberOfSubModelParts(), 4);

        std::size_t id = 1;
        for(auto i_SubModelPart = r_model_part.SubModelPartsBegin() ; i_SubModelPart != r_model_part.SubModelPartsEnd() ; i_SubModelPart++){
            i_SubModelPart->CreateNewNode(id++, 0.00,0.00,0.00);
        }

        KRATOS_CHECK_EQUAL(r_model_part.NumberOfNodes(), 4);
        KRATOS_CHECK_EQUAL(r_model_part.GetSubModelPart("Inlet1").NumberOfNodes(), 1);
        KRATOS_CHECK_EQUAL(r_model_part.GetSubModelPart("Outlet").NumberOfNodes(), 1);

        const auto& r_const_ref = r_model_part;
        auto r_smp_names = r_const_ref.GetSubModelPartNames();
        std::sort(r_smp_names.begin(), r_smp_names.end());
        const std::vector<std::string> r_smp_ref_names {"AnotherOutlet", "Inlet1", "Inlet2", "Outlet"};

        for (std::size_t i=0; i<r_smp_names.size(); ++i) {
            KRATOS_CHECK_EQUAL(r_smp_names[i], r_smp_ref_names[i]);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartAddNodalSolutionStepVariable, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        r_model_part.AddNodalSolutionStepVariable(VELOCITY);

        r_model_part.CreateNewNode(123, 0.00,0.00,0.00);

        KRATOS_CHECK_EXCEPTION_IS_THROWN(r_model_part.AddNodalSolutionStepVariable(PRESSURE),
            "Error: Attempting to add the variable \"PRESSURE\" to the model part with name \"Main\" which is not empty");

        r_model_part.AddNodalSolutionStepVariable(VELOCITY); // Adding the same Variable twice is fine bcs it wont do anything
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartHasNodalSolutionStepVariable, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");


        r_model_part.AddNodalSolutionStepVariable(VELOCITY);

        KRATOS_CHECK(r_model_part.HasNodalSolutionStepVariable(VELOCITY));
        KRATOS_CHECK_IS_FALSE(r_model_part.HasNodalSolutionStepVariable(PRESSURE));
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartFullName, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        ModelPart& r_sub_model_part = r_model_part.CreateSubModelPart("SubModelPart");
        ModelPart& r_sub_sub_model_part = r_sub_model_part.CreateSubModelPart("SubSubModelPart");

        KRATOS_CHECK_STRING_EQUAL(r_model_part.FullName(), "Main");
        KRATOS_CHECK_STRING_EQUAL(r_sub_model_part.FullName(), "Main.SubModelPart");
        KRATOS_CHECK_STRING_EQUAL(r_sub_sub_model_part.FullName(), "Main.SubModelPart.SubSubModelPart");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartEmptyName, KratosCoreFastSuite)
    {
        Model current_model;

        // Constructor with name
        KRATOS_CHECK_EXCEPTION_IS_THROWN(current_model.CreateModelPart(""),
            "Error: Please don't use empty names (\"\") when creating a ModelPart");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartBaseCreation, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Check results
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 6);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 4);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 4);
        KRATOS_CHECK(r_model_part.NumberOfGeometries() == 8);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartAddGeometries, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // Fill model part
        GenerateGenericModelPart(r_model_part);

        // Generate submodelpart (I)
        auto& r_sub_model_1 = r_model_part.CreateSubModelPart("Sub1");
        KRATOS_CHECK(r_sub_model_1.NumberOfGeometries() == 0);

        // Copy one
        r_sub_model_1.AddGeometries(r_model_part.GeometriesBegin(), r_model_part.GeometriesEnd());

        // Check results
        KRATOS_CHECK(r_sub_model_1.NumberOfGeometries() == 8);

        // Generate submodelpart (II)
        auto& r_sub_model_2 = r_model_part.CreateSubModelPart("Sub2");
        KRATOS_CHECK(r_sub_model_2.NumberOfGeometries() == 0);

        // Copy one
        const std::vector<std::size_t> indexes = {1,2,3,4};
        r_sub_model_2.AddGeometries(indexes);

        // Check results
        KRATOS_CHECK(r_sub_model_2.NumberOfGeometries() == 4);
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
        KRATOS_CHECK_DOUBLE_EQUAL(nearest, 4.0);
        KRATOS_CHECK_DOUBLE_EQUAL(r_table_model_part.GetValue(2.1), 4.2);
        KRATOS_CHECK_DOUBLE_EQUAL(r_table_model_part(2.1), 4.2);
        KRATOS_CHECK_DOUBLE_EQUAL(r_table_model_part.GetDerivative(2.1), 2.0);

        auto& r_data = r_table_model_part.Data();
        KRATOS_CHECK_EQUAL(r_data.size(), 6);

        ModelPart& r_sub_model_part = r_model_part.CreateSubModelPart("Sub");

        r_sub_model_part.AddTable(1, p_table);

        Table<double>& r_table_sub_model_part = r_sub_model_part.GetTable(1);

        nearest = (r_table_sub_model_part.GetNearestRow(2.1))[0];
        KRATOS_CHECK_DOUBLE_EQUAL(nearest, 4.0);
        KRATOS_CHECK_DOUBLE_EQUAL(r_table_sub_model_part.GetValue(2.1), 4.2);
        KRATOS_CHECK_DOUBLE_EQUAL(r_table_sub_model_part(2.1), 4.2);
        KRATOS_CHECK_DOUBLE_EQUAL(r_table_sub_model_part.GetDerivative(2.1), 2.0);

        r_data = r_table_sub_model_part.Data();
        KRATOS_CHECK_EQUAL(r_data.size(), 6);
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
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 6);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 2);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 4);
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
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 4);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 2);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 2);
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
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 6);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 2);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 4);
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
        KRATOS_CHECK(r_model_part.NumberOfNodes() == 3);
        KRATOS_CHECK(r_model_part.NumberOfConditions() == 2);
        KRATOS_CHECK(r_model_part.NumberOfElements() == 2);
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

        KRATOS_CHECK_EQUAL(r_model_part.NumberOfProperties(), 1);

        // Call method
        auto aux_util = AuxiliarModelPartUtilities(r_model_part);
        aux_util.EnsureModelPartOwnsProperties(false);

        // Check results
        KRATOS_CHECK_EQUAL(r_model_part.NumberOfProperties(), 2);

        aux_util.EnsureModelPartOwnsProperties(true);

        // Check results
        KRATOS_CHECK_EQUAL(r_model_part.NumberOfProperties(), 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartGetSubModelPart, KratosCoreFastSuite)
    {
        Model model;

        auto& model_part = model.CreateModelPart("Main");
        const auto& r_const_model_part = model_part;

        // Checking SubModelPart
        model_part.CreateSubModelPart("Inlet1");

        ModelPart& smp = model_part.GetSubModelPart("Inlet1");
        KRATOS_CHECK_EQUAL("Inlet1", smp.Name());

        const ModelPart& c_smp = r_const_model_part.GetSubModelPart("Inlet1");
        KRATOS_CHECK_EQUAL("Inlet1", c_smp.Name());

        KRATOS_CHECK_EXCEPTION_IS_THROWN(model_part.GetSubModelPart("Random"),
            "Error: There is no sub model part with name \"Random\" in model part \"Main\"\nThe following sub model parts are available:");

        KRATOS_CHECK_EXCEPTION_IS_THROWN(r_const_model_part.GetSubModelPart("Random"),
            "Error: There is no sub model part with name \"Random\" in model part \"Main\"\nThe following sub model parts are available:");

        // Checking SubSubModelPart
        auto& ssmp = smp.CreateSubModelPart("sub_inlet");

        KRATOS_CHECK(model_part.HasSubModelPart("Inlet1"));
        KRATOS_CHECK(model_part.HasSubModelPart("Inlet1.sub_inlet"));
        KRATOS_CHECK(r_const_model_part.HasSubModelPart("Inlet1"));
        KRATOS_CHECK(r_const_model_part.HasSubModelPart("Inlet1.sub_inlet"));

        KRATOS_CHECK_EQUAL("sub_inlet", model_part.GetSubModelPart("Inlet1.sub_inlet").Name());
        KRATOS_CHECK_EQUAL("sub_inlet", r_const_model_part.GetSubModelPart("Inlet1.sub_inlet").Name());

        KRATOS_CHECK_EXCEPTION_IS_THROWN(model_part.GetSubModelPart("Inlet1.random_sub_inlet"),
            "Error: There is no sub model part with name \"random_sub_inlet\" in model part \"Main.Inlet1\"\nThe following sub model parts are available:\n\t\"sub_inlet\"");

        KRATOS_CHECK_EXCEPTION_IS_THROWN(r_const_model_part.GetSubModelPart("Inlet1.random_sub_inlet"),
            "Error: There is no sub model part with name \"random_sub_inlet\" in model part \"Main.Inlet1\"\nThe following sub model parts are available:\n\t\"sub_inlet\"");

        // Checking SubSubSubModelPart
        ssmp.CreateSubModelPart("tiny_inlet");

        KRATOS_CHECK(model_part.HasSubModelPart("Inlet1.sub_inlet.tiny_inlet"));
        KRATOS_CHECK(r_const_model_part.HasSubModelPart("Inlet1.sub_inlet.tiny_inlet"));

        KRATOS_CHECK_EQUAL("tiny_inlet", model_part.GetSubModelPart("Inlet1.sub_inlet.tiny_inlet").Name());
        KRATOS_CHECK_EQUAL("tiny_inlet", r_const_model_part.GetSubModelPart("Inlet1.sub_inlet.tiny_inlet").Name());

        KRATOS_CHECK_EXCEPTION_IS_THROWN(model_part.GetSubModelPart("Inlet1.sub_inlet.big_inlet"),
            "Error: There is no sub model part with name \"big_inlet\" in model part \"Main.Inlet1.sub_inlet\"\nThe following sub model parts are available:\n\t\"tiny_inlet\"");

        KRATOS_CHECK_EXCEPTION_IS_THROWN(r_const_model_part.GetSubModelPart("Inlet1.sub_inlet.big_inlet"),
            "Error: There is no sub model part with name \"big_inlet\" in model part \"Main.Inlet1.sub_inlet\"\nThe following sub model parts are available:\n\t\"tiny_inlet\"");
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartHasSubModelPart, KratosCoreFastSuite)
    {
        Model model;

        auto& model_part = model.CreateModelPart("Main");

        model_part.CreateSubModelPart("Inlet1");

        KRATOS_CHECK(model_part.HasSubModelPart("Inlet1"));

        ModelPart& smp = model_part.GetSubModelPart("Inlet1");
        auto& ssmp = smp.CreateSubModelPart("SubInlet");
        ssmp.CreateSubModelPart("SuperSubInlet");

        KRATOS_CHECK(model_part.HasSubModelPart("Inlet1.SubInlet"));
        KRATOS_CHECK(model_part.HasSubModelPart("Inlet1.SubInlet.SuperSubInlet"));
        KRATOS_CHECK(smp.HasSubModelPart("SubInlet"));
        KRATOS_CHECK(smp.HasSubModelPart("SubInlet.SuperSubInlet"));
        KRATOS_CHECK(ssmp.HasSubModelPart("SuperSubInlet"));

        KRATOS_CHECK_IS_FALSE(model_part.HasSubModelPart("SubInlet"));
        KRATOS_CHECK_IS_FALSE(model_part.HasSubModelPart("SuperSubInlet"));
        KRATOS_CHECK_IS_FALSE(smp.HasSubModelPart("SuperSubInlet"));

        KRATOS_CHECK_IS_FALSE(model_part.HasSubModelPart("Random"));
    }

    KRATOS_TEST_CASE_IN_SUITE(ModelPartCreateSubModelPart, KratosCoreFastSuite)
    {
        Model model;

        auto& model_part = model.CreateModelPart("Main");

        // Checking SubModelPart
        auto& smp = model_part.CreateSubModelPart("Inlet1");
        KRATOS_CHECK_EQUAL("Inlet1", smp.Name());

        // Checking SubSubModelPart
        // here Inlet1 exists already
        auto& ssmp = model_part.CreateSubModelPart("Inlet1.sub_inlet");

        KRATOS_CHECK(model_part.HasSubModelPart("Inlet1"));
        KRATOS_CHECK(model_part.HasSubModelPart("Inlet1.sub_inlet"));
        KRATOS_CHECK(smp.HasSubModelPart("sub_inlet"));

        // here InletCustom does NOT exists yet
        auto& ssmp_2 = model_part.CreateSubModelPart("InletCustom.ccc_sub_inlet");

        KRATOS_CHECK(model_part.HasSubModelPart("InletCustom"));
        KRATOS_CHECK(model_part.HasSubModelPart("InletCustom.ccc_sub_inlet"));
        KRATOS_CHECK(model_part.GetSubModelPart("InletCustom").HasSubModelPart("ccc_sub_inlet"));
        KRATOS_CHECK_EQUAL(ssmp_2.Name(), "ccc_sub_inlet");

        // Checking SubSubSubModelPart
        auto& sssmp = model_part.CreateSubModelPart("Inlet1.sub_inlet.aabbcc");

        KRATOS_CHECK(model_part.HasSubModelPart("Inlet1"));
        KRATOS_CHECK(model_part.HasSubModelPart("Inlet1.sub_inlet"));
        KRATOS_CHECK(model_part.HasSubModelPart("Inlet1.sub_inlet.aabbcc"));
        KRATOS_CHECK(ssmp.HasSubModelPart("aabbcc"));
        KRATOS_CHECK_EQUAL(sssmp.Name(), "aabbcc");

        // here nothing exists yet
        auto& sssmp_2 = model_part.CreateSubModelPart("Fancy.xxx_sub_inlet.uztr");

        KRATOS_CHECK(model_part.HasSubModelPart("Fancy"));
        KRATOS_CHECK(model_part.HasSubModelPart("Fancy.xxx_sub_inlet"));
        KRATOS_CHECK(model_part.HasSubModelPart("Fancy.xxx_sub_inlet.uztr"));
        KRATOS_CHECK(model_part.GetSubModelPart("Fancy").HasSubModelPart("xxx_sub_inlet"));
        KRATOS_CHECK(model_part.GetSubModelPart("Fancy").HasSubModelPart("xxx_sub_inlet.uztr"));
        KRATOS_CHECK_EQUAL(sssmp_2.Name(), "uztr");
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
        KRATOS_CHECK_IS_FALSE(ssmp_2.HasSubModelPart("ccc_sub_inlet"));

        model_part.RemoveSubModelPart("Inlet1.sub_inlet.sub_sub_inlet");
        KRATOS_CHECK_IS_FALSE(smp.HasSubModelPart("sub_inlet.sub_sub_inlet"));
        KRATOS_CHECK(smp.HasSubModelPart("sub_inlet"));
        KRATOS_CHECK(model_part.HasSubModelPart("Inlet1.sub_inlet"));

        model_part.RemoveSubModelPart("test1");
        KRATOS_CHECK_IS_FALSE(model_part.HasSubModelPart("test1"));
        KRATOS_CHECK_IS_FALSE(model_part.HasSubModelPart("test1.InletCustom"));
        KRATOS_CHECK_IS_FALSE(model_part.HasSubModelPart("test1.InletCustom"));

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            model_part.RemoveSubModelPart("Inlet1.sub_inlet.sub_sub_inlet.test"),
            "Error: There is no sub model part with name \"sub_sub_inlet\" in model part \"Main.Inlet1.sub_inlet\"");
    }
  }  // namespace Testing.
}  // namespace Kratos.
