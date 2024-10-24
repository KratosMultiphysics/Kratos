//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//
//


// Project includes
#include "testing/testing.h"
#include "containers/model.h"


namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ModelGetModelPart, KratosCoreFastSuite)
{
    Model model;

    auto& model_part = model.CreateModelPart("Main");

    model_part.CreateSubModelPart("Inlet1");


    KRATOS_EXPECT_EQ(model.GetModelPart("Main").Name(), model_part.Name());

    ModelPart& smp = model_part.GetSubModelPart("Inlet1");
    KRATOS_EXPECT_EQ(model.GetModelPart("Main.Inlet1").Name(), smp.Name());

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(model.GetModelPart("Main.Random"),
        "Error: There is no sub model part with name \"Random\" in model part \"Main\"\nThe following sub model parts are available:");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(model.GetModelPart("Inlet1"),
        "Error: The ModelPart named : \"Inlet1\" was not found as root-ModelPart. The total input string was \"Inlet1\"");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        model.GetModelPart("Maiiiiin"),
        "Error: The ModelPart named : \"Maiiiiin\" was not found as root-ModelPart. The total input string was \"Maiiiiin\"");
}

KRATOS_TEST_CASE_IN_SUITE(ModelDataValueContainer, KratosCoreFastSuite)
{
    Model model;
    model.GetDataValueContainer().SetValue(DENSITY,1.0);
    KRATOS_EXPECT_TRUE(model.GetDataValueContainer().Has(DENSITY));
    KRATOS_EXPECT_FALSE(model.GetDataValueContainer().Has(TEMPERATURE));
    KRATOS_EXPECT_DOUBLE_EQ(model.GetDataValueContainer().GetValue(DENSITY),1.0);
}

KRATOS_TEST_CASE_IN_SUITE(ModelHasModelPart, KratosCoreFastSuite)
{
    Model model;

    auto& model_part = model.CreateModelPart("Main");

    model_part.CreateSubModelPart("Inlet1");


    KRATOS_EXPECT_TRUE(model.HasModelPart("Main"));
    KRATOS_EXPECT_TRUE(model.HasModelPart("Main.Inlet1"));

    KRATOS_EXPECT_FALSE(model.HasModelPart("Inlet1"));

    ModelPart& smp = model_part.GetSubModelPart("Inlet1");
    smp.CreateSubModelPart("SubInlet");

    KRATOS_EXPECT_TRUE(model.HasModelPart("Main.Inlet1.SubInlet"));
    KRATOS_EXPECT_FALSE(model.HasModelPart("Main.SubInlet"));

    KRATOS_EXPECT_FALSE(model.HasModelPart("Random"));
}

KRATOS_TEST_CASE_IN_SUITE(ModelDeleteModelPart, KratosCoreFastSuite)
{
    Model model;

    auto& model_part = model.CreateModelPart("Main");
    model_part.CreateSubModelPart("Inlet1");
    model_part.GetSubModelPart("Inlet1").CreateSubModelPart("SubSub");

    KRATOS_EXPECT_TRUE(model.HasModelPart("Main"));
    KRATOS_EXPECT_TRUE(model.HasModelPart("Main.Inlet1"));
    KRATOS_EXPECT_TRUE(model.HasModelPart("Main.Inlet1.SubSub"));

    model.DeleteModelPart("Main");

    KRATOS_EXPECT_FALSE(model.HasModelPart("Main"));
    KRATOS_EXPECT_FALSE(model.HasModelPart("Main.Inlet1"));
    KRATOS_EXPECT_FALSE(model.HasModelPart("Main.Inlet1.SubSub"));
}

KRATOS_TEST_CASE_IN_SUITE(ModelDeleteSubModelPart, KratosCoreFastSuite)
{
    Model model;

    auto& model_part = model.CreateModelPart("Main");
    model_part.CreateSubModelPart("Inlet1");

    KRATOS_EXPECT_TRUE(model.HasModelPart("Main"));
    KRATOS_EXPECT_TRUE(model.HasModelPart("Main.Inlet1"));

    model.DeleteModelPart("Main.Inlet1");

    KRATOS_EXPECT_TRUE(model.HasModelPart("Main"));
    KRATOS_EXPECT_FALSE(model.HasModelPart("Main.Inlet1"));
}

KRATOS_TEST_CASE_IN_SUITE(ModelGetModel, KratosCoreFastSuite)
{
    Model model;

    auto& model_part = model.CreateModelPart("Main");
    model_part.CreateSubModelPart("Inlet1");
    model_part.GetSubModelPart("Inlet1").CreateSubModelPart("SubSub");

    KRATOS_EXPECT_TRUE(&model == &model_part.GetModel());
}

KRATOS_TEST_CASE_IN_SUITE(ModelGetModelPartNames, KratosCoreFastSuite)
{
    Model model;

    auto& model_part = model.CreateModelPart("Main");
    model_part.CreateSubModelPart("Inlet1");
    model_part.GetSubModelPart("Inlet1").CreateSubModelPart("SubSub");

    std::vector<std::string> model_part_names = model.GetModelPartNames();

    KRATOS_EXPECT_TRUE(std::find(model_part_names.begin(), model_part_names.end(), "Main") != model_part_names.end());
    KRATOS_EXPECT_TRUE(std::find(model_part_names.begin(), model_part_names.end(), "Main.Inlet1") != model_part_names.end());
    KRATOS_EXPECT_TRUE(std::find(model_part_names.begin(), model_part_names.end(), "Main.Inlet1.SubSub") != model_part_names.end());
}

KRATOS_TEST_CASE_IN_SUITE(ModelCreateModelPart, KratosCoreFastSuite)
{
    Model model;

    auto& model_part = model.CreateModelPart("Main");
    KRATOS_EXPECT_EQ("Main", model_part.Name());

    // Checking SubModelPart
    auto& smp = model.CreateModelPart("Main.Inlet1");
    KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1"));
    KRATOS_EXPECT_EQ("Inlet1", smp.Name());

    // Checking SubSubModelPart
    // here Inlet1 exists already
    auto& ssmp = model.CreateModelPart("Main.Inlet1.sub_inlet");
    KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("Inlet1.sub_inlet"));
    KRATOS_EXPECT_EQ("sub_inlet", ssmp.Name());

    // here InletCustom does NOT exists yet
    auto& ssmp_2 = model.CreateModelPart("Main.InletCustom.ccc_sub_inlet");
    KRATOS_EXPECT_TRUE(model_part.HasSubModelPart("InletCustom.ccc_sub_inlet"));
    KRATOS_EXPECT_EQ("ccc_sub_inlet", ssmp_2.Name());

    // here nothing exists yet
    auto& smp_sub_outlet = model.CreateModelPart("Main2.Outlet.sub_outlet");
    auto& model_part_outlet = smp_sub_outlet.GetRootModelPart();
    KRATOS_EXPECT_EQ("Main2", model_part_outlet.Name());
    KRATOS_EXPECT_EQ("sub_outlet", smp_sub_outlet.Name());
    KRATOS_EXPECT_TRUE(model_part_outlet.HasSubModelPart("Outlet.sub_outlet"));
}

}   // namespace Testing
}  // namespace Kratos.
