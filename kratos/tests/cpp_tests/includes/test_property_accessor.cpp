//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/properties.h"
#include "includes/table_accessor.h"
#include "geometries/quadrilateral_2d_4.h"
#include "tests/test_utilities/test_element.h"
#include "tests/test_utilities/test_constitutive_law.h"
#include "includes/serializer.h"

namespace Kratos::Testing 
{

/**
* Checks the correct work of the Has methods
*/
KRATOS_TEST_CASE_IN_SUITE(PropertyAccessorSimpleProperties, KratosCoreFastSuite)
{
    Model current_model;
    auto &this_model_part = current_model.CreateModelPart("ModelPart",1);
    auto p_prop = this_model_part.CreateNewProperties(0);
    p_prop->SetValue(YOUNG_MODULUS, 2.1e11);

    const double initial_E = ((*p_prop)[YOUNG_MODULUS]);
    KRATOS_EXPECT_NEAR(2.1e11, initial_E, 1.0e-8);

    // custom accessor that returns 2.0
    class CustomAccessor
        : public  Accessor 
    {
        public:
        double GetValue(
            const Variable<double> &rVariable,
            const Properties &rProperties,
            const GeometryType &rGeometry,
            const Vector &rShapeFunctionVector,
            const ProcessInfo &rProcessInfo) const override
        {
            return mValue * rProperties[rVariable];
        }

        Accessor::UniquePointer Clone() const override
        {
            return Kratos::make_unique<CustomAccessor>(*this);
        }

        private:
            const double mValue = 2.0;
    };

    auto& r_process_info = this_model_part.GetProcessInfo();

    auto p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    auto p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);

    std::vector<Node::Pointer> geom(4);
    geom[0] = p_node_1;
    geom[1] = p_node_2;
    geom[2] = p_node_3;
    geom[3] = p_node_4;
    auto pgeom = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{geom});

    auto p_elem = Kratos::make_intrusive<TestElement>(0, pgeom, p_prop, TestElement::ResidualType::LINEAR);

    p_prop->SetAccessor(YOUNG_MODULUS, std::make_unique<CustomAccessor>());
    KRATOS_EXPECT_TRUE(p_prop->HasAccessor(YOUNG_MODULUS))

    Vector N;
    const double modified_E = p_prop->GetValue(YOUNG_MODULUS, *pgeom, N, r_process_info);
    KRATOS_EXPECT_NEAR(2.1e11 * 2.0,  modified_E, 1.0e-8);

    const auto& r_accessor = p_prop->GetAccessor(YOUNG_MODULUS);
    const double modified_E_from_acc = r_accessor.GetValue(YOUNG_MODULUS, *p_prop, *pgeom, N, r_process_info);
    KRATOS_EXPECT_NEAR(2.1e11 * 2.0, modified_E_from_acc, 1.0e-8);
}

/**
* Checks the correct work of the TableAccessor
*/
KRATOS_TEST_CASE_IN_SUITE(TableAccessorSimpleProperties, KratosCoreFastSuite)
{
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

        r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);

        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);

        std::vector<Node::Pointer> geom(4);
        geom[0] = p_node_1;
        geom[1] = p_node_2;
        geom[2] = p_node_3;
        geom[3] = p_node_4;

        for (auto& r_node : r_model_part.Nodes()){
            r_node.AddDof(TEMPERATURE);
        }

        auto p_geom = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{geom});
        Vector N = ZeroVector(4);
        N[0] = 0.1;
        N[0] = 0.2;
        N[0] = 0.3;
        N[0] = 0.4;

        KRATOS_EXPECT_EQ(2.0e6, (*p_elem_prop)[YOUNG_MODULUS]);
        KRATOS_EXPECT_EQ(2.0e6, (*p_elem_prop).GetValue(YOUNG_MODULUS));
        KRATOS_EXPECT_EQ(2.0e6, (*p_elem_prop).GetValue(YOUNG_MODULUS, *p_geom, N, r_model_part.GetProcessInfo()));
        KRATOS_EXPECT_EQ(false, (*p_elem_prop).HasAccessor(YOUNG_MODULUS));

        p_node_1->GetSolutionStepValue(TEMPERATURE) = 25.0;
        p_node_2->GetSolutionStepValue(TEMPERATURE) = 30.0;
        p_node_3->GetSolutionStepValue(TEMPERATURE) = 35.0;
        p_node_4->GetSolutionStepValue(TEMPERATURE) = 40.0;

        Table<double> T_E_table;
        T_E_table.PushBack(0.0,   2.0e6);
        T_E_table.PushBack(25.0,  1.0e6);
        T_E_table.PushBack(50.0,  0.5e6);
        T_E_table.PushBack(200.0, 0.25e6);

        p_elem_prop->SetTable(TEMPERATURE, YOUNG_MODULUS, T_E_table);
        KRATOS_EXPECT_EQ(true, (*p_elem_prop).HasTable(TEMPERATURE, YOUNG_MODULUS));

        TableAccessor E_table_accessor = TableAccessor(TEMPERATURE, "node_historical");
        p_elem_prop->SetAccessor(YOUNG_MODULUS, E_table_accessor.Clone());
        KRATOS_EXPECT_EQ(true, (*p_elem_prop).HasAccessor(YOUNG_MODULUS));
        KRATOS_EXPECT_EQ(1.6e6, (*p_elem_prop).GetValue(YOUNG_MODULUS, *p_geom, N, r_model_part.GetProcessInfo()));

        Table<double> T_NU_table;
        T_NU_table.PushBack(0.0,   0.3);
        T_NU_table.PushBack(25.0,  0.4);
        T_NU_table.PushBack(50.0,  0.41);
        T_NU_table.PushBack(200.0, 0.43);
        p_elem_prop->SetTable(TEMPERATURE, POISSON_RATIO, T_NU_table);
        TableAccessor nu_table_accessor = TableAccessor(TEMPERATURE); // using the default nodal_historical
        p_elem_prop->SetAccessor(POISSON_RATIO, nu_table_accessor.Clone());

        KRATOS_EXPECT_EQ(true, (*p_elem_prop).HasAccessor(POISSON_RATIO));
        KRATOS_EXPECT_EQ(true, (*p_elem_prop).HasAccessor(YOUNG_MODULUS));
        KRATOS_EXPECT_EQ(1.6e6, (*p_elem_prop).GetValue(YOUNG_MODULUS, *p_geom, N, r_model_part.GetProcessInfo()));
        KRATOS_EXPECT_EQ(0.34, (*p_elem_prop).GetValue(POISSON_RATIO, *p_geom, N, r_model_part.GetProcessInfo()));
}

/**
* Checks the correct work of the TableAccessor when using ProcessInfo variables (TIME)
*/
KRATOS_TEST_CASE_IN_SUITE(TableAccessorSimplePropertiesProcessInfo, KratosCoreFastSuite)
{
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
        r_model_part.GetProcessInfo().SetValue(TIME, 0.0);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0);

        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);

        std::vector<Node::Pointer> geom(4);
        geom[0] = p_node_1;
        geom[1] = p_node_2;
        geom[2] = p_node_3;
        geom[3] = p_node_4;

        auto p_geom = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{geom});
        Vector N = ZeroVector(4);
        N[0] = 0.1;
        N[0] = 0.2;
        N[0] = 0.3;
        N[0] = 0.4;

        KRATOS_EXPECT_EQ(2.0, (*p_elem_prop)[YOUNG_MODULUS]);
        KRATOS_EXPECT_EQ(2.0, (*p_elem_prop).GetValue(YOUNG_MODULUS));
        KRATOS_EXPECT_EQ(2.0, (*p_elem_prop).GetValue(YOUNG_MODULUS, *p_geom, N, r_model_part.GetProcessInfo()));
        KRATOS_EXPECT_EQ(false, (*p_elem_prop).HasAccessor(YOUNG_MODULUS));

        Table<double> Time_E_table;
        Time_E_table.PushBack(0.0,   2.0);
        Time_E_table.PushBack(1.0,   1.0);

        p_elem_prop->SetTable(TIME, YOUNG_MODULUS, Time_E_table);
        KRATOS_EXPECT_EQ(true, (*p_elem_prop).HasTable(TIME, YOUNG_MODULUS));

        TableAccessor E_table_accessor = TableAccessor(TIME, "process_info");
        p_elem_prop->SetAccessor(YOUNG_MODULUS, E_table_accessor.Clone());
        KRATOS_EXPECT_EQ(true, (*p_elem_prop).HasAccessor(YOUNG_MODULUS));
        KRATOS_EXPECT_EQ(2.0, (*p_elem_prop).GetValue(YOUNG_MODULUS, *p_geom, N, r_model_part.GetProcessInfo()));

        r_model_part.GetProcessInfo().SetValue(TIME, 0.5);
        KRATOS_EXPECT_EQ(1.5, (*p_elem_prop).GetValue(YOUNG_MODULUS, *p_geom, N, r_model_part.GetProcessInfo()));
}

KRATOS_TEST_CASE_IN_SUITE(TableTableAccessorSerialization, KratosCoreFastSuite)
{
    StreamSerializer serializer;
    TableAccessor table_accessor = TableAccessor(TEMPERATURE, "node_historical");

    serializer.save("table_accessor_info", table_accessor);

    TableAccessor table_accessor_loaded;

    // Variable<double> *p_var_loaded;
    serializer.load("table_accessor_info", table_accessor_loaded);

    KRATOS_EXPECT_EQ(TEMPERATURE.Key(), table_accessor_loaded.GetInputVariable().Key());
}

}  // namespace Kratos::Testing.
