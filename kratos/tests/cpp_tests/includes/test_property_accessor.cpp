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
#include "geometries/quadrilateral_2d_4.h"
#include "tests/cpp_tests/auxiliar_files_for_cpp_unnitest/test_element.h"
#include "tests/cpp_tests/auxiliar_files_for_cpp_unnitest/test_constitutive_law.h"

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
    KRATOS_CHECK_NEAR(2.1e11, initial_E, 1.0e-8);

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
    KRATOS_CHECK(p_prop->HasAccessor(YOUNG_MODULUS))

    Vector N;
    const double modified_E = p_prop->GetValue(YOUNG_MODULUS, *pgeom, N, r_process_info);
    KRATOS_CHECK_NEAR(2.1e11 * 2.0,  modified_E, 1.0e-8);

    const auto& r_accessor = p_prop->GetAccessor(YOUNG_MODULUS);
    const double modified_E_from_acc = r_accessor.GetValue(YOUNG_MODULUS, *p_prop, *pgeom, N, r_process_info);
    KRATOS_CHECK_NEAR(2.1e11 * 2.0, modified_E_from_acc, 1.0e-8);
}

/**
* Checks the correct work of the TableAccessor
*/
KRATOS_TEST_CASE_IN_SUITE(AAAATableAccessorSimpleProperties, KratosCoreFastSuite)
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
        N[0] = 1.0;
        N[0] = 2.0;
        N[0] = 3.0;
        N[0] = 4.0;

        KRATOS_CHECK_EQUAL(2.0e6, (*p_elem_prop)[YOUNG_MODULUS]);
        KRATOS_CHECK_EQUAL(2.0e6, (*p_elem_prop).GetValue(YOUNG_MODULUS));
        KRATOS_CHECK_EQUAL(2.0e6, (*p_elem_prop).GetValue(YOUNG_MODULUS, *p_geom, N, r_model_part.GetProcessInfo()));
        KRATOS_CHECK_EQUAL(false, (*p_elem_prop).HasAccessor(YOUNG_MODULUS));

        p_node_1->GetSolutionStepValue(TEMPERATURE) = 25.0;
        p_node_2->GetSolutionStepValue(TEMPERATURE) = 30.0;
        p_node_3->GetSolutionStepValue(TEMPERATURE) = 35.0;
        p_node_4->GetSolutionStepValue(TEMPERATURE) = 40.0;

        // we create a table relating TEMPERATURE with YOUNG_MODULUS

        // auto table_accessor = 

}

}  // namespace Kratos::Testing.
