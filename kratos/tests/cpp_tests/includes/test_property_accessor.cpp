//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
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
#include"tests/cpp_tests/auxiliar_files_for_cpp_unnitest/test_element.h"
#include"tests/cpp_tests/auxiliar_files_for_cpp_unnitest/test_constitutive_law.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3> NodeType;

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
            KRATOS_CHECK_LESS_EQUAL(2.1e11 - initial_E, 1.0e-8);

            // custom accessor that returns 2.0
            class CustomAccessor
                : public  Accessor 
            {
                public:
                double GetProperty(
                    const Variable<double> &rVariable,
                    const Properties &rProperties,
                    const GeometryType &rGeometry,
                    const Vector &rShapeFunctionVector,
                    const ProcessInfo &rProcessInfo) override
                {
                    return mValue * rProperties[YOUNG_MODULUS];
                }

                Accessor::Pointer Clone() const override
                {
                    return Kratos::make_shared<CustomAccessor>(*this);
                }

                private:
                    const double mValue = 2.0;
            };

            auto& r_process_info = this_model_part.GetProcessInfo();
            TestConstitutiveLaw r_clone_cl = TestConstitutiveLaw();
            p_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

            auto p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            auto p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            auto p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            auto p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);

            std::vector<NodeType::Pointer> geom(4);
            geom[0] = p_node_1;
            geom[1] = p_node_2;
            geom[2] = p_node_3;
            geom[3] = p_node_4;
            auto pgeom = Kratos::make_shared<Quadrilateral2D4<NodeType>>(PointerVector<NodeType>{geom});

            auto p_elem = Kratos::make_intrusive<TestElement>(0, pgeom, p_prop, TestElement::ResidualType::LINEAR);
            p_elem->Initialize(r_process_info);
            this_model_part.AddElement(p_elem);


            CustomAccessor custom_accessor = CustomAccessor();
            p_prop->SetAccessor(YOUNG_MODULUS, custom_accessor.Clone());

            Vector N;
            const double modified_E = p_prop->GetValue(YOUNG_MODULUS, *pgeom, N, r_process_info);
            KRATOS_CHECK_LESS_EQUAL(2.1e11 * 2.0 - modified_E, 1.0e-8);
        }


    } // namespace Testing
}  // namespace Kratos.
