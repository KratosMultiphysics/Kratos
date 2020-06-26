//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/tetrahedra_3d_4.h"
#include "testing/testing.h"
#include "includes/model_part.h"

/* Utilities */
#include "utilities/cpp_tests_utilities.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_fast_point_locator_conditions.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3> NodeType;
        
        /** 
        * Checks the correct work of the binbased fast point locator
        * Test triangle 
        */
        KRATOS_TEST_CASE_IN_SUITE(BinBasedFastPointLocator1, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N");
            
            // We create the locator
            auto point_locator = BinBasedFastPointLocator<2>(this_model_part);
            point_locator.UpdateSearchDatabase();

            array_1d<double, 3> coordinates(3, 0.0);
            coordinates[0] = 0.5;
            coordinates[1] = 0.5;
            Vector shape_functions;
            Element::Pointer p_element;
            bool is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, 1000, 5.0e-2);

            Vector ref_shape_functions(3);
            ref_shape_functions[0] = 0.5;
            ref_shape_functions[1] = 0.0;
            ref_shape_functions[2] = 0.5;

            const double tolerance = 1.0e-16;
            KRATOS_CHECK(is_found);
            KRATOS_CHECK_EQUAL(p_element->Id(), this_model_part.pGetElement(1)->Id());
            KRATOS_CHECK_LESS_EQUAL(norm_2(shape_functions - ref_shape_functions), tolerance);

            coordinates[0] = -0.5;
            coordinates[1] = -0.5;
            is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, 1000, 5.0e-2);
            KRATOS_CHECK_IS_FALSE(is_found);
        }
        
        /** 
        * Checks the correct work of the binbased fast point locator
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(BinBasedFastPointLocator2, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            CppTestsUtilities::Create3DGeometry(this_model_part, "Element3D4N");
            
            // We create the locator
            auto point_locator = BinBasedFastPointLocator<3>(this_model_part);
            point_locator.UpdateSearchDatabase();

            array_1d<double, 3> coordinates(3, 0.0);
            coordinates[0] = 0.5;
            coordinates[1] = 0.5;
            coordinates[2] = 0.5;
            Vector shape_functions;
            Element::Pointer p_element;
            bool is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, 1000, 5.0e-2);

            Vector ref_shape_functions(4);
            ref_shape_functions[0] = 0.0;
            ref_shape_functions[1] = 0.5;
            ref_shape_functions[2] = 0.0;
            ref_shape_functions[3] = 0.5;

            const double tolerance = 1.0e-16;
            KRATOS_CHECK(is_found);
            KRATOS_CHECK_EQUAL(p_element->Id(), this_model_part.pGetElement(4)->Id());
            KRATOS_CHECK_LESS_EQUAL(norm_2(shape_functions - ref_shape_functions), tolerance);

            coordinates[0] = -0.5;
            coordinates[1] = -0.5;
            coordinates[2] = -0.5;
            is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, 1000, 5.0e-2);
            KRATOS_CHECK_IS_FALSE(is_found);
        }
        
        /** 
        * Checks the correct work of the binbased fast point locator
        * Test triangle for conditions
        */

        KRATOS_TEST_CASE_IN_SUITE(BinBasedFastPointLocator3, KratosCoreFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("test_model_part",2);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            CppTestsUtilities::Create2DGeometry(this_model_part, "SurfaceCondition3D3N", true, false);
            
            // We create the locator
            auto point_locator = BinBasedFastPointLocatorConditions<3>(this_model_part);
            point_locator.UpdateSearchDatabase();

            array_1d<double, 3> coordinates(3, 0.0);
            coordinates[0] = 0.5;
            coordinates[1] = 0.5;
            Vector shape_functions;
            Condition::Pointer p_condition;
            bool is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_condition, 1000, 5.0e-2);

            Vector ref_shape_functions(3);
            ref_shape_functions[0] = 0.5;
            ref_shape_functions[1] = 0.0;
            ref_shape_functions[2] = 0.5;

            const double tolerance = 1.0e-16;
            KRATOS_CHECK(is_found);
            KRATOS_CHECK_EQUAL(p_condition->Id(), this_model_part.pGetCondition(1)->Id());
            KRATOS_CHECK_LESS_EQUAL(norm_2(shape_functions - ref_shape_functions), tolerance);

            coordinates[0] = -0.5;
            coordinates[1] = -0.5;
            is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_condition, 1000, 5.0e-2);
            KRATOS_CHECK_IS_FALSE(is_found);
        }
    } // namespace Testing
}  // namespace Kratos.
