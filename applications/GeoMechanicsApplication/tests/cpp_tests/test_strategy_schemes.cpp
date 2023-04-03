// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#pragma once

// Project includes
#include "testing/testing.h"

#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "custom_elements/U_Pw_small_strain_element.hpp"

// test.Initialize();
// test.InitializeSolutionStep();
// test.Check();
// test.Predict();
// test.InitializeNonLinIteration();
// test.FinalizeNonLinIteration();
// test.GetDofList();
// test.GetDofList();
// test.EquationId();
// test.EquationId();
// test.FinalizeSolutionStep();
// test.FinalizeSolutionStepActiveEntities();
// test.Update();
// test.CalculateSystemContributions();
// test.CalculateSystemContributions();
// test.CalculateRHSContribution();
// test.CalculateRHSContribution();
// test.CalculateLHSContribution();
// test.CalculateLHSContribution();

namespace Kratos
{
    namespace Testing
    {
        KRATOS_TEST_CASE_IN_SUITE(StrategySchemesInit, KratosGeoMechanicsFastSuite)
        {
            using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
            using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
            using TSystemMatrixType = NewmarkQuasistaticUPwScheme< SparseSpaceType, LocalSpaceType >::TSystemMatrixType;
            using TSystemVectorType = NewmarkQuasistaticUPwScheme< SparseSpaceType, LocalSpaceType >::TSystemVectorType;

            Kratos::VariablesList::Pointer variable_list;
            Kratos::Model model;
            Kratos::ModelPart& model_part = model.CreateModelPart("TestPart", 5);
            

            auto test = NewmarkQuasistaticUPwScheme< SparseSpaceType, LocalSpaceType >(1.0, 1.0, 1.0);

            // auto model_part = Kratos::ModelPart(variable_list, model);

            test.Initialize(model_part);


            TSystemMatrixType A;
            TSystemVectorType Dx;
            TSystemVectorType b;
            test.InitializeSolutionStep(model_part, A, Dx, b);

            KRATOS_CHECK_EQUAL(true, true);
        }
    
        KRATOS_TEST_CASE_IN_SUITE(StrategySchemesCheck, KratosGeoMechanicsFastSuite)
        {
            using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
            using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
            using TSystemMatrixType = NewmarkQuasistaticUPwScheme< SparseSpaceType, LocalSpaceType >::TSystemMatrixType;
            using TSystemVectorType = NewmarkQuasistaticUPwScheme< SparseSpaceType, LocalSpaceType >::TSystemVectorType;

            Kratos::VariablesList::Pointer variable_list;
            Kratos::Model model;
            Kratos::ModelPart& model_part = model.CreateModelPart("TestPart", 5);
        
            auto test = NewmarkQuasistaticUPwScheme< SparseSpaceType, LocalSpaceType >(1.0, 1.0, 1.0);

            auto node = Kratos::make_intrusive<Kratos::Node<3>>();
            model_part.AddNode(node);

            KRATOS_CHECK_EXCEPTION_IS_THROWN(test.Check(model_part), "DISPLACEMENT variable is not allocated");

            KRATOS_CHECK_EQUAL(true, true);
        }
    }
}