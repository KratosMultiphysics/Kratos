//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \ 
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/process_info.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ProcessInfoGetPreviousTimeStepInfo, KratosCoreFastSuite)
{
    //An empty history must report the intended errors for any number of steps asked,
    //including (and especially) StepsBefore > 1, which recursively walks the history
    //buffer: the null check must happen BEFORE the recursive call is made, otherwise
    //a null mpPreviousTimeStepInfo is dereferenced (UB/crash).
    ProcessInfo info;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(info.pGetPreviousTimeStepInfo(0), "Steps before = 0");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(info.pGetPreviousTimeStepInfo(1), "No previous time step exists.");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(info.pGetPreviousTimeStepInfo(2), "No previous time step exists.");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(info.pGetPreviousTimeStepInfo(3), "No previous time step exists.");

    //Same behaviour through the const overload
    const ProcessInfo& r_const_info = info;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_const_info.pGetPreviousTimeStepInfo(2), "No previous time step exists.");

    //A history of one time step: the stored step is returned and asking for more
    //steps still raises the error instead of crashing
    info.SetAsTimeStepInfo();
    info.CreateSolutionStepInfo(1);
    KRATOS_EXPECT_NE(nullptr, info.pGetPreviousTimeStepInfo(1));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(info.pGetPreviousTimeStepInfo(2), "No previous time step exists.");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(info.pGetPreviousTimeStepInfo(5), "No previous time step exists.");

    //A history of two time steps: valid requests return, overflows throw
    info.SetAsTimeStepInfo();
    info.CreateSolutionStepInfo(2);
    ProcessInfo::Pointer p_first = info.pGetPreviousTimeStepInfo(1);
    ProcessInfo::Pointer p_second = info.pGetPreviousTimeStepInfo(2);
    KRATOS_EXPECT_NE(nullptr, p_first);
    KRATOS_EXPECT_NE(nullptr, p_second);
    KRATOS_EXPECT_NE(p_first.get(), p_second.get());
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(info.pGetPreviousTimeStepInfo(3), "No previous time step exists.");
}

}   // namespace Testing
}  // namespace Kratos.
