// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "boost/numeric/ublas/assignment.hpp"
#include "custom_elements/membrane_elements/membrane_element.hpp"
#include "custom_strategies/schemes/load_stepping_scheme.hpp"
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;

class MockElementForLoadSteppingScheme : public Element
{
public:
    void Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rVariable == INTERNAL_FORCES_VECTOR) Output = mInternalForces;
        if (rVariable == EXTERNAL_FORCES_VECTOR) Output = mExternalForces;
    }

    void SetInternalForces(const Vector& rInternalForces) { mInternalForces = rInternalForces; }

    void SetExternalForces(const Vector& rExternalForces) { mExternalForces = rExternalForces; }

private:
    Vector mInternalForces;
    Vector mExternalForces;
};

KRATOS_TEST_CASE_IN_SUITE(LoadSteppingSchemeRHSAtStartOfStageIsEqualToCurrentInternalForcesAtStartMinusInternalForcesAtStart,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    auto element = Kratos::make_intrusive<MockElementForLoadSteppingScheme>();
    element->SetId(1);

    ProcessInfo CurrentProcessInfo;
    CurrentProcessInfo[TIME]       = 0.0;
    CurrentProcessInfo[START_TIME] = 0.0;
    CurrentProcessInfo[END_TIME]   = 1.0;
    std::vector<std::size_t> EquationId;

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddElement(element);
    CompressedMatrix A;
    Vector           Dx;
    Vector           b;

    Vector internal_forces_at_start_of_stage(4);
    internal_forces_at_start_of_stage <<= 1.0, 2.0, 3.0, 4.0; // Change later to new creation function
    Vector external_forces_at_start_of_stage = internal_forces_at_start_of_stage * 2;
    element->SetInternalForces(internal_forces_at_start_of_stage);
    element->SetExternalForces(external_forces_at_start_of_stage);
    scheme.InitializeSolutionStep(model_part, A, Dx, b);
    Vector actual_right_hand_side;

    Vector difference{ScalarVector(4, 0.1)};
    Vector current_internal_forces = internal_forces_at_start_of_stage + difference;
    element->SetInternalForces(current_internal_forces);
    scheme.CalculateRHSContribution(*element, actual_right_hand_side, EquationId, CurrentProcessInfo);

    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, difference, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(LoadSteppingSchemeRHSAtEndOfStageIsEqualToCurrentInternalForcesPlusExternalForcesAtStart,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    auto element = Kratos::make_intrusive<MockElementForLoadSteppingScheme>();
    element->SetId(1);

    ProcessInfo CurrentProcessInfo;
    CurrentProcessInfo[TIME]       = 1.0;
    CurrentProcessInfo[START_TIME] = 0.0;
    CurrentProcessInfo[END_TIME]   = 1.0;
    std::vector<std::size_t> EquationId;

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddElement(element);
    CompressedMatrix A;
    Vector           Dx;
    Vector           b;

    Vector internal_forces_at_start_of_stage(4);
    internal_forces_at_start_of_stage <<= 1.0, 2.0, 3.0, 4.0; // Change later to new creation function

    Vector external_forces_at_start_of_stage = internal_forces_at_start_of_stage * 2;
    element->SetInternalForces(internal_forces_at_start_of_stage);
    element->SetExternalForces(external_forces_at_start_of_stage);
    scheme.InitializeSolutionStep(model_part, A, Dx, b);
    Vector actual_right_hand_side;

    Vector difference{ScalarVector(4, 0.1)};
    Vector current_internal_forces = internal_forces_at_start_of_stage + difference;
    element->SetInternalForces(current_internal_forces);
    scheme.CalculateRHSContribution(*element, actual_right_hand_side, EquationId, CurrentProcessInfo);

    auto expected = current_internal_forces + external_forces_at_start_of_stage;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(LoadSteppingSchemeRHSAtHalfWayIsWeightedSum, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    auto element = Kratos::make_intrusive<MockElementForLoadSteppingScheme>();
    element->SetId(1);

    ProcessInfo CurrentProcessInfo;
    CurrentProcessInfo[TIME]       = 0.5;
    CurrentProcessInfo[START_TIME] = 0.0;
    CurrentProcessInfo[END_TIME]   = 1.0;
    std::vector<std::size_t> EquationId;

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddElement(element);
    CompressedMatrix A;
    Vector           Dx;
    Vector           b;

    Vector internal_forces_at_start_of_stage(4);
    internal_forces_at_start_of_stage <<= -1.0, -2.0, -3.0, -4.0; // Change later to new creation function

    Vector external_forces_at_start_of_stage(4);
    external_forces_at_start_of_stage <<= 5.0, 6.0, 7.0, 8.0;
    element->SetInternalForces(internal_forces_at_start_of_stage);
    element->SetExternalForces(external_forces_at_start_of_stage);
    scheme.InitializeSolutionStep(model_part, A, Dx, b);
    Vector actual_right_hand_side;

    Vector current_internal_forces(4);
    current_internal_forces <<= -2.0, -3.0, -4.0, -5.0;
    element->SetInternalForces(current_internal_forces);
    scheme.CalculateRHSContribution(*element, actual_right_hand_side, EquationId, CurrentProcessInfo);

    auto expected = Vector{ScalarVector(4, 1.0)};
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected, 1e-6);
}

} // namespace Kratos::Testing