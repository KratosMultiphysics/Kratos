//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License: BSD License

#include <array>
#include <cmath>

#include "containers/model.h"
#include "includes/cfd_variables.h"
#include "includes/stream_serializer.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_elements/d_vms.h"
#include "custom_elements/data_containers/time_integrated_qs_vms/time_integrated_qs_vms_data.h"
#include "custom_utilities/qsvms_gap_element_utility.h"
#include "utilities/gap_quadrature_utility.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos::Testing
{
namespace
{
ModelPart& CreateDVMSQ1Part(Model& rModel)
{
    auto& r_part = rModel.CreateModelPart("Fluid", 3);
    for (const auto* p_variable : {&VELOCITY, &MESH_VELOCITY, &ACCELERATION, &BODY_FORCE, &ADVPROJ}) {
        r_part.AddNodalSolutionStepVariable(*p_variable);
    }
    r_part.AddNodalSolutionStepVariable(PRESSURE);
    r_part.AddNodalSolutionStepVariable(DIVPROJ);
    r_part.AddNodalSolutionStepVariable(NODAL_AREA);
    auto& r_info = r_part.GetProcessInfo();
    r_info[DELTA_TIME] = 0.1;
    r_info[DYNAMIC_TAU] = 1.0;
    r_info[OSS_SWITCH] = 0;
    r_info[STEP] = 1;
    Vector bdf(3);
    bdf[0] = 15.0; bdf[1] = -20.0; bdf[2] = 5.0;
    r_info[BDF_COEFFICIENTS] = bdf;
    auto p_properties = r_part.CreateNewProperties(1);
    p_properties->SetValue(DENSITY, 1.0);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 0.1);
    p_properties->SetValue(CONSTITUTIVE_LAW, Kratos::make_shared<Newtonian2DLaw>());
    r_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    r_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    for (auto& r_node : r_part.Nodes()) {
        r_node.AddDof(VELOCITY_X);
        r_node.AddDof(VELOCITY_Y);
        r_node.AddDof(PRESSURE);
        r_node.FastGetSolutionStepValue(PRESSURE) = 0.2 * r_node.X();
        r_node.FastGetSolutionStepValue(BODY_FORCE)[0] = 0.4;
        for (unsigned int step = 0; step < 3; ++step) {
            auto& r_velocity = r_node.FastGetSolutionStepValue(VELOCITY, step);
            r_velocity[0] = 0.2 + 0.1 * r_node.Y() - 0.03 * step;
            r_velocity[1] = -0.1 * r_node.X() + 0.01 * step;
        }
    }
    return r_part;
}

void SetBDFAcceleration(ModelPart& rPart)
{
    const auto& r_bdf = rPart.GetProcessInfo()[BDF_COEFFICIENTS];
    for (auto& r_node : rPart.Nodes()) {
        auto& r_acceleration = r_node.FastGetSolutionStepValue(ACCELERATION);
        for (unsigned int d = 0; d < 2; ++d) {
            r_acceleration[d] = 0.0;
            for (unsigned int step = 0; step < 3; ++step) {
                r_acceleration[d] += r_bdf[step] * r_node.FastGetSolutionStepValue(VELOCITY, step)[d];
            }
        }
    }
}
}

KRATOS_TEST_CASE_IN_SUITE(DVMS2D4N_BDFAssembly, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& r_part = CreateDVMSQ1Part(model);
    auto p_standard = r_part.CreateNewElement("DVMS2D4N", 1, {1, 2, 3, 4}, r_part.pGetProperties(1));
    auto p_integrated = r_part.CreateNewElement("TimeIntegratedDVMS2D4N", 2, {1, 2, 3, 4}, r_part.pGetProperties(1));
    auto& r_info = r_part.GetProcessInfo();
    p_standard->Initialize(r_info);
    p_integrated->Initialize(r_info);
    KRATOS_EXPECT_EQ(p_integrated->GetGeometry().IntegrationPointsNumber(p_integrated->GetIntegrationMethod()), 4);
    KRATOS_EXPECT_EQ(p_integrated->Check(r_info), 0);
    for (unsigned int order = 1; order <= 2; ++order) {
        auto& r_bdf = r_info[BDF_COEFFICIENTS];
        r_bdf[0] = order == 1 ? 10.0 : 15.0;
        r_bdf[1] = order == 1 ? -10.0 : -20.0;
        r_bdf[2] = order == 1 ? 0.0 : 5.0;
        SetBDFAcceleration(r_part);
        p_standard->InitializeNonLinearIteration(r_info);
        p_integrated->InitializeNonLinearIteration(r_info);
        Matrix velocity_lhs, mass, lhs, lhs_only;
        Vector velocity_rhs, rhs, rhs_only;
        p_standard->CalculateLocalVelocityContribution(velocity_lhs, velocity_rhs, r_info);
        p_standard->CalculateMassMatrix(mass, r_info);
        p_integrated->CalculateLocalSystem(lhs, rhs, r_info);
        p_integrated->CalculateLeftHandSide(lhs_only, r_info);
        p_integrated->CalculateRightHandSide(rhs_only, r_info);
        KRATOS_EXPECT_EQ(lhs.size1(), 12);
        KRATOS_EXPECT_EQ(rhs.size(), 12);
        Matrix expected_lhs = velocity_lhs + r_bdf[0] * mass;
        Vector acceleration = ZeroVector(12);
        for (unsigned int i = 0; i < 4; ++i) {
            const auto& r_acceleration = p_standard->GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);
            acceleration[3*i] = r_acceleration[0];
            acceleration[3*i+1] = r_acceleration[1];
        }
        Vector expected_rhs = velocity_rhs - prod(mass, acceleration);
        KRATOS_EXPECT_MATRIX_NEAR(lhs, expected_lhs, 1e-10);
        KRATOS_EXPECT_VECTOR_NEAR(rhs, expected_rhs, 1e-10);
        KRATOS_EXPECT_MATRIX_NEAR(lhs_only, lhs, 1e-10);
        KRATOS_EXPECT_VECTOR_NEAR(rhs_only, rhs, 1e-10);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DVMS2D4N_SubscaleCommitAndRestart, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& r_part = CreateDVMSQ1Part(model);
    auto p_element = r_part.CreateNewElement("TimeIntegratedDVMS2D4N", 1, {1, 2, 3, 4}, r_part.pGetProperties(1));
    auto& r_info = r_part.GetProcessInfo();
    p_element->Initialize(r_info);
    p_element->InitializeNonLinearIteration(r_info);
    Matrix before_lhs, after_lhs;
    Vector before, after;
    p_element->CalculateLocalSystem(before_lhs, before, r_info);
    p_element->InitializeNonLinearIteration(r_info);
    p_element->CalculateLocalSystem(after_lhs, after, r_info);
    KRATOS_EXPECT_VECTOR_NEAR(before, after, 1e-10);

    p_element->FinalizeSolutionStep(r_info);
    p_element->CalculateLocalSystem(after_lhs, after, r_info);
    double history_effect = 0.0;
    for (std::size_t i = 0; i < before.size(); ++i) {
        history_effect += std::abs(after[i] - before[i]);
    }
    KRATOS_EXPECT_GT(history_effect, 1e-12);
    before = after;
    p_element->FinalizeSolutionStep(r_info);
    p_element->CalculateLocalSystem(after_lhs, after, r_info);
    KRATOS_EXPECT_VECTOR_NEAR(before, after, 1e-12);

    StreamSerializer serializer;
    serializer.save("Model", model);
    Model restored_model;
    serializer.load("Model", restored_model);
    auto& r_restored = restored_model.GetModelPart("Fluid");
    auto& r_restored_element = r_restored.GetElement(1);
    r_restored_element.Initialize(r_restored.GetProcessInfo());
    r_restored_element.FinalizeSolutionStep(r_restored.GetProcessInfo());
    p_element->InitializeNonLinearIteration(r_info);
    r_restored_element.InitializeNonLinearIteration(r_restored.GetProcessInfo());
    p_element->CalculateLocalSystem(before_lhs, before, r_info);
    r_restored_element.CalculateLocalSystem(after_lhs, after, r_restored.GetProcessInfo());
    KRATOS_EXPECT_VECTOR_NEAR(before, after, 1e-10);

    // Advancing the private predictor does not commit the accepted element.
    r_restored.GetProcessInfo()[STEP] = 2;
    r_restored_element.FinalizeSolutionStep(r_restored.GetProcessInfo());
    p_element->CalculateLocalSystem(after_lhs, after, r_info);
    KRATOS_EXPECT_VECTOR_NEAR(before, after, 1e-12);

    // Reload an independent copy for the uninterrupted/restarted comparison.
    StreamSerializer restart_serializer;
    restart_serializer.save("Model", model);
    Model restart_model;
    restart_serializer.load("Model", restart_model);
    auto& r_restart_part = restart_model.GetModelPart("Fluid");
    auto& r_restart_element = r_restart_part.GetElement(1);
    r_restart_element.Initialize(r_restart_part.GetProcessInfo());
    r_info[STEP] = 2;
    r_restart_part.GetProcessInfo()[STEP] = 2;
    p_element->InitializeNonLinearIteration(r_info);
    r_restart_element.InitializeNonLinearIteration(r_restart_part.GetProcessInfo());
    p_element->FinalizeSolutionStep(r_info);
    r_restart_element.FinalizeSolutionStep(r_restart_part.GetProcessInfo());
    p_element->CalculateLocalSystem(before_lhs, before, r_info);
    r_restart_element.CalculateLocalSystem(after_lhs, after, r_restart_part.GetProcessInfo());
    KRATOS_EXPECT_VECTOR_NEAR(before, after, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(DVMS2D4N_UniformTransientManufacturedSolution, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& r_part = CreateDVMSQ1Part(model);
    auto p_element = r_part.CreateNewElement("TimeIntegratedDVMS2D4N", 1, {1, 2, 3, 4}, r_part.pGetProperties(1));
    auto& r_info = r_part.GetProcessInfo();
    p_element->Initialize(r_info);
    for (unsigned int order = 1; order <= 2; ++order) {
        auto& r_bdf = r_info[BDF_COEFFICIENTS];
        r_bdf[0] = order == 1 ? 10.0 : 15.0;
        r_bdf[1] = order == 1 ? -10.0 : -20.0;
        r_bdf[2] = order == 1 ? 0.0 : 5.0;
        for (auto& r_node : r_part.Nodes()) {
            r_node.FastGetSolutionStepValue(PRESSURE) = 0.0;
            for (unsigned int step = 0; step < 3; ++step) {
                auto& r_velocity = r_node.FastGetSolutionStepValue(VELOCITY, step);
                r_velocity[0] = 0.4 - 0.02 * step;
                r_velocity[1] = 0.0;
            }
        }
        SetBDFAcceleration(r_part);
        for (auto& r_node : r_part.Nodes()) {
            r_node.FastGetSolutionStepValue(BODY_FORCE) = r_node.FastGetSolutionStepValue(ACCELERATION);
            // The element must calculate BDF acceleration independently of this field.
            r_node.FastGetSolutionStepValue(ACCELERATION).clear();
        }
        p_element->InitializeNonLinearIteration(r_info);
        Matrix lhs;
        Vector rhs;
        p_element->CalculateLocalSystem(lhs, rhs, r_info);
        for (const double value : rhs) {
            KRATOS_EXPECT_NEAR(value, 0.0, 1e-10);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DVMS2D4N_LinearVelocityPressurePatch, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& r_part = CreateDVMSQ1Part(model);
    r_part.CreateNewNode(5, 0.5, 0.0, 0.0);
    r_part.CreateNewNode(6, 1.0, 0.5, 0.0);
    r_part.CreateNewNode(7, 0.5, 1.0, 0.0);
    r_part.CreateNewNode(8, 0.0, 0.5, 0.0);
    r_part.CreateNewNode(9, 0.5, 0.5, 0.0);
    for (auto& r_node : r_part.Nodes()) {
        if (r_node.Id() > 4) {
            r_node.AddDof(VELOCITY_X);
            r_node.AddDof(VELOCITY_Y);
            r_node.AddDof(PRESSURE);
        }
        r_node.FastGetSolutionStepValue(PRESSURE) = r_node.X();
        r_node.FastGetSolutionStepValue(BODY_FORCE)[0] = 1.0;
        r_node.FastGetSolutionStepValue(BODY_FORCE)[1] = 0.0;
        r_node.FastGetSolutionStepValue(ACCELERATION).clear();
        for (unsigned int step = 0; step < 3; ++step) {
            r_node.FastGetSolutionStepValue(VELOCITY, step).clear();
            r_node.FastGetSolutionStepValue(VELOCITY, step)[0] = r_node.Y();
        }
    }
    const std::array<std::array<ModelPart::IndexType, 4>, 4> cells{{
        {{1, 5, 9, 8}}, {{5, 2, 6, 9}}, {{9, 6, 3, 7}}, {{8, 9, 7, 4}}}};
    auto& r_info = r_part.GetProcessInfo();
    double center_x = 0.0;
    double center_y = 0.0;
    double center_pressure = 0.0;
    for (std::size_t e = 0; e < cells.size(); ++e) {
        const auto& r_cell = cells[e];
        auto p_element = r_part.CreateNewElement("TimeIntegratedDVMS2D4N", e + 1,
            {r_cell[0], r_cell[1], r_cell[2], r_cell[3]}, r_part.pGetProperties(1));
        p_element->Initialize(r_info);
        p_element->InitializeNonLinearIteration(r_info);
        Matrix lhs;
        Vector rhs;
        p_element->CalculateLocalSystem(lhs, rhs, r_info);
        for (unsigned int i = 0; i < 4; ++i) {
            if (r_cell[i] == 9) {
                center_x += rhs[3*i];
                center_y += rhs[3*i+1];
                center_pressure += rhs[3*i+2];
            }
        }
    }
    KRATOS_EXPECT_NEAR(center_x, 0.0, 1e-10);
    KRATOS_EXPECT_NEAR(center_y, 0.0, 1e-10);
    KRATOS_EXPECT_NEAR(center_pressure, 0.0, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(DVMS2D4N_InvalidTimeData, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& r_part = CreateDVMSQ1Part(model);
    auto p_element = r_part.CreateNewElement("TimeIntegratedDVMS2D4N", 1, {1, 2, 3, 4}, r_part.pGetProperties(1));
    auto& r_info = r_part.GetProcessInfo();
    r_info[DELTA_TIME] = 0.0;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(r_info), "positive DELTA_TIME");
    r_info[DELTA_TIME] = 0.1;
    r_info[BDF_COEFFICIENTS] = ZeroVector(2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Check(r_info), "three BDF coefficients");
}

KRATOS_TEST_CASE_IN_SUITE(DVMS2D4N_GapQuadratureAndContinuitySource, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& r_part = CreateDVMSQ1Part(model);
    auto p_parent = r_part.CreateNewElement("TimeIntegratedDVMS2D4N", 1, {1, 2, 3, 4}, r_part.pGetProperties(1));
    auto& r_info = r_part.GetProcessInfo();
    p_parent->Initialize(r_info);
    const double a = 0.5 - 0.5 / std::sqrt(3.0);
    const double b = 0.5 + 0.5 / std::sqrt(3.0);
    Geometry<Node>::IntegrationPointsArrayType points;
    points.emplace_back(a, a, 0.0, 0.25);
    points.emplace_back(b, a, 0.0, 0.25);
    points.emplace_back(b, b, 0.0, 0.25);
    points.emplace_back(a, b, 0.0, 0.25);
    auto p_geometry = GapQuadratureUtility::Create(p_parent->GetGeometry(), points);
    auto p_gap = QSVMSGapElementUtility::Create(2, *p_parent, p_geometry);
    using DVMSQ1 = DVMS<TimeIntegratedQSVMSData<2,4>>;
    KRATOS_EXPECT_TRUE(dynamic_cast<DVMSQ1*>(p_gap.get()) != nullptr);
    p_gap->Initialize(r_info);
    Matrix parent_lhs, gap_lhs;
    Vector parent_rhs, gap_rhs, source_rhs;
    p_parent->CalculateLocalSystem(parent_lhs, parent_rhs, r_info);
    p_gap->CalculateLocalSystem(gap_lhs, gap_rhs, r_info);
    KRATOS_EXPECT_MATRIX_NEAR(parent_lhs, gap_lhs, 1e-10);
    KRATOS_EXPECT_VECTOR_NEAR(parent_rhs, gap_rhs, 1e-10);
    p_gap->SetValue(CONTINUITY_SOURCE, 2.0);
    p_gap->CalculateLocalSystem(gap_lhs, source_rhs, r_info);
    double continuity = 0.0;
    double stabilization = 0.0;
    for (unsigned int i = 0; i < 4; ++i) {
        continuity += source_rhs[3*i+2] - gap_rhs[3*i+2];
        for (unsigned int d = 0; d < 2; ++d) {
            stabilization += std::abs(source_rhs[3*i+d] - gap_rhs[3*i+d]);
        }
    }
    KRATOS_EXPECT_NEAR(continuity, 2.0, 1e-10);
    KRATOS_EXPECT_GT(stabilization, 1e-12);
    auto* p_dvms = dynamic_cast<DVMSQ1*>(p_gap.get());
    const auto source_contribution = p_dvms->CalculateContinuitySourceRHS(r_info);
    for (unsigned int i = 0; i < 12; ++i) {
        KRATOS_EXPECT_NEAR(source_contribution[i], source_rhs[i] - gap_rhs[i], 1e-10);
    }
    Matrix previous_source = ZeroMatrix(4, 3);
    for (unsigned int g = 0; g < 4; ++g) previous_source(g, 2) = 0.75;
    p_dvms->SetAcceptedSubscaleState(previous_source, -1);
    p_gap->CalculateLocalSystem(gap_lhs, source_rhs, r_info);
    const auto with_previous = p_dvms->CalculateContinuitySourceRHS(r_info);
    p_gap->SetValue(CONTINUITY_SOURCE, 0.0);
    previous_source.clear();
    p_dvms->SetAcceptedSubscaleState(previous_source, -1);
    p_gap->CalculateLocalSystem(gap_lhs, gap_rhs, r_info);
    for (unsigned int i = 0; i < 12; ++i) {
        KRATOS_EXPECT_NEAR(with_previous[i], source_rhs[i] - gap_rhs[i], 1e-10);
    }
    // Restore the current positive source for the sign check below.
    p_gap->SetValue(CONTINUITY_SOURCE, 2.0);
    p_gap->CalculateLocalSystem(gap_lhs, source_rhs, r_info);
    p_gap->SetValue(CONTINUITY_SOURCE, -2.0);
    Vector negative_rhs;
    p_gap->CalculateLocalSystem(gap_lhs, negative_rhs, r_info);
    for (unsigned int i = 0; i < 12; ++i) {
        KRATOS_EXPECT_NEAR(source_rhs[i] - gap_rhs[i], gap_rhs[i] - negative_rhs[i], 1e-10);
    }
}

} // namespace Kratos::Testing
