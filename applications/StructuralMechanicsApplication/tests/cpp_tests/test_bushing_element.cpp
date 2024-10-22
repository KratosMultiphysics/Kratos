// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//


// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/table.h"
#include "includes/variables.h"
#include "structural_mechanics_application_variables.h"
#include "structural_mechanics_fast_suite.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {
namespace {

void CreateBushingElementModelPart(
    ModelPart& rModelPart,
    const bool IsLinearStiffnessUsed)
{
    KRATOS_TRY;
    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    r_process_info[DOMAIN_SIZE] = 3;
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(ROTATION);

    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.5, 0.0, 0.0);

    for (auto& r_node : rModelPart.Nodes()) {
        auto& disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        disp[0] =  r_node.Id() / 4.0;
        disp[1] = -r_node.Id() / 4.0;
        disp[2] =  r_node.Id() * 0.1;

        auto& rot = r_node.FastGetSolutionStepValue(ROTATION);
        rot[0] = -r_node.Id() * 0.3 / 4.0;
        rot[1] = -r_node.Id() * 0.5 / 4.0;
        rot[2] =  r_node.Id() * 0.7 / 4.0;
    }

    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
        r_node.AddDof(ROTATION_X);
        r_node.AddDof(ROTATION_Y);
        r_node.AddDof(ROTATION_Z);
    }

    auto p_prop = rModelPart.CreateNewProperties(1);
    std::vector<IndexType> node_ids({1, 2});
    rModelPart.CreateNewElement("BushingElement3D2N", 1, node_ids, p_prop);

    if (IsLinearStiffnessUsed) {
        // constant stiffness based
        (*p_prop)[NODAL_DISPLACEMENT_STIFFNESS] = array_1d<double, 3>({572.0, 2960.0, 700.0});
        (*p_prop)[NODAL_ROTATIONAL_STIFFNESS] = array_1d<double, 3>({1.908e6, 3.223e5, 5.873e5});
    } else {
        // tabular data based
        Table<double> table_x_fx;
        auto& data_x_fx = table_x_fx.Data();
        data_x_fx.push_back(std::make_pair<double, std::array<double, 1UL>>(-1.0, std::array<double, 1UL>({-2e+4})));
        data_x_fx.push_back(std::make_pair<double, std::array<double, 1UL>>(-0.75,std::array<double, 1UL>({-1e+3})));
        data_x_fx.push_back(std::make_pair<double, std::array<double, 1UL>>(-0.25,std::array<double, 1UL>({-1e+1})));
        data_x_fx.push_back(std::make_pair<double, std::array<double, 1UL>>(0.00,std::array<double, 1UL>({0.0})));
        data_x_fx.push_back(std::make_pair<double, std::array<double, 1UL>>(0.25,std::array<double, 1UL>({1e+1})));
        data_x_fx.push_back(std::make_pair<double, std::array<double, 1UL>>(0.75,std::array<double, 1UL>({1e+3})));
        data_x_fx.push_back(std::make_pair<double, std::array<double, 1UL>>(1.0,std::array<double, 1UL>({2e+4})));
        p_prop->SetTable(DISPLACEMENT_X, FORCE_X, table_x_fx);

        Table<double> table_y_fy;
        auto& data_y_fy = table_y_fy.Data();
        data_y_fy.push_back(std::make_pair<double, std::array<double, 1UL>>(-1.0, std::array<double, 1UL>({-2e+4 * 2})));
        data_y_fy.push_back(std::make_pair<double, std::array<double, 1UL>>(-0.75,std::array<double, 1UL>({-1e+3 * 2})));
        data_y_fy.push_back(std::make_pair<double, std::array<double, 1UL>>(-0.25,std::array<double, 1UL>({-1e+1 * 2})));
        data_y_fy.push_back(std::make_pair<double, std::array<double, 1UL>>(0.00,std::array<double, 1UL>({0.0 * 2})));
        data_y_fy.push_back(std::make_pair<double, std::array<double, 1UL>>(0.25,std::array<double, 1UL>({1e+1 * 2})));
        data_y_fy.push_back(std::make_pair<double, std::array<double, 1UL>>(0.75,std::array<double, 1UL>({1e+3 * 2})));
        data_y_fy.push_back(std::make_pair<double, std::array<double, 1UL>>(1.0,std::array<double, 1UL>({2e+4 * 2})));
        p_prop->SetTable(DISPLACEMENT_Y, FORCE_Y, table_y_fy);

        Table<double> table_z_fz;
        auto& data_z_fz = table_z_fz.Data();
        data_z_fz.push_back(std::make_pair<double, std::array<double, 1UL>>(-1.0, std::array<double, 1UL>({-2e+4 * 1.5})));
        data_z_fz.push_back(std::make_pair<double, std::array<double, 1UL>>(-0.75,std::array<double, 1UL>({-1e+3 * 1.5})));
        data_z_fz.push_back(std::make_pair<double, std::array<double, 1UL>>(-0.25,std::array<double, 1UL>({-1e+1 * 1.5})));
        data_z_fz.push_back(std::make_pair<double, std::array<double, 1UL>>(0.00,std::array<double, 1UL>({0.0 * 1.5})));
        data_z_fz.push_back(std::make_pair<double, std::array<double, 1UL>>(0.25,std::array<double, 1UL>({1e+1 * 1.5})));
        data_z_fz.push_back(std::make_pair<double, std::array<double, 1UL>>(0.75,std::array<double, 1UL>({1e+3 * 1.5})));
        data_z_fz.push_back(std::make_pair<double, std::array<double, 1UL>>(1.0,std::array<double, 1UL>({2e+4 * 1.5})));
        p_prop->SetTable(DISPLACEMENT_Z, FORCE_Z, table_z_fz);

        Table<double> table_rx_mx;
        auto& data_rx_mx = table_rx_mx.Data();
        data_rx_mx.push_back(std::make_pair<double, std::array<double, 1UL>>(-1.0, std::array<double, 1UL>({-2e+4 * 0.7})));
        data_rx_mx.push_back(std::make_pair<double, std::array<double, 1UL>>(-0.75,std::array<double, 1UL>({-1e+3 * 0.7})));
        data_rx_mx.push_back(std::make_pair<double, std::array<double, 1UL>>(-0.25,std::array<double, 1UL>({-1e+1 * 0.7})));
        data_rx_mx.push_back(std::make_pair<double, std::array<double, 1UL>>(0.00,std::array<double, 1UL>({0.0 * 0.7})));
        data_rx_mx.push_back(std::make_pair<double, std::array<double, 1UL>>(0.25,std::array<double, 1UL>({1e+1 * 0.7})));
        data_rx_mx.push_back(std::make_pair<double, std::array<double, 1UL>>(0.75,std::array<double, 1UL>({1e+3 * 0.7})));
        data_rx_mx.push_back(std::make_pair<double, std::array<double, 1UL>>(1.0,std::array<double, 1UL>({2e+4 * 0.7})));
        p_prop->SetTable(ROTATION_X, MOMENT_X, table_rx_mx);

        Table<double> table_ry_my;
        auto& data_ry_my = table_ry_my.Data();
        data_ry_my.push_back(std::make_pair<double, std::array<double, 1UL>>(-1.0, std::array<double, 1UL>({-2e+4 * 2 * 0.7})));
        data_ry_my.push_back(std::make_pair<double, std::array<double, 1UL>>(-0.75,std::array<double, 1UL>({-1e+3 * 2 * 0.7})));
        data_ry_my.push_back(std::make_pair<double, std::array<double, 1UL>>(-0.25,std::array<double, 1UL>({-1e+1 * 2 * 0.7})));
        data_ry_my.push_back(std::make_pair<double, std::array<double, 1UL>>(0.00,std::array<double, 1UL>({0.0 * 2 * 0.7})));
        data_ry_my.push_back(std::make_pair<double, std::array<double, 1UL>>(0.25,std::array<double, 1UL>({1e+1 * 2 * 0.7})));
        data_ry_my.push_back(std::make_pair<double, std::array<double, 1UL>>(0.75,std::array<double, 1UL>({1e+3 * 2 * 0.7})));
        data_ry_my.push_back(std::make_pair<double, std::array<double, 1UL>>(1.0,std::array<double, 1UL>({2e+4 * 2 * 0.7})));
        p_prop->SetTable(ROTATION_Y, MOMENT_Y, table_ry_my);

        Table<double> table_rz_mz;
        auto& data_rz_mz = table_rz_mz.Data();
        data_rz_mz.push_back(std::make_pair<double, std::array<double, 1UL>>(-1.0, std::array<double, 1UL>({-2e+4 * 1.5 * 0.7})));
        data_rz_mz.push_back(std::make_pair<double, std::array<double, 1UL>>(-0.75,std::array<double, 1UL>({-1e+3 * 1.5 * 0.7})));
        data_rz_mz.push_back(std::make_pair<double, std::array<double, 1UL>>(-0.25,std::array<double, 1UL>({-1e+1 * 1.5 * 0.7})));
        data_rz_mz.push_back(std::make_pair<double, std::array<double, 1UL>>(0.00,std::array<double, 1UL>({0.0 * 1.5 * 0.7})));
        data_rz_mz.push_back(std::make_pair<double, std::array<double, 1UL>>(0.25,std::array<double, 1UL>({1e+1 * 1.5 * 0.7})));
        data_rz_mz.push_back(std::make_pair<double, std::array<double, 1UL>>(0.75,std::array<double, 1UL>({1e+3 * 1.5 * 0.7})));
        data_rz_mz.push_back(std::make_pair<double, std::array<double, 1UL>>(1.0,std::array<double, 1UL>({2e+4 * 1.5 * 0.7})));
        p_prop->SetTable(ROTATION_Z, MOMENT_Z, table_rz_mz);
    }

    rModelPart.GetElement(1).Initialize(r_process_info);
    rModelPart.GetElement(1).Check(r_process_info);

    KRATOS_CATCH("CreateBushingElementModelPart");
}

} // anonymous namespace


KRATOS_TEST_CASE_IN_SUITE(BushingElement_LHS_Constant, KratosStructuralMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("test");
    CreateBushingElementModelPart(r_model_part, true);

    Matrix lhs;
    r_model_part.GetElement(1).CalculateLeftHandSide(lhs, r_model_part.GetProcessInfo());

    Matrix ref_lhs(12, 12, 0.0);
    ref_lhs(0, 0) = 5.720000000000e+02;
    ref_lhs(0, 6) = -5.720000000000e+02;
    ref_lhs(1, 1) = 2.960000000000e+03;
    ref_lhs(1, 5) = 2.220000000000e+03;
    ref_lhs(1, 7) = -2.960000000000e+03;
    ref_lhs(1, 11) = 2.220000000000e+03;
    ref_lhs(2, 2) = 7.000000000000e+02;
    ref_lhs(2, 4) = -5.250000000000e+02;
    ref_lhs(2, 8) = -7.000000000000e+02;
    ref_lhs(2, 10) = -5.250000000000e+02;
    ref_lhs(3, 3) = 1.908000000000e+06;
    ref_lhs(3, 9) = -1.908000000000e+06;
    ref_lhs(4, 2) = -5.250000000000e+02;
    ref_lhs(4, 4) = 3.226937500000e+05;
    ref_lhs(4, 8) = 5.250000000000e+02;
    ref_lhs(4, 10) = -3.219062500000e+05;
    ref_lhs(5, 1) = 2.220000000000e+03;
    ref_lhs(5, 5) = 5.889650000000e+05;
    ref_lhs(5, 7) = -2.220000000000e+03;
    ref_lhs(5, 11) = -5.856350000000e+05;
    ref_lhs(6, 0) = -5.720000000000e+02;
    ref_lhs(6, 6) = 5.720000000000e+02;
    ref_lhs(7, 1) = -2.960000000000e+03;
    ref_lhs(7, 5) = -2.220000000000e+03;
    ref_lhs(7, 7) = 2.960000000000e+03;
    ref_lhs(7, 11) = -2.220000000000e+03;
    ref_lhs(8, 2) = -7.000000000000e+02;
    ref_lhs(8, 4) = 5.250000000000e+02;
    ref_lhs(8, 8) = 7.000000000000e+02;
    ref_lhs(8, 10) = 5.250000000000e+02;
    ref_lhs(9, 3) = -1.908000000000e+06;
    ref_lhs(9, 9) = 1.908000000000e+06;
    ref_lhs(10, 2) = -5.250000000000e+02;
    ref_lhs(10, 4) = -3.219062500000e+05;
    ref_lhs(10, 8) = 5.250000000000e+02;
    ref_lhs(10, 10) = 3.226937500000e+05;
    ref_lhs(11, 1) = 2.220000000000e+03;
    ref_lhs(11, 5) = -5.856350000000e+05;
    ref_lhs(11, 7) = -2.220000000000e+03;
    ref_lhs(11, 11) = 5.889650000000e+05;

    KRATOS_EXPECT_MATRIX_NEAR(ref_lhs, lhs, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(BushingElement_LHS_NonConstant, KratosStructuralMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("test");
    CreateBushingElementModelPart(r_model_part, false);

    Matrix lhs;
    r_model_part.GetElement(1).CalculateLeftHandSide(lhs, r_model_part.GetProcessInfo());

    Matrix ref_lhs(12, 12, 0.0);
    ref_lhs(0, 0) = 4.0000000000000000e+01;
    ref_lhs(0, 6) = -4.0000000000000000e+01;
    ref_lhs(1, 1) = 1.4965079365079362e+03;
    ref_lhs(1, 5) = 1.1223809523809521e+03;
    ref_lhs(1, 7) = -1.4965079365079362e+03;
    ref_lhs(1, 11) = 1.1223809523809521e+03;
    ref_lhs(2, 2) = 1.1400000000000000e+05;
    ref_lhs(2, 4) = -8.5500000000000000e+04;
    ref_lhs(2, 8) = -1.1400000000000000e+05;
    ref_lhs(2, 10) = -8.5500000000000000e+04;
    ref_lhs(3, 3) = 2.8000000000000000e+01;
    ref_lhs(3, 9) = -2.8000000000000000e+01;
    ref_lhs(4, 2) = -8.5500000000000000e+04;
    ref_lhs(4, 4) = 6.4181000000000000e+04;
    ref_lhs(4, 8) = 8.5500000000000000e+04;
    ref_lhs(4, 10) = 6.4069000000000000e+04;
    ref_lhs(5, 1) = 1.1223809523809521e+03;
    ref_lhs(5, 5) = 8.8378571428571399e+02;
    ref_lhs(5, 7) = -1.1223809523809521e+03;
    ref_lhs(5, 11) = 7.9978571428571399e+02;
    ref_lhs(6, 0) = -4.0000000000000000e+01;
    ref_lhs(6, 6) = 4.0000000000000000e+01;
    ref_lhs(7, 1) = -1.4965079365079362e+03;
    ref_lhs(7, 5) = -1.1223809523809521e+03;
    ref_lhs(7, 7) = 1.4965079365079362e+03;
    ref_lhs(7, 11) = -1.1223809523809521e+03;
    ref_lhs(8, 2) = -1.1400000000000000e+05;
    ref_lhs(8, 4) = 8.5500000000000000e+04;
    ref_lhs(8, 8) = 1.1400000000000000e+05;
    ref_lhs(8, 10) = 8.5500000000000000e+04;
    ref_lhs(9, 3) = -2.8000000000000000e+01;
    ref_lhs(9, 9) = 2.8000000000000000e+01;
    ref_lhs(10, 2) = -8.5500000000000000e+04;
    ref_lhs(10, 4) = 6.4069000000000000e+04;
    ref_lhs(10, 8) = 8.5500000000000000e+04;
    ref_lhs(10, 10) = 6.4181000000000000e+04;
    ref_lhs(11, 1) = 1.1223809523809521e+03;
    ref_lhs(11, 5) = 7.9978571428571399e+02;
    ref_lhs(11, 7) = -1.1223809523809521e+03;
    ref_lhs(11, 11) = 8.8378571428571399e+02;

    KRATOS_EXPECT_MATRIX_NEAR(ref_lhs, lhs, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(BushingElement_RHS_Constant, KratosStructuralMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("test");
    CreateBushingElementModelPart(r_model_part, true);

    Vector rhs, u;
    Matrix lhs;
    r_model_part.GetElement(1).GetValuesVector(u);
    r_model_part.GetElement(1).CalculateRightHandSide(rhs, r_model_part.GetProcessInfo());
    r_model_part.GetElement(1).CalculateLeftHandSide(lhs, r_model_part.GetProcessInfo());

    Vector ref_rhs = -prod(lhs, u);

    KRATOS_EXPECT_VECTOR_NEAR(ref_rhs, rhs, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(BushingElement_RHS_NonConstant, KratosStructuralMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("test");
    CreateBushingElementModelPart(r_model_part, false);

    Vector rhs, u;
    Matrix lhs;
    r_model_part.GetElement(1).GetValuesVector(u);
    r_model_part.GetElement(1).CalculateRightHandSide(rhs, r_model_part.GetProcessInfo());
    r_model_part.GetElement(1).CalculateLeftHandSide(lhs, r_model_part.GetProcessInfo());

    Vector ref_rhs = -prod(lhs, u);

    KRATOS_EXPECT_VECTOR_NEAR(ref_rhs, rhs, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(BushingElement_CalculateLocalSystem_Constant, KratosStructuralMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("test");
    CreateBushingElementModelPart(r_model_part, true);

    Vector rhs, ref_rhs;
    Matrix lhs, ref_lhs;
    r_model_part.GetElement(1).CalculateRightHandSide(ref_rhs, r_model_part.GetProcessInfo());
    r_model_part.GetElement(1).CalculateLeftHandSide(ref_lhs, r_model_part.GetProcessInfo());
    r_model_part.GetElement(1).CalculateLocalSystem(lhs, rhs, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_VECTOR_NEAR(ref_rhs, rhs, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(ref_lhs, lhs, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(BushingElement_CalculateLocalSystem_NonConstant, KratosStructuralMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("test");
    CreateBushingElementModelPart(r_model_part, false);

    Vector rhs, ref_rhs;
    Matrix lhs, ref_lhs;
    r_model_part.GetElement(1).CalculateRightHandSide(ref_rhs, r_model_part.GetProcessInfo());
    r_model_part.GetElement(1).CalculateLeftHandSide(ref_lhs, r_model_part.GetProcessInfo());
    r_model_part.GetElement(1).CalculateLocalSystem(lhs, rhs, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_VECTOR_NEAR(ref_rhs, rhs, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(ref_lhs, lhs, 1e-12);
}

} // namespace Testing
} // namespace Kratos
