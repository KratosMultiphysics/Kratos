// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//


// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "structural_mechanics_fast_suite.h"
#include "containers/model.h"
#include "includes/variables.h"
#include "includes/debug_helpers.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/linear_plane_strain.h"


namespace Kratos {
namespace Testing {
namespace {

void CreateShellTestModelPart(std::string const& rElementName, ModelPart& rModelPart)
{
    KRATOS_TRY;
    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    const Element& r_elem = KratosComponents<Element>::Get(rElementName);
    r_process_info[DOMAIN_SIZE] = 3;
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(ROTATION);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(THICKNESS);

    Matrix coordinates;
    r_elem.GetGeometry().PointsLocalCoordinates(coordinates);
    for (std::size_t i = 0; i < r_elem.GetGeometry().PointsNumber(); ++i) {
        rModelPart.CreateNewNode(i + 1, coordinates(i, 0), coordinates(i, 1), 0.0);
    }

    std::vector<ModelPart::IndexType> node_ids(r_elem.GetGeometry().PointsNumber());
    for (std::size_t i = 0; i < r_elem.GetGeometry().PointsNumber(); ++i) {
        node_ids.at(i) = i + 1;
    }

    auto p_prop = rModelPart.CreateNewProperties(1);
    rModelPart.CreateNewElement(rElementName, 1, node_ids, p_prop);
    rModelPart.SetBufferSize(2);

    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
        r_node.AddDof(ROTATION_X);
        r_node.AddDof(ROTATION_Y);
        r_node.AddDof(ROTATION_Z);
    }

    (*p_prop)[CONSTITUTIVE_LAW] = LinearPlaneStrain::Pointer(new LinearPlaneStrain());
    (*p_prop)[DENSITY] = 1000.0;
    (*p_prop)[YOUNG_MODULUS] = 1400000.0;
    (*p_prop)[THICKNESS] = 0.1;
    (*p_prop)[POISSON_RATIO] = 0.3;
    (*p_prop)[RAYLEIGH_ALPHA] = 0.02;
    (*p_prop)[RAYLEIGH_BETA] = 0.03;

    rModelPart.GetElement(1).Initialize(r_process_info);
    rModelPart.GetElement(1).Check(r_process_info);
    rModelPart.GetElement(1).InitializeSolutionStep(r_process_info);
    rModelPart.GetElement(1).InitializeNonLinearIteration(r_process_info);

    KRATOS_CATCH("CreateShellTestModelPart");
}

void ConductShellDampingMatrixTest(std::string const& rElementName, const Matrix& rRefMatrix)
{
    KRATOS_TRY;

    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("test");
    CreateShellTestModelPart(rElementName, test_model_part);
    auto p_elem = test_model_part.pGetElement(1);
    auto& r_process_info = test_model_part.GetProcessInfo();
    r_process_info[COMPUTE_LUMPED_MASS_MATRIX] = true;

    Matrix lhs;
    p_elem->CalculateDampingMatrix(lhs, r_process_info);

    // for writing the reference values
    // int counter = 0;
    // std::cout << std::endl<< std::endl << "    ";
    // for (int i=0;i<rRefMatrix.size1();++i) {
    //     for (int j=0;j<rRefMatrix.size2();++j) {
    //         const double val = lhs(i,j);
    //         if (std::abs(val) > 1e-10) {
    //             std::cout << std::setprecision(13) << "ref_damping_matrix("<<i<<","<<j<< ")=" << val << "; ";
    //             if (++counter%3 == 0) { std::cout << std::endl << "    "; }
    //         }
    //     }
    // }

    KRATOS_EXPECT_MATRIX_NEAR(lhs, rRefMatrix, 1e-8);

    KRATOS_CATCH("ConductShellDampingMatrixTest");
}

void ConductShellMassMatrixTest(std::string const& rElementName, const Matrix& rRefMatrix, const bool ComputeLumpedMassMatrix)
{
    KRATOS_TRY;

    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("test");
    CreateShellTestModelPart(rElementName, test_model_part);
    auto p_elem = test_model_part.pGetElement(1);
    auto& r_process_info = test_model_part.GetProcessInfo();
    r_process_info[COMPUTE_LUMPED_MASS_MATRIX] = ComputeLumpedMassMatrix;

    Matrix lhs;
    p_elem->CalculateMassMatrix(lhs, r_process_info);

    // for writing the reference values
    // int counter = 0;
    // std::cout << std::endl<< std::endl << "    ";
    // for (int i=0;i<rRefMatrix.size1();++i) {
    //     for (int j=0;j<rRefMatrix.size2();++j) {
    //         const double val = lhs(i,j);
    //         if (std::abs(val) > 1e-10) {
    //             std::cout << std::setprecision(13) << "ref_mass_matrix("<<i<<","<<j<< ")=" << val << "; ";
    //             if (++counter%3 == 0) { std::cout << std::endl << "    "; }
    //         }
    //     }
    // }

    KRATOS_EXPECT_MATRIX_NEAR(lhs, rRefMatrix, 1e-8);

    KRATOS_CATCH("ConductShellMassMatrixTest");
}

}

KRATOS_TEST_CASE_IN_SUITE(ShellElementCorotational_3N_LumpedMassMatrix, KratosStructuralMechanicsFastSuite)
{
    Matrix ref_mass_matrix(18, 18);
    ref_mass_matrix = ZeroMatrix(18,18);
    const double ref_val = 16.6666666666666666666666;
    for (std::size_t i=0; i<3; ++i) {
        std::size_t idx = i*6;
        ref_mass_matrix(idx, idx) = ref_val;
        ref_mass_matrix(idx+1, idx+1) = ref_val;
        ref_mass_matrix(idx+2, idx+2) = ref_val;
    }

    ConductShellMassMatrixTest("ShellThickElementCorotational3D3N", ref_mass_matrix, true);
    ConductShellMassMatrixTest("ShellThinElementCorotational3D3N", ref_mass_matrix, true);
}

KRATOS_TEST_CASE_IN_SUITE(ShellElementCorotational_4N_LumpedMassMatrix, KratosStructuralMechanicsFastSuite)
{
    Matrix ref_mass_matrix(24, 24);
    ref_mass_matrix = ZeroMatrix(24,24);
    const double ref_val = 100;
    for (std::size_t i=0; i<4; ++i) {
        std::size_t idx = i*6;
        ref_mass_matrix(idx, idx) = ref_val;
        ref_mass_matrix(idx+1, idx+1) = ref_val;
        ref_mass_matrix(idx+2, idx+2) = ref_val;
    }

    ConductShellMassMatrixTest("ShellThickElementCorotational3D4N", ref_mass_matrix, true);
    ConductShellMassMatrixTest("ShellThinElementCorotational3D4N", ref_mass_matrix, true);
}

KRATOS_TEST_CASE_IN_SUITE(ShellElementCorotational_3N_ConsistentMassMatrix, KratosStructuralMechanicsFastSuite)
{
    Matrix ref_mass_matrix(18, 18);
    ref_mass_matrix = ZeroMatrix(18,18);

    ref_mass_matrix(0,0)=8.333333333333; ref_mass_matrix(0,6)=4.166666666667; ref_mass_matrix(0,12)=4.166666666667;
    ref_mass_matrix(1,1)=8.333333333333; ref_mass_matrix(1,7)=4.166666666667; ref_mass_matrix(1,13)=4.166666666667;
    ref_mass_matrix(2,2)=8.333333333333; ref_mass_matrix(2,8)=4.166666666667; ref_mass_matrix(2,14)=4.166666666667;
    ref_mass_matrix(3,3)=0.006944444444444; ref_mass_matrix(3,9)=0.003472222222222; ref_mass_matrix(3,15)=0.003472222222222;
    ref_mass_matrix(4,4)=0.006944444444444; ref_mass_matrix(4,10)=0.003472222222222; ref_mass_matrix(4,16)=0.003472222222222;
    ref_mass_matrix(5,5)=0.006944444444444; ref_mass_matrix(5,11)=0.003472222222222; ref_mass_matrix(5,17)=0.003472222222222;
    ref_mass_matrix(6,0)=4.166666666667; ref_mass_matrix(6,6)=8.333333333333; ref_mass_matrix(6,12)=4.166666666667;
    ref_mass_matrix(7,1)=4.166666666667; ref_mass_matrix(7,7)=8.333333333333; ref_mass_matrix(7,13)=4.166666666667;
    ref_mass_matrix(8,2)=4.166666666667; ref_mass_matrix(8,8)=8.333333333333; ref_mass_matrix(8,14)=4.166666666667;
    ref_mass_matrix(9,3)=0.003472222222222; ref_mass_matrix(9,9)=0.006944444444444; ref_mass_matrix(9,15)=0.003472222222222;
    ref_mass_matrix(10,4)=0.003472222222222; ref_mass_matrix(10,10)=0.006944444444444; ref_mass_matrix(10,16)=0.003472222222222;
    ref_mass_matrix(11,5)=0.003472222222222; ref_mass_matrix(11,11)=0.006944444444444; ref_mass_matrix(11,17)=0.003472222222222;
    ref_mass_matrix(12,0)=4.166666666667; ref_mass_matrix(12,6)=4.166666666667; ref_mass_matrix(12,12)=8.333333333333;
    ref_mass_matrix(13,1)=4.166666666667; ref_mass_matrix(13,7)=4.166666666667; ref_mass_matrix(13,13)=8.333333333333;
    ref_mass_matrix(14,2)=4.166666666667; ref_mass_matrix(14,8)=4.166666666667; ref_mass_matrix(14,14)=8.333333333333;
    ref_mass_matrix(15,3)=0.003472222222222; ref_mass_matrix(15,9)=0.003472222222222; ref_mass_matrix(15,15)=0.006944444444444;
    ref_mass_matrix(16,4)=0.003472222222222; ref_mass_matrix(16,10)=0.003472222222222; ref_mass_matrix(16,16)=0.006944444444444;
    ref_mass_matrix(17,5)=0.003472222222222; ref_mass_matrix(17,11)=0.003472222222222; ref_mass_matrix(17,17)=0.006944444444444;

    ConductShellMassMatrixTest("ShellThickElementCorotational3D3N", ref_mass_matrix, false);
    ConductShellMassMatrixTest("ShellThinElementCorotational3D3N", ref_mass_matrix, false);
}

KRATOS_TEST_CASE_IN_SUITE(ShellElementCorotational_4N_ConsistentMassMatrix, KratosStructuralMechanicsFastSuite)
{
    Matrix ref_mass_matrix(24, 24);
    ref_mass_matrix = ZeroMatrix(24,24);

    ref_mass_matrix(0,0)=44.44444444444; ref_mass_matrix(0,6)=22.22222222222; ref_mass_matrix(0,12)=11.11111111111;
    ref_mass_matrix(0,18)=22.22222222222; ref_mass_matrix(1,1)=44.44444444444; ref_mass_matrix(1,7)=22.22222222222;
    ref_mass_matrix(1,13)=11.11111111111; ref_mass_matrix(1,19)=22.22222222222; ref_mass_matrix(2,2)=44.44444444444;
    ref_mass_matrix(2,8)=22.22222222222; ref_mass_matrix(2,14)=11.11111111111; ref_mass_matrix(2,20)=22.22222222222;
    ref_mass_matrix(3,3)=0.03703703703704; ref_mass_matrix(3,9)=0.01851851851852; ref_mass_matrix(3,15)=0.009259259259259;
    ref_mass_matrix(3,21)=0.01851851851852; ref_mass_matrix(4,4)=0.03703703703704; ref_mass_matrix(4,10)=0.01851851851852;
    ref_mass_matrix(4,16)=0.009259259259259; ref_mass_matrix(4,22)=0.01851851851852; ref_mass_matrix(5,5)=0.03703703703704;
    ref_mass_matrix(5,11)=0.01851851851852; ref_mass_matrix(5,17)=0.009259259259259; ref_mass_matrix(5,23)=0.01851851851852;
    ref_mass_matrix(6,0)=22.22222222222; ref_mass_matrix(6,6)=44.44444444444; ref_mass_matrix(6,12)=22.22222222222;
    ref_mass_matrix(6,18)=11.11111111111; ref_mass_matrix(7,1)=22.22222222222; ref_mass_matrix(7,7)=44.44444444444;
    ref_mass_matrix(7,13)=22.22222222222; ref_mass_matrix(7,19)=11.11111111111; ref_mass_matrix(8,2)=22.22222222222;
    ref_mass_matrix(8,8)=44.44444444444; ref_mass_matrix(8,14)=22.22222222222; ref_mass_matrix(8,20)=11.11111111111;
    ref_mass_matrix(9,3)=0.01851851851852; ref_mass_matrix(9,9)=0.03703703703704; ref_mass_matrix(9,15)=0.01851851851852;
    ref_mass_matrix(9,21)=0.009259259259259; ref_mass_matrix(10,4)=0.01851851851852; ref_mass_matrix(10,10)=0.03703703703704;
    ref_mass_matrix(10,16)=0.01851851851852; ref_mass_matrix(10,22)=0.009259259259259; ref_mass_matrix(11,5)=0.01851851851852;
    ref_mass_matrix(11,11)=0.03703703703704; ref_mass_matrix(11,17)=0.01851851851852; ref_mass_matrix(11,23)=0.009259259259259;
    ref_mass_matrix(12,0)=11.11111111111; ref_mass_matrix(12,6)=22.22222222222; ref_mass_matrix(12,12)=44.44444444444;
    ref_mass_matrix(12,18)=22.22222222222; ref_mass_matrix(13,1)=11.11111111111; ref_mass_matrix(13,7)=22.22222222222;
    ref_mass_matrix(13,13)=44.44444444444; ref_mass_matrix(13,19)=22.22222222222; ref_mass_matrix(14,2)=11.11111111111;
    ref_mass_matrix(14,8)=22.22222222222; ref_mass_matrix(14,14)=44.44444444444; ref_mass_matrix(14,20)=22.22222222222;
    ref_mass_matrix(15,3)=0.009259259259259; ref_mass_matrix(15,9)=0.01851851851852; ref_mass_matrix(15,15)=0.03703703703704;
    ref_mass_matrix(15,21)=0.01851851851852; ref_mass_matrix(16,4)=0.009259259259259; ref_mass_matrix(16,10)=0.01851851851852;
    ref_mass_matrix(16,16)=0.03703703703704; ref_mass_matrix(16,22)=0.01851851851852; ref_mass_matrix(17,5)=0.009259259259259;
    ref_mass_matrix(17,11)=0.01851851851852; ref_mass_matrix(17,17)=0.03703703703704; ref_mass_matrix(17,23)=0.01851851851852;
    ref_mass_matrix(18,0)=22.22222222222; ref_mass_matrix(18,6)=11.11111111111; ref_mass_matrix(18,12)=22.22222222222;
    ref_mass_matrix(18,18)=44.44444444444; ref_mass_matrix(19,1)=22.22222222222; ref_mass_matrix(19,7)=11.11111111111;
    ref_mass_matrix(19,13)=22.22222222222; ref_mass_matrix(19,19)=44.44444444444; ref_mass_matrix(20,2)=22.22222222222;
    ref_mass_matrix(20,8)=11.11111111111; ref_mass_matrix(20,14)=22.22222222222; ref_mass_matrix(20,20)=44.44444444444;
    ref_mass_matrix(21,3)=0.01851851851852; ref_mass_matrix(21,9)=0.009259259259259; ref_mass_matrix(21,15)=0.01851851851852;
    ref_mass_matrix(21,21)=0.03703703703704; ref_mass_matrix(22,4)=0.01851851851852; ref_mass_matrix(22,10)=0.009259259259259;
    ref_mass_matrix(22,16)=0.01851851851852; ref_mass_matrix(22,22)=0.03703703703704; ref_mass_matrix(23,5)=0.01851851851852;
    ref_mass_matrix(23,11)=0.009259259259259; ref_mass_matrix(23,17)=0.01851851851852; ref_mass_matrix(23,23)=0.03703703703704;

    ConductShellMassMatrixTest("ShellThickElementCorotational3D4N", ref_mass_matrix, false);
    ConductShellMassMatrixTest("ShellThinElementCorotational3D4N", ref_mass_matrix, false);
}

KRATOS_TEST_CASE_IN_SUITE(ShellThickElementCorotational3D3N_DampingMatrix, KratosStructuralMechanicsFastSuite)
{
    Matrix ref_damping_matrix(18, 18);
    ref_damping_matrix = ZeroMatrix(18,18);

    ref_damping_matrix(0,0)=3637.677406812; ref_damping_matrix(0,1)=2016.502080368; ref_damping_matrix(0,2)=-0.002727320768932;
    ref_damping_matrix(0,5)=-1.818216573204; ref_damping_matrix(0,6)=-2826.923076923; ref_damping_matrix(0,7)=-804.9649831734;
    ref_damping_matrix(0,11)=-1.818216573204; ref_damping_matrix(0,12)=-810.4196322112; ref_damping_matrix(0,13)=-1211.538461538;
    ref_damping_matrix(0,17)=-1.818216573204; ref_damping_matrix(1,0)=2016.502080368; ref_damping_matrix(1,1)=3637.677406812;
    ref_damping_matrix(1,2)=0.002727320768932; ref_damping_matrix(1,5)=1.818216573203; ref_damping_matrix(1,6)=-1211.538461538;
    ref_damping_matrix(1,7)=-810.4196322112; ref_damping_matrix(1,11)=1.818216573203; ref_damping_matrix(1,12)=-804.9649831734;
    ref_damping_matrix(1,13)=-2826.923076923; ref_damping_matrix(1,17)=1.818216573203; ref_damping_matrix(2,0)=-0.002727320768932;
    ref_damping_matrix(2,1)=0.002727320768932; ref_damping_matrix(2,2)=64.43590016185; ref_damping_matrix(2,3)=16.02564102564;
    ref_damping_matrix(2,4)=-16.02564102564; ref_damping_matrix(2,5)=0.001817304966422; ref_damping_matrix(2,7)=-0.002725957108889;
    ref_damping_matrix(2,8)=-32.05128205128; ref_damping_matrix(2,10)=-16.02564102564; ref_damping_matrix(2,11)=0.001817304966422;
    ref_damping_matrix(2,12)=0.002725957108889; ref_damping_matrix(2,14)=-32.05128205128; ref_damping_matrix(2,15)=16.02564102564;
    ref_damping_matrix(2,17)=0.001817304966422; ref_damping_matrix(3,2)=16.02564102564; ref_damping_matrix(3,3)=11.04166666667;
    ref_damping_matrix(3,4)=-1.682692307692; ref_damping_matrix(3,9)=-0.6730769230769; ref_damping_matrix(3,10)=1.009615384615;
    ref_damping_matrix(3,14)=-16.02564102564; ref_damping_matrix(3,15)=5.657051282051; ref_damping_matrix(3,16)=0.6730769230769;
    ref_damping_matrix(4,2)=-16.02564102564; ref_damping_matrix(4,3)=-1.682692307692; ref_damping_matrix(4,4)=11.04166666667;
    ref_damping_matrix(4,8)=16.02564102564; ref_damping_matrix(4,9)=0.6730769230769; ref_damping_matrix(4,10)=5.657051282051;
    ref_damping_matrix(4,15)=1.009615384615; ref_damping_matrix(4,16)=-0.6730769230769; ref_damping_matrix(5,0)=-1.818216573204;
    ref_damping_matrix(5,1)=1.818216573203; ref_damping_matrix(5,2)=0.001817304966422; ref_damping_matrix(5,5)=3.634615384615;
    ref_damping_matrix(5,7)=-1.817307465144; ref_damping_matrix(5,12)=1.817307465144; ref_damping_matrix(6,0)=-2826.923076923;
    ref_damping_matrix(6,1)=-1211.538461538; ref_damping_matrix(6,6)=2827.256410256; ref_damping_matrix(6,13)=1211.538461538;
    ref_damping_matrix(7,0)=-804.9649831734; ref_damping_matrix(7,1)=-810.4196322112; ref_damping_matrix(7,2)=-0.002725957108889;
    ref_damping_matrix(7,5)=-1.817307465144; ref_damping_matrix(7,7)=810.7516018826; ref_damping_matrix(7,11)=-1.817307465144;
    ref_damping_matrix(7,12)=804.9663468353; ref_damping_matrix(7,17)=-1.817307465144; ref_damping_matrix(8,2)=-32.05128205128;
    ref_damping_matrix(8,4)=16.02564102564; ref_damping_matrix(8,8)=32.38461538462; ref_damping_matrix(8,10)=16.02564102564;
    ref_damping_matrix(9,3)=-0.6730769230769; ref_damping_matrix(9,4)=0.6730769230769; ref_damping_matrix(9,9)=0.6730769230769;
    ref_damping_matrix(9,16)=-0.6730769230769; ref_damping_matrix(10,2)=-16.02564102564; ref_damping_matrix(10,3)=1.009615384615;
    ref_damping_matrix(10,4)=5.657051282051; ref_damping_matrix(10,8)=16.02564102564; ref_damping_matrix(10,10)=10.36858974359;
    ref_damping_matrix(10,15)=-1.009615384615; ref_damping_matrix(11,0)=-1.818216573204; ref_damping_matrix(11,1)=1.818216573203;
    ref_damping_matrix(11,2)=0.001817304966422; ref_damping_matrix(11,7)=-1.817307465144; ref_damping_matrix(11,11)=3.634615384615;
    ref_damping_matrix(11,12)=1.817307465144; ref_damping_matrix(12,0)=-810.4196322112; ref_damping_matrix(12,1)=-804.9649831734;
    ref_damping_matrix(12,2)=0.002725957108889; ref_damping_matrix(12,5)=1.817307465144; ref_damping_matrix(12,7)=804.9663468353;
    ref_damping_matrix(12,11)=1.817307465144; ref_damping_matrix(12,12)=810.7516018826; ref_damping_matrix(12,17)=1.817307465144;
    ref_damping_matrix(13,0)=-1211.538461538; ref_damping_matrix(13,1)=-2826.923076923; ref_damping_matrix(13,6)=1211.538461538;
    ref_damping_matrix(13,13)=2827.256410256; ref_damping_matrix(14,2)=-32.05128205128; ref_damping_matrix(14,3)=-16.02564102564;
    ref_damping_matrix(14,14)=32.38461538462; ref_damping_matrix(14,15)=-16.02564102564; ref_damping_matrix(15,2)=16.02564102564;
    ref_damping_matrix(15,3)=5.657051282051; ref_damping_matrix(15,4)=1.009615384615; ref_damping_matrix(15,10)=-1.009615384615;
    ref_damping_matrix(15,14)=-16.02564102564; ref_damping_matrix(15,15)=10.36858974359; ref_damping_matrix(16,3)=0.6730769230769;
    ref_damping_matrix(16,4)=-0.6730769230769; ref_damping_matrix(16,9)=-0.6730769230769; ref_damping_matrix(16,16)=0.6730769230769;
    ref_damping_matrix(17,0)=-1.818216573204; ref_damping_matrix(17,1)=1.818216573203; ref_damping_matrix(17,2)=0.001817304966422;
    ref_damping_matrix(17,7)=-1.817307465144; ref_damping_matrix(17,12)=1.817307465144; ref_damping_matrix(17,17)=3.634615384615;

    ConductShellDampingMatrixTest("ShellThickElementCorotational3D3N", ref_damping_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(ShellThinElementCorotational3D3N_DampingMatrix, KratosStructuralMechanicsFastSuite)
{
    Matrix ref_damping_matrix(18, 18);
    ref_damping_matrix = ZeroMatrix(18,18);

    ref_damping_matrix(0,0)=4203.324786325; ref_damping_matrix(0,1)=1450.854700855; ref_damping_matrix(0,5)=-957.264957265;
    ref_damping_matrix(0,6)=-2826.923076923; ref_damping_matrix(0,7)=-239.3162393162; ref_damping_matrix(0,11)=818.9102564103;
    ref_damping_matrix(0,12)=-1376.068376068; ref_damping_matrix(0,13)=-1211.538461538; ref_damping_matrix(0,17)=-998.3974358974;
    ref_damping_matrix(1,0)=1450.854700855; ref_damping_matrix(1,1)=4203.324786325; ref_damping_matrix(1,5)=957.264957265;
    ref_damping_matrix(1,6)=-1211.538461538; ref_damping_matrix(1,7)=-1376.068376068; ref_damping_matrix(1,11)=998.3974358974;
    ref_damping_matrix(1,12)=-239.3162393162; ref_damping_matrix(1,13)=-2826.923076923; ref_damping_matrix(1,17)=-818.9102564103;
    ref_damping_matrix(2,2)=46.77564102564; ref_damping_matrix(2,3)=8.918269230769; ref_damping_matrix(2,4)=-8.918269230769;
    ref_damping_matrix(2,8)=-23.22115384615; ref_damping_matrix(2,9)=2.692307692308; ref_damping_matrix(2,10)=-11.61057692308;
    ref_damping_matrix(2,14)=-23.22115384615; ref_damping_matrix(2,15)=11.61057692308; ref_damping_matrix(2,16)=-2.692307692308;
    ref_damping_matrix(3,2)=8.918269230769; ref_damping_matrix(3,3)=6.225961538462; ref_damping_matrix(3,4)=-1.598557692308;
    ref_damping_matrix(3,8)=0.4206730769231; ref_damping_matrix(3,9)=0.7992788461538; ref_damping_matrix(3,10)=1.219951923077;
    ref_damping_matrix(3,14)=-9.338942307692; ref_damping_matrix(3,15)=2.313701923077; ref_damping_matrix(3,16)=0.7992788461538;
    ref_damping_matrix(4,2)=-8.918269230769; ref_damping_matrix(4,3)=-1.598557692308; ref_damping_matrix(4,4)=6.225961538462;
    ref_damping_matrix(4,8)=9.338942307692; ref_damping_matrix(4,9)=0.7992788461538; ref_damping_matrix(4,10)=2.313701923077;
    ref_damping_matrix(4,14)=-0.4206730769231; ref_damping_matrix(4,15)=1.219951923077; ref_damping_matrix(4,16)=0.7992788461538;
    ref_damping_matrix(5,0)=-957.264957265; ref_damping_matrix(5,1)=957.264957265; ref_damping_matrix(5,5)=822.6495726496;
    ref_damping_matrix(5,6)=403.8461538462; ref_damping_matrix(5,7)=-553.4188034188; ref_damping_matrix(5,11)=142.094017094;
    ref_damping_matrix(5,12)=553.4188034188; ref_damping_matrix(5,13)=-403.8461538462; ref_damping_matrix(5,17)=142.094017094;
    ref_damping_matrix(6,0)=-2826.923076923; ref_damping_matrix(6,1)=-1211.538461538; ref_damping_matrix(6,5)=403.8461538462;
    ref_damping_matrix(6,6)=2827.256410256; ref_damping_matrix(6,11)=-706.7307692308; ref_damping_matrix(6,13)=1211.538461538;
    ref_damping_matrix(6,17)=302.8846153846; ref_damping_matrix(7,0)=-239.3162393162; ref_damping_matrix(7,1)=-1376.068376068;
    ref_damping_matrix(7,5)=-553.4188034188; ref_damping_matrix(7,7)=1376.401709402; ref_damping_matrix(7,11)=-695.5128205128;
    ref_damping_matrix(7,12)=239.3162393162; ref_damping_matrix(7,17)=112.1794871795; ref_damping_matrix(8,2)=-23.22115384615;
    ref_damping_matrix(8,3)=0.4206730769231; ref_damping_matrix(8,4)=9.338942307692; ref_damping_matrix(8,8)=26.07852564103;
    ref_damping_matrix(8,9)=2.1875; ref_damping_matrix(8,10)=11.52644230769; ref_damping_matrix(8,14)=-2.524038461538;
    ref_damping_matrix(8,15)=-0.08413461538462; ref_damping_matrix(8,16)=4.879807692308; ref_damping_matrix(9,2)=2.692307692308;
    ref_damping_matrix(9,3)=0.7992788461538; ref_damping_matrix(9,4)=0.7992788461538; ref_damping_matrix(9,8)=2.1875;
    ref_damping_matrix(9,9)=2.313701923077; ref_damping_matrix(9,10)=0.4206730769231; ref_damping_matrix(9,14)=-4.879807692308;
    ref_damping_matrix(9,15)=1.766826923077; ref_damping_matrix(9,16)=0.9675480769231; ref_damping_matrix(10,2)=-11.61057692308;
    ref_damping_matrix(10,3)=1.219951923077; ref_damping_matrix(10,4)=2.313701923077; ref_damping_matrix(10,8)=11.52644230769;
    ref_damping_matrix(10,9)=0.4206730769231; ref_damping_matrix(10,10)=7.445913461538; ref_damping_matrix(10,14)=0.08413461538462;
    ref_damping_matrix(10,15)=-1.724759615385; ref_damping_matrix(10,16)=1.766826923077; ref_damping_matrix(11,0)=818.9102564103;
    ref_damping_matrix(11,1)=998.3974358974; ref_damping_matrix(11,5)=142.094017094; ref_damping_matrix(11,6)=-706.7307692308;
    ref_damping_matrix(11,7)=-695.5128205128; ref_damping_matrix(11,11)=580.5288461538; ref_damping_matrix(11,12)=-112.1794871795;
    ref_damping_matrix(11,13)=-302.8846153846; ref_damping_matrix(11,17)=-139.2895299145; ref_damping_matrix(12,0)=-1376.068376068;
    ref_damping_matrix(12,1)=-239.3162393162; ref_damping_matrix(12,5)=553.4188034188; ref_damping_matrix(12,7)=239.3162393162;
    ref_damping_matrix(12,11)=-112.1794871795; ref_damping_matrix(12,12)=1376.401709402; ref_damping_matrix(12,17)=695.5128205128;
    ref_damping_matrix(13,0)=-1211.538461538; ref_damping_matrix(13,1)=-2826.923076923; ref_damping_matrix(13,5)=-403.8461538462;
    ref_damping_matrix(13,6)=1211.538461538; ref_damping_matrix(13,11)=-302.8846153846; ref_damping_matrix(13,13)=2827.256410256;
    ref_damping_matrix(13,17)=706.7307692308; ref_damping_matrix(14,2)=-23.22115384615; ref_damping_matrix(14,3)=-9.338942307692;
    ref_damping_matrix(14,4)=-0.4206730769231; ref_damping_matrix(14,8)=-2.524038461538; ref_damping_matrix(14,9)=-4.879807692308;
    ref_damping_matrix(14,10)=0.08413461538462; ref_damping_matrix(14,14)=26.07852564103; ref_damping_matrix(14,15)=-11.52644230769;
    ref_damping_matrix(14,16)=-2.1875; ref_damping_matrix(15,2)=11.61057692308; ref_damping_matrix(15,3)=2.313701923077;
    ref_damping_matrix(15,4)=1.219951923077; ref_damping_matrix(15,8)=-0.08413461538462; ref_damping_matrix(15,9)=1.766826923077;
    ref_damping_matrix(15,10)=-1.724759615385; ref_damping_matrix(15,14)=-11.52644230769; ref_damping_matrix(15,15)=7.445913461538;
    ref_damping_matrix(15,16)=0.4206730769231; ref_damping_matrix(16,2)=-2.692307692308; ref_damping_matrix(16,3)=0.7992788461538;
    ref_damping_matrix(16,4)=0.7992788461538; ref_damping_matrix(16,8)=4.879807692308; ref_damping_matrix(16,9)=0.9675480769231;
    ref_damping_matrix(16,10)=1.766826923077; ref_damping_matrix(16,14)=-2.1875; ref_damping_matrix(16,15)=0.4206730769231;
    ref_damping_matrix(16,16)=2.313701923077; ref_damping_matrix(17,0)=-998.3974358974; ref_damping_matrix(17,1)=-818.9102564103;
    ref_damping_matrix(17,5)=142.094017094; ref_damping_matrix(17,6)=302.8846153846; ref_damping_matrix(17,7)=112.1794871795;
    ref_damping_matrix(17,11)=-139.2895299145; ref_damping_matrix(17,12)=695.5128205128; ref_damping_matrix(17,13)=706.7307692308;
    ref_damping_matrix(17,17)=580.5288461538;

    ConductShellDampingMatrixTest("ShellThinElementCorotational3D3N", ref_damping_matrix);
}


KRATOS_TEST_CASE_IN_SUITE(ShellThickElementCorotational3D4N_DampingMatrix, KratosStructuralMechanicsFastSuite)
{
    Matrix ref_damping_matrix(24, 24);
    ref_damping_matrix = ZeroMatrix(24,24);

    ref_damping_matrix(0,0)=2338.538461538; ref_damping_matrix(0,1)=908.6538461538; ref_damping_matrix(0,5)=-269.2307692308;
    ref_damping_matrix(0,6)=-1326.923076923; ref_damping_matrix(0,7)=302.8846153846; ref_damping_matrix(0,11)=-134.6153846154;
    ref_damping_matrix(0,12)=-1500; ref_damping_matrix(0,13)=-908.6538461538; ref_damping_matrix(0,17)=-134.6153846154;
    ref_damping_matrix(0,18)=490.3846153846; ref_damping_matrix(0,19)=-302.8846153846; ref_damping_matrix(0,23)=-269.2307692308;
    ref_damping_matrix(1,0)=908.6538461538; ref_damping_matrix(1,1)=2338.538461538; ref_damping_matrix(1,5)=269.2307692308;
    ref_damping_matrix(1,6)=-302.8846153846; ref_damping_matrix(1,7)=490.3846153846; ref_damping_matrix(1,11)=269.2307692308;
    ref_damping_matrix(1,12)=-908.6538461538; ref_damping_matrix(1,13)=-1500; ref_damping_matrix(1,17)=134.6153846154;
    ref_damping_matrix(1,18)=302.8846153846; ref_damping_matrix(1,19)=-1326.923076923; ref_damping_matrix(1,23)=134.6153846154;
    ref_damping_matrix(2,2)=23.88868042527; ref_damping_matrix(2,3)=10.94434021263; ref_damping_matrix(2,4)=-10.94434021263;
    ref_damping_matrix(2,8)=-5.472170106316; ref_damping_matrix(2,9)=5.472170106316; ref_damping_matrix(2,10)=-10.94434021263;
    ref_damping_matrix(2,14)=-10.94434021263; ref_damping_matrix(2,15)=5.472170106316; ref_damping_matrix(2,16)=-5.472170106316;
    ref_damping_matrix(2,20)=-5.472170106316; ref_damping_matrix(2,21)=10.94434021263; ref_damping_matrix(2,22)=-5.472170106316;
    ref_damping_matrix(3,2)=10.94434021263; ref_damping_matrix(3,3)=12.96357098186; ref_damping_matrix(3,4)=-0.8413461538462;
    ref_damping_matrix(3,8)=5.472170106316; ref_damping_matrix(3,9)=5.808708567855; ref_damping_matrix(3,10)=0.1682692307692;
    ref_damping_matrix(3,14)=-5.472170106316; ref_damping_matrix(3,15)=4.462554721701; ref_damping_matrix(3,16)=0.8413461538462;
    ref_damping_matrix(3,20)=-10.94434021263; ref_damping_matrix(3,21)=9.598186366479; ref_damping_matrix(3,22)=-0.1682692307692;
    ref_damping_matrix(4,2)=-10.94434021263; ref_damping_matrix(4,3)=-0.8413461538462; ref_damping_matrix(4,4)=12.96357098186;
    ref_damping_matrix(4,8)=10.94434021263; ref_damping_matrix(4,9)=-0.1682692307692; ref_damping_matrix(4,10)=9.598186366479;
    ref_damping_matrix(4,14)=5.472170106316; ref_damping_matrix(4,15)=0.8413461538462; ref_damping_matrix(4,16)=4.462554721701;
    ref_damping_matrix(4,20)=-5.472170106316; ref_damping_matrix(4,21)=0.1682692307692; ref_damping_matrix(4,22)=5.808708567855;
    ref_damping_matrix(5,0)=-269.2307692308; ref_damping_matrix(5,1)=269.2307692308; ref_damping_matrix(5,5)=717.9487179487;
    ref_damping_matrix(5,6)=-134.6153846154; ref_damping_matrix(5,7)=-269.2307692308; ref_damping_matrix(5,11)=358.9743589744;
    ref_damping_matrix(5,12)=134.6153846154; ref_damping_matrix(5,13)=-134.6153846154; ref_damping_matrix(5,17)=179.4871794872;
    ref_damping_matrix(5,18)=269.2307692308; ref_damping_matrix(5,19)=134.6153846154; ref_damping_matrix(5,23)=358.9743589744;
    ref_damping_matrix(6,0)=-1326.923076923; ref_damping_matrix(6,1)=-302.8846153846; ref_damping_matrix(6,5)=-134.6153846154;
    ref_damping_matrix(6,6)=2338.538461538; ref_damping_matrix(6,7)=-908.6538461538; ref_damping_matrix(6,11)=-269.2307692308;
    ref_damping_matrix(6,12)=490.3846153846; ref_damping_matrix(6,13)=302.8846153846; ref_damping_matrix(6,17)=-269.2307692308;
    ref_damping_matrix(6,18)=-1500; ref_damping_matrix(6,19)=908.6538461538; ref_damping_matrix(6,23)=-134.6153846154;
    ref_damping_matrix(7,0)=302.8846153846; ref_damping_matrix(7,1)=490.3846153846; ref_damping_matrix(7,5)=-269.2307692308;
    ref_damping_matrix(7,6)=-908.6538461538; ref_damping_matrix(7,7)=2338.538461538; ref_damping_matrix(7,11)=-269.2307692308;
    ref_damping_matrix(7,12)=-302.8846153846; ref_damping_matrix(7,13)=-1326.923076923; ref_damping_matrix(7,17)=-134.6153846154;
    ref_damping_matrix(7,18)=908.6538461538; ref_damping_matrix(7,19)=-1500; ref_damping_matrix(7,23)=-134.6153846154;
    ref_damping_matrix(8,2)=-5.472170106316; ref_damping_matrix(8,3)=5.472170106316; ref_damping_matrix(8,4)=10.94434021263;
    ref_damping_matrix(8,8)=23.88868042527; ref_damping_matrix(8,9)=10.94434021263; ref_damping_matrix(8,10)=10.94434021263;
    ref_damping_matrix(8,14)=-5.472170106316; ref_damping_matrix(8,15)=10.94434021263; ref_damping_matrix(8,16)=5.472170106316;
    ref_damping_matrix(8,20)=-10.94434021263; ref_damping_matrix(8,21)=5.472170106316; ref_damping_matrix(8,22)=5.472170106316;
    ref_damping_matrix(9,2)=5.472170106316; ref_damping_matrix(9,3)=5.808708567855; ref_damping_matrix(9,4)=-0.1682692307692;
    ref_damping_matrix(9,8)=10.94434021263; ref_damping_matrix(9,9)=12.96357098186; ref_damping_matrix(9,10)=0.8413461538462;
    ref_damping_matrix(9,14)=-10.94434021263; ref_damping_matrix(9,15)=9.598186366479; ref_damping_matrix(9,16)=0.1682692307692;
    ref_damping_matrix(9,20)=-5.472170106316; ref_damping_matrix(9,21)=4.462554721701; ref_damping_matrix(9,22)=-0.8413461538462;
    ref_damping_matrix(10,2)=-10.94434021263; ref_damping_matrix(10,3)=0.1682692307692; ref_damping_matrix(10,4)=9.598186366479;
    ref_damping_matrix(10,8)=10.94434021263; ref_damping_matrix(10,9)=0.8413461538462; ref_damping_matrix(10,10)=12.96357098186;
    ref_damping_matrix(10,14)=5.472170106316; ref_damping_matrix(10,15)=-0.1682692307692; ref_damping_matrix(10,16)=5.808708567855;
    ref_damping_matrix(10,20)=-5.472170106316; ref_damping_matrix(10,21)=-0.8413461538462; ref_damping_matrix(10,22)=4.462554721701;
    ref_damping_matrix(11,0)=-134.6153846154; ref_damping_matrix(11,1)=269.2307692308; ref_damping_matrix(11,5)=358.9743589744;
    ref_damping_matrix(11,6)=-269.2307692308; ref_damping_matrix(11,7)=-269.2307692308; ref_damping_matrix(11,11)=717.9487179487;
    ref_damping_matrix(11,12)=269.2307692308; ref_damping_matrix(11,13)=-134.6153846154; ref_damping_matrix(11,17)=358.9743589744;
    ref_damping_matrix(11,18)=134.6153846154; ref_damping_matrix(11,19)=134.6153846154; ref_damping_matrix(11,23)=179.4871794872;
    ref_damping_matrix(12,0)=-1500; ref_damping_matrix(12,1)=-908.6538461538; ref_damping_matrix(12,5)=134.6153846154;
    ref_damping_matrix(12,6)=490.3846153846; ref_damping_matrix(12,7)=-302.8846153846; ref_damping_matrix(12,11)=269.2307692308;
    ref_damping_matrix(12,12)=2338.538461538; ref_damping_matrix(12,13)=908.6538461538; ref_damping_matrix(12,17)=269.2307692308;
    ref_damping_matrix(12,18)=-1326.923076923; ref_damping_matrix(12,19)=302.8846153846; ref_damping_matrix(12,23)=134.6153846154;
    ref_damping_matrix(13,0)=-908.6538461538; ref_damping_matrix(13,1)=-1500; ref_damping_matrix(13,5)=-134.6153846154;
    ref_damping_matrix(13,6)=302.8846153846; ref_damping_matrix(13,7)=-1326.923076923; ref_damping_matrix(13,11)=-134.6153846154;
    ref_damping_matrix(13,12)=908.6538461538; ref_damping_matrix(13,13)=2338.538461538; ref_damping_matrix(13,17)=-269.2307692308;
    ref_damping_matrix(13,18)=-302.8846153846; ref_damping_matrix(13,19)=490.3846153846; ref_damping_matrix(13,23)=-269.2307692308;
    ref_damping_matrix(14,2)=-10.94434021263; ref_damping_matrix(14,3)=-5.472170106316; ref_damping_matrix(14,4)=5.472170106316;
    ref_damping_matrix(14,8)=-5.472170106316; ref_damping_matrix(14,9)=-10.94434021263; ref_damping_matrix(14,10)=5.472170106316;
    ref_damping_matrix(14,14)=23.88868042527; ref_damping_matrix(14,15)=-10.94434021263; ref_damping_matrix(14,16)=10.94434021263;
    ref_damping_matrix(14,20)=-5.472170106316; ref_damping_matrix(14,21)=-5.472170106316; ref_damping_matrix(14,22)=10.94434021263;
    ref_damping_matrix(15,2)=5.472170106316; ref_damping_matrix(15,3)=4.462554721701; ref_damping_matrix(15,4)=0.8413461538462;
    ref_damping_matrix(15,8)=10.94434021263; ref_damping_matrix(15,9)=9.598186366479; ref_damping_matrix(15,10)=-0.1682692307692;
    ref_damping_matrix(15,14)=-10.94434021263; ref_damping_matrix(15,15)=12.96357098186; ref_damping_matrix(15,16)=-0.8413461538462;
    ref_damping_matrix(15,20)=-5.472170106316; ref_damping_matrix(15,21)=5.808708567855; ref_damping_matrix(15,22)=0.1682692307692;
    ref_damping_matrix(16,2)=-5.472170106316; ref_damping_matrix(16,3)=0.8413461538462; ref_damping_matrix(16,4)=4.462554721701;
    ref_damping_matrix(16,8)=5.472170106316; ref_damping_matrix(16,9)=0.1682692307692; ref_damping_matrix(16,10)=5.808708567855;
    ref_damping_matrix(16,14)=10.94434021263; ref_damping_matrix(16,15)=-0.8413461538462; ref_damping_matrix(16,16)=12.96357098186;
    ref_damping_matrix(16,20)=-10.94434021263; ref_damping_matrix(16,21)=-0.1682692307692; ref_damping_matrix(16,22)=9.598186366479;
    ref_damping_matrix(17,0)=-134.6153846154; ref_damping_matrix(17,1)=134.6153846154; ref_damping_matrix(17,5)=179.4871794872;
    ref_damping_matrix(17,6)=-269.2307692308; ref_damping_matrix(17,7)=-134.6153846154; ref_damping_matrix(17,11)=358.9743589744;
    ref_damping_matrix(17,12)=269.2307692308; ref_damping_matrix(17,13)=-269.2307692308; ref_damping_matrix(17,17)=717.9487179487;
    ref_damping_matrix(17,18)=134.6153846154; ref_damping_matrix(17,19)=269.2307692308; ref_damping_matrix(17,23)=358.9743589744;
    ref_damping_matrix(18,0)=490.3846153846; ref_damping_matrix(18,1)=302.8846153846; ref_damping_matrix(18,5)=269.2307692308;
    ref_damping_matrix(18,6)=-1500; ref_damping_matrix(18,7)=908.6538461538; ref_damping_matrix(18,11)=134.6153846154;
    ref_damping_matrix(18,12)=-1326.923076923; ref_damping_matrix(18,13)=-302.8846153846; ref_damping_matrix(18,17)=134.6153846154;
    ref_damping_matrix(18,18)=2338.538461538; ref_damping_matrix(18,19)=-908.6538461538; ref_damping_matrix(18,23)=269.2307692308;
    ref_damping_matrix(19,0)=-302.8846153846; ref_damping_matrix(19,1)=-1326.923076923; ref_damping_matrix(19,5)=134.6153846154;
    ref_damping_matrix(19,6)=908.6538461538; ref_damping_matrix(19,7)=-1500; ref_damping_matrix(19,11)=134.6153846154;
    ref_damping_matrix(19,12)=302.8846153846; ref_damping_matrix(19,13)=490.3846153846; ref_damping_matrix(19,17)=269.2307692308;
    ref_damping_matrix(19,18)=-908.6538461538; ref_damping_matrix(19,19)=2338.538461538; ref_damping_matrix(19,23)=269.2307692308;
    ref_damping_matrix(20,2)=-5.472170106316; ref_damping_matrix(20,3)=-10.94434021263; ref_damping_matrix(20,4)=-5.472170106316;
    ref_damping_matrix(20,8)=-10.94434021263; ref_damping_matrix(20,9)=-5.472170106316; ref_damping_matrix(20,10)=-5.472170106316;
    ref_damping_matrix(20,14)=-5.472170106316; ref_damping_matrix(20,15)=-5.472170106316; ref_damping_matrix(20,16)=-10.94434021263;
    ref_damping_matrix(20,20)=23.88868042527; ref_damping_matrix(20,21)=-10.94434021263; ref_damping_matrix(20,22)=-10.94434021263;
    ref_damping_matrix(21,2)=10.94434021263; ref_damping_matrix(21,3)=9.598186366479; ref_damping_matrix(21,4)=0.1682692307692;
    ref_damping_matrix(21,8)=5.472170106316; ref_damping_matrix(21,9)=4.462554721701; ref_damping_matrix(21,10)=-0.8413461538462;
    ref_damping_matrix(21,14)=-5.472170106316; ref_damping_matrix(21,15)=5.808708567855; ref_damping_matrix(21,16)=-0.1682692307692;
    ref_damping_matrix(21,20)=-10.94434021263; ref_damping_matrix(21,21)=12.96357098186; ref_damping_matrix(21,22)=0.8413461538462;
    ref_damping_matrix(22,2)=-5.472170106316; ref_damping_matrix(22,3)=-0.1682692307692; ref_damping_matrix(22,4)=5.808708567855;
    ref_damping_matrix(22,8)=5.472170106316; ref_damping_matrix(22,9)=-0.8413461538462; ref_damping_matrix(22,10)=4.462554721701;
    ref_damping_matrix(22,14)=10.94434021263; ref_damping_matrix(22,15)=0.1682692307692; ref_damping_matrix(22,16)=9.598186366479;
    ref_damping_matrix(22,20)=-10.94434021263; ref_damping_matrix(22,21)=0.8413461538462; ref_damping_matrix(22,22)=12.96357098186;
    ref_damping_matrix(23,0)=-269.2307692308; ref_damping_matrix(23,1)=134.6153846154; ref_damping_matrix(23,5)=358.9743589744;
    ref_damping_matrix(23,6)=-134.6153846154; ref_damping_matrix(23,7)=-134.6153846154; ref_damping_matrix(23,11)=179.4871794872;
    ref_damping_matrix(23,12)=134.6153846154; ref_damping_matrix(23,13)=-269.2307692308; ref_damping_matrix(23,17)=358.9743589744;
    ref_damping_matrix(23,18)=269.2307692308; ref_damping_matrix(23,19)=269.2307692308; ref_damping_matrix(23,23)=717.9487179487;

    ConductShellDampingMatrixTest("ShellThickElementCorotational3D4N", ref_damping_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(ShellThinElementCorotational3D4N_DampingMatrix, KratosStructuralMechanicsFastSuite)
{
    Matrix ref_damping_matrix(24, 24);
    ref_damping_matrix = ZeroMatrix(24,24);

    ref_damping_matrix(0,0)=2138.346153846; ref_damping_matrix(0,1)=908.6538461538; ref_damping_matrix(0,5)=-759.2307692308;
    ref_damping_matrix(0,6)=-1126.730769231; ref_damping_matrix(0,7)=302.8846153846; ref_damping_matrix(0,11)=355.3846153846;
    ref_damping_matrix(0,12)=-1700.192307692; ref_damping_matrix(0,13)=-908.6538461538; ref_damping_matrix(0,17)=-452.3076923077;
    ref_damping_matrix(0,18)=690.5769230769; ref_damping_matrix(0,19)=-302.8846153846; ref_damping_matrix(0,23)=48.46153846154;
    ref_damping_matrix(1,0)=908.6538461538; ref_damping_matrix(1,1)=2138.346153846; ref_damping_matrix(1,5)=759.2307692308;
    ref_damping_matrix(1,6)=-302.8846153846; ref_damping_matrix(1,7)=690.5769230769; ref_damping_matrix(1,11)=-48.46153846154;
    ref_damping_matrix(1,12)=-908.6538461538; ref_damping_matrix(1,13)=-1700.192307692; ref_damping_matrix(1,17)=452.3076923077;
    ref_damping_matrix(1,18)=302.8846153846; ref_damping_matrix(1,19)=-1126.730769231; ref_damping_matrix(1,23)=-355.3846153846;
    ref_damping_matrix(2,2)=13.77884615385; ref_damping_matrix(2,3)=5.721153846154; ref_damping_matrix(2,4)=-5.721153846154;
    ref_damping_matrix(2,8)=-4.711538461538; ref_damping_matrix(2,9)=1.346153846154; ref_damping_matrix(2,10)=-4.711538461538;
    ref_damping_matrix(2,14)=-2.355769230769; ref_damping_matrix(2,15)=2.355769230769; ref_damping_matrix(2,16)=-2.355769230769;
    ref_damping_matrix(2,20)=-4.711538461538; ref_damping_matrix(2,21)=4.711538461538; ref_damping_matrix(2,22)=-1.346153846154;
    ref_damping_matrix(3,2)=5.721153846154; ref_damping_matrix(3,3)=6.394230769231; ref_damping_matrix(3,4)=-2.019230769231;
    ref_damping_matrix(3,8)=1.346153846154; ref_damping_matrix(3,9)=3.028846153846; ref_damping_matrix(3,14)=-2.355769230769;
    ref_damping_matrix(3,15)=1.682692307692; ref_damping_matrix(3,20)=-4.711538461538; ref_damping_matrix(3,21)=3.028846153846;
    ref_damping_matrix(4,2)=-5.721153846154; ref_damping_matrix(4,3)=-2.019230769231; ref_damping_matrix(4,4)=6.394230769231;
    ref_damping_matrix(4,8)=4.711538461538; ref_damping_matrix(4,10)=3.028846153846; ref_damping_matrix(4,14)=2.355769230769;
    ref_damping_matrix(4,16)=1.682692307692; ref_damping_matrix(4,20)=-1.346153846154; ref_damping_matrix(4,22)=3.028846153846;
    ref_damping_matrix(5,0)=-759.2307692308; ref_damping_matrix(5,1)=759.2307692308; ref_damping_matrix(5,5)=1233.076923077;
    ref_damping_matrix(5,6)=355.3846153846; ref_damping_matrix(5,7)=48.46153846154; ref_damping_matrix(5,12)=452.3076923077;
    ref_damping_matrix(5,13)=-452.3076923077; ref_damping_matrix(5,17)=382.3076923077; ref_damping_matrix(5,18)=-48.46153846154;
    ref_damping_matrix(5,19)=-355.3846153846; ref_damping_matrix(6,0)=-1126.730769231; ref_damping_matrix(6,1)=-302.8846153846;
    ref_damping_matrix(6,5)=355.3846153846; ref_damping_matrix(6,6)=2138.346153846; ref_damping_matrix(6,7)=-908.6538461538;
    ref_damping_matrix(6,11)=-759.2307692308; ref_damping_matrix(6,12)=690.5769230769; ref_damping_matrix(6,13)=302.8846153846;
    ref_damping_matrix(6,17)=48.46153846154; ref_damping_matrix(6,18)=-1700.192307692; ref_damping_matrix(6,19)=908.6538461538;
    ref_damping_matrix(6,23)=-452.3076923077; ref_damping_matrix(7,0)=302.8846153846; ref_damping_matrix(7,1)=690.5769230769;
    ref_damping_matrix(7,5)=48.46153846154; ref_damping_matrix(7,6)=-908.6538461538; ref_damping_matrix(7,7)=2138.346153846;
    ref_damping_matrix(7,11)=-759.2307692308; ref_damping_matrix(7,12)=-302.8846153846; ref_damping_matrix(7,13)=-1126.730769231;
    ref_damping_matrix(7,17)=355.3846153846; ref_damping_matrix(7,18)=908.6538461538; ref_damping_matrix(7,19)=-1700.192307692;
    ref_damping_matrix(7,23)=-452.3076923077; ref_damping_matrix(8,2)=-4.711538461538; ref_damping_matrix(8,3)=1.346153846154;
    ref_damping_matrix(8,4)=4.711538461538; ref_damping_matrix(8,8)=13.77884615385; ref_damping_matrix(8,9)=5.721153846154;
    ref_damping_matrix(8,10)=5.721153846154; ref_damping_matrix(8,14)=-4.711538461538; ref_damping_matrix(8,15)=4.711538461538;
    ref_damping_matrix(8,16)=1.346153846154; ref_damping_matrix(8,20)=-2.355769230769; ref_damping_matrix(8,21)=2.355769230769;
    ref_damping_matrix(8,22)=2.355769230769; ref_damping_matrix(9,2)=1.346153846154; ref_damping_matrix(9,3)=3.028846153846;
    ref_damping_matrix(9,8)=5.721153846154; ref_damping_matrix(9,9)=6.394230769231; ref_damping_matrix(9,10)=2.019230769231;
    ref_damping_matrix(9,14)=-4.711538461538; ref_damping_matrix(9,15)=3.028846153846; ref_damping_matrix(9,20)=-2.355769230769;
    ref_damping_matrix(9,21)=1.682692307692; ref_damping_matrix(10,2)=-4.711538461538; ref_damping_matrix(10,4)=3.028846153846;
    ref_damping_matrix(10,8)=5.721153846154; ref_damping_matrix(10,9)=2.019230769231; ref_damping_matrix(10,10)=6.394230769231;
    ref_damping_matrix(10,14)=1.346153846154; ref_damping_matrix(10,16)=3.028846153846; ref_damping_matrix(10,20)=-2.355769230769;
    ref_damping_matrix(10,22)=1.682692307692; ref_damping_matrix(11,0)=355.3846153846; ref_damping_matrix(11,1)=-48.46153846154;
    ref_damping_matrix(11,6)=-759.2307692308; ref_damping_matrix(11,7)=-759.2307692308; ref_damping_matrix(11,11)=1233.076923077;
    ref_damping_matrix(11,12)=-48.46153846154; ref_damping_matrix(11,13)=355.3846153846; ref_damping_matrix(11,18)=452.3076923077;
    ref_damping_matrix(11,19)=452.3076923077; ref_damping_matrix(11,23)=382.3076923077; ref_damping_matrix(12,0)=-1700.192307692;
    ref_damping_matrix(12,1)=-908.6538461538; ref_damping_matrix(12,5)=452.3076923077; ref_damping_matrix(12,6)=690.5769230769;
    ref_damping_matrix(12,7)=-302.8846153846; ref_damping_matrix(12,11)=-48.46153846154; ref_damping_matrix(12,12)=2138.346153846;
    ref_damping_matrix(12,13)=908.6538461538; ref_damping_matrix(12,17)=759.2307692308; ref_damping_matrix(12,18)=-1126.730769231;
    ref_damping_matrix(12,19)=302.8846153846; ref_damping_matrix(12,23)=-355.3846153846; ref_damping_matrix(13,0)=-908.6538461538;
    ref_damping_matrix(13,1)=-1700.192307692; ref_damping_matrix(13,5)=-452.3076923077; ref_damping_matrix(13,6)=302.8846153846;
    ref_damping_matrix(13,7)=-1126.730769231; ref_damping_matrix(13,11)=355.3846153846; ref_damping_matrix(13,12)=908.6538461538;
    ref_damping_matrix(13,13)=2138.346153846; ref_damping_matrix(13,17)=-759.2307692308; ref_damping_matrix(13,18)=-302.8846153846;
    ref_damping_matrix(13,19)=690.5769230769; ref_damping_matrix(13,23)=48.46153846154; ref_damping_matrix(14,2)=-2.355769230769;
    ref_damping_matrix(14,3)=-2.355769230769; ref_damping_matrix(14,4)=2.355769230769; ref_damping_matrix(14,8)=-4.711538461538;
    ref_damping_matrix(14,9)=-4.711538461538; ref_damping_matrix(14,10)=1.346153846154; ref_damping_matrix(14,14)=13.77884615385;
    ref_damping_matrix(14,15)=-5.721153846154; ref_damping_matrix(14,16)=5.721153846154; ref_damping_matrix(14,20)=-4.711538461538;
    ref_damping_matrix(14,21)=-1.346153846154; ref_damping_matrix(14,22)=4.711538461538; ref_damping_matrix(15,2)=2.355769230769;
    ref_damping_matrix(15,3)=1.682692307692; ref_damping_matrix(15,8)=4.711538461538; ref_damping_matrix(15,9)=3.028846153846;
    ref_damping_matrix(15,14)=-5.721153846154; ref_damping_matrix(15,15)=6.394230769231; ref_damping_matrix(15,16)=-2.019230769231;
    ref_damping_matrix(15,20)=-1.346153846154; ref_damping_matrix(15,21)=3.028846153846; ref_damping_matrix(16,2)=-2.355769230769;
    ref_damping_matrix(16,4)=1.682692307692; ref_damping_matrix(16,8)=1.346153846154; ref_damping_matrix(16,10)=3.028846153846;
    ref_damping_matrix(16,14)=5.721153846154; ref_damping_matrix(16,15)=-2.019230769231; ref_damping_matrix(16,16)=6.394230769231;
    ref_damping_matrix(16,20)=-4.711538461538; ref_damping_matrix(16,22)=3.028846153846; ref_damping_matrix(17,0)=-452.3076923077;
    ref_damping_matrix(17,1)=452.3076923077; ref_damping_matrix(17,5)=382.3076923077; ref_damping_matrix(17,6)=48.46153846154;
    ref_damping_matrix(17,7)=355.3846153846; ref_damping_matrix(17,12)=759.2307692308; ref_damping_matrix(17,13)=-759.2307692308;
    ref_damping_matrix(17,17)=1233.076923077; ref_damping_matrix(17,18)=-355.3846153846; ref_damping_matrix(17,19)=-48.46153846154;
    ref_damping_matrix(18,0)=690.5769230769; ref_damping_matrix(18,1)=302.8846153846; ref_damping_matrix(18,5)=-48.46153846154;
    ref_damping_matrix(18,6)=-1700.192307692; ref_damping_matrix(18,7)=908.6538461538; ref_damping_matrix(18,11)=452.3076923077;
    ref_damping_matrix(18,12)=-1126.730769231; ref_damping_matrix(18,13)=-302.8846153846; ref_damping_matrix(18,17)=-355.3846153846;
    ref_damping_matrix(18,18)=2138.346153846; ref_damping_matrix(18,19)=-908.6538461538; ref_damping_matrix(18,23)=759.2307692308;
    ref_damping_matrix(19,0)=-302.8846153846; ref_damping_matrix(19,1)=-1126.730769231; ref_damping_matrix(19,5)=-355.3846153846;
    ref_damping_matrix(19,6)=908.6538461538; ref_damping_matrix(19,7)=-1700.192307692; ref_damping_matrix(19,11)=452.3076923077;
    ref_damping_matrix(19,12)=302.8846153846; ref_damping_matrix(19,13)=690.5769230769; ref_damping_matrix(19,17)=-48.46153846154;
    ref_damping_matrix(19,18)=-908.6538461538; ref_damping_matrix(19,19)=2138.346153846; ref_damping_matrix(19,23)=759.2307692308;
    ref_damping_matrix(20,2)=-4.711538461538; ref_damping_matrix(20,3)=-4.711538461538; ref_damping_matrix(20,4)=-1.346153846154;
    ref_damping_matrix(20,8)=-2.355769230769; ref_damping_matrix(20,9)=-2.355769230769; ref_damping_matrix(20,10)=-2.355769230769;
    ref_damping_matrix(20,14)=-4.711538461538; ref_damping_matrix(20,15)=-1.346153846154; ref_damping_matrix(20,16)=-4.711538461538;
    ref_damping_matrix(20,20)=13.77884615385; ref_damping_matrix(20,21)=-5.721153846154; ref_damping_matrix(20,22)=-5.721153846154;
    ref_damping_matrix(21,2)=4.711538461538; ref_damping_matrix(21,3)=3.028846153846; ref_damping_matrix(21,8)=2.355769230769;
    ref_damping_matrix(21,9)=1.682692307692; ref_damping_matrix(21,14)=-1.346153846154; ref_damping_matrix(21,15)=3.028846153846;
    ref_damping_matrix(21,20)=-5.721153846154; ref_damping_matrix(21,21)=6.394230769231; ref_damping_matrix(21,22)=2.019230769231;
    ref_damping_matrix(22,2)=-1.346153846154; ref_damping_matrix(22,4)=3.028846153846; ref_damping_matrix(22,8)=2.355769230769;
    ref_damping_matrix(22,10)=1.682692307692; ref_damping_matrix(22,14)=4.711538461538; ref_damping_matrix(22,16)=3.028846153846;
    ref_damping_matrix(22,20)=-5.721153846154; ref_damping_matrix(22,21)=2.019230769231; ref_damping_matrix(22,22)=6.394230769231;
    ref_damping_matrix(23,0)=48.46153846154; ref_damping_matrix(23,1)=-355.3846153846; ref_damping_matrix(23,6)=-452.3076923077;
    ref_damping_matrix(23,7)=-452.3076923077; ref_damping_matrix(23,11)=382.3076923077; ref_damping_matrix(23,12)=-355.3846153846;
    ref_damping_matrix(23,13)=48.46153846154; ref_damping_matrix(23,18)=759.2307692308; ref_damping_matrix(23,19)=759.2307692308;
    ref_damping_matrix(23,23)=1233.076923077;

    ConductShellDampingMatrixTest("ShellThinElementCorotational3D4N", ref_damping_matrix);
}

}
}
