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


namespace Kratos {
namespace Testing {
namespace {

void CreateMassElementTestModelPart(std::string const& rElementName, ModelPart& rModelPart)
{
    KRATOS_TRY;
    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    const Element& r_elem = KratosComponents<Element>::Get(rElementName);
    r_process_info[DOMAIN_SIZE] = 3;
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);

    std::vector<ModelPart::IndexType> node_ids(r_elem.GetGeometry().PointsNumber());
    for (std::size_t i = 0; i < r_elem.GetGeometry().PointsNumber(); ++i) {
        node_ids.at(i) = i + 1;
        const array_1d<double,3> vol_acc {1,2,3};
        noalias(rModelPart.Nodes()[i+1].FastGetSolutionStepValue(VOLUME_ACCELERATION)) = vol_acc;
    }

    auto p_prop = rModelPart.CreateNewProperties(1);
    rModelPart.CreateNewElement(rElementName, 1, node_ids, p_prop);

    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
    }

    (*p_prop)[DENSITY] = 1000.0;
    (*p_prop)[THICKNESS] = 0.1;
    (*p_prop)[CROSS_AREA] = 0.3;

    rModelPart.GetElement(1).Initialize(r_process_info);
    rModelPart.GetElement(1).Check(r_process_info);

    KRATOS_CATCH("CreateMassElementTestModelPart");
}

void ConductMassElementRHSTest(std::string const& rElementName, const std::vector<double>& rRefRHS)
{
    KRATOS_TRY;

    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("test");
    CreateMassElementTestModelPart(rElementName, test_model_part);
    const auto& r_process_info = test_model_part.GetProcessInfo();
    auto p_elem = test_model_part.pGetElement(1);

    Vector rhs;
    p_elem->CalculateRightHandSide(rhs, r_process_info);

    KRATOS_EXPECT_VECTOR_NEAR(rhs, rRefRHS, 1e-10);

    KRATOS_CATCH("ConductMassElementRHSTest");
}

void ConductMassElementMassMatrixTest(std::string const& rElementName, const Matrix& rRefMatrix)
{
    KRATOS_TRY;

    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("test");
    CreateMassElementTestModelPart(rElementName, test_model_part);
    auto p_elem = test_model_part.pGetElement(1);
    const auto& r_process_info = test_model_part.GetProcessInfo();

    Matrix lhs;
    p_elem->CalculateMassMatrix(lhs, r_process_info);

    KRATOS_EXPECT_MATRIX_NEAR(lhs, rRefMatrix, 1e-8);

    KRATOS_CATCH("ConductMassElementMassMatrixTest");
}

} // anonymous namespace


KRATOS_TEST_CASE_IN_SUITE(MassElement_2N_RHS, KratosStructuralMechanicsFastSuite)
{
    std::vector<double> ref_rhs {150,300,450,150,300,450};

    ConductMassElementRHSTest("LineMassElement3D2N", ref_rhs);
}

KRATOS_TEST_CASE_IN_SUITE(MassElement_2N_Massmatrix, KratosStructuralMechanicsFastSuite)
{
    Matrix ref_mass_matrix(6,6);
    noalias(ref_mass_matrix) = IdentityMatrix(6,6);
    ref_mass_matrix *= 150;

    ConductMassElementMassMatrixTest("LineMassElement3D2N", ref_mass_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(MassElement_3N_RHS, KratosStructuralMechanicsFastSuite)
{
    std::vector<double> ref_rhs {16.666666666667,33.3333333333,50,16.666666666667,33.3333333333,50,16.666666666667,33.3333333333,50};

    ConductMassElementRHSTest("SurfaceMassElement3D3N", ref_rhs);
}

KRATOS_TEST_CASE_IN_SUITE(MassElement_3N_Massmatrix, KratosStructuralMechanicsFastSuite)
{
    Matrix ref_mass_matrix(9,9);
    noalias(ref_mass_matrix) = IdentityMatrix(9,9);
    ref_mass_matrix *= 16.666666666667;

    ConductMassElementMassMatrixTest("SurfaceMassElement3D3N", ref_mass_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(MassElement_4N_RHS, KratosStructuralMechanicsFastSuite)
{
    std::vector<double> ref_rhs {25,50,75,25,50,75,25,50,75,25,50,75};

    ConductMassElementRHSTest("SurfaceMassElement3D4N", ref_rhs);
}

KRATOS_TEST_CASE_IN_SUITE(MassElement_4N_Massmatrix, KratosStructuralMechanicsFastSuite)
{
    Matrix ref_mass_matrix(12,12);
    noalias(ref_mass_matrix) = IdentityMatrix(12,12);
    ref_mass_matrix *= 25;

    ConductMassElementMassMatrixTest("SurfaceMassElement3D4N", ref_mass_matrix);
}

}
}
