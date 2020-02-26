// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Yaman Zendaki , Hoang-Giang Bui
//

// System includes

// External includes

// Project includes
//#include "includes/gid_io.h"
#include "testing/testing.h"
#include "containers/model.h"
#include "spaces/ublas_space.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/variable_utils.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"

namespace Kratos {
namespace Testing {

typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
typedef LinearSolver<SparseSpaceType,LocalSpaceType> LinearSolverType;
typedef SkylineLUFactorizationSolver<SparseSpaceType,  LocalSpaceType > SkylineLUFactorizationSolverType;
typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
typedef ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverType;
typedef Scheme< SparseSpaceType, LocalSpaceType >  SchemeType;
typedef ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType> ResidualBasedIncrementalUpdateStaticSchemeType;
typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;
typedef ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedLinearStrategyType;

KRATOS_TEST_CASE_IN_SUITE(PatchTestMPCPlateTension, KratosStructuralMechanicsFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Main");
    model_part.SetBufferSize(2);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(REACTION);
    //Material properties
    auto p_Properties = model_part.CreateNewProperties(0);
    p_Properties->SetValue(YOUNG_MODULUS, 210e3);
    p_Properties->SetValue(POISSON_RATIO, 0.3);
    p_Properties->SetValue(DENSITY, 1.0);
    auto p_constitutive_law = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw").Clone();
    p_Properties->SetValue(CONSTITUTIVE_LAW,p_constitutive_law);
    //Create Nodes
    auto pnode1 = model_part.CreateNewNode(1, 0.00,0.00,0.00);
    auto pnode2 = model_part.CreateNewNode(2, 1.00,0.00,0.00);
    auto pnode3 = model_part.CreateNewNode(3, 1.00,1.00,0.00);
    auto pnode4 = model_part.CreateNewNode(4, 0.00,1.00,0.00);
    auto pnode5 = model_part.CreateNewNode(5, 0.50,0.00,0.00);
    auto pnode6 = model_part.CreateNewNode(6, 1.00,0.50,0.00);
    auto pnode7 = model_part.CreateNewNode(7, 0.50,1.00,0.00);
    auto pnode8 = model_part.CreateNewNode(8, 0.50,0.50,0.00);
    //Create Elements
    model_part.CreateNewElement("SmallDisplacementElement2D4N", 1, {{1, 5, 7, 4}}, p_Properties);
    model_part.CreateNewElement("SmallDisplacementElement2D4N", 2, {{5, 2, 6, 8}}, p_Properties);
    model_part.CreateNewElement("SmallDisplacementElement2D4N", 3, {{8, 6, 3, 7}}, p_Properties);
    //Add DOF
    VariableUtils().AddDofWithReaction(DISPLACEMENT_X, REACTION_X, model_part);
    VariableUtils().AddDofWithReaction(DISPLACEMENT_Y, REACTION_Y, model_part);
    //Apply BC
    pnode1->Fix(DISPLACEMENT_X);
    pnode1->Fix(DISPLACEMENT_Y);
    pnode2->Fix(DISPLACEMENT_Y);
    pnode4->Fix(DISPLACEMENT_X);
    //Creat Load
    auto p_condition1 = model_part.CreateNewCondition("LineLoadCondition2D2N", 1, std::vector<IndexType>({2,6}), p_Properties);
    auto p_condition2 = model_part.CreateNewCondition("LineLoadCondition2D2N", 2, std::vector<IndexType>({6,3}), p_Properties);
    array_1d<double, 3> load;
    load[0] = 1000.0; load[1] = 0.0; load[2] = 0.0;
    p_condition1->SetValue(LINE_LOAD, load);
    p_condition2->SetValue(LINE_LOAD, load);
    //Apply Constraints
    std::vector< Dof<double>::Pointer >slave_dofs(2);
    slave_dofs[0] = pnode8->pGetDof(DISPLACEMENT_X);
    slave_dofs[1] = pnode8->pGetDof(DISPLACEMENT_Y);
    std::vector< Dof<double>::Pointer >master_dofs(4);
    master_dofs[0] = pnode7->pGetDof(DISPLACEMENT_X);
    master_dofs[1] = pnode5->pGetDof(DISPLACEMENT_X);
    master_dofs[2] = pnode7->pGetDof(DISPLACEMENT_Y);
    master_dofs[3] = pnode5->pGetDof(DISPLACEMENT_Y);
    Matrix relation_matrix = ZeroMatrix(2,4);
    relation_matrix(0,0) = 0.5;
    relation_matrix(0,1) = 0.5;
    relation_matrix(1,2) = 0.5;
    relation_matrix(1,3) = 0.5;
    Vector constant_vector = ZeroVector(2);
    model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, master_dofs, slave_dofs, relation_matrix, constant_vector);
    //Setting up System
    SchemeType::Pointer p_scheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticSchemeType>();
    LinearSolverType::Pointer p_solver = Kratos::make_shared<SkylineLUFactorizationSolverType>();
    BuilderAndSolverType::Pointer p_builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolverType>(p_solver);
    SolvingStrategyType::Pointer p_strategy = Kratos::make_shared<ResidualBasedLinearStrategyType>(model_part, p_scheme, p_solver, p_builder_and_solver, true);
    p_strategy->SetEchoLevel(0);
    p_strategy->Solve();
    //Output to GID
    // GidIO<> gid_io("TEST_MPC_Plate_Tension", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
    // const double label = static_cast<double>(nl_iter);
    // gid_io.InitializeMesh(label);
    // gid_io.WriteMesh(model_part.GetMesh());
    // gid_io.FinalizeMesh();
    // gid_io.InitializeResults(label, model_part.GetMesh());
    // gid_io.WriteNodalResults(DISPLACEMENT, model_part.Nodes(), label, 0);
    // gid_io.WriteNodalResults(REACTION, model_part.Nodes(), label, 0);
    // gid_io.PrintOnGaussPoints(STRESSES, model_part , label ,0);
    //Check Resault
    // std::cout << "DISPLACEMENTS" << std::endl;
    // std::cout << "u1x=" << pnode1->GetSolutionStepValue(DISPLACEMENT_X) << std::endl;
    // std::cout << "u2x=" << pnode2->GetSolutionStepValue(DISPLACEMENT_X) << std::endl;
    // std::cout << "u3x=" << pnode3->GetSolutionStepValue(DISPLACEMENT_X) << std::endl;
    // std::cout << "u4x=" << pnode4->GetSolutionStepValue(DISPLACEMENT_X) << std::endl;
    // std::cout << "u5x=" << pnode5->GetSolutionStepValue(DISPLACEMENT_X) << std::endl;
    // std::cout << "u6x=" << pnode6->GetSolutionStepValue(DISPLACEMENT_X) << std::endl;
    // std::cout << "u7x=" << pnode7->GetSolutionStepValue(DISPLACEMENT_X) << std::endl;
    // std::cout << "u8x=" << pnode8->GetSolutionStepValue(DISPLACEMENT_X) << std::endl;
    // std::cout << "u1y=" << pnode1->GetSolutionStepValue(DISPLACEMENT_Y) << std::endl;
    // std::cout << "u2y=" << pnode2->GetSolutionStepValue(DISPLACEMENT_Y) << std::endl;
    // std::cout << "u3y=" << pnode3->GetSolutionStepValue(DISPLACEMENT_Y) << std::endl;
    // std::cout << "u4y=" << pnode4->GetSolutionStepValue(DISPLACEMENT_Y) << std::endl;
    // std::cout << "u5y=" << pnode5->GetSolutionStepValue(DISPLACEMENT_Y) << std::endl;
    // std::cout << "u6y=" << pnode6->GetSolutionStepValue(DISPLACEMENT_Y) << std::endl;
    // std::cout << "u7y=" << pnode7->GetSolutionStepValue(DISPLACEMENT_Y) << std::endl;
    // std::cout << "u8y=" << pnode8->GetSolutionStepValue(DISPLACEMENT_Y) << std::endl;

    KRATOS_CHECK_GREATER(std::abs(pnode8->GetSolutionStepValue(DISPLACEMENT_X)),0);
    KRATOS_CHECK_DOUBLE_EQUAL(pnode8->GetSolutionStepValue(DISPLACEMENT_X),pnode5->GetSolutionStepValue(DISPLACEMENT_X));
    KRATOS_CHECK_DOUBLE_EQUAL(pnode8->GetSolutionStepValue(DISPLACEMENT_X),pnode7->GetSolutionStepValue(DISPLACEMENT_X));
    KRATOS_CHECK_DOUBLE_EQUAL(pnode3->GetSolutionStepValue(DISPLACEMENT_X),pnode6->GetSolutionStepValue(DISPLACEMENT_X));
    KRATOS_CHECK_DOUBLE_EQUAL(pnode6->GetSolutionStepValue(DISPLACEMENT_X),pnode2->GetSolutionStepValue(DISPLACEMENT_X));

    KRATOS_CHECK_DOUBLE_EQUAL(pnode4->GetSolutionStepValue(DISPLACEMENT_Y),pnode7->GetSolutionStepValue(DISPLACEMENT_Y));
    KRATOS_CHECK_DOUBLE_EQUAL(pnode7->GetSolutionStepValue(DISPLACEMENT_Y),pnode3->GetSolutionStepValue(DISPLACEMENT_Y));
    KRATOS_CHECK_DOUBLE_EQUAL(pnode8->GetSolutionStepValue(DISPLACEMENT_Y),pnode6->GetSolutionStepValue(DISPLACEMENT_Y));
}

} // namespace Testing

} // namespace Kratos
