//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes
#include <limits>

/* External includes */

/* Project includes */
#include "testing/testing.h"
#include "containers/model.h"
#include "geometries/point_3d.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/auxiliar_files_for_cpp_unnitest/test_element.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/explicit_solving_strategy_runge_kutta_4.h"

namespace Kratos
{
namespace Testing
{
    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef ExplicitBuilder< SparseSpaceType, LocalSpaceType > ExplicitBuilderType;
    typedef ExplicitSolvingStrategyRungeKutta4<SparseSpaceType, LocalSpaceType> ExplicitSolvingStrategyRK4Type;

    class AuxiliaryExplicitStrategiesTestElement : public Element
    {
    public:
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AuxiliaryExplicitStrategiesTestElement);

        AuxiliaryExplicitStrategiesTestElement() = default;

        AuxiliaryExplicitStrategiesTestElement(
            IndexType NewId,
            GeometryType::Pointer pGeometry) : Element(NewId, pGeometry) {};

        AuxiliaryExplicitStrategiesTestElement(AuxiliaryExplicitStrategiesTestElement const &rOther) = delete;

        void GetDofList(
            DofsVectorType& rElementalDofList,
            const ProcessInfo& rCurrentProcessInfo) const override
        {
            if (rElementalDofList.size() != 1) {
                rElementalDofList.resize(1);
            }
            rElementalDofList[0] = GetGeometry()[0].pGetDof(TEMPERATURE);
        }

        void CalculateMassMatrix(
            MatrixType& rMassMatrix,
            const ProcessInfo& rCurrentProcessInfo) override
        {
            if (rMassMatrix.size1() != 1 || rMassMatrix.size2() != 1) {
                rMassMatrix.resize(1, 1, false);
            }

            rMassMatrix(0,0) = 1.0;
        }

        void AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo) override
        {
            auto& r_node = GetGeometry()[0];
            const auto aux = r_node.FastGetSolutionStepValue(TEMPERATURE);
            r_node.FastGetSolutionStepValue(REACTION_FLUX) = 37.5 - 3.5 * aux;
        }

    };


    void GenerateTestExplicitStrategiesModelPart(ModelPart& rModelPart)
    {
        // Model part settings
        rModelPart.SetBufferSize(2);
        rModelPart.AddNodalSolutionStepVariable(TEMPERATURE);
        rModelPart.AddNodalSolutionStepVariable(REACTION_FLUX);
        rModelPart.GetNodalSolutionStepVariablesList().AddDof(&TEMPERATURE, &REACTION_FLUX);

        // Create the auxiliary element geometry
        auto p_node = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        rModelPart.AddNode(p_node);
        p_node->AddDof(TEMPERATURE, REACTION_FLUX);

        // Create and add the auxiliary test element
        typename GeometryType::PointsArrayType points_vect;
        points_vect.push_back(p_node);
        auto p_geom = Kratos::make_shared<Point3D<NodeType>>(p_node);
        auto p_elem = Kratos::make_intrusive<AuxiliaryExplicitStrategiesTestElement>(1, p_geom);
        rModelPart.AddElement(p_elem);
    }

    /**
     * Checks if the Linear strategy performs correctly the resolution of the system
     */

    KRATOS_TEST_CASE_IN_SUITE(ExplicitSolvingStrategyRungeKutta4, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart");

        // Set the test model part
        const double delta_time = 1.5; // Set time step
        GenerateTestExplicitStrategiesModelPart(r_model_part); // Create the geometry

        // Create the RK4 explicit strategy
        const bool move_mesh_flag = false;
        const unsigned int rebuild_level = 0;
        auto p_explicit_bs = Kratos::make_shared<ExplicitBuilderType>();
        auto p_explicit_strategy = Kratos::make_unique<ExplicitSolvingStrategyRK4Type>(
            r_model_part,
            p_explicit_bs,
            move_mesh_flag,
            rebuild_level);

        // Solve and check the two first RK4 steps
        p_explicit_strategy->Initialize();
        auto p_test_node = r_model_part.pGetNode(1);
        p_test_node->FastGetSolutionStepValue(TEMPERATURE) = 50.0; // Set initial solution

        // 1st step
        r_model_part.CloneTimeStep(delta_time);
        p_explicit_strategy->InitializeSolutionStep();
        p_explicit_strategy->SolveSolutionStep();
        p_explicit_strategy->FinalizeSolutionStep();
        KRATOS_CHECK_NEAR(p_test_node->FastGetSolutionStepValue(TEMPERATURE), 681.238, 1.0e-3);
        // 2nd step
        r_model_part.CloneTimeStep(2.0 * delta_time);
        p_explicit_strategy->InitializeSolutionStep();
        p_explicit_strategy->SolveSolutionStep();
        p_explicit_strategy->FinalizeSolutionStep();
        KRATOS_CHECK_NEAR(p_test_node->FastGetSolutionStepValue(TEMPERATURE), 11455.1, 1.0e-1);
    }

} // namespace Testing
} // namespace Kratos
