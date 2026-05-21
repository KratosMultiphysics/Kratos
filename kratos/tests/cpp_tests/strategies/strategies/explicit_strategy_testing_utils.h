//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla, Eduard Gómez
//

#if !defined(KRATOS_EXPLICIT_STRATEGY_TEST_UTILS_H_INCLUDED)
#define KRATOS_EXPLICIT_STRATEGY_TEST_UTILS_H_INCLUDED

// System includes

/* External includes */

/* Project includes */
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"
#include "testing/testing.h"
#include "containers/model.h"
#include "geometries/point_3d.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "tests/test_utilities/test_element.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/explicit_solving_strategy.h"
#include "utilities/math_utils.h"

namespace Kratos
{

namespace ExplicitStrategyTestingUtils
{

typedef Node NodeType;
typedef Geometry<NodeType> GeometryType;
typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
typedef ExplicitBuilder< SparseSpaceType, LocalSpaceType > ExplicitBuilderType;

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

    void CalculateLumpedMassVector(
        VectorType& rLumpedMassVector,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        rLumpedMassVector = VectorType(1);
        rLumpedMassVector(0) = 1.0;
    }

    void AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo) override
    {
        auto& r_node = GetGeometry()[0];
        const auto aux = r_node.FastGetSolutionStepValue(TEMPERATURE);
        r_node.FastGetSolutionStepValue(REACTION_FLUX) = 37.5 - 3.5 * aux;
    }

    void EquationIdVector(
        EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        if (rEquationIdVector.size() != 1) {
            rEquationIdVector.resize(1);
        }
        rEquationIdVector[0] = GetGeometry()[0].GetDof(TEMPERATURE).EquationId();
    }

};

// Fills the modelpart with nodes and dofs
inline void GenerateTestExplicitStrategiesModelPart(ModelPart& rModelPart)
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

// Solves the problem with the given explicit strategy
template<typename TStrategyType>
inline void Solve(
    ModelPart& rModelPart,
    Node& rTestNode,
    const double time,
    const unsigned int n_steps)
{
    KRATOS_TRY

    rModelPart.CloneTimeStep(0.0);

    // Create the explicit strategy
    const bool move_mesh_flag = false;
    const unsigned int rebuild_level = 0;
    auto p_explicit_bs = Kratos::make_shared<ExplicitBuilderType>();
    auto p_explicit_strategy = Kratos::make_unique<TStrategyType>(
        rModelPart,
        p_explicit_bs,
        move_mesh_flag,
        rebuild_level);

    p_explicit_strategy->Initialize();
    rTestNode.FastGetSolutionStepValue(TEMPERATURE) = 50.0; // Set initial solution

    const double delta_time = time / n_steps;
    for(unsigned int step=1; step<=n_steps; ++step)
    {
        rModelPart.CloneTimeStep(step * delta_time);
        p_explicit_strategy->InitializeSolutionStep();
        p_explicit_strategy->SolveSolutionStep();
        p_explicit_strategy->FinalizeSolutionStep();
    }

    KRATOS_CATCH("");
}

/**
 * Fits a function to a set of value pairs (x, y):
 *
 *  log(y) ~ p0 + p1*log(x)
 *
 * it fits a linear polynomial with coefficients p0 and p1
 * with least squares.
 *
 * @param rXValues The x values
 * @param rYValues The y values
 *
 * @return The slope of the fitted polynomial
 */

template<unsigned int TDataSize>
inline double LogFittingSlope(
    const std::array<double, TDataSize>& rX,
    const std::array<double, TDataSize>& rY)
{
    BoundedMatrix<double, TDataSize, 2> A;
    array_1d<double, TDataSize> y;

    for(unsigned int i=0; i<TDataSize; ++i){
        A(i, 0) = 1.0;
        A(i, 1) = std::log(rX[i]);
        y(i)    = std::log(rY[i]);
    }

    const BoundedMatrix<double, 2, 2> AtA = prod(trans(A), A);
    const array_1d<double, 2> AtY = prod(trans(A), y);

    Vector p;
    MathUtils<double>::Solve(AtA, p, AtY);

    return p[1]; // slope
}

/**
 * Runs a whole problem with the given strategy and
 * checks that the deviation against the analytical solution
 * stays below a tolerance.
 */
template<typename TStrategyType>
inline void RunTest(const double tolerance)
{
    KRATOS_TRY

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart");

    // Set the test model part
    constexpr double delta_time = 0.01; // Set time step. Must be small to preserve stability for low order methods
    GenerateTestExplicitStrategiesModelPart(r_model_part); // Create the geometry

    // Create the explicit strategy
    const bool move_mesh_flag = false;
    const unsigned int rebuild_level = 0;
    auto p_explicit_bs = Kratos::make_shared<ExplicitBuilderType>();
    auto p_explicit_strategy = Kratos::make_unique<TStrategyType>(
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

    double analytical_1 = (37.5 / 3.5) + (50 - 37.5/3.5) * std::exp(- 3.5 * delta_time);
    KRATOS_EXPECT_NEAR(p_test_node->FastGetSolutionStepValue(TEMPERATURE), analytical_1, tolerance);
    // 2nd step
    r_model_part.CloneTimeStep(2.0 * delta_time);
    p_explicit_strategy->InitializeSolutionStep();
    p_explicit_strategy->SolveSolutionStep();
    p_explicit_strategy->FinalizeSolutionStep();

    double analytical_2 = (37.5 / 3.5) + (50 - 37.5/3.5) * std::exp(- 3.5 * 2.0 * delta_time);
    KRATOS_EXPECT_NEAR(p_test_node->FastGetSolutionStepValue(TEMPERATURE), analytical_2, tolerance);

    KRATOS_CATCH("");
}

/**
 * Runs a whole problem with the given strategy and
 * checks that the deviation against the analytical solution.
 *
 * This process is repeated with diferent time-step sizes in order
 * to determine the convergence rate of the explicit strategy. This result
 * is then checked againts the expected one.
 */
template<typename TStrategyType>
inline void ConvergenceTest(const unsigned int ExpectedOrder)
{
    KRATOS_TRY

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart");

    // Set the test model part
    GenerateTestExplicitStrategiesModelPart(r_model_part); // Create the geometry
    auto& r_test_node = r_model_part.GetNode(1);

    constexpr double total_time = 1.0;

    constexpr unsigned int n_testpoints = 4;
    std::array<unsigned int, n_testpoints> n_steps = {32, 46, 68, 100}; // evenly distributed (in logspace) bewtween 10^1.5 and 10^2.0
    std::array<double, n_testpoints> delta_time;
    std::array<double, n_testpoints> error;

    for(unsigned int i=0; i<n_testpoints; ++i)
    {
        Solve<TStrategyType>(r_model_part, r_test_node, total_time, n_steps[i]);
        delta_time[i] = total_time / n_steps[i];

        const double result = r_test_node.FastGetSolutionStepValue(TEMPERATURE);
        const double analytical = (37.5 / 3.5) + (50 - 37.5/3.5) * std::exp(- 3.5 * n_steps[i] * delta_time[i]);

        error[i] = std::abs(result - analytical);
    }

    const double convergence_rate = LogFittingSlope<n_testpoints>(delta_time, error);

    KRATOS_EXPECT_NEAR(convergence_rate, ExpectedOrder, 0.1);

    KRATOS_CATCH("");
}

} // namespace ExplicitStrategyTestingUtils
} // namespace Kratos

#endif // KRATOS_EXPLICIT_STRATEGY_TEST_UTILS_H_INCLUDED