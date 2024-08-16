// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Juan Carlos Alvarado Morones
//

#pragma once

// System includes
#include <filesystem>

// External includes

// Project includes
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
template <class TSparseSpace,
class TDenseSpace,  // = DenseSpace<double>,
class TLinearSolver // = LinearSolver<TSparseSpace,TDenseSpace>
>
class InverseFormingStrategy
: public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(InverseFormingStrategy);

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    ///@}

    ///@name Life Cycle
    ///@{

    /**
    * Constructor.
    */

   // Constructor with Builder and Solver
   InverseFormingStrategy(
    ModelPart& rModelPart,
    typename TSchemeType::Pointer pScheme,
    typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
    typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
    ModelPart& rInverseFormingPart,
    int MaxIterations = 30,
    bool CalculateReactions = false,
    bool ReformDofSetAtEachStep = false,
    bool MoveMeshFlag = false)
   : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme,
   pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag),
   mrInverseFormingPart(rInverseFormingPart)
   {}

   // Destructor
    ~InverseFormingStrategy() = default;

private:
    void UpdateDatabase(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        const bool MoveMesh) override
    {
        // call Update from scheme to apply DofUpdate of solution
        typename TSchemeType::Pointer p_scheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        p_scheme->Update(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), A, Dx, b);

        // maybe displacement needs to be changed/implemented adjusted here so that solution is visible. Maybe relevant 'mrInverseFormingModelPart' in loop
        for (auto& r_node : mrInverseFormingPart.Nodes()) {
            // Updating reference
            const array_1d<double, 3>& disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3>& disp_non_historical = r_node.GetValue(DISPLACEMENT);
            // KRATOS_WATCH(r_node.FastGetSolutionStepValue(DISPLACEMENT))
            // KRATOS_WATCH(r_node.GetValue(DISPLACEMENT))

            disp_non_historical = disp_non_historical + disp;
            // KRATOS_WATCH(r_node.GetValue(DISPLACEMENT))
            // KRATOS_INFO("InverseFormingStrategy") << "Hello, inside UpdateDatabase function" << std::endl;
            // std::cout << "disp: " << disp << std::endl;

            array_1d<double, 3> CurrentCoords;
            CurrentCoords[0] = r_node.X();
            CurrentCoords[1] = r_node.Y();
            CurrentCoords[2] = r_node.Z();
            // std::cout << "CurrentCoords: " << CurrentCoords << std::endl;

            array_1d<double, 3> InitCoords = (CurrentCoords - disp);
            // std::cout << "InitCoords: " << InitCoords << std::endl;
            r_node.SetInitialPosition(InitCoords[0], InitCoords[1], InitCoords[2]);
            // r_node.GetInitialPosition() += disp;
            // double X0_new = r_node.X0() + disp[0];
            // double Y0_new = r_node.Y0() + disp[1];
            // double Z0_new = r_node.Z0() + disp[2];
            // r_node.SetInitialPosition(X0_new, Y0_new, Z0_new);

            // r_node.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
        }
    }

    ModelPart& mrInverseFormingPart;
    
};  /* Class InverseFormingStrategy */
}   /* namespace Kratos. */