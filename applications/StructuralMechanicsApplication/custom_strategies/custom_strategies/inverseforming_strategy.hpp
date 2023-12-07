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
class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
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
    ModelPart& model_part,
    typename TSchemeType::Pointer pScheme,
    typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
    typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
    ModelPart& rInverseFormingModelPart,    // maybe not needed (?)
    const std::string& rPrintingFormat,     // may not be needed, used in formfinding for output files format
    int MaxIterations = 30,
    bool CalculateReactions = false,
    bool ReformDofSetAtEachStep = false,
    bool MoveMeshFlag = false
   )
   : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
   pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag),
   mrInverseFormingModelPart(rInverseFormingModelPart), // maybe not needed (?)
   mPrintingFormat(rPrintingFormat)                     // may not be needed
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
        BaseType::UpdateDatabase(A, Dx, b, MoveMesh);
        for (auto& r_node : mrInverseFormingModelPart.Nodes()) {
            // Updating reference
            const array_1d<double, 3>& disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3>& disp_non_historical = r_node.GetValue(DISPLACEMENT);

            disp_non_historical = disp_non_historical + disp;
            r_node.GetInitialPosition() += disp;

            r_node.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
        }
    }

    ModelPart& mrInverseFormingModelPart;
    std::string mPrintingFormat;
    
};  /* Class InverseFormingStrategy */
}   /* namespace Kratos. */