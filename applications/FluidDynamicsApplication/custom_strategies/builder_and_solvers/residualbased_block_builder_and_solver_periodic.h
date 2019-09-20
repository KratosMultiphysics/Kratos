//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#ifndef KRATOS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_PERIODIC_H
#define KRATOS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_PERIODIC_H

/* System includes */
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "includes/model_part.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/** Variant of ResidualBasedBlockBuilderAndSolver for problems with periodic boundary conditions.
 * @see PeriodicCondition
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedBlockBuilderAndSolverPeriodic
    : public ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverPeriodic);


    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ResidualBasedBlockBuilderAndSolverPeriodic(typename TLinearSolver::Pointer pNewLinearSystemSolver,
                                               const Variable<int>& PeriodicVariable ):
        ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver),
        mPeriodicIdVar(PeriodicVariable)
    {}

    /** Destructor.
     */
    ~ResidualBasedBlockBuilderAndSolverPeriodic() override
    {
    }


    /*@} */
    /**@name Operators */
    /*@{ */


    /*@} */
    /**@name Operations */
    /*@{ */


    void SetUpSystem(ModelPart &r_model_part) override
    {
        // Assign an Equation Id to all non-duplicate nodes
        unsigned int EqId = 0;

        // Modify edge/corner periodic condtions so that one of the nodes gets EquationId
        for (ModelPart::ConditionIterator itCond = r_model_part.ConditionsBegin(); itCond != r_model_part.ConditionsEnd(); itCond++) {
            auto& r_geometry = itCond->GetGeometry();
            if (itCond->Is(PERIODIC) && (r_geometry.PointsNumber() == 4 || r_geometry.PointsNumber() == 8)) {
                r_geometry[0].GetSolutionStepValue(mPeriodicIdVar) = 0;
            }
        }

        for (typename DofsArrayType::iterator itDof = BaseType::mDofSet.begin(); itDof != BaseType::mDofSet.end(); ++itDof)
        {
            if ( itDof->GetSolutionStepValue(mPeriodicIdVar) < static_cast<int>(itDof->Id()) )
                itDof->SetEquationId(EqId++);
        }

        // Reset edge/corner periodic condtions
        for (ModelPart::ConditionIterator itCond = r_model_part.ConditionsBegin(); itCond != r_model_part.ConditionsEnd(); itCond++) {
            auto& r_geometry = itCond->GetGeometry();
            if (itCond->Is(PERIODIC) && (r_geometry.PointsNumber() == 4 || r_geometry.PointsNumber() == 8)) {
                r_geometry[0].GetSolutionStepValue(mPeriodicIdVar) = r_geometry[1].GetSolutionStepValue(mPeriodicIdVar);
            }
        }

        // Copy Equation Id to duplicate nodes.
        for (ModelPart::ConditionIterator itCond = r_model_part.ConditionsBegin(); itCond != r_model_part.ConditionsEnd(); itCond++)
        {
            // PeriodicCondition always have exactly 2 nodes
            ModelPart::ConditionType::GeometryType& rGeom = itCond->GetGeometry();
            if (itCond->Is(PERIODIC)) {
                if (rGeom.PointsNumber() == 2) {
                    int Node0 = rGeom[0].Id();
                    int Node0Pair = rGeom[0].FastGetSolutionStepValue(mPeriodicIdVar);

                    int Node1 = rGeom[1].Id();
                    int Node1Pair = rGeom[1].FastGetSolutionStepValue(mPeriodicIdVar);

                    // If the nodes are marked as a periodic pair (this is to avoid acting on two-noded conditions that are not PeriodicCondition)
                    if ( ( Node0 == Node1Pair ) && ( Node1 == Node0Pair ) ) {
                        if ( Node0 < Node0Pair ) // If Node0 is the one with lower Id (the one that does not have an EquationId yet)
                            CopyEquationId(rGeom[1],rGeom[0]);
                        else
                            CopyEquationId(rGeom[0],rGeom[1]);
                    }
                }
                else {
                    for (unsigned int i = 1; i < rGeom.PointsNumber(); i++) {
                        CopyEquationId(rGeom[0],rGeom[i]);
                    }
                }
            }
        }

        BaseType::mEquationSystemSize = EqId;
    }

    void ApplyDirichletConditions(typename TSchemeType::Pointer pScheme,
                                          ModelPart &r_model_part,
                                          TSystemMatrixType &A,
                                          TSystemVectorType &Dx,
                                          TSystemVectorType &b) override
    {
        double* Avalues = A.value_data().begin();
        std::size_t* Arow_indices = A.index1_data().begin();
        std::size_t* Acol_indices = A.index2_data().begin();

        for (typename DofsArrayType::iterator itDof = BaseType::mDofSet.begin(); itDof != BaseType::mDofSet.end(); ++itDof)
        {
            if (itDof->IsFixed())
            {
                std::size_t RowId = itDof->EquationId();
                std::size_t RowBegin = Arow_indices[RowId];
                std::size_t RowEnd = Arow_indices[RowId+1];

                for (std::size_t k = RowBegin; k != RowEnd; k++)
                {
                    if ( Acol_indices[k] == RowId )
                        Avalues[k] = 1.0;
                    else
                        Avalues[k] = 0.0;
                }

                b[RowId] = 0.0;
            }
        }
    }

    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */



    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    const Variable<int>& mPeriodicIdVar;

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */

    /// Duplicate EquationIds to the second node on each periodic pair
    void CopyEquationId(ModelPart::NodeType& rOrigin,
                        ModelPart::NodeType& rDest)
    {
        for (auto itDof = rOrigin.GetDofs().begin(); itDof != rOrigin.GetDofs().end(); itDof++)
            rDest.pGetDof( itDof->GetVariable() )->SetEquationId( itDof->EquationId() );
    }

    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class ResidualBasedBlockBuilderAndSolverPeriodic */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif // KRATOS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_PERIODIC_H
