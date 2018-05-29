#ifndef KRATOS_RESIDUALBASED_ELIMINATION_BUILDER_AND_SOLVER_PERIODIC_NORMAL_ONLY_H
#define KRATOS_RESIDUALBASED_ELIMINATION_BUILDER_AND_SOLVER_PERIODIC_NORMAL_ONLY_H

/* System includes */
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "includes/model_part.h"
#include "fluid_rve_lagrange_multipliers_application.h"

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

/** Variant of ResidualBasedEliminationBuilderAndSolver for problems with periodic boundary conditions.
 * @see PeriodicCondition
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedEliminationBuilderAndSolverPeriodicNormalOnly
    : public ResidualBasedEliminationBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedEliminationBuilderAndSolverPeriodicNormalOnly);


    typedef ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

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
    ResidualBasedEliminationBuilderAndSolverPeriodicNormalOnly(typename TLinearSolver::Pointer pNewLinearSystemSolver,
                                               const Variable<int>& PeriodicVariable ):
        ResidualBasedEliminationBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver),
        mPeriodicIdVar(PeriodicVariable)
    {}

    /** Destructor.
     */
    virtual ~ResidualBasedEliminationBuilderAndSolverPeriodicNormalOnly()
    {
    }


    /*@} */
    /**@name Operators */
    /*@{ */


    /*@} */
    /**@name Operations */
    /*@{ */


    virtual void SetUpSystem(ModelPart &r_model_part)
    {
        // Assign an Equation Id to all non-duplicate nodes
        
        KRATOS_WATCH("setting up normal periodic system")
        
        unsigned int EqId = 0;

        for (typename DofsArrayType::iterator itDof = BaseType::mDofSet.begin(); itDof != BaseType::mDofSet.end(); ++itDof)
        {
			if(itDof->GetVariable()==VELOCITY_X)
			{
				if ( ( itDof->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) < static_cast<int>(itDof->Id()) && itDof->GetSolutionStepValue(NODE_PAIR_X_COMPONENT_ANTIPERIODIC)==0 )  || (   itDof->GetSolutionStepValue(NODE_PAIR_X_COMPONENT_ANTIPERIODIC) < static_cast<int>(itDof->Id()) && itDof->GetSolutionStepValue(NODE_PAIR_X_COMPONENT)==0 ) )
					itDof->SetEquationId(EqId++);
			}
			else if(itDof->GetVariable()==VELOCITY_Y)
			{
				if ( ( itDof->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) < static_cast<int>(itDof->Id()) && itDof->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC)==0 ) || ( itDof->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC) < static_cast<int>(itDof->Id()) &&  itDof->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT)==0 ) )
					itDof->SetEquationId(EqId++);
			}
			else if(itDof->GetVariable()==PRESSURE)
			{
				if ( ( itDof->GetSolutionStepValue(NODE_PAIR_PRESSURE) < static_cast<int>(itDof->Id())   ))
					itDof->SetEquationId(EqId++);
			}
			else
			{
				itDof->SetEquationId(EqId++);
			}
        }

        // Copy Equation Id to duplicate nodes.
        for (ModelPart::ConditionIterator itCond = r_model_part.ConditionsBegin(); itCond != r_model_part.ConditionsEnd(); itCond++)
        {
            // PeriodicCondition always have exactly 2 nodes
            ModelPart::ConditionType::GeometryType& rGeom = itCond->GetGeometry();
            if (rGeom.PointsNumber() == 2)
            {
                const int Node0 = rGeom[0].Id();
                const int Node0Pair_x = rGeom[0].FastGetSolutionStepValue(NODE_PAIR_X_COMPONENT);
                const int Node0Pair_y = rGeom[0].FastGetSolutionStepValue(NODE_PAIR_Y_COMPONENT);
                const int Node0Pair_pressure = rGeom[0].FastGetSolutionStepValue(NODE_PAIR_PRESSURE);
                const int Node0Pair_x_antiperiodic = rGeom[0].FastGetSolutionStepValue(NODE_PAIR_X_COMPONENT_ANTIPERIODIC);
                const int Node0Pair_y_antiperiodic = rGeom[0].FastGetSolutionStepValue(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC);

                const int Node1 = rGeom[1].Id();
                const int Node1Pair_x = rGeom[1].FastGetSolutionStepValue(NODE_PAIR_X_COMPONENT);
                const int Node1Pair_y = rGeom[1].FastGetSolutionStepValue(NODE_PAIR_Y_COMPONENT);
                const int Node1Pair_pressure = rGeom[1].FastGetSolutionStepValue(NODE_PAIR_PRESSURE);
				const int Node1Pair_x_antiperiodic = rGeom[1].FastGetSolutionStepValue(NODE_PAIR_X_COMPONENT_ANTIPERIODIC);
                const int Node1Pair_y_antiperiodic = rGeom[1].FastGetSolutionStepValue(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC);

                // If the nodes are marked as a periodic pair (this is to avoid acting on two-noded conditions that are not PeriodicCondition)
                //X periodicity
                if ( ( Node0 == Node1Pair_x ) && ( Node1 == Node0Pair_x ) )
                {
                    if ( Node0 < Node0Pair_x ) // If Node0 is the one with lower Id (the one that does not have an EquationId yet)
                        CopyEquationId_x(rGeom[1],rGeom[0]);
                    else
                        CopyEquationId_x(rGeom[0],rGeom[1]);
                }
                //X anti-periodicity
                if ( ( Node0 == Node1Pair_x_antiperiodic ) && ( Node1 == Node0Pair_x_antiperiodic ) )
                {
                    if ( Node0 < Node0Pair_x_antiperiodic ) // If Node0 is the one with lower Id (the one that does not have an EquationId yet)
                        CopyEquationId_x(rGeom[1],rGeom[0]);
                    else
                        CopyEquationId_x(rGeom[0],rGeom[1]);
                }

                //Y periodicity
                if ( ( Node0 == Node1Pair_y ) && ( Node1 == Node0Pair_y ) )
                {
                    if ( Node0 < Node0Pair_y ) // If Node0 is the one with lower Id (the one that does not have an EquationId yet)
                        CopyEquationId_y(rGeom[1],rGeom[0]);
                    else
                        CopyEquationId_y(rGeom[0],rGeom[1]);
                }
                 //Y anti-periodicity
                if ( ( Node0 == Node1Pair_y_antiperiodic ) && ( Node1 == Node0Pair_y_antiperiodic ) )
                {
                    if ( Node0 < Node0Pair_y_antiperiodic ) // If Node0 is the one with lower Id (the one that does not have an EquationId yet)
                        CopyEquationId_y(rGeom[1],rGeom[0]);
                    else
                        CopyEquationId_y(rGeom[0],rGeom[1]);
                }
                
                if ( ( Node0 == Node1Pair_pressure ) && ( Node1 == Node0Pair_pressure ) )
                {
                    if ( Node0 < Node0Pair_pressure ) // If Node0 is the one with lower Id (the one that does not have an EquationId yet)
                        CopyEquationId_pressure(rGeom[1],rGeom[0]);
                    else
                        CopyEquationId_pressure(rGeom[0],rGeom[1]);
                }
            }
        }

        BaseType::mEquationSystemSize = EqId;
    }

    virtual void ApplyDirichletConditions(typename TSchemeType::Pointer pScheme,
                                          ModelPart &r_model_part,
                                          TSystemMatrixType &A,
                                          TSystemVectorType &Dx,
                                          TSystemVectorType &b)
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
    void CopyEquationId_x(ModelPart::NodeType& rOrigin,
                        ModelPart::NodeType& rDest)
    {
			//copying velocity dof of correct component
			rDest.pGetDof( VELOCITY_X )->SetEquationId( rOrigin.pGetDof(VELOCITY_X)->EquationId() );
			
			//if we are not a corner node, we also copy the pressure dof
            //if(rDest.FastGetSolutionStepValue(NODE_PAIR_Y_COMPONENT)==0)
			//	rDest.pGetDof( PRESSURE )->SetEquationId( PRESSURE );			
    }
    void CopyEquationId_y(ModelPart::NodeType& rOrigin,
                        ModelPart::NodeType& rDest)
    {
			//copying velocity dof of correct component
			rDest.pGetDof( VELOCITY_Y )->SetEquationId(  rOrigin.pGetDof(VELOCITY_Y)->EquationId() );
			//if we are not a corner node, we also copy the pressure dof
            //if(rDest.FastGetSolutionStepValue(NODE_PAIR_X_COMPONENT)==0)
			//	rDest.pGetDof( PRESSURE )->SetEquationId( PRESSURE );			
    }
    void CopyEquationId_pressure(ModelPart::NodeType& rOrigin,
                        ModelPart::NodeType& rDest)
    {
			//copying velocity dof of correct component
			rDest.pGetDof( PRESSURE )->SetEquationId(  rOrigin.pGetDof(PRESSURE)->EquationId() );
			//if we are not a corner node, we also copy the pressure dof
            //if(rDest.FastGetSolutionStepValue(NODE_PAIR_X_COMPONENT)==0)
			//	rDest.pGetDof( PRESSURE )->SetEquationId( PRESSURE );			
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

}; /* Class ResidualBasedEliminationBuilderAndSolverPeriodicNormalOnly */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif // KRATOS_RESIDUALBASED_ELIMINATION_BUILDER_AND_SOLVER_PERIODIC_NORMAL_ONLY_H
