//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author                      $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_BASED_NEWTON_RAPHSON_LINE_SEARCH_IMPLEX_STRATEGY )
#define  KRATOS_RESIDUAL_BASED_NEWTON_RAPHSON_LINE_SEARCH_IMPLEX_STRATEGY

/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/strategies/residual_based_newton_raphson_line_search_strategy.hpp"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "custom_utilities/line_search_calculation_utilities.hpp"  // is the solid mechanics and not the Pfem

#include "solid_mechanics_application_variables.h"

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

   /// Short class definition.

   /**   Detail class definition.

     \URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

     \URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

     \URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

     \URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


     \URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

     \URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

     \URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

     \URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


    */
   template<class TSparseSpace,
      class TDenseSpace, // = DenseSpace<double>,
      class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
         class ResidualBasedNewtonRaphsonLineSearchImplexStrategy
         : public ResidualBasedNewtonRaphsonLineSearchStrategy < TSparseSpace, TDenseSpace, TLinearSolver>
           //: public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
         {
            public:
               /**@name Type Definitions */
               /*@{ */
               typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

               /** Counted pointer of ClassName */

	       KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedNewtonRaphsonLineSearchImplexStrategy );

               typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

               typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

               typedef typename BaseType::TDataType TDataType;

               typedef TSparseSpace SparseSpaceType;

               typedef typename BaseType::TSchemeType TSchemeType;

               //typedef typename BaseType::DofSetType DofSetType;

               typedef typename BaseType::DofsArrayType DofsArrayType;

               typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

               typedef typename BaseType::TSystemVectorType TSystemVectorType;

               typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

               typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

               typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
               typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;


               /*@} */
               /**@name Life Cycle
                */
               /*@{ */

               /** Constructor.
                */
               ResidualBasedNewtonRaphsonLineSearchImplexStrategy(
                     ModelPart& model_part,
                     typename TSchemeType::Pointer pScheme,
                     typename TLinearSolver::Pointer pNewLinearSolver,
                     typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                     int MaxIterations = 30,
                     bool CalculateReactions = false,
                     bool ReformDofSetAtEachStep = false,
                     bool MoveMeshFlag = false
                     )
                  //: SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
                  : ResidualBasedNewtonRaphsonLineSearchStrategy < TSparseSpace, TDenseSpace, TLinearSolver> ( model_part, pScheme, pNewLinearSolver, 
                        pNewConvergenceCriteria, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
            {
               KRATOS_TRY

      KRATOS_CATCH( "" )
            }

               ResidualBasedNewtonRaphsonLineSearchImplexStrategy(
                     ModelPart& model_part,
                     typename TSchemeType::Pointer pScheme,
                     typename TLinearSolver::Pointer pNewLinearSolver,
                     typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                     typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                     int MaxIterations = 30,
                     bool CalculateReactions = false,
                     bool ReformDofSetAtEachStep = false,
                     bool MoveMeshFlag = false
                     )
                  //: SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag)
                  : ResidualBasedNewtonRaphsonLineSearchStrategy < TSparseSpace, TDenseSpace, TLinearSolver> ( model_part, pScheme, pNewLinearSolver, 
                        pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
            {
               KRATOS_TRY

      KRATOS_CATCH( "" )
            }

               /** Destructor.
                */
               virtual ~ResidualBasedNewtonRaphsonLineSearchImplexStrategy()
               {
               }


               //Set and Get Scheme ... containing Builder, Update and other
               virtual bool  GetImplexSetToConstitutiveEquations()
               {
                  if ( mImplexFlag == true) {
                     if ( this->GetModelPart().GetProcessInfo()[IMPLEX] == 1) {
                        return true;
                     }
                  }
                  return false; 
               }

               //**********************************************************************
               //**********************************************************************

               void InitializeSolutionStep()
               {
                  KRATOS_TRY

      mImplexFlag = true;
                  if ( mImplexFlag == true) {
                     this->GetModelPart().GetProcessInfo()[IMPLEX] = 1;
                  }

                  ResidualBasedNewtonRaphsonLineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::InitializeSolutionStep();

                  KRATOS_CATCH( "" )

               }


               virtual void FinalizeSolutionStep( )
               {
                  KRATOS_TRY

      if ( mImplexFlag == true) {
         this->GetModelPart().GetProcessInfo()[IMPLEX] = 0;
      }

                  if ( ResidualBasedNewtonRaphsonLineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::GetCalculateReactionsFlag() == true)
                     //if ( this->GetCalculateReactionsFlag() == true) 
                  {
                     typename TSchemeType::Pointer pScheme = this->GetScheme();
                     typename TBuilderAndSolverType::Pointer pBuilderAndSolver = this->GetBuilderAndSolver();
                     TSystemMatrixType& mA  = this->GetSystemMatrix();
                     TSystemVectorType& mDx = this->GetSystemDx();
                     TSystemVectorType& mb  = this->GetSystemb();

                     (pBuilderAndSolver)->CalculateReactions( pScheme, BaseType::GetModelPart(), mA, mDx, mb);
                  }

                  ResidualBasedNewtonRaphsonLineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::FinalizeSolutionStep();


                  KRATOS_CATCH( "" )
               }

               /*@} */
               /**@name Operators
                */
               /*@{ */

               /*@} */
               /**@name Operations */
               /*@{ */


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

            private:
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

            protected:
               /**@name Static Member Variables */
               /*@{ */


               /*@} */
               /**@name Member Variables */
               /*@{ */


               bool mImplexFlag; // ho he posat perquè sóc així però només entra aquí si hi ha implex... :'P
               //flag to allow to not finalize the solution step, so the historical variables are not updated



               /*@} */
               /**@name Private Operators*/
               /*@{ */
               //**********************************************************************
               //**********************************************************************


               /*@} */
               /**@name Private Operations*/
               /*@{ */


               /*@} */
               /**@name Private  Access */
               /*@{ */


               /*@} */
               /**@name Private Inquiry */
               /*@{ */


               /*@} */
               /**@name Un accessible methods */
               /*@{ */

               /** Copy constructor.
                */
               ResidualBasedNewtonRaphsonLineSearchImplexStrategy(const ResidualBasedNewtonRaphsonLineSearchImplexStrategy& Other)
               {
               };


               /*@} */

         }; /* Class ResidualBasedNewtonRaphsonLineSearchStrategy */

   /*@} */

   /**@name Type Definitions */
   /*@{ */


   /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_NEWTON_RAPHSON_LINE_SEARCH_IMPLEX_STRATEGY  defined */

