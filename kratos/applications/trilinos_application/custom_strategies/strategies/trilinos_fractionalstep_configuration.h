#if !defined(KRATOS_RESIDUALBASED_TRILINOS_FRACTIONALSTEP_CONFIGURATION )
#define  KRATOS_RESIDUALBASED_TRILINOS_FRACTIONALSTEP_CONFIGURATION


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"

#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/builder_and_solvers/trilinos_residualbased_elimination_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_residualbased_elimination_builder_and_solver_componentwise_slip.h"


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
	class TDenseSpace,
	class TLinearSolver
	>
	class TrilinosFractionalStepConfiguration : public SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>
	{
	public:
		/**@name Type Definitions */
		/*@{ */

		/** Counted pointer of ClassName */

		/*@} */
		/**@name Life Cycle
		*/
		/*@{ */

		/** Constructor.
		*/
                TrilinosFractionalStepConfiguration(Epetra_MpiComm&  Comm,
                                            ModelPart& model_part,
                                            typename TLinearSolver::Pointer pNewVelocityLinearSolver,
                                            typename TLinearSolver::Pointer pNewPressureLinearSolver,
                                            unsigned int mDomainSize,
                                            unsigned int laplacian_form,
                                            bool use_dt_in_stabilization
                                            )
                        :  SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, mDomainSize)
                {
                        bool CalculateReactions = false;
                        bool CalculateNormDxFlag = true;
                        bool ReformDofAtEachIteration = false;

                        //definitions for trilinos
                        int guess_row_size;
                        if(mDomainSize == 2) guess_row_size = 15;
                        else guess_row_size = 40;
                        
                        //computation of the fractional vel velocity (first step)
                        //3 dimensional case
			typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;
			typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

                        //initializing fractional velocity solution step
			typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
			typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
				( new TrilinosResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,  TDenseSpace >() );

                        //CONSTRUCTION OF VELOCITY
			typedef typename Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > VarComponent;			
			BuilderSolverTypePointer vel_build =
			BuilderSolverTypePointer(	
			new TrilinosResidualBasedEliminationBuilderAndSolverComponentwiseSplit<TSparseSpace,TDenseSpace,TLinearSolver,VarComponent>(Comm,guess_row_size,pNewVelocityLinearSolver,this->mDomainSize,FRACT_VEL_X,FRACT_VEL_Y,FRACT_VEL_Z) );
			this->mpfracvel_strategy = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >(model_part,pscheme,pNewVelocityLinearSolver,vel_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
			this->mpfracvel_strategy->SetEchoLevel(1);

// 			BuilderSolverTypePointer vel_y_build	= BuilderSolverTypePointer(	new TrilinosResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(Comm,guess_row_size,pNewVelocityLinearSolver) );
// 			this->mpfracvel_y_strategy = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >(model_part,pscheme,pNewVelocityLinearSolver,vel_y_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
// 			this->mpfracvel_y_strategy->SetEchoLevel(1);
//  
// 			if(this->mDomainSize == 3)
// 			{
// 				BuilderSolverTypePointer vel_z_build = BuilderSolverTypePointer(	new TrilinosResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(Comm,guess_row_size,pNewVelocityLinearSolver) );
// 				this->mpfracvel_z_strategy = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >(model_part,pscheme,pNewVelocityLinearSolver,vel_z_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
// 				this->mpfracvel_z_strategy->SetEchoLevel(1);
// 			}
                        
                        
                        
                        //CONSTRUCTION OF PRESSURE SOLVER
			if( laplacian_form == 1) //laplacian form
			{
				std::cout << "standard laplacian form" << std::endl;
//				use_dt_in_stabilization = false;
//				this->mppressurestep = typename BaseType::Pointer(
//					new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >
//					(model_part,pscheme,pNewPressureLinearSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
//				this->mppressurestep->SetEchoLevel(1);

				BuilderSolverTypePointer pressure_builderandsolver = BuilderSolverTypePointer(	new TrilinosResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(Comm,guess_row_size,pNewPressureLinearSolver) );
				this->mppressurestep = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >(model_part,pscheme,pNewVelocityLinearSolver,pressure_builderandsolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
				this->mppressurestep->SetEchoLevel(1);

			}
			else
                        {
                            KRATOS_ERROR(std::logic_error,"trilinos discrete laplacian not yet implemented","")
                        }

                }

		/** Destructor.
		*/

		/*@} */
		/**@name Operators
		*/
		/*@{ */


                typename SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer pGetStrategy(const std::string& strategy_name )
                {
                    KRATOS_TRY

                    if(strategy_name == std::string("vel_strategy"))
                        return mpfracvel_strategy;
//                     else if(strategy_name == std::string("vel_y_strategy"))
//                          return mpfracvel_y_strategy;
//                     else if(strategy_name == std::string("vel_z_strategy"))
//                          return mpfracvel_z_strategy;
                    else if(strategy_name == std::string("pressure_strategy"))
                         return mppressurestep;
                    else
                        KRATOS_ERROR(std::invalid_argument,"trying to get an inexisting strategy","");
                    
                    KRATOS_CATCH("")
                }

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
                typename SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer mpfracvel_strategy;
/*		typename SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer mpfracvel_y_strategy;
		typename SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer mpfracvel_z_strategy;*/
		typename SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer mppressurestep;

                /*@{ */
		//this funcion is needed to ensure that all the memory is allocated correctly


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


		/*@} */

	}; /* Class FractionalStepStrategy */

	/*@} */

	/**@name Type Definitions */
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_TRILINOS_FRACTIONALSTEP_CONFIGURATION  defined */
