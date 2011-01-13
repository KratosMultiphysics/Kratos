/* *********************************************************   
*          
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2007-03-06 10:30:32 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY_NONLINEAR )
#define  KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY_NONLINEAR


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "convection_diffusion_application.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "includes/convection_diffusion_settings.h"


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
	class ResidualBasedConvectionDiffusionStrategyNonLinear 
		: public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
	{
	public:
		/**@name Type Definitions */       
		/*@{ */

		/** Counted pointer of ClassName */
		//typedef boost::shared_ptr< ResidualBasedConvectionDiffusionStrategyNonLinear<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
		KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedConvectionDiffusionStrategyNonLinear );
		//typedef boost::shared_ptr< ResidualBasedConvectionDiffusionStrategyNonLinear<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;

		typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

		typedef typename BaseType::TDataType TDataType;

		//typedef typename BaseType::DofSetType DofSetType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;



		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */

		/** Constructor.
		*/


		ResidualBasedConvectionDiffusionStrategyNonLinear(
			ModelPart& model_part, 
			typename TLinearSolver::Pointer pNewLinearSolver,
			bool ReformDofAtEachIteration = true,
			int time_order = 2,
			int max_iter = 30,
			double toll = 1e-6
			)
			: SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part,false)
		{
			KRATOS_TRY

			mtime_order = time_order;
			mOldDt = 0.00;
			mtoll = toll;
			mmax_iter = max_iter;
			
			ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
			ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
			const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();

			//initializing fractional velocity solution step
			typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
			typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
				( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,  TDenseSpace >() );

			bool CalculateReactions = false;
			bool CalculateNormDxFlag = true;



			typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

			BuilderSolverTypePointer componentwise_build = BuilderSolverTypePointer(new	ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace,TDenseSpace,TLinearSolver,Variable<double> > (pNewLinearSolver,rUnknownVar) );
			mstep1 = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > 				(model_part,pscheme,pNewLinearSolver,componentwise_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
			mstep1->SetEchoLevel(2);

			KRATOS_CATCH("")
		}



		/** Destructor.
		*/
		virtual ~ResidualBasedConvectionDiffusionStrategyNonLinear() {}

		/** Destructor.
		*/

		//*********************************************************************************
		//**********************************************************************
		double Solve()
		{
			KRATOS_TRY

			//calculate the BDF coefficients
			ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
			double Dt = rCurrentProcessInfo[DELTA_TIME];
			
			if(mOldDt == 0.00) //needed for the first step
				mOldDt = Dt;
			if(mtime_order == 2)
			{
				if(BaseType::GetModelPart().GetBufferSize() < 3)
					KRATOS_ERROR(std::logic_error,"insufficient buffer size for BDF2","")
										
				double dt_old = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];
				
				double rho = dt_old/Dt;
				double coeff = 1.0/(Dt*rho*rho+Dt*rho);

				rCurrentProcessInfo[BDF_COEFFICIENTS].resize(3);
				Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
				BDFcoeffs[0] =	coeff*(rho*rho+2.0*rho);	//coefficient for step n+1
				BDFcoeffs[1] =	-coeff*(rho*rho+2.0*rho+1.0);//coefficient for step n
				BDFcoeffs[2] =	coeff;
			}
			else
			{
				rCurrentProcessInfo[BDF_COEFFICIENTS].resize(2);
				Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
				BDFcoeffs[0] =	1.0 / Dt;	//coefficient for step n+1
				BDFcoeffs[1] =	-1.0 / Dt;//coefficient for step n
			}



			unsigned int iter = 0;
			double ratio;
			bool is_converged = false;
			double dT_norm=0.0;
			double T_norm=0.0;

			while(iter++ < mmax_iter && is_converged == false)
			{
				rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
				dT_norm = mstep1->Solve();
				T_norm = CalculateTemperatureNorm();
				CalculateProjection();

				
				ratio = 1.00;
				if( T_norm != 0.00)
					ratio = dT_norm/T_norm;
				else
				{
					std::cout << "T_norm = " << T_norm << " dT_norm = " << dT_norm << std::endl;
				}

				if(dT_norm < 1e-11)
					ratio = 0; //converged
			

				if(ratio < mtoll)
					is_converged = true;

				std::cout <<" iter = " << iter << " ratio = " << ratio <<std::endl;

			}
			


			return dT_norm;
			KRATOS_CATCH("")
		}



		//******************************************************************************************************
		//******************************************************************************************************
		//calculation of projection 
		double CalculateTemperatureNorm()
		{
			KRATOS_TRY;

			double norm = 0.00;

			ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
			ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
			const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();


			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
				{
					norm += pow( i->FastGetSolutionStepValue(rUnknownVar) , 2);
				}

			return sqrt(norm);
							
			KRATOS_CATCH("")
		}

		//******************************************************************************************************
		//******************************************************************************************************
		//calculation of projection 
		void CalculateProjection()
		{
			KRATOS_TRY;

			ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

			//first of all set to zero the nodal variables to be updated nodally
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				if( (i->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0)
				{
					(i)->FastGetSolutionStepValue(TEMP_CONV_PROJ) = 0.00;
					(i)->FastGetSolutionStepValue(NODAL_AREA) = 0.00;
				}
				else
				{
					(i)->FastGetSolutionStepValue(NODAL_AREA) = 1.00;
				}
			}

			//add the elemental contributions for the calculation of the velocity
			//and the determination of the nodal area
			rCurrentProcessInfo[FRACTIONAL_STEP] = 2;
			for(ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin() ; 
				i != BaseType::GetModelPart().ElementsEnd() ; ++i)
			{
				(i)->InitializeSolutionStep(rCurrentProcessInfo);
			}

			//solve nodally for the velocity
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				double& conv_proj = (i)->FastGetSolutionStepValue(TEMP_CONV_PROJ);
				double temp = 1.00 / (i)->FastGetSolutionStepValue(NODAL_AREA);
				conv_proj *= temp; 
			}

			KRATOS_CATCH("")
		}

		virtual void SetEchoLevel(int Level) 
		{
			mstep1->SetEchoLevel(Level);
		}

		virtual void Clear() 
		{
			mstep1->Clear();
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
		typename BaseType::Pointer mstep1;
		double mOldDt;
		int mtime_order;
		unsigned int mmax_iter;
		double mtoll;





		/*@} */
		/**@name Private Operators*/
		/*@{ */

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
		ResidualBasedConvectionDiffusionStrategyNonLinear(const ResidualBasedConvectionDiffusionStrategyNonLinear& Other);


		/*@} */   

	}; /* Class ResidualBasedConvectionDiffusionStrategyNonLinear */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY_NONLINEAR  defined */

