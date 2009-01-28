/* *********************************************************   
*          
*   Last Modified by:    $Author: antonia $
*   Date:                $Date: 2008-05-13 14:09:23 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_ND_FLUID_STRATEGY2 )
#define  KRATOS_RESIDUALBASED_ND_FLUID_STRATEGY2


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "incompressible_fluid_application.h"

#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "custom_strategies/builder_and_solvers/residualbased_elimination_ND_discretelaplacian_builder_and_solver.h"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"


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
	class ResidualBasedNDFluidStrategy 
		: public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
	{
	public:
		/**@name Type Defintubo_12DL.flavia.resitions */       
		/*@{ */

		/** Counted pointer of ClassName */
		typedef boost::shared_ptr< ResidualBasedNDFluidStrategy<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;

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
		ResidualBasedNDFluidStrategy(
			ModelPart& model_part, 
			typename TLinearSolver::Pointer pNewVelocityLinearSolver,
			typename TLinearSolver::Pointer pNewPressureLinearSolver,
			bool CalculateReactions = false,
			bool ReformDofAtEachIteration = true,
			bool CalculateNormDxFlag = true,
			double velocity_toll = 0.01,
			double pressure_toll = 0.01,
			int MaxVelocityIterations = 3,
			int MaxPressureIterations = 1,
			unsigned int time_order  = 2,
			//unsigned int prediction_order  = 2,
			unsigned int domain_size = 2,
			unsigned int laplacian_form = 2, //1 = laplacian, 2 = discrete laplacian
			bool predictor_corrector = false
			)
			: SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part,false)
		{
			KRATOS_TRY
//std::cout << "SONO QUI" << std::endl;
			this->mvelocity_toll = velocity_toll;
			this->mpressure_toll = pressure_toll;
			this->mMaxVelIterations = MaxVelocityIterations;
			this->mMaxPressIterations = MaxPressureIterations;
			this->mtime_order = time_order;
			this->mprediction_order = time_order;
			this->mdomain_size = domain_size;
			this->mpredictor_corrector = predictor_corrector;
			this->mReformDofAtEachIteration = ReformDofAtEachIteration;
			this->proj_is_initialized = false;

			//the system will be cleared at the end!
			ReformDofAtEachIteration = false;

				//initializing fractional velocity solution step
				typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
			typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
				( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,  TDenseSpace >() );

			//bool CalculateReactions = false;
			//bool ReformDofAtEachIteration = true;
			//bool CalculateNormDxFlag = true;

			//computation of the fractional vel velocity (first step)
				//3 dimensional case
			typedef typename Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > VarComponent;

			typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

			BuilderSolverTypePointer vel_x_build =BuilderSolverTypePointer(				new	ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace,TDenseSpace,TLinearSolver, VarComponent>(pNewVelocityLinearSolver,FRACT_VEL_X) );
			
			this->mfracvel_x_strategy = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >				(model_part,pscheme,pNewVelocityLinearSolver,vel_x_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
			this->mfracvel_x_strategy->SetEchoLevel(2);

			BuilderSolverTypePointer vel_y_build	= BuilderSolverTypePointer(new	ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace,TDenseSpace,TLinearSolver,VarComponent> (pNewVelocityLinearSolver,FRACT_VEL_Y) );
			this->mfracvel_y_strategy = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > 				(model_part,pscheme,pNewVelocityLinearSolver,vel_y_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
			this->mfracvel_y_strategy->SetEchoLevel(2);
 
			if(this->mdomain_size == 3)
			{
				BuilderSolverTypePointer vel_z_build	= BuilderSolverTypePointer(
					new	ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace,TDenseSpace,TLinearSolver,VarComponent>(pNewVelocityLinearSolver,FRACT_VEL_Z) );
				this->mfracvel_z_strategy = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >					(model_part,pscheme,pNewVelocityLinearSolver,vel_z_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
				this->mfracvel_z_strategy->SetEchoLevel(2);

			}

			if( laplacian_form == 1) //laplacian form
			{
				std::cout << "standard laplacian form" << std::endl;
				
				this->mpressurestep = typename BaseType::Pointer(
					new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >
					(model_part,pscheme,pNewPressureLinearSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
				this->mpressurestep->SetEchoLevel(2); 
			}
			else //discrete laplacian form
			{
				std::cout << "discrete laplacian form" << std::endl;
				BuilderSolverTypePointer discretebuild;
				if(domain_size == 2)
				{
				//2 dimensional case
				discretebuild = BuilderSolverTypePointer(
					new	ResidualBasedEliminationNDDiscreteLaplacianBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver, 2>(pNewPressureLinearSolver)
					);
				}
				else if (domain_size == 3)
				{
				//3 dimensional case
				discretebuild = BuilderSolverTypePointer(
					new	ResidualBasedEliminationNDDiscreteLaplacianBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver, 3>(pNewPressureLinearSolver)
					);
				}

				this->mpressurestep = typename BaseType::Pointer(
					new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >
					(model_part,pscheme,pNewPressureLinearSolver,discretebuild,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
				this->mpressurestep->SetEchoLevel(2);

			}

			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				if(i->IsFixed(VELOCITY_X))
					(i)->Fix(FRACT_VEL_X);
				if(i->IsFixed(VELOCITY_Y))
					(i)->Fix(FRACT_VEL_Y);
				if(i->IsFixed(VELOCITY_Z))
					(i)->Fix(FRACT_VEL_Z);
			}

			this->m_step = 1;

			this->mOldDt  =0.00;

			KRATOS_CATCH("")
		}



		/** Destructor.
		*/
		virtual ~ResidualBasedNDFluidStrategy() {}

		/** Destructor.
		*/

		//*********************************************************************************
		//**********************************************************************
		double Solve()
		{
			KRATOS_TRY

			//assign the correct fractional step coefficients (BDF_COEFFICIENTS..)
			InitializeFractionalStep(this->m_step, this->mtime_order);
			double Dp_norm;

			//predicting the velocity 
			PredictVelocity(this->m_step,this->mprediction_order);

			//initialize projections at the first steps
			InitializeProjections(this->m_step, this->proj_is_initialized);

			//Assign Velocity To Fract Step Velocity and Node Area to Zero
			AssignInitialStepValues();

			if(this->m_step <= this->mtime_order)
				Dp_norm = IterativeSolve();
			else
			{
				if(this->mpredictor_corrector == false) //standard fractional step
					Dp_norm = FracStepSolution();
				else  //iterative solution
					Dp_norm = IterativeSolve();
			}
			
			if(this->mReformDofAtEachIteration == true )
				this->Clear();

			this->m_step += 1;
			this->mOldDt =  BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME];

			return Dp_norm;
			KRATOS_CATCH("")
		}


		//*********************************************************************************
		//**********************************************************************

		double IterativeSolve()
		{
			KRATOS_TRY


			double Dp_norm = 1.00;
			int iteration = 0;
 
			int MaxPressureIterations = this->mMaxPressIterations;

			while(	Dp_norm >= this->mpressure_toll && iteration++ < MaxPressureIterations  ) 
			{
/*				double p_norm = 0.0;
				for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
						i != BaseType::GetModelPart().NodesEnd() ; ++i)
				{
					//setting the old value of the pressure to the current one
					const double& p = (i)->FastGetSolutionStepValue(PRESSURE);
					(i)->FastGetSolutionStepValue(PRESSURE,1) = p;

					p_norm+=p*p;
				}	
				p_norm = sqrt(p_norm);*/
				double p_norm = SavePressureIteration();				

				Dp_norm = FracStepSolution();

				Dp_norm /= p_norm;

				std::cout << "it = " << iteration << " Pressure Variation Norm = " << Dp_norm << std::endl;

			}

			if(this->mReformDofAtEachIteration == true )
				this->Clear();

			return Dp_norm;
			KRATOS_CATCH("")
		}


		//*********************************************************************************
		//**********************************************************************

		double SavePressureIteration()
		{
			KRATOS_TRY


			double p_norm = 0.0;
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				//setting the old value of the pressure to the current one
				const double& p = (i)->FastGetSolutionStepValue(PRESSURE);
				//(i)->FastGetSolutionStepValue(PRESSURE,1) = p;
				(i)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = p;

				p_norm+=p*p;
			}	
			p_norm = sqrt(p_norm);
			KRATOS_WATCH(p_norm);	

			return p_norm;
			KRATOS_CATCH("")
		}

		//******************************************************************************************************
		//******************************************************************************************************
		//explicit correction for velocities and eventually stabilization terms
		double FracStepSolution()
		{
			KRATOS_TRY
				
			//setting the fractional velocity to the value of the velocity
			AssignInitialStepValues();

			//solve first step for fractional step velocities
			boost::timer step1time;
			this->SolveStep1(this->mvelocity_toll, this->mMaxVelIterations);
			std::cout << "step1 time " << step1time.elapsed() << std::endl;

		
// 			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
// 				i != BaseType::GetModelPart().NodesEnd() ; ++i)
// 			{
// 				std::cout << "fractional velocity " << i->Id()<< " " << i->FastGetSolutionStepValue(FRACT_VEL) << std::endl;
// 			}

			//solve for pressures (and recalculate the nodal area)
			boost::timer step2time;
			double Dp_norm = this->SolveStep2();
			std::cout << "pressure calculation time " << step2time.elapsed() << std::endl;

			this->ActOnLonelyNodes();

			//calculate projection terms
			boost::timer projection_time;
			this->SolveStep3();
			std::cout << "projection calculation time " << projection_time.elapsed() << std::endl;

			//correct velocities
			boost::timer vel_time;
			this->SolveStep4();
			std::cout << "velocity correction time " << vel_time.elapsed() << std::endl;
	
			

			return Dp_norm;

			KRATOS_CATCH("")
		}

		//******************************************************************************************************
		//******************************************************************************************************
		//explicit correction for velocities and eventually stabilization terms
		void SolveStep4()
		{
			KRATOS_TRY;

			ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
			array_1d<double,3> zero = ZeroVector(3);
			Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

			//first of all set to zero the nodal variables to be updated nodally
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
 				double eps = (i)->FastGetSolutionStepValue(POROSITY);
				array_1d<double,3>& fract_v = (i)->FastGetSolutionStepValue(FRACT_VEL);
				fract_v *= (i)->FastGetSolutionStepValue(NODAL_MASS) * BDFcoeffs[0];
				fract_v /= eps;
			}

			//add the elemental contributions for the calculation of the velocity
			//and the determination of the nodal area
			rCurrentProcessInfo[FRACTIONAL_STEP] = 6;
			for(ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin() ; 
				i != BaseType::GetModelPart().ElementsEnd() ; ++i)
			{
				(i)->InitializeSolutionStep(rCurrentProcessInfo);
			}



			//solve nodally for the velocity
			if(this->mdomain_size == 2)
			{
				for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
					i != BaseType::GetModelPart().NodesEnd() ; ++i)
				{
					array_1d<double,3>& v = (i)->FastGetSolutionStepValue(VELOCITY);
					const array_1d<double,3>& fract_v = (i)->FastGetSolutionStepValue(FRACT_VEL);
					double A = (i)->FastGetSolutionStepValue(NODAL_MASS);
 					double eps = (i)->FastGetSolutionStepValue(POROSITY);
					
					double temp = (1.0/BDFcoeffs[0]) / A;
					temp *= eps;
					if(!i->IsFixed(VELOCITY_X))
					{
						v[0] = fract_v[0]*temp;
					}
					if(!i->IsFixed(VELOCITY_Y))
					{
						v[1] = fract_v[1]*temp;
					}				
				}
			}
			else if (this->mdomain_size == 3)
			{
				for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
					i != BaseType::GetModelPart().NodesEnd() ; ++i)
				{
					array_1d<double,3>& v = (i)->FastGetSolutionStepValue(VELOCITY);
					const array_1d<double,3>& fract_v = (i)->FastGetSolutionStepValue(FRACT_VEL);
					double A = (i)->FastGetSolutionStepValue(NODAL_MASS);
 					double eps =(i)->FastGetSolutionStepValue(POROSITY);	

					double temp = (1.0/BDFcoeffs[0]) / A;
					temp *= eps;
					if(!i->IsFixed(VELOCITY_X))
					{
						v[0] = fract_v[0]*temp;
					}
					if(!i->IsFixed(VELOCITY_Y))
					{
						v[1] = fract_v[1]*temp;
					}				
					if(!i->IsFixed(VELOCITY_Z))
					{
						v[2] = fract_v[2]*temp;
					}		
				}			
			}



			KRATOS_CATCH("");
		}

		//******************************************************************************************************
		//******************************************************************************************************
		//sets the BDF coefficients to the correct value
		void InitializeFractionalStep(const int step,const int time_order)
		{
			KRATOS_TRY;

			//calculate the BDF coefficients
			ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
			double Dt = rCurrentProcessInfo[DELTA_TIME];
			
			if(this->mOldDt == 0.00) //needed for the first step
				this->mOldDt = Dt;

			if(time_order == 2 && step > time_order)
			{
				if(BaseType::GetModelPart().GetBufferSize() < 3)
					KRATOS_ERROR(std::logic_error,"insufficient buffer size for BDF2","")

				rCurrentProcessInfo[BDF_COEFFICIENTS].resize(3,false);
				Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
				BDFcoeffs[0] =	1.5 / Dt;	//coefficient for step n+1
				BDFcoeffs[1] =	-2.0 / Dt;//coefficient for step n
				BDFcoeffs[2] =	0.5 / Dt;//coefficient for step n-1
			}
			else
			{
				rCurrentProcessInfo[BDF_COEFFICIENTS].resize(2,false);
				Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
				BDFcoeffs[0] =	1.0 / Dt;	//coefficient for step n+1
				BDFcoeffs[1] =	-1.0 / Dt;//coefficient for step n
			}

			
			KRATOS_CATCH("");
		}

		//******************************************************************************************************
		//******************************************************************************************************
		void InitializeProjections(int step, bool proj_is_initialized)
		{
			if(step <= 2 && proj_is_initialized==false) 
			{ 
				for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
					i != BaseType::GetModelPart().NodesEnd() ; ++i)
				{
					noalias( i->FastGetSolutionStepValue(PRESS_PROJ) ) = i->FastGetSolutionStepValue(BODY_FORCE);
				}
				proj_is_initialized = true;
			}
		}

		//******************************************************************************************************
		//******************************************************************************************************
		//predict values for the fractional step velocities
		//and set to zero the nodal mass
		void AssignInitialStepValues()
		{
			KRATOS_TRY
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				//predicting the values for the fluid velocity
				array_1d<double,3>& v = (i)->FastGetSolutionStepValue(VELOCITY);
				array_1d<double,3>& fracv = (i)->FastGetSolutionStepValue(FRACT_VEL);
				noalias(fracv) = v;

				//setting the old pressure iteration to the value of the pressure
				(i)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = (i)->FastGetSolutionStepValue(PRESSURE);;


				//resetting the nodal area
				double area = 0.00;
				(i)->FastGetSolutionStepValue(NODAL_MASS) = area;
			}
			KRATOS_CATCH("");
		}

		//******************************************************************************************************
		//******************************************************************************************************
		void PredictVelocity(int step,int prediction_order)
		{
			KRATOS_TRY
			if(prediction_order == 2)
			{
				if(BaseType::GetModelPart().GetBufferSize() < 3)
					KRATOS_ERROR(std::logic_error,"insufficient buffer size for second order prediction","")
			}

			if(prediction_order == 2 && step > 2)
			{
				//second order prediction for the velocity
				for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
				{
					array_1d<double,3>& vel = i->FastGetSolutionStepValue(VELOCITY); 
					const array_1d<double,3>& v1 = i->FastGetSolutionStepValue(VELOCITY,1);
					const array_1d<double,3>& v2 = i->FastGetSolutionStepValue(VELOCITY,2);
					if(!i->IsFixed(VELOCITY_X))
						vel[0] = 2.00*v1[0] - v2[0];
					if(!i->IsFixed(VELOCITY_Y))
						vel[1] = 2.00*v1[1] - v2[1];
					if(!i->IsFixed(VELOCITY_Z))
						vel[2] = 2.00*v1[2] - v2[2];
				}
			}
			KRATOS_CATCH("");
		}
		//******************************************************************************************************
		//******************************************************************************************************
		//calculation of projection 
		void SolveStep1(double velocity_toll, int MaxIterations)
		{
			KRATOS_TRY;
			array_1d<double,3> normDx = ZeroVector(3);

			bool is_converged = false;
			double iteration = 1;

			//solve for fractional step velocities
			while(	is_converged == false && iteration++<MaxIterations  ) 
			{
				//perform one iteration over the fractional step velocity
				FractionalVelocityIteration(normDx);

				is_converged = ConvergenceCheck(normDx,velocity_toll);
			}

			if (is_converged == false)
				std::cout << "ATTENTION: convergence NOT achieved" << std::endl;

			//clear if needed
			if(mReformDofAtEachIteration == true && mpredictor_corrector == false )
			{
				this->mfracvel_x_strategy->Clear();
				this->mfracvel_y_strategy->Clear();
				if(this->mdomain_size == 3)
					this->mfracvel_z_strategy->Clear();
			}

			KRATOS_CATCH("");
		}

		//******************************************************************************************************
		//******************************************************************************************************
		double SolveStep2() 
		{
			KRATOS_TRY;
			 BaseType::GetModelPart().GetProcessInfo()[FRACTIONAL_STEP] = 4;
			return mpressurestep->Solve();
			KRATOS_CATCH("");
		}

		//******************************************************************************************************
		//******************************************************************************************************
		//calculation of projection 
		void SolveStep3()
		{
			KRATOS_TRY;

			ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
			array_1d<double,3> zero = ZeroVector(3);

			//first of all set to zero the nodal variables to be updated nodally
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				array_1d<double,3>& press_proj = (i)->FastGetSolutionStepValue(PRESS_PROJ);
				array_1d<double,3>& conv_proj = (i)->FastGetSolutionStepValue(CONV_PROJ);
				noalias(press_proj) = zero;
				noalias(conv_proj) = zero;
			}

			//add the elemental contributions for the calculation of the velocity
			//and the determination of the nodal area
			rCurrentProcessInfo[FRACTIONAL_STEP] = 5;
			for(ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin() ; 
				i != BaseType::GetModelPart().ElementsEnd() ; ++i)
			{
				(i)->InitializeSolutionStep(rCurrentProcessInfo);
			}

			//solve nodally for the velocity
			double temp;
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				array_1d<double,3>& press_proj = (i)->FastGetSolutionStepValue(PRESS_PROJ);
				array_1d<double,3>& conv_proj = (i)->FastGetSolutionStepValue(CONV_PROJ);
				double A = (i)->FastGetSolutionStepValue(NODAL_MASS);

				temp = 1.00 / A;
				press_proj *= temp; 
				conv_proj *= temp; 
				
			}

			KRATOS_CATCH("");
		}

		//******************************************************************************************************
		//******************************************************************************************************
		//calculation of projection 
		void ActOnLonelyNodes()
		{
			KRATOS_TRY;

			//array_1d<double,3> body_force = BaseType::GetModelPart().ElementsBegin()->GetProperties()[BODY_FORCE];
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				double& A = (i)->FastGetSolutionStepValue(NODAL_MASS);

				//the area is zero on lonely nodes, in this case set it to 1.00
				if(	A <= 1e-12)
				{
					A = 1.0;
				}
			}


			KRATOS_CATCH("");
		}

		//******************************************************************************************************
		//******************************************************************************************************
		void ApplyFractionalVelocityFixity()
		{
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				if(i->IsFixed(VELOCITY_X))
					(i)->Fix(FRACT_VEL_X);
				else
					(i)->Free(FRACT_VEL_X);

				if(i->IsFixed(VELOCITY_Y))
					(i)->Fix(FRACT_VEL_Y);
				else
					(i)->Free(FRACT_VEL_Y);

				if(i->IsFixed(VELOCITY_Z))
					(i)->Fix(FRACT_VEL_Z);
				else
					(i)->Free(FRACT_VEL_Z);
			}
		}


		//******************************************************************************************************
		//******************************************************************************************************
		//explicit correction for velocities and eventually stabilization terms
		bool ConvergenceCheck(const array_1d<double,3>& normDx, double toll)
		{
			KRATOS_TRY;
			double norm_vx = 0.00;
			double norm_vy = 0.00;
			double norm_vz = 0.00;


			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				const array_1d<double,3>& v = (i)->FastGetSolutionStepValue(FRACT_VEL);

				norm_vx += v[0]*v[0];
				norm_vy += v[1]*v[1];
				norm_vz += v[2]*v[2];
			}

			norm_vx = sqrt(norm_vx);
			norm_vy = sqrt(norm_vy);
			norm_vz = sqrt(norm_vz);

			if(norm_vx == 0.0) norm_vx  = 1.00;
			if(norm_vy == 0.0) norm_vy  = 1.00;
			if(norm_vz == 0.0) norm_vz  = 1.00;

			double ratio_x = normDx[0]/norm_vx;
			double ratio_y = normDx[1]/norm_vy;
			double ratio_z = normDx[2]/norm_vz;

			std::cout << "ratio_X = " << ratio_x << " ratio_Y = " << ratio_y << " ratio_Z = " << ratio_z << std::endl; 


			if(ratio_x < toll && ratio_y < toll && ratio_z < toll)
			{
				std::cout << "convergence achieved" << std::endl;
				return true;
			}
			
			return false;

			
			KRATOS_CATCH("");
		}

		//******************************************************************************************************
		//******************************************************************************************************
		void FractionalVelocityIteration(array_1d<double,3>& normDx)
		{	
			KRATOS_TRY	

			ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
			
			rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
			normDx[0] = mfracvel_x_strategy->Solve();
			KRATOS_WATCH(normDx[0]);

			rCurrentProcessInfo[FRACTIONAL_STEP] = 2;
			normDx[1] = mfracvel_y_strategy->Solve();
			KRATOS_WATCH(normDx[1]);

			if(mdomain_size == 3)
			{
				rCurrentProcessInfo[FRACTIONAL_STEP] = 3;
				normDx[2] = mfracvel_z_strategy->Solve();
			}
			KRATOS_CATCH("");
		}

		//******************************************************************************************************
		//******************************************************************************************************
		void VelocityPrediction(array_1d<double,3>& normDx)
		{	
			KRATOS_TRY	
			
			KRATOS_CATCH("");
		}



		//******************************************************************************************************
		//******************************************************************************************************
		virtual void SetEchoLevel(int Level) 
		{
			mfracvel_x_strategy->SetEchoLevel(Level);
			mfracvel_y_strategy->SetEchoLevel(Level);
			if(mdomain_size == 3)
				mfracvel_z_strategy->SetEchoLevel(Level);
			//
			mpressurestep->SetEchoLevel(Level);
		}

		//******************************************************************************************************
		//******************************************************************************************************
		virtual void Clear() 
		{
			KRATOS_WATCH("ResidualBasedNDFluidStrategy Clear Function called");
			mfracvel_x_strategy->Clear();
			mfracvel_y_strategy->Clear();
			if(mdomain_size == 3)
				mfracvel_z_strategy->Clear();
			mpressurestep->Clear();
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
		typename BaseType::Pointer mfracvel_x_strategy;
		typename BaseType::Pointer mfracvel_y_strategy;
		typename BaseType::Pointer mfracvel_z_strategy;
		typename BaseType::Pointer mpressurestep;

		double mvelocity_toll;
		double mpressure_toll;
		int mMaxVelIterations;
		int mMaxPressIterations;
		unsigned int mtime_order;
		unsigned int mprediction_order;
		bool mpredictor_corrector;
		bool mReformDofAtEachIteration;

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
		unsigned int m_step;
		int mdomain_size;
		double mOldDt;
		bool proj_is_initialized;


		/*@} */
		/**@name Private Operators*/
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
		ResidualBasedNDFluidStrategy(const ResidualBasedNDFluidStrategy& Other);


		/*@} */   

	}; /* Class ResidualBasedNDFluidStrategy */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_ND_FLUID_STRATEGY2  defined */

