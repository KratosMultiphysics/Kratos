/* *********************************************************   
*          
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2009-01-13 15:39:56 $
*   Revision:            $Revision: 1.14 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_FRACTIONALSTEP_STRATEGY )
#define  KRATOS_RESIDUALBASED_FRACTIONALSTEP_STRATEGY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
//#include "incompressible_fluid_application.h"
#include "custom_strategies/strategies/solver_configuration.h"
#include "utilities/geometry_utilities.h"

#ifdef _OPENMP
    #include "omp.h"
#endif

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
	class FractionalStepStrategy
		: public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
	{
	public:
		/**@name Type Definitions */       
		/*@{ */

		/** Counted pointer of ClassName */
		typedef boost::shared_ptr< FractionalStepStrategy<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;

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
		FractionalStepStrategy(
			ModelPart& model_part,
                        SolverConfiguration<TSparseSpace,TDenseSpace,TLinearSolver>& solver_config,
                        bool ReformDofAtEachIteration = true,
			double velocity_toll = 0.01,
			double pressure_toll = 0.01,
			int MaxVelocityIterations = 3,
			int MaxPressureIterations = 1,
			unsigned int time_order  = 2,
			unsigned int domain_size = 2,
			bool predictor_corrector = false
			)
			: SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part,false)
		{
			KRATOS_TRY

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
                                
			//veryfying that the model part has all the variables needed
			if (model_part.NodesBegin()->SolutionStepsDataHas(FRACT_VEL)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----FRACT_VEL---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(VELOCITY)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----VELOCITY---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(MESH_VELOCITY)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----MESH_VELOCITY---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(PRESSURE)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----PRESSURE---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(PRESSURE_OLD_IT)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----PRESSURE_OLD_IT---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(PRESS_PROJ)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----PRESS_PROJ---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(CONV_PROJ)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----CONV_PROJ---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(NODAL_MASS)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----NODAL_MASS---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(BODY_FORCE)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----BODY_FORCE---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(DENSITY)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----DENSITY---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(VISCOSITY)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----VISCOSITY---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(IS_STRUCTURE)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----IS_STRUCTURE---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(EXTERNAL_PRESSURE)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----EXTERNAL_PRESSURE---- variable!!!!!! ERROR","");
			if (model_part.NodesBegin()->SolutionStepsDataHas(IS_INTERFACE)==false )
				KRATOS_ERROR(std::logic_error,"Add  ----IS_INTERFACE---- variable!!!!!! ERROR","");


                        if(this->mdomain_size != solver_config.GetDomainSize())
                            KRATOS_ERROR(std::logic_error,"domain size not coinciding","")

                        this->mpfracvel_x_strategy = solver_config.pGetStrategy(std::string("vel_x_strategy"));
                        this->mpfracvel_y_strategy = solver_config.pGetStrategy(std::string("vel_y_strategy"));
                        this->mpfracvel_z_strategy = solver_config.pGetStrategy(std::string("vel_z_strategy"));
                        this->mppressurestep = solver_config.pGetStrategy(std::string("pressure_strategy"));



			
			
			
			
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
		virtual ~FractionalStepStrategy() {}

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

                        int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
                        while(	Dp_norm >= this->mpressure_toll && iteration++ < MaxPressureIterations  )
			{
				double p_norm = SavePressureIteration();				

				Dp_norm = FracStepSolution();

				if(fabs(p_norm) > 1e-10)
					Dp_norm /= p_norm;
				else
					Dp_norm = 1.0;

				if(rank == 0)  std::cout << "it = " << iteration << " Pressure Variation Norm = " << Dp_norm << std::endl;

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


			double local_p_norm = 0.0;
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				//setting the old value of the pressure to the current one
				const double& p = (i)->FastGetSolutionStepValue(PRESSURE);
				(i)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = p;

				local_p_norm+=p*p;
			}

                        double p_norm = BaseType::GetModelPart().GetCommunicator().SumAll(local_p_norm);
                        
                        //TODO: prepare for parallelization
			p_norm = sqrt(p_norm);
				

			return p_norm;
			KRATOS_CATCH("")
		}

		//******************************************************************************************************
		//******************************************************************************************************
		//explicit correction for velocities and eventually stabilization terms
		double FracStepSolution()
		{
			KRATOS_TRY

                        int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
				
			//setting the fractional velocity to the value of the velocity
			AssignInitialStepValues();

			//solve first step for fractional step velocities
			boost::timer step1time;
			this->SolveStep1(this->mvelocity_toll, this->mMaxVelIterations);
			if(rank == 0) std::cout << "step1 time " << step1time.elapsed() << std::endl;

			//solve for pressures (and recalculate the nodal area)
			boost::timer step2time;
			double Dp_norm = this->SolveStep2();
			if(rank == 0) std::cout << "pressure calculation time " << step2time.elapsed() << std::endl;


			this->ActOnLonelyNodes();

			//calculate projection terms
			boost::timer projection_time;
			this->SolveStep3();
			if(rank == 0) std::cout << "projection calculation time " << projection_time.elapsed() << std::endl;

//			if(mdomain_size == 2)
//  				this->SolveStep2_Mp();


			//correct velocities
			boost::timer vel_time;
			this->SolveStep4();
			if(rank == 0) std::cout << "velocity correction time " << vel_time.elapsed() << std::endl;

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
//			for(ModelPart::NodeIterator i = BaseType::GetModelPart().GetCommunicator().LocalMesh().NodesBegin() ;
//				i != BaseType::GetModelPart().GetCommunicator().LocalMesh().NodesEnd() ; ++i)
			{
				array_1d<double,3>& fract_v = (i)->FastGetSolutionStepValue(FRACT_VEL);
				fract_v *= (i)->FastGetSolutionStepValue(NODAL_MASS) * BDFcoeffs[0];
			}

//                        //set to zero fract_v on ghost nodes --> does nothing on serial version
                        for(ModelPart::NodeIterator i = BaseType::GetModelPart().GetCommunicator().GhostMesh().NodesBegin() ;
				i != BaseType::GetModelPart().GetCommunicator().GhostMesh().NodesEnd() ; ++i)
			{
				array_1d<double,3>& fract_v = (i)->FastGetSolutionStepValue(FRACT_VEL);
                                noalias(fract_v) = ZeroVector(3);
			}

			//add the elemental contributions for the calculation of the velocity
			//and the determination of the nodal area
			rCurrentProcessInfo[FRACTIONAL_STEP] = 6;
			for(ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin() ; 
				i != BaseType::GetModelPart().ElementsEnd() ; ++i)
			{
				(i)->InitializeSolutionStep(rCurrentProcessInfo);
			}

                        BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(FRACT_VEL);


			//solve nodally for the velocity
			if(this->mdomain_size == 2)
			{
				for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
					i != BaseType::GetModelPart().NodesEnd() ; ++i)
				{
					array_1d<double,3>& v = (i)->FastGetSolutionStepValue(VELOCITY);
					const array_1d<double,3>& fract_v = (i)->FastGetSolutionStepValue(FRACT_VEL);
					double A = (i)->FastGetSolutionStepValue(NODAL_MASS);
					
					double temp = (1.0/BDFcoeffs[0]) / A;
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
					
					double temp = (1.0/BDFcoeffs[0]) / A;
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
                        int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

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
				if(rank == 0)  std::cout << "ATTENTION: convergence NOT achieved" << std::endl;

			//clear if needed
			if(mReformDofAtEachIteration == true && mpredictor_corrector == false )
			{
				this->mpfracvel_x_strategy->Clear();
				this->mpfracvel_y_strategy->Clear();
				if(this->mdomain_size == 3)
					this->mpfracvel_z_strategy->Clear();
			}

			KRATOS_CATCH("");
		}

		//******************************************************************************************************
		//******************************************************************************************************
		double SolveStep2() 
		{
			KRATOS_TRY;
			 BaseType::GetModelPart().GetProcessInfo()[FRACTIONAL_STEP] = 4;
			return mppressurestep->Solve();
			KRATOS_CATCH("");
		}

		//******************************************************************************************************
		//******************************************************************************************************
		//calculation of projection 
		void SolveStep3()
		{
		  //			KRATOS_TRY;

#ifdef _OPENMP
                        int number_of_threads = omp_get_max_threads();
#else
                        int number_of_threads = 1;
#endif

                        vector<unsigned int> partition;
			CreatePartition(number_of_threads, BaseType::GetModelPart().Nodes().size(), partition);

                        #pragma omp parallel for schedule(static,1)
                        for(int k=0; k<number_of_threads; k++)
                        {
                            ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
                            ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k+1];

//                            ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
                            array_1d<double,3> zero = ZeroVector(3);

                            //first of all set to zero the nodal variables to be updated nodally
                            for(ModelPart::NodeIterator i = it_begin ; i != it_end ; ++i)
                            {
                                    array_1d<double,3>& press_proj = (i)->FastGetSolutionStepValue(PRESS_PROJ);
                                    array_1d<double,3>& conv_proj = (i)->FastGetSolutionStepValue(CONV_PROJ);
                                    noalias(press_proj) = zero;
                                    noalias(conv_proj) = zero;
                            }
                        }

			//add the elemental contributions for the calculation of the velocity
			//and the determination of the nodal area
                        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
			rCurrentProcessInfo[FRACTIONAL_STEP] = 5;

                        vector<unsigned int> elem_partition;
			CreatePartition(number_of_threads, BaseType::GetModelPart().Elements().size(), elem_partition);

                        #pragma omp parallel for schedule(static,1)
                        for(int k=0; k<number_of_threads; k++)
                        {
                            ModelPart::ElementIterator it_begin = BaseType::GetModelPart().ElementsBegin() + elem_partition[k];
                            ModelPart::ElementIterator it_end = BaseType::GetModelPart().ElementsBegin() + elem_partition[k+1];
                            for(ModelPart::ElementIterator i = it_begin ; i!=it_end ; ++i)
                            {
                                    (i)->InitializeSolutionStep(BaseType::GetModelPart().GetProcessInfo());
                            }
                        }

                        BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(NODAL_MASS);
                        BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(PRESS_PROJ);
                        BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(CONV_PROJ);

			//solve nodally for the velocity
                        #pragma omp parallel for schedule(static,1)
                        for(int k=0; k<number_of_threads; k++)
                        {
                            ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
                            ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k+1];

                            for(ModelPart::NodeIterator i =it_begin ; i != it_end ; ++i)
                            {
                                    array_1d<double,3>& press_proj = (i)->FastGetSolutionStepValue(PRESS_PROJ);
                                    array_1d<double,3>& conv_proj = (i)->FastGetSolutionStepValue(CONV_PROJ);
                                    double A = (i)->FastGetSolutionStepValue(NODAL_MASS);

                                    double temp = 1.00 / A;
                                    press_proj *= temp;
                                    conv_proj *= temp;

                            }
                        }

		
			//	KRATOS_CATCH("");
		}

		//******************************************************************************************************
		//******************************************************************************************************
		//correct pressure taking in account viscosity (should accelerate the convergence of the predictor corrector) 
		void SolveStep2_Mp()
		{
			KRATOS_TRY;

			ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
			array_1d<double,3> zero = ZeroVector(3);
			Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				//setting the old value of the pressure to the current one
				double& p = (i)->FastGetSolutionStepValue(PRESSURE);
				(i)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = p;
				p = 0.0;
				(i)->FastGetSolutionStepValue(NODAL_MASS) = 0.0;


			}

			//const Vector& BDFcoeffs = BaseType::GetModelPart().GetProcessInfo()[BDF_COEFFICIENTS];

			//calculate divergence of the fractional velocity element by element
			if(this->mdomain_size == 2)
			{
				array_1d<double, 3> N;
				array_1d<double, 2> temp, proj_aux;
				array_1d<double, 3> pressures;
				array_1d<double, 3> elemental_stabilization;
				boost::numeric::ublas::bounded_matrix <double, 3,2> DN_DX;
				array_1d<double,3> vg;
				for(ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin() ; 
					i != BaseType::GetModelPart().ElementsEnd() ; ++i)
				{
					Geometry< Node<3> >& geom = i->GetGeometry();
					double volume;

					//calculate derivatives 
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);

					//getting fractional velocities on nodes
					const array_1d<double,3>& fv0 = geom[0].FastGetSolutionStepValue(FRACT_VEL);
					const array_1d<double,3>& fv1 = geom[1].FastGetSolutionStepValue(FRACT_VEL);
					const array_1d<double,3>& fv2 = geom[2].FastGetSolutionStepValue(FRACT_VEL);

					pressures[0] = geom[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
					pressures[1] = geom[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
					pressures[2] = geom[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);

					//calculate N
					double Gaux;
					Gaux =  DN_DX(0,0)*fv0[0] + DN_DX(0,1)*fv0[1];
					Gaux += DN_DX(1,0)*fv1[0] + DN_DX(1,1)*fv1[1];
					Gaux += DN_DX(2,0)*fv2[0] + DN_DX(2,1)*fv2[1];

					double density = geom[0].FastGetSolutionStepValue(DENSITY);
					double nu = geom[0].FastGetSolutionStepValue(VISCOSITY);
					for(int ii = 1; ii<3; ii++)
					{
						density += geom[ii].FastGetSolutionStepValue(DENSITY);
						nu += geom[ii].FastGetSolutionStepValue(VISCOSITY);
					}
					density *= 0.333333333333333333333;
					nu *= 0.333333333333333333333;

					noalias(vg) = ZeroVector(3);; //velocity on the gauss points
					for(int kk = 0; kk<3; kk++)
					{
						//adding the elemental contribution to the nodal volume
						geom[kk].FastGetSolutionStepValue(NODAL_MASS) += density*volume*N[kk];
	
						//calculating the velocity on the gauss point
						noalias(vg) +=	N[kk]* ( geom[kk].FastGetSolutionStepValue(FRACT_VEL) - 
							geom[kk].FastGetSolutionStepValue(MESH_VELOCITY) );
					}

					//calculating the stabilization
					double dt_contrib_to_tau = 0.0;
					if( muse_dt_in_stabilization == true)
						dt_contrib_to_tau = 1.0/BDFcoeffs[0];
					double	h = sqrt(2.0*volume);
					double c1 = 4.00;
					double c2 = 2.00;
					double norm_u = norm_2(vg);
					double tau = 1.00 / (dt_contrib_to_tau + c1*nu/(h*h) + c2*norm_u/h );

					//calculating stabilization laplacian LHS
					noalias(temp) = prod(trans(DN_DX),pressures);
					noalias(elemental_stabilization) =  -tau/density * prod(DN_DX,temp);
	
					const array_1d<double,3>& proj_temp = geom[0].FastGetSolutionStepValue(PRESS_PROJ);
					for(int iii = 0; iii<2; iii++)
						proj_aux[iii] = N[0]*proj_temp[iii];
					for(int kk = 1; kk<3; kk++)
					{
						const array_1d<double,3>& proj_temp = geom[kk].FastGetSolutionStepValue(PRESS_PROJ);
						for(int iii = 0; iii<2; iii++)
							proj_aux[iii] += N[kk]*proj_temp[iii];
					}
					proj_aux *= tau;
					noalias(elemental_stabilization) += prod(DN_DX , proj_aux);

					geom[0].FastGetSolutionStepValue(PRESSURE) += (elemental_stabilization[0] - Gaux*N[0])*(volume*density);
					geom[1].FastGetSolutionStepValue(PRESSURE) += (elemental_stabilization[1] - Gaux*N[1])*(volume*density);
					geom[2].FastGetSolutionStepValue(PRESSURE) += (elemental_stabilization[2] - Gaux*N[2])*(volume*density);

					double nodal_vol = volume/3.0;
					geom[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_vol * density;
					geom[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_vol * density;
					geom[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_vol * density;
				}
			}
			else
			{
				KRATOS_ERROR(std::logic_error,"not yet implemented the 3D","");
			
			}

                        BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(NODAL_MASS);
                        BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(PRESSURE);
		
			//correct pressure
			//double temp;
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				double& p = (i)->FastGetSolutionStepValue(PRESSURE);
				const double p_old_it = (i)->FastGetSolutionStepValue(PRESSURE_OLD_IT);
				const double rho = (i)->FastGetSolutionStepValue(DENSITY);
				const double nu = (i)->FastGetSolutionStepValue(VISCOSITY);
				const double mass = (i)->FastGetSolutionStepValue(NODAL_MASS);
				
				if(mass >1e-12)
				{
   					p *= nu*rho/mass;
   					p += p_old_it; 
				}
				else
				{
					(i)->FastGetSolutionStepValue(NODAL_MASS) = 1.0;
					p = p_old_it;
				}
				
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
#ifdef _OPENMP
                        int number_of_threads = omp_get_max_threads();
#else
                        int number_of_threads = 1;
#endif

                        vector<unsigned int> partition;
			CreatePartition(number_of_threads, BaseType::GetModelPart().Nodes().size(), partition);

//                        #pragma omp parallel for schedule(static,1)
                        for(int k=0; k<number_of_threads; k++)
                        {
                            ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
                            ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k+1];
                                for(ModelPart::NodeIterator i = it_begin ; i != it_end ; ++i)
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

                        BaseType::GetModelPart().GetCommunicator().SumAll(norm_vx);
                        BaseType::GetModelPart().GetCommunicator().SumAll(norm_vy);
                        BaseType::GetModelPart().GetCommunicator().SumAll(norm_vz);

			norm_vx = sqrt(norm_vx);
			norm_vy = sqrt(norm_vy);
			norm_vz = sqrt(norm_vz);

			if(norm_vx == 0.0) norm_vx  = 1.00;
			if(norm_vy == 0.0) norm_vy  = 1.00;
			if(norm_vz == 0.0) norm_vz  = 1.00;

			double ratio_x = normDx[0]/norm_vx;
			double ratio_y = normDx[1]/norm_vy;
			double ratio_z = normDx[2]/norm_vz;

                        int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
                        if(rank == 0)  std::cout << "ratio_x = " << ratio_x << " ratio_Y = " << ratio_y << " ratio_Z = " << ratio_z << std::endl;


			if(ratio_x < toll && ratio_y < toll && ratio_z < toll)
			{
				if(rank == 0)  std::cout << "convergence achieved" << std::endl;
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
			normDx[0] = mpfracvel_x_strategy->Solve();

			rCurrentProcessInfo[FRACTIONAL_STEP] = 2;
			normDx[1] = mpfracvel_y_strategy->Solve();

			if(mdomain_size == 3)
			{
				rCurrentProcessInfo[FRACTIONAL_STEP] = 3;
				normDx[2] = mpfracvel_z_strategy->Solve();
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
			mpfracvel_x_strategy->SetEchoLevel(Level);
			mpfracvel_y_strategy->SetEchoLevel(Level);
			if(mdomain_size == 3)
				mpfracvel_z_strategy->SetEchoLevel(Level);
			//
			mppressurestep->SetEchoLevel(Level);
		}

		//******************************************************************************************************
		//******************************************************************************************************
		virtual void Clear() 
		{
                        int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
                        if(rank == 0)  KRATOS_WATCH("FractionalStepStrategy Clear Function called");
			mpfracvel_x_strategy->Clear();
			mpfracvel_y_strategy->Clear();
			if(mdomain_size == 3)
				mpfracvel_z_strategy->Clear();
			mppressurestep->Clear();
		}

		virtual double GetStageResidualNorm(unsigned int step)
		{
			if(step == 1)
				return mpfracvel_x_strategy->GetResidualNorm();
			if(step == 2)
				return mpfracvel_y_strategy->GetResidualNorm();
			if(step == 3)
				return mpfracvel_z_strategy->GetResidualNorm();
			if(step == 4)
				return mppressurestep->GetResidualNorm();
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
		typename BaseType::Pointer mpfracvel_x_strategy;
		typename BaseType::Pointer mpfracvel_y_strategy;
		typename BaseType::Pointer mpfracvel_z_strategy;
		typename BaseType::Pointer mppressurestep;

		double mvelocity_toll;
		double mpressure_toll;
		int mMaxVelIterations;
		int mMaxPressIterations;
		unsigned int mtime_order;
		unsigned int mprediction_order;
		bool mpredictor_corrector;
		bool mReformDofAtEachIteration;

		bool muse_dt_in_stabilization;


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
		unsigned int mdomain_size;
		double mOldDt;
		bool proj_is_initialized;

                //******************************************************************************************
		//******************************************************************************************
		inline void CreatePartition(unsigned int number_of_threads,const int number_of_rows, vector<unsigned int>& partitions)
		{
			partitions.resize(number_of_threads+1);
			int partition_size = number_of_rows / number_of_threads;
			partitions[0] = 0;
			partitions[number_of_threads] = number_of_rows;
			for(unsigned int i = 1; i<number_of_threads; i++)
			   partitions[i] = partitions[i-1] + partition_size ;
		}


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
		FractionalStepStrategy(const FractionalStepStrategy& Other);


		/*@} */   

	}; /* Class FractionalStepStrategy */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_FRACTIONALSTEP_STRATEGY  defined */

