/* *********************************************************   
*          
*   Last Modified by:    $Author: jmarti $
*   Date:                $Date: 2009-01-23 14:34:00 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_FLUID_STRATEGY2 ) /////CAMBIARRRR
#define  KRATOS_RESIDUALBASED_FLUID_STRATEGY2


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "Qcompressible_fluid_application.h"

#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
//#include "custom_strategies/builder_and_solvers/residualbased_elimination_discretelaplacian_builder_and_solver.h"

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
	class ResidualBasedFluidStrategy 
		: public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
	{
	public:
		/**@name Type Definitions */       
		/*@{ */

		/** Counted pointer of ClassName */
		typedef boost::shared_ptr< ResidualBasedFluidStrategy<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;

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
		ResidualBasedFluidStrategy(
			ModelPart& model_part, 
			typename TLinearSolver::Pointer pNewVelocityLinearSolver,
			typename TLinearSolver::Pointer pNewPressureLinearSolver,
			bool CalculateReactions = false,
			bool ReformDofAtEachIteration = true,
			bool CalculateNormDxFlag = false,
			double velocity_toll = 0.01,
			double pressure_toll = 0.01,
			int MaxVelocityIterations = 10,
			int MaxPressureIterations = 4,
			unsigned int time_order  = 2,
			unsigned int domain_size = 2,
			unsigned int laplacian_form = 1, //1 = laplacian, 2 = discrete laplacian
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
			ReformDofAtEachIteration = false;

				//initializing fractional velocity solution step
				typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
			typename SchemeType::Pointer pscheme = typename SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,  TDenseSpace >() );

			typedef typename Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > VarComponent;
			
			typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;
			
			BuilderSolverTypePointer vel_x_build =BuilderSolverTypePointer(	new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(pNewVelocityLinearSolver) );
			
			
			this->mfracvel_strategy = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >				(model_part,pscheme,pNewVelocityLinearSolver,vel_x_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
			this->mfracvel_strategy->SetEchoLevel(2);
			
			
			if( laplacian_form == 1) //laplacian form
			  {
				std::cout << "standard laplacian form" << std::endl;
				
					    
		
	  }
			else if( laplacian_form == 2) //discrete laplacian form
			  {
			    std::cout << "discrete laplacian form" << std::endl;
			    
			    if(domain_size == 2)
			      {
				
			      }
			    else if (domain_size == 3)
				  {
				    
				  }
			    
			   
			    
			  }
			else if( laplacian_form == 3) //discrete laplacian form - stabilized only with dt
			  			
			    
			  {  
			    
			  }
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
			    i != BaseType::GetModelPart().NodesEnd() ; ++i)
			  {
			    if(i->IsFixed(VELOCITY_X)){
			      (i)->Fix(FRACT_VEL_X);
			    }
			    if(i->IsFixed(VELOCITY_Y)){
			      (i)->Fix(FRACT_VEL_Y);
			    }
			    if(i->IsFixed(VELOCITY_Z))
			      (i)->Fix(FRACT_VEL_Z);
			  }
			
			this->m_step = 1;
			
			this->mOldDt  =0.00;
			
			KRATOS_CATCH("")
			  }
		
		
		
		/** Destructor.
		 */
		virtual ~ResidualBasedFluidStrategy() {}
		
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
		  
		  double Dt = 0.1;
		  
		  		  
		  AssignInitialValues();
		  //ApplyFractionalVelocityFixity();

		  //PredictVelocity(1,1);
		  Dp_norm = IterativeSolve();
		  

		  if(this->mReformDofAtEachIteration == true )
		    this->Clear();
		  
		  this->m_step += 1;
		  this->mOldDt =  BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME];
		  
		  //Nodemove();
		  
		  
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
			
			
			    
			    
			    			
			    Dp_norm = FracStepSolution();
			    
			    			
			//if(this->mReformDofAtEachIteration == true )
			//  this->Clear();
			
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
		      
		      (i)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = p;
		      
		      p_norm+=p*p;
		    }	
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
				
			//setting the fractional velocity to the value of the velocity
			  AssignInitialStepValues();

			//solve first step for fractional step velocities
			boost::timer step1time;
			this->SolveStep1(this->mvelocity_toll, this->mMaxVelIterations);
			std::cout << "step1 time " << step1time.elapsed() << std::endl;
			

			//solve for pressures (and recalculate the nodal area)
			boost::timer step2time;
			double Dp_norm =0.001;//=this->SolveStep2();
			std::cout << "pressure calculation time " << step2time.elapsed() << std::endl;
			
			/*boost::timer vel_time;
			this->SolveStep4();
			std::cout << "velocity correction time " << vel_time.elapsed() << std::endl;
	*/
			//this->SolveStep3();
			//double Dp_norm=0.001;

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
			double Dt = rCurrentProcessInfo[DELTA_TIME];
		

			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
			i != BaseType::GetModelPart().NodesEnd() ; ++i)
		      	{
			(i)->FastGetSolutionStepValue(NODAL_MAUX) = 0;
			(i)->FastGetSolutionStepValue(NODAL_PAUX) = 0;

			double pp=(i)->FastGetSolutionStepValue(NODAL_MAUX);
			}

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
				double& p = (i)->FastGetSolutionStepValue(PRESSURE);
				const int f= (i)->FastGetSolutionStepValue(FLAG_VARIABLE); 
				const int fi= (i)->FastGetSolutionStepValue(IS_INTERFACE); 
				if(A <= 1e-8 ) //hay q mirar esto
				{
				    p=0; 		
				}
				else{	
					double B = (i)->FastGetSolutionStepValue(NODAL_PRESS);//
	      				if(A <= 1e-8){p=0;}
					else{ p=B/A;}

				}				

				
				/*if(!i->IsFixed(VELOCITY_X))
					{
					  v[0] = fract_v[0];
					}
				if(!i->IsFixed(VELOCITY_Y))
				  {
				    v[1] = fract_v[1];
				  }*/				
			      }
			}
			
			KRATOS_CATCH("");
		}
		
		//******************************************************************************************************

		void SolveStep3()
		{
			KRATOS_TRY;

			ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
			array_1d<double,3> zero = ZeroVector(3);

			//first of all set to zero the nodal variables to be updated nodally
			

			//add the elemental contributions for the calculation of the velocity
			//and the determination of the nodal area
			
			double T= rCurrentProcessInfo[TIME];	
			  
		
		  	if(T>0.085){



				rCurrentProcessInfo[FRACTIONAL_STEP] = 5;
				for(ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin() ;i != BaseType::GetModelPart().ElementsEnd() ; ++i)
				{
					(i)->InitializeSolutionStep(rCurrentProcessInfo);
				}
				if(this->mdomain_size == 2)
				{
					for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; i != BaseType::GetModelPart().NodesEnd() ; ++i)
					{
	
						double A1 = (i)->FastGetSolutionStepValue(NODAL_MASSAUX);	
				/*****/
						double& d = (i)->FastGetSolutionStepValue(DENSITY);
				/*****/
						const int f= (i)->FastGetSolutionStepValue(FLAG_VARIABLE); 
						const int fi= (i)->FastGetSolutionStepValue(IS_INTERFACE); 
				
				//double temp = (1.0 / A);
			//d=1.3;

						if( A1 <= 1e-12 ) //hay q mirar esto
						{
						d = 1.3;		
						}
						else	
						{
							if(f==1)
							{	
   								d = 1.3;//NO TIENE QUE HACER NADA
							}
							if(f==2)
							{
								double B1 = (i)->FastGetSolutionStepValue(NODAL_DENSITYAUX);//
								d = B1/A1;
							}	
						}			
					}
				}


			}

		KRATOS_CATCH("");
	}

		//******************************************************************************************************
		//******************************************************************************************************



		//******************************************************************************************************
		//explicit correction for velocities and eventually stabilization terms
		void Nodemove()
		{
		  KRATOS_TRY;
		  
		  ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
		  array_1d<double,3> zero = ZeroVector(3);
		  Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
		  double Dt = rCurrentProcessInfo[DELTA_TIME];
		  
	 	double T= rCurrentProcessInfo[TIME];	
			  
		
		  if(T>0.085){



		rCurrentProcessInfo[FRACTIONAL_STEP] = 5;
		for(ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin() ;i != BaseType::GetModelPart().ElementsEnd() ; ++i)
		{
		(i)->InitializeSolutionStep(rCurrentProcessInfo);
		}
		if(this->mdomain_size == 2)
		{
		for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{

			double A1 = (i)->FastGetSolutionStepValue(NODAL_MASSAUX);	
				/*****/
			double& d = (i)->FastGetSolutionStepValue(DENSITY);
				/*****/
			const int f= (i)->FastGetSolutionStepValue(FLAG_VARIABLE); 
			const int fi= (i)->FastGetSolutionStepValue(IS_INTERFACE); 
				
				//double temp = (1.0 / A);
			d=1.3;

			if( A1 <= 1e-12 ) //hay q mirar esto
			{
			d = 1.3;		
			}
			else	
			{
			if(f==1)
			{	
   				d = 1.3;//NO TIENE QUE HACER NADA
			}
			if(f==2)
			{
				double B1 = (i)->FastGetSolutionStepValue(NODAL_DENSITYAUX);//
				d = B1/A1;
			}	
			}			
			}
			}

}
		  /*if(this->mdomain_size == 2){
		    
		    
		    for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; i != BaseType::GetModelPart().NodesEnd() ; ++i)
		      {
			const array_1d<double,3>& v = (i)->FastGetSolutionStepValue(VELOCITY);
			array_1d<double,3>& disp = (i)->FastGetSolutionStepValue(DISPLACEMENT);			  
			const array_1d<double,3>& disp_old = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
			const array_1d<double,3>& v_old = (i)->FastGetSolutionStepValue(VELOCITY,1);
			
			disp = disp_old+ v*Dt ; 
			
			i->X() = i->X0() + disp[0];
			i->Y() = i->Y0() + disp[1];
		      }
		  }
		  
		  else{
		    for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; i != BaseType::GetModelPart().NodesEnd() ; ++i)
		      {
			const array_1d<double,3>& v = (i)->FastGetSolutionStepValue(VELOCITY);
			array_1d<double,3>& disp = (i)->FastGetSolutionStepValue(DISPLACEMENT);			  
			const array_1d<double,3>& disp_old = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
			const array_1d<double,3>& v_old = (i)->FastGetSolutionStepValue(VELOCITY,1);
			
			
			noalias(disp) = disp_old+ v*Dt ; 
			
			
			i->X() = i->X0() + disp[0];
			i->Y() = i->Y0() + disp[1];
			i->Z() = i->Z0() + disp[2];
		      }
		    
		    
		    
		  }		*/  
		  
		  
		  KRATOS_CATCH("");
		}
		//******************************************************************************************************
		//******************************************************************************************************
		//sets the BDF coefficients to the correct value
		void InitializeFractionalStep(const int step,const int time_order)
		{
		   ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
		  double T= rCurrentProcessInfo[TIME];	
		  KRATOS_TRY;
		  
		  //calculate the BDF coefficients
		 
		  double Dt = rCurrentProcessInfo[DELTA_TIME];
		
		  if(T>0.085){
			double pp=0;
		}
			
		  if(this->mOldDt == 0.00) //needed for the first step
		    this->mOldDt = Dt;
		  
		  
		  rCurrentProcessInfo[BDF_COEFFICIENTS].resize(2,false);
		  Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
		  BDFcoeffs[0] =	1.0 / Dt;	//coefficient for step n+1
		  BDFcoeffs[1] =	-1.0 / Dt;//coefficient for step n
		  /*}*/
		  
		  
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
		void AssignInitialValues()
		{
		  KRATOS_TRY
		    for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
			i != BaseType::GetModelPart().NodesEnd() ; ++i)
		      {
			//predicting the values for the fluid velocity
			array_1d<double,3>& v = (i)->FastGetSolutionStepValue(VELOCITY);
			array_1d<double,3>& fracv = (i)->FastGetSolutionStepValue(FRACT_VEL);
			array_1d<double,3> zero = ZeroVector(3);			
			//noalias(fracv) = v; 
			
			
			(i)->FastGetSolutionStepValue(PRESSURE)=0;////nuevo
			(i)->FastGetSolutionStepValue(PRESSUREAUX)=0;////nuevo
			
			//resetting the nodal area
			double area = 0.00;
			(i)->FastGetSolutionStepValue(NODAL_MASS) = area;
			(i)->FastGetSolutionStepValue(NODAL_MASSAUX) = area;
			(i)->FastGetSolutionStepValue(NODAL_DENSITYAUX) = area;
			int fi= (i)->FastGetSolutionStepValue(IS_STRUCTURE);//(IS_FREE_SURFACE)
			/*if(fi){

 			array_1d<double,3>& fract_v = (i)->FastGetSolutionStepValue(FRACT_VEL);
			array_1d<double,3>& vel = (i)->FastGetSolutionStepValue(VELOCITY);//ZeroVector(3);
			fract_v[0]=0.05;
			vel[0]=0.05;
			(i)->Fix(VELOCITY_X);
			(i)->Fix(VELOCITY_Y);
			(i)->Fix(VELOCITY_Z);
			(i)->Fix(FRACT_VEL_X);
			(i)->Fix(FRACT_VEL_Y);
			(i)->Fix(FRACT_VEL_Z);


			}*/
		      }
		  KRATOS_CATCH("");
		}
		
		//******************************************************************************************************
		//******************************************************************************************************
		void AssignInitialStepValues()
		{
		  KRATOS_TRY
		    for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
			i != BaseType::GetModelPart().NodesEnd() ; ++i)
		      {
			
			//(i)->FastGetSolutionStepValue(PRESSURE)=0;////nuevo	
			//(i)->FastGetSolutionStepValue(PRESSUREAUX)=0;////nuevo	
			//resetting the nodal area
			double area = 0.00;
			array_1d<double,3> desp = ZeroVector(3);
			(i)->FastGetSolutionStepValue(NODAL_MASS) = area;
			/*****/	
			(i)->FastGetSolutionStepValue(NODAL_MASSAUX) = area;
			/*****/	
			(i)->FastGetSolutionStepValue(DESP) = desp;
			
			double press = 0.00;
			(i)->FastGetSolutionStepValue(NODAL_PRESS) = press;
				/*****/		
			(i)->FastGetSolutionStepValue(NODAL_PRESSAUX) = press;
			/*****/		
			//double pp=(i)->FastGetSolutionStepValue(NODAL_MASSAUX);
			//KRATOS_WATCH(pp);////nuevo
			
		      }
		  KRATOS_CATCH("");
		}
		
		//******************************************************************************************************
		//******************************************************************************************************
		void AssignInitialDesp()
		{
		  KRATOS_TRY
		    for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
			i != BaseType::GetModelPart().NodesEnd() ; ++i)
		      {
 			ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
			double Dt = rCurrentProcessInfo[DELTA_TIME];
			const array_1d<double,3>& fract_v = (i)->FastGetSolutionStepValue(FRACT_VEL);
			array_1d<double,3>& desp = (i)->FastGetSolutionStepValue(DESP);//ZeroVector(3);
			noalias(desp)=fract_v*Dt;

			
			
		      }
		  KRATOS_CATCH("");
		}
		
		//******************************************************************************************************
		//******************************************************************************************************
		void PredictVelocity(int step,int prediction_order)
		{
		  KRATOS_TRY

		      for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ;   i != BaseType::GetModelPart().NodesEnd() ; ++i)
				{
				  array_1d<double,3>& vel = i->FastGetSolutionStepValue(VELOCITY); 
				  array_1d<double,3>& vel1 = i->FastGetSolutionStepValue(FRACT_VEL);
				  	
				  if(!i->IsFixed(VELOCITY_X)) {vel1[0] = 0.5;};// i->Fix(FRACT_VEL_X);}
				  if(!i->IsFixed(VELOCITY_Y)) {vel1[1] = vel[1];} //i->Fix(FRACT_VEL_Y);}
				  if(!i->IsFixed(VELOCITY_Z)) {vel1[2] = vel[2];}// i->Fix(FRACT_VEL_X);}
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
		  
	          FractionalVelocityIteration(normDx);
		  //clear if needed
		  if(mReformDofAtEachIteration == true && mpredictor_corrector == false )
		    {
		      this->mfracvel_strategy->Clear();
		      
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
			    if(i->IsFixed(VELOCITY_X)){
				std::cout << "VLOCIDADDDDDDDDDDDDD " ;
				std::cout << "VLOCIDADDDDDDDDDDDDD " ;
				std::cout << "VLOCIDADDDDDDDDDDDDD " ;
				std::cout << "VLOCIDADDDDDDDDDDDDD " ;
			      (i)->Fix(FRACT_VEL_X); }
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
			
			double norm_tot=0.00;

			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				const array_1d<double,3>& v = (i)->FastGetSolutionStepValue(FRACT_VEL);

				norm_tot += v[0]*v[0];
				norm_tot += v[1]*v[1];
				norm_tot += v[2]*v[2];


			}

			
			norm_tot=sqrt(norm_tot);
			if(norm_tot==0.0) norm_tot=1.00;
			double ratio=normDx[0]/norm_tot;
			
			

			/*std::cout << "ratio_x = " << ratio_x << " ratio_Y = " << ratio_y << " ratio_Z = " << ratio_z << std::endl; */
			std::cout << "toll = " << toll << " ratio = " << ratio << std::endl; 
			std::cout << "norm_tot = " << norm_tot <<  std::endl; 
			/*if(ratio_x < toll && ratio_y < toll && ratio_z < toll)
			{
				std::cout << "convergence achieved" << std::endl;
				return true;
				}*/
			if(ratio < toll)
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
			//AssignInitialStepValues();
			rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
			normDx[0] = mfracvel_strategy->Solve();/*mfracvel_x_strategy->Solve();*/
								 
			normDx[1] = normDx[0];/*mfracvel_y_strategy->Solve();*/
			normDx[2] = normDx[0];
						if(mdomain_size == 3)
			{
				
				normDx[2] = normDx[0];
				}


			//printing all of the nodes with negative Id

if(rCurrentProcessInfo[TIME] > 0.08)
{

			void* p_temp;					 /*rCurrentProcessInfo[FRACTIONAL_STEP] = 2;*/

	
			//looping on all elements and then on the geometry of all of the nodes
			for(ModelPart::ElementsContainerType::const_iterator ie = BaseType::GetModelPart().ElementsBegin();ie != BaseType::GetModelPart().ElementsEnd(); ie++)	
			{
				Geometry< Node<3> >& geom = ie->GetGeometry();

				for(int i=0; i<geom.size();i++)
				{
					Node<3>::DofsContainerType& reference_dofs = geom[i].GetDofs();

					//generating the dofs
					for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)
					{
						if( (iii->GetVariable() == FRACT_VEL_X && iii->EquationId() > 100000)  || (iii->GetVariable() == FRACT_VEL_Y && iii->EquationId() > 100000)  )
						{	
							KRATOS_WATCH(geom[i].Id());
							KRATOS_WATCH(&(geom[i]) );					
							KRATOS_WATCH( & (*(BaseType::GetModelPart().Nodes().find(geom[i].Id())) ));

							KRATOS_WATCH( &(*(geom[i].pGetDof(FRACT_VEL_X) ))     );
							KRATOS_WATCH((&(*((BaseType::GetModelPart().Nodes().find(geom[i].Id()))->pGetDof(FRACT_VEL_X)  ))) );

// 							KRATOS_WATCH( geom[i].pGetDof(FRACT_VEL_X)->GetVariable().Key()     );
// 							KRATOS_WATCH( geom[i].pGetDof(FRACT_VEL_X)->Id()     );	
// 							KRATOS_WATCH( (BaseType::GetModelPart().Nodes().find(geom[i].Id()))->pGetDof(FRACT_VEL_X)->GetVariable().Key() );
// 							KRATOS_WATCH( (BaseType::GetModelPart().Nodes().find(geom[i].Id()))->pGetDof(FRACT_VEL_X)->Id() );
						}
					}
	/*				if (geom[i].Id() == 16)
{
std::cout << "Elements16 ";
						KRATOS_WATCH(geom[i]);
						KRATOS_WATCH(&(geom[i]));
						KRATOS_WATCH(p_temp);
}*/
				}

			}
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

		  mfracvel_strategy->SetEchoLevel(Level);
		
		}

		//******************************************************************************************************
		//******************************************************************************************************
		virtual void Clear() 
		{
		  mfracvel_strategy->Clear();
		  //mpressurestep->Clear();
		  
		
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
		typename BaseType::Pointer mfracvel_strategy;
		/*typename BaseType::Pointer mfracvel_x_strategy;
		typename BaseType::Pointer mfracvel_y_strategy;
		typename BaseType::Pointer mfracvel_z_strategy;*/
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
		ResidualBasedFluidStrategy(const ResidualBasedFluidStrategy& Other);


		/*@} */   

	}; /* Class ResidualBasedFluidStrategy */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_FLUID_STRATEGY2  defined */

