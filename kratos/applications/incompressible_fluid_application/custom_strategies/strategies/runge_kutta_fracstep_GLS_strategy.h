/* *********************************************************   
*          
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2008-11-10 14:23:32 $
*   Revision:            $Revision: 1.12 $
*
* ***********************************************************/


#if !defined(KRATOS_RUNGE_KUTTA_GLS_STRATEGY)
#define  KRATOS_RUNGE_KUTTA_GLS_STRATEGY


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "incompressible_fluid_application.h"

//#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
//#include "custom_strategies/builder_and_solvers/residualbased_elimination_discretelaplacian_builder_and_solver.h"
//#include "custom_strategies/builder_and_solvers/residualbased_elimination_discretelaplacian_builder_and_solver_flexiblefsi.h"

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
	template<unsigned int TDim, class TSparseSpace,
	class TDenseSpace,
	class TLinearSolver
	>
	class RungeKuttaFracStepStrategy 
		: public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
	{
	public:
		/**@name Type Definitions */       
		/*@{ */

		/** Counted pointer of ClassName */
		typedef std::vector<unsigned int> IndicesVectorType;
		
		typedef boost::shared_ptr< RungeKuttaFracStepStrategy<TDim, TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;

		typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

		typedef typename BaseType::TDataType TDataType;

		//typedef typename BaseType::DofSetType DofSetType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

		typedef Node<3> PointType;
		
		typedef Node<3>::Pointer PointPointerType;

		typedef std::vector<PointType::Pointer>           PointVector;

		typedef PointVector::iterator PointIterator;



		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */

		/** Constructor.
		*/
		RungeKuttaFracStepStrategy(
			ModelPart& model_part, 
			typename TLinearSolver::Pointer pNewPressureLinearSolver,
			bool CalculateReactions = false,
			bool ReformDofAtEachIteration = true,
			bool CalculateNormDxFlag = true
			//double velocity_toll = 0.01,
			//double pressure_toll = 0.01,
			//int MaxVelocityIterations = 3,
			//int MaxPressureIterations = 1,
			//unsigned int time_order  = 2,
			//unsigned int prediction_order  = 2,
			//unsigned int domain_size = 2
			//unsigned int laplacian_form = 2, //1 = laplacian, 2 = discrete laplacian
			//bool predictor_corrector = false
			)
			: SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part,false)
		{
			KRATOS_TRY
//std::cout << "SONO QUI" << std::endl;
			//this->mvelocity_toll = velocity_toll;
			//this->mpressure_toll = pressure_toll;
			//this->mMaxVelIterations = MaxVelocityIterations;
			//this->mMaxPressIterations = MaxPressureIterations;
			//this->mtime_order = time_order;
			//this->mprediction_order = time_order;
			//this->mdomain_size = domain_size;

			//this->mpredictor_corrector = predictor_corrector;
			this->mReformDofAtEachIteration = ReformDofAtEachIteration;

			//this->proj_is_initialized = false;

			//the system will be cleared at the end!
			//ReformDofAtEachIteration = false;

				//initializing fractional velocity solution step
				typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
			typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
				( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,  TDenseSpace >() );

			//commented the 3 lines below

			/*
			bool CalculateReactions = false;
			bool ReformDofAtEachIteration = true;
			bool CalculateNormDxFlag = true;
			*/


			//computation of the fractional vel velocity (first step)
				//3 dimensional case
			//typedef typename Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > > VarComponent;

			//typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

			
			//std::cout << "standard laplacian form" << std::endl;
			//mmin_conv_vel_norm = 0.0;
			this->mpressurestep = typename BaseType::Pointer(
				new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >
				(model_part,pscheme,pNewPressureLinearSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
			this->mpressurestep->SetEchoLevel(2); 
			
			//identify nodes, weher slip shall be imposed .. store them in a list
			mSlipBoundaryList.clear();
			//ModelPart& model_part=GetModelPart();
			for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
			{
			//FLAG_VAR = 1 is SLIP
			if (it->FastGetSolutionStepValue(FLAG_VARIABLE)==3.0)		

			mSlipBoundaryList.push_back(*(it.base()));
			//mSlipBoundaryList.push_back(it);
			}

			for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
			{
				if(it->pGetDof(VELOCITY_X)->IsFixed() == true)
				{
					mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_X) );
					mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_X)->GetSolutionStepValue() );
				}

				if(it->pGetDof(VELOCITY_Y)->IsFixed() == true)
				{
					mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_Y) );
					mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_Y)->GetSolutionStepValue() );
				}

				if(it->pGetDof(VELOCITY_Z)->IsFixed() == true)
				{
					mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_Z) );
					mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_Z)->GetSolutionStepValue() );
				}
			}
			
			

			Element & ref_el = model_part.Elements().front();
			Geometry<Node<3> >::Pointer p_null_geom=Geometry< Node<3> >::Pointer(new Geometry< Node<3> >);

			//int id=1;
			if (TDim==2)
			{
			Fluid2DGLS_expl el(1, p_null_geom);

			if (typeid(ref_el) != typeid(el))
				KRATOS_ERROR(std::logic_error,  "Incompressible Runge Kutta Strategy requires utilization of Fluid2DGLS_expl elements " , "");
			}

			if (TDim==3)
			{
			Fluid3DGLS_expl el(1, p_null_geom);

			if (typeid(ref_el) != typeid(el))
				KRATOS_ERROR(std::logic_error,  "Incompressible Runge Kutta Strategy requires utilization of Fluid3DGLS_expl elements " , "");
			}

			KRATOS_CATCH("")
		}



		/** Destructor.
		*/
		virtual ~RungeKuttaFracStepStrategy() {}

		/** Destructor.
		*/

		//*********************************************************************************
		//**********************************************************************
		double Solve()
		{
		KRATOS_WATCH("Solve of Runge Kutta GLS Frac Step Strategy")
		//we estimate the time step for the explicit time integration schem estability
		//ComputeTimeStep(0.8);
		Timer::Start("SolveStep1");
		SolveStep1();
		Timer::Stop("SolveStep1");
		KRATOS_WATCH("First frac step completed")
		//we write now the beginning of step pressure to OLD_IT to use it in the second frac step
		
		SavePressureIt();
		Timer::Start("SolveStep2");
		double Dp_norm = this->SolveStep2();
		Timer::Stop("SolveStep2");
		//if(this->mReformDofAtEachIteration == true )
		//		this->Clear();
		Timer::Start("SolveStep3");
		SolveStep3();
		Timer::Stop("SolveStep3");

		//return 0.0;//
		return Dp_norm;
		
		}
				
		//*********************************************************************************
		//**********************************************************************
		void ApplyVelocityBoundaryConditions(DofsArrayType& mFixedVelocityDofSet,std::vector<double>& mFixedVelocityDofValues)
		{	
			KRATOS_TRY
				
			unsigned int i=0;
			for(typename DofsArrayType::iterator i_dof = mFixedVelocityDofSet.begin() ; i_dof != mFixedVelocityDofSet.end() ; ++i_dof)
			{
				i_dof->GetSolutionStepValue() = mFixedVelocityDofValues[i];
				i++;
			}
		
			KRATOS_CATCH("")
		}
		//*********************************************************************************
		//**********************************************************************
		void SetToZero( Variable<array_1d<double,3> >& rVariable, ModelPart::NodesContainerType& rNodes)
		{
			KRATOS_TRY
			array_1d<double,3> zero = ZeroVector(3);
			for(ModelPart::NodesContainerType::iterator i = rNodes.begin(); i!=rNodes.end(); i++)
					noalias(i->FastGetSolutionStepValue(rVariable)) = zero;
			KRATOS_CATCH("")
		}
		//*********************************************************************************
		//**********************************************************************

		
		void SetToZero( Variable<  double >& rVariable, ModelPart::NodesContainerType& rNodes)
		{
			KRATOS_TRY
			for(ModelPart::NodesContainerType::iterator i = rNodes.begin(); i!=rNodes.end(); i++)
					i->FastGetSolutionStepValue(rVariable) = 0.0;
			KRATOS_CATCH("")
		}
		//************************************************************************
		//************************************************************************
		
		//*************************************************************************
		//*************************************************************************
		void SolveStep1()		
		{
			KRATOS_TRY
			//vector that we shall use to store temporary results in the Runge-Kutta context
			array_1d<double,3> aux;	

			ModelPart& model_part=BaseType::GetModelPart();

			//getting delta time 
			//ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
			double& delta_t = model_part.GetProcessInfo()[DELTA_TIME];
			//KRATOS_WATCH(delta_t)
			//Write the mass factors (that corresponds to the assembly of global lumped mass matrix... instead of building and
			//storing this matris, we write the mass factors to the nodes in the var MASS_FACTORS
			for(ModelPart::NodeIterator in = BaseType::GetModelPart().NodesBegin() ; 
					in != BaseType::GetModelPart().NodesEnd() ; ++in)
			{
				//reset the mass factors to Zero			
				in->FastGetSolutionStepValue(NODAL_MASS)=0.0;
			}
			double dummy;
			ProcessInfo& proc_info = model_part.GetProcessInfo();

			for(ModelPart::ElementIterator im = model_part.ElementsBegin() ; 
				im != model_part.ElementsEnd() ; ++im)
			{
				//note that the lumped mass factors are saved nodally, and are equal for x, y, and z velocity				
				im->Calculate(NODAL_MASS, dummy, proc_info);				
			}				

			//////////////////////////////////////////////////////////////////////////////////////////////////
			//												//
			//			SAVING VELOCITY B.C.'s							//
			//												//
			//////////////////////////////////////////////////////////////////////////////////////////////////
			//if the DOFsets are reformed at every step, we have to find the Dirichlet B.C.s every time again (think of coupling with Lag)
			//otherwise do nothing
			if (this->mReformDofAtEachIteration==true)
			{
				mFixedVelocityDofSet.clear();
				mFixedVelocityDofValues.clear();

				for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
				{
					if(it->pGetDof(VELOCITY_X)->IsFixed() == true)
					{
						mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_X) );
						mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_X)->GetSolutionStepValue() );
					}

					if(it->pGetDof(VELOCITY_Y)->IsFixed() == true)
					{
						mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_Y) );
						mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_Y)->GetSolutionStepValue() );
					}

					if(it->pGetDof(VELOCITY_Z)->IsFixed() == true)
					{
						mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_Z) );
						mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_Z)->GetSolutionStepValue() );
					}
				}
			}
			
			//perform the First Fractional Step (using Runge-Kutta for finding u_tilda)
			//important is to apply the boundary conditions upon the intermediate velocity
		
			//set WORK = VELOCITY of the old step
			for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
			{ noalias(it->FastGetSolutionStepValue(AUX_VECTOR)) = it->FastGetSolutionStepValue(VELOCITY,1);	}
			
			//and set the initial guess for the velocity (for computation of teh first residual) to the one of teh old step
			//for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
			//{ noalias(it->FastGetSolutionStepValue(VELOCITY)) = it->FastGetSolutionStepValue(VELOCITY,1);	}


			//reset the RHS
			SetToZero(RHS_VECTOR,model_part.Nodes());
			
			array_1d<double,3> Frac_Step_Switch; //switch variable, that needs to be passed to the fct Calculate (of the element)
						// which decides weather 1st or last Frac Step(correction of vel) should be performed
						//if Frac_Step_Switch==1 - 1st  if 2-last, otherwise - ERROR

			Frac_Step_Switch[0]=1.0;
			Frac_Step_Switch[1]=1.0;
			Frac_Step_Switch[2]=1.0;
			//this variable is an array_1d just because of the structure of the function calculate....[1] and [2] aren't important

			//ProcessInfo& proc_info = model_part.GetProcessInfo();

			//loop over elements calculating the Right Hand Side, that is stored directly to the node.. this is done by fct Calculate
			for(ModelPart::ElementIterator im = model_part.ElementsBegin() ; im != model_part.ElementsEnd() ; ++im)
			{
			//compute the momentum residual, add it to the RHS_VECTOR on nodes
			im->Calculate(VELOCITY, Frac_Step_Switch, proc_info);
			}
			//////////////////////////////////////////
			//	first step of Runge Kutta	//
			//////////////////////////////////////////
			//KRATOS_WATCH("RUNGE KUTTA 1st STEP")
			double one_sixt = 0.166666666666667; 
			for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
			{
				//AUX_VECTOR = AUX_VECTOR + delta_T/6 * 1/NODAL_MASS * RHS
				//VELOCITY = VELOCITY_old + delta_T/2 * 1/NODAL_MASS * RHS
				noalias(aux) = delta_t/(it->FastGetSolutionStepValue(NODAL_MASS)) * it->FastGetSolutionStepValue(RHS_VECTOR);
				noalias(it->FastGetSolutionStepValue(AUX_VECTOR)) += one_sixt * aux;
				noalias(it->FastGetSolutionStepValue(VELOCITY)) += 0.5 * aux;			
			}
			ApplyVelocityBoundaryConditions(mFixedVelocityDofSet,mFixedVelocityDofValues);
			//apply the slip BC only if there are some slip BCs identified 	
			if (mSlipBoundaryList.size()!=0)
				ApplySlipBC();

			SetToZero(RHS_VECTOR,model_part.Nodes());

			/////////////////////////////////	
			//  second step of Runge Kutta //
			////////////////////////////////
			//KRATOS_WATCH("RUNGE KUTTA 2st STEP")
			//...now the residual will be computed with the intermediate velocity, that was computed at 1st step of R-K
			//loop over elements calculating the Right Hand Side, that is stored directly to the node
			for(ModelPart::ElementIterator im = model_part.ElementsBegin() ; im != model_part.ElementsEnd() ; ++im)
			{
			im->Calculate(VELOCITY, Frac_Step_Switch, proc_info);
			}
		
			double one_third = 0.33333333333333333333333333; 
			for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
			{
				//AUX_VECTOR = AUX_VECTOR + delta_T/3 * 1/NODAL_MASS * RHS
				//VELOCITY = VELOCITY_old + delta_T/2 * 1/NODAL_MASS * RHS
				noalias(aux) = delta_t/(it->FastGetSolutionStepValue(NODAL_MASS)) * it->FastGetSolutionStepValue(RHS_VECTOR);
				noalias(it->FastGetSolutionStepValue(AUX_VECTOR)) += one_third * aux;
				noalias(it->FastGetSolutionStepValue(VELOCITY)) += 0.5 * aux;		
			}
			ApplyVelocityBoundaryConditions(mFixedVelocityDofSet,mFixedVelocityDofValues);
			//apply the slip BC only if there are some slip BCs identified 	
			if (mSlipBoundaryList.size()!=0)
				ApplySlipBC();

			SetToZero(RHS_VECTOR,model_part.Nodes());

			///////////////////////////////
			//third step of Runge Kutta  //
			///////////////////////////////
			//KRATOS_WATCH("RUNGE KUTTA 3rd STEP")
			//...now the residual will be computed with the intermediate velocity, that was computed at 2nd step of R-K
			//loop over elements calculating the Right Hand Side, that is stored directly to the node
			for(ModelPart::ElementIterator im = model_part.ElementsBegin() ; im != model_part.ElementsEnd() ; ++im)
			{
			im->Calculate(VELOCITY, Frac_Step_Switch, proc_info);				
			}
			for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
			{
				//AUX_VECTOR = AUX_VECTOR + delta_T/3 * 1/NODAL_MASS * RHS
				//VELOCITY = VELOCITY_old + delta_T* 1/NODAL_MASS * RHS
				noalias(aux) = delta_t/(it->FastGetSolutionStepValue(NODAL_MASS)) * it->FastGetSolutionStepValue(RHS_VECTOR);
				noalias(it->FastGetSolutionStepValue(AUX_VECTOR)) += one_third * aux;
				noalias(it->FastGetSolutionStepValue(VELOCITY)) += aux;		
			}
			ApplyVelocityBoundaryConditions(mFixedVelocityDofSet,mFixedVelocityDofValues);
			//apply the slip BC only if there are some slip BCs identified 	
			if (mSlipBoundaryList.size()!=0)
				ApplySlipBC();

			SetToZero(RHS_VECTOR,model_part.Nodes());
			/////////////////////////////
			//last step of Runge Kutta //
			////////////////////////////
			//KRATOS_WATCH("RUNGE KUTTA LAst STEP")
			//...now the residual will be computed with the intermediate velocity, that was computed at 3rd step of R-K
			//loop over elements calculating the Right Hand Side, that is stored directly at the node
			for(ModelPart::ElementIterator im = model_part.ElementsBegin() ; im != model_part.ElementsEnd() ; ++im)
			{
			im->Calculate(VELOCITY, Frac_Step_Switch, proc_info);
			}
			for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
			{
				//AUX_VECTOR = AUX_VECTOR + delta_T/6 * 1/NODAL_MASS * RHS
				//VELOCITY = AUX_VSetToZero_VectorVarECTOR 
				noalias(aux) = delta_t/(it->FastGetSolutionStepValue(NODAL_MASS)) * it->FastGetSolutionStepValue(RHS_VECTOR);

				//KRATOS_WATCH(it->FastGetSolutionStepValue(RHS_VECTOR))

				noalias(it->FastGetSolutionStepValue(AUX_VECTOR)) += one_sixt * aux;

				//KRATOS_WATCH(it->FastGetSolutionStepValue(AUX_VECTOR))

				noalias(it->FastGetSolutionStepValue(VELOCITY)) = it->FastGetSolutionStepValue(AUX_VECTOR);			

				//KRATOS_WATCH(it->FastGetSolutionStepValue(VELOCITY))
			}
			ApplyVelocityBoundaryConditions(mFixedVelocityDofSet,mFixedVelocityDofValues);
			//apply the slip BC only if there are some slip BCs identified 	
			if (mSlipBoundaryList.size()!=0)
				ApplySlipBC();


			KRATOS_WATCH("FINISHED STAGE1 OF FRACTIONAL STEP")
			
			KRATOS_CATCH("")
		}

		//*********************************************************************************
		//**********************************************************************
		//solve the pressure equation
		double SolveStep2()
		{
		KRATOS_TRY
		KRATOS_WATCH("Second stage of Frac Step")
		
		//solves the system that is assembled within "calculateLocalSystem" of the element
		return mpressurestep->Solve();	
		KRATOS_CATCH("");
		}
		
	
		
		//******************************************************************************************************
		//******************************************************************************************************
		
		void SolveStep3()
		{
		KRATOS_TRY
					
		ModelPart& model_part=BaseType::GetModelPart();
		
		const double dt = model_part.GetProcessInfo()[DELTA_TIME];
				
		//set to zero AUX vector
		SetToZero(AUX_VECTOR,model_part.Nodes());
		for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
		{
			it->FastGetSolutionStepValue(AUX_VECTOR)=ZeroVector(3);
		}		

		if (TDim==2)
		{
			//allocation of work space
			boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
			array_1d<double,3> N;
			//array_1d<double,3> aux0, aux1, aux2; //this are sized to 3 even in 2D!!		
	//		double lumping_factor = 0.33333333333333;
		
		
			//calculate the velocity correction and store it in AUX_VECTOR
			for (typename ModelPart::ElementsContainerType::iterator it=model_part.ElementsBegin(); it!=model_part.ElementsEnd(); ++it)
			{
				//get the list of nodes of the element
				Geometry< Node<3> >& geom = it->GetGeometry();

				double volume;
				GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);			
						
				array_1d<double,3> pres_inc;
				pres_inc[0] = geom[0].FastGetSolutionStepValue(PRESSURE,1)-geom[0].FastGetSolutionStepValue(PRESSURE);
				pres_inc[1] = geom[1].FastGetSolutionStepValue(PRESSURE,1)-geom[1].FastGetSolutionStepValue(PRESSURE);
				pres_inc[2] = geom[2].FastGetSolutionStepValue(PRESSURE,1)-geom[2].FastGetSolutionStepValue(PRESSURE);
				
				//Riccardo's modification: multiply the G(p_n+1-p_n) by 1/2
				pres_inc*=0.5;
				//KRATOS_WATCH(pres_inc)

				//Gradient operator G:
				boost::numeric::ublas::bounded_matrix<double,6,2> shape_func = ZeroMatrix(6, 2);
				boost::numeric::ublas::bounded_matrix<double,6,3> G = ZeroMatrix(6,3);
				for (int ii = 0; ii< 3; ii++)
				    {
					int column = ii*2;				
					shape_func(column,0) = N[ii];
					shape_func(column + 1, 1) = shape_func(column,0);
				    }
				noalias(G)=prod(shape_func, trans(DN_DX));
				G*=volume;

				array_1d<double,6> aaa;
				noalias(aaa) = prod(G,pres_inc);

				array_1d<double,3> aux;
				aux[0]=aaa[0];
				aux[1]=aaa[1];			
				//z-component is zero
				aux[2]=0.0;

				geom[0].FastGetSolutionStepValue(AUX_VECTOR) += aux;
				//reusing aux for the second node 
				aux[0]=aaa[2];
				aux[1]=aaa[3];			
				//z-component is zero
				geom[1].FastGetSolutionStepValue(AUX_VECTOR) += aux;
				//reusing aux for the third node
				aux[0]=aaa[4];
				aux[1]=aaa[5];			
				geom[2].FastGetSolutionStepValue(AUX_VECTOR) += aux;
				//for(unsigned int i=0;i<3;i++)
				//  geom[i].FastGetSolutionStepValue(AUX_VECTOR) += aux0;
			}
		}
		if (TDim==3)
			{
			KRATOS_WATCH("Last step in 3D")
			
			array_1d<double,4> pres_inc;
			boost::numeric::ublas::bounded_matrix<double,12,3> shape_func = ZeroMatrix(12, 3);
			boost::numeric::ublas::bounded_matrix<double,12,4> G = ZeroMatrix(12,4);
			boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;
			array_1d<double,4> N;
			//array_1d<double,3> aux0, aux1, aux2, aux3; //this are sized to 3 even in 2D!!		
			for (typename ModelPart::ElementsContainerType::iterator it=model_part.ElementsBegin(); it!=model_part.ElementsEnd(); ++it)
			{
				Geometry< Node<3> >& geom = it->GetGeometry();

				pres_inc[0] = geom[0].FastGetSolutionStepValue(PRESSURE,1)-geom[0].FastGetSolutionStepValue(PRESSURE);
				pres_inc[1] = geom[1].FastGetSolutionStepValue(PRESSURE,1)-geom[1].FastGetSolutionStepValue(PRESSURE);
				pres_inc[2] = geom[2].FastGetSolutionStepValue(PRESSURE,1)-geom[2].FastGetSolutionStepValue(PRESSURE);
				pres_inc[3] = geom[3].FastGetSolutionStepValue(PRESSURE,1)-geom[3].FastGetSolutionStepValue(PRESSURE);
				
				//Riccardo's modification: multiply the G(p_n+1-p_n) by 1/2
				pres_inc*=0.5;
	
				double volume;
				GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);		

				//Gradient operator G:
			
				for (int ii = 0; ii< 4; ii++)
				    {
					int column = ii*3;				
					shape_func(column,0) = N[ii];
					shape_func(column + 1, 1) = shape_func(column,0);
					shape_func(column + 2, 2) = shape_func(column,0);
				    }
				noalias(G)=prod(shape_func, trans(DN_DX));
				G*=volume;

				array_1d<double,12> aaa;
				noalias(aaa) = prod(G,pres_inc);

				array_1d<double,3> aux;
				aux[0]=aaa[0];
				aux[1]=aaa[1];			
				aux[2]=aaa[2];			

				geom[0].FastGetSolutionStepValue(AUX_VECTOR) += aux;
				//reusing aux for the second node 
				aux[0]=aaa[3];
				aux[1]=aaa[4];
				aux[2]=aaa[5];						
				//z-component is zero
				geom[1].FastGetSolutionStepValue(AUX_VECTOR) += aux;
				//reusing aux for the third node
				aux[0]=aaa[6];
				aux[1]=aaa[7];
				aux[2]=aaa[8];						
				geom[2].FastGetSolutionStepValue(AUX_VECTOR) += aux;

				aux[0]=aaa[9];
				aux[1]=aaa[10];
				aux[2]=aaa[11];						
				geom[3].FastGetSolutionStepValue(AUX_VECTOR) += aux;

				//for(unsigned int i=0;i<3;i++)
				//  geom[i].FastGetSolutionStepValue(AUX_VECTOR) += aux0;
			}
			}
		//correct the velocities
		for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
		{
			//VELOCITY = VELOCITY + dt * Minv * AUX_VECTOR
			double dt_Minv = dt / it->FastGetSolutionStepValue(NODAL_MASS);
			array_1d<double,3>& temp = it->FastGetSolutionStepValue(AUX_VECTOR);		
			if(!it->IsFixed(VELOCITY_X))
				{
				it->FastGetSolutionStepValue(VELOCITY_X)+=dt_Minv*temp[0];
				}
			if(!it->IsFixed(VELOCITY_Y))
				{
				it->FastGetSolutionStepValue(VELOCITY_Y)+=dt_Minv*temp[1];						
				}
			if(!it->IsFixed(VELOCITY_Z))
				{
				it->FastGetSolutionStepValue(VELOCITY_Z)+=dt_Minv*temp[2];						
				}
			
			
		}
		
		KRATOS_CATCH("");
		}
		
		void SavePressureIt()
		{
			KRATOS_TRY
			
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
				//setting the old value of the pressure to the current one
				const double& p = (i)->FastGetSolutionStepValue(PRESSURE);
				//(i)->FastGetSolutionStepValue(PRESSURE,1) = p;
				(i)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = p;

			}	
			KRATOS_CATCH("")
		}
		//************************************
		//************************************
		void ApplySlipBC()
		{
		KRATOS_TRY
		
		
		for (PointIterator it=mSlipBoundaryList.begin(); it!=mSlipBoundaryList.end(); ++it)
		{
		//KRATOS_WATCH("slip node")
		array_1d<double, 3> normal = (*it)->FastGetSolutionStepValue(NORMAL);
		double length = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
		if (length==0)
			{
			KRATOS_ERROR(std::logic_error,  "TO apply SLIP you should calculate normals first! Dont forget to assign Condition2D/3D resp for that " , "");
			}
		normal*=1.0/length;
		array_1d<double, 3> normal_comp_vec;
		//CHECK IF NORMAL IS NORMALIZED (divided by the length)
		//double length = ...
		array_1d<double, 3> vel = (*it)->FastGetSolutionStepValue(VELOCITY);
		double normal_comp;
		normal_comp=inner_prod(normal, vel);
		normal_comp_vec = normal_comp*normal;
		(*it)->FastGetSolutionStepValue(VELOCITY)-=normal_comp_vec;
		}					
		KRATOS_CATCH("")
		}


		//******************************************************************************************************
		//******************************************************************************************************
		virtual void SetEchoLevel(int Level) 
		{
			//mfracvel_x_strategy->SetEchoLevel(Level);
			//mfracvel_y_strategy->SetEchoLevel(Level);
			//if(mdomain_size == 3)
			//	mfracvel_z_strategy->SetEchoLevel(Level);
			//
			mpressurestep->SetEchoLevel(Level);
		}

		//******************************************************************************************************
		//******************************************************************************************************
		virtual void Clear() 
		{
			KRATOS_WATCH("RungeKuttaFractStepGLSStrategy Clear Function called");
			//mfracvel_x_strategy->Clear();
			//mfracvel_y_strategy->Clear();
			//if(mdomain_size == 3)
			//	mfracvel_z_strategy->Clear();
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
		
		//typename BaseType::Pointer mfracvel_x_strategy;
		//typename BaseType::Pointer mfracvel_y_strategy;
		//typename BaseType::Pointer mfracvel_z_strategy;
		typename BaseType::Pointer mpressurestep;

		//double mvelocity_toll;
		double mpressure_toll;
		//unsigned int mdomain_size;
		//int mMaxVelIterations;
		//int mMaxPressIterations;
		//unsigned int mtime_order;
		//unsigned int mprediction_order;
		//bool mpredictor_corrector;
		bool mReformDofAtEachIteration;

		//double mmin_conv_vel_norm;


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
		/*		
		unsigned int m_step;
		int mdomain_size;
		double mOldDt;
		*/
		//ModelPart& model_part;
	
		DofsArrayType mFixedVelocityDofSet;
		std::vector<double> mFixedVelocityDofValues;
		PointVector mSlipBoundaryList;
		//bool proj_is_initialized;


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
		RungeKuttaFracStepStrategy(const RungeKuttaFracStepStrategy& Other);


		/*@} */   

	}; /* Class RungeKuttaFracStepStrategy */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RUNGE_KUTTA_GLS_STRATEGY  defined */

