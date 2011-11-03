/* *********************************************************   
*          
*   Last Modified by:    $Author: jmarti $
*   Date:                $Date: 2008-11-10 14:23:32 $
*   Revision:            $Revision: 1.12 $
*
* ***********************************************************/


#if !defined(KRATOS_GLS_STRATEGY)
#define  KRATOS_GLS_STRATEGY


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
	class FracStepStrategy 
		: public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
	{
	public:
		/**@name Type Definitions */       
		/*@{ */

		/** Counted pointer of ClassName */
		typedef std::vector<unsigned int> IndicesVectorType;
		
		typedef boost::shared_ptr< FracStepStrategy<TDim, TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;

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
		FracStepStrategy(
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
				KRATOS_ERROR(std::logic_error,  "Incompressible Strategy requires utilization of Fluid2DGLS_expl elements " , "");
			}

			if (TDim==3)
			{
			Fluid3DGLS_expl el(1, p_null_geom);

			if (typeid(ref_el) != typeid(el))
				KRATOS_ERROR(std::logic_error,  "Incompressible Strategy requires utilization of Fluid3DGLS_expl elements " , "");
			}

			KRATOS_CATCH("")
		}



		/** Destructor.
		*/
		virtual ~FracStepStrategy() {}

		/** Destructor.
		*/

		//*********************************************************************************
		//**********************************************************************
		double Solve()
		{
		KRATOS_WATCH("Solve of GLS Frac Step Strategy")
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

		SaveAccelerations();
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
			//vector that we shall use to store temporary results in the context
			array_1d<double,3> aux;	
			array_1d<double,3> aux1;	
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
			
			array_1d<double,3> zero = ZeroVector(3);
			//set WORK = VELOCITY of the old step
			for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
			{ 
			noalias(it->FastGetSolutionStepValue(AUX_VECTOR)) = it->FastGetSolutionStepValue(VELOCITY,1);	
			noalias(it->FastGetSolutionStepValue(FORCE)) =    zero;		
			}
			
			
			//reset the RHS
			SetToZero(RHS_VECTOR,model_part.Nodes());
			
			array_1d<double,3> Frac_Step_Switch; 

			Frac_Step_Switch[0]=1.0;
			Frac_Step_Switch[1]=1.0;
			Frac_Step_Switch[2]=1.0;


			//loop over elements calculating the Right Hand Side, that is stored directly to the node.. this is done by fct Calculate
			for(ModelPart::ElementIterator im = model_part.ElementsBegin() ; im != model_part.ElementsEnd() ; ++im)
			{
			//compute the momentum residual, add it to the RHS_VECTOR on nodes
			im->Calculate(VELOCITY, Frac_Step_Switch, proc_info);
			}
			
			
			for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)	{	
			if(it->pGetDof(VELOCITY_X)->IsFixed() == true)
				{
				noalias(it->FastGetSolutionStepValue(FORCE)) =    zero;	
				}
			}

		

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
			it->FastGetSolutionStepValue(FORCE)=ZeroVector(3);
		}		

		if (TDim==2)
		{
			//allocation of work space
			boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
			array_1d<double,3> N;
		
		
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
				geom[0].FastGetSolutionStepValue(FORCE) += aux;
				//reusing aux for the second node 
				aux[0]=aaa[2];
				aux[1]=aaa[3];			
				//z-component is zero
				geom[1].FastGetSolutionStepValue(AUX_VECTOR) += aux;
				geom[1].FastGetSolutionStepValue(FORCE) += aux;
				//reusing aux for the third node
				aux[0]=aaa[4];
				aux[1]=aaa[5];			
				geom[2].FastGetSolutionStepValue(AUX_VECTOR) += aux;
				geom[2].FastGetSolutionStepValue(FORCE) += aux;
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
			double dt_Minv = (dt / 2.00) / it->FastGetSolutionStepValue(NODAL_MASS);
			array_1d<double,3>& temp = it->FastGetSolutionStepValue(AUX_VECTOR);
			array_1d<double,3>& force_temp = it->FastGetSolutionStepValue(FORCE);
			force_temp *=(1.0/ it->FastGetSolutionStepValue(NODAL_MASS));
			//KRATOS_WATCH(force_temp);

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
		void SaveAccelerations()
		{
		KRATOS_TRY
		array_1d<double, 3> acc=ZeroVector(3);
		ModelPart& model_part=BaseType::GetModelPart();
		const double dt = model_part.GetProcessInfo()[DELTA_TIME];
		for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
				i != BaseType::GetModelPart().NodesEnd() ; ++i)
			{
			acc=(i)->FastGetSolutionStepValue(VELOCITY)-(i)->FastGetSolutionStepValue(VELOCITY,1);
			(i)->FastGetSolutionStepValue(ACCELERATION)=acc/dt;
			}
		KRATOS_CATCH("")
		}
		///////////////////////////////////////////
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
			mpressurestep->SetEchoLevel(Level);
		}

		//******************************************************************************************************
		//******************************************************************************************************
		virtual void Clear() 
		{
			KRATOS_WATCH("FractStepGLSStrategy Clear Function called");
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
		FracStepStrategy(const FracStepStrategy& Other);


		/*@} */   

	}; /* Class FracStepStrategy */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_GLS_STRATEGY  defined */


