//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Pavel Ryzhakov and Julio Marti 
//


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
//#include "incompressible_fluid_application.h"
#include "ULF_application.h"
#include "includes/cfd_variables.h"
#include "includes/mat_variables.h"
#include "includes/dem_variables.h"


#include <stdio.h>      /* printf */
#include <math.h>       /* cos */

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

#define UPDATE_MASS //chek if the element is correct

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

      KRATOS_CLASS_POINTER_DEFINITION(  RungeKuttaFracStepStrategy );

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
		  mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_X).get() );
		  mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_X)->GetSolutionStepValue() );
		}

	      if(it->pGetDof(VELOCITY_Y)->IsFixed() == true)
		{
		  mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_Y).get() );
		  mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_Y)->GetSolutionStepValue() );
		}

	      if(it->pGetDof(VELOCITY_Z)->IsFixed() == true)
		{
		  mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_Z).get() );
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
                KRATOS_THROW_ERROR(std::logic_error,  "Incompressible Runge Kutta Strategy requires utilization of Fluid2DGLS_expl elements " , "");
	    }

	  if (TDim==3)
	    {
	      KRATOS_THROW_ERROR(std::logic_error,  "not Runge Kutta Strategy for 3D problems " , "");
	      //            Fluid3DGLS_expl el(1, p_null_geom);

	      //            if (typeid(ref_el) != typeid(el))
	      //                KRATOS_THROW_ERROR(std::logic_error,  "Incompressible Runge Kutta Strategy requires utilization of Fluid3DGLS_expl elements " , "");
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
	  KRATOS_WATCH("Started RK step1");
        Timer::Start("SolveStep1");
	SolveStep1();
        Timer::Stop("SolveStep1");
	KRATOS_WATCH("Finished RK step1");
	
	Timer::Start("SolveStep2");
        double Dp_norm = this->SolveStep2();
	Timer::Stop("SolveStep2");
	KRATOS_WATCH("Finished RK step2");
        
	if(this->mReformDofAtEachIteration == true )
	  this->Clear();
        Timer::Start("SolveStep3");
	SolveStep3();
        Timer::Stop("SolveStep3");
	KRATOS_WATCH("Finished RK step3");
	
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
      void SolveStep1()
      {
	KRATOS_TRY
	  //vector that we shall use to store temporary results in the Runge-Kutta context
	  array_1d<double,3> aux;
	array_1d<double,3> aux_res;
      
	ModelPart& model_part=BaseType::GetModelPart();
      
	//getting delta time
	//ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
	double delta_t = model_part.GetProcessInfo()[DELTA_TIME];
	double dummy;
	ProcessInfo& proc_info = model_part.GetProcessInfo();
	bool inverted_element=true;

	//inverted_element=CalculateLumpedMassAux();

	double time_ccc = model_part.GetProcessInfo()[TIME];
	double time_cal=time_ccc - delta_t;
      
	for(ModelPart::NodesContainerType::iterator in = model_part.NodesBegin() ;  in != model_part.NodesEnd() ; ++in)
          {
	    in->FastGetSolutionStepValue(IS_INTERFACE)=0.0; 
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
		    mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_X).get() );
		    mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_X)->GetSolutionStepValue() );
		  }

		if(it->pGetDof(VELOCITY_Y)->IsFixed() == true)
		  {
		    mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_Y).get() );
		    mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_Y)->GetSolutionStepValue() );
		  }

		if(it->pGetDof(VELOCITY_Z)->IsFixed() == true)
		  {
		    mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_Z).get() );
		    mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_Z)->GetSolutionStepValue() );
		  }
	      }
	  }

	//perform the First Fractional Step (using Runge-Kutta for finding u_tilda)
	//important is to apply the boundary conditions upon the intermediate velocity
	//set WORK = VELOCITY of the old step

	bool inverted=false;
	while(inverted_element==true )
	  {	
	    inverted_element=false;
	  
	    for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	      {
		noalias(it->FastGetSolutionStepValue(VELOCITY_OLD)) = it->FastGetSolutionStepValue(VELOCITY,1);
		noalias(it->FastGetSolutionStepValue(CONV_PROJ)) = ZeroVector(3);
		noalias(it->FastGetSolutionStepValue(CONV_PROJ,1)) = ZeroVector(3);
		it->FastGetSolutionStepValue(PRESSURE) = it->FastGetSolutionStepValue(PRESSURE,1);	
		it->FastGetSolutionStepValue(VELOCITY) = it->FastGetSolutionStepValue(VELOCITY,1);
		it->FastGetSolutionStepValue(NODAL_MASS)=0.0;		
	      }


	  
	    for(ModelPart::ElementIterator im = model_part.ElementsBegin() ;
		im != model_part.ElementsEnd() ; ++im)
	      {
		//note that the lumped mass factors are saved nodally, and are equal for x, y, and z velocity
		im->Calculate(NODAL_MASS, dummy, proc_info);
	      }
	  
	    //reset the RHS
	    SetToZero(VELOCITY_OLD_OLD,model_part.Nodes());
	  
	    array_1d<double,3> Frac_Step_Switch; //switch variable, that needs to be passed to the fct Calculate (of the element)
	    // which decides weather 1st or last Frac Step(correction of vel) should be performed
	    //if Frac_Step_Switch==1 - 1st  if 2-last, otherwise - ERROR

	    Frac_Step_Switch[0]=1.0;
	    Frac_Step_Switch[1]=1.0;
	    Frac_Step_Switch[2]=1.0;
	    //this variable is an array_1d just because of the structure of the function calculate....[1] and [2] aren't important

	    //ProcessInfo& proc_info = model_part.GetProcessInfo();

#if defined( UPDATE_MASS)
	    inverted=CalculateLumpedMassAux();
	    if(inverted==true ) inverted_element=true;
#endif

	    double cal_time=time_cal;


	    //loop over elements calculating the Right Hand Side, that is stored directly to the node.. this is done by fct Calculate
	    for(ModelPart::ElementIterator im = model_part.ElementsBegin() ; im != model_part.ElementsEnd() ; ++im)
	      {
		//compute the momentum residual, add it to the VELOCITY_OLD_OLD on nodes
		im->Calculate(VELOCITY, Frac_Step_Switch, proc_info);
	      }
	    //////////////////////////////////////////
	    //	first step of Runge Kutta	//
	    //////////////////////////////////////////
	    //KRATOS_WATCH("RUNGE KUTTA 1st STEP")
	    double one_sixt = 0.166666666666667;
	    for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	      {
		noalias(it->FastGetSolutionStepValue(CONV_PROJ))= one_sixt * delta_t * it->FastGetSolutionStepValue(VELOCITY,1);

		noalias(aux) = delta_t/(it->FastGetSolutionStepValue(NODAL_MASS)) * it->FastGetSolutionStepValue(VELOCITY_OLD_OLD);
		noalias(it->FastGetSolutionStepValue(VELOCITY_OLD)) += one_sixt * aux;
		noalias(it->FastGetSolutionStepValue(VELOCITY)) = 0.5 * aux + it->FastGetSolutionStepValue(VELOCITY,1); //V_beta1
	       
	      }
	    ApplyVelocityBoundaryConditions(mFixedVelocityDofSet,mFixedVelocityDofValues);
	    //apply the slip BC only if there are some slip BCs identified
	    if (mSlipBoundaryList.size()!=0)
	      ApplySlipBC();

	    SetToZero(VELOCITY_OLD_OLD,model_part.Nodes());

	    //updating position
	    for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	      {
		it->FastGetSolutionStepValue(ADVPROJ)= it->Coordinates();


		if(it->FastGetSolutionStepValue(IS_STRUCTURE)==0.0)
		  {
		    noalias(it->Coordinates())+= 0.5 * delta_t * it->FastGetSolutionStepValue(VELOCITY,1);  //move with vn

		  }
	      }
	    
	    /////////////////////////////////
	    //  second step of Runge Kutta //
	    ////////////////////////////////
	    //KRATOS_WATCH("RUNGE KUTTA 2st STEP")
	    
	    //...now the residual will be computed with the intermediate velocity, that was computed at 1st step of R-K
	    //loop over elements calculating the Right Hand Side, that is stored directly to the node
	    
	    //LETS TRY NOT TO RECCOMPUTE MASS MATRIX	  
#if defined( UPDATE_MASS)
	    inverted=CalculateLumpedMassAux();
	    if(inverted==true ) inverted_element=true;
	    
#endif
	    
	    cal_time=time_cal + 0.5 * delta_t;
	    
	    for(ModelPart::ElementIterator im = model_part.ElementsBegin() ; im != model_part.ElementsEnd() ; ++im)
	      {
		im->Calculate(VELOCITY, Frac_Step_Switch, proc_info);
	      }
	    
	    //updating position_ x_beta2=xn+0.5*dt*v_beta1
	    for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	      {
		//resetting coordinates to x_n
		it->Coordinates()=it->FastGetSolutionStepValue(ADVPROJ);
		if(it->FastGetSolutionStepValue(IS_STRUCTURE)==0.0)
		  {
		    //x_beta2=x_n+v_beta1
		    noalias(it->Coordinates())+= 0.5 * delta_t * it->FastGetSolutionStepValue(VELOCITY);  
		   
		  }
	      }

	    double one_third = 0.33333333333333333333333333;

	    for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	      {
	      
		noalias(it->FastGetSolutionStepValue(CONV_PROJ))+= one_third * delta_t * it->FastGetSolutionStepValue(VELOCITY);	//v_betta2
		noalias(aux) = delta_t/(it->FastGetSolutionStepValue(NODAL_MASS)) * it->FastGetSolutionStepValue(VELOCITY_OLD_OLD);
		noalias(it->FastGetSolutionStepValue(VELOCITY_OLD)) += one_third * aux;
		//computing v_beta2=v_n+0.5dt*r
		noalias(it->FastGetSolutionStepValue(VELOCITY)) = 0.5 * aux  + it->FastGetSolutionStepValue(VELOCITY,1); //v beta2
	      }
	    
	    ApplyVelocityBoundaryConditions(mFixedVelocityDofSet,mFixedVelocityDofValues);
	    //apply the slip BC only if there are some slip BCs identified
	    if (mSlipBoundaryList.size()!=0)
	      ApplySlipBC();
	    
	    SetToZero(VELOCITY_OLD_OLD,model_part.Nodes());
	    
	    //LETS TRY NOT TO RECCOMPUTE MASS MATRIX	  
#if defined( UPDATE_MASS)
	    inverted=CalculateLumpedMassAux();
	    if(inverted==true ) inverted_element=true;
#endif	

	    ///////////////////////////////
	    //third step of Runge Kutta  //
	    ///////////////////////////////
	    //KRATOS_WATCH("RUNGE KUTTA 3rd STEP")
	    //...now the residual will be computed with the intermediate velocity, that was computed at 2nd step of R-K
	    //loop over elements calculating the Right Hand Side, that is stored directly to the node

	    cal_time=time_cal + 0.5 * delta_t;
	    
	    for(ModelPart::ElementIterator im = model_part.ElementsBegin() ; im != model_part.ElementsEnd() ; ++im)
	      {
		im->Calculate(VELOCITY, Frac_Step_Switch, proc_info);
	      }
	    
	    for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	      {

		if(it->FastGetSolutionStepValue(IS_STRUCTURE)==0.0)
		  {
		    //resetting coords to x_n
		    noalias(it->Coordinates())= it->FastGetSolutionStepValue(ADVPROJ);
		    //x_beta3=xn+dt*vbeta2
		    noalias(it->Coordinates())+= delta_t * it->FastGetSolutionStepValue(VELOCITY);
		  }
	      }
	    
	    for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	      {
		noalias(it->FastGetSolutionStepValue(CONV_PROJ))+= one_third * delta_t * it->FastGetSolutionStepValue(VELOCITY);
		noalias(aux) = delta_t/(it->FastGetSolutionStepValue(NODAL_MASS)) * it->FastGetSolutionStepValue(VELOCITY_OLD_OLD);
		noalias(it->FastGetSolutionStepValue(VELOCITY_OLD)) += one_third * aux;
		//v_beta3=vn+dt*r3
		noalias(it->FastGetSolutionStepValue(VELOCITY)) = aux + it->FastGetSolutionStepValue(VELOCITY,1);//aux;
		//NEW VERSION OF VELbeta3 - see blue text in the paper
	      }
	    ApplyVelocityBoundaryConditions(mFixedVelocityDofSet,mFixedVelocityDofValues);
	    //apply the slip BC only if there are some slip BCs identified
	    if (mSlipBoundaryList.size()!=0)
	      ApplySlipBC();
	    
	    SetToZero(VELOCITY_OLD_OLD,model_part.Nodes());
	    
	    cal_time=time_cal + delta_t;
	    
	    /////////////////////////////
	    //last step of Runge Kutta //
	    ////////////////////////////
	    
#if defined( UPDATE_MASS)
	    inverted=CalculateLumpedMassAux();
	    if(inverted==true ) inverted_element=true;
#endif
	    
	    //KRATOS_WATCH("RUNGE KUTTA LAst STEP")
	    //...now the residual will be computed with the intermediate velocity, that was computed at 3rd step of R-K
	    //loop over elements calculating the Right Hand Side, that is stored directly at the node
	    for(ModelPart::ElementIterator im = model_part.ElementsBegin() ; im != model_part.ElementsEnd() ; ++im)
	      {
		im->Calculate(VELOCITY, Frac_Step_Switch, proc_info);
	      }
	    
	    for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	      {
		
		noalias(it->FastGetSolutionStepValue(CONV_PROJ))+= one_sixt * delta_t * it->FastGetSolutionStepValue(VELOCITY); 
		//VELOCITY_OLD = VELOCITY_OLD + delta_T/6 * 1/NODAL_MASS * RHS
		//VELOCITY = AUX_VSetToZero_VectorVarECTOR
		noalias(aux) = delta_t/(it->FastGetSolutionStepValue(NODAL_MASS)) * it->FastGetSolutionStepValue(VELOCITY_OLD_OLD);
		noalias(it->FastGetSolutionStepValue(VELOCITY_OLD)) += one_sixt * aux;
	      }
	    
	    for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	      {
		noalias(it->Coordinates())=it->FastGetSolutionStepValue(ADVPROJ);
		if(it->FastGetSolutionStepValue(IS_STRUCTURE)==0.0)
		  {
		    noalias(it->Coordinates())+= it->FastGetSolutionStepValue(CONV_PROJ);
		  }
	      }
	    
	    for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	      {
		noalias(it->FastGetSolutionStepValue(VELOCITY)) = it->FastGetSolutionStepValue(VELOCITY_OLD); //  V tilde
	      }
	    ApplyVelocityBoundaryConditions(mFixedVelocityDofSet,mFixedVelocityDofValues);
	    //apply the slip BC only if there are some slip BCs identified
	    if (mSlipBoundaryList.size()!=0)
	      ApplySlipBC();
	    
#if defined( UPDATE_MASS)
	    inverted=CalculateLumpedMassAux();
	    if(inverted==true ) inverted_element=true;
#endif
	    
	    if(inverted_element==true) { 
	      KRATOS_WATCH("AN ELEMENT HAS BEEN INVERTEDDDDDDDDD ");
	      KRATOS_WATCH("INSIDE -REDUCE TIME STEP- process")		
		inverted_element=false;
	    }
	  }
	KRATOS_WATCH("FINISHED STAGE1 OF FRACTIONAL STEP")
	  KRATOS_CATCH("")
	  }
      
      //******************************************************************************************************
      //******************************************************************************************************
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
        SetToZero(VELOCITY_OLD,model_part.Nodes());
        for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	  {
            it->FastGetSolutionStepValue(VELOCITY_OLD)=ZeroVector(3);
	    it->FastGetSolutionStepValue(CONV_PROJ)=it->FastGetSolutionStepValue(VELOCITY);

	  }
        double one_sixth = 0.166666666666667;
	
        if (TDim==2)
	  {
            //allocation of work space
            boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
            array_1d<double,3> N;
            //array_1d<double,3> aux0, aux1, aux2; //this are sized to 3 even in 2D!!
            //double lumping_factor = 0.33333333333333;
	    
	    
            //calculate the velocity correction and store it in VELOCITY_OLD
            for (typename ModelPart::ElementsContainerType::iterator it=model_part.ElementsBegin(); it!=model_part.ElementsEnd(); ++it)
	      {
                //get the list of nodes of the element
                Geometry< Node<3> >& geom = it->GetGeometry();
                
                double volume;
                GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);
		// OLD VERSION PAVEL'S IMPLEMENTATION
                array_1d<double,3> pres_inc;
                pres_inc[0] = geom[0].FastGetSolutionStepValue(PRESSURE)-geom[0].FastGetSolutionStepValue(PRESSURE,1);
                pres_inc[1] = geom[1].FastGetSolutionStepValue(PRESSURE)-geom[1].FastGetSolutionStepValue(PRESSURE,1);
		pres_inc[2] = geom[2].FastGetSolutionStepValue(PRESSURE)-geom[2].FastGetSolutionStepValue(PRESSURE,1);

		pres_inc*=one_sixth;

                array_1d<double,3> aux=ZeroVector(2);

 		double p_avg = N[0]*pres_inc[0]
		  + N[1]*pres_inc[1]
		  + N[2]*pres_inc[2];
		
		p_avg *= volume;

		//KRATOS_WATCH(p_avg)

		array_1d<double,6> aaa=ZeroVector(6);

	        aaa[0] += DN_DX(0, 0) * p_avg;
        	aaa[1] += DN_DX(0, 1) * p_avg;
	        aaa[2] += DN_DX(1, 0) * p_avg;
        	aaa[3] += DN_DX(1, 1) * p_avg;
	        aaa[4] += DN_DX(2, 0) * p_avg;
        	aaa[5] += DN_DX(2, 1) * p_avg;

		//KRATOS_WATCH(aaa)

                geom[0].FastGetSolutionStepValue(VELOCITY_OLD_X) += aaa[0];
                geom[0].FastGetSolutionStepValue(VELOCITY_OLD_Y) += aaa[1];
                
                geom[1].FastGetSolutionStepValue(VELOCITY_OLD_X) += aaa[2];
                geom[1].FastGetSolutionStepValue(VELOCITY_OLD_Y) += aaa[3];

		geom[2].FastGetSolutionStepValue(VELOCITY_OLD_X) += aaa[4];
                geom[2].FastGetSolutionStepValue(VELOCITY_OLD_Y) += aaa[5];
                //reusing aux for the third node
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

                geom[0].FastGetSolutionStepValue(VELOCITY_OLD) += aux;
                //reusing aux for the second node
                aux[0]=aaa[3];
                aux[1]=aaa[4];
                aux[2]=aaa[5];
                //z-component is zero
                geom[1].FastGetSolutionStepValue(VELOCITY_OLD) += aux;
                //reusing aux for the third node
                aux[0]=aaa[6];
                aux[1]=aaa[7];
                aux[2]=aaa[8];
                geom[2].FastGetSolutionStepValue(VELOCITY_OLD) += aux;

                aux[0]=aaa[9];
                aux[1]=aaa[10];
                aux[2]=aaa[11];
                geom[3].FastGetSolutionStepValue(VELOCITY_OLD) += aux;

                //for(unsigned int i=0;i<3;i++)
                //  geom[i].FastGetSolutionStepValue(VELOCITY_OLD) += aux0;
	      }
	  }
        //correct the velocities
        for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	  {
            //VELOCITY = VELOCITY + dt * Minv * VELOCITY_OLD
            double dt_Minv = dt / it->FastGetSolutionStepValue(NODAL_MASS);
            array_1d<double,3>& temp = it->FastGetSolutionStepValue(VELOCITY_OLD);
            if(!it->IsFixed(VELOCITY_X))
	      {
                it->FastGetSolutionStepValue(VELOCITY_X)+= dt_Minv*temp[0];

	      }
            if(!it->IsFixed(VELOCITY_Y))
	      {
                it->FastGetSolutionStepValue(VELOCITY_Y)+= dt_Minv*temp[1];
	      }
            if(!it->IsFixed(VELOCITY_Z))
	      {
                it->FastGetSolutionStepValue(VELOCITY_Z)+= dt_Minv*temp[2];
	      }
	
	  }

        KRATOS_CATCH("");
      }
      //************************************
      //************************************
      void SavePressureIt()
      {
        KRATOS_TRY
	  
	  for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ;
	      i != BaseType::GetModelPart().NodesEnd() ; ++i)
	    {
	      //setting the old value of the pressure to the current one
	      const double& p = (i)->FastGetSolutionStepValue(PRESSURE);
	      
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
		  KRATOS_THROW_ERROR(std::logic_error,  "TO apply SLIP you should calculate normals first! Dont forget to assign Condition2D/3D resp for that " , "");
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
      //************************************
      //************************************
      void AssembleMassMatrices(TSystemMatrixType& Mconsistent, TSystemVectorType& mMdiagInv,  ModelPart& r_model_part)
      {
	//first we assemble the diagonal mass matrix
	KRATOS_TRY
	  //KRATOS_WATCH("BUILDING MASS MATRICES ")
	  boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
	array_1d<double,TDim+1> N;
	array_1d<unsigned int ,TDim+1> local_indices;
	//array_1d<double,TDim+1> rhs_contribution;
	double Volume;
	double temp;
	//getting the dof position
	unsigned int dof_position = (r_model_part.NodesBegin())->GetDofPosition(DISPLACEMENT_X);
      
	double aaa = 1.0/(TDim+1.0);
      
	for(ModelPart::ElementsContainerType::iterator i = r_model_part.ElementsBegin(); i!=r_model_part.ElementsEnd(); i++)
	  {
	  
	    Geometry< Node<3> >& geom = i->GetGeometry();
	    //counting number of structural nodes
	    unsigned int str_nr=0;
	    //for (int k = 0;k<TDim+1;k++)
	    for (unsigned int k = 0; k<geom.size(); k++)
	      {
		str_nr+=int(i->GetGeometry()[k].FastGetSolutionStepValue(IS_STRUCTURE));
	      }
	    //we do not do anything for the elements of the structure (all nodes are IS_STR)
	    if (geom.size()!=str_nr)
	      {
	      
		GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
		if (Volume<0)
		  Volume*=-1.0;
	      
		//finiding local indices
		//for(int ii = 0; ii<TDim+1; ii++)
		for(unsigned int ii = 0; ii<geom.size(); ii++)
		  {
		    local_indices[ii] = geom[ii].GetDof(DISPLACEMENT_X,dof_position).EquationId();
		  }
	      
		temp = Volume*aaa;
		for(unsigned int row = 0; row<TDim+1; row++)
		  {
                    unsigned int row_index = local_indices[row] / (TDim);
                    mMdiagInv[row_index] += temp;
		  }
	      }
	      
	  }
	//KRATOS_WATCH(mMdiagInv)
	//inverting the mass matrix
	for(unsigned int i = 0; i<TSparseSpace::Size(mMdiagInv); i++)
	  {
	    if (mMdiagInv[i]>1e-26)
	      mMdiagInv[i] = 1.0/mMdiagInv[i];
	    else  //if (mMdiagInv[i]==0.0)
	      {
	      
		//KRATOS_WATCH(mMdiagInv[i])
		//KRATOS_THROW_ERROR(std::logic_error,"something is wrong with the mass matrix entry - ZERO!!!","")
		mMdiagInv[i] = 1000000000000.0;
	      
		//KRATOS_WATCH(mMdiagInv[i])
		//KRATOS_THROW_ERROR(std::logic_error,"Zero ELEMENT VOLUMEE!!!!!!!!!!!!!!","")
		//mMdiagInv[i] = 0.0;
	      
	      }
	  }
	
	//KRATOS_WATCH(mMdiagInv)
	//AND NOW WE BUILD THE CONSISTENT MASS MATRIX
	
	for(ModelPart::ElementsContainerType::iterator i = r_model_part.ElementsBegin(); i!=r_model_part.ElementsEnd(); i++)
	  {
	  
	    Geometry< Node<3> >& geom = i->GetGeometry();
	    unsigned int str_nr=0;
	    for (unsigned int k = 0; k<i->GetGeometry().size(); k++)
	      {
		str_nr+=(unsigned int)(i->GetGeometry()[k].FastGetSolutionStepValue(IS_STRUCTURE));
	      }
	      
	    if (geom.size()!=str_nr)
	      {
	      
		GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
		if (Volume<0)
		  Volume*=-1.0;
		//finiding local indices
		//for(int ii = 0; ii<TDim+1; ii++)
		for(unsigned int ii = 0; ii<geom.size(); ii++)
		  {
		    local_indices[ii] = geom[ii].GetDof(DISPLACEMENT_X,dof_position).EquationId();
		  }
	      
		temp = Volume*aaa;
		//element mass matrix has a shape:
		//			2 1 1
		//  A/12.0* 		1 2 1	in 2D
		//			1 1 2
		//
		//			and
		//
		//		    2 1 1 1
		//	V/20.0*     1 2 1 1		in 3D
		//		    1 1 2 1
		//		    1 1 1 2
	      
		//nothing should be added in case of mmassembrane
		for(unsigned int row = 0; row<TDim+1; row++)
		  {
		    unsigned int row_index = local_indices[row] / (TDim); //pressure is a scalar=>matrix size is Tdim times smaller than for vector
		    for(unsigned int col = 0; col<TDim+1; col++)
		      {
			unsigned int col_index = local_indices[col] /(TDim);
			if (row_index==col_index)
			  {
			    //Mconsistent(row_index,col_index) += temp * 2.0;
			    if (TDim==2)
			      Mconsistent(row_index,col_index) += 0.25*temp * 2.0;
			    else if (TDim==3)
			      Mconsistent(row_index,col_index) += 0.2*temp * 2.0;
			  }
			else
			  {
			  
			    //Mconsistent(row_index,col_index) += temp ;
			    if (TDim==2)
			      Mconsistent(row_index,col_index) += 0.25*temp ;
			    else if (TDim==3)
			      Mconsistent(row_index,col_index) += 0.2*temp;
			  
			  }
		      }
		  
		  
		  }
	      }
	  }
      
	//	KRATOS_WATCH("FINISHED BUILDING MASS MATRICES ")
	KRATOS_CATCH("")
	  
	  }
      
      //************************************
      //************************************
      bool CalculateLumpedMassAux()
      {
	  
        KRATOS_TRY

	  ModelPart& model_part=BaseType::GetModelPart();

	bool inv_element=false;
	unsigned int bad_elements_number=0;
	//////////////////////////////////////////////////////////////////////////

	for(ModelPart::ElementsContainerType::iterator im = model_part.ElementsBegin() ;  im != model_part.ElementsEnd() ; ++im)
	  {
	    //get the geometry
	    Geometry< Node<3> >& geom = im->GetGeometry();
	
	    double Area = GeometryUtils::CalculateVolume2D(geom);
	
	    if(Area<0.0) 
	      {
		inv_element=true;

		//WE MARK THE NODES OF THE INVERTED ELEMENTS
		im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE)=1.0;
		im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE)=1.0;
		im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE)=1.0;
		bad_elements_number++;
		     
	      }
	  }
	KRATOS_WATCH("Number of bad elements!")
	  KRATOS_WATCH(bad_elements_number)
	  unsigned int bad_nodes=0;
	for(ModelPart::NodesContainerType::iterator in = model_part.NodesBegin() ;  in != model_part.NodesEnd() ; ++in)
	  {
	    if (in->FastGetSolutionStepValue(IS_INTERFACE)==1.0)
	      bad_nodes++;
	  }
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	while (bad_nodes>=1)
	  {
	
	    for(ModelPart::NodesContainerType::iterator in = model_part.NodesBegin() ;  in != model_part.NodesEnd() ; ++in)
	      {
		//of thats a node of a bad element...
		if (in->FastGetSolutionStepValue(IS_INTERFACE)==1.0)
		  {
		    array_1d<double,3> v = ZeroVector(3);
		    double p=0.0; 
		    unsigned int counter=0;
		
		    for (WeakPointerVector< Node < 3 > >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
			 i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
		      {
			//If the a node is a node of a "good neighbor"
			if (i->FastGetSolutionStepValue(IS_INTERFACE)==0.0)
			  {
			
			    v+=i->FastGetSolutionStepValue(VELOCITY);
			    p+=i->FastGetSolutionStepValue(PRESSURE);

			    counter++;
			  }
		

		      }
		    //if the node has at least one "good" neighbor
		    if (counter>=1)
		      {
			//averaging
			
			v/=counter;
			p/=counter;
			in->FastGetSolutionStepValue(VELOCITY)=v;
			in->FastGetSolutionStepValue(PRESSURE)=p;
			in->FastGetSolutionStepValue(PRESSURE,1)=p;
			in->FastGetSolutionStepValue(IS_INTERFACE)=0.0;					
			bad_nodes+=-1;

					
		      }
		    else{
		      if(in->FastGetSolutionStepValue(IS_STRUCTURE)==1.0) 
			{
			}
		      else			 
			{ 
			  in->Set(TO_ERASE,true);
			  KRATOS_WATCH("Counter is zero.. no clean node " );	
			  KRATOS_WATCH(in->X());
			  KRATOS_WATCH(in->Y());
			  bad_nodes+=-1;
			}
		    }
		    //closing if  
		  }
		//closing for
	      }
	    //closing while
	  }
	KRATOS_WATCH("Finished While Loop, number of bad nodes is")
	  KRATOS_WATCH(bad_nodes)
	  if (bad_nodes==0)
	    inv_element=false;
	return inv_element;
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

