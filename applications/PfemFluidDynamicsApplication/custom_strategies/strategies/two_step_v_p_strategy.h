//
//   Project Name:        KratosPFEMFluidDynamicsApplication $
//   Last modified by:    $Author:                   AFranci $
//   Date:                $Date:                January 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#ifndef KRATOS_TWO_STEP_V_P_STRATEGY_H
#define KRATOS_TWO_STEP_V_P_STRATEGY_H

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "utilities/openmp_utils.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_utilities/modeler_utilities.hpp"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"

#include "custom_utilities/solver_settings.h"

#include "custom_strategies/strategies/gauss_seidel_linear_strategy.h"

#include "pfem_fluid_dynamics_application_variables.h"


#include <stdio.h>      
#include <math.h>     


namespace Kratos {

///@addtogroup PFEMFluidDynamicsApplication
///@{

///@name Kratos Globals
///@{


///@}
///@name Type Definitions
///@{

///@}


///@name  Enum's
///@{


///@}
///@name  Functions
///@{



///@}
///@name Kratos Classes
///@{

template<class TSparseSpace,
class TDenseSpace,
class TLinearSolver
>
class TwoStepVPStrategy : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of TwoStepVPStrategy
    typedef boost::shared_ptr< TwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver> > Pointer;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TDataType TDataType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

    typedef TwoStepVPSolverSettings<TSparseSpace,TDenseSpace,TLinearSolver> SolverSettingsType;

    ///@}
    ///@name Life Cycle
    ///@{


    TwoStepVPStrategy(ModelPart& rModelPart,
               SolverSettingsType& rSolverConfig):
        BaseType(rModelPart)
    {
        InitializeStrategy(rSolverConfig);
    }

    TwoStepVPStrategy(ModelPart& rModelPart,
               /*SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& rSolverConfig,*/
               typename TLinearSolver::Pointer pVelocityLinearSolver,
               typename TLinearSolver::Pointer pPressureLinearSolver,
               bool ReformDofSet = true,
               double VelTol = 0.0001,
               double PresTol = 0.0001,
               int MaxPressureIterations = 1,// Only for predictor-corrector
               unsigned int TimeOrder = 2, 
               unsigned int DomainSize = 2):
        BaseType(rModelPart), // Move Mesh flag, pass as input?
        mVelocityTolerance(VelTol),
        mPressureTolerance(PresTol),
        mMaxPressureIter(MaxPressureIterations),
        mDomainSize(DomainSize),
        mTimeOrder(TimeOrder),
        mReformDofSet(ReformDofSet)
    {
        KRATOS_TRY;

        BaseType::SetEchoLevel(1);

        // Check that input parameters are reasonable and sufficient.
        this->Check();

        bool CalculateNormDxFlag = true;

        bool ReformDofAtEachIteration = false; // DofSet modifiaction is managed by the fractional step strategy, auxiliary strategies should not modify the DofSet directly.

        // Additional Typedefs
        //typedef typename Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3 > > > VarComponent;
        typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
        typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

        //initializing fractional velocity solution step
        typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
        typename SchemeType::Pointer pScheme;
 
	typename SchemeType::Pointer Temp = typename SchemeType::Pointer(new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());
	pScheme.swap(Temp);

        //CONSTRUCTION OF VELOCITY
        BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (pVelocityLinearSolver));

        this->mpMomentumStrategy = typename BaseType::Pointer(new GaussSeidelLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (rModelPart, pScheme, pVelocityLinearSolver, vel_build, ReformDofAtEachIteration, CalculateNormDxFlag));

        this->mpMomentumStrategy->SetEchoLevel( BaseType::GetEchoLevel() );

	vel_build->SetCalculateReactionsFlag(false);
	
        BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pPressureLinearSolver, PRESSURE));

	this->mpPressureStrategy = typename BaseType::Pointer(new GaussSeidelLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (rModelPart, pScheme, pPressureLinearSolver, pressure_build, ReformDofAtEachIteration, CalculateNormDxFlag));

        this->mpPressureStrategy->SetEchoLevel( BaseType::GetEchoLevel() );

	pressure_build->SetCalculateReactionsFlag(false);
	
        KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~TwoStepVPStrategy(){}

    virtual int Check()
    {
        KRATOS_TRY;

        // Check elements and conditions in the model part
        int ierr = BaseType::Check();
        if (ierr != 0) return ierr;

        if(DELTA_TIME.Key() == 0)
            KRATOS_THROW_ERROR(std::runtime_error,"DELTA_TIME Key is 0. Check that the application was correctly registered.","");
        if(BDF_COEFFICIENTS.Key() == 0)
            KRATOS_THROW_ERROR(std::runtime_error,"BDF_COEFFICIENTS Key is 0. Check that the application was correctly registered.","");

        ModelPart& rModelPart = BaseType::GetModelPart();

        if ( mTimeOrder == 2 && rModelPart.GetBufferSize() < 3 )
            KRATOS_THROW_ERROR(std::invalid_argument,"Buffer size too small for fractional step strategy (BDF2), needed 3, got ",rModelPart.GetBufferSize());
        if ( mTimeOrder == 1 && rModelPart.GetBufferSize() < 2 )
            KRATOS_THROW_ERROR(std::invalid_argument,"Buffer size too small for fractional step strategy (Backward Euler), needed 2, got ",rModelPart.GetBufferSize());

        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        for ( ModelPart::ElementIterator itEl = rModelPart.ElementsBegin(); itEl != rModelPart.ElementsEnd(); ++itEl )
        {
            ierr = itEl->Check(rCurrentProcessInfo);
            if (ierr != 0) break;
        }

        for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond)
        {
            ierr = itCond->Check(rCurrentProcessInfo);
            if (ierr != 0) break;
        }

        return ierr;

        KRATOS_CATCH("");
    }

    virtual double Solve()
    {
      // Initialize BDF2 coefficients
      ModelPart& rModelPart = BaseType::GetModelPart();
      this->SetTimeCoefficients(rModelPart.GetProcessInfo());

      if ( BaseType::GetEchoLevel() > 1)
	std::cout << "Solve in two_step_vp strategy "  << std::endl;

      double NormDp = 0.0;

      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
      double currentTime = rCurrentProcessInfo[TIME];
      double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      bool timeIntervalChanged=  rCurrentProcessInfo[TIME_INTERVAL_CHANGED];
 
      unsigned int maxNonLinearIterations=mMaxPressureIter;

      if(timeIntervalChanged==true && currentTime>10*timeInterval ){
	maxNonLinearIterations*=2;
      }
      if(currentTime<10*timeInterval){
	if ( BaseType::GetEchoLevel() > 1)
	  std::cout << "within the first 10 time steps, I consider the given iteration number x3"<< std::endl;
	maxNonLinearIterations*=3;
      }
      if(currentTime<20*timeInterval && currentTime>=10*timeInterval){
	if ( BaseType::GetEchoLevel() > 1)
	  std::cout << "within the second 10 time steps, I consider the given iteration number x2"<< std::endl;
	maxNonLinearIterations*=2;
      }
      bool momentumConverged = true;
      bool continuityConverged = false;
      bool fixedTimeStep=false;
      boost::timer solve_step_time;
      // Iterative solution for pressure
      /* unsigned int timeStep = rCurrentProcessInfo[STEP]; */
      /* if(timeStep==1){ */
      /* 	unsigned int iter=0; */
      /* 	continuityConverged = this->SolveContinuityIteration(iter,maxNonLinearIterations); */
      /* }else if(timeStep==2){ */
      /* 	unsigned int iter=0; */
      /* 	momentumConverged = this->SolveMomentumIteration(iter,maxNonLinearIterations,fixedTimeStep); */
      /* }else{ */

      for(unsigned int it = 0; it < maxNonLinearIterations; ++it)
	{
	  if ( BaseType::GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)
	    std::cout << "----- > iteration: " << it << std::endl;


	  momentumConverged = this->SolveMomentumIteration(it,maxNonLinearIterations,fixedTimeStep);

	  this->CalculateDisplacements();
	  BaseType::MoveMesh();
	  
	  if( fixedTimeStep==false){
	    continuityConverged = this->SolveContinuityIteration(it,maxNonLinearIterations);
	  }
	  if(it==maxNonLinearIterations-1 || ((continuityConverged && momentumConverged) && it>2)){
	    this->UpdateStressStrain();
	  }
	  if ( (continuityConverged && momentumConverged) && it>2)
	    {
	      rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE,false);
	      rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE,false);
	      if ( BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
		std::cout << "V-P strategy converged in " << it+1 << " iterations." << std::endl;
	      break;
	    }

	}

   
      /* } */

      if (!continuityConverged && !momentumConverged && BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
	std::cout << "Convergence tolerance not reached." << std::endl;

      std::cout << "solve_step_time : " << solve_step_time.elapsed() << std::endl;

      if (mReformDofSet)
	this->Clear();

      return NormDp;
    }

    virtual void FinalizeSolutionStep(){
      /* this->UpdateStressStrain(); */
    }

    virtual void InitializeSolutionStep(){
      /* this->ResetActiveToSlivers(); */
      TimeIntervalDurationControl();
    }


    void TimeIntervalDurationControl()
    {


      KRATOS_TRY

      ModelPart& rModelPart = BaseType::GetModelPart();
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

      const double initialTimeInterval = rCurrentProcessInfo[INITIAL_DELTA_TIME];	
      const double currentTimeInterval = rCurrentProcessInfo[CURRENT_DELTA_TIME];
      double updatedTime = rCurrentProcessInfo[TIME];
      double updatedTimeInterval = rCurrentProcessInfo[DELTA_TIME];
      double deltaTimeToNewMilestone=initialTimeInterval;
      double minimumTimeInterval=initialTimeInterval*0.0001;

      rCurrentProcessInfo.SetValue(PREVIOUS_DELTA_TIME,currentTimeInterval);
      rCurrentProcessInfo.SetValue(TIME_INTERVAL_CHANGED,false);
	    
      bool milestoneTimeReached=true;
      bool increaseTimeInterval=true;
      bool timeIntervalReduced=false;

      double tolerance=0.0001;
      updatedTime -= initialTimeInterval;
      unsigned int previousMilestoneStep=updatedTime/initialTimeInterval;
      deltaTimeToNewMilestone=initialTimeInterval*(previousMilestoneStep+1)-updatedTime;

      updatedTimeInterval =currentTimeInterval;
      
      bool badVelocityConvergence=rCurrentProcessInfo[BAD_VELOCITY_CONVERGENCE];
      bool badPressureConvergence=rCurrentProcessInfo[BAD_PRESSURE_CONVERGENCE];
 
      if(updatedTimeInterval<2.0*minimumTimeInterval && BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0){
	std::cout<<"ATTENTION! time step much smaller than initial time step, I'll not reduce it"<<std::endl;
      }
      if(badVelocityConvergence==true && updatedTimeInterval>(2.0*minimumTimeInterval)){
	updatedTimeInterval *=0.5;
	/* std::cout<<"reducing time step (bad convergence at the previous step)"<<updatedTimeInterval<<std::endl; */
	rCurrentProcessInfo.SetValue(TIME_INTERVAL_CHANGED,true);
	timeIntervalReduced=true;
      }

      if(deltaTimeToNewMilestone<(1.0+tolerance)*updatedTimeInterval && deltaTimeToNewMilestone>initialTimeInterval*tolerance){
	rCurrentProcessInfo.SetValue(DELTA_TIME,deltaTimeToNewMilestone);
	if(deltaTimeToNewMilestone<0.75*updatedTimeInterval){
	  timeIntervalReduced=true;
	  rCurrentProcessInfo.SetValue(TIME_INTERVAL_CHANGED,true);
	}
	updatedTimeInterval =deltaTimeToNewMilestone;
	milestoneTimeReached=true;
      }else{
	milestoneTimeReached=false;
	rCurrentProcessInfo.SetValue(DELTA_TIME,updatedTimeInterval);
      }

      if(timeIntervalReduced==false){
	if(updatedTimeInterval>(2.0*minimumTimeInterval)){
	  
	  CheckNodalConditionForTimeStepReduction(updatedTimeInterval,increaseTimeInterval,timeIntervalReduced);

	  if(timeIntervalReduced==false){

	    CheckElementalConditionForTimeStepReduction(increaseTimeInterval);

	  }
	}

	if(increaseTimeInterval==true && initialTimeInterval>(1.0+tolerance)*updatedTimeInterval && badPressureConvergence==false && badVelocityConvergence==false ){
	  IncreaseTimeInterval(updatedTimeInterval,deltaTimeToNewMilestone,tolerance,increaseTimeInterval);
	}
	else{
	  increaseTimeInterval=false;
	}

      }

      double newTimeInterval = rCurrentProcessInfo[DELTA_TIME];
      double milestoneGap=fabs(newTimeInterval-deltaTimeToNewMilestone);
      if(milestoneGap<0.49*newTimeInterval && milestoneTimeReached==false){
	/* std::cout<<"the milestone is very close, I add "<<milestoneGap<<" to "<<newTimeInterval<<std::endl;*/
	newTimeInterval+=milestoneGap;
	rCurrentProcessInfo.SetValue(DELTA_TIME,newTimeInterval);
	milestoneTimeReached=true;
      }

      updatedTime+=newTimeInterval;
      rCurrentProcessInfo.SetValue(TIME,updatedTime);
      rCurrentProcessInfo.SetValue(CURRENT_DELTA_TIME,newTimeInterval);

      /* if(newTimeInterval<(0.49*currentTimeInterval)){ */
      /* std::cout<<"ATTENTION! new time step is more than 2 times smaller than the previous one"<<std::endl; */
      /* } */

      if(increaseTimeInterval==false && milestoneTimeReached==true && fabs(newTimeInterval-initialTimeInterval)>tolerance && !(deltaTimeToNewMilestone>newTimeInterval*(1.0+tolerance))){
	rCurrentProcessInfo.SetValue(CURRENT_DELTA_TIME,currentTimeInterval);
      }


      if (newTimeInterval<initialTimeInterval){
	std::cout<<"current time "<<updatedTime<<" time step: new  "<<newTimeInterval<<" previous "<<currentTimeInterval<<" initial  "<<initialTimeInterval<<"\n"<<std::endl;
      }

    
      KRATOS_CATCH("");


    }


    void CheckNodalConditionForTimeStepReduction(double updatedTimeInterval,
						 bool &increaseTimeInterval,
						 bool &timeIntervalReduced)
    {

      ModelPart& rModelPart = BaseType::GetModelPart();
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

#pragma omp parallel 
	  {
	    ModelPart::NodeIterator NodeBegin;
	    ModelPart::NodeIterator NodeEnd;
	    OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodeBegin,NodeEnd);

	    for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
	      {
		if(itNode->IsNot(TO_ERASE) && itNode->IsNot(ISOLATED) && itNode->IsNot(SOLID)){
		  const array_1d<double,3> &Vel = itNode->FastGetSolutionStepValue(VELOCITY);
		  double NormVelNode=0;
		  for (unsigned int d = 0; d < 3; ++d){
		    NormVelNode+=Vel[d] * Vel[d];
		  }
		  double motionInStep=sqrt(NormVelNode)*updatedTimeInterval;
		  double unsafetyFactor=0;
		  WeakPointerVector< Node < 3 > >& neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
		  for (WeakPointerVector< Node <3> >::iterator nn = neighb_nodes.begin();nn != neighb_nodes.end(); nn++)
		    {
		      array_1d<double,3>  CoorNeighDifference=itNode->Coordinates()-nn->Coordinates();
		      double squaredDistance=0;
		      for (unsigned int d = 0; d < 3; ++d){
			squaredDistance+=CoorNeighDifference[d]*CoorNeighDifference[d];
		      }
		      double nodeDistance=sqrt(squaredDistance);
		      double tempUnsafetyFactor=motionInStep/nodeDistance;
		      if(tempUnsafetyFactor>unsafetyFactor){
			unsafetyFactor=tempUnsafetyFactor;
		      }
		    }

		  if(unsafetyFactor>0.35){
		    increaseTimeInterval=false;
		    if(unsafetyFactor>1.0){
		      double temporaryTimeInterval = rCurrentProcessInfo[DELTA_TIME];
		      double reducedTimeInterval=0.5*updatedTimeInterval;
		      if(reducedTimeInterval<temporaryTimeInterval){
			rCurrentProcessInfo.SetValue(DELTA_TIME,reducedTimeInterval);
			/* std::cout<<"reducing time step (nodal criterion)"<<reducedTimeInterval<<std::endl; */
			rCurrentProcessInfo.SetValue(TIME_INTERVAL_CHANGED,true);
			timeIntervalReduced=true;
			break;
		      }
		    }
		  }
		}
	      }
	  }

    }


    void CheckElementalConditionForTimeStepReduction(bool &increaseTimeInterval)
    {

      ModelPart& rModelPart = BaseType::GetModelPart();
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

#pragma omp parallel
      {
	ModelPart::ElementIterator ElemBegin;
	ModelPart::ElementIterator ElemEnd;
	OpenMPUtils::PartitionedIterators(rModelPart.Elements(),ElemBegin,ElemEnd);
	for ( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem )
	  {
	    double temporaryTimeInterval=rCurrentProcessInfo[DELTA_TIME];
	    double currentElementalArea =  0;
	    const unsigned int dimension = (itElem)->GetGeometry().WorkingSpaceDimension();
	    if(dimension==2){
	      currentElementalArea =  (itElem)->GetGeometry().Area();
	      Geometry<Node<3> >  updatedElementCoordinates;
	      bool solidElement=false;
	      for(unsigned int i=0; i<itElem->GetGeometry().size(); i++)
		{
		  if(itElem->GetGeometry()[i].Is(SOLID) || itElem->GetGeometry()[i].Is(TO_ERASE) || itElem->IsNot(ACTIVE)){
		    solidElement=true;
		  }
	
		  const array_1d<double,3> &Vel = itElem->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
		  Point<3> updatedNodalCoordinates=itElem->GetGeometry()[i].Coordinates()+Vel*temporaryTimeInterval;
		  updatedElementCoordinates.push_back(Node<3>::Pointer(new Node<3>(i,updatedNodalCoordinates.X(),updatedNodalCoordinates.Y(),updatedNodalCoordinates.Z())));
		}

	      double newArea=0;
	      if(itElem->GetGeometry().size()==3){
		Triangle2D3<Node<3> > myGeometry(updatedElementCoordinates);
		newArea=myGeometry.Area();
	      }else if(itElem->GetGeometry().size()==6){
		Triangle2D6<Node<3> > myGeometry(updatedElementCoordinates);
		newArea=myGeometry.Area();
	      }else{
		std::cout<<"GEOMETRY NOT DEFINED"<<std::endl;
	      }

	      if(solidElement==true){
		newArea=currentElementalArea;
	      }

	      if(newArea<0.001*currentElementalArea && currentElementalArea>0){
		double reducedTimeInterval=0.5*temporaryTimeInterval;
	      
		if(reducedTimeInterval<temporaryTimeInterval){
		  rCurrentProcessInfo.SetValue(DELTA_TIME,reducedTimeInterval);
		  /* std::cout<<"reducing time step (elemental inversion)"<<reducedTimeInterval<<std::endl; */
		  rCurrentProcessInfo.SetValue(TIME_INTERVAL_CHANGED,true);
		  increaseTimeInterval=false;
		  break;
		}
	      }else{
		Geometry<Node<3> >  updatedEnlargedElementCoordinates;

		for(unsigned int i=0; i<itElem->GetGeometry().size(); i++)
		  {
		    const array_1d<double,3> &Vel = itElem->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
		    Point<3> updatedNodalCoordinates=itElem->GetGeometry()[i].Coordinates()+Vel*temporaryTimeInterval*2.5;
		    updatedEnlargedElementCoordinates.push_back(Node<3>::Pointer(new Node<3>(i,updatedNodalCoordinates.X(),updatedNodalCoordinates.Y(),updatedNodalCoordinates.Z())));

		  }

		if(itElem->GetGeometry().size()==3){
		  Triangle2D3<Node<3> > myGeometry(updatedEnlargedElementCoordinates);
		  newArea=myGeometry.Area();
		}else if(itElem->GetGeometry().size()==6){
		  Triangle2D6<Node<3> > myGeometry(updatedEnlargedElementCoordinates);
		  newArea=myGeometry.Area();
		}else{
		  std::cout<<"GEOMETRY NOT DEFINED"<<std::endl;
		}

		if(newArea<0.001*currentElementalArea && currentElementalArea>0){
		  increaseTimeInterval=false;
		  /* std::cout<<"I'll not reduce the time step but I'll not allow to increase it"<<std::endl; */
		}

	      }
	    }
	    else if(dimension==3){
	      double currentElementalVolume =  (itElem)->GetGeometry().Volume();
	      Geometry<Node<3> >  updatedElementCoordinates;
	      bool solidElement=false;
	      for(unsigned int i=0; i<itElem->GetGeometry().size(); i++)
		{
		  if(itElem->GetGeometry()[i].Is(SOLID)){
		    solidElement=true;
		  }
		  const array_1d<double,3> &Vel = itElem->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
		  Point<3> updatedNodalCoordinates=itElem->GetGeometry()[i].Coordinates()+Vel*temporaryTimeInterval;
		  updatedElementCoordinates.push_back(Node<3>::Pointer(new Node<3>(i,updatedNodalCoordinates.X(),updatedNodalCoordinates.Y(),updatedNodalCoordinates.Z())));
		}

	      double newVolume=0;
	      if(itElem->GetGeometry().size()==4){
		Tetrahedra3D4<Node<3> > myGeometry(updatedElementCoordinates);
		newVolume=myGeometry.Volume();
	      }else if(itElem->GetGeometry().size()==10){
		Tetrahedra3D10<Node<3> > myGeometry(updatedElementCoordinates);
		newVolume=myGeometry.Volume();
	      }else{
		std::cout<<"GEOMETRY NOT DEFINED"<<std::endl;
	      }

	      if(solidElement==true){
		newVolume=currentElementalVolume;
	      }

	      if(newVolume<0.001*currentElementalVolume && currentElementalVolume>0){
		double reducedTimeInterval=0.5*temporaryTimeInterval;
	      
		if(reducedTimeInterval<temporaryTimeInterval){
		  rCurrentProcessInfo.SetValue(DELTA_TIME,reducedTimeInterval);
		  /* std::cout<<"reducing time step (elemental inversion)"<<reducedTimeInterval<<std::endl; */
		  rCurrentProcessInfo.SetValue(TIME_INTERVAL_CHANGED,true);
		  increaseTimeInterval=false;
		  break;
		}
	      }else{
		Geometry<Node<3> >  updatedEnlargedElementCoordinates;

		for(unsigned int i=0; i<itElem->GetGeometry().size(); i++)
		  {
		    const array_1d<double,3> &Vel = itElem->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
		    Point<3> updatedNodalCoordinates=itElem->GetGeometry()[i].Coordinates()+Vel*temporaryTimeInterval*2.5;
		    updatedEnlargedElementCoordinates.push_back(Node<3>::Pointer(new Node<3>(i,updatedNodalCoordinates.X(),updatedNodalCoordinates.Y(),updatedNodalCoordinates.Z())));
		  }

		if(itElem->GetGeometry().size()==4){
		  Tetrahedra3D4<Node<3> > myGeometry(updatedEnlargedElementCoordinates);
		  newVolume=myGeometry.Volume();
		}else if(itElem->GetGeometry().size()==10){
		  Tetrahedra3D10<Node<3> > myGeometry(updatedEnlargedElementCoordinates);
		  newVolume=myGeometry.Volume();
		}else{
		  std::cout<<"GEOMETRY NOT DEFINED"<<std::endl;
		}

		if(newVolume<0.001*currentElementalVolume && currentElementalVolume>0){
		  increaseTimeInterval=false;
		  /* std::cout<<"I'll not reduce the time step but I'll not allow to increase it"<<std::endl; */
		}


	      }



	    }
		
	  }

      }
    }

    void IncreaseTimeInterval(double updatedTimeInterval,
			      double deltaTimeToNewMilestone,
			      double tolerance,
			      bool &increaseTimeInterval)
    {
      ModelPart& rModelPart = BaseType::GetModelPart();
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();    
      double increasedTimeInterval=2.0*updatedTimeInterval;
      if(increasedTimeInterval<deltaTimeToNewMilestone*(1.0+tolerance)){
	rCurrentProcessInfo.SetValue(DELTA_TIME,increasedTimeInterval);
	/* std::cout<<"increasing time step "<<increasedTimeInterval<<" previous one="<<updatedTimeInterval<<std::endl; */
	rCurrentProcessInfo.SetValue(TIME_INTERVAL_CHANGED,true);
      }else{
	increaseTimeInterval=false;
      }
    }


    void CalculatePressureVelocity()
    {
      ModelPart& rModelPart = BaseType::GetModelPart();
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
      const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      unsigned int timeStep = rCurrentProcessInfo[STEP];

      for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
	   i != rModelPart.NodesEnd(); ++i)
        {
	  if(timeStep==1){
	    (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0)=0;
	    (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1)=0;
	  }else{
	    double  & CurrentPressure      = (i)->FastGetSolutionStepValue(PRESSURE, 0);
	    double  & PreviousPressure     = (i)->FastGetSolutionStepValue(PRESSURE, 1);
	    double  & CurrentPressureVelocity  = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0);
	    CurrentPressureVelocity = (CurrentPressure-PreviousPressure)/timeInterval;
	  }
       
        }
    }

    void CalculatePressureAcceleration()
    {
      ModelPart& rModelPart = BaseType::GetModelPart();
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
      const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      unsigned int timeStep = rCurrentProcessInfo[STEP];

      for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
	   i != rModelPart.NodesEnd(); ++i)
        {
	  if(timeStep==1){
	    (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0)=0;
	    (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 1)=0;
	  }else{
	    double & CurrentPressureVelocity      = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 0);
	    double & PreviousPressureVelocity     = (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY, 1);
	    double & CurrentPressureAcceleration  = (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION, 0);
	    CurrentPressureAcceleration = (CurrentPressureVelocity-PreviousPressureVelocity)/timeInterval;

	  }
        }
    }

    void CalculateAccelerations()
    {
      ModelPart& rModelPart = BaseType::GetModelPart();
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
      Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
      
      for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
	   i != rModelPart.NodesEnd(); ++i)
        {

	  array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY, 0);
	  array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);

	  array_1d<double, 3 > & CurrentAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
	  array_1d<double, 3 > & PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

	  if((i)->IsNot(ISOLATED) || (i)->Is(SOLID)){

	    UpdateAccelerations (CurrentAcceleration, CurrentVelocity, PreviousAcceleration, PreviousVelocity,BDFcoeffs);


	  }else {
	    (i)->FastGetSolutionStepValue(PRESSURE,0) = 0.0; 
	    (i)->FastGetSolutionStepValue(PRESSURE,1) = 0.0; 
	    (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY,0) = 0.0; 
	    (i)->FastGetSolutionStepValue(PRESSURE_VELOCITY,1) = 0.0; 
	    (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION,0) = 0.0; 
	    (i)->FastGetSolutionStepValue(PRESSURE_ACCELERATION,1) = 0.0; 
	    if((i)->SolutionStepsDataHas(VOLUME_ACCELERATION)){
	      array_1d<double, 3 >& VolumeAcceleration = (i)->FastGetSolutionStepValue(VOLUME_ACCELERATION);
	      (i)->FastGetSolutionStepValue(ACCELERATION,0) = VolumeAcceleration;
	      (i)->FastGetSolutionStepValue(VELOCITY,0) += VolumeAcceleration*rCurrentProcessInfo[DELTA_TIME];
	    }
	  }


        }
    }

    inline void UpdateAccelerations(array_1d<double, 3 > & CurrentAcceleration,
				    const array_1d<double, 3 > & CurrentVelocity,
				    array_1d<double, 3 > & PreviousAcceleration,
				    const array_1d<double, 3 > & PreviousVelocity,
				    Vector& BDFcoeffs)
    {
      /* noalias(PreviousAcceleration)=CurrentAcceleration; */
      noalias(CurrentAcceleration) = -BDFcoeffs[1]*(CurrentVelocity-PreviousVelocity) - PreviousAcceleration ;
      // std::cout<<"rBDFCoeffs[0] is "<<rBDFCoeffs[0]<<std::endl;//3/(2*delta_t)
      // std::cout<<"rBDFCoeffs[1] is "<<rBDFCoeffs[1]<<std::endl;//-2/(delta_t)
      // std::cout<<"rBDFCoeffs[2] is "<<rBDFCoeffs[2]<<std::endl;//1/(2*delta_t)
    }

    void CalculateDisplacements()
    {
      ModelPart& rModelPart = BaseType::GetModelPart();
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
      const double TimeStep = rCurrentProcessInfo[DELTA_TIME];
      
      for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
	   i != rModelPart.NodesEnd(); ++i)
        {

	  array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY, 0);
	  array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);

	  array_1d<double, 3 > & CurrentDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT, 0);
	  array_1d<double, 3 > & PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
	  
	  if( i->IsFixed(DISPLACEMENT_X) == false )
	    CurrentDisplacement[0] = 0.5* TimeStep *(CurrentVelocity[0]+PreviousVelocity[0]) + PreviousDisplacement[0];	  

	  if( i->IsFixed(DISPLACEMENT_Y) == false )
	    CurrentDisplacement[1] = 0.5* TimeStep *(CurrentVelocity[1]+PreviousVelocity[1]) + PreviousDisplacement[1];

	  if( i->IsFixed(DISPLACEMENT_Z) == false )
	    CurrentDisplacement[2] = 0.5* TimeStep *(CurrentVelocity[2]+PreviousVelocity[2]) + PreviousDisplacement[2];

        }
    }

    void ResetActiveToSlivers()
    {

#pragma omp parallel
      {
	ModelPart& rModelPart = BaseType::GetModelPart();
	ModelPart::ElementIterator ElemBegin;
	ModelPart::ElementIterator ElemEnd;
	ModelerUtilities ModelerUtils;
	OpenMPUtils::PartitionedIterators(rModelPart.Elements(),ElemBegin,ElemEnd);
	for ( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem )
	  {
	    double ElementalVolume =  0;
	    const unsigned int dimension = (itElem)->GetGeometry().WorkingSpaceDimension();
	    if(dimension==2){
	      ElementalVolume =  (itElem)->GetGeometry().Area();
	    }else if(dimension==3){
	      ElementalVolume =  (itElem)->GetGeometry().Volume();
	    }else{
	      ElementalVolume = 0;
	    }
	    double ModelPartVolume=ModelerUtils.ComputeModelPartVolume(rModelPart);
	    double CriticalVolume=0.01*ModelPartVolume/double(rModelPart.Elements().size());

	    if(ElementalVolume<CriticalVolume && ElementalVolume>0){
	      (itElem)->Set(ACTIVE,false);
	      std::cout<<"RESET ACTIVE FOR THIS SLIVER! \t";
	      std::cout<<"its volume is "<<ElementalVolume<<" vs CriticalVolume "<<CriticalVolume<<std::endl;
	    }else{
	      (itElem)->Set(ACTIVE,true);
	    }
	   

	  }

      }
        
    }


   void UpdateStressStrain()
   {
     ModelPart& rModelPart = BaseType::GetModelPart();
     ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

#pragma omp parallel
     {
       ModelPart::ElementIterator ElemBegin;
       ModelPart::ElementIterator ElemEnd;
       OpenMPUtils::PartitionedIterators(rModelPart.Elements(),ElemBegin,ElemEnd);

       for ( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem )
	 {
	   /* itElem-> InitializeElementStrainStressState(); */
	   itElem-> InitializeSolutionStep(rCurrentProcessInfo);
	 }

     }

     this->CalculateAccelerations(); 
     this->CalculatePressureVelocity();
     this->CalculatePressureAcceleration();

   }

    virtual void Clear()
    {
        mpMomentumStrategy->Clear();
        mpPressureStrategy->Clear();
    }


    ///@}
    ///@name Access
    ///@{

    virtual void SetEchoLevel(int Level)
    {
        BaseType::SetEchoLevel(Level);
        int StrategyLevel = Level > 0 ? Level - 1 : 0;
        mpMomentumStrategy->SetEchoLevel(StrategyLevel);
        mpPressureStrategy->SetEchoLevel(StrategyLevel);
    }


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "TwoStepVPStrategy" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "TwoStepVPStrategy";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:

    ///@name Protected Life Cycle
    ///@{


    ///@}
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    /// Calculate the coefficients for time iteration.
    /**
     * @param rCurrentProcessInfo ProcessInfo instance from the fluid ModelPart. Must contain DELTA_TIME and BDF_COEFFICIENTS variables.
     */
    void SetTimeCoefficients(ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if (mTimeOrder == 2)
        {
            //calculate the BDF coefficients
            double Dt = rCurrentProcessInfo[DELTA_TIME];
            double OldDt = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

            double Rho = OldDt / Dt;
            double TimeCoeff = 1.0 / (Dt * Rho * Rho + Dt * Rho);

            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            BDFcoeffs.resize(3, false);

            BDFcoeffs[0] = TimeCoeff * (Rho * Rho + 2.0 * Rho); //coefficient for step n+1 (3/2Dt if Dt is constant)
            BDFcoeffs[1] = -TimeCoeff * (Rho * Rho + 2.0 * Rho + 1.0); //coefficient for step n (-4/2Dt if Dt is constant)
            BDFcoeffs[2] = TimeCoeff; //coefficient for step n-1 (1/2Dt if Dt is constant)
        }
        else if (mTimeOrder == 1)
        {
            double Dt = rCurrentProcessInfo[DELTA_TIME];
            double TimeCoeff = 1.0 / Dt;

            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            BDFcoeffs.resize(2, false);

            BDFcoeffs[0] = TimeCoeff; //coefficient for step n+1 (1/Dt)
            BDFcoeffs[1] = -TimeCoeff; //coefficient for step n (-1/Dt)
        }

        KRATOS_CATCH("");
    }

    bool SolveMomentumIteration(unsigned int it,unsigned int maxIt, bool & fixedTimeStep)
    {
      ModelPart& rModelPart = BaseType::GetModelPart();
      int Rank = rModelPart.GetCommunicator().MyPID();
      bool ConvergedMomentum = false;
      double NormDv = 0;
      fixedTimeStep=false;
      // build momentum system and solve for fractional step velocity increment
      rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);

      std::cout<<"-------- m o m e n t u m   e q u a t i o n s ----------"<<std::endl;
      if(it==0){
	mpMomentumStrategy->InitializeSolutionStep();
      }
      /* else{ */
      /* 	NormDv = mpMomentumStrategy->Solve(); */
      /* } */
      NormDv = mpMomentumStrategy->Solve(); 
	  
      if (BaseType::GetEchoLevel() > 1 && Rank == 0)
	std::cout<<"-------------- s o l v e d ! ------------------"<<std::endl;

      double DvErrorNorm = 0; 
      ConvergedMomentum = this->CheckVelocityConvergence(NormDv,DvErrorNorm);
      // Check convergence
      if(it==maxIt-1){
	fixedTimeStep=this->FixTimeStepMomentum(DvErrorNorm);
      }

      if (!ConvergedMomentum && BaseType::GetEchoLevel() > 0 && Rank == 0)
	std::cout << "Momentum equations did not reach the convergence tolerance." << std::endl;

      return ConvergedMomentum;
    }


    bool SolveContinuityIteration(unsigned int it,unsigned int maxIt)
    {
      ModelPart& rModelPart = BaseType::GetModelPart();
      int Rank = rModelPart.GetCommunicator().MyPID();
      bool ConvergedContinuity = false;
      double NormDp = 0;

      // 2. Pressure solution 
      rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,5);

      std::cout<<"          -------- c o n t i n u i t y   e q u a t i o n ----------"<<std::endl;
 
      if(it==0){
	mpPressureStrategy->InitializeSolutionStep();
      }
      /* else{ */
      /* 	NormDp = mpPressureStrategy->Solve(); */
      /* } */
      NormDp = mpPressureStrategy->Solve();

      if (BaseType::GetEchoLevel() > 0 && Rank == 0)
	std::cout << "The norm of pressure is: " << NormDp << std::endl;

      double DpErrorNorm = 0; 
      ConvergedContinuity = this->CheckPressureConvergence(NormDp,DpErrorNorm);

      // Check convergence
      if(it==maxIt-1){
      	ConvergedContinuity=this->FixTimeStepContinuity(DpErrorNorm);
      }

      if (!ConvergedContinuity && BaseType::GetEchoLevel() > 0 && Rank == 0)
	std::cout << "Continuity equation did not reach the convergence tolerance." << std::endl;

      return ConvergedContinuity;
    }


    bool CheckVelocityConvergence(const double NormDv, double& errorNormDv)
    {
        ModelPart& rModelPart = BaseType::GetModelPart();

        double NormV = 0.00;
        errorNormDv = 0;

#pragma omp parallel reduction(+:NormV)
        {
            ModelPart::NodeIterator NodeBegin;
            ModelPart::NodeIterator NodeEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodeBegin,NodeEnd);

            for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
            {
                const array_1d<double,3> &Vel = itNode->FastGetSolutionStepValue(VELOCITY);

		double NormVelNode=0;

                for (unsigned int d = 0; d < 3; ++d){
		  NormVelNode+=Vel[d] * Vel[d];
		  NormV += Vel[d] * Vel[d];
		}

            }
        }

        BaseType::GetModelPart().GetCommunicator().SumAll(NormV);

        NormV = sqrt(NormV);

        if (NormV == 0.0) NormV = 1.00;

	errorNormDv = NormDv / NormV;
	
        if ( BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0){
	  std::cout << "The norm of velocity increment is: " << NormDv << std::endl;
	  std::cout << "The norm of velocity is: " << NormV << std::endl;
	  std::cout << "Velocity error: " << errorNormDv << "mVelocityTolerance: " << mVelocityTolerance<< std::endl;
	}else{
	  std::cout<<"Velocity error: "<< errorNormDv <<" velTol: " << mVelocityTolerance<< std::endl;
	}
	
        if (errorNormDv < mVelocityTolerance)
        {
            return true;
        }
        else{
            return false;
	}
    }



    bool CheckPressureConvergence(const double NormDp, double& errorNormDp)
    {
        ModelPart& rModelPart = BaseType::GetModelPart();

        double NormP = 0.00;
        errorNormDp = 0;

#pragma omp parallel reduction(+:NormP)
        {
            ModelPart::NodeIterator NodeBegin;
            ModelPart::NodeIterator NodeEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodeBegin,NodeEnd);

            for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
            {
                const double Pr = itNode->FastGetSolutionStepValue(PRESSURE);
                NormP += Pr * Pr;
            }
        }

        BaseType::GetModelPart().GetCommunicator().SumAll(NormP);

        NormP = sqrt(NormP);

        if (NormP == 0.0) NormP = 1.00;

        errorNormDp = NormDp / NormP;

        if ( BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0){
	  std::cout << "         The norm of pressure increment is: " << NormDp << std::endl;
	  std::cout << "         The norm of pressure is: " << NormP << std::endl;
            std::cout << "         Pressure error: " <<errorNormDp  << std::endl;
	}else{
            std::cout<<"         Pressure error: "<<errorNormDp <<" presTol: "<<mPressureTolerance << std::endl;
	}

        if ( errorNormDp< mPressureTolerance)
        {
            return true;
        }
        else
            return false;
    }

    bool FixTimeStepMomentum(const double DvErrorNorm)
    {
      ModelPart& rModelPart = BaseType::GetModelPart();
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
      double currentTime = rCurrentProcessInfo[TIME];
      double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      double minTolerance=0.005;
      bool fixedTimeStep=false;
      if(currentTime<10*timeInterval){
	minTolerance=10;
      }

      bool isItNan=false;
      isItNan=std::isnan(DvErrorNorm);
      bool isItInf=false;
      isItInf=std::isinf(DvErrorNorm);
      if((DvErrorNorm>minTolerance || (DvErrorNorm<0 && DvErrorNorm>0) || (DvErrorNorm!=DvErrorNorm) || isItNan==true || isItInf==true) && DvErrorNorm!=0 && DvErrorNorm!=1){
	rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE,true);
	std::cout << "NOT GOOD CONVERGENCE!!! I'll reduce the next time interval"<<DvErrorNorm<< std::endl;
	minTolerance=0.05;
	if(DvErrorNorm>minTolerance){
	  std::cout<< "BAD CONVERGENCE!!! I GO AHEAD WITH THE PREVIOUS VELOCITY AND PRESSURE FIELDS"<<DvErrorNorm<< std::endl;
	  fixedTimeStep=true;
#pragma omp parallel 
	  {
	    ModelPart::NodeIterator NodeBegin;
	    ModelPart::NodeIterator NodeEnd;
	    OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodeBegin,NodeEnd);

	    for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
	      {
		itNode->FastGetSolutionStepValue(VELOCITY,0)=itNode->FastGetSolutionStepValue(VELOCITY,1);
		itNode->FastGetSolutionStepValue(PRESSURE,0)=itNode->FastGetSolutionStepValue(PRESSURE,1);
		itNode->FastGetSolutionStepValue(ACCELERATION,0)=itNode->FastGetSolutionStepValue(ACCELERATION,1);
	      }
	  }
	}

      }else{
	rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE,false);
      }
      return fixedTimeStep;
    }

   bool FixTimeStepContinuity(const double DvErrorNorm)
    {
      ModelPart& rModelPart = BaseType::GetModelPart();
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
      double currentTime = rCurrentProcessInfo[TIME];
      double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      double minTolerance=0.005;
      bool fixedTimeStep=false;
      if(currentTime<10*timeInterval){
	minTolerance=10;
      }

      bool isItNan=false;
      isItNan=std::isnan(DvErrorNorm);
      bool isItInf=false;
      isItInf=std::isinf(DvErrorNorm);
      if((DvErrorNorm>minTolerance || (DvErrorNorm<0 && DvErrorNorm>0) || (DvErrorNorm!=DvErrorNorm) || isItNan==true || isItInf==true) && DvErrorNorm!=0 && DvErrorNorm!=1){
	fixedTimeStep=true;
	rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE,true);
      }else{
	rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE,false);
      }
      return fixedTimeStep;
    }




    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    double mVelocityTolerance;

    double mPressureTolerance;

    unsigned int mMaxPressureIter;

    unsigned int mDomainSize;

    unsigned int mTimeOrder;

    bool mReformDofSet;

    // Fractional step index.
    /*  1 : Momentum step (calculate fractional step velocity)
      * 2-3 : Unused (reserved for componentwise calculation of frac step velocity)
      * 4 : Pressure step
      * 5 : Computation of projections
      * 6 : End of step velocity
      */
//    unsigned int mStepId;

    /// Scheme for the solution of the momentum equation
    StrategyPointerType mpMomentumStrategy;

    /// Scheme for the solution of the mass equation
    StrategyPointerType mpPressureStrategy;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    void InitializeStrategy(SolverSettingsType& rSolverConfig)
    {
        KRATOS_TRY;
        
        mTimeOrder = rSolverConfig.GetTimeOrder();
        
        // Check that input parameters are reasonable and sufficient.
        this->Check();

        //ModelPart& rModelPart = this->GetModelPart();

        mDomainSize = rSolverConfig.GetDomainSize();

        mReformDofSet = rSolverConfig.GetReformDofSet();

        BaseType::SetEchoLevel(rSolverConfig.GetEchoLevel());

        // Initialize strategies for each step
        bool HaveVelStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Velocity,mpMomentumStrategy);

        if (HaveVelStrategy)
        {
            rSolverConfig.FindTolerance(SolverSettingsType::Velocity,mVelocityTolerance);
            /* rSolverConfig.FindMaxIter(SolverSettingsType::Velocity,mMaxVelocityIter); */
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,"TwoStepVPStrategy error: No Velocity strategy defined in FractionalStepSettings","");
        }

        bool HavePressStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Pressure,mpPressureStrategy);

        if (HavePressStrategy)
        {
            rSolverConfig.FindTolerance(SolverSettingsType::Pressure,mPressureTolerance);
            rSolverConfig.FindMaxIter(SolverSettingsType::Pressure,mMaxPressureIter);
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,"TwoStepVPStrategy error: No Pressure strategy defined in FractionalStepSettings","");
        }

        // Check input parameters
        this->Check();

        KRATOS_CATCH("");
    }


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TwoStepVPStrategy& operator=(TwoStepVPStrategy const& rOther){}

    /// Copy constructor.
    TwoStepVPStrategy(TwoStepVPStrategy const& rOther){}


    ///@}

}; /// Class TwoStepVPStrategy

///@}
///@name Type Definitions
///@{


///@}

///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_TWO_STEP_V_P_STRATEGY_H
