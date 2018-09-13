//
//   Project Name:        KratosPFEMFluidDynamicsApplication $
//   Last modified by:    $Author:                   AFranci $
//   Date:                $Date:                   June 2018 $
//   Revision:            $Revision:                     0.0 $
//
//

#ifndef KRATOS_NODAL_TWO_STEP_V_P_STRATEGY_H
#define KRATOS_NODAL_TWO_STEP_V_P_STRATEGY_H

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"
#include "utilities/openmp_utils.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/builders_and_solvers/nodal_residualbased_elimination_builder_and_solver.h"
#include "custom_strategies/builders_and_solvers/nodal_residualbased_elimination_builder_and_solver_continuity.h"
#include "custom_strategies/builders_and_solvers/nodal_residualbased_block_builder_and_solver.h"

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
    class NodalTwoStepVPStrategy : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
    {
    public:
      ///@name Type Definitions
      ///@{
      KRATOS_CLASS_POINTER_DEFINITION(NodalTwoStepVPStrategy);

      /// Counted pointer of NodalTwoStepVPStrategy
      //typedef boost::shared_ptr< NodalTwoStepVPStrategy<TSparseSpace, TDenseSpace, TLinearSolver> > Pointer;

      typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

      typedef typename BaseType::TDataType TDataType;
    
      /// Node type (default is: Node<3>)
      typedef Node <3> NodeType;
    
      /// Geometry type (using with given NodeType)
      typedef Geometry<NodeType> GeometryType;

      typedef std::size_t SizeType;

      //typedef typename BaseType::DofSetType DofSetType;

      typedef typename BaseType::DofsArrayType DofsArrayType;

      typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

      typedef typename BaseType::TSystemVectorType TSystemVectorType;

      typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

      typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

      typedef typename BaseType::ElementsArrayType ElementsArrayType;

      typedef typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

      typedef TwoStepVPSolverSettings<TSparseSpace,TDenseSpace,TLinearSolver> SolverSettingsType;

      typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;


      ///@}
      ///@name Life Cycle
      ///@{


    NodalTwoStepVPStrategy(ModelPart& rModelPart,
		      SolverSettingsType& rSolverConfig):
      BaseType(rModelPart)
      {
        InitializeStrategy(rSolverConfig);
      }

    NodalTwoStepVPStrategy(ModelPart& rModelPart,
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
	  /* typename SchemeType::Pointer Temp = typename SchemeType::Pointer(new IncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ()); */
	  pScheme.swap(Temp);

	  //CONSTRUCTION OF VELOCITY
	  BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new NodalResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (pVelocityLinearSolver));
	  /* BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (pVelocityLinearSolver)); */

	  this->mpMomentumStrategy = typename BaseType::Pointer(new GaussSeidelLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (rModelPart, pScheme, pVelocityLinearSolver, vel_build, ReformDofAtEachIteration, CalculateNormDxFlag));

	  this->mpMomentumStrategy->SetEchoLevel( BaseType::GetEchoLevel() );

	  vel_build->SetCalculateReactionsFlag(false);
	
	  /* BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pPressureLinearSolver, PRESSURE)); */
	  /* BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver >(pPressureLinearSolver)); */
	  BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new NodalResidualBasedEliminationBuilderAndSolverContinuity<TSparseSpace, TDenseSpace, TLinearSolver > (pPressureLinearSolver));
	  /* BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new NodalResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver >(pPressureLinearSolver)); */

	  this->mpPressureStrategy = typename BaseType::Pointer(new GaussSeidelLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (rModelPart, pScheme, pPressureLinearSolver, pressure_build, ReformDofAtEachIteration, CalculateNormDxFlag));

	  this->mpPressureStrategy->SetEchoLevel( BaseType::GetEchoLevel() );

	  pressure_build->SetCalculateReactionsFlag(false);
	
	  KRATOS_CATCH("");
	}

      /// Destructor.
      virtual ~NodalTwoStepVPStrategy(){}

      int Check() override
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

      double Solve() override
      {
	// Initialize BDF2 coefficients
	ModelPart& rModelPart = BaseType::GetModelPart();
	this->SetTimeCoefficients(rModelPart.GetProcessInfo());
	double NormDp = 0.0;
	ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
	double currentTime = rCurrentProcessInfo[TIME];
	double timeInterval = rCurrentProcessInfo[DELTA_TIME];
	bool timeIntervalChanged=  rCurrentProcessInfo[TIME_INTERVAL_CHANGED];
 
	unsigned int maxNonLinearIterations=mMaxPressureIter;
	/* if ( BaseType::GetEchoLevel() > 1) */
	/* 	std::cout << "Solve with two_step_vp strategy "  << std::endl; */


	std::cout << "\n                   Solve with nodally_integrated_two_step_vp strategy at t="<< currentTime<<"s"<<std::endl;


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
	/* boost::timer solve_step_time; */
      
	this->InitializeSolutionStep();
	for(unsigned int it = 0; it < maxNonLinearIterations; ++it)
	  {
	  
	    if ( BaseType::GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)
	      std::cout << "----- > iteration: " << it << std::endl;

	    if(it==0){
	      this->ComputeNodalVolume();
	      this->InitializeNonLinearIterations();

	      this->CalcNodalStrainsAndStresses();

	    }else{
	      this->CalcNodalStrainsAndStresses();
	    }

	    momentumConverged = this->SolveMomentumIteration(it,maxNonLinearIterations,fixedTimeStep);

	    this->UpdateTopology(rModelPart, BaseType::GetEchoLevel());

	    this->ComputeNodalVolume();
	    this->InitializeNonLinearIterations();
	    this->CalcNodalStrains();
	    	     	      
	    if( fixedTimeStep==false){
	      continuityConverged = this->SolveContinuityIteration(it,maxNonLinearIterations);
	    }
		
	    if(it==maxNonLinearIterations-1 || ((continuityConverged && momentumConverged) && it>2)){
	      this->CalculateAccelerations();
	    }
	    if ( (continuityConverged && momentumConverged) && it>2)
	      {
		rCurrentProcessInfo.SetValue(BAD_VELOCITY_CONVERGENCE,false);
		rCurrentProcessInfo.SetValue(BAD_PRESSURE_CONVERGENCE,false);
		/* if ( BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0) */
		std::cout << "nodal V-P strategy converged in " << it+1 << " iterations." << std::endl;
		break;
	      }

	  }

   
	/* } */

	if (!continuityConverged && !momentumConverged && BaseType::GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
	  std::cout << "Convergence tolerance not reached." << std::endl;

	/* std::cout << "solve_step_time : " << solve_step_time.elapsed() << std::endl; */

	if (mReformDofSet)
	  this->Clear();

	return NormDp;
      }

      void FinalizeSolutionStep() override
      {
	/* this->UpdateStressStrain(); */
      }

      void Initialize() override
      {
	std::cout << "                 Initialize!!! " << std::endl;

	ModelPart& rModelPart = BaseType::GetModelPart();
	const unsigned int dimension =  rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
	unsigned int sizeStrains=3*(dimension-1);


#pragma omp parallel
	{
	  ModelPart::NodeIterator NodesBegin;
	  ModelPart::NodeIterator NodesEnd;
	  OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

	  for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
	    {
	      WeakPointerVector< Node < 3 > >& neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
	      unsigned int neighbourNodes=neighb_nodes.size();
	      unsigned int sizeSDFNeigh=neighbourNodes*dimension;
	      
	      if(itNode->SolutionStepsDataHas(NODAL_CAUCHY_STRESS)){
		Vector& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS);
		if(rNodalStress.size() != sizeStrains)
		  rNodalStress.resize(sizeStrains,false);
		noalias(rNodalStress) = ZeroVector(sizeStrains);
	      }else{
		std::cout<<"THIS node does not have NODAL_CAUCHY_STRESS... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_DEVIATORIC_CAUCHY_STRESS)){
		Vector& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS);
		if(rNodalStress.size() != sizeStrains)
		  rNodalStress.resize(sizeStrains,false);
		noalias(rNodalStress) = ZeroVector(sizeStrains);
	      }else{
		std::cout<<"THIS node does not have NODAL_DEVIATORIC_CAUCHY_STRESS... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_VOLUME)){
		itNode->FastGetSolutionStepValue(NODAL_VOLUME)=0;
	      }else{
		std::cout<<"THIS node does not have NODAL_VOLUME... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_MEAN_MESH_SIZE)){
		itNode->FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE)=0;
	      }else{
		std::cout<<"THIS node does not have NODAL_MEAN_MESH_SIZE... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      if(itNode->SolutionStepsDataHas(NODAL_FREESURFACE_AREA)){
		itNode->FastGetSolutionStepValue(NODAL_FREESURFACE_AREA)=0;
	      }else{
		std::cout<<"THIS node does not have NODAL_FREESURFACE_AREA... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_SFD_NEIGHBOURS)){
		Vector& rNodalSFDneighbours=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
		if(rNodalSFDneighbours.size() != sizeSDFNeigh)
		  rNodalSFDneighbours.resize(sizeSDFNeigh,false);
		noalias(rNodalSFDneighbours)=ZeroVector(sizeSDFNeigh);
	      }else{
		std::cout<<"THIS node does not have NODAL_SFD_NEIGHBOURS... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_SPATIAL_DEF_RATE)){
		Vector& rSpatialDefRate=itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
		if(rSpatialDefRate.size() != sizeStrains)
		  rSpatialDefRate.resize(sizeStrains,false);
		noalias(rSpatialDefRate)=ZeroVector(sizeStrains);
	      }else{
		std::cout<<"THIS node does not have NODAL_SPATIAL_DEF_RATE... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_DEFORMATION_GRAD)){
		Matrix& rFgrad=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
		if(rFgrad.size1() != dimension)
		  rFgrad.resize(dimension,dimension,false);
		noalias(rFgrad)=ZeroMatrix(dimension,dimension);
	      }else{
		std::cout<<"THIS node does not have NODAL_DEFORMATION_GRAD... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_DEFORMATION_GRAD_VEL)){
		Matrix& rFgradVel=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
		if(rFgradVel.size1() != dimension)
		  rFgradVel.resize(dimension,dimension,false);
		noalias(rFgradVel)=ZeroMatrix(dimension,dimension);
	      }else{
		std::cout<<"THIS node does not have NODAL_DEFORMATION_GRAD_VEL... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }

	    
	    }

	}
      
      }

       
      void ComputeNodalVolume()
      {
	/* std::cout << "ComputeNodalVolume: " << std::endl; */

	ModelPart& rModelPart = BaseType::GetModelPart();
	ElementsArrayType& pElements       = rModelPart.Elements();
	const unsigned int dimension =  rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
	
#ifdef _OPENMP
	int number_of_threads = omp_get_max_threads();
#else
	int number_of_threads = 1;
#endif

	vector<unsigned int> element_partition;
	OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

#pragma omp parallel
	{
	  int k = OpenMPUtils::ThisThread();
	  typename ElementsArrayType::iterator ElemBegin = pElements.begin() + element_partition[k];
	  typename ElementsArrayType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];

	  for (typename ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++)  //MSI: To be parallelized
	    {
	      Element::GeometryType& geometry = itElem->GetGeometry();
     
	      double elementalVolume=0; 

	      if(dimension==2){
		elementalVolume=geometry.Area()/3.0;
	      }else if(dimension==3){
		elementalVolume=geometry.Volume()*0.25;
	      }
	      // index = 0;
	      for (unsigned int i = 0; i <geometry.size(); i++)
		{
		  // index = i*dimension;

		  double& mass = geometry(i)->FastGetSolutionStepValue(NODAL_VOLUME);
		  mass+=elementalVolume;

		}

	    }
	}
      }

      void InitializeSolutionStep() override
      {
	this->FillNodalSFDVector();
      }  
   
      void FillNodalSFDVector()
      {

	ModelPart& rModelPart = BaseType::GetModelPart();
	
#pragma omp parallel
	{
	  ModelPart::NodeIterator NodesBegin;
	  ModelPart::NodeIterator NodesEnd;
	  OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

	  for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
	    {

	      if(itNode->Is(FLUID) && itNode->IsNot(SOLID)){
		/* unsigned int idNode=itNode->Id(); */
		InitializeNodalVariablesForRemeshedDomain(itNode);
	      }
	      
	      WeakPointerVector< Node < 3 > >& neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
	      unsigned int neighbourNodes=neighb_nodes.size()+1;
	      
	      Vector& rNodeOrderedNeighbours=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
	      if(rNodeOrderedNeighbours.size() != neighbourNodes)
		rNodeOrderedNeighbours.resize(neighbourNodes,false);
	      noalias(rNodeOrderedNeighbours)=ZeroVector(neighbourNodes);
	      rNodeOrderedNeighbours[0]=itNode->Id();
	      if(neighbourNodes>1){
		for(unsigned int k = 0; k< neighbourNodes-1; k++){
		  rNodeOrderedNeighbours[k+1]=neighb_nodes[k].Id();
		}
	      }	      
	      
	      /* if(itNode->SolutionStepsDataHas(NODAL_SFD_NEIGHBOURS_PRESSURE)){ */
	      /* 	Vector& rNodeOrderedNeighboursPressure=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_PRESSURE); */
	      /* 	if(rNodeOrderedNeighboursPressure.size() != neighbourNodes) */
	      /* 	  rNodeOrderedNeighboursPressure.resize(neighbourNodes,false); */
	      /* 	noalias(rNodeOrderedNeighboursPressure)=ZeroVector(neighbourNodes); */
	      /* }else{ */
	      /* 	std::cout<<"THIS node does not have NODAL_SFD_NEIGHBOURS_PRESSURE... "<<itNode->X()<<" "<<itNode->Y()<<std::endl; */
	      /* } */

	      /* if(itNode->SolutionStepsDataHas(NODAL_SFD_NEIGHBOURS_PRESSURE_ID)){ */
	      /* 	Vector& rNodeOrderedNeighboursPressureId=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_PRESSURE_ID); */
	      /* 	if(rNodeOrderedNeighboursPressureId.size() != neighbourNodes) */
	      /* 	  rNodeOrderedNeighboursPressureId.resize(neighbourNodes,false); */
	      /* 	noalias(rNodeOrderedNeighboursPressureId)=ZeroVector(neighbourNodes); */
	      /* 	const unsigned int xpos = itNode->GetDofPosition(PRESSURE); */
	      /* 	rNodeOrderedNeighboursPressureId[0]=itNode->GetDof(PRESSURE,xpos).EquationId(); */
	      /* 	if(neighbourNodes>1){ */
	      /* 	  for(unsigned int k = 0; k< neighbourNodes-1; k++){ */
	      /* 	    rNodeOrderedNeighboursPressureId[k+1]=neighb_nodes[k].GetDof(PRESSURE,xpos).EquationId(); */
	      /* 	  } */
	      /* 	} */
	      /* }else{ */
	      /* 	std::cout<<"THIS node does not have NODAL_SFD_NEIGHBOURS_PRESSURE_ID... "<<itNode->X()<<" "<<itNode->Y()<<std::endl; */
	      /* } */

	      /* 	      if(itNode->SolutionStepsDataHas(NODAL_SFD_NEIGHBOURS_VELOCITY_X_ID)){ */
	      /* 	Vector& rNodeOrderedNeighboursVelocityXId=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_VELOCITY_X_ID); */
	      /* 	Vector& rNodeOrderedNeighboursVelocityYId=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_VELOCITY_Y_ID); */
	      /* 	Vector& rNodeOrderedNeighboursVelocityZId=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_VELOCITY_Z_ID); */
	      /* 	if(rNodeOrderedNeighboursVelocityXId.size() != neighbourNodes){ */
	      /* 	  rNodeOrderedNeighboursVelocityXId.resize(neighbourNodes,false); */
	      /* 	  rNodeOrderedNeighboursVelocityYId.resize(neighbourNodes,false); */
	      /* 	  rNodeOrderedNeighboursVelocityZId.resize(neighbourNodes,false); */
	      /* 	} */
	      /* 	/\* noalias(rNodeOrderedNeighboursVelocityXId)=ZeroVector(neighbourNodes); *\/ */
	      /* 	const unsigned int xpos = itNode->GetDofPosition(VELOCITY_X); */
	      /* 	rNodeOrderedNeighboursVelocityXId[0]=itNode->GetDof(VELOCITY_X,xpos).EquationId(); */
	      /* 	rNodeOrderedNeighboursVelocityYId[0]=itNode->GetDof(VELOCITY_Y,xpos+1).EquationId(); */
	      /* 	rNodeOrderedNeighboursVelocityZId[0]=itNode->GetDof(VELOCITY_Z,xpos+2).EquationId(); */

	      /* 	if(neighbourNodes>1){ */
	      /* 	  for(unsigned int k = 0; k< neighbourNodes-1; k++){ */
	      /* 	    rNodeOrderedNeighboursVelocityXId[k+1]=neighb_nodes[k].GetDof(VELOCITY_X,xpos).EquationId(); */
	      /* 	    rNodeOrderedNeighboursVelocityYId[k+1]=neighb_nodes[k].GetDof(VELOCITY_Y,xpos+1).EquationId(); */
	      /* 	    rNodeOrderedNeighboursVelocityZId[k+1]=neighb_nodes[k].GetDof(VELOCITY_Z,xpos+2).EquationId(); */
	      /* 	  } */
	      /* 	} */
	      /* }else{ */
	      /* 	std::cout<<"THIS node does not have NODAL_SFD_NEIGHBOURS_VELOCITY_X_ID... "<<itNode->X()<<" "<<itNode->Y()<<std::endl; */
	      /* } */
		      

	    }
	}
      }


      void PrintVectors()
      {

	ModelPart& rModelPart = BaseType::GetModelPart();
	
#pragma omp parallel
	{
	  ModelPart::NodeIterator NodesBegin;
	  ModelPart::NodeIterator NodesEnd;
	  OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

	  for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
	    {
	      std::cout<<"THIS node "<<itNode->Id()<<" x="<<itNode->X()<<" y="<<itNode->Y();
	      std::cout<<"has neigh "<<itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER)<<std::endl;
	      std::cout<<"has area "<<itNode->FastGetSolutionStepValue(NODAL_VOLUME)<<std::endl;
	      std::cout<<"has h "<<itNode->FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE)<<std::endl;
	      std::cout<<"and shape derivatives "<<itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)<<std::endl;
	      std::cout<<"\n ";
	    }
	}
      }



      void ResetAndCheckNodalVariables()
      {

	ModelPart& rModelPart = BaseType::GetModelPart();
	const unsigned int dimension =  rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
	unsigned int sizeStrains=3*(dimension-1);
	
#pragma omp parallel
	{
	  ModelPart::NodeIterator NodesBegin;
	  ModelPart::NodeIterator NodesEnd;
	  OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

	  for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
	    {
	      WeakPointerVector< Node < 3 > >& neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
	      unsigned int neighbourNodes=neighb_nodes.size()+1;
	      unsigned int sizeSDFNeigh=neighbourNodes*dimension;
	      if(itNode->SolutionStepsDataHas(NODAL_VOLUME)){
		itNode->FastGetSolutionStepValue(NODAL_VOLUME)=0;
	      }else{
		std::cout<<"THIS node does not have NODAL_VOLUME... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_MEAN_MESH_SIZE)){
		itNode->FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE)=0;
	      }else{
		std::cout<<"THIS node does not have NODAL_MEAN_MESH_SIZE... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_FREESURFACE_AREA)){
		itNode->FastGetSolutionStepValue(NODAL_FREESURFACE_AREA)=0;
	      }else{
		std::cout<<"THIS node does not have NODAL_FREESURFACE_AREA... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_SFD_NEIGHBOURS)){
		Vector& rNodalSFDneighbours=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
		if(rNodalSFDneighbours.size() != sizeSDFNeigh)
		  rNodalSFDneighbours.resize(sizeSDFNeigh,false);
		noalias(rNodalSFDneighbours)=ZeroVector(sizeSDFNeigh);
	      }else{
		std::cout<<"THIS node does not have NODAL_SFD_NEIGHBOURS... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_SPATIAL_DEF_RATE)){
		Vector& rSpatialDefRate=itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
		if(rSpatialDefRate.size() != sizeStrains)
		  rSpatialDefRate.resize(sizeStrains,false);
		noalias(rSpatialDefRate)=ZeroVector(sizeStrains);
	      }else{
		std::cout<<"THIS node does not have NODAL_SPATIAL_DEF_RATE... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_VOLUMETRIC_DEF_RATE)){
		itNode->FastGetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE)=0;
	      }else{
		std::cout<<"THIS node does not have NODAL_VOLUMETRIC_DEF_RATE... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_DEFORMATION_GRAD)){
		Matrix& rFgrad=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
		if(rFgrad.size1() != dimension)
		  rFgrad.resize(dimension,dimension,false);
		noalias(rFgrad)=ZeroMatrix(dimension,dimension);
	      }else{
		std::cout<<"THIS node does not have NODAL_DEFORMATION_GRAD... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	      
	      if(itNode->SolutionStepsDataHas(NODAL_DEFORMATION_GRAD_VEL)){
		Matrix& rFgradVel=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
		if(rFgradVel.size1() != dimension)
		  rFgradVel.resize(dimension,dimension,false);
		noalias(rFgradVel)=ZeroMatrix(dimension,dimension);
	      }else{
		std::cout<<"THIS node does not have NODAL_DEFORMATION_GRAD_VEL... "<<itNode->X()<<" "<<itNode->Y()<<std::endl;
	      }
	    }
	}
      }

      void ResetNodalVariables()
      {

	ModelPart& rModelPart = BaseType::GetModelPart();
	const unsigned int dimension =  rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
	unsigned int sizeStrains=3*(dimension-1);
	
#pragma omp parallel
	{
	  ModelPart::NodeIterator NodesBegin;
	  ModelPart::NodeIterator NodesEnd;
	  OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

	  for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
	    {
	      unsigned int neighbourNodes=itNode->GetValue(NEIGHBOUR_NODES).size()+1;
	      unsigned int sizeSDFNeigh=neighbourNodes*dimension;
	      
	      Vector& rNodalSFDneighbours=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
	      noalias(rNodalSFDneighbours)=ZeroVector(sizeSDFNeigh);

	      Vector& rSpatialDefRate=itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
	      noalias(rSpatialDefRate)=ZeroVector(sizeStrains);
	
	      Matrix& rFgrad=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
	      noalias(rFgrad)=ZeroMatrix(dimension,dimension);

	      Matrix& rFgradVel=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
	      noalias(rFgradVel)=ZeroMatrix(dimension,dimension);
	      
	      itNode->FastGetSolutionStepValue(NODAL_VOLUME)=0;
	   
	      itNode->FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE)=0;
	
	      itNode->FastGetSolutionStepValue(NODAL_FREESURFACE_AREA)=0;

	      itNode->FastGetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE)=0;

	    }
	}
      }

      
      void InitializeNodalVariablesForRemeshedDomain(ModelPart::NodeIterator itNode)
      {

	ModelPart& rModelPart = BaseType::GetModelPart();
	const unsigned int dimension =  rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
	unsigned int sizeStrains=3*(dimension-1);
	WeakPointerVector< Node < 3 > >& neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
	unsigned int neighbourNodes=neighb_nodes.size()+1;
	unsigned int sizeSDFNeigh=neighbourNodes*dimension;

	Vector& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS);
	if(rNodalStress.size() != sizeStrains)
	  rNodalStress.resize(sizeStrains,false);
	noalias(rNodalStress) = ZeroVector(sizeStrains);
	
	Vector& rNodalDevStress = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS);
	if(rNodalDevStress.size() != sizeStrains)
	  rNodalDevStress.resize(sizeStrains,false);
	noalias(rNodalDevStress) = ZeroVector(sizeStrains);
	
	Vector& rNodalSFDneighbours=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
	if(rNodalSFDneighbours.size() != sizeSDFNeigh)
	  rNodalSFDneighbours.resize(sizeSDFNeigh,false);
	noalias(rNodalSFDneighbours)=ZeroVector(sizeSDFNeigh);

	Vector& rSpatialDefRate=itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
	if(rSpatialDefRate.size() != sizeStrains)
	  rSpatialDefRate.resize(sizeStrains,false);
	noalias(rSpatialDefRate)=ZeroVector(sizeStrains);
	
	Matrix& rFgrad=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
	if(rFgrad.size1() != dimension)
	  rFgrad.resize(dimension,dimension,false);
	noalias(rFgrad)=ZeroMatrix(dimension,dimension);

	Matrix& rFgradVel=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
	if(rFgradVel.size1() != dimension)
	  rFgradVel.resize(dimension,dimension,false);
	noalias(rFgradVel)=ZeroMatrix(dimension,dimension);

	itNode->FastGetSolutionStepValue(NODAL_VOLUME)=0;
	
	itNode->FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE)=0;
	
	itNode->FastGetSolutionStepValue(NODAL_FREESURFACE_AREA)=0;
	
	itNode->FastGetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE)=0;
	

      }


	
    void InitializeNonLinearIterations()
    {

      ModelPart& rModelPart = BaseType::GetModelPart();
      ElementsArrayType& pElements       = rModelPart.Elements();
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

#ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
#else
      int number_of_threads = 1;
#endif

      vector<unsigned int> element_partition;
      OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

#pragma omp parallel
      {
	int k = OpenMPUtils::ThisThread();
	typename ElementsArrayType::iterator ElemBegin = pElements.begin() + element_partition[k];
	typename ElementsArrayType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];

	for (typename ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++)  //MSI: To be parallelized
	  {
	    itElem->InitializeNonLinearIteration(rCurrentProcessInfo);
	  }
      }
    }



      
  void CalcNodalStrainsAndStresses()
  {

    ModelPart& rModelPart = BaseType::GetModelPart();

    const unsigned int dimension =  rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
    const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
#pragma omp parallel
    {
      ModelPart::NodeIterator NodesBegin;
      ModelPart::NodeIterator NodesEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

      for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
	{

	  if(itNode->Is(SOLID) || itNode->Is(FLUID)){
   
	    Vector nodalSFDneighboursId=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
	    double nodalVolume=itNode->FastGetSolutionStepValue(NODAL_VOLUME);
	    if(nodalVolume>0){
	      double deviatoricCoeff=0;
	      double volumetricCoeff=0;
	      if(itNode->Is(SOLID)){
		double youngModulus=itNode->FastGetSolutionStepValue(YOUNG_MODULUS);
		double poissonRatio=itNode->FastGetSolutionStepValue(POISSON_RATIO);
		deviatoricCoeff = timeInterval*youngModulus/(1.0+poissonRatio)*0.5;
		volumetricCoeff = timeInterval*poissonRatio*youngModulus/((1.0+poissonRatio)*(1.0-2.0*poissonRatio)) + 2.0*deviatoricCoeff/3.0;
	      }else if(itNode->Is(FLUID)){
		deviatoricCoeff = itNode->FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
		volumetricCoeff = timeInterval*itNode->FastGetSolutionStepValue(BULK_MODULUS);
	      }
	      double currFirstLame=volumetricCoeff - 2.0*deviatoricCoeff/3.0;

	      Vector& rNodalSFDneigh = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
	    
	      double theta=0.5;
	      
	      this->ComputeAndStoreNodalDeformationGradient(itNode, nodalSFDneighboursId, rNodalSFDneigh, theta);

	      Matrix Fgrad=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
	      Matrix FgradVel=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
	      double detFgrad=1.0;
	      Matrix InvFgrad=ZeroMatrix(dimension,dimension);
	      Matrix SpatialVelocityGrad=ZeroMatrix(dimension,dimension);
	      //Inverse

	      MathUtils<double>::InvertMatrix(Fgrad,InvFgrad,detFgrad);
	  
	      //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
	      SpatialVelocityGrad=prod(FgradVel,InvFgrad);

	      if(dimension==2){
		itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0]=SpatialVelocityGrad(0,0);
		itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1]=SpatialVelocityGrad(1,1);
		itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2]=0.5*(SpatialVelocityGrad(1,0)+SpatialVelocityGrad(0,1));

		double DefX=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0];
		double DefY=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1];
		double DefXY=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2];

		double DefVol=DefX+DefY;
	      
		itNode->GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE)=DefVol;

		double nodalSigmaTot_xx= currFirstLame*DefVol + 2.0*deviatoricCoeff*DefX;
		double nodalSigmaTot_yy= currFirstLame*DefVol + 2.0*deviatoricCoeff*DefY;
		double nodalSigmaTot_xy= 2.0*deviatoricCoeff*DefXY;

		double nodalSigmaDev_xx= 2.0*deviatoricCoeff*(DefX - DefVol/3.0);
		double nodalSigmaDev_yy= 2.0*deviatoricCoeff*(DefY - DefVol/3.0);
		double nodalSigmaDev_xy= 2.0*deviatoricCoeff*DefXY;
		
		if(itNode->Is(SOLID)){
		  nodalSigmaTot_xx+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[0];
		  nodalSigmaTot_yy+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[1];
		  nodalSigmaTot_xy+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[2];

		  nodalSigmaDev_xx+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[0];
		  nodalSigmaDev_yy+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[1];
		  nodalSigmaDev_xy+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[2];
		}
 
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[0]=nodalSigmaTot_xx;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[1]=nodalSigmaTot_yy;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[2]=nodalSigmaTot_xy;

		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[0]=nodalSigmaDev_xx;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[1]=nodalSigmaDev_yy;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[2]=nodalSigmaDev_xy;
	    
	      }else if (dimension==3){
	   
		itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0]=SpatialVelocityGrad(0,0);
		itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1]=SpatialVelocityGrad(1,1);
		itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2]=SpatialVelocityGrad(2,2);
		itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[3]=0.5*(SpatialVelocityGrad(1,0)+SpatialVelocityGrad(0,1));
		itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[4]=0.5*(SpatialVelocityGrad(2,0)+SpatialVelocityGrad(0,2));
		itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[5]=0.5*(SpatialVelocityGrad(2,1)+SpatialVelocityGrad(1,2));

		double DefX=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0];
		double DefY=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1];
		double DefZ=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2];
		double DefXY=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[3];
		double DefXZ=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[4];
		double DefYZ=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[5];

		double DefVol=DefX+DefY+DefZ;

		itNode->GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE)=DefVol;

		double nodalSigmaTot_xx= currFirstLame*DefVol + 2.0*deviatoricCoeff*DefX;
		double nodalSigmaTot_yy= currFirstLame*DefVol + 2.0*deviatoricCoeff*DefY;
		double nodalSigmaTot_zz= currFirstLame*DefVol + 2.0*deviatoricCoeff*DefZ;
		double nodalSigmaTot_xy= 2.0*deviatoricCoeff*DefXY;
		double nodalSigmaTot_xz= 2.0*deviatoricCoeff*DefXZ;
		double nodalSigmaTot_yz= 2.0*deviatoricCoeff*DefYZ;

		double nodalSigmaDev_xx= 2.0*deviatoricCoeff*(DefX - DefVol/3.0);
		double nodalSigmaDev_yy= 2.0*deviatoricCoeff*(DefY - DefVol/3.0);
		double nodalSigmaDev_zz= 2.0*deviatoricCoeff*(DefZ - DefVol/3.0);
		double nodalSigmaDev_xy= 2.0*deviatoricCoeff*DefXY;
		double nodalSigmaDev_xz= 2.0*deviatoricCoeff*DefXZ;
		double nodalSigmaDev_yz= 2.0*deviatoricCoeff*DefYZ;

		if(itNode->Is(SOLID)){
		  nodalSigmaTot_xx+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[0];
		  nodalSigmaTot_yy+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[1];
		  nodalSigmaTot_zz+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[2];
		  nodalSigmaTot_xy+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[3];
		  nodalSigmaTot_xz+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[4];
		  nodalSigmaTot_yz+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[5];  

		  nodalSigmaDev_xx+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[0];
		  nodalSigmaDev_yy+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[1];
		  nodalSigmaDev_zz+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[2];
		  nodalSigmaDev_xy+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[3];
		  nodalSigmaDev_xz+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[4];
		  nodalSigmaDev_yz+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[5];
		}
		  
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[0]=nodalSigmaTot_xx;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[1]=nodalSigmaTot_yy;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[2]=nodalSigmaTot_zz;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[3]=nodalSigmaTot_xy;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[4]=nodalSigmaTot_xz;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[5]=nodalSigmaTot_yz;
		
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[0]=nodalSigmaDev_xx;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[1]=nodalSigmaDev_yy;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[2]=nodalSigmaDev_zz;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[3]=nodalSigmaDev_xy;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[4]=nodalSigmaDev_xz;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[5]=nodalSigmaDev_yz;
	
	      }
	    }
	    else{
	      /* unsigned int  idNode=itNode->Id(); */
	      InitializeNodalVariablesForRemeshedDomain(itNode);
	     
	    }
	  }

	}
    }

    /* std::cout << "Calc Nodal Strains And Stresses DONE " << std::endl; */

  }


  
  void CalcNodalStresses()
  {

    /* std::cout << "Calc Nodal Strains And Stresses " << std::endl; */

    ModelPart& rModelPart = BaseType::GetModelPart();

    const unsigned int dimension =  rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
    ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
    const double timeInterval = rCurrentProcessInfo[DELTA_TIME];
#pragma omp parallel
    {
      ModelPart::NodeIterator NodesBegin;
      ModelPart::NodeIterator NodesEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

      for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
	{

	  if(itNode->Is(SOLID) || itNode->Is(FLUID)){
   	    double nodalVolume=itNode->FastGetSolutionStepValue(NODAL_VOLUME);
	    if(nodalVolume>0){
	      double deviatoricCoeff=0;
	      double volumetricCoeff=0;
	      if(itNode->Is(SOLID)){
		double youngModulus=itNode->FastGetSolutionStepValue(YOUNG_MODULUS);
		double poissonRatio=itNode->FastGetSolutionStepValue(POISSON_RATIO);
		deviatoricCoeff = timeInterval*youngModulus/(1.0+poissonRatio)*0.5;
		volumetricCoeff = timeInterval*poissonRatio*youngModulus/((1.0+poissonRatio)*(1.0-2.0*poissonRatio)) + 2.0*deviatoricCoeff/3.0;
	      }else if(itNode->Is(FLUID)){
		deviatoricCoeff = itNode->FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
		volumetricCoeff = timeInterval*itNode->FastGetSolutionStepValue(BULK_MODULUS);
	      }
	      double currFirstLame=volumetricCoeff - 2.0*deviatoricCoeff/3.0;
	      
	      if(dimension==2){

		double DefX=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0];
		double DefY=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1];
		double DefXY=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2];

		double DefVol=itNode->GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE);

		double nodalSigmaTot_xx= currFirstLame*DefVol + 2.0*deviatoricCoeff*DefX;
		double nodalSigmaTot_yy= currFirstLame*DefVol + 2.0*deviatoricCoeff*DefY;
		double nodalSigmaTot_xy= 2.0*deviatoricCoeff*DefXY;

		double nodalSigmaDev_xx= 2.0*deviatoricCoeff*(DefX - DefVol/3.0);
		double nodalSigmaDev_yy= 2.0*deviatoricCoeff*(DefY - DefVol/3.0);
		double nodalSigmaDev_xy= 2.0*deviatoricCoeff*DefXY;
		
		if(itNode->Is(SOLID)){
		  nodalSigmaTot_xx+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[0];
		  nodalSigmaTot_yy+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[1];
		  nodalSigmaTot_xy+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[2];

		  nodalSigmaDev_xx+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[0];
		  nodalSigmaDev_yy+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[1];
		  nodalSigmaDev_xy+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[2];
		}
 
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[0]=nodalSigmaTot_xx;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[1]=nodalSigmaTot_yy;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[2]=nodalSigmaTot_xy;

		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[0]=nodalSigmaDev_xx;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[1]=nodalSigmaDev_yy;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[2]=nodalSigmaDev_xy;
	    
	      }else if (dimension==3){
	   
		double DefX=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0];
		double DefY=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1];
		double DefZ=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2];
		double DefXY=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[3];
		double DefXZ=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[4];
		double DefYZ=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[5];

		double DefVol=itNode->GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE);

		double nodalSigmaTot_xx= currFirstLame*DefVol + 2.0*deviatoricCoeff*DefX;
		double nodalSigmaTot_yy= currFirstLame*DefVol + 2.0*deviatoricCoeff*DefY;
		double nodalSigmaTot_zz= currFirstLame*DefVol + 2.0*deviatoricCoeff*DefZ;
		double nodalSigmaTot_xy= 2.0*deviatoricCoeff*DefXY;
		double nodalSigmaTot_xz= 2.0*deviatoricCoeff*DefXZ;
		double nodalSigmaTot_yz= 2.0*deviatoricCoeff*DefYZ;

		double nodalSigmaDev_xx= 2.0*deviatoricCoeff*(DefX - DefVol/3.0);
		double nodalSigmaDev_yy= 2.0*deviatoricCoeff*(DefY - DefVol/3.0);
		double nodalSigmaDev_zz= 2.0*deviatoricCoeff*(DefZ - DefVol/3.0);
		double nodalSigmaDev_xy= 2.0*deviatoricCoeff*DefXY;
		double nodalSigmaDev_xz= 2.0*deviatoricCoeff*DefXZ;
		double nodalSigmaDev_yz= 2.0*deviatoricCoeff*DefYZ;

		if(itNode->Is(SOLID)){
		  nodalSigmaTot_xx+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[0];
		  nodalSigmaTot_yy+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[1];
		  nodalSigmaTot_zz+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[2];
		  nodalSigmaTot_xy+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[3];
		  nodalSigmaTot_xz+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[4];
		  nodalSigmaTot_yz+=itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,1)[5];  

		  nodalSigmaDev_xx+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[0];
		  nodalSigmaDev_yy+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[1];
		  nodalSigmaDev_zz+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[2];
		  nodalSigmaDev_xy+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[3];
		  nodalSigmaDev_xz+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[4];
		  nodalSigmaDev_yz+=itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,1)[5];
		}
		  
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[0]=nodalSigmaTot_xx;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[1]=nodalSigmaTot_yy;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[2]=nodalSigmaTot_zz;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[3]=nodalSigmaTot_xy;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[4]=nodalSigmaTot_xz;
		itNode->GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[5]=nodalSigmaTot_yz;
		
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[0]=nodalSigmaDev_xx;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[1]=nodalSigmaDev_yy;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[2]=nodalSigmaDev_zz;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[3]=nodalSigmaDev_xy;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[4]=nodalSigmaDev_xz;
		itNode->GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[5]=nodalSigmaDev_yz;
	
	      }
	    }
	    else{
	      /* unsigned int  idNode=itNode->Id(); */
	      InitializeNodalVariablesForRemeshedDomain(itNode);
	     
	    }
	  }

	}
    }

    /* std::cout << "Calc Nodal Strains And Stresses DONE " << std::endl; */

  }


      
  void CalcNodalStrains()
  {

    /* std::cout << "Calc Nodal Strains  " << std::endl; */
    ModelPart& rModelPart = BaseType::GetModelPart();

    const unsigned int dimension =  rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

#pragma omp parallel
    {
      ModelPart::NodeIterator NodesBegin;
      ModelPart::NodeIterator NodesEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

      for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
	{

	  if(itNode->Is(SOLID) || itNode->Is(FLUID)){
   
	    Vector nodalSFDneighboursId=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);

	    Vector& rNodalSFDneigh = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
	    
	    double theta=1.0;
	    
	    this->ComputeAndStoreNodalDeformationGradient(itNode, nodalSFDneighboursId, rNodalSFDneigh, theta);
	  
	    Matrix Fgrad=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
	    Matrix FgradVel=itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
	    double detFgrad=1.0;
	    Matrix InvFgrad=ZeroMatrix(dimension,dimension);
	    Matrix SpatialVelocityGrad=ZeroMatrix(dimension,dimension);
	    //Inverse

	    
	    MathUtils<double>::InvertMatrix(Fgrad,InvFgrad,detFgrad);
	  
	    //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
	    SpatialVelocityGrad=prod(FgradVel,InvFgrad);

	    if(dimension==2){
	      itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0]=SpatialVelocityGrad(0,0);
	      itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1]=SpatialVelocityGrad(1,1);
	      itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2]=0.5*(SpatialVelocityGrad(1,0)+SpatialVelocityGrad(0,1));

	      double DefX=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0];
	      double DefY=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1];

	      double DefVol=DefX+DefY;
	      
	      itNode->GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE)=DefVol;

	    
	    }else if (dimension==3){

	   
	      itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0]=SpatialVelocityGrad(0,0);
	      itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1]=SpatialVelocityGrad(1,1);
	      itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2]=SpatialVelocityGrad(2,2);
	      itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[3]=0.5*(SpatialVelocityGrad(1,0)+SpatialVelocityGrad(0,1));
	      itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[4]=0.5*(SpatialVelocityGrad(2,0)+SpatialVelocityGrad(0,2));
	      itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[5]=0.5*(SpatialVelocityGrad(2,1)+SpatialVelocityGrad(1,2));

	      double DefX=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[0];
	      double DefY=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[1];
	      double DefZ=itNode->GetSolutionStepValue(NODAL_SPATIAL_DEF_RATE)[2];

	      double DefVol=DefX+DefY+DefZ;

	      itNode->GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE)=DefVol;
	    
	    }

	  }

	}
    }

    /* std::cout << "Calc Nodal Strains And Stresses DONE " << std::endl; */

  }

  void ComputeAndStoreNodalDeformationGradient(ModelPart::NodeIterator itNode, Vector nodalSFDneighboursId, Vector rNodalSFDneigh, double theta){
    
    KRATOS_TRY;

    ModelPart& rModelPart = BaseType::GetModelPart();
    const unsigned int dimension =  rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
    /* unsigned int idThisNode=nodalSFDneighboursId[0]; */
    const unsigned int neighSize = nodalSFDneighboursId.size();
    Matrix Fgrad=ZeroMatrix(dimension,dimension);
    Matrix FgradVel=ZeroMatrix(dimension,dimension);
    WeakPointerVector< Node < 3 > >& neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);

    if(dimension==2){
      
      double dNdXi=rNodalSFDneigh[0];
      double dNdYi=rNodalSFDneigh[1];

      double VelocityX = itNode->FastGetSolutionStepValue(VELOCITY_X,0)*theta +	itNode->FastGetSolutionStepValue(VELOCITY_X,1)*(1-theta);
      double VelocityY = itNode->FastGetSolutionStepValue(VELOCITY_Y,0)*theta +	itNode->FastGetSolutionStepValue(VELOCITY_Y,1)*(1-theta);

      Fgrad(0,0)+=dNdXi*itNode->X();
      Fgrad(0,1)+=dNdYi*itNode->X();
      Fgrad(1,0)+=dNdXi*itNode->Y();
      Fgrad(1,1)+=dNdYi*itNode->Y();

      FgradVel(0,0)+=dNdXi*VelocityX;
      FgradVel(0,1)+=dNdYi*VelocityX;
      FgradVel(1,0)+=dNdXi*VelocityY;
      FgradVel(1,1)+=dNdYi*VelocityY;

      unsigned int firstRow=2;

      if(neighSize>0){
	for (unsigned int i = 0; i< neighSize-1; i++)
	  {
	    dNdXi=rNodalSFDneigh[firstRow];
	    dNdYi=rNodalSFDneigh[firstRow+1];

	    VelocityX = neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_X,0)*theta + neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_X,1)*(1-theta);
	    VelocityY = neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_Y,0)*theta + neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_Y,1)*(1-theta);

	    Fgrad(0,0)+=dNdXi*neighb_nodes[i].X();
	    Fgrad(0,1)+=dNdYi*neighb_nodes[i].X();
	    Fgrad(1,0)+=dNdXi*neighb_nodes[i].Y();
	    Fgrad(1,1)+=dNdYi*neighb_nodes[i].Y();
		    
	    FgradVel(0,0)+=dNdXi*VelocityX;
	    FgradVel(0,1)+=dNdYi*VelocityX;
	    FgradVel(1,0)+=dNdXi*VelocityY;
	    FgradVel(1,1)+=dNdYi*VelocityY;
		  
	    firstRow+=2;
		  
	  }
      }

    }else{


      double dNdXi=rNodalSFDneigh[0];
      double dNdYi=rNodalSFDneigh[1];
      double dNdZi=rNodalSFDneigh[2];
		  
      double VelocityX = itNode->FastGetSolutionStepValue(VELOCITY_X,0)*theta +	itNode->FastGetSolutionStepValue(VELOCITY_X,1)*(1-theta);
      double VelocityY = itNode->FastGetSolutionStepValue(VELOCITY_Y,0)*theta +	itNode->FastGetSolutionStepValue(VELOCITY_Y,1)*(1-theta);
      double VelocityZ = itNode->FastGetSolutionStepValue(VELOCITY_Z,0)*theta +	itNode->FastGetSolutionStepValue(VELOCITY_Z,1)*(1-theta);
		    
      Fgrad(0,0)+=dNdXi*itNode->X();
      Fgrad(0,1)+=dNdYi*itNode->X();
      Fgrad(0,2)+=dNdZi*itNode->X();
		  
      Fgrad(1,0)+=dNdXi*itNode->Y();
      Fgrad(1,1)+=dNdYi*itNode->Y();
      Fgrad(1,2)+=dNdZi*itNode->Y();
		  
      Fgrad(2,0)+=dNdXi*itNode->Z();
      Fgrad(2,1)+=dNdYi*itNode->Z();
      Fgrad(2,2)+=dNdZi*itNode->Z();

		    
      FgradVel(0,0)+=dNdXi*VelocityX;
      FgradVel(0,1)+=dNdYi*VelocityX;
      FgradVel(0,2)+=dNdZi*VelocityX;
	  
      FgradVel(1,0)+=dNdXi*VelocityY;
      FgradVel(1,1)+=dNdYi*VelocityY;
      FgradVel(1,2)+=dNdZi*VelocityY;
		  
      FgradVel(2,0)+=dNdXi*VelocityZ;
      FgradVel(2,1)+=dNdYi*VelocityZ;
      FgradVel(2,2)+=dNdZi*VelocityZ;
		  	  
      unsigned int firstRow=3;

      if(neighSize>0){
	for (unsigned int i = 0; i< neighSize-1; i++)
	  {
	  
	    dNdXi=rNodalSFDneigh[firstRow];
	    dNdYi=rNodalSFDneigh[firstRow+1];
	    dNdZi=rNodalSFDneigh[firstRow+2];

	    VelocityX = neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_X,0)*theta + neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_X,1)*(1-theta);
	    VelocityY = neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_Y,0)*theta + neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_Y,1)*(1-theta);
	    VelocityZ = neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_Z,0)*theta + neighb_nodes[i].FastGetSolutionStepValue(VELOCITY_Z,1)*(1-theta);
		  		    
	    Fgrad(0,0)+=dNdXi*neighb_nodes[i].X();
	    Fgrad(0,1)+=dNdYi*neighb_nodes[i].X();
	    Fgrad(0,2)+=dNdZi*neighb_nodes[i].X();
		  
	    Fgrad(1,0)+=dNdXi*neighb_nodes[i].Y();
	    Fgrad(1,1)+=dNdYi*neighb_nodes[i].Y();
	    Fgrad(1,2)+=dNdZi*neighb_nodes[i].Y();
		  
	    Fgrad(2,0)+=dNdXi*neighb_nodes[i].Z();
	    Fgrad(2,1)+=dNdYi*neighb_nodes[i].Z();
	    Fgrad(2,2)+=dNdZi*neighb_nodes[i].Z();

		    
	    FgradVel(0,0)+=dNdXi*VelocityX;
	    FgradVel(0,1)+=dNdYi*VelocityX;
	    FgradVel(0,2)+=dNdZi*VelocityX;
	  
	    FgradVel(1,0)+=dNdXi*VelocityY;
	    FgradVel(1,1)+=dNdYi*VelocityY;
	    FgradVel(1,2)+=dNdZi*VelocityY;
		  
	    FgradVel(2,0)+=dNdXi*VelocityZ;
	    FgradVel(2,1)+=dNdYi*VelocityZ;
	    FgradVel(2,2)+=dNdZi*VelocityZ;
		  
	    firstRow+=3;

		  
	  }
      }
	      

    }

    itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD)=Fgrad;
    itNode->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL)=FgradVel;
    KRATOS_CATCH("");
  }



    
    void UpdateTopology(ModelPart& rModelPart, unsigned int echoLevel)
    {
      KRATOS_TRY;
      
      /* this->CalculateDisplacements(); */
      this->CalculateDisplacementsAndResetNodalVariables();
      BaseType::MoveMesh();
      BoundaryNormalsCalculationUtilities BoundaryComputation;
      BoundaryComputation.CalculateWeightedBoundaryNormals(rModelPart, echoLevel);
      
      KRATOS_CATCH("");
  
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

      for (ModelPart::NodeIterator i = rModelPart.NodesBegin(); i != rModelPart.NodesEnd(); ++i)
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
      
      for (ModelPart::NodeIterator i = rModelPart.NodesBegin(); i != rModelPart.NodesEnd(); ++i)
        {

	  array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY, 0);
	  array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);

	  array_1d<double, 3 > & CurrentAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
	  array_1d<double, 3 > & PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

	  /* if((i)->IsNot(ISOLATED) || (i)->Is(SOLID)){ */
	  if((i)->IsNot(ISOLATED) && (i)->IsNot(RIGID)){
	    UpdateAccelerations (CurrentAcceleration, CurrentVelocity, PreviousAcceleration, PreviousVelocity,BDFcoeffs);
	  }else if((i)->Is(RIGID)){
	    array_1d<double, 3>  Zeros(3,0.0);
	    (i)->FastGetSolutionStepValue(ACCELERATION,0) = Zeros;
	    (i)->FastGetSolutionStepValue(ACCELERATION,1) = Zeros;
	  }else {
	    (i)->FastGetSolutionStepValue(NODAL_VOLUME) = 0.0;
	    (i)->FastGetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE) = 0.0;
	    (i)->FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE) = 0.0;
	    (i)->FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) = 0.0;
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
      
      for (ModelPart::NodeIterator i = rModelPart.NodesBegin(); i != rModelPart.NodesEnd(); ++i)
        {

	  array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY, 0);
	  array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);

	  array_1d<double, 3 > & CurrentDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT, 0);
	  array_1d<double, 3 > & PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
	  
	  /* if( i->IsFixed(DISPLACEMENT_X) == false ) */
	    CurrentDisplacement[0] = 0.5* TimeStep *(CurrentVelocity[0]+PreviousVelocity[0]) + PreviousDisplacement[0];	  

	  /* if( i->IsFixed(DISPLACEMENT_Y) == false ) */
	    CurrentDisplacement[1] = 0.5* TimeStep *(CurrentVelocity[1]+PreviousVelocity[1]) + PreviousDisplacement[1];

	  /* if( i->IsFixed(DISPLACEMENT_Z) == false ) */
	    CurrentDisplacement[2] = 0.5* TimeStep *(CurrentVelocity[2]+PreviousVelocity[2]) + PreviousDisplacement[2];

        }
    }

    
    void CalculateDisplacementsAndResetNodalVariables()
    {
      ModelPart& rModelPart = BaseType::GetModelPart();
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
      const double TimeStep = rCurrentProcessInfo[DELTA_TIME];
      const unsigned int dimension =  rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
      unsigned int sizeStrains=3*(dimension-1);
 
#pragma omp parallel
      {
	ModelPart::NodeIterator NodesBegin;
	ModelPart::NodeIterator NodesEnd;
	OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

	for (ModelPart::NodeIterator i = NodesBegin; i != NodesEnd; ++i)
	  {
	    array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY, 0);
	    array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);

	    array_1d<double, 3 > & CurrentDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT, 0);
	    array_1d<double, 3 > & PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
	  
	    /* if( i->IsFixed(DISPLACEMENT_X) == false ) */
	    CurrentDisplacement[0] = 0.5* TimeStep *(CurrentVelocity[0]+PreviousVelocity[0]) + PreviousDisplacement[0];	  

	    /* if( i->IsFixed(DISPLACEMENT_Y) == false ) */
	    CurrentDisplacement[1] = 0.5* TimeStep *(CurrentVelocity[1]+PreviousVelocity[1]) + PreviousDisplacement[1];

	    /* if( i->IsFixed(DISPLACEMENT_Z) == false ) */
	    if(dimension==3){
	      CurrentDisplacement[2] = 0.5* TimeStep *(CurrentVelocity[2]+PreviousVelocity[2]) + PreviousDisplacement[2];
	    }

	    ///// reset Nodal variables //////
	    unsigned int neighbourNodes=i->GetValue(NEIGHBOUR_NODES).size()+1;
	    unsigned int sizeSDFNeigh=neighbourNodes*dimension;
	    i->FastGetSolutionStepValue(NODAL_VOLUME)=0;
	   
	    i->FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE)=0;
	
	    i->FastGetSolutionStepValue(NODAL_FREESURFACE_AREA)=0;

	    i->FastGetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE)=0;

	    Vector& rNodalSFDneighbours=i->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS);
	    noalias(rNodalSFDneighbours)=ZeroVector(sizeSDFNeigh);

	    Vector& rSpatialDefRate=i->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
	    noalias(rSpatialDefRate)=ZeroVector(sizeStrains);
	
	    Matrix& rFgrad=i->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD);
	    noalias(rFgrad)=ZeroMatrix(dimension,dimension);

	    Matrix& rFgradVel=i->FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL);
	    noalias(rFgradVel)=ZeroMatrix(dimension,dimension);
	      

	  }
      }
    }

  

   void UpdatePressureAccelerations()
   {
     this->CalculateAccelerations(); 
     this->CalculatePressureVelocity();
     this->CalculatePressureAcceleration();

   }

   void Clear() override
    {
        mpMomentumStrategy->Clear();
        mpPressureStrategy->Clear();
    }


    ///@}
    ///@name Access
    ///@{

    void SetEchoLevel(int Level) override
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
        buffer << "NodalTwoStepVPStrategy" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "NodalTwoStepVPStrategy";}

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

      /* std::cout<<"---- m o m e n t u m   e q u a t i o n s ----"<<std::endl; */
      if(it==0){
	mpMomentumStrategy->InitializeSolutionStep();
	/* this->SetNeighboursVelocityId(); */
      }
      
      NormDv = mpMomentumStrategy->Solve(); 
	  
      if (BaseType::GetEchoLevel() > 1 && Rank == 0)
	std::cout<<"-------------- s o l v e d ! ------------------"<<std::endl;

      double DvErrorNorm = 0; 
      ConvergedMomentum = this->CheckVelocityConvergence(NormDv,DvErrorNorm);
      // Check convergence

      if(it==maxIt-1){
	std::cout<<"         iteration("<<it<<") Final Velocity error: "<< DvErrorNorm <<" velTol: " << mVelocityTolerance<< std::endl;
	fixedTimeStep=this->FixTimeStepMomentum(DvErrorNorm);
      }else{
	std::cout<<"iteration("<<it<<")  Velocity error: "<< DvErrorNorm <<" velTol: " << mVelocityTolerance<< std::endl;
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

      /* std::cout<<"     ---- c o n t i n u i t y   e q u a t i o n ----"<<std::endl; */
 
      if(it==0){
	mpPressureStrategy->InitializeSolutionStep();
	/* this->SetNeighboursPressureId(); */
      }

      NormDp = mpPressureStrategy->Solve();

      if (BaseType::GetEchoLevel() > 0 && Rank == 0)
	std::cout << "The norm of pressure is: " << NormDp << std::endl;

      double DpErrorNorm = 0; 
      ConvergedContinuity = this->CheckPressureConvergence(NormDp,DpErrorNorm);

      // Check convergence
      if(it==maxIt-1){
	std::cout<<"                  iteration("<<it<<") Final Pressure error: "<<DpErrorNorm <<" presTol: "<<mPressureTolerance << std::endl;
      	ConvergedContinuity=this->FixTimeStepContinuity(DpErrorNorm);
      }else{
	std::cout<<"iteration("<<it<<") Pressure error: "<<DpErrorNorm <<" presTol: "<<mPressureTolerance << std::endl;
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
	}
	/* else{ */
	/*   std::cout<<"Velocity error: "<< errorNormDv <<" velTol: " << mVelocityTolerance<< std::endl; */
	/* } */
	
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
	}
	/* else{ */
        /*     std::cout<<"         Pressure error: "<<errorNormDp <<" presTol: "<<mPressureTolerance << std::endl; */
	/* } */

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
      double minTolerance=0.01;
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
            KRATOS_THROW_ERROR(std::runtime_error,"NodalTwoStepVPStrategy error: No Velocity strategy defined in FractionalStepSettings","");
        }

        bool HavePressStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Pressure,mpPressureStrategy);

        if (HavePressStrategy)
        {
            rSolverConfig.FindTolerance(SolverSettingsType::Pressure,mPressureTolerance);
            rSolverConfig.FindMaxIter(SolverSettingsType::Pressure,mMaxPressureIter);
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,"NodalTwoStepVPStrategy error: No Pressure strategy defined in FractionalStepSettings","");
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
    NodalTwoStepVPStrategy& operator=(NodalTwoStepVPStrategy const& rOther){}

    /// Copy constructor.
    NodalTwoStepVPStrategy(NodalTwoStepVPStrategy const& rOther){}


    ///@}

}; /// Class NodalTwoStepVPStrategy

///@}
///@name Type Definitions
///@{


///@}

///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_NODAL_TWO_STEP_V_P_STRATEGY_H
