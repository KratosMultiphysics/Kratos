//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:       BSD License
//                Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi, Alessandro Franci
//
//

#if !defined(KRATOS_NODAL_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_CONTINUITY )
#define  KRATOS_NODAL_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_CONTINUITY

/* System includes */
#include <set>
#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */
// #define USE_GOOGLE_HASH
#ifdef USE_GOOGLE_HASH
#include "sparsehash/dense_hash_set" //included in external libraries
#else
#include <unordered_set>
#endif

/* Project includes */
#include "utilities/timer.h"
#include "includes/define.h"
#include "includes/key_hash.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"

#include "pfem_fluid_dynamics_application_variables.h"


namespace Kratos
{

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

  /**
   * @class NodalResidualBasedEliminationBuilderAndSolverContinuity
   * @ingroup KratosCore
   * @brief Current class provides an implementation for standard builder and solving operations.
   * @details The RHS is constituted by the unbalanced loads (residual)
   * Degrees of freedom are reordered putting the restrained degrees of freedom at
   * the end of the system ordered in reverse order with respect to the DofSet.
   * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
   * this information.
   * Calculation of the reactions involves a cost very similiar to the calculation of the total residual
   * @author Riccardo Rossi
   */
  template<class TSparseSpace,
    class TDenseSpace, //= DenseSpace<double>,
    class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
    >
    class NodalResidualBasedEliminationBuilderAndSolverContinuity
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
    {
    public:
      ///@name Type Definitions
      ///@{
      KRATOS_CLASS_POINTER_DEFINITION(NodalResidualBasedEliminationBuilderAndSolverContinuity);

      typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

      typedef typename BaseType::TSchemeType TSchemeType;

      typedef typename BaseType::TDataType TDataType;

      typedef typename BaseType::DofsArrayType DofsArrayType;

      typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

      typedef typename BaseType::TSystemVectorType TSystemVectorType;

      typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

      typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

      typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
      typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

      typedef Node<3> NodeType;

      typedef typename BaseType::NodesArrayType NodesArrayType;
      typedef typename BaseType::ElementsArrayType ElementsArrayType;
      typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

      typedef typename BaseType::ElementsContainerType ElementsContainerType;

      typedef Vector VectorType;


      ///@}
      ///@name Life Cycle
      ///@{

      /** Constructor.
       */
    NodalResidualBasedEliminationBuilderAndSolverContinuity(
				typename TLinearSolver::Pointer pNewLinearSystemSolver)
      : BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
	{
	  //         KRATOS_INFO("NodalResidualBasedEliminationBuilderAndSolverContinuity") << "Using the standard builder and solver " << std::endl;
	}

      /** Destructor.
       */
      ~NodalResidualBasedEliminationBuilderAndSolverContinuity() override
	{
	}

      ///@}
      ///@name Operators
      ///@{

      ///@}
      ///@name Operations
      ///@{
  
     
      void BuildNodally(
			typename TSchemeType::Pointer pScheme,
			ModelPart& rModelPart,
			TSystemMatrixType& A,
			TSystemVectorType& b)
      {
	KRATOS_TRY

	  KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

	/* std::cout<<"Building LHS and RHS of Momentum Equation Nodally"<<std::endl; */
     
        //contributions to the system
	LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
	LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

	//vector containing the localization in the system of the different terms
	Element::EquationIdVectorType EquationId;

	ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

	const unsigned int dimension =  rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
	const double timeInterval = CurrentProcessInfo[DELTA_TIME];

	double deviatoricCoeff=1.0;
	double volumetricCoeff=1.0;
	double pressure=0;
	double deltaPressure=0;
	double meanMeshSize=0;
	double characteristicLength=0;
	double density=0;
	double nodalVelocity=0;		 
	double tauStab=0;
	double dNdXi=0;
	double dNdYi=0;
	double dNdZi=0;
	double dNdXj=0;
	double dNdYj=0;
	double dNdZj=0;
	bool stabilizationNeeded=false;
	unsigned int firstRow=0;
	unsigned int firstCol=0;
	
	/* #pragma omp parallel */
	{
	  ModelPart::NodeIterator NodesBegin;
	  ModelPart::NodeIterator NodesEnd;
	  OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

	  for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
	    {

	      WeakPointerVector< Node < 3 > >& neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
	      const unsigned int neighSize = neighb_nodes.size() +1 ;

	      if(neighSize>1){
		
		const double nodalVolume=itNode->FastGetSolutionStepValue(NODAL_VOLUME);

		LHS_Contribution= ZeroMatrix(neighSize,neighSize);
		RHS_Contribution= ZeroVector(neighSize);
		
		if (EquationId.size() != neighSize)
		  EquationId.resize(neighSize, false);

		if(itNode->Is(SOLID)){
		
		  deviatoricCoeff = timeInterval*itNode->FastGetSolutionStepValue(YOUNG_MODULUS)/(1.0+itNode->FastGetSolutionStepValue(POISSON_RATIO))*0.5;
		  volumetricCoeff = timeInterval*itNode->FastGetSolutionStepValue(POISSON_RATIO)*itNode->FastGetSolutionStepValue(YOUNG_MODULUS)/((1.0+itNode->FastGetSolutionStepValue(POISSON_RATIO))*(1.0-2.0*itNode->FastGetSolutionStepValue(POISSON_RATIO))) + 2.0*deviatoricCoeff/3.0;
		}
		else if(itNode->Is(FLUID)){
		  deviatoricCoeff = itNode->FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
		  if(deviatoricCoeff>1.0){
		    deviatoricCoeff=1.0;
		  }
		  volumetricCoeff = timeInterval*itNode->FastGetSolutionStepValue(BULK_MODULUS);
		}

		deltaPressure=itNode->FastGetSolutionStepValue(PRESSURE,0)-itNode->FastGetSolutionStepValue(PRESSURE,1);

		LHS_Contribution(0,0)+= nodalVolume/volumetricCoeff;
		RHS_Contribution[0]  += (-deltaPressure/volumetricCoeff +  itNode->GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE))*nodalVolume;

		const unsigned int xpos = itNode->GetDofPosition(PRESSURE);
		EquationId[0]=itNode->GetDof(PRESSURE,xpos).EquationId();

		stabilizationNeeded=false;
		if((itNode->Is(FLUID) || (itNode->Is(SOLID) && itNode->FastGetSolutionStepValue(POISSON_RATIO)>0.49))){
		  stabilizationNeeded=true;
		}else{

		  for (unsigned int i = 0; i< neighb_nodes.size(); i++)
		    {
		      EquationId[i+1]=neighb_nodes[i].GetDof(PRESSURE,xpos).EquationId();
		    }

		}
		if(stabilizationNeeded==true){
		  
		  firstRow=0;
		  firstCol=0;
		  meanMeshSize=itNode->FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
		  characteristicLength=2.0*meanMeshSize;
		  density=itNode->FastGetSolutionStepValue(DENSITY);
		
		  /* double tauStab=1.0/(8.0*deviatoricCoeff/(meanMeshSize*meanMeshSize)+2.0*density/timeInterval); */
		  
		  if(dimension==2){
		    nodalVelocity= sqrt(itNode->FastGetSolutionStepValue(VELOCITY_X)*itNode->FastGetSolutionStepValue(VELOCITY_X) +
					itNode->FastGetSolutionStepValue(VELOCITY_Y)*itNode->FastGetSolutionStepValue(VELOCITY_Y));
		  }else if(dimension==3){
		    nodalVelocity=sqrt(itNode->FastGetSolutionStepValue(VELOCITY_X)*itNode->FastGetSolutionStepValue(VELOCITY_X) +
				       itNode->FastGetSolutionStepValue(VELOCITY_Y)*itNode->FastGetSolutionStepValue(VELOCITY_Y) +
				       itNode->FastGetSolutionStepValue(VELOCITY_Z)*itNode->FastGetSolutionStepValue(VELOCITY_Z));
		  }
		  
		  tauStab= 1.0 * (characteristicLength * characteristicLength * timeInterval) / ( density * nodalVelocity * timeInterval * characteristicLength + density * characteristicLength * characteristicLength +  8.0 * deviatoricCoeff * timeInterval );

		  
		  /* std::cout<<"tauStab= "<<tauStab<<std::endl; */
		  
		  LHS_Contribution(0,0)+= +nodalVolume*tauStab*density/(volumetricCoeff*timeInterval);
		  RHS_Contribution[0]  += -nodalVolume*tauStab*density/(volumetricCoeff*timeInterval)*(deltaPressure-itNode->FastGetSolutionStepValue(PRESSURE_VELOCITY,0)*timeInterval);
		  
		  
		  if(itNode->Is(FREE_SURFACE)){
		    
		    /* LHS_Contribution(0,0) += + 2.0 * tauStab * nodalFreesurfaceArea / meanMeshSize; */
		    /* RHS_Contribution[0]   += - 2.0 * tauStab * nodalFreesurfaceArea / meanMeshSize * itNode->FastGetSolutionStepValue(PRESSURE,0); */

		    if(itNode->IsNot(RIGID)){
		      LHS_Contribution(0,0) += + 4.0 * 2.0 * tauStab * nodalVolume /(meanMeshSize*meanMeshSize);
		      RHS_Contribution[0]   += - 4.0 * 2.0 * tauStab * nodalVolume /(meanMeshSize*meanMeshSize) * itNode->FastGetSolutionStepValue(PRESSURE,0);
		    }
		    else{
		      LHS_Contribution(0,0) += + 4.0 * tauStab * nodalVolume /(meanMeshSize*meanMeshSize);
		      RHS_Contribution[0]   += - 4.0 * tauStab * nodalVolume /(meanMeshSize*meanMeshSize) * itNode->FastGetSolutionStepValue(PRESSURE,0);
		    }
		    const array_1d<double, 3> &Normal    = itNode->FastGetSolutionStepValue(NORMAL);
		    Vector& SpatialDefRate=itNode->FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
		    array_1d<double, 3> nodalAcceleration=  0.5*(itNode->FastGetSolutionStepValue(VELOCITY,0)-itNode->FastGetSolutionStepValue(VELOCITY,1))/timeInterval - itNode->FastGetSolutionStepValue(ACCELERATION,1);
		    /* nodalAcceleration=  (itNode->FastGetSolutionStepValue(VELOCITY,0)-itNode->FastGetSolutionStepValue(VELOCITY,1))/timeInterval; */

		    double nodalNormalAcceleration=0;
		    double nodalNormalProjDefRate=0;
		    if(dimension==2){
		      nodalNormalProjDefRate=Normal[0]*SpatialDefRate[0]*Normal[0] + Normal[1]*SpatialDefRate[1]*Normal[1] + 2*Normal[0]*SpatialDefRate[2]*Normal[1];
		      /* nodalNormalAcceleration=Normal[0]*itNode->FastGetSolutionStepValue(ACCELERATION_X,1) + Normal[1]*itNode->FastGetSolutionStepValue(ACCELERATION_Y,1); */
		      nodalNormalAcceleration=(0.5*(itNode->FastGetSolutionStepValue(VELOCITY_X,0)-itNode->FastGetSolutionStepValue(VELOCITY_X,1))/timeInterval+0.5*itNode->FastGetSolutionStepValue(ACCELERATION_X,1))*Normal[0] +
		      	(0.5*(itNode->FastGetSolutionStepValue(VELOCITY_Y,0)-itNode->FastGetSolutionStepValue(VELOCITY_Y,1))/timeInterval+0.5*itNode->FastGetSolutionStepValue(ACCELERATION_Y,1))*Normal[1];
		      nodalNormalAcceleration=Normal[0]*nodalAcceleration[0] + Normal[1]*nodalAcceleration[1];
		    }else if(dimension==3){
		      nodalNormalProjDefRate=Normal[0]*SpatialDefRate[0]*Normal[0] + Normal[1]*SpatialDefRate[1]*Normal[1] + Normal[2]*SpatialDefRate[2]*Normal[2] +
		  	2*Normal[0]*SpatialDefRate[3]*Normal[1] + 2*Normal[0]*SpatialDefRate[4]*Normal[2] + 2*Normal[1]*SpatialDefRate[5]*Normal[2];
		      
		      /* nodalNormalAcceleration=Normal[0]*itNode->FastGetSolutionStepValue(ACCELERATION_X) + Normal[1]*itNode->FastGetSolutionStepValue(ACCELERATION_Y) + Normal[2]*itNode->FastGetSolutionStepValue(ACCELERATION_Z); */
		      /* nodalNormalAcceleration=Normal[0]*nodalAcceleration[0] + Normal[1]*nodalAcceleration[1] + Normal[2]*nodalAcceleration[2]; */
		    }
		    /* RHS_Contribution[0]  += tauStab * (density*nodalNormalAcceleration - 4.0*deviatoricCoeff*nodalNormalProjDefRate/meanMeshSize) * nodalFreesurfaceArea; */
		    double accelerationContribution=2.0*density*nodalNormalAcceleration/meanMeshSize;
		    
		    double deviatoricContribution=8.0*deviatoricCoeff*nodalNormalProjDefRate/(meanMeshSize*meanMeshSize);

		    /* accelerationContribution=0; */
		    /* deviatoricContribution=0; */
		    if(itNode->IsNot(RIGID)){
		      RHS_Contribution[0]  += 2.0* tauStab * (accelerationContribution - deviatoricContribution) * nodalVolume;
		    }else{
		      RHS_Contribution[0]  += 1.0* tauStab * (accelerationContribution - deviatoricContribution) * nodalVolume;
		    }
		  }

		   
		  array_1d<double, 3 >& VolumeAcceleration = itNode->FastGetSolutionStepValue(VOLUME_ACCELERATION);


		  for (unsigned int i = 0; i< neighSize; i++)
		    {
		
		      dNdXi=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstCol];
		      dNdYi=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstCol+1];
		      
		      if(dimension==2){
		      	RHS_Contribution[i]  += - tauStab * density * (dNdXi* VolumeAcceleration[0] + dNdYi* VolumeAcceleration[1]) * nodalVolume;
		      }
		      else if(dimension==3){
		      	dNdZi=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstCol+2];
		      	RHS_Contribution[i]  += - tauStab * density * (dNdXi* VolumeAcceleration[0] + dNdYi* VolumeAcceleration[1] + dNdZi* VolumeAcceleration[2]) * nodalVolume;
		      }
		      			
		      if(i<neighb_nodes.size()){
		  	// i==0 of EquationIs has been already filled with the master node (that is not included in neighb_nodes)
		      	EquationId[i+1]=neighb_nodes[i].GetDof(PRESSURE,xpos).EquationId();
		  	// at i==0 density and volume acceleration are taken from the master node
		      	density= neighb_nodes[i].FastGetSolutionStepValue(DENSITY);
		      	VolumeAcceleration = neighb_nodes[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
		      }
		      		      
		      firstRow=0;

		      for (unsigned int j = 0; j< neighSize; j++)
		      	{
		      	  dNdXj=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstRow];
		      	  dNdYj=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstRow+1];
		      	  if(j!=0){
		      	    pressure=neighb_nodes[j-1].FastGetSolutionStepValue(PRESSURE,0);
		      	  }else{
		      	    pressure=itNode->FastGetSolutionStepValue(PRESSURE,0);
		      	  }
		      	  if(dimension==2){
		      	    ////////////////// Laplacian term for LHS
		      	    LHS_Contribution(i,j)+= + tauStab * (dNdXi*dNdXj + dNdYi*dNdYj) * nodalVolume;
		      	    ////////////////// Laplacian term L_ij*P_j for RHS
		      	    RHS_Contribution[i]  += - tauStab * (dNdXi*dNdXj + dNdYi*dNdYj)  * nodalVolume * pressure;
		      	  }
		      	  else if(dimension==3){
		      	    dNdZj=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstRow+2];
		      	    ////////////////// Laplacian term for LHS
		      	    LHS_Contribution(i,j) += + tauStab * (dNdXi*dNdXj + dNdYi*dNdYj + dNdZi*dNdZj) * nodalVolume;
		      	    ////////////////// Laplacian term L_ij*P_j for RHS
		      	    RHS_Contribution[i]   += - tauStab * (dNdXi*dNdXj + dNdYi*dNdYj + dNdZi*dNdZj) * nodalVolume * pressure;
		      	  }
		      	  firstRow+=dimension;
		      	}
		      
		      firstCol+=dimension;

		    }

	
		}

#ifdef _OPENMP
		Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
#else
		Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
#endif

	      }	      

	    }

	}


	KRATOS_CATCH("")

	  }

 
      /**
       * @brief This is a call to the linear system solver
       * @param A The LHS matrix
       * @param Dx The Unknowns vector
       * @param b The RHS vector
       */
      void SystemSolve(
		       TSystemMatrixType& A,
		       TSystemVectorType& Dx,
		       TSystemVectorType& b
		       ) override
      {
        KRATOS_TRY

	  double norm_b;
        if (TSparseSpace::Size(b) != 0)
	  norm_b = TSparseSpace::TwoNorm(b);
        else
	  norm_b = 0.00;

        if (norm_b != 0.00)
	  {
            //do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
	  }
        else
	  TSparseSpace::SetToZero(Dx);

        // Prints informations about the current time
        KRATOS_INFO_IF("NodalResidualBasedEliminationBuilderAndSolverContinuity", this->GetEchoLevel() > 1) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")

	  }

      /**
       *@brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
       * @param A The LHS matrix
       * @param Dx The Unknowns vector
       * @param b The RHS vector
       * @param rModelPart The model part of the problem to solve
       */
      void SystemSolveWithPhysics(
				  TSystemMatrixType& A,
				  TSystemVectorType& Dx,
				  TSystemVectorType& b,
				  ModelPart& rModelPart
				  )
      {
        KRATOS_TRY

	  double norm_b;
        if (TSparseSpace::Size(b) != 0)
	  norm_b = TSparseSpace::TwoNorm(b);
        else
	  norm_b = 0.00;

        if (norm_b != 0.00)
	  {
            //provide physical data as needed
            if(BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded() )
	      BaseType::mpLinearSystemSolver->ProvideAdditionalData(A, Dx, b, BaseType::mDofSet, rModelPart);

            //do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
	  }
        else
	  {
            TSparseSpace::SetToZero(Dx);
            KRATOS_WARNING_IF("NodalResidualBasedEliminationBuilderAndSolverContinuity", rModelPart.GetCommunicator().MyPID() == 0) << "ATTENTION! setting the RHS to zero!" << std::endl;
	  }

        // Prints informations about the current time
        KRATOS_INFO_IF("NodalResidualBasedEliminationBuilderAndSolverContinuity", this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")

	  }

      /**
       * @brief Function to perform the building and solving phase at the same time.
       * @details It is ideally the fastest and safer function to use when it is possible to solve
       * just after building
       * @param pScheme The integration scheme considered
       * @param rModelPart The model part of the problem to solve
       * @param A The LHS matrix
       * @param Dx The Unknowns vector
       * @param b The RHS vector
       */
      void BuildAndSolve(
			 typename TSchemeType::Pointer pScheme,
			 ModelPart& rModelPart,
			 TSystemMatrixType& A,
			 TSystemVectorType& Dx,
			 TSystemVectorType& b) override
      {
        KRATOS_TRY

	  Timer::Start("Build");

	/* boost::timer c_build_time; */
	BuildNodally(pScheme, rModelPart, A, b);
	/* std::cout << "CONTINUITY EQ: build_time : " << c_build_time.elapsed() << std::endl; */


        Timer::Stop("Build");

	//         ApplyPointLoads(pScheme,rModelPart,b);

        // Does nothing...dirichlet conditions are naturally dealt with in defining the residual
        ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        /* const double start_solve = OpenMPUtils::GetCurrentTime(); */
        Timer::Start("Solve");

	/* boost::timer c_solve_time; */
	SystemSolveWithPhysics(A, Dx, b, rModelPart);
	/* std::cout << "CONTINUITY EQ: solve_time : " << c_solve_time.elapsed() << std::endl; */

        Timer::Stop("Solve");
        /* const double stop_solve = OpenMPUtils::GetCurrentTime(); */

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        KRATOS_CATCH("")
	  }


      /**
       * @brief Builds the list of the DofSets involved in the problem by "asking" to each element
       * and condition its Dofs.
       * @details The list of dofs is stores insde the BuilderAndSolver as it is closely connected to the
       * way the matrix and RHS are built
       * @param pScheme The integration scheme considered
       * @param rModelPart The model part of the problem to solve
       */
      void SetUpDofSet(
		       typename TSchemeType::Pointer pScheme,
		       ModelPart& rModelPart
		       ) override
      {
        KRATOS_TRY;

        KRATOS_INFO_IF("NodalResidualBasedEliminationBuilderAndSolverContinuity", this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0) << "Setting up the dofs" << std::endl;

        //Gets the array of elements from the modeler
        ElementsArrayType& pElements = rModelPart.Elements();
        const int nelements = static_cast<int>(pElements.size());

        Element::DofsVectorType ElementalDofList;

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        unsigned int nthreads = OpenMPUtils::GetNumThreads();

	//         typedef boost::fast_pool_allocator< NodeType::DofType::Pointer > allocator_type;
	//         typedef std::unordered_set < NodeType::DofType::Pointer,
	//             DofPointerHasher,
	//             DofPointerComparor,
	//             allocator_type    >  set_type;

#ifdef USE_GOOGLE_HASH
        typedef google::dense_hash_set < NodeType::DofType::Pointer, DofPointerHasher>  set_type;
#else
        typedef std::unordered_set < NodeType::DofType::Pointer, DofPointerHasher>  set_type;
#endif
	//


        std::vector<set_type> dofs_aux_list(nthreads);
	//         std::vector<allocator_type> allocators(nthreads);

        for (int i = 0; i < static_cast<int>(nthreads); i++)
	  {
#ifdef USE_GOOGLE_HASH
            dofs_aux_list[i].set_empty_key(NodeType::DofType::Pointer());
#else
	    //             dofs_aux_list[i] = set_type( allocators[i]);
            dofs_aux_list[i].reserve(nelements);
#endif
	  }

#pragma omp parallel for firstprivate(nelements, ElementalDofList)
        for (int i = 0; i < static_cast<int>(nelements); i++)
	  {
            typename ElementsArrayType::iterator it = pElements.begin() + i;
            const unsigned int this_thread_id = OpenMPUtils::ThisThread();

            // gets list of Dof involved on every element
            pScheme->GetElementalDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);

            dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
	  }

        ConditionsArrayType& pConditions = rModelPart.Conditions();
        const int nconditions = static_cast<int>(pConditions.size());
#pragma omp parallel for firstprivate(nconditions, ElementalDofList)
        for (int i = 0; i < nconditions; i++)
	  {
            typename ConditionsArrayType::iterator it = pConditions.begin() + i;
            const unsigned int this_thread_id = OpenMPUtils::ThisThread();

            // gets list of Dof involved on every element
            pScheme->GetConditionDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);
            dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
	  }

        //here we do a reduction in a tree so to have everything on thread 0
        unsigned int old_max = nthreads;
        unsigned int new_max = ceil(0.5*static_cast<double>(old_max));
        while (new_max >= 1 && new_max != old_max)
	  {
	    //          //just for debugging
	    //          std::cout << "old_max" << old_max << " new_max:" << new_max << std::endl;
	    //          for (int i = 0; i < new_max; i++)
	    //          {
	    //             if (i + new_max < old_max)
	    //             {
	    //                std::cout << i << " - " << i + new_max << std::endl;
	    //             }
	    //          }
	    //          std::cout << "********************" << std::endl;

#pragma omp parallel for
            for (int i = 0; i < static_cast<int>(new_max); i++)
	      {
                if (i + new_max < old_max)
		  {
                    dofs_aux_list[i].insert(dofs_aux_list[i + new_max].begin(), dofs_aux_list[i + new_max].end());
                    dofs_aux_list[i + new_max].clear();
		  }
	      }

            old_max = new_max;
            new_max = ceil(0.5*static_cast<double>(old_max));
	  }

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        Doftemp.reserve(dofs_aux_list[0].size());
        for (auto it = dofs_aux_list[0].begin(); it != dofs_aux_list[0].end(); it++)
	  {
            Doftemp.push_back(it->get());
	  }
        Doftemp.Sort();

        BaseType::mDofSet = Doftemp;

        // Throws an execption if there are no Degrees of freedom involved in the analysis
        KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;

        BaseType::mDofSetIsInitialized = true;

        KRATOS_INFO_IF("NodalResidualBasedEliminationBuilderAndSolverContinuity", this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Finished setting up the dofs" << std::endl;

#ifdef _OPENMP
        if (mlock_array.size() != 0)
	  {
            for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
	      omp_destroy_lock(&mlock_array[i]);
	  }

        mlock_array.resize(BaseType::mDofSet.size());

        for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
	  omp_init_lock(&mlock_array[i]);
#endif

        // If reactions are to be calculated, we check if all the dofs have reactions defined
        // This is tobe done only in debug mode
#ifdef KRATOS_DEBUG
        if(BaseType::GetCalculateReactionsFlag())
	  {
            for(auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
	      {
		KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction()) << "Reaction variable not set for the following : " <<std::endl
								 << "Node : "<<dof_iterator->Id()<< std::endl
								 << "Dof : "<<(*dof_iterator)<<std::endl<<"Not possible to calculate reactions."<<std::endl;
	      }
	  }
#endif

        KRATOS_CATCH("");
      }

      /**
       * @brief Organises the dofset in order to speed up the building phase
       * @param rModelPart The model part of the problem to solve
       */
      void SetUpSystem(
		       ModelPart& rModelPart
		       ) override
      {
        // Set equation id for degrees of freedom
        // the free degrees of freedom are positioned at the beginning of the system,
        // while the fixed one are at the end (in opposite order).
        //
        // that means that if the EquationId is greater than "mEquationSystemSize"
        // the pointed degree of freedom is restrained
        //
        int free_id = 0;
        int fix_id = BaseType::mDofSet.size();

        for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
	  if (dof_iterator->IsFixed())
	    dof_iterator->SetEquationId(--fix_id);
	  else
	    dof_iterator->SetEquationId(free_id++);

        BaseType::mEquationSystemSize = fix_id;

      }

      //**************************************************************************
      //**************************************************************************

      void ResizeAndInitializeVectors(
				      typename TSchemeType::Pointer pScheme,
				      TSystemMatrixPointerType& pA,
				      TSystemVectorPointerType& pDx,
				      TSystemVectorPointerType& pb,
				      ModelPart& rModelPart
				      ) override
      {
        KRATOS_TRY

	  /* boost::timer c_contruct_matrix; */

	  if (pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
	    {
	      TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
	      pA.swap(pNewA);
	    }
        if (pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
	  {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
            pDx.swap(pNewDx);
	  }
        if (pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
	  {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
            pb.swap(pNewb);
	  }
        if (BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
	  {
            TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(0));
            BaseType::mpReactionsVector.swap(pNewReactionsVector);
	  }

        TSystemMatrixType& A = *pA;
        TSystemVectorType& Dx = *pDx;
        TSystemVectorType& b = *pb;

        //resizing the system vectors and matrix
        if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
	  {
            A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            ConstructMatrixStructure(pScheme, A, rModelPart);
	  }
        else
	  {
            if (A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
	      {
                /* KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW"); */
                KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permited."<<std::endl;
                A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
                ConstructMatrixStructure(pScheme, A, rModelPart);
	      }
	  }
        if (Dx.size() != BaseType::mEquationSystemSize)
	  Dx.resize(BaseType::mEquationSystemSize, false);
        if (b.size() != BaseType::mEquationSystemSize)
	  b.resize(BaseType::mEquationSystemSize, false);

        //if needed resize the vector for the calculation of reactions
        if (BaseType::mCalculateReactionsFlag == true)
	  {
            unsigned int ReactionsVectorSize = BaseType::mDofSet.size();
            if (BaseType::mpReactionsVector->size() != ReactionsVectorSize)
	      BaseType::mpReactionsVector->resize(ReactionsVectorSize, false);
	  }
	/* std::cout << "CONTINUITY EQ: contruct_matrix : " << c_contruct_matrix.elapsed() << std::endl; */
	
        KRATOS_CATCH("")

	  }

      //**************************************************************************
      //**************************************************************************


      /**
       * @brief Applies the dirichlet conditions. This operation may be very heavy or completely
       * unexpensive depending on the implementation choosen and on how the System Matrix is built.
       * @details For explanation of how it works for a particular implementation the user
       * should refer to the particular Builder And Solver choosen
       * @param pScheme The integration scheme considered
       * @param rModelPart The model part of the problem to solve
       * @param A The LHS matrix
       * @param Dx The Unknowns vector
       * @param b The RHS vector
       */
      void ApplyDirichletConditions(
				    typename TSchemeType::Pointer pScheme,
				    ModelPart& rModelPart,
				    TSystemMatrixType& A,
				    TSystemVectorType& Dx,
				    TSystemVectorType& b) override
      {
      }

      /**
       * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
       */
      void Clear() override
      {
        this->mDofSet = DofsArrayType();

        if (this->mpReactionsVector != NULL)
	  TSparseSpace::Clear((this->mpReactionsVector));
	//             this->mReactionsVector = TSystemVectorType();

        this->mpLinearSystemSolver->Clear();

        KRATOS_INFO_IF("NodalResidualBasedEliminationBuilderAndSolverContinuity", this->GetEchoLevel() > 1) << "Clear Function called" << std::endl;
      }

      /**
       * @brief This function is designed to be called once to perform all the checks needed
       * on the input provided. Checks can be "expensive" as the function is designed
       * to catch user's errors.
       * @param rModelPart The model part of the problem to solve
       * @return 0 all ok
       */
      int Check(ModelPart& rModelPart) override
      {
        KRATOS_TRY

	  return 0;
        KRATOS_CATCH("");
      }

      ///@}
      ///@name Access
      ///@{

      ///@}
      ///@name Inquiry
      ///@{

      ///@}
      ///@name Friends
      ///@{

      ///@}

    protected:
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

      void Assemble(
		    TSystemMatrixType& A,
		    TSystemVectorType& b,
		    const LocalSystemMatrixType& LHS_Contribution,
		    const LocalSystemVectorType& RHS_Contribution,
		    const Element::EquationIdVectorType& EquationId
#ifdef _OPENMP
		    ,std::vector< omp_lock_t >& lock_array
#endif
		    )
      {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
	  {
            unsigned int i_global = EquationId[i_local];

            if (i_global < BaseType::mEquationSystemSize)
	      {
#ifdef _OPENMP
                omp_set_lock(&lock_array[i_global]);
#endif
                b[i_global] += RHS_Contribution(i_local);
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
		  {
                    unsigned int j_global = EquationId[j_local];
                    if (j_global < BaseType::mEquationSystemSize)
		      {
                        A(i_global, j_global) += LHS_Contribution(i_local, j_local);
		      }
		  }
#ifdef _OPENMP
                omp_unset_lock(&lock_array[i_global]);
#endif

	      }
            //note that assembly on fixed rows is not performed here
	  }
      }


      //**************************************************************************


      virtual void ConstructMatrixStructure(
					    typename TSchemeType::Pointer pScheme,
					    TSystemMatrixType& A,
					    ModelPart& rModelPart)
      {
	//filling with zero the matrix (creating the structure)
	Timer::Start("MatrixStructure");

	ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

	// Getting the array of the conditions
	const int nconditions = static_cast<int>(rModelPart.Conditions().size());
	ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

	const std::size_t equation_size = BaseType::mEquationSystemSize;

#ifdef USE_GOOGLE_HASH
	std::vector<google::dense_hash_set<std::size_t> > indices(equation_size);
	const std::size_t empty_key = 2 * equation_size + 10;
#else
	std::vector<std::unordered_set<std::size_t> > indices(equation_size);
#endif

#pragma omp parallel for firstprivate(equation_size)
	for (int iii = 0; iii < static_cast<int>(equation_size); iii++)
	  {
#ifdef USE_GOOGLE_HASH
	    indices[iii].set_empty_key(empty_key);
#else
	    indices[iii].reserve(40);
#endif
	  }


	Element::EquationIdVectorType EquationId;

	ModelPart::NodeIterator NodesBegin;
	ModelPart::NodeIterator NodesEnd;
	OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

	for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
	  {
	    WeakPointerVector< Node < 3 > >& neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
	    const unsigned int neighSize =itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER).size();
	    if (EquationId.size() != neighSize)
	      EquationId.resize(neighSize, false);

	    /* const unsigned int xpos = itNode->GetDofPosition(VELOCITY_X); */
	    const unsigned int xpos = itNode->GetDofPosition(PRESSURE);
	    EquationId[0]=itNode->GetDof(PRESSURE,xpos).EquationId();

	    /* Vector& rNodeOrderedNeighboursPressureId=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_PRESSURE_ID); */
	  
	    
	    for (unsigned int i = 0; i< neighb_nodes.size(); i++)
	      {
		EquationId[i+1]=neighb_nodes[i].GetDof(PRESSURE,xpos).EquationId();
		/* unsigned int idNode=itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER)[i]; */
		/* EquationId[i]=rModelPart.Nodes()[idNode].GetDof(PRESSURE,xpos).EquationId(); */
		/* rNodeOrderedNeighboursPressureId[i]=EquationId[i]; */
	      }
	    


	    
	    for (std::size_t i = 0; i < EquationId.size(); i++)
	      {
		if (EquationId[i] < BaseType::mEquationSystemSize)
		  {
#ifdef _OPENMP
		    omp_set_lock(&mlock_array[EquationId[i]]);
#endif
		    auto& row_indices = indices[EquationId[i]];
		    for (auto it = EquationId.begin(); it != EquationId.end(); it++)
		      {
			if (*it < BaseType::mEquationSystemSize)
			  row_indices.insert(*it);
		      }
#ifdef _OPENMP
		    omp_unset_lock(&mlock_array[EquationId[i]]);
#endif
		  }
	      }
	

	  }

	Element::EquationIdVectorType ids(3, 0);


#pragma omp parallel for firstprivate(nconditions, ids)
	for (int iii = 0; iii<nconditions; iii++)
	  {
	    typename ConditionsArrayType::iterator i_condition = cond_begin + iii;
	    pScheme->Condition_EquationId( *(i_condition.base()) , ids, CurrentProcessInfo);
	    for (std::size_t i = 0; i < ids.size(); i++)
	      {
		if (ids[i] < BaseType::mEquationSystemSize)
		  {
#ifdef _OPENMP
		    omp_set_lock(&mlock_array[ids[i]]);
#endif
		    auto& row_indices = indices[ids[i]];
		    for (auto it = ids.begin(); it != ids.end(); it++)
		      {
			if (*it < BaseType::mEquationSystemSize)
			  row_indices.insert(*it);
		      }
#ifdef _OPENMP
		    omp_unset_lock(&mlock_array[ids[i]]);
#endif
		  }
	      }
	  }

	//count the row sizes
	unsigned int nnz = 0;
	for (unsigned int i = 0; i < indices.size(); i++)
	  nnz += indices[i].size();

	A = boost::numeric::ublas::compressed_matrix<double>(indices.size(), indices.size(), nnz);

	double* Avalues = A.value_data().begin();
	std::size_t* Arow_indices = A.index1_data().begin();
	std::size_t* Acol_indices = A.index2_data().begin();

	//filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
	Arow_indices[0] = 0;
	for (int i = 0; i < static_cast<int>(A.size1()); i++)
	  Arow_indices[i + 1] = Arow_indices[i] + indices[i].size();



#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(A.size1()); i++)
	  {
	    const unsigned int row_begin = Arow_indices[i];
	    const unsigned int row_end = Arow_indices[i + 1];
	    unsigned int k = row_begin;
	    for (auto it = indices[i].begin(); it != indices[i].end(); it++)
	      {
		Acol_indices[k] = *it;
		Avalues[k] = 0.0;
		k++;
	      }

	    std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

	  }

	A.set_filled(indices.size() + 1, nnz);

	Timer::Stop("MatrixStructure");

	/* std::cout<<"..... ConstructMatrixStructure for Momentum equation DONE"<<std::endl; */

      }




      void AssembleLHS(
		       TSystemMatrixType& A,
		       LocalSystemMatrixType& LHS_Contribution,
		       Element::EquationIdVectorType& EquationId
		       )
      {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
	  {
            unsigned int i_global = EquationId[i_local];
            if (i_global < BaseType::mEquationSystemSize)
	      {
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
		  {
                    unsigned int j_global = EquationId[j_local];
                    if (j_global < BaseType::mEquationSystemSize)
		      A(i_global, j_global) += LHS_Contribution(i_local, j_local);
		  }
	      }
	  }
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
#ifdef _OPENMP
      std::vector< omp_lock_t > mlock_array;
#endif
      ///@}
      ///@name Private Operators
      ///@{

      ///@}
      ///@name Private Operations
      ///@{

      inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate)
      {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();
        while (i != endit && (*i) != candidate)
	  {
            i++;
	  }
        if (i == endit)
	  {
            v.push_back(candidate);
	  }

      }

      void AssembleRHS(
		       TSystemVectorType& b,
		       const LocalSystemVectorType& RHS_Contribution,
		       const Element::EquationIdVectorType& EquationId
		       )
      {
        unsigned int local_size = RHS_Contribution.size();

        if (BaseType::mCalculateReactionsFlag == false)
	  {
            for (unsigned int i_local = 0; i_local < local_size; i_local++)
	      {
                const unsigned int i_global = EquationId[i_local];

                if (i_global < BaseType::mEquationSystemSize) //free dof
		  {
                    // ASSEMBLING THE SYSTEM VECTOR
                    double& b_value = b[i_global];
                    const double& rhs_value = RHS_Contribution[i_local];

#pragma omp atomic
                    b_value += rhs_value;
		  }
	      }
	  }
        else
	  {
            TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
            for (unsigned int i_local = 0; i_local < local_size; i_local++)
	      {
                const unsigned int i_global = EquationId[i_local];

                if (i_global < BaseType::mEquationSystemSize) //free dof
		  {
                    // ASSEMBLING THE SYSTEM VECTOR
                    double& b_value = b[i_global];
                    const double& rhs_value = RHS_Contribution[i_local];

#pragma omp atomic
                    b_value += rhs_value;
		  }
                else //fixed dof
		  {
                    double& b_value = ReactionsVector[i_global - BaseType::mEquationSystemSize];
                    const double& rhs_value = RHS_Contribution[i_local];

#pragma omp atomic
                    b_value += rhs_value;
		  }
	      }
	  }
      }

      //**************************************************************************

      void AssembleLHS_CompleteOnFreeRows(
					  TSystemMatrixType& A,
					  LocalSystemMatrixType& LHS_Contribution,
					  Element::EquationIdVectorType& EquationId
					  )
      {
        unsigned int local_size = LHS_Contribution.size1();
        for (unsigned int i_local = 0; i_local < local_size; i_local++)
	  {
            unsigned int i_global = EquationId[i_local];
            if (i_global < BaseType::mEquationSystemSize)
	      {
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
		  {
                    int j_global = EquationId[j_local];
                    A(i_global, j_global) += LHS_Contribution(i_local, j_local);
		  }
	      }
	  }
      }

      ///@}
      ///@name Private Operations
      ///@{

      ///@}
      ///@name Private  Access
      ///@{

      ///@}
      ///@name Private Inquiry
      ///@{

      ///@}
      ///@name Un accessible methods
      ///@{

      ///@}

    }; /* Class NodalResidualBasedEliminationBuilderAndSolverContinuity */

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}

} /* namespace Kratos.*/

#endif /* KRATOS_NODAL_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER  defined */
