//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti.
//

#if !defined(KRATOS_GLS_STRATEGY)
#define  KRATOS_GLS_STRATEGY

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */


#include "utilities/geometry_utilities.h"


#include "pfem_2_application.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"
#include "utilities/openmp_utils.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"

//#include "custom_utilities/solver_settings.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define QCOMP

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
    class FracStepStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    {
    public:
      /**@name Type Definitions */
      /*@{ */

      /** Counted pointer of ClassName */
      KRATOS_CLASS_POINTER_DEFINITION(  FracStepStrategy );

      typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

      typedef typename BaseType::TDataType TDataType;

      //typedef typename BaseType::DofSetType DofSetType;

      typedef typename BaseType::DofsArrayType DofsArrayType;

      typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

      typedef typename BaseType::TSystemVectorType TSystemVectorType;

      typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

      typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

      typedef OpenMPUtils::PartitionVector PartitionVector;


      /*@} */
      /**@name Life Cycle
       */
      /*@{ */

      /**
       * Constructor of the FracStepStrategy. Implements the solutions strategy for a Navier Stokes solver
       * using the fractional step approach. Prepared for both openmp parallelism and mpi parallelism. The function
       * also calls internally the "Check" function to verify that the input is complete
       * @param model_part - contains Nodes, elements, etc.
       * @param solver_config - auxiliary file to ease the configuration. Prescribes the linear solvers and builiding
       *        strategies to be used in defining the current composite solver.
       *        @see FractionalStepConfiguration for OpenMP setting or
       *        @see TrilinosFractionalStepConfiguration (in the Trilinos application) for the MPI version
       * @param ReformDofAtEachIteration - if set to true the graph of the matrix is recomputed at each iteration
       * @param velocity_toll - tolerance used in the velocity convergence check
       * @param pressure_toll - pressure tolerance in finalizing the predictor corrector strategy
       * @param MaxVelocityIterations - maximum number of iterations of the velocity solver
       * @param MaxPressureIterations - max number of iteration for the predictor corrector strategy
       * @param time_order - 1=BDF1 , 2=BDF2
       * @param domain_size 2=2D, 3=3D
       * @param predictor_corrector - true->for predictor corrector, false->standard Fractional Step (default = false)
       */
    FracStepStrategy( ModelPart& model_part, typename TLinearSolver::Pointer pNewVelocityLinearSolver,typename TLinearSolver::Pointer pNewPressureLinearSolver,
		      bool ReformDofAtEachIteration = true,
		      double velocity_toll = 0.01,
		      double pressure_toll = 0.01,
		      int MaxVelocityIterations = 3,
		      int MaxPressureIterations = 1,
		      unsigned int time_order = 2,
		      unsigned int domain_size = 2,
		      bool predictor_corrector = false
		      )
      : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, false)//, msolver_config(solver_config)
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
	  this->mecho_level = 1;


	  bool CalculateReactions = false;
	  bool CalculateNormDxFlag = true;
	  bool ReformDofAtEachIteration = false;

	  //computation of the fractional vel velocity (first step)
	  //3 dimensional case
	  //typedef typename Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3 > > > VarComponent;
	  typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
	  typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

	  //initializing fractional velocity solution step
	  typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
	  typename SchemeType::Pointer pscheme = typename SchemeType::Pointer(new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());

	  BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (pNewVelocityLinearSolver));


	  this->mpfracvel_strategy = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (model_part, pscheme, pNewVelocityLinearSolver, vel_build, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag));

	  this->mpfracvel_strategy->SetEchoLevel(1);


	  BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pNewPressureLinearSolver, PRESSURE));

	  this->mppressurestep = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (model_part, pscheme,pNewPressureLinearSolver, pressure_build, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag));
	  this->mppressurestep->SetEchoLevel(2);

	  this->m_step = 1;

	  mHasSlipProcess = false;

	  KRATOS_CATCH("")
	    }

      /** Destructor.
       */
      virtual ~FracStepStrategy()
	{
	}

      /** Destructor.
       */

      double Solve() override
      {
	KRATOS_TRY
	  Timer time;
	Timer::Start("Solve_strategy");

#if defined(QCOMP)
	double Dp_norm;
	Dp_norm = IterativeSolve();
#else
	//multifluids
	AssignInitialStepValues();

	double Dp_norm = 1.00;
	//int iteration = 0;
	//int MaxPressureIterations = this->mMaxPressIterations;
	Dp_norm = IterativeSolve();
#endif
	//this->Clear();
	this->m_step += 1;
	return Dp_norm;
	KRATOS_CATCH("")
	  }

      double SolvePressure()
      {
	KRATOS_TRY

	//ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
	this->SolveStep7(); //pold=pn+1
	double Dp_norm = this->SolveStep2();
      	return Dp_norm;

	KRATOS_CATCH("")
	  }

      double IterativeSolve()
      {
	KRATOS_TRY
	  Timer time;
	Timer::Start("Solve_ambos");

	double Dp_norm = 1.00;
	ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
	//KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
	rCurrentProcessInfo[VISCOSITY] = 1.0;

	for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();i != BaseType::GetModelPart().NodesEnd(); ++i)
	  {
	    i->FastGetSolutionStepValue(VELOCITY_X,1) =  i->FastGetSolutionStepValue(VELOCITY_X);
	    i->FastGetSolutionStepValue(VELOCITY_Y,1) =  i->FastGetSolutionStepValue(VELOCITY_Y);
	    i->FastGetSolutionStepValue(VELOCITY_Z,1) =  i->FastGetSolutionStepValue(VELOCITY_Z);
	  }

#if defined(QCOMP)
	this->SolveStep1(this->mvelocity_toll, this->mMaxVelIterations);
#else
	this->SolveStep3();
#endif
        //double p_norm=0.0;
#if defined(QCOMP)

	//polimero
	this->SolveStepaux();
	//int MaxPressureIterations = this->mMaxPressIterations;
	//int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
	//double p_norm = SavePressureIteration();
	Dp_norm = 1.0;
	//Timer::Stop("Solve_ambos");
	//KRATOS_WATCH(time)
#else
	int iteration = 0;
	while (  iteration++ < 3)
	  {
	    Dp_norm = SolvePressure();
 	    double p_norm = SavePressureIteration();
	    if (fabs(p_norm) > 1e-10){
	      Dp_norm /= p_norm;
	    }
	    else
	      Dp_norm = 1.0;
	    this->SolveStep4();
	  }
#endif
	this->Clear();
	return Dp_norm;
	KRATOS_CATCH("")
	  }

      /**
       * copies PRESSURE->PRESSURE_OLD_IT
       * @return the norm of the pressure vector
       */
      double SavePressureIteration()
      {
        KRATOS_TRY

	  double local_p_norm = 0.0;
        for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
	     i != BaseType::GetModelPart().NodesEnd(); ++i)
	  {
            //setting the old value of the pressure to the current one
            const double& p = (i)->FastGetSolutionStepValue(PRESSURE);
	    local_p_norm += p*p;
	  }

        double p_norm = local_p_norm;

        //TODO: prepare for parallelization
        p_norm = sqrt(p_norm);

        return p_norm;
        KRATOS_CATCH("")
	  }



      void AssignInitialStepValues()
      {
	KRATOS_TRY

	  ModelPart& model_part=BaseType::GetModelPart();

	const double dt = model_part.GetProcessInfo()[DELTA_TIME];

	for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();i != BaseType::GetModelPart().NodesEnd(); ++i)
	  {
	    (i)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = 0.0;
	    (i)->FastGetSolutionStepValue(PRESSURE) = 0.0;
	    (i)->FastGetSolutionStepValue(PRESSURE,1) = 0.0;
	  }
	KRATOS_CATCH("");
      }


      /**
       * this function performs the iterative solution of the non-linear velocity problem in the first step
     * of the fractional step procedure
     * @param velocity_toll - tolerance used in the velocity convergence check
     * @param MaxIterations - max number of iterations
     */
      void SolveStep1(double velocity_toll, int MaxIterations)
      {
	KRATOS_TRY;

	Timer time;
	Timer::Start("SolveStep1");

	int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        double normDx = 0.0;

        bool is_converged = false;
	int iteration = 0;
        //double iteration = 1;
	//ModelPart& model_part=BaseType::GetModelPart();

	while (is_converged == false && iteration++<3)
	  {
	    //perform one iteration over the fractional step velocity
	    normDx = FractionalVelocityIteration();
            is_converged = ConvergenceCheck(normDx, velocity_toll);
	  }
        if (is_converged == false)
	  if (rank == 0) std::cout << "ATTENTION: convergence NOT achieved" << std::endl;

	KRATOS_CATCH("");

      }

      double FractionalVelocityIteration()
      {
	KRATOS_TRY
	  ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
	rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
	double normDx = mpfracvel_strategy->Solve();
	return normDx;
	KRATOS_CATCH("");
      }

      void SolveStep4()
      {
        KRATOS_TRY;
        Timer time;
        Timer::Start("paso_4");

        array_1d<double, 3 > zero = ZeroVector(3);

//#ifdef _OPENMP
//        int number_of_threads = omp_get_max_threads();
//#else
//        int number_of_threads = 1;
//#endif

        //ModelPart& model_part=BaseType::GetModelPart();

        //double dt = model_part.GetProcessInfo()[DELTA_TIME];
	//dt=0.005;
	for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();i != BaseType::GetModelPart().NodesEnd(); ++i)
	  {
	    array_1d<double, 3 > zero = ZeroVector(3);
	    i->FastGetSolutionStepValue(FORCE)=ZeroVector(3);
	    double & nodal_mass = (i)->FastGetSolutionStepValue(NODAL_MASS);
	    nodal_mass = 0.0;
	  }

 	ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
	rCurrentProcessInfo[FRACTIONAL_STEP] = 6;
	for (ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin(); i != BaseType::GetModelPart().ElementsEnd(); ++i)
	  {
            (i)->InitializeSolutionStep(BaseType::GetModelPart().GetProcessInfo());
	  }

	for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();i != BaseType::GetModelPart().NodesEnd(); ++i)
	  {
	    array_1d<double,3>& force_temp = i->FastGetSolutionStepValue(FORCE);
	    double A = (i)->FastGetSolutionStepValue(NODAL_MASS);
	    if(A<0.0000000000000001){
	      A=1.0;
	    }

	    double dt_Minv = 0.005  / A ;
	    //dt_Minv=1.0;
	    force_temp *= dt_Minv;

	  //KRATOS_WATCH(force_temp);
	    if(!i->IsFixed(VELOCITY_X)) //FRACT_VEL_X
	      {
		i->FastGetSolutionStepValue(VELOCITY_X) += force_temp[0] ;
	      }
	    if(!i->IsFixed(VELOCITY_Y))
	      {
		i->FastGetSolutionStepValue(VELOCITY_Y) +=force_temp[1];
	      }
	    if(!i->IsFixed(VELOCITY_Z))
	      {
		i->FastGetSolutionStepValue(VELOCITY_Z) +=force_temp[2] ;
	      }
	    if(i->IsFixed(VELOCITY_X))
	      {
		i->FastGetSolutionStepValue(VELOCITY_X)=0.0; //i->FastGetSolutionStepValue(VELOCITY_X,1);
	      }
	    if(i->IsFixed(VELOCITY_Y))
	      {
		i->FastGetSolutionStepValue(VELOCITY_Y)= 0.0; //i->FastGetSolutionStepValue(VELOCITY_Y,1) ;
	      }
	    if(i->IsFixed(VELOCITY_Z))
	      {
		i->FastGetSolutionStepValue(VELOCITY_Z)=0.0; //i->FastGetSolutionStepValue(VELOCITY_Z,1) ;
	      }
	  }
	KRATOS_CATCH("");
    }

      void SolveStep7()
      {
	KRATOS_TRY;
	//  ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
	array_1d<double, 3 > zero = ZeroVector(3);
	//Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

#ifdef _OPENMP
	int number_of_threads = omp_get_max_threads();
#else
	int number_of_threads = 1;
#endif

	//ModelPart& model_part=BaseType::GetModelPart();
	//const double dt = model_part.GetProcessInfo()[DELTA_TIME];
	vector<unsigned int> partition;
	CreatePartition(number_of_threads, BaseType::GetModelPart().Nodes().size(), partition);

#pragma omp parallel for schedule(static,1)
      for (int k = 0; k < number_of_threads; k++)
        {
            ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
            ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k + 1];

            array_1d<double, 3 > zero = ZeroVector(3);
            for (typename ModelPart::NodesContainerType::iterator it=it_begin; it!=it_end; ++it)
	      {
		it->FastGetSolutionStepValue(PRESSURE_OLD_IT)=it->FastGetSolutionStepValue(PRESSURE);
	      }
        }
      KRATOS_CATCH("");
      }

      double SolveStep2()
      {
	KRATOS_TRY;
	Timer::Start("Presion");
	BaseType::GetModelPart().GetProcessInfo()[FRACTIONAL_STEP] = 4;
	return mppressurestep->Solve();
	Timer::Stop("Presion");
	//KRATOS_WATCH(*time)
	//mppressurestep->Clear();
	KRATOS_CATCH("");

    }

      void SolveStep3()
      {
	KRATOS_TRY

	  ModelPart& model_part=BaseType::GetModelPart();
	const double dt = model_part.GetProcessInfo()[DELTA_TIME];
        Timer time;
        Timer::Start("paso_3");

#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        vector<unsigned int> partition;
        CreatePartition(number_of_threads, BaseType::GetModelPart().Nodes().size(), partition);



#pragma omp parallel for schedule(static,1)
        for (int k = 0; k < number_of_threads; k++)
	  {
            ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
            ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k + 1];


            array_1d<double, 3 > zero = ZeroVector(3);


            for (ModelPart::NodeIterator i = it_begin; i != it_end; ++i)
	      {
                double & nodal_mass = (i)->FastGetSolutionStepValue(NODAL_MASS);
                nodal_mass = 0.0;
	        noalias(i->FastGetSolutionStepValue(FORCE)) =    zero;
                //double & nodal_area = (i)->FastGetSolutionStepValue(NODAL_AREA);
                //nodal_area = 0.0;
	      }
	  }
        array_1d<double,3> zero = ZeroVector(3);

        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        rCurrentProcessInfo[FRACTIONAL_STEP] = 5;

        vector<unsigned int> elem_partition;
        CreatePartition(number_of_threads, BaseType::GetModelPart().Elements().size(), elem_partition);

#pragma omp parallel for schedule(static,1)
        for (int k = 0; k < number_of_threads; k++)
	  {
            ModelPart::ElementIterator it_begin = BaseType::GetModelPart().ElementsBegin() + elem_partition[k];
            ModelPart::ElementIterator it_end = BaseType::GetModelPart().ElementsBegin() + elem_partition[k + 1];
            for (ModelPart::ElementIterator i = it_begin; i != it_end; ++i)
	      {
                (i)->InitializeSolutionStep(BaseType::GetModelPart().GetProcessInfo());
	      }
	  }
#pragma omp parallel for schedule(static,1)
        for (int k = 0; k < number_of_threads; k++)
	  {
            ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
            ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k + 1];

            for (ModelPart::NodeIterator i = it_begin; i != it_end; ++i)
	      {
                array_1d<double,3>& force_temp = i->FastGetSolutionStepValue(FORCE);

                force_temp *=(1.0/ i->FastGetSolutionStepValue(NODAL_MASS));

		//array_1d<double,3>& vel = i->FastGetSolutionStepValue(VELOCITY);
 		i->FastGetSolutionStepValue(VELOCITY) = i->FastGetSolutionStepValue(VELOCITY,1) + dt * force_temp;
	      }
	  }

        KRATOS_CATCH("");
      }

      void SolveStepaux()
      {
	KRATOS_TRY

	  Timer time;
	Timer::Start("SolveStepaux");

	//ModelPart& model_part=BaseType::GetModelPart();
	//const double dt = model_part.GetProcessInfo()[DELTA_TIME];

#ifdef _OPENMP
	int number_of_threads = omp_get_max_threads();
#else
	int number_of_threads = 1;
#endif

	//number_of_threads = 1;
	vector<unsigned int> partition;
	CreatePartition(number_of_threads, BaseType::GetModelPart().Nodes().size(), partition);

#pragma omp parallel for schedule(static,1)
	for (int k = 0; k < number_of_threads; k++)
	  {
            ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
            ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k + 1];

            array_1d<double, 3 > zero = ZeroVector(3);

            for (ModelPart::NodeIterator i = it_begin; i != it_end; ++i)
	      {
                double & nodal_mass = (i)->FastGetSolutionStepValue(NODAL_MASS);
                nodal_mass = 0.0;
		i->FastGetSolutionStepValue(PRESSUREAUX)=0.0;
                i->FastGetSolutionStepValue(PRESSURE)=0.0;
	      }
	  }


	ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
	rCurrentProcessInfo[FRACTIONAL_STEP] = 7;

	vector<unsigned int> elem_partition;
	CreatePartition(number_of_threads, BaseType::GetModelPart().Elements().size(), elem_partition);
#pragma omp parallel for schedule(static,1)
	for (int k = 0; k < number_of_threads; k++)
	  {
	    ModelPart::ElementIterator it_begin = BaseType::GetModelPart().ElementsBegin() + elem_partition[k];
	    ModelPart::ElementIterator it_end = BaseType::GetModelPart().ElementsBegin() + elem_partition[k + 1];
            for (ModelPart::ElementIterator i = it_begin; i != it_end; ++i)
	      {
                (i)->InitializeSolutionStep(BaseType::GetModelPart().GetProcessInfo());
	      }
	  }


#pragma omp parallel for schedule(static,1)
	for (int k = 0; k < number_of_threads; k++)
	  {
	    ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
	    ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k + 1];

            for (ModelPart::NodeIterator i = it_begin; i != it_end; ++i)
	      {
		if(i->FastGetSolutionStepValue(NODAL_MASS)==0.0)
		  {
		    i->FastGetSolutionStepValue(PRESSURE)=0.0;
		  }
		else
		  {
		    //if()
		    i->FastGetSolutionStepValue(PRESSURE)=i->FastGetSolutionStepValue(PRESSUREAUX) * (1.0/ i->FastGetSolutionStepValue(NODAL_MASS));
		  }

	      }
	  }


        KRATOS_CATCH("");


      }


    /**
     * implements the convergence check for the velocities
     * convergence is considered achieved when normDx/norm(v) is less than tol
     * @param normDx norm of the VELOCITY correction
     * @param toll tolerance accepted
     * @return true if converged
     */
      bool ConvergenceCheck(const double& normDx, double tol)
      {
	KRATOS_TRY;
	double norm_v = 0.00;


	for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
	     i != BaseType::GetModelPart().NodesEnd(); ++i)
	  {
	    const array_1d<double, 3 > & v = (i)->FastGetSolutionStepValue(VELOCITY);

	    norm_v += v[0] * v[0];
	    norm_v += v[1] * v[1];
	    norm_v += v[2] * v[2];
	  }

	//BaseType::GetModelPart().GetCommunicator().SumAll(norm_v);

        double norm_v1 = sqrt(norm_v);

        if (norm_v1 == 0.0) norm_v1 = 1.00;

        double ratio = normDx / norm_v1;

        int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
        if (rank == 0) std::cout << "velocity ratio = " << ratio << std::endl;


        if (ratio < tol)
	  {
            if (rank == 0) std::cout << "convergence achieved" << std::endl;
            return true;
	  }

        return false;


        KRATOS_CATCH("");
      }


      void Compute()
      {
        KRATOS_TRY

#ifdef _OPENMP
	  int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        array_1d<double,3> aux;
        array_1d<double,3> aux1;

        ModelPart& model_part=BaseType::GetModelPart();
        const double dt = model_part.GetProcessInfo()[DELTA_TIME];


        vector<unsigned int> partition;
        CreatePartition(number_of_threads, BaseType::GetModelPart().Nodes().size(), partition);



#pragma omp parallel for schedule(static,1)
        for (int k = 0; k < number_of_threads; k++)
	  {
            ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
            ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k + 1];


            array_1d<double, 3 > zero = ZeroVector(3);

            //first of all set to zero the nodal variables to be updated nodally
            for (ModelPart::NodeIterator i = it_begin; i != it_end; ++i)
	      {
                noalias(i->FastGetSolutionStepValue(FORCE)) =    zero;
                (i)->FastGetSolutionStepValue(NODAL_MASS)=0.0;
	      }
	  }

        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        rCurrentProcessInfo[FRACTIONAL_STEP] = 6;


        vector<unsigned int> elem_partition;
        CreatePartition(number_of_threads, BaseType::GetModelPart().Elements().size(), elem_partition);

#pragma omp parallel for schedule(static,1)
        for (int k = 0; k < number_of_threads; k++)
	  {

            ModelPart::ElementIterator it_begin = BaseType::GetModelPart().ElementsBegin() + elem_partition[k];
            ModelPart::ElementIterator it_end = BaseType::GetModelPart().ElementsBegin() + elem_partition[k + 1];


            for (ModelPart::ElementIterator i = it_begin; i != it_end; ++i)
	      {

                (i)->InitializeSolutionStep(BaseType::GetModelPart().GetProcessInfo());
	      }

	  }


#pragma omp parallel for schedule(static,1)
        for (int k = 0; k < number_of_threads; k++)
	  {
            ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
            ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k + 1];


            array_1d<double, 3 > zero = ZeroVector(3);

            //first of all set to zero the nodal variables to be updated nodally
            for (ModelPart::NodeIterator i = it_begin; i != it_end; ++i)
	      {
                array_1d<double,3>& force_temp = i->FastGetSolutionStepValue(FORCE);
                force_temp -=i->FastGetSolutionStepValue(NODAL_MASS) * ((i)->FastGetSolutionStepValue(VELOCITY)-(i)->FastGetSolutionStepValue(VELOCITY,1))/dt;
	      }
	  }
        KRATOS_CATCH("");
      }


      /**
       *
       * @param Level
       */
      virtual void SetEchoLevel(int Level) override
      {
        mecho_level = Level;
        mpfracvel_strategy->SetEchoLevel(Level);
        mppressurestep->SetEchoLevel(Level);
      }


      virtual void Clear() override
      {
        int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
        if (rank == 0) KRATOS_WATCH("FracStepStrategy Clear Function called");
        mpfracvel_strategy->Clear();
        mppressurestep->Clear();
      }

      virtual double GetStageResidualNorm(unsigned int step)
      {
        if (step <= 3)
	  return mpfracvel_strategy->GetResidualNorm();
        if (step == 4)
	  return mppressurestep->GetResidualNorm();
        else
	  return 0.0;
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
      typename BaseType::Pointer mpfracvel_strategy;
      typename BaseType::Pointer mppressurestep;

      double mvelocity_toll;
      double mpressure_toll;
      int mMaxVelIterations;
      int mMaxPressIterations;
      unsigned int mtime_order;
      unsigned int mprediction_order;
      bool mpredictor_corrector;
      bool mReformDofAtEachIteration;
      int mecho_level;

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
      bool proj_is_initialized;

      //GenerateSlipConditionProcess::Pointer mpSlipProcess;
      bool mHasSlipProcess;

      std::vector< Process::Pointer > mInitializeIterationProcesses;

      inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
      {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (unsigned int i = 1; i < number_of_threads; i++)
	  partitions[i] = partitions[i - 1] + partition_size;
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
      FracStepStrategy(const FracStepStrategy& Other);


      /*@} */

    }; /* Class FracStepStrategy */

  /*@} */

  /**@name Type Definitions */
  /*@{ */


  /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_FRACTIONALSTEP_STRATEGY  defined */
