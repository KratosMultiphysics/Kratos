/* *********************************************************   
          
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
//#include "incompressible_fluid_application.h"
#include "custom_strategies/strategies/solver_configuration.h"
#include "utilities/geometry_utilities.h"
#include "custom_processes/generate_slip_condition_process.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "incompressible_fluid_application.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

#ifdef _OPENMP
#include "omp.h"
#endif






//#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
//#include "custom_strategies/builder_and_solvers/residualbased_elimination_discretelaplacian_builder_and_solver.h"
//#include "custom_strategies/builder_and_solvers/residualbased_elimination_discretelaplacian_builder_and_solver_flexiblefsi.h"





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
        typedef boost::shared_ptr< FracStepStrategy<TSparseSpace, TDenseSpace, TLinearSolver> > Pointer;

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
        FracStepStrategy(
                ModelPart& model_part,
                SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& solver_config,
                bool ReformDofAtEachIteration = true,
                double velocity_toll = 0.01,
                double pressure_toll = 0.01,
                int MaxVelocityIterations = 3,
                int MaxPressureIterations = 1,
                unsigned int time_order = 2,
                unsigned int domain_size = 2,
                bool predictor_corrector = false
                )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, false), msolver_config(solver_config)
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

            //performs checks to verify the quality of the input
            //this->Check();

            //initialize strategy
            this->mpfracvel_strategy = solver_config.pGetStrategy(std::string("vel_strategy"));
            this->mppressurestep = solver_config.pGetStrategy(std::string("pressure_strategy"));

            //fix fractional_velocities as needed
            for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                    i != BaseType::GetModelPart().NodesEnd(); ++i)
            {
                if (i->IsFixed(VELOCITY_X))
                    (i)->Fix(FRACT_VEL_X);
                if (i->IsFixed(VELOCITY_Y))
                    (i)->Fix(FRACT_VEL_Y);
                if (i->IsFixed(VELOCITY_Z))
                    (i)->Fix(FRACT_VEL_Z);
            }

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

        //*********************************************************************************
        //**********************************************************************
       double Solve()
        {
            KRATOS_TRY
	    Timer time;
	    Timer::Start("FractionalStep");

            double Dp_norm;

	    Dp_norm = FracStepSolution();

            if (this->mReformDofAtEachIteration == true)
                this->Clear();

            this->m_step += 1;
	    Timer::Stop("FractionalStep");
	    KRATOS_WATCH(time)
            return Dp_norm;

            KRATOS_CATCH("")
        }


        //*********************************************************************************
         //******************************************************************************************************
        double FracStepSolution()
        {
            KRATOS_TRY

            double Dp_norm = this->SolveStep2();

            return Dp_norm;

            KRATOS_CATCH("")
        }

        //******************************************************************************************************
        //******************************************************************************************************


        void SolveStep4()
        {
	  KRATOS_TRY;
	  Timer time;
	  Timer::Start("paso_4");
	//  ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
	  array_1d<double, 3 > zero = ZeroVector(3);
	  //Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
	  
#ifdef _OPENMP
	  int number_of_threads = omp_get_max_threads();
#else
	  int number_of_threads = 1;
#endif
	    
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
		for (typename ModelPart::NodesContainerType::iterator it=it_begin; it!=it_end; ++it)
		  {
		    it->FastGetSolutionStepValue(FORCE)=ZeroVector(3);
		    array_1d<double, 3 > & press_proj = (it)->FastGetSolutionStepValue(PRESS_PROJ);
                    noalias(press_proj) = ZeroVector(3);
		  }		
	      }
	 
	  
	  
	  vector<unsigned int> elem_partition;
	  CreatePartition(number_of_threads, BaseType::GetModelPart().Elements().size(), elem_partition);
	  
#pragma omp parallel for schedule(static,1)
	  for (int k = 0; k < number_of_threads; k++)
	    {
	      ModelPart::ElementIterator it_begin = BaseType::GetModelPart().ElementsBegin() + elem_partition[k];
	      ModelPart::ElementIterator it_end = BaseType::GetModelPart().ElementsBegin() + elem_partition[k + 1];
	      for (ModelPart::ElementIterator i = it_begin; i != it_end; ++i)
		
		{
		  //get the list of nodes of the element
		  boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
		  array_1d<double,3> N;
		  Geometry< Node<3> >& geom = i->GetGeometry();
		  array_1d<double, 2 > vel_gauss;
		  double volume;
		  GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);			
		  
		  array_1d<double,3> aux;
		  double p0 = geom[0].FastGetSolutionStepValue(PRESSURE);
		  double p0old = geom[0].FastGetSolutionStepValue(PRESSURE,1);
		  
		  double p1 = geom[1].FastGetSolutionStepValue(PRESSURE);
		  double p1old = geom[1].FastGetSolutionStepValue(PRESSURE,1);
		  
		  double p2 = geom[2].FastGetSolutionStepValue(PRESSURE);
		  double p2old = geom[2].FastGetSolutionStepValue(PRESSURE,1);
		    
		  double p_avg = N[0]*(p0 - p0old) + N[1]*(p1 - p1old) + N[2]*(p2 - p2old);
		  p_avg *= volume;
		  aux [0] = DN_DX(0, 0) * p_avg;
		  aux [1] = DN_DX(0, 1) * p_avg;
		  aux [2] =0.0;
		  
		  
		  geom[0].SetLock();
		  geom[0].FastGetSolutionStepValue(FORCE) += aux;
		  geom[0].UnSetLock();


		  aux[0] = DN_DX(1, 0) * p_avg;
		  aux[1] = DN_DX(1, 1) * p_avg;
		  
		  geom[1].SetLock();
		  geom[1].FastGetSolutionStepValue(FORCE) += aux;
		  geom[1].UnSetLock();
		  
		  aux[0] = DN_DX(2, 0) * p_avg;
		  aux[1] = DN_DX(2, 1) * p_avg;
		  
		  geom[2].SetLock();
		  geom[2].FastGetSolutionStepValue(FORCE) += aux;
		  geom[2].UnSetLock();

		  ////////////
                  array_1d<double, 3 > & press_proj0 = geom[0].FastGetSolutionStepValue(PRESS_PROJ);

                  array_1d<double, 3 > & press_proj1 = geom[1].FastGetSolutionStepValue(PRESS_PROJ);

                  array_1d<double, 3 > & press_proj2 = geom[2].FastGetSolutionStepValue(PRESS_PROJ);


		    //calculation of the pressure gradient (saved in vel_gauss)
		    //note that here we calculate it "strong"
		    vel_gauss[0] = DN_DX(0, 0)*(p0) + DN_DX(1, 0)*(p1) + DN_DX(2, 0)*(p2);
		    vel_gauss[1] = DN_DX(0, 1)*(p0) + DN_DX(1, 1)*(p1) + DN_DX(2, 1)*(p2);
		    vel_gauss *= volume;

		    //press_proj += G*p
		    geom[0].SetLock();
		    press_proj0[0] += N[0] * vel_gauss[0];
		    press_proj0[1] += N[0] * vel_gauss[1];
		    geom[0].UnSetLock();

		    geom[1].SetLock();
		    press_proj1[0] += N[1] * vel_gauss[0];
		    press_proj1[1] += N[1] * vel_gauss[1];
		    geom[1].UnSetLock();

		    geom[2].SetLock();
		    press_proj2[0] += N[2] * vel_gauss[0];
		    press_proj2[1] += N[2] * vel_gauss[1];
		    geom[2].UnSetLock();
		    ////////////////

		}
	    }
	  //correct the velocities
	  
#pragma omp parallel for schedule(static,1)
	  for (int k = 0; k < number_of_threads; k++)
	    {
	      ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
	      ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k + 1];
	      
	    for (typename ModelPart::NodesContainerType::iterator it=it_begin; it!=it_end; ++it)
	      {
		double dt_Minv = (dt / 2.00) / it->FastGetSolutionStepValue(NODAL_MASS);
		array_1d<double,3>& force_temp = it->FastGetSolutionStepValue(FORCE);
		force_temp *= dt_Minv;

		array_1d<double, 3 > & press_proj = (it)->FastGetSolutionStepValue(PRESS_PROJ);
                //double A = (it)->FastGetSolutionStepValue(NODAL_MASS);
		double A = (it)->FastGetSolutionStepValue(NODAL_AREA);

		//double density_inverse = 1.0 / it->FastGetSolutionStepValue(DENSITY);

                double temp = 1.00 / A;
                press_proj *= temp;


		//press_proj *= density_inverse ;
		


		if(!it->IsFixed(VELOCITY_X))
		  {
		    it->FastGetSolutionStepValue(VELOCITY_X)+=force_temp[0] ; //* density_inverse;
		  }
		if(!it->IsFixed(VELOCITY_Y))
		  {
		    it->FastGetSolutionStepValue(VELOCITY_Y)+=force_temp[1] ; //* density_inverse;						
		  }
		if(!it->IsFixed(VELOCITY_Z))
		  {
		    it->FastGetSolutionStepValue(VELOCITY_Z)+=force_temp[2] ; //* density_inverse;						
		  }
	      }
	    }
	  
	  Timer::Stop("paso_4");
	  KRATOS_WATCH(time)
	    
	    KRATOS_CATCH("");
        }
	
	//******************************************************************************************************
        //******************************************************************************************************
        /**
         * solution of the pressure. Implements the second step of the fractional step
         * @return norm of the pressure variation vector
         */
        double SolveStep2()
        {
	  KRATOS_TRY;
	  Timer::Start("Presion");
	  BaseType::GetModelPart().GetProcessInfo()[FRACTIONAL_STEP] = 4;
	  return mppressurestep->Solve();
	  Timer::Stop("Presion");
	  //KRATOS_WATCH(*time)
	    
	    KRATOS_CATCH("");

        }
	
        //******************************************************************************************************
        //******************************************************************************************************
	
        /**
         * calculation of the projections. Needed for OSS
         */
        void SolveStep3()
        {
            KRATOS_TRY
	      
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
		
                //first of all set to zero the nodal variables to be updated nodally
                for (ModelPart::NodeIterator i = it_begin; i != it_end; ++i)
		  {
		    double & nodal_mass = (i)->FastGetSolutionStepValue(NODAL_MASS);
		    nodal_mass = 0.0;
		    double & particle_mass = (i)->FastGetSolutionStepValue(PARTICLE_MASS);
		    particle_mass = 0.0;

 		    double & nodal_area = (i)->FastGetSolutionStepValue(NODAL_AREA);
		    nodal_area = 0.0;
		  }
	      }
 
    
	    array_1d<double,3> zero = ZeroVector(3);
	    //set WORK = VELOCITY of the old step
	    
	    
	    for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin(); i != BaseType::GetModelPart().NodesEnd(); ++i){
	      noalias(i->FastGetSolutionStepValue(FORCE)) =    zero;
	      noalias(i->FastGetSolutionStepValue(ANGULAR_ACCELERATION)) =    zero;
	      noalias(i->FastGetSolutionStepValue(ACCELERATION)) =    zero;		

	    }
	    
            //add the elemental contributions for the calculation of the velocity
            //and the determination of the nodal area
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
                    //array_1d<double, 3 > & press_proj = (i)->FastGetSolutionStepValue(PRESS_PROJ);
                    //array_1d<double, 3 > & conv_proj = (i)->FastGetSolutionStepValue(CONV_PROJ);
                    //double A = (i)->FastGetSolutionStepValue(NODAL_MASS);
		    //array_1d<double, 3 > & v = (i)->FastGetSolutionStepValue(VELOCITY);

		    array_1d<double,3>& force_temp = i->FastGetSolutionStepValue(FORCE);

		    force_temp *=(1.0/ i->FastGetSolutionStepValue(NODAL_MASS));
		    //force_temp *=(1.0/ i->FastGetSolutionStepValue(PARTICLE_MASS));

		    //array_1d<double,3>& force_temp_g = i->FastGetSolutionStepValue(ANGULAR_ACCELERATION);
 		    //force_temp_g *=(1.0/ i->FastGetSolutionStepValue(NODAL_MASS));
		    //force_temp_g *=(1.0/ i->FastGetSolutionStepValue(PARTICLE_MASS));

		    array_1d<double,3>& acc_temp = i->FastGetSolutionStepValue(ACCELERATION);

			

		    acc_temp =force_temp;
		    if(i->IsFixed(VELOCITY_X) == true){
		      noalias(i->FastGetSolutionStepValue(FORCE)) =    zero;
		      noalias(i->FastGetSolutionStepValue(ANGULAR_ACCELERATION)) =    zero;

		      noalias(i->FastGetSolutionStepValue(ACCELERATION)) =    zero;
		      
		    }
		  }
	      }
	    
	    Timer::Stop("paso_3");
	    KRATOS_WATCH(time)
	      
	      KRATOS_CATCH("");
	    
	    
        }
	
        //******************************************************************************************************
        //******************************************************************************************************
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
        virtual void SetEchoLevel(int Level)
        {
	  mecho_level = Level;
	  mpfracvel_strategy->SetEchoLevel(Level);
	  mppressurestep->SetEchoLevel(Level);
        }
	
        //******************************************************************************************************
        //******************************************************************************************************
	
        virtual void Clear()
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

        GenerateSlipConditionProcess::Pointer mpSlipProcess;
        bool mHasSlipProcess;

        std::vector< Process::Pointer > mInitializeIterationProcesses;

        SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>& msolver_config;

        //******************************************************************************************
        //******************************************************************************************

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






