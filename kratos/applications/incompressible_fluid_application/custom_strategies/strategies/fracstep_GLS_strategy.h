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
            ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
            array_1d<double, 3 > zero = ZeroVector(3);
            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];


	    ModelPart& model_part=BaseType::GetModelPart();
		
	    const double dt = model_part.GetProcessInfo()[DELTA_TIME];
				
	   
	    for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	      {
		//it->FastGetSolutionStepValue(AUX_VECTOR)=ZeroVector(3);
		it->FastGetSolutionStepValue(FORCE)=ZeroVector(3);
	      }		
	    
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
		
		//geom[0].FastGetSolutionStepValue(AUX_VECTOR) += aux;
		geom[0].FastGetSolutionStepValue(FORCE) += aux;
		//reusing aux for the second node 
		aux[0]=aaa[2];
		aux[1]=aaa[3];			
		//z-component is zero
		//geom[1].FastGetSolutionStepValue(AUX_VECTOR) += aux;
		geom[1].FastGetSolutionStepValue(FORCE) += aux;
		//reusing aux for the third node
		aux[0]=aaa[4];
		aux[1]=aaa[5];			
		//geom[2].FastGetSolutionStepValue(AUX_VECTOR) += aux;
		geom[2].FastGetSolutionStepValue(FORCE) += aux;
	      }
	    
	    //correct the velocities
	    for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
	      {
		//VELOCITY = VELOCITY + dt * Minv * AUX_VECTOR
		double dt_Minv = (dt / 2.00) / it->FastGetSolutionStepValue(NODAL_MASS);
		//array_1d<double,3>& temp = it->FastGetSolutionStepValue(AUX_VECTOR);
		array_1d<double,3>& force_temp = it->FastGetSolutionStepValue(FORCE);
		force_temp *= dt_Minv;//(1.0/ it->FastGetSolutionStepValue(NODAL_MASS));
		//KRATOS_WATCH(force_temp);
		
		if(!it->IsFixed(VELOCITY_X))
		  {
				it->FastGetSolutionStepValue(VELOCITY_X)+=force_temp[0];
		  }
		if(!it->IsFixed(VELOCITY_Y))
		  {
		    it->FastGetSolutionStepValue(VELOCITY_Y)+=force_temp[1];						
		  }
		if(!it->IsFixed(VELOCITY_Z))
		  {
		    it->FastGetSolutionStepValue(VELOCITY_Z)+=force_temp[2];						
				}
		
		
	      }
	    //Timer::Stop("Ultimo_paso");
	
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
	  KRATOS_WATCH(time)
	    
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
                    /*array_1d<double, 3 > & press_proj = (i)->FastGetSolutionStepValue(PRESS_PROJ);
		      array_1d<double, 3 > & conv_proj = (i)->FastGetSolutionStepValue(CONV_PROJ);
		      noalias(press_proj) = zero;
		      noalias(conv_proj) = zero;*/
		    double & nodal_mass = (i)->FastGetSolutionStepValue(NODAL_MASS);
		    nodal_mass = 0.0;
		  }
	      }
	    
	    
	    array_1d<double,3> zero = ZeroVector(3);
	    //set WORK = VELOCITY of the old step
	    
	    
	    for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin(); i != BaseType::GetModelPart().NodesEnd(); ++i){
	      //noalias(i->FastGetSolutionStepValue(AUX_VECTOR)) = i->FastGetSolutionStepValue(VELOCITY,1);	
	      noalias(i->FastGetSolutionStepValue(FORCE)) =    zero;		
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
	    
	    
            /*BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(NODAL_MASS);
            BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(PRESS_PROJ);
            BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(CONV_PROJ);*/
	    
            //solve nodally for the velocity
	    
#pragma omp parallel for schedule(static,1)
            for (int k = 0; k < number_of_threads; k++)
	      {
                ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
                ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k + 1];
		
                for (ModelPart::NodeIterator i = it_begin; i != it_end; ++i)
		  {
                    //array_1d<double, 3 > & press_proj = (i)->FastGetSolutionStepValue(PRESS_PROJ);
                    //array_1d<double, 3 > & conv_proj = (i)->FastGetSolutionStepValue(CONV_PROJ);
                    double A = (i)->FastGetSolutionStepValue(NODAL_MASS);
		    //array_1d<double, 3 > & v = (i)->FastGetSolutionStepValue(VELOCITY);
		    array_1d<double,3>& force_temp = i->FastGetSolutionStepValue(FORCE);
		    force_temp *=(1.0/ i->FastGetSolutionStepValue(NODAL_MASS));
		    if(i->IsFixed(VELOCITY_X) == true){
		      noalias(i->FastGetSolutionStepValue(FORCE)) =    zero;
		      //noalias(i->FastGetSolutionStepValue(VELOCITY)) =    zero;
		      
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
			array_1d<double,3> aux;	
			array_1d<double,3> aux1;

			ModelPart& model_part=BaseType::GetModelPart();
			const double dt = model_part.GetProcessInfo()[DELTA_TIME];


			array_1d<double,3> zero = ZeroVector(3);

			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin(); i != BaseType::GetModelPart().NodesEnd(); ++i)
			{ 
			noalias(i->FastGetSolutionStepValue(FORCE)) =    zero;		
			}
			
			
			ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
            		rCurrentProcessInfo[FRACTIONAL_STEP] = 5;
			
			//loop over elements calculating the Right Hand Side, that is stored directly to the node.. this is done by fct Calculate
			for(ModelPart::ElementIterator im = model_part.ElementsBegin() ; im != model_part.ElementsEnd() ; ++im)
			{
			//compute the momentum residual, add it to the RHS_VECTOR on nodes
 			(im)->InitializeSolutionStep(BaseType::GetModelPart().GetProcessInfo());
			}
			
			
			for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin(); i != BaseType::GetModelPart().NodesEnd(); ++i)
			{	

			array_1d<double,3>& force_temp = i->FastGetSolutionStepValue(FORCE);
			//KRATOS_WATCH(it->FastGetSolutionStepValue(NODAL_MASS))
			force_temp -=i->FastGetSolutionStepValue(NODAL_MASS) * ((i)->FastGetSolutionStepValue(VELOCITY)-(i)->FastGetSolutionStepValue(VELOCITY,1))/dt;
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






