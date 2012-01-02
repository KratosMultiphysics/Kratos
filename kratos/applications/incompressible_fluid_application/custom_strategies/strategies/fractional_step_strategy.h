/* *********************************************************   
 *
 *   Last Modified by:    $Author: rrossi $
 *   Date:                $Date: 2009-01-13 15:39:56 $
 *   Revision:            $Revision: 1.14 $
 *
 * ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_FRACTIONALSTEP_STRATEGY )
#define  KRATOS_RESIDUALBASED_FRACTIONALSTEP_STRATEGY


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

#ifdef _OPENMP
#include "omp.h"
#endif

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

    /// FractionalStep Solver

    /**   This class implements a compound strategy that allows the solution of the
     * Navier Stokes problem using a multi-step strategy.
     * The solver is based on the Algebraic splitting proposed by Ramon Codina.
     * The solver is designed to work EXCLUSIVELY with @see Fluid2D and @see Fluid3D elements.
     * Although the solver is implicit, and theoretically no time step
     * is needed for stability, it is typically needed to provide a restriction to allow
     * the convergence of the convection non-lienarity
     */
    template<class TSparseSpace,
    class TDenseSpace,
    class TLinearSolver
    >
    class FractionalStepStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    {
    public:
        /**@name Type Definitions */
        /*@{ */

        /** Counted pointer of ClassName */
        typedef boost::shared_ptr< FractionalStepStrategy<TSparseSpace, TDenseSpace, TLinearSolver> > Pointer;

        typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

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

        /**
         * Constructor of the FractionalStepStrategy. Implements the solutions strategy for a Navier Stokes solver
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
        FractionalStepStrategy(
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
            this->Check();

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
        virtual ~FractionalStepStrategy()
        {
        }

        /** Destructor.
         */

        //*********************************************************************************
        //**********************************************************************

        /**
         * This function performs the compound solve of the problem, that is devolves a solution
         * ready to be advanced in time
         * @return returns the norm of the pressure correction vector
         */
        double Solve()
        {
            KRATOS_TRY
	Timer::Start("solve");
            //assign the correct fractional step coefficients (BDF_COEFFICIENTS..)
            InitializeFractionalStep(this->m_step, this->mtime_order);
            double Dp_norm;

            //predicting the velocity
            PredictVelocity(this->m_step, this->mprediction_order);

            //initialize projections at the first steps
            InitializeProjections(this->m_step, this->proj_is_initialized);

            //Assign Velocity To Fract Step Velocity and Node Area to Zero
            AssignInitialStepValues();

            if (this->m_step <= this->mtime_order)
                Dp_norm = IterativeSolve();
            else
            {
                if (this->mpredictor_corrector == false) //standard fractional step
                    Dp_norm = FracStepSolution();
                else //iterative solution
                    Dp_norm = IterativeSolve();
            }

            if (this->mReformDofAtEachIteration == true)
                this->Clear();

            this->m_step += 1;
	Timer::Stop("solve");
            return Dp_norm;
            KRATOS_CATCH("")
        }


        //*********************************************************************************
        //**********************************************************************

        /**
         * internal function that implements predictor corrector strategy
         * @return  the norm of the pressure correction vector
         */
        double IterativeSolve()
        {
            KRATOS_TRY


                    double Dp_norm = 1.00;
            int iteration = 0;

            int MaxPressureIterations = this->mMaxPressIterations;

            int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
            while (Dp_norm >= this->mpressure_toll && iteration++ < MaxPressureIterations)
            {
                double p_norm = SavePressureIteration();

                Dp_norm = FracStepSolution();

                if (fabs(p_norm) > 1e-10)
                    Dp_norm /= p_norm;
                else
                    Dp_norm = 1.0;

                if (rank == 0) std::cout << "it = " << iteration << " Pressure Variation Norm = " << Dp_norm << std::endl;

            }

            if (this->mReformDofAtEachIteration == true)
                this->Clear();

            return Dp_norm;
            KRATOS_CATCH("")
        }

        //*********************************************************************************
        //**********************************************************************

        /**
         * utility to activate an auxiliary process for the computation of the slip condition
         * @param pSlipProcess
         */
        void SetSlipProcess(GenerateSlipConditionProcess::Pointer pSlipProcess)
        {
            mpSlipProcess = pSlipProcess;
            mHasSlipProcess = true;
        }

        //*********************************************************************************
        //**********************************************************************

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
                (i)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = p;

                local_p_norm += p*p;
            }

            double p_norm = BaseType::GetModelPart().GetCommunicator().SumAll(local_p_norm);

            //TODO: prepare for parallelization
            p_norm = sqrt(p_norm);

            return p_norm;
            KRATOS_CATCH("")
        }

        //******************************************************************************************************
        //******************************************************************************************************
        //explicit correction for velocities and eventually stabilization terms

        /**
         * routine that implements the complete solution step
         *          * @return  returns the norm of the pressure correction vector
         */
        double FracStepSolution()
        {
            KRATOS_TRY

                    int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

            //setting the fractional velocity to the value of the velocity
            AssignInitialStepValues();

            //solve first step for fractional step velocities
            boost::timer step1time;
            this->SolveStep1(this->mvelocity_toll, this->mMaxVelIterations);
            if (rank == 0) std::cout << "step1 time " << step1time.elapsed() << std::endl;

            //solve for pressures (and recalculate the nodal area)
            boost::timer step2time;
            double Dp_norm = this->SolveStep2();
            if (rank == 0) std::cout << "pressure calculation time " << step2time.elapsed() << std::endl;


            this->ActOnLonelyNodes();

            //calculate projection terms
            boost::timer projection_time;
            this->SolveStep3();
            if (rank == 0) std::cout << "projection calculation time " << projection_time.elapsed() << std::endl;

            //			if(mdomain_size == 2)
            //  				this->SolveStep2_Mp();


            //correct velocities
            boost::timer vel_time;
            this->SolveStep4();
            if (rank == 0) std::cout << "velocity correction time " << vel_time.elapsed() << std::endl;

            return Dp_norm;

            KRATOS_CATCH("")
        }

        //******************************************************************************************************
        //******************************************************************************************************

        /**
         * Implements the third step of the FractionalStep scheme, in which the velocity is corrected
         * according to the pressure correction derived from the second step
         */
        void SolveStep4()
        {
            KRATOS_TRY;

            ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
            array_1d<double, 3 > zero = ZeroVector(3);
            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];


            //first of all set to zero the nodal variables to be updated nodally
            for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                    i != BaseType::GetModelPart().NodesEnd(); ++i)
            {
                array_1d<double, 3 > & fract_v = (i)->FastGetSolutionStepValue(FRACT_VEL);
                fract_v *= (i)->FastGetSolutionStepValue(NODAL_MASS) * BDFcoeffs[0];
            }

            //set to zero fract_v on ghost nodes --> does nothing on serial version
            for (ModelPart::NodeIterator i = BaseType::GetModelPart().GetCommunicator().GhostMesh().NodesBegin();
                    i != BaseType::GetModelPart().GetCommunicator().GhostMesh().NodesEnd(); ++i)
            {
                array_1d<double, 3 > & fract_v = (i)->FastGetSolutionStepValue(FRACT_VEL);
                noalias(fract_v) = ZeroVector(3);
            }

            //add the elemental contributions for the calculation of the velocity
            //and the determination of the nodal area
            rCurrentProcessInfo[FRACTIONAL_STEP] = 6;
            for (ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin();
                    i != BaseType::GetModelPart().ElementsEnd(); ++i)
            {
                (i)->InitializeSolutionStep(rCurrentProcessInfo);
            }

            BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(FRACT_VEL);


            //solve nodally for the velocity
            if (this->mdomain_size == 2)
            {
                for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                        i != BaseType::GetModelPart().NodesEnd(); ++i)
                {
                    array_1d<double, 3 > & v = (i)->FastGetSolutionStepValue(VELOCITY);
                    const array_1d<double, 3 > & fract_v = (i)->FastGetSolutionStepValue(FRACT_VEL);
                    double A = (i)->FastGetSolutionStepValue(NODAL_MASS);

                    double temp = (1.0 / BDFcoeffs[0]) / A;
                    if (!i->IsFixed(VELOCITY_X))
                    {
                        v[0] = fract_v[0] * temp;
                    }
                    if (!i->IsFixed(VELOCITY_Y))
                    {
                        v[1] = fract_v[1] * temp;
                    }
                }
            } else if (this->mdomain_size == 3)
            {
                for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                        i != BaseType::GetModelPart().NodesEnd(); ++i)
                {
                    array_1d<double, 3 > & v = (i)->FastGetSolutionStepValue(VELOCITY);
                    const array_1d<double, 3 > & fract_v = (i)->FastGetSolutionStepValue(FRACT_VEL);
                    double A = (i)->FastGetSolutionStepValue(NODAL_MASS);

                    double temp = (1.0 / BDFcoeffs[0]) / A;
                    if (!i->IsFixed(VELOCITY_X))
                    {
                        v[0] = fract_v[0] * temp;
                    }
                    if (!i->IsFixed(VELOCITY_Y))
                    {
                        v[1] = fract_v[1] * temp;
                    }
                    if (!i->IsFixed(VELOCITY_Z))
                    {
                        v[2] = fract_v[2] * temp;
                    }
                }
            }

            //if we have slip condition apply it
            if (mHasSlipProcess == true)
            {
                mpSlipProcess->SetNormalVelocityToZero(VELOCITY);
            }



            KRATOS_CATCH("");
        }

        //******************************************************************************************************
        //******************************************************************************************************
        //sets the BDF coefficients to the correct value

        /**
         * computes the time coefficients
         * @param step number (the first time the system is solved this is 1)
         * @param time_order time accuracy of the solver. 1-BDF1-first order, 2-BDF2-second order
         */
        void InitializeFractionalStep(const int step, const int time_order)
        {
            KRATOS_TRY;

            //calculate the BDF coefficients
            ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
            double Dt = rCurrentProcessInfo[DELTA_TIME];

            if (time_order == 2 && step > time_order)
            {
                if (BaseType::GetModelPart().GetBufferSize() < 3)
                    KRATOS_ERROR(std::logic_error, "insufficient buffer size for BDF2", "")

                    double dt_old = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

                double rho = dt_old / Dt;
                double coeff = 1.0 / (Dt * rho * rho + Dt * rho);

                rCurrentProcessInfo[BDF_COEFFICIENTS].resize(3, false);
                Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
                BDFcoeffs[0] = coeff * (rho * rho + 2.0 * rho); //coefficient for step n+1
                BDFcoeffs[1] = -coeff * (rho * rho + 2.0 * rho + 1.0); //coefficient for step n
                BDFcoeffs[2] = coeff;
            } else
            {
                rCurrentProcessInfo[BDF_COEFFICIENTS].resize(2, false);
                Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
                BDFcoeffs[0] = 1.0 / Dt; //coefficient for step n+1
                BDFcoeffs[1] = -1.0 / Dt; //coefficient for step n
            }


            KRATOS_CATCH("");
        }

        //******************************************************************************************************
        //******************************************************************************************************

        /**
         * helper function that initialized the pressure projection to the body force
         * this is important at the first step, since the pressure does not yet have a sensible
         * distribution
         * @param step - solution step (1 at the first solution step, then incremented)
         * @param proj_is_initialized flag that controls if this function was called or not
         */
        void InitializeProjections(int step, bool proj_is_initialized)
        {
            if (step <= 2 && proj_is_initialized == false)
            {
                for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                        i != BaseType::GetModelPart().NodesEnd(); ++i)
                {
                    noalias(i->FastGetSolutionStepValue(PRESS_PROJ)) = i->FastGetSolutionStepValue(BODY_FORCE);
                }
                proj_is_initialized = true;
            }
        }

        //******************************************************************************************************
        //******************************************************************************************************
        //predict values for the fractional step velocities
        //and set to zero the nodal mass

        void AssignInitialStepValues()
        {
            KRATOS_TRY
            for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                    i != BaseType::GetModelPart().NodesEnd(); ++i)
            {
                //predicting the values for the fluid velocity
                array_1d<double, 3 > & v = (i)->FastGetSolutionStepValue(VELOCITY);
                array_1d<double, 3 > & fracv = (i)->FastGetSolutionStepValue(FRACT_VEL);
                noalias(fracv) = v;

                //setting the old pressure iteration to the value of the pressure
                (i)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = (i)->FastGetSolutionStepValue(PRESSURE);

                //resetting the nodal area
                double area = 0.00;
                (i)->FastGetSolutionStepValue(NODAL_MASS) = area;
            }
            KRATOS_CATCH("");
        }

        //******************************************************************************************************
        //******************************************************************************************************

        void PredictVelocity(int step, int prediction_order)
        {
            KRATOS_TRY
            if (prediction_order == 2)
            {
                if (BaseType::GetModelPart().GetBufferSize() < 3)
                    KRATOS_ERROR(std::logic_error, "insufficient buffer size for second order prediction", "")
                }

            if (prediction_order == 2 && step > 2)
            {
                //second order prediction for the velocity
                for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                        i != BaseType::GetModelPart().NodesEnd(); ++i)
                {
                    array_1d<double, 3 > & vel = i->FastGetSolutionStepValue(VELOCITY);
                    const array_1d<double, 3 > & v1 = i->FastGetSolutionStepValue(VELOCITY, 1);
                    const array_1d<double, 3 > & v2 = i->FastGetSolutionStepValue(VELOCITY, 2);
                    if (!i->IsFixed(VELOCITY_X))
                        vel[0] = 2.00 * v1[0] - v2[0];
                    if (!i->IsFixed(VELOCITY_Y))
                        vel[1] = 2.00 * v1[1] - v2[1];
                    if (!i->IsFixed(VELOCITY_Z))
                        vel[2] = 2.00 * v1[2] - v2[2];
                }
            }
            KRATOS_CATCH("");
        }
        //******************************************************************************************************
        //******************************************************************************************************

        /**
         * this function performs the iterative solution of the non-linear velocity problem in the first step
         * of the fractional step procedure
         * @param velocity_toll - tolerance used in the velocity convergence check
         * @param MaxIterations - max number of iterations
         */
        void SolveStep1(double velocity_toll, int MaxIterations)
        {
            KRATOS_TRY;
            int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

            double normDx = 0.0;

            bool is_converged = false;
            double iteration = 1;

            if (mHasSlipProcess == true)
            {
                mpSlipProcess->SetNormalVelocityToZero(FRACT_VEL);
            }

            //solve for fractional step velocities
            while (is_converged == false && iteration++<MaxIterations)
            {
                //execute initialize iteration processes;
                for(unsigned int i=0; i<mInitializeIterationProcesses.size(); i++)
                    mInitializeIterationProcesses[i]->Execute();

                //perform one iteration over the fractional step velocity
                normDx = FractionalVelocityIteration();
                is_converged = ConvergenceCheck(normDx, velocity_toll);


            }

            if (mHasSlipProcess == true)
            {
                mpSlipProcess->SetNormalVelocityToZero(FRACT_VEL);
                // 			    mpSlipProcess->ApplyEdgeConstraints(FRACT_VEL);
            }

            if (is_converged == false)
                if (rank == 0) std::cout << "ATTENTION: convergence NOT achieved" << std::endl;

            //clear if needed
            if (mReformDofAtEachIteration == true && mpredictor_corrector == false)
            {
                this->mpfracvel_strategy->Clear();
            }

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
            BaseType::GetModelPart().GetProcessInfo()[FRACTIONAL_STEP] = 4;
            return mppressurestep->Solve();
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

                //                            ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
                array_1d<double, 3 > zero = ZeroVector(3);

                //first of all set to zero the nodal variables to be updated nodally
                for (ModelPart::NodeIterator i = it_begin; i != it_end; ++i)
                {
                    array_1d<double, 3 > & press_proj = (i)->FastGetSolutionStepValue(PRESS_PROJ);
                    array_1d<double, 3 > & conv_proj = (i)->FastGetSolutionStepValue(CONV_PROJ);
                    noalias(press_proj) = zero;
                    noalias(conv_proj) = zero;
                }
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

            BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(NODAL_MASS);
            BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(PRESS_PROJ);
            BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(CONV_PROJ);

            //solve nodally for the velocity
#pragma omp parallel for schedule(static,1)
            for (int k = 0; k < number_of_threads; k++)
            {
                ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
                ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k + 1];

                for (ModelPart::NodeIterator i = it_begin; i != it_end; ++i)
                {
                    array_1d<double, 3 > & press_proj = (i)->FastGetSolutionStepValue(PRESS_PROJ);
                    array_1d<double, 3 > & conv_proj = (i)->FastGetSolutionStepValue(CONV_PROJ);
                    double A = (i)->FastGetSolutionStepValue(NODAL_MASS);

                    double temp = 1.00 / A;
                    press_proj *= temp;
                    conv_proj *= temp;

                }
            }


            KRATOS_CATCH("");
        }

        //******************************************************************************************************
        //******************************************************************************************************
        //correct pressure taking in account viscosity (should accelerate the convergence of the predictor corrector)

        void SolveStep2_Mp()
        {
            KRATOS_TRY;

            ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
            array_1d<double, 3 > zero = ZeroVector(3);
            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

            for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                    i != BaseType::GetModelPart().NodesEnd(); ++i)
            {
                //setting the old value of the pressure to the current one
                double& p = (i)->FastGetSolutionStepValue(PRESSURE);
                (i)->FastGetSolutionStepValue(PRESSURE_OLD_IT) = p;
                p = 0.0;
                (i)->FastGetSolutionStepValue(NODAL_MASS) = 0.0;


            }

            //const Vector& BDFcoeffs = BaseType::GetModelPart().GetProcessInfo()[BDF_COEFFICIENTS];

            //calculate divergence of the fractional velocity element by element
            if (this->mdomain_size == 2)
            {
                array_1d<double, 3 > N;
                array_1d<double, 2 > temp, proj_aux;
                array_1d<double, 3 > pressures;
                array_1d<double, 3 > elemental_stabilization;
                boost::numeric::ublas::bounded_matrix <double, 3, 2 > DN_DX;
                array_1d<double, 3 > vg;
                for (ModelPart::ElementIterator i = BaseType::GetModelPart().ElementsBegin();
                        i != BaseType::GetModelPart().ElementsEnd(); ++i)
                {
                    Geometry< Node < 3 > >& geom = i->GetGeometry();
                    double volume;

                    //calculate derivatives
                    GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);

                    //getting fractional velocities on nodes
                    const array_1d<double, 3 > & fv0 = geom[0].FastGetSolutionStepValue(FRACT_VEL);
                    const array_1d<double, 3 > & fv1 = geom[1].FastGetSolutionStepValue(FRACT_VEL);
                    const array_1d<double, 3 > & fv2 = geom[2].FastGetSolutionStepValue(FRACT_VEL);

                    pressures[0] = geom[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
                    pressures[1] = geom[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
                    pressures[2] = geom[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);

                    //calculate N
                    double Gaux;
                    Gaux = DN_DX(0, 0) * fv0[0] + DN_DX(0, 1) * fv0[1];
                    Gaux += DN_DX(1, 0) * fv1[0] + DN_DX(1, 1) * fv1[1];
                    Gaux += DN_DX(2, 0) * fv2[0] + DN_DX(2, 1) * fv2[1];

                    double density = geom[0].FastGetSolutionStepValue(DENSITY);
                    double nu = geom[0].FastGetSolutionStepValue(VISCOSITY);
                    for (int ii = 1; ii < 3; ii++)
                    {
                        density += geom[ii].FastGetSolutionStepValue(DENSITY);
                        nu += geom[ii].FastGetSolutionStepValue(VISCOSITY);
                    }
                    density *= 0.333333333333333333333;
                    nu *= 0.333333333333333333333;

                    noalias(vg) = ZeroVector(3);
                    ; //velocity on the gauss points
                    for (int kk = 0; kk < 3; kk++)
                    {
                        //adding the elemental contribution to the nodal volume
                        geom[kk].FastGetSolutionStepValue(NODAL_MASS) += density * volume * N[kk];

                        //calculating the velocity on the gauss point
                        noalias(vg) += N[kk]* (geom[kk].FastGetSolutionStepValue(FRACT_VEL) -
                                geom[kk].FastGetSolutionStepValue(MESH_VELOCITY));
                    }

                    //calculating the stabilization
                    double dt_contrib_to_tau = 0.0;
                    if (muse_dt_in_stabilization == true)
                        dt_contrib_to_tau = 1.0 / BDFcoeffs[0];
                    double h = sqrt(2.0 * volume);
                    double c1 = 4.00;
                    double c2 = 2.00;
                    double norm_u = norm_2(vg);
                    double tau = 1.00 / (dt_contrib_to_tau + c1 * nu / (h * h) + c2 * norm_u / h);

                    //calculating stabilization laplacian LHS
                    noalias(temp) = prod(trans(DN_DX), pressures);
                    noalias(elemental_stabilization) = -tau / density * prod(DN_DX, temp);

                    const array_1d<double, 3 > & proj_temp = geom[0].FastGetSolutionStepValue(PRESS_PROJ);
                    for (int iii = 0; iii < 2; iii++)
                        proj_aux[iii] = N[0] * proj_temp[iii];
                    for (int kk = 1; kk < 3; kk++)
                    {
                        const array_1d<double, 3 > & proj_temp = geom[kk].FastGetSolutionStepValue(PRESS_PROJ);
                        for (int iii = 0; iii < 2; iii++)
                            proj_aux[iii] += N[kk] * proj_temp[iii];
                    }
                    proj_aux *= tau;
                    noalias(elemental_stabilization) += prod(DN_DX, proj_aux);

                    geom[0].FastGetSolutionStepValue(PRESSURE) += (elemental_stabilization[0] - Gaux * N[0])*(volume * density);
                    geom[1].FastGetSolutionStepValue(PRESSURE) += (elemental_stabilization[1] - Gaux * N[1])*(volume * density);
                    geom[2].FastGetSolutionStepValue(PRESSURE) += (elemental_stabilization[2] - Gaux * N[2])*(volume * density);

                    double nodal_vol = volume / 3.0;
                    geom[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_vol * density;
                    geom[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_vol * density;
                    geom[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_vol * density;
                }
            } else
            {
                KRATOS_ERROR(std::logic_error, "not yet implemented the 3D", "");

            }

            BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(NODAL_MASS);
            BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(PRESSURE);

            //correct pressure
            //double temp;
            for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                    i != BaseType::GetModelPart().NodesEnd(); ++i)
            {
                double& p = (i)->FastGetSolutionStepValue(PRESSURE);
                const double p_old_it = (i)->FastGetSolutionStepValue(PRESSURE_OLD_IT);
                const double rho = (i)->FastGetSolutionStepValue(DENSITY);
                const double nu = (i)->FastGetSolutionStepValue(VISCOSITY);
                const double mass = (i)->FastGetSolutionStepValue(NODAL_MASS);

                if (mass > 1e-12)
                {
                    p *= nu * rho / mass;
                    p += p_old_it;
                } else
                {
                    (i)->FastGetSolutionStepValue(NODAL_MASS) = 1.0;
                    p = p_old_it;
                }

            }

            KRATOS_CATCH("");
        }

        //******************************************************************************************************
        //******************************************************************************************************
        //calculation of projection

        void ActOnLonelyNodes()
        {
            KRATOS_TRY;

            for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                    i != BaseType::GetModelPart().NodesEnd(); ++i)
            {
                double& A = (i)->FastGetSolutionStepValue(NODAL_MASS);

                //the area is zero on lonely nodes, in this case set it to 1.00
                if (A <= 1e-12)
                {
                    A = 1.0;
                }
            }


            KRATOS_CATCH("");
        }

        //******************************************************************************************************
        //******************************************************************************************************

        /**
         * Utility to apply correctly the fractional velocity fixity, following variation of the fixity of velocity
         */
        void ApplyFractionalVelocityFixity()
        {
#ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
#else
            int number_of_threads = 1;
#endif

            vector<unsigned int> partition;
            CreatePartition(number_of_threads, BaseType::GetModelPart().Nodes().size(), partition);

            //                        #pragma omp parallel for schedule(static,1)
            for (int k = 0; k < number_of_threads; k++)
            {
                ModelPart::NodeIterator it_begin = BaseType::GetModelPart().NodesBegin() + partition[k];
                ModelPart::NodeIterator it_end = BaseType::GetModelPart().NodesBegin() + partition[k + 1];
                for (ModelPart::NodeIterator i = it_begin; i != it_end; ++i)
                {
                    if (i->IsFixed(VELOCITY_X))
                    {
                        (i)->FastGetSolutionStepValue(FRACT_VEL_X) = (i)->FastGetSolutionStepValue(VELOCITY_X);
                        (i)->Fix(FRACT_VEL_X);
                    }
                    else
                    {
                        (i)->Free(FRACT_VEL_X);
                    }

                    if (i->IsFixed(VELOCITY_Y))
                    {
                        (i)->FastGetSolutionStepValue(FRACT_VEL_Y) = (i)->FastGetSolutionStepValue(VELOCITY_Y);
                        (i)->Fix(FRACT_VEL_Y);
                    }
                    else
                        (i)->Free(FRACT_VEL_Y);

                    if (i->IsFixed(VELOCITY_Z))
                    {
                        (i)->FastGetSolutionStepValue(FRACT_VEL_Z) = (i)->FastGetSolutionStepValue(VELOCITY_Z);
                        (i)->Fix(FRACT_VEL_Z);
                    }
                    else
                        (i)->Free(FRACT_VEL_Z);
                }
            }
        }


        //******************************************************************************************************
        //******************************************************************************************************

        /**
         * implements the convergence check for the velocities
         * convergence is considered achieved when normDx/norm(v) is less than tol
         * @param normDx norm of the velocity correction
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
                const array_1d<double, 3 > & v = (i)->FastGetSolutionStepValue(FRACT_VEL);

                norm_v += v[0] * v[0];
                norm_v += v[1] * v[1];
                norm_v += v[2] * v[2];
            }

            BaseType::GetModelPart().GetCommunicator().SumAll(norm_v);

            norm_v = sqrt(norm_v);

            if (norm_v == 0.0) norm_v = 1.00;

            double ratio = normDx / norm_v;

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

        //******************************************************************************************************
        //******************************************************************************************************

        /**
         * solves one single iteration of the step1 of the FractionalStep algorithm
         * @param normDx
         */
        double FractionalVelocityIteration()
        {
            KRATOS_TRY

            ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

            rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
            //if we have slip condition apply it
            double normDx = mpfracvel_strategy->Solve();
	    return normDx;

            KRATOS_CATCH("");
        }

        //******************************************************************************************************
        //******************************************************************************************************

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
            if (rank == 0) KRATOS_WATCH("FractionalStepStrategy Clear Function called");
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

        /**
         * This function checks the input extensively to verify that no common error is present
         */
        virtual int Check()
        {
            KRATOS_TRY

                    //veryfying that the model part has all the variables needed
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(FRACT_VEL) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----FRACT_VEL---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----VELOCITY---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(MESH_VELOCITY) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----MESH_VELOCITY---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----PRESSURE---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(PRESSURE_OLD_IT) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----PRESSURE_OLD_IT---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(PRESS_PROJ) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----PRESS_PROJ---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(CONV_PROJ) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----CONV_PROJ---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(NODAL_MASS) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----NODAL_MASS---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(BODY_FORCE) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----BODY_FORCE---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(DENSITY) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----DENSITY---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(VISCOSITY) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----VISCOSITY---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(IS_STRUCTURE) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----IS_STRUCTURE---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(EXTERNAL_PRESSURE) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----EXTERNAL_PRESSURE---- variable!!!!!! ERROR", "");
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(IS_INTERFACE) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----IS_INTERFACE---- variable!!!!!! ERROR", "");

            //check that the domain size is correctly prescribed
            if (this->mdomain_size != msolver_config.GetDomainSize())
                KRATOS_ERROR(std::logic_error, "domain size not coinciding", "")

                //verify buffer size
                if (BaseType::GetModelPart().GetBufferSize() < mtime_order + 1)
                    KRATOS_ERROR(std::logic_error, "insufficient buffer size. Buffer size should be >= time_order+1", "");

            //check that, in the 2D case, the xy plane is used.
            if (this->mdomain_size == 2)
            {
                double zmin = BaseType::GetModelPart().NodesBegin()->Z();
                double zmax = zmin;
                for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
                        i != BaseType::GetModelPart().NodesEnd(); ++i)
                {
                    if (i->Z() < zmin) zmin = i->Z();
                    else if (i->Z() > zmax) zmax = i->Z();
                }
                if (fabs(zmax - zmin) > 1e-20)
                    KRATOS_ERROR(std::logic_error, "2D model is not in the XY plane!", "")
                }

            //verify element type, check that the Id is non zero, and calls the Check function for all of the elements
            if (this->mdomain_size == 2)
            {
                const char ElementName[] = "Fluid2D";
                Element const& ref_el = KratosComponents<Element>::Get(ElementName);

                const char ElementName2[] = "Fluid2Dlevelset";
                Element const& ref_el2 = KratosComponents<Element>::Get(ElementName2);

                for (ModelPart::ElementsContainerType::iterator it = BaseType::GetModelPart().ElementsBegin();
                        it != BaseType::GetModelPart().ElementsEnd(); ++it)
                {
                    if (it->Id() < 1)
                        KRATOS_ERROR(std::logic_error, "Element Id can not be lesser than 1 (0 is not allowed as Id)", "");
                    if (typeid (ref_el) != typeid (*it) && typeid (ref_el2) != typeid (*it))
                    {
                        std::cout << "wrong element found --> " << it->Id() << std::endl;
                        KRATOS_ERROR(std::logic_error, "Fractional step strategy requires Fluid2D element for the 2D case", "");
                    }
                    it->Check(BaseType::GetModelPart().GetProcessInfo());
                }
            } else
            {
                const char ElementName[] = "Fluid3D";
                Element const& ref_el = KratosComponents<Element>::Get(ElementName);

                for (ModelPart::ElementsContainerType::iterator it = BaseType::GetModelPart().ElementsBegin();
                        it != BaseType::GetModelPart().ElementsEnd(); ++it)
                {
                    if (it->Id() < 1)
                        KRATOS_ERROR(std::logic_error, "Element Id can not be lesser than 1 (0 is not allowed as Id)", "");
                    if (typeid (ref_el) != typeid (*it))
                    {
                        std::cout << "wrong element found --> " << it->Id() << std::endl;
                        KRATOS_ERROR(std::logic_error, "Fractional step strategy requires Fluid3D element for the 3D case", "");
                    }
                    it->Check(BaseType::GetModelPart().GetProcessInfo());
                }

            }

            return 0;

            //verify
            KRATOS_CATCH("")
        }

                //******************************************************************************************
        //******************************************************************************************
        /**
         * This function computes the reactions and stores them in "ReactionVar"
         * @param rReactionVar variable used in storing the reactions
         * @param change_sign
         */
        void ComputeReactions(Variable<array_1d<double,3> >& rReactionVar)
        {
            KRATOS_TRY

            //check if the variable used is existing in the model part
            if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(rReactionVar) == false)
                KRATOS_ERROR(std::logic_error, "ReactionVar does not exist! please Add  ----rReactionVar---- variable!!!!!! ERROR", "");

            InitializeFractionalStep(this->m_step, this->mtime_order);

            ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();

            //set reactions to zero
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(BaseType::GetModelPart().Nodes().size()); i++)
            {
                ModelPart::NodesContainerType::iterator it = BaseType::GetModelPart().NodesBegin() + i;
                it->FastGetSolutionStepValue(FRACT_VEL) = it->FastGetSolutionStepValue(VELOCITY);
                it->FastGetSolutionStepValue(PRESSURE_OLD_IT) = it->FastGetSolutionStepValue(PRESSURE);
                it->FastGetSolutionStepValue(rReactionVar) = ZeroVector(3);
            }

            for (unsigned int component = 0; component != this->mdomain_size; component++)
            {
                rCurrentProcessInfo[FRACTIONAL_STEP] = component+1;

                Vector rhs(this->mdomain_size+1);
                Matrix lhs(this->mdomain_size+1,this->mdomain_size+1);

                #pragma omp parallel for firstprivate(rhs,lhs)
                for(int i=0; i<static_cast<int>(BaseType::GetModelPart().Elements().size()); i++)
                {
                    ModelPart::ElementsContainerType::iterator it = BaseType::GetModelPart().ElementsBegin()+i;

                    it->CalculateLocalSystem(lhs,rhs,rCurrentProcessInfo);

                    //now sum contributions where needed
                    Geometry<Node<3> >& geom = it->GetGeometry();
                    for(unsigned int k=0; k<rhs.size(); k++)
                    {
                        array_1d<double,3>& react = geom[k].FastGetSolutionStepValue(rReactionVar);

                        #pragma omp atomic
                        react[component] += rhs[k];
                    }
                }
            }

            BaseType::GetModelPart().GetCommunicator().AssembleCurrentData(rReactionVar);

            KRATOS_CATCH("");
        }


         //******************************************************************************************
        void AddInitializeIterationProcess(Process::Pointer pnew_process)
        {
            KRATOS_TRY
            mInitializeIterationProcesses.push_back(pnew_process);
            KRATOS_CATCH("");
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
        FractionalStepStrategy(const FractionalStepStrategy& Other);


        /*@} */

    }; /* Class FractionalStepStrategy */

    /*@} */

    /**@name Type Definitions */
    /*@{ */


    /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_FRACTIONALSTEP_STRATEGY  defined */

