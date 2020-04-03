//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson (thanks Klaus Sautter)
//
//


#if !defined(KRATOS_MPM_EXPLICIT_STRATEGY )
#define  KRATOS_MPM_EXPLICIT_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "includes/kratos_flags.h"

#include "solving_strategies/strategies/solving_strategy.h"

// Application includes
#include "particle_mechanics_application_variables.h"

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
/**
 * @class MPMExplicitStrategy
 * @ingroup ParticleMechanicsApplciation
 * @brief This strategy is used for the explicit time integration
 * @author Peter Wilson (based on the work of Klaus B Sautter)
 */
    template<class TSparseSpace,
        class TDenseSpace,
        class TLinearSolver
    >
        class MPMExplicitStrategy
        : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    {
    public:
        /**@name Type Definitions */
        /*@{ */
        /** Counted pointer of ClassName */
        KRATOS_CLASS_POINTER_DEFINITION(MPMExplicitStrategy);

        typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

        typedef typename BaseType::TDataType TDataType;

        typedef TSparseSpace SparseSpaceType;

        typedef typename BaseType::TSchemeType TSchemeType;

        typedef typename BaseType::DofsArrayType DofsArrayType;

        typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

        typedef typename BaseType::TSystemVectorType TSystemVectorType;

        typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

        typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

        typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
        typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

        typedef typename BaseType::ElementsArrayType ElementsArrayType;

        typedef typename BaseType::NodesArrayType NodesArrayType;

        typedef typename BaseType::ConditionsArrayType ConditionsArrayType;


        /*@} */
        /**@name Life Cycle
         */
         /*@{ */
        MPMExplicitStrategy(
            ModelPart& model_part,
            typename TSchemeType::Pointer pScheme,
            bool CalculateReactions = false,
            bool ReformDofSetAtEachStep = false,
            bool MoveMeshFlag = false
            )
            : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, MoveMeshFlag),
            mpScheme(pScheme),
            mReformDofSetAtEachStep(ReformDofSetAtEachStep),
            mCalculateReactionsFlag(CalculateReactions)
        {
            KRATOS_TRY

            mKeepSystemConstantDuringIterations = false;

            // Set flags to start correcty the calculations
            mSolutionStepIsInitialized = false;

            mInitializeWasPerformed = false;

            mFinalizeSolutionStep = true;

            // Set EchoLevel to the default value (only time is displayed)
            SetEchoLevel(1);

            // By default the matrices are rebuilt at each solution step
            this->SetRebuildLevel(1);

            KRATOS_CATCH("")
        }

        /** Destructor.
         */
        virtual ~MPMExplicitStrategy()
        {
        }

        /** Destructor.
         */

         //Set and Get Scheme ... containing Builder, Update and other

        void SetScheme(typename TSchemeType::Pointer pScheme)
        {
            mpScheme = pScheme;
        };

        typename TSchemeType::Pointer GetScheme()
        {
            return mpScheme;
        };


        void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
        {
            mInitializeWasPerformed = InitializePerformedFlag;
        }

        bool GetInitializePerformedFlag()
        {
            return mInitializeWasPerformed;
        }

        void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
        {
            mCalculateReactionsFlag = CalculateReactionsFlag;
        }

        bool GetCalculateReactionsFlag()
        {
            return mCalculateReactionsFlag;
        }

        bool GetReformDofSetAtEachStepFlag()
        {
            return mReformDofSetAtEachStep;
        }

        void SetFinalizeSolutionStepFlag(bool FinalizeSolutionStepFlag = true)
        {
            mFinalizeSolutionStep = FinalizeSolutionStepFlag;
        }

        bool GetFinalizeSolutionStepFlag()
        {
            return mFinalizeSolutionStep;
        }

        //level of echo for the solving strategy
        // 0 -> mute... no echo at all
        // 1 -> printing time and basic informations
        // 2 -> printing linear solver data
        // 3 -> Print of debug informations:
        //		Echo of stiffness matrix, Dx, b...

        void SetEchoLevel(int Level) override
        {
            BaseType::mEchoLevel = Level;
        }

        //*********************************************************************************
        /**OPERATIONS ACCESSIBLE FROM THE INPUT:*/
        //*********************************************************************************
        /**
        Initialize members
         */
         //**********************************************************************


        void Initialize() override
        {
            KRATOS_TRY

                typename TSchemeType::Pointer pScheme = GetScheme();

            // OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
            // if the operations needed were already performed this does nothing
            if (mInitializeWasPerformed == false)
            {
                KRATOS_INFO_IF("MPM_Explicit_Strategy", this->GetEchoLevel() > 1) << "Initializing solving strategy" << std::endl;
                KRATOS_ERROR_IF(mInitializeWasPerformed == true) << "Initialization was already performed " << mInitializeWasPerformed << std::endl;

                // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
                KRATOS_INFO_IF("MPM_Explicit_Strategy", this->GetEchoLevel() > 1) << "Initializing scheme" << std::endl;
                if (pScheme->SchemeIsInitialized() == false)
                    pScheme->Initialize(BaseType::GetModelPart());

                // Initialize The Elements - OPERATIONS TO BE DONE ONCE
                KRATOS_INFO_IF("MPM_Explicit_Strategy", this->GetEchoLevel() > 1) << "Initializing elements" << std::endl;
                if (pScheme->ElementsAreInitialized() == false)
                    pScheme->InitializeElements(BaseType::GetModelPart());

                // Initialize The Conditions - OPERATIONS TO BE DONE ONCE
                KRATOS_INFO_IF("MPM_Explicit_Strategy", this->GetEchoLevel() > 1) << "Initializing conditions" << std::endl;
                if (pScheme->ConditionsAreInitialized() == false)
                    pScheme->InitializeConditions(BaseType::GetModelPart());

                mInitializeWasPerformed = true;
            }

            // Prints informations about the current time
            if (this->GetEchoLevel() == 2 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
            {
                KRATOS_INFO("MPM_Explicit_Strategy") << "CurrentTime = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;
            }

            KRATOS_CATCH("")
        }

        void InitializeSolutionStep() override
        {
            KRATOS_TRY

                // Initialize solution step
                if (mSolutionStepIsInitialized == false)
                {
                    typename TSchemeType::Pointer pScheme = GetScheme();

                    ModelPart& r_model_part = BaseType::GetModelPart();

                    Vector dummy_vector;
                    TSystemVectorType rDx = TSystemVectorType();
                    TSystemVectorType rb = TSystemVectorType();


                    TSystemMatrixType mA = TSystemMatrixType();
                    TSystemVectorType mDx = TSystemVectorType();
                    TSystemVectorType mb = TSystemVectorType();

                    // Initial operations ... things that are constant over the Solution Step
                    pScheme->InitializeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
                }

            mSolutionStepIsInitialized = true;

            KRATOS_INFO_IF("MPM_Explicit_Strategy", this->GetEchoLevel() >= 3) << "Initialize Solution Step in strategy finished" << std::endl;

            KRATOS_CATCH("")
        }


        //*********************************************************************************
        /**
        the problem of interest is solved
         */
         //**********************************************************************
        bool SolveSolutionStep() override
        {
            typename TSchemeType::Pointer pScheme = GetScheme();
            DofsArrayType dof_set_dummy;
            TSystemMatrixType mA = TSystemMatrixType();
            TSystemVectorType mDx = TSystemVectorType();
            TSystemVectorType mb = TSystemVectorType();

            pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            // Compute residual forces on the model part
            this->CalculateAndAddRHS(pScheme, BaseType::GetModelPart());

            pScheme->Update(BaseType::GetModelPart(), dof_set_dummy, mA, mDx, mb);

            //pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb); // TODO check if we need this

            // Calculate reactions if required
            if (mCalculateReactionsFlag) {
                CalculateReactions(pScheme, BaseType::GetModelPart(), mA, mDx, mb);
            }

            return true;
        }

        //*********************************************************************************

        //**********************************************************************
        //**********************************************************************

        void Clear() override
        {
            KRATOS_TRY
                GetScheme()->Clear();

            KRATOS_CATCH("")
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

        void SetKeepSystemConstantDuringIterations(bool value)
        {
            mKeepSystemConstantDuringIterations = value;
        }

        bool GetKeepSystemConstantDuringIterations()
        {
            return mKeepSystemConstantDuringIterations;
        }


        /*@} */
        /**@name Inquiry */
        /*@{ */


        /*@} */
        /**@name Friends */
        /*@{ */


        /*@} */

    private:
        /**@name Protected static Member Variables */
        /*@{ */


        /*@} */
        /**@name Protected member Variables */
        /*@{ */


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

    protected:
        /**@name Static Member Variables */
        /*@{ */


        /*@} */
        /**@name Member Variables */
        /*@{ */

        typename TSchemeType::Pointer mpScheme;



        /**
        Flag telling if it is needed to reform the DofSet at each
        solution step or if it is possible to form it just once
        - true  => reforme at each time step
        - false => form just one (more efficient)

        Default = false
         */
        bool mReformDofSetAtEachStep;

        /**
        Flag telling if it is needed or not to compute the reactions

        default = true
         */
        bool mCalculateReactionsFlag;

        bool mSolutionStepIsInitialized;

        bool mInitializeWasPerformed;

        //flag to allow keeping system matrix constant during iterations
        bool mKeepSystemConstantDuringIterations;

        //flag to allow to not finalize the solution step, so the historical variables are not updated
        bool mFinalizeSolutionStep;

        /*@} */
        /**@name Private Operators*/
        /*@{ */
        //**********************************************************************
        //**********************************************************************




        //**********************************************************************
        //**********************************************************************
        void FinalizeSolutionStep() override
        {
            KRATOS_TRY
                typename TSchemeType::Pointer pScheme = GetScheme();
            TSystemMatrixType mA = TSystemMatrixType();
            TSystemVectorType mDx = TSystemVectorType();
            TSystemVectorType mb = TSystemVectorType();


            /*Finalization of the solution step,
            operations to be done after achieving convergence, for example the
            Final Residual Vector (mb) has to be saved in there
            to avoid error accumulation*/
            if (mFinalizeSolutionStep)
            {
                KRATOS_INFO_IF("MPM_Explicit_Strategy", this->GetEchoLevel() >= 3) << "Calling FinalizeSolutionStep" << std::endl;

                pScheme->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
                if (BaseType::MoveMeshFlag()) BaseType::MoveMesh();
            }

            // Cleaning memory after the solution
            pScheme->Clean();

            // Reset flags for next step
            mSolutionStepIsInitialized = false;
            KRATOS_CATCH("")
        }

        /**
         * function to perform expensive checks.
         * It is designed to be called ONCE to verify that the input is correct.
         */
        int Check() override
        {
            KRATOS_TRY

                BaseType::Check();
            GetScheme()->Check(BaseType::GetModelPart());
            return 0;

            KRATOS_CATCH("")
        }


        void CalculateAndAddRHS(
            typename TSchemeType::Pointer pScheme,
            ModelPart& rModelPart
            )
        {
            KRATOS_TRY

            ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
            ConditionsArrayType& r_conditions = rModelPart.Conditions();
            ElementsArrayType& r_elements = rModelPart.Elements();

            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
            Element::EquationIdVectorType equation_id_vector_dummy; // Dummy
            
            #pragma omp parallel for firstprivate(RHS_Contribution, equation_id_vector_dummy), schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_conditions.size()); ++i) {
                auto it_cond = r_conditions.begin() + i;
                pScheme->Condition_Calculate_RHS_Contribution((*it_cond.base()), RHS_Contribution, equation_id_vector_dummy, r_current_process_info);
            }

            #pragma omp parallel for firstprivate(RHS_Contribution, equation_id_vector_dummy), schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
                auto it_elem = r_elements.begin() + i;
                pScheme->Calculate_RHS_Contribution((*it_elem.base()), RHS_Contribution, equation_id_vector_dummy, r_current_process_info);
            }

            KRATOS_CATCH("")
        }


        void CalculateReactions(
            typename TSchemeType::Pointer pScheme,
            ModelPart& rModelPart,
            TSystemMatrixType& rA,
            TSystemVectorType& rDx,
            TSystemVectorType& rb
            )
        {
            // We iterate over the nodes
            auto& r_nodes = rModelPart.Nodes();

            // If we consider rotation dofs
            const bool has_dof_for_rot_z = (r_nodes.begin())->HasDofFor(ROTATION_Z);

            // Auxiliar values
            const array_1d<double, 3> zero_array = ZeroVector(3);
            array_1d<double, 3> force_residual = ZeroVector(3);
            array_1d<double, 3> moment_residual = ZeroVector(3);

            // Getting
            const auto it_node_begin = r_nodes.begin();
            const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);
            const IndexType rotppos = it_node_begin->GetDofPosition(ROTATION_X);

            // Iterating nodes
#pragma omp parallel for firstprivate(force_residual, moment_residual), schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                auto it_node = it_node_begin + i;

                noalias(force_residual) = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);
                if (has_dof_for_rot_z) {
                    noalias(moment_residual) = it_node->FastGetSolutionStepValue(MOMENT_RESIDUAL);
                }
                else {
                    noalias(moment_residual) = zero_array;
                }

                if (it_node->GetDof(DISPLACEMENT_X, disppos).IsFixed()) {
                    double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_X);
                    r_reaction = force_residual[0];
                }
                if (it_node->GetDof(DISPLACEMENT_Y, disppos + 1).IsFixed()) {
                    double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_Y);
                    r_reaction = force_residual[1];
                }
                if (it_node->GetDof(DISPLACEMENT_Z, disppos + 2).IsFixed()) {
                    double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_Z);
                    r_reaction = force_residual[2];
                }
                if (has_dof_for_rot_z) {
                    if (it_node->GetDof(ROTATION_X, rotppos).IsFixed()) {
                        double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_MOMENT_X);
                        r_reaction = moment_residual[0];
                    }
                    if (it_node->GetDof(ROTATION_Y, rotppos + 1).IsFixed()) {
                        double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_MOMENT_Y);
                        r_reaction = moment_residual[1];
                    }
                    if (it_node->GetDof(ROTATION_Z, rotppos + 2).IsFixed()) {
                        double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_MOMENT_Z);
                        r_reaction = moment_residual[2];
                    }
                }
            }
        }



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
        MPMExplicitStrategy(const MPMExplicitStrategy& Other)
        {
        };


        /*@} */

    }; /* Class MPMExplicitStrategy */

    /*@} */

    /**@name Type Definitions */
    /*@{ */


    /*@} */

}; /* namespace Kratos.*/

#endif /* KRATOS_MPM_EXPLICIT_STRATEGY  defined */
