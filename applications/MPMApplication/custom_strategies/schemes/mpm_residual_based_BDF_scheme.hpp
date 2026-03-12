//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


#if !defined(KRATOS_MPM_RESIDUAL_BASED_BDF_SCHEME )
#define      KRATOS_MPM_RESIDUAL_BASED_BDF_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "containers/array_1d.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residual_based_implicit_time_scheme.h"
//#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "custom_utilities/mpm_boundary_rotation_utility.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_base_condition.h"
#include "utilities/parallel_utilities.h"
#include "utilities/coordinate_transformation_utilities.h"


namespace Kratos
{

/**
 * @class MPMResidualBasedBDFScheme
 * @ingroup KratosMPM
 * @brief BDF1 integration scheme (for linear and nonlinear dynamic problems) for velocities adjusted for Material Point Method
 * @details BDF1 integration scheme (for linear and nonlinear dynamic problems) for velocities adjusted for Material Point Method
 */
template<class TSparseSpace,  class TDenseSpace >
class MPMResidualBasedBDFScheme
    : public ResidualBasedImplicitTimeScheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( MPMResidualBasedBDFScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                                      BaseType;

    typedef ResidualBasedImplicitTimeScheme<TSparseSpace,TDenseSpace>     ImplicitBaseType;

    //typedef ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace> BossakBaseType;

    typedef typename BaseType::TDataType                                   TDataType;

    typedef typename BaseType::DofsArrayType                           DofsArrayType;

    typedef typename Element::DofsVectorType                                DofsVectorType;

    typedef typename BaseType::TSystemMatrixType                   TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                   TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType           LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType           LocalSystemMatrixType;

    typedef ModelPart::ElementsContainerType                             ElementsArrayType;

    typedef ModelPart::ConditionsContainerType                         ConditionsArrayType;

    typedef typename BaseType::Pointer                                     BaseTypePointer;

    typedef CoordinateTransformationUtils<LocalSystemMatrixType, LocalSystemVectorType, double> RotationToolType;
    
    //typedef typename RotationToolType::UniquePointer RotationToolPointerType;

    //using BossakBaseType::mBossak; // !! non necessario

    using ImplicitBaseType::mMatrix;

    /**
     * Constructor.
     * The MPM bdf scheme
     */
    /**
     * @brief Constructor.
     * @detail The MPM bdf method
     * @param rGridModelPart The background grid model part
     * @param DomainSize Domain size or space dimension of the problems
     * @param BlockSize Block size of DOF
     */
    MPMResidualBasedBDFScheme(ModelPart& rGridModelPart, unsigned int DomainSize,
        unsigned int BlockSize, bool IsDynamic = true)
                :ResidualBasedImplicitTimeScheme<TSparseSpace,TDenseSpace>(),
                mGridModelPart(rGridModelPart), mRotationTool(DomainSize, BlockSize)
    {
        // To distinguish quasi-static and dynamic
        mIsDynamic = IsDynamic;

        // For Rotation Utility
        mDomainSize = DomainSize;
        mBlockSize  = BlockSize;

        // For friction-related info
        mFrictionIsActive = mGridModelPart.GetProcessInfo()[FRICTION_ACTIVE];
    }

    /**
     * @brief Copy Constructor.
     */
     MPMResidualBasedBDFScheme(MPMResidualBasedBDFScheme& rOther)
        :ImplicitBaseType(rOther)
        ,mGridModelPart(rOther.mGridModelPart)
        ,mFrictionIsActive(rOther.mFrictionIsActive)
        ,mRotationTool(rOther.mDomainSize,rOther.mBlockSize)
    {
    }

    /**
     * @brief Clone method
     */
    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new MPMResidualBasedBDFScheme(*this) );
    }

    /** Destructor.
     */
    virtual ~MPMResidualBasedBDFScheme
    () {}


    //***************************************************************************
    //***************************************************************************

    void Initialize(ModelPart& rModelPart) override {
        MPMParticleBaseCondition::SetRotationUtility(&mRotationTool);

        BaseType::Initialize(rModelPart);
    }


    /**
     * @brief Performing the update of the solution
     * @details Incremental update within newton iteration. It updates the state variables at the end of the time step u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param rModelPart The model of the problem to solve
     * @param rDofSet Set of all primary variables
     * @param rA LHS matrix
     * @param rDx incremental update of primary variables
     * @param rb RHS Vector
     */
    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb ) override
    {
        KRATOS_TRY

        // Rotate the current displacement to the modified coordinate system since rDx is currently at the modified coordinate system
        mRotationTool.RotateDisplacements(rModelPart); // !! il funzionamento è analogo indipendentemente dalla variabile; vedi mpm_boundary_rotation_utility, riga 268

        // Update of displacement (by DOF)
        mpDofUpdater->UpdateDofs(rDofSet, rDx);

        // Rotate the displacement back to the original coordinate system to calculate the velocity and acceleration
        mRotationTool.RecoverDisplacements(rModelPart); // !!

        // Updating time derivatives (nodally for efficiency)
        block_for_each(rModelPart.Nodes(), [&](Node& rNode)
        {
            // In MPM the delta_displacement is the current displacement as the previous_displacement is always reset
            const array_1d<double, 3 > & r_current_velocity = rNode.FastGetSolutionStepValue(VELOCITY); // !!
            const array_1d<double, 3 > & r_previous_velocity = rNode.FastGetSolutionStepValue(VELOCITY, 1);
            
            array_1d<double, 3 > & r_current_displacement = rNode.FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3 > & r_previous_displacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);

            r_current_displacement = r_previous_displacement + r_current_velocity * mDeltaTime;
            //KRATOS_WATCH(r_previous_velocity);
            //KRATOS_WATCH(r_previous_displacement);
            array_1d<double, 3>& r_current_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION);
            

            if (mIsDynamic){
                r_current_acceleration = (r_current_velocity) / mDeltaTime;
            }
        });

        KRATOS_CATCH( "" )
    }

    /**
     * @brief Performing the prediction of the solution
     * @details It predicts the solution for the current step x = xold + vold * Dt
     * @param rModelPart The model of the problem to solve
     * @param rDofSet set of all primary variables
     * @param rA LHS matrix
     * @param rDx Incremental update of primary variables
     * @param rb RHS Vector
     */
    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        KRATOS_TRY;

        // Start of TEMP: These are moved here to comply to pr #13432
        // Particle to Grid mapping for elements is done here because predict needs the velocity field (PR #13432)        
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const double delta_time = r_current_process_info[DELTA_TIME];

        // Loop over the grid nodes performed to clear all nodal information
        block_for_each(rModelPart.Nodes(), [&](Node& rNode)
		{
            // Variables to be cleaned
            double & r_nodal_mass     = rNode.FastGetSolutionStepValue(NODAL_MASS);
            array_1d<double, 3 > & r_nodal_momentum = rNode.FastGetSolutionStepValue(NODAL_MOMENTUM);
            array_1d<double, 3 > & r_nodal_inertia  = rNode.FastGetSolutionStepValue(NODAL_INERTIA);

            array_1d<double, 3 > & r_nodal_displacement = rNode.FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3 > & r_nodal_velocity     = rNode.FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & r_nodal_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION,1);

            double & r_nodal_old_pressure = rNode.FastGetSolutionStepValue(PRESSURE,1);
            double & r_nodal_pressure = rNode.FastGetSolutionStepValue(PRESSURE);

            // Clear
            r_nodal_mass = 0.0;
            r_nodal_momentum.clear();
            r_nodal_inertia.clear();

            r_nodal_displacement.clear();
            r_nodal_velocity.clear();
            r_nodal_acceleration.clear();
            r_nodal_old_pressure = 0.0;
            r_nodal_pressure = 0.0;

            // Other additional variables
            if (rNode.SolutionStepsDataHas(NODAL_AREA)){
                double & r_nodal_area = rNode.FastGetSolutionStepValue(NODAL_AREA);
                r_nodal_area          = 0.0;
            }
            if(rNode.SolutionStepsDataHas(NODAL_MPRESSURE)) {
                double & r_nodal_mpressure = rNode.FastGetSolutionStepValue(NODAL_MPRESSURE);
                r_nodal_mpressure = 0.0;
            }
            
            /*
            // friction-related
            if(mFrictionIsActive){
                rNode.FastGetSolutionStepValue(STICK_FORCE).clear();
                rNode.FastGetSolutionStepValue(FRICTION_STATE) = mRotationTool.GetSlidingState();
                rNode.SetValue(FRICTION_ASSIGNED, false);
            }
            */
		});

        // Extrapolate from Material Point Elements and Conditions (P2G Mapping)
        const auto &r_elements_array = rModelPart.Elements();
        const std::size_t n_elems = r_elements_array.size();
        IndexPartition<std::size_t>(n_elems).for_each([&](std::size_t i_elem) {
            auto it_elem = r_elements_array.begin() + i_elem;

            it_elem->AddExplicitContribution(r_current_process_info);
        });

        // Assign nodal variables after extrapolation
        block_for_each(rModelPart.Nodes(), [&](Node& rNode)
        {
            const double & r_nodal_mass = rNode.FastGetSolutionStepValue(NODAL_MASS);

            if (r_nodal_mass > std::numeric_limits<double>::epsilon())
            {
                const array_1d<double, 3 > & r_nodal_momentum   = rNode.FastGetSolutionStepValue(NODAL_MOMENTUM);
                const array_1d<double, 3 > & r_nodal_inertia    = rNode.FastGetSolutionStepValue(NODAL_INERTIA);

                array_1d<double, 3 > & r_previous_velocity     = rNode.FastGetSolutionStepValue(VELOCITY,1);
                array_1d<double, 3 > & r_previous_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION,1);
                double & r_previous_pressure = rNode.FastGetSolutionStepValue(PRESSURE,1);

                double delta_nodal_pressure = 0.0;

                // For mixed formulation
                if (rNode.HasDofFor(PRESSURE) && rNode.SolutionStepsDataHas(NODAL_MPRESSURE))
                {
                    double & nodal_mpressure = rNode.FastGetSolutionStepValue(NODAL_MPRESSURE);
                    delta_nodal_pressure = nodal_mpressure/r_nodal_mass;
                }

                const array_1d<double, 3 > delta_nodal_velocity = r_nodal_momentum/r_nodal_mass;
                const array_1d<double, 3 > delta_nodal_acceleration = r_nodal_inertia/r_nodal_mass;

                r_previous_velocity += delta_nodal_velocity;
                r_previous_acceleration += delta_nodal_acceleration;

                r_previous_pressure += delta_nodal_pressure;
                
                /*
                // mark nodes which have non-zero momentum in the 1st timestep s.t. these nodes can have
                // an initial friction state of SLIDING instead of STICK
                if(mFrictionIsActive){
                    const bool has_initial_momentum = (mGridModelPart.GetProcessInfo()[STEP] ==  1 && norm_2(r_nodal_momentum) > std::numeric_limits<double>::epsilon());
                    rNode.SetValue(HAS_INITIAL_MOMENTUM, has_initial_momentum);
                }*/
            }
        });
        // End of TEMP: These are moved here to comply to pr #13432

		block_for_each(rModelPart.Nodes(), [&](Node& rNode)
		{
            const array_1d<double, 3 > & r_previous_displacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
			const array_1d<double, 3 > & r_previous_velocity     = rNode.FastGetSolutionStepValue(VELOCITY, 1);
            const array_1d<double, 3 > & r_previous_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

            array_1d<double, 3 > & r_current_velocity  = rNode.FastGetSolutionStepValue(VELOCITY);

            // Displacement prediction for implicit MPM
            if (!(rNode.pGetDof(VELOCITY_X)->IsFixed()))
                {
                //KRATOS_INFO("Velocity x  fixed ---------------------------------------------");
                r_current_velocity[0] = 0.0;
                }
            else
                r_current_velocity[0]  = r_previous_velocity[0]; // !!
                //r_current_velocity[0] = 0.0;

            if (!(rNode.pGetDof(VELOCITY_Y)->IsFixed()))
                {
                //KRATOS_INFO("Velocity y  fixed ---------------------------------------------");
                r_current_velocity[1] = 0.0;
                }
            else
                r_current_velocity[1]  = r_previous_velocity[1];
                //r_current_velocity[1] = 0.0;

            if (rNode.HasDofFor(VELOCITY_Z))
            {
                if (!(rNode.pGetDof(VELOCITY_Z)->IsFixed()))
                    r_current_velocity[2] = 0.0;
                else
                    r_current_velocity[2]  = r_previous_velocity[2];
                    //r_current_velocity[2] = 0.0;
            }

            // Pressure prediction for implicit MPM
            if (rNode.HasDofFor(PRESSURE))
            {
                double& r_current_pressure        = rNode.FastGetSolutionStepValue(PRESSURE);
                const double& r_previous_pressure = rNode.FastGetSolutionStepValue(PRESSURE, 1);

                if (!(rNode.pGetDof(PRESSURE))->IsFixed())
                    r_current_pressure = r_previous_pressure;
            }

            // Updating time derivatives
            //array_1d<double, 3 > & current_velocity       = rNode.FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & current_acceleration   = rNode.FastGetSolutionStepValue(ACCELERATION);

            if (mIsDynamic){
                current_acceleration = (r_current_velocity) / mDeltaTime; // !! predizione di accelerazione nulla
            }

		});

        KRATOS_CATCH( "" );
    }

    void FinalizeNonLinIteration(ModelPart &rModelPart, TSystemMatrixType &rA, TSystemVectorType &rDx,
                                   TSystemVectorType &rb) override {

        // Special treatment of particle based dirichlet conditions to calculate the reaction forces at the boundary particles
        // ***
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        
        // clear any nodal reaction values
        ClearReactionVariable();
        
        // Calculating as an intermediate step the nodal reaction forces due to the boundary particles
        block_for_each(rModelPart.Conditions(), std::vector<bool>(), [&r_current_process_info](Condition& rCondition, auto& r_dummy)
        {  
            rCondition.CalculateOnIntegrationPoints(MPC_CALCULATE_NODAL_REACTIONS, r_dummy, r_current_process_info);
        });
        
        // Calculating the reaction forces at the boundary particles due to the nodal reaction forces
        block_for_each(rModelPart.Conditions(), std::vector<bool>(), [&r_current_process_info](Condition& rCondition, auto& r_dummy)
        {  
            rCondition.CalculateOnIntegrationPoints(MPC_CALCULATE_CONTACT_FORCE, r_dummy, r_current_process_info);
        });  

        // clear nodal reaction values again
        ClearReactionVariable();    
        
        // *** 

        ImplicitBaseType::FinalizeNonLinIteration(rModelPart, rA, rDx, rb);

    }

    void InitializeNonLinIteration(ModelPart &rModelPart, TSystemMatrixType &rA, TSystemVectorType &rDx,
                                 TSystemVectorType &rb) override {

        ImplicitBaseType::InitializeNonLinIteration(rModelPart, rA, rDx, rb);

        // clear nodal reaction values again
        ClearReactionVariable();
        /*
        if(mFrictionIsActive) {
            mRotationTool.ComputeFrictionAndResetFlags(rModelPart);
        }
        */
    }

    /**
     * @brief It initializes time step solution for MPM simulations.
     * @details The initialize solution step here also perform the following procedures:
     * 1. Loop over the grid nodes performed to clear all historical nodal information
     * 2. Assign a new nodal variables by means of extrapolation from material points
     * @param rModelPart The model of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        KRATOS_TRY
        
        // Particle to Grid mapping for elements is moved to predict because it needs the velocity field (PR #13432)
        ImplicitBaseType::InitializeSolutionStep(rModelPart,rA,rDx,rb);

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        mDeltaTime = r_current_process_info[DELTA_TIME];

        KRATOS_CATCH( "" )
    }

    /**
     * @brief Function called once at the end of a solution step, after convergence is reached if an iterative process is needed
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        ImplicitBaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);
        /*
        if(mFrictionIsActive) {
            block_for_each(mGridModelPart.Nodes(), [&](Node& rNode)
            {
                const Node& rConstNode = rNode; // const Node reference to avoid issues with previously unset GetValue()
                if( mRotationTool.IsConformingSlip(rConstNode) && rConstNode.GetValue(FRICTION_COEFFICIENT) > 0 )
                    rNode.FastGetSolutionStepValue(REACTION).clear();       
            });
            
            mRotationTool.ComputeFrictionAndResetFlags(rModelPart);
        }

        block_for_each(mGridModelPart.Nodes(), [&](Node& rNode) {
            const Node& rConstNode = rNode; // const Node reference to avoid issues with previously unset GetValue()

            // rotate forces stored in REACTION to global coordinates on conforming boundaries
            if (mRotationTool.IsConformingSlip(rConstNode) ) {
                mRotationTool.RotateVector(rNode.FastGetSolutionStepValue(REACTION), rConstNode, true);
            }
        });
        */

        //!! Qui va calcolato il displacement
        
    }

    void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& M,
        const ProcessInfo& rCurrentProcessInfo
        ) 
    {
        // Adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
            noalias(LHS_Contribution) += M * 1 / mDeltaTime; // * 1 / deltat
        // !! Controllare che la mass matrix sia sommata nel posto giusto
    }

    void AddDynamicsToRHS(
        Element& rElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& M,
        const ProcessInfo& rCurrentProcessInfo
        ) 
    {
        const auto& r_const_elem_ref = rElement;
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (M.size1() != 0) {
            LocalSystemVectorType vel; // !! come funziona questo in parallelo?
            rElement.GetFirstDerivativesVector(vel);
            //KRATOS_INFO("---------------RHS:");
            //KRATOS_WATCH(vel);
            //KRATOS_WATCH(mDeltaTime);
            vel *= - 1 / mDeltaTime; // !! negativo

            LocalSystemVectorType acc;
            //rElement.GetSecondDerivativesVector(acc);
            //KRATOS_WATCH(acc);
            acc += vel;

            RHS_Contribution -= prod(M, vel);
        }
    }

    void AddDynamicsToRHSConditions(
        Condition& rCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& M,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        const auto& r_const_cond_ref = rCondition;
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (M.size1() != 0) {
            LocalSystemVectorType vel; // !! come funziona questo in parallelo?
            rCondition.GetFirstDerivativesVector(vel);
            vel *= - 1 / mDeltaTime; // !! negativo

            LocalSystemVectorType acc;
            //rCondition.GetSecondDerivativesVector(acc);

            acc += vel;

            noalias(RHS_Contribution) -= prod(M, acc);
        }
    }



    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @param rCurrentElement The element to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY
        // !! questa funzione è presa da residual_based_bossak_displacement_scheme.hpp e opera in parallelo; qui funziona?
        const IndexType this_thread = OpenMPUtils::ThisThread();
        const auto& rConstElemRef = rCurrentElement;
        rCurrentElement.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);
        rConstElemRef.EquationIdVector(EquationId,rCurrentProcessInfo);

        if(mIsDynamic)
        {
            // !! 
            rCurrentElement.CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);
            AddDynamicsToLHS(LHS_Contribution, mMatrix.M[this_thread], rCurrentProcessInfo);
            //AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMatrix.M[this_thread], rCurrentProcessInfo);
        }

        /*
        // Rotate contributions (to match coordinates for slip conditions)
        if(!mRotationTool.IsParticleBasedSlip(rCurrentElement.GetGeometry())){
            // prevent rotation in case of particle-based slip (handled by condition itself)
            mRotationTool.Rotate(LHS_Contribution, RHS_Contribution, rCurrentElement.GetGeometry());
            mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentElement.GetGeometry());
        }
        */

        KRATOS_CATCH( "" )
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rCurrentElement The element to compute
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const IndexType this_thread = OpenMPUtils::ThisThread();
        const auto& r_const_elem_ref = rCurrentElement;

        // Basic operations for the element considered
        rCurrentElement.CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);
        r_const_elem_ref.EquationIdVector(EquationId,rCurrentProcessInfo);

        if(mIsDynamic)
        {
            // !!
            rCurrentElement.CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);
            //AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMatrix.M[this_thread], rCurrentProcessInfo);
        }

        /*
        // Rotate contributions (to match coordinates for slip conditions)
        if(!mRotationTool.IsParticleBasedSlip(rCurrentElement.GetGeometry())){
            // prevent rotation in case of particle-based slip (handled by condition itself)
            mRotationTool.Rotate(RHS_Contribution, rCurrentElement.GetGeometry());
            mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentElement.GetGeometry());
        }
        */

        KRATOS_CATCH( "" )
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCurrentCondition The condition to compute
     * @param rLHSContribution The LHS matrix contribution
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {

        KRATOS_TRY

        const IndexType this_thread = OpenMPUtils::ThisThread();
        const auto& r_const_cond_ref = rCurrentCondition;

        rCurrentCondition.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);
        r_const_cond_ref.EquationIdVector(EquationId,rCurrentProcessInfo);

        if(mIsDynamic)
        {
            // !! Da cambiare
            rCurrentCondition.CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);
            AddDynamicsToLHS(LHS_Contribution, mMatrix.M[this_thread], rCurrentProcessInfo);
            //AddDynamicsToRHSConditions(rCurrentCondition, RHS_Contribution, mMatrix.M[this_thread], rCurrentProcessInfo);
        }
        /*
        // Rotate contributions (to match coordinates for slip conditions)
        if(!mRotationTool.IsParticleBasedSlip(rCurrentCondition.GetGeometry())){
            // prevent rotation in case of particle-based slip (handled by condition itself)
            mRotationTool.Rotate(LHS_Contribution, RHS_Contribution, rCurrentCondition.GetGeometry());
            mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentCondition.GetGeometry());
        }
        */

        KRATOS_CATCH( "" )
    }

    /**
     * @brief Functions that calculates the RHS of a "condition" object
     * @param rCurrentCondition The condition to compute
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const IndexType this_thread = OpenMPUtils::ThisThread();
        const auto& r_const_cond_ref = rCurrentCondition;
        rCurrentCondition.CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);
        r_const_cond_ref.EquationIdVector(EquationId,rCurrentProcessInfo);

        if(mIsDynamic)
        {
            rCurrentCondition.CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);
            //AddDynamicsToRHSConditions(rCurrentCondition, RHS_Contribution, mMatrix.M[this_thread], rCurrentProcessInfo);
        }
        /*
        // Rotate contributions (to match coordinates for slip conditions)
        if(!mRotationTool.IsParticleBasedSlip(rCurrentCondition.GetGeometry())){
            // prevent rotation in case of particle-based slip (handled by condition itself)
            mRotationTool.Rotate(RHS_Contribution, rCurrentCondition.GetGeometry());
            mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentCondition.GetGeometry());
        }
        */

        KRATOS_CATCH( "" )
    }

protected:
    // MPM Background Grid
    ModelPart& mGridModelPart;

    // To distinguish quasi-static and dynamic
    bool mIsDynamic;

    // Identifies cases where friction is active in at least 1 slip boundary
    bool mFrictionIsActive;

    double mDeltaTime;

    // For Rotation Utility
    unsigned int mDomainSize;
    unsigned int mBlockSize;
    MPMBoundaryRotationUtility<LocalSystemMatrixType,LocalSystemVectorType> mRotationTool;

    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();;

    void ClearReactionVariable() const
    {
        block_for_each(mGridModelPart.Nodes(), [&](Node& rNode)
        {
            rNode.FastGetSolutionStepValue(REACTION).clear();
        });
    }

}; /* Class MPMResidualBasedBossakScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_MPM_RESIDUAL_BASED_BOSSAK_SCHEME defined */

