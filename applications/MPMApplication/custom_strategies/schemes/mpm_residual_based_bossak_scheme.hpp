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


#if !defined(KRATOS_MPM_RESIDUAL_BASED_BOSSAK_SCHEME )
#define      KRATOS_MPM_RESIDUAL_BASED_BOSSAK_SCHEME

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
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "custom_utilities/mpm_boundary_rotation_utility.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_base_condition.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

/**
 * @class MPMResidualBasedBossakScheme
 * @ingroup KratosMPM
 * @brief Bossak integration scheme (for linear and nonlinear dynamic problems) for displacements adjusted for Material Point Method
 * @details This is an implicit scheme based of the Bossak algorithm for displacements suitable for quasi-static and dynamic problems.
 * Furthermore, this scheme has been adjusted for mixed formulation where pressure is also solved as one of the DoFs.
 * The parameter Alpha of Bossak introduces damping, the value of Bossak is from 0 to -0.5 (negative)
 * Implementation according to: "An alpha modification of Newmark's method; W.L. Wood, M. Bossak, O.C. Zienkiewicz;
 * Numerical Methods in Engineering; 1980"
 * MPM implementation according to: "Implicit time integration for the material point method"; J. E. Guilkey, J. A. Weiss
 * The parameter Alpha of Bossak introduces damping, the value of Bossak is from 0 to -0.5 (negative)
 */
template<class TSparseSpace,  class TDenseSpace >
class MPMResidualBasedBossakScheme
    : public ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( MPMResidualBasedBossakScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                                      BaseType;

    typedef ResidualBasedImplicitTimeScheme<TSparseSpace,TDenseSpace>     ImplicitBaseType;

    typedef ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace> BossakBaseType;

    typedef typename BossakBaseType::TDataType                                   TDataType;

    typedef typename BossakBaseType::DofsArrayType                           DofsArrayType;

    typedef typename Element::DofsVectorType                                DofsVectorType;

    typedef typename BossakBaseType::TSystemMatrixType                   TSystemMatrixType;

    typedef typename BossakBaseType::TSystemVectorType                   TSystemVectorType;

    typedef typename BossakBaseType::LocalSystemVectorType           LocalSystemVectorType;

    typedef typename BossakBaseType::LocalSystemMatrixType           LocalSystemMatrixType;

    typedef ModelPart::ElementsContainerType                             ElementsArrayType;

    typedef ModelPart::ConditionsContainerType                         ConditionsArrayType;

    typedef typename BaseType::Pointer                                     BaseTypePointer;

    using BossakBaseType::mpDofUpdater;

    using BossakBaseType::mBossak;

    using ImplicitBaseType::mMatrix;

    /**
     * Constructor.
     * The MPM bossak scheme
     */
    /**
     * @brief Constructor.
     * @detail The MPM bossak method
     * @param rGridModelPart The background grid model part
     * @param DomainSize Domain size or space dimension of the problems
     * @param BlockSize Block size of DOF
     * @param Alpha is the Bossak parameter. Default value is 0, which is the Newmark method
     * @param NewmarkBeta the Newmark parameter. Default value is 0.25, for mean constant acceleration.
     */
    MPMResidualBasedBossakScheme(ModelPart& rGridModelPart, unsigned int DomainSize,
        unsigned int BlockSize, double Alpha = 0.0,
        double NewmarkBeta = 0.25, bool IsDynamic = true)
                :ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>(Alpha, NewmarkBeta),
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
     MPMResidualBasedBossakScheme(MPMResidualBasedBossakScheme& rOther)
        :BossakBaseType(rOther)
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
        return BaseTypePointer( new MPMResidualBasedBossakScheme(*this) );
    }

    /** Destructor.
     */
    virtual ~MPMResidualBasedBossakScheme
    () {}


    //***************************************************************************
    //***************************************************************************

    void Initialize(ModelPart& rModelPart) override {
        MPMParticleBaseCondition::SetRotationUtility(&mRotationTool);

        BossakBaseType::Initialize(rModelPart);
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
        mRotationTool.RotateDisplacements(rModelPart);

        // Update of displacement (by DOF)
        mpDofUpdater->UpdateDofs(rDofSet, rDx);

        // Rotate the displacement back to the original coordinate system to calculate the velocity and acceleration
        mRotationTool.RecoverDisplacements(rModelPart);

        // Updating time derivatives (nodally for efficiency)
        block_for_each(rModelPart.Nodes(), [&](Node& rNode)
        {
            // In MPM the delta_displacement is the current displacement as the previous_displacement is always reset
            const array_1d<double, 3 > & r_delta_displacement = rNode.FastGetSolutionStepValue(DISPLACEMENT);

            array_1d<double, 3>& r_current_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& r_previous_velocity = rNode.FastGetSolutionStepValue(VELOCITY, 1);

            array_1d<double, 3>& r_current_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION);
            const array_1d<double, 3>& r_previous_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

            if (mIsDynamic){
                BossakBaseType::UpdateVelocity(r_current_velocity, r_delta_displacement, r_previous_velocity, r_previous_acceleration);
                BossakBaseType::UpdateAcceleration(r_current_acceleration, r_delta_displacement, r_previous_velocity, r_previous_acceleration);
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

		block_for_each(rModelPart.Nodes(), [&](Node& rNode)
		{
            const array_1d<double, 3 > & r_previous_displacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
			const array_1d<double, 3 > & r_previous_velocity     = rNode.FastGetSolutionStepValue(VELOCITY, 1);
            const array_1d<double, 3 > & r_previous_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

            array_1d<double, 3 > & r_current_displacement  = rNode.FastGetSolutionStepValue(DISPLACEMENT);

            // Displacement prediction for implicit MPM
            if (!(rNode.pGetDof(DISPLACEMENT_X)->IsFixed()))
                r_current_displacement[0] = 0.0;
            else
                r_current_displacement[0]  = r_previous_displacement[0];

            if (!(rNode.pGetDof(DISPLACEMENT_Y)->IsFixed()))
                r_current_displacement[1] = 0.0;
            else
                r_current_displacement[1]  = r_previous_displacement[1];

            if (rNode.HasDofFor(DISPLACEMENT_Z))
            {
                if (!(rNode.pGetDof(DISPLACEMENT_Z)->IsFixed()))
                    r_current_displacement[2] = 0.0;
                else
                    r_current_displacement[2]  = r_previous_displacement[2];
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
            array_1d<double, 3 > & current_velocity       = rNode.FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & current_acceleration   = rNode.FastGetSolutionStepValue(ACCELERATION);

            if (mIsDynamic){
                BossakBaseType::UpdateVelocity(current_velocity, r_current_displacement, r_previous_velocity, r_previous_acceleration);
                BossakBaseType::UpdateAcceleration (current_acceleration, r_current_displacement, r_previous_velocity, r_previous_acceleration);
            }

		});

        KRATOS_CATCH( "" );
    }

    void FinalizeNonLinIteration(ModelPart &rModelPart, TSystemMatrixType &rA, TSystemVectorType &rDx,
                                   TSystemVectorType &rb) override {

        // clear any nodal reaction values for conforming friction (needed in current penalty-based formulation)
        ClearConformingFrictionReaction();

        BossakBaseType::FinalizeNonLinIteration(rModelPart, rA, rDx, rb);

        // modify reaction forces for material point particle slip conditions (Penalty)
        mRotationTool.CalculateReactionForces(mGridModelPart);

        if(mFrictionIsActive) {
            mRotationTool.ComputeFrictionAndResetFlags(rModelPart);
        }
    }

    void InitializeNonLinIteration(ModelPart &rModelPart, TSystemMatrixType &rA, TSystemVectorType &rDx,
                                 TSystemVectorType &rb) override {

        BossakBaseType::InitializeNonLinIteration(rModelPart, rA, rDx, rb);

        if(mFrictionIsActive) {
            mRotationTool.ComputeFrictionAndResetFlags(rModelPart);
        }
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

        // Loop over the grid nodes performed to clear all nodal information
        block_for_each(rModelPart.Nodes(), [&](Node& rNode)
		{
            // Variables to be cleaned
            double & r_nodal_mass     = rNode.FastGetSolutionStepValue(NODAL_MASS);
            array_1d<double, 3 > & r_nodal_momentum = rNode.FastGetSolutionStepValue(NODAL_MOMENTUM);
            array_1d<double, 3 > & r_nodal_inertia  = rNode.FastGetSolutionStepValue(NODAL_INERTIA);

            array_1d<double, 3 > & r_nodal_displacement = rNode.FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 > & r_nodal_velocity     = rNode.FastGetSolutionStepValue(VELOCITY,1);
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

            // friction-related
            if(mFrictionIsActive){
                rNode.FastGetSolutionStepValue(STICK_FORCE).clear();
                rNode.FastGetSolutionStepValue(FRICTION_STATE) = mRotationTool.GetSlidingState();
                rNode.SetValue(FRICTION_ASSIGNED, false);
            }
		});

        // Extrapolate from Material Point Elements and Conditions
        ImplicitBaseType::InitializeSolutionStep(rModelPart,rA,rDx,rb);

        // Assign nodal variables after extrapolation
        block_for_each(rModelPart.Nodes(), [&](Node& rNode)
        {
            const double & r_nodal_mass = rNode.FastGetSolutionStepValue(NODAL_MASS);

            if (r_nodal_mass > std::numeric_limits<double>::epsilon())
            {
                const array_1d<double, 3 > & r_nodal_momentum   = rNode.FastGetSolutionStepValue(NODAL_MOMENTUM);
                const array_1d<double, 3 > & r_nodal_inertia    = rNode.FastGetSolutionStepValue(NODAL_INERTIA);

                array_1d<double, 3 > & r_nodal_velocity     = rNode.FastGetSolutionStepValue(VELOCITY,1);
                array_1d<double, 3 > & r_nodal_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION,1);
                double & r_nodal_pressure = rNode.FastGetSolutionStepValue(PRESSURE,1);

                double delta_nodal_pressure = 0.0;

                // For mixed formulation
                if (rNode.HasDofFor(PRESSURE) && rNode.SolutionStepsDataHas(NODAL_MPRESSURE))
                {
                    double & nodal_mpressure = rNode.FastGetSolutionStepValue(NODAL_MPRESSURE);
                    delta_nodal_pressure = nodal_mpressure/r_nodal_mass;
                }

                const array_1d<double, 3 > delta_nodal_velocity = r_nodal_momentum/r_nodal_mass;
                const array_1d<double, 3 > delta_nodal_acceleration = r_nodal_inertia/r_nodal_mass;

                r_nodal_velocity += delta_nodal_velocity;
                r_nodal_acceleration += delta_nodal_acceleration;

                r_nodal_pressure += delta_nodal_pressure;

                // mark nodes which have non-zero momentum in the 1st timestep s.t. these nodes can have
                // an initial friction state of SLIDING instead of STICK
                if(mFrictionIsActive){
                    const bool has_initial_momentum = (mGridModelPart.GetProcessInfo()[STEP] ==  1 && norm_2(r_nodal_momentum) > std::numeric_limits<double>::epsilon());
                    rNode.SetValue(HAS_INITIAL_MOMENTUM, has_initial_momentum);
                }
            }
        });

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const double delta_time = r_current_process_info[DELTA_TIME];

        // Initializing Bossak constants
        mBossak.c0 = ( 1.0 / (mBossak.beta * delta_time * delta_time) );
        mBossak.c1 = ( mBossak.gamma / (mBossak.beta * delta_time) );
        mBossak.c2 = ( 1.0 / (mBossak.beta * delta_time) );
        mBossak.c3 = ( 0.5 / (mBossak.beta) - 1.0 );
        mBossak.c4 = ( (mBossak.gamma / mBossak.beta) - 1.0  );
        mBossak.c5 = ( delta_time * 0.5 * ( ( mBossak.gamma / mBossak.beta ) - 2.0 ) );

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
        BossakBaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        if(mFrictionIsActive) {
            block_for_each(mGridModelPart.Nodes(), [&](Node& rNode) {
                const Node& rConstNode = rNode; // const Node reference to avoid issues with previously unset GetValue()

                const double mu = rConstNode.GetValue(FRICTION_COEFFICIENT);

                // rotate friction forces stored in REACTION to global coordinates on conforming boundaries
                if (mRotationTool.IsConformingSlip(rNode) && mu > 0) {
                    mRotationTool.RotateVector(rNode.FastGetSolutionStepValue(REACTION), rNode, true);
                }
            });
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

        const IndexType this_thread = OpenMPUtils::ThisThread();
        const auto& rConstElemRef = rCurrentElement;
        rCurrentElement.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);
        rConstElemRef.EquationIdVector(EquationId,rCurrentProcessInfo);

        if(mIsDynamic)
        {
            rCurrentElement.CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);
            rCurrentElement.CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);
            BossakBaseType::AddDynamicsToLHS(LHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);
            BossakBaseType::AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);
        }

        // Rotate contributions (to match coordinates for slip conditions)
        if(!mRotationTool.IsParticleBasedSlip(rCurrentElement.GetGeometry())){
            // prevent rotation in case of particle-based slip (handled by condition itself)
            mRotationTool.Rotate(LHS_Contribution, RHS_Contribution, rCurrentElement.GetGeometry());
            mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentElement.GetGeometry());
        }


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
            rCurrentElement.CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);
            rCurrentElement.CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);
            BossakBaseType::AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);
        }

        // Rotate contributions (to match coordinates for slip conditions)
        if(!mRotationTool.IsParticleBasedSlip(rCurrentElement.GetGeometry())){
            // prevent rotation in case of particle-based slip (handled by condition itself)
            mRotationTool.Rotate(RHS_Contribution, rCurrentElement.GetGeometry());
            mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentElement.GetGeometry());
        }

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
            rCurrentCondition.CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);
            rCurrentCondition.CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);
            BossakBaseType::AddDynamicsToLHS(LHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);
            BossakBaseType::AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);
        }

        // Rotate contributions (to match coordinates for slip conditions)
        if(!mRotationTool.IsParticleBasedSlip(rCurrentCondition.GetGeometry())){
            // prevent rotation in case of particle-based slip (handled by condition itself)
            mRotationTool.Rotate(LHS_Contribution, RHS_Contribution, rCurrentCondition.GetGeometry());
            mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentCondition.GetGeometry());
        }


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
            rCurrentCondition.CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);
            BossakBaseType::AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);
        }

        // Rotate contributions (to match coordinates for slip conditions)
        if(!mRotationTool.IsParticleBasedSlip(rCurrentCondition.GetGeometry())){
            // prevent rotation in case of particle-based slip (handled by condition itself)
            mRotationTool.Rotate(RHS_Contribution, rCurrentCondition.GetGeometry());
            mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentCondition.GetGeometry());
        }


        KRATOS_CATCH( "" )
    }

protected:

    // MPM Background Grid
    ModelPart& mGridModelPart;

    // To distinguish quasi-static and dynamic
    bool mIsDynamic;

    // Identifies cases where friction is active in at least 1 slip boundary
    bool mFrictionIsActive;

    // For Rotation Utility
    unsigned int mDomainSize;
    unsigned int mBlockSize;
    MPMBoundaryRotationUtility<LocalSystemMatrixType,LocalSystemVectorType> mRotationTool;

    void ClearConformingFrictionReaction() const
    {
        block_for_each(mGridModelPart.Nodes(), [&](Node& rNode)
        {
            const Node& rConstNode = rNode; // const Node reference to avoid issues with previously unset GetValue()

            if( mRotationTool.IsConformingSlip(rConstNode) && rConstNode.GetValue(FRICTION_COEFFICIENT) > 0 )
                rNode.FastGetSolutionStepValue(REACTION).clear();
        });
    }

}; /* Class MPMResidualBasedBossakScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_MPM_RESIDUAL_BASED_BOSSAK_SCHEME defined */

