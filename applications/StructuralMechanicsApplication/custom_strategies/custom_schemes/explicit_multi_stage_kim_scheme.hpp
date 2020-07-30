//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Klaus B Sautter
//  Based on : "An accurate two‐stage explicit time integration scheme for structural dynamics and various dynamic problems" - Wooram Kim
//

#if !defined(KRATOS_EXPLICIT_MULTI_STAGE_KIM_SCHEME_HPP_INCLUDED)
#define KRATOS_EXPLICIT_MULTI_STAGE_KIM_SCHEME_HPP_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/scheme.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/explicit_integration_utilities.h"

namespace Kratos {

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
 * @class ExplicitMultiStageKimScheme
 * @ingroup StructuralMechanicsApplciation
 * @brief An explicit multi stage scheme
 * @details "An accurate two stage explicit time integration scheme for
 * structural dynamics and various dynamic problems" - Wooram Kim
 * DOI: 10.1002/nme.6098
 * @author Klaus B Sautter
 */
template <class TSparseSpace,
          class TDenseSpace //= DenseSpace<double>
          >
class ExplicitMultiStageKimScheme
    : public Scheme<TSparseSpace, TDenseSpace> {

public:
    ///@name Type Definitions
    ///@{

    /// The definition of the base type
    typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

    /// Some definitions related with the base class
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    /// The arrays of elements and nodes
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesArrayType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Definition fo the node iterator
    typedef typename ModelPart::NodeIterator NodeIterator;

    /// The definition of the numerical limit
    static constexpr double numerical_limit = std::numeric_limits<double>::epsilon();

    typedef Variable<array_1d<double,3>> ArrayVarType;
    typedef Variable<double> DoubleVarType;
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double,3>>> VarComponentType;

    typedef array_1d<double, 3> Double3DArray;

    /// Counted pointer of ExplicitMultiStageKimScheme
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitMultiStageKimScheme);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details The ExplicitMultiStageKimScheme method
     * @param MaximumDeltaTime The maximum delta time to be considered
     * @param DeltaTimeFraction The delta ttime fraction
     * @param DeltaTimePredictionLevel The prediction level
     */
    ExplicitMultiStageKimScheme(
        const double DeltaTimeFraction
        )
        : Scheme<TSparseSpace, TDenseSpace>()
    {
        mDeltaTime.Fraction = DeltaTimeFraction;
    }

    /**
     * @brief Constructor with parameters
     * @details The ExplicitMultiStageKimScheme method
     * @param rParameters The parameters containing the configuration parameters
     * @warning time_step_prediction_level should be an integer
     */
    ExplicitMultiStageKimScheme(Parameters rParameters =  Parameters(R"({})"))
        : Scheme<TSparseSpace, TDenseSpace>()
    {
        Parameters default_parameters = Parameters(R"(
        {
            "fraction_delta_time"        : 0.333333333333333333333333333333333333
        })" );

        rParameters.ValidateAndAssignDefaults(default_parameters);
        mDeltaTime.Fraction = rParameters["fraction_delta_time"].GetDouble();
    }

    /** Destructor.
    */
    virtual ~ExplicitMultiStageKimScheme() {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided.
     * @details Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model of the problem to solve
     * @return Zero means  all ok
     */
    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        BaseType::Check(rModelPart);

        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2) << "Insufficient buffer size for Central Difference Scheme. It has to be > 2" << std::endl;

        KRATOS_ERROR_IF_NOT(rModelPart.GetProcessInfo().Has(DOMAIN_SIZE)) << "DOMAIN_SIZE not defined on ProcessInfo. Please define" << std::endl;

        KRATOS_ERROR_IF((mDeltaTime.Fraction <= 0.0) || (mDeltaTime.Fraction >= 1.0)) << "fraction_delta_time must be >0 && <1 !" << std::endl;

        KRATOS_ERROR_IF_NOT(rModelPart.NumberOfNodes()>0) << "model part contains 0 nodes" << std::endl;

        return 0;

        KRATOS_CATCH("");
    }

    /**
     * @brief This is the place to initialize the Scheme. This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model of the problem to solve
     */
    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Preparing the time values for the first step (where time = initial_time +
        // dt)
        mTime.Delta = r_current_process_info[DELTA_TIME];
        mTime.MidStep = mTime.Delta * mDeltaTime.Fraction;

        /// Working in 2D/3D (the definition of DOMAIN_SIZE is check in the Check method)
        const SizeType dim = r_current_process_info[DOMAIN_SIZE];

        // Initialize scheme
        if (!BaseType::SchemeIsInitialized())
            InitializeExplicitScheme(rModelPart, dim);
        else
            SchemeCustomInitialization(rModelPart, dim);

        BaseType::SetSchemeIsInitialized();

        KRATOS_CATCH("")
    }

    /**
     * @brief It initializes time step solution. Only for reasons if the time step solution is restarted
     * @param rModelPart The model of the problem to solve
     * @param rA LHS matrix
     * @param rDx Incremental update of primary variables
     * @param rb RHS Vector
     * @todo I cannot find the formula for the higher orders with variable time step. I tried to deduce by myself but the result was very unstable
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY
        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);
        InitializeResidual(rModelPart);
        KRATOS_CATCH("")
    }

    /**
     * @brief It initializes the non-linear iteration
     * @param rModelPart The model of the problem to solve
     * @param rA LHS matrix
     * @param rDx Incremental update of primary variables
     * @param rb RHS Vector
     * @todo I cannot find the formula for the higher orders with variable time step. I tried to deduce by myself but the result was very unstable
     */
    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY;

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        const auto it_elem_begin = rModelPart.ElementsBegin();
        #pragma omp parallel for schedule(guided,512)
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); ++i) {
            auto it_elem = it_elem_begin + i;
            it_elem->InitializeNonLinearIteration(r_current_process_info);
        }

        const auto it_cond_begin = rModelPart.ConditionsBegin();
        #pragma omp parallel for schedule(guided,512)
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); ++i) {
            auto it_elem = it_cond_begin + i;
            it_elem->InitializeNonLinearIteration(r_current_process_info);
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief This method initializes the residual in the nodes of the model part
     * @param rModelPart The model of the problem to solve
     */
    void InitializeResidual(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // The array of nodes
        NodesArrayType& r_nodes = rModelPart.Nodes();

        // Auxiliar values
        const Double3DArray zero_array = ZeroVector(3);
        // Initializing the variables
        VariableUtils().SetVariable(FORCE_RESIDUAL, zero_array,r_nodes);
        const bool has_dof_for_rot_z = (r_nodes.begin())->HasDofFor(ROTATION_Z);
        if (has_dof_for_rot_z)
            VariableUtils().SetVariable(MOMENT_RESIDUAL,zero_array,r_nodes);

        KRATOS_CATCH("")
    }

    /**
     * @brief This method initializes some rutines related with the explicit scheme
     * @param rModelPart The model of the problem to solve
     * @param DomainSize The current dimention of the problem
     */
    void InitializeExplicitScheme(
        ModelPart& rModelPart,
        const SizeType DomainSize = 3
        )
    {
        KRATOS_TRY

        /// The array of ndoes
        NodesArrayType& r_nodes = rModelPart.Nodes();

        // The first iterator of the array of nodes
        const auto it_node_begin = rModelPart.NodesBegin();

        /// Initialise the database of the nodes
        const Double3DArray zero_array = ZeroVector(3);
        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            auto it_node = (it_node_begin + i);
            it_node->SetValue(NODAL_MASS, 0.0);
            Double3DArray& r_fractional_acceleration = it_node->FastGetSolutionStepValue(FRACTIONAL_ACCELERATION);
            r_fractional_acceleration  = ZeroVector(3);
        }
        const bool has_dof_for_rot_z = it_node_begin->HasDofFor(ROTATION_Z);
        if (has_dof_for_rot_z) {
            #pragma omp parallel for schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                auto it_node = (it_node_begin + i);
                it_node->SetValue(NODAL_INERTIA, zero_array);
                Double3DArray& r_fractional_acceleration = it_node->FastGetSolutionStepValue(FRACTIONAL_ANGULAR_ACCELERATION);
                r_fractional_acceleration  = ZeroVector(3);
            }
        }

        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            auto it_node = (it_node_begin + i);

            Double3DArray& r_current_residual = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);

            for (IndexType j = 0; j < DomainSize; j++) {
                r_current_residual[j] = 0.0;}

            if (has_dof_for_rot_z) {
                Double3DArray& r_current_residual_moment = it_node->FastGetSolutionStepValue(MOMENT_RESIDUAL);

                const IndexType initial_j = DomainSize == 3 ? 0 : 2; // We do this because in 2D only the rotation Z is needed, then we start with 2, instead of 0
                for (IndexType j = initial_j; j < 3; j++) {
                    r_current_residual_moment[j] = 0.0;}
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Performing the update of the solution
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
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        // The current process info
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // The array of nodes
        NodesArrayType& r_nodes = rModelPart.Nodes();

        /// Working in 2D/3D (the definition of DOMAIN_SIZE is check in the Check method)
        const SizeType dim = r_current_process_info[DOMAIN_SIZE];

        // The iterator of the first node
        const auto it_node_begin = rModelPart.NodesBegin();
        const bool has_dof_for_rot_z = it_node_begin->HasDofFor(ROTATION_Z);


        // Getting dof position
        const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);
        const IndexType rotppos = has_dof_for_rot_z ? it_node_begin->GetDofPosition(ROTATION_X) : 0;

        // _____________________________________________________________________
        // ______________________________ STAGE 3 ______________________________
        // _____________________________________________________________________

        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            // Current step information "N+1" (before step update).
            this->UpdateAccelerationStage(
                it_node_begin + i, disppos,ACCELERATION,VELOCITY,
                NODAL_DISPLACEMENT_DAMPING,FORCE_RESIDUAL,NODAL_MASS,
                DISPLACEMENT_X,DISPLACEMENT_Y,DISPLACEMENT_Z,dim);
        } // for Node parallel

        if (has_dof_for_rot_z){
            #pragma omp parallel for schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                // Current step information "N+1" (before step update).
                this->UpdateAccelerationStage(
                    it_node_begin + i, rotppos,ANGULAR_ACCELERATION,
                    ANGULAR_VELOCITY,NODAL_ROTATION_DAMPING,MOMENT_RESIDUAL,
                    NODAL_INERTIA,ROTATION_X,ROTATION_Y,ROTATION_Z,dim);
            } // for Node parallel
        }


        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            // Current step information "N+1" (before step update).
            this->UpdateDegreesOfFreedomStage3(it_node_begin + i,
                VELOCITY,DISPLACEMENT,ACCELERATION,FRACTIONAL_ACCELERATION,dim);
        } // for Node parallel

        if (has_dof_for_rot_z){
            #pragma omp parallel for schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                // Current step information "N+1" (before step update).
                this->UpdateDegreesOfFreedomStage3(it_node_begin + i,
                    ANGULAR_VELOCITY,ROTATION,ANGULAR_ACCELERATION,
                    FRACTIONAL_ANGULAR_ACCELERATION,dim);
            } // for Node parallel
        }

        KRATOS_CATCH("")
    }


    template <typename TObjectType>
    void UpdateAccelerationStage(
        NodeIterator itCurrentNode,
        const IndexType FieldPosition,
        const ArrayVarType& rAccelerationVariable,
        const ArrayVarType& rVelocityVariable,
        const TObjectType& rDampingVariable,
        const ArrayVarType& rResidualVariable,
        const TObjectType& rIntertiaVariable,
        const VarComponentType& rFixVariable1,
        const VarComponentType& rFixVariable2,
        const VarComponentType& rFixVariable3,
        const SizeType DomainSize = 3
        )
    {
        Double3DArray& r_acceleration = itCurrentNode->FastGetSolutionStepValue(rAccelerationVariable);
        Double3DArray& r_current_velocity = itCurrentNode->FastGetSolutionStepValue(rVelocityVariable);
        const Double3DArray& r_current_residual = itCurrentNode->FastGetSolutionStepValue(rResidualVariable);

        if (std::is_same<TObjectType,DoubleVarType>::value) {
            const double& r_nodal_inertia = itCurrentNode->GetValue(NODAL_MASS);
            const double& r_nodal_damping = itCurrentNode->GetValue(NODAL_DISPLACEMENT_DAMPING);

            // Solution of the explicit equation:
            if (r_nodal_inertia > numerical_limit)
            {
                noalias(r_acceleration) = (r_current_residual - r_nodal_damping * r_current_velocity) / r_nodal_inertia;
            }
            else
                noalias(r_acceleration) = ZeroVector(3);
        }
        else if (std::is_same<TObjectType,ArrayVarType>::value) {
            const Double3DArray& r_nodal_inertia = itCurrentNode->GetValue(NODAL_INERTIA);
            const Double3DArray& r_nodal_rotational_damping = itCurrentNode->GetValue(NODAL_ROTATION_DAMPING);

            const IndexType initial_k = DomainSize == 3 ? 0 : 2; // We do this because in 2D only the rotation Z is needed, then we start with 2, instead of 0
            for (IndexType kk = initial_k; kk < 3; ++kk) {
                if (r_nodal_inertia[kk] > numerical_limit)
                    r_acceleration[kk] = (r_current_residual[kk] - r_nodal_rotational_damping[kk] * r_current_velocity[kk]) / r_nodal_inertia[kk];
                else
                    r_acceleration[kk] = 0.0;
            }
        }
        else KRATOS_ERROR << "cannot handle input variable in \"explicit_multi_stage_kim_scheme.hpp\"" << std::endl;

        std::array<bool, 3> fix_field = {false, false, false};
        fix_field[0] = (itCurrentNode->GetDof(rFixVariable1, FieldPosition).IsFixed());
        fix_field[1] = (itCurrentNode->GetDof(rFixVariable2, FieldPosition + 1).IsFixed());
        if (DomainSize == 3)
            fix_field[2] = (itCurrentNode->GetDof(rFixVariable3, FieldPosition + 2).IsFixed());

        for (IndexType j = 0; j < DomainSize; j++) {
            if (fix_field[j]) {
                r_acceleration[j] = 0.0;
            }
        }
    }

    /**
     * @brief This method updates the translation DoF
     * @param itCurrentNode The iterator of the current node
     * @param DisplacementPosition The position of the displacement dof on the database
     * @param DomainSize The current dimention of the problem
     */
    void UpdateDegreesOfFreedomStage1(
        NodeIterator itCurrentNode,
        const ArrayVarType& rVelocityVariable,
        const ArrayVarType& rDisplacementVariable,
        const ArrayVarType& rAccelerationVariable,
        const SizeType DomainSize = 3
        )
    {
        Double3DArray& r_current_velocity = itCurrentNode->FastGetSolutionStepValue(rVelocityVariable);
        Double3DArray& r_current_displacement = itCurrentNode->FastGetSolutionStepValue(rDisplacementVariable);

        const Double3DArray& r_previous_displacement = itCurrentNode->FastGetSolutionStepValue(rDisplacementVariable, 1);
        const Double3DArray& r_previous_velocity = itCurrentNode->FastGetSolutionStepValue(rVelocityVariable, 1);
        const Double3DArray& r_previous_acceleration = itCurrentNode->FastGetSolutionStepValue(rAccelerationVariable, 1);

        for (IndexType j = 0; j < DomainSize; j++) {

            r_current_displacement[j] = r_previous_displacement[j] + mTime.MidStep * r_previous_velocity[j];
            r_current_displacement[j] += 0.50 * mTime.MidStep * mTime.MidStep * r_previous_acceleration[j];

            r_current_velocity[j] =  r_previous_velocity[j] + mTime.MidStep * r_previous_acceleration[j];
        } // for DomainSize
    }

    void UpdateDegreesOfFreedomStage2(
        NodeIterator itCurrentNode,
        const ArrayVarType& rVelocityVariable,
        const ArrayVarType& rDisplacementVariable,
        const ArrayVarType& rAccelerationVariable,
        const ArrayVarType& rFractionalAccelerationVariable,
        const SizeType DomainSize = 3
        )
    {
        Double3DArray& r_current_velocity = itCurrentNode->FastGetSolutionStepValue(rVelocityVariable);
        Double3DArray& r_current_displacement = itCurrentNode->FastGetSolutionStepValue(rDisplacementVariable);
        Double3DArray& r_fractional_acceleration = itCurrentNode->FastGetSolutionStepValue(rFractionalAccelerationVariable);

        const Double3DArray& r_previous_displacement = itCurrentNode->FastGetSolutionStepValue(rDisplacementVariable, 1);
        const Double3DArray& r_previous_velocity = itCurrentNode->FastGetSolutionStepValue(rVelocityVariable, 1);
        const Double3DArray& r_previous_acceleration = itCurrentNode->FastGetSolutionStepValue(rAccelerationVariable, 1);


        for (IndexType j = 0; j < DomainSize; j++) {

            r_current_displacement[j] = r_previous_displacement[j] + mTime.Delta * r_previous_velocity[j];
            r_current_displacement[j] += 0.50 * mTime.Delta * mTime.Delta * r_previous_acceleration[j] * ((3.0*mDeltaTime.Fraction-1.0)/(3.0*mDeltaTime.Fraction));
            r_current_displacement[j] += 0.50 * mTime.Delta * mTime.Delta * r_fractional_acceleration[j] * (1.0/(3.0*mDeltaTime.Fraction));


            r_current_velocity[j] = r_previous_velocity[j];
            r_current_velocity[j] += mTime.Delta * r_previous_acceleration[j] * ((2.0*mDeltaTime.Fraction-1.0)/(2.0*mDeltaTime.Fraction));
            r_current_velocity[j] += mTime.Delta * r_fractional_acceleration[j] * ((1.0)/(2.0*mDeltaTime.Fraction));
        } // for DomainSize
    }

    void UpdateDegreesOfFreedomStage3(
        NodeIterator itCurrentNode,
        const ArrayVarType& rVelocityVariable,
        const ArrayVarType& rDisplacementVariable,
        const ArrayVarType& rAccelerationVariable,
        const ArrayVarType& rFractionalAccelerationVariable,
        const SizeType DomainSize = 3
        )
    {
        Double3DArray& r_current_velocity = itCurrentNode->FastGetSolutionStepValue(rVelocityVariable);
        Double3DArray& r_current_displacement = itCurrentNode->FastGetSolutionStepValue(rDisplacementVariable);
        Double3DArray& r_current_acceleration = itCurrentNode->FastGetSolutionStepValue(rAccelerationVariable);
        Double3DArray& r_fractional_acceleration = itCurrentNode->FastGetSolutionStepValue(rFractionalAccelerationVariable);

        const Double3DArray& r_previous_displacement = itCurrentNode->FastGetSolutionStepValue(rDisplacementVariable, 1);
        const Double3DArray& r_previous_velocity = itCurrentNode->FastGetSolutionStepValue(rVelocityVariable, 1);
        const Double3DArray& r_previous_acceleration = itCurrentNode->FastGetSolutionStepValue(rAccelerationVariable, 1);


        for (IndexType j = 0; j < DomainSize; j++) {

            r_current_displacement[j] = r_previous_displacement[j] + mTime.Delta * r_previous_velocity[j];
            r_current_displacement[j] += 0.50 * mTime.Delta * mTime.Delta * r_previous_acceleration[j] * ((4.0*mDeltaTime.Fraction-1.0)/(6.0*mDeltaTime.Fraction));
            r_current_displacement[j] -= 0.50 * mTime.Delta * mTime.Delta * r_fractional_acceleration[j] * (1.0/(6.0*mDeltaTime.Fraction*(mDeltaTime.Fraction-1.0)));
            r_current_displacement[j] += 0.50 * mTime.Delta * mTime.Delta * r_current_acceleration[j] * ((2.0*mDeltaTime.Fraction-1.0)/(6.0*(mDeltaTime.Fraction-1.0)));


            r_current_velocity[j] = r_previous_velocity[j];
            r_current_velocity[j] += mTime.Delta * r_previous_acceleration[j] * ((3.0*mDeltaTime.Fraction-1.0)/(6.0*mDeltaTime.Fraction));
            r_current_velocity[j] -= mTime.Delta * r_fractional_acceleration[j] * ((1.0)/(6.0*mDeltaTime.Fraction*(mDeltaTime.Fraction-1.0)));
            r_current_velocity[j] += mTime.Delta * r_current_acceleration[j] * ((3.0*mDeltaTime.Fraction-2.0)/(6.0*(mDeltaTime.Fraction-1.0)));
        } // for DomainSize
    }

    /**
     * @brief This method performs some custom operations to initialize the scheme
     * @param rModelPart The model of the problem to solve
     * @param DomainSize The current dimention of the problem
     */
    virtual void SchemeCustomInitialization(
        ModelPart& rModelPart,
        const SizeType DomainSize = 3
        )
    {
        KRATOS_TRY

        // The array containing the nodes
        NodesArrayType& r_nodes = rModelPart.Nodes();

        // The fisrt node interator
        const auto it_node_begin = rModelPart.NodesBegin();

        // Auxiliar zero array
        const Double3DArray zero_array = ZeroVector(3);

        // If we consider the rotation DoF
        const bool has_dof_for_rot_z = it_node_begin->HasDofFor(ROTATION_Z);

        // Getting dof position
        const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);
        const IndexType rotppos = has_dof_for_rot_z ? it_node_begin->GetDofPosition(ROTATION_X) : 0;

        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            // Current step information "N+1" (before step update).
            auto it_node = it_node_begin + i;

            const double nodal_mass = it_node->GetValue(NODAL_MASS);
            const Double3DArray& r_current_residual = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);
            Double3DArray& r_current_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION);

            // Solution of the explicit equation:
            if (nodal_mass > numerical_limit) {
                r_current_acceleration = r_current_residual / nodal_mass;
            } else {
                r_current_acceleration = zero_array;
            }

            std::array<bool, 3> fix_displacements = {false, false, false};

            fix_displacements[0] = (it_node->GetDof(DISPLACEMENT_X, disppos).IsFixed());
            fix_displacements[1] = (it_node->GetDof(DISPLACEMENT_Y, disppos + 1).IsFixed());
            if (DomainSize == 3)
                fix_displacements[2] = (it_node->GetDof(DISPLACEMENT_Z, disppos + 2).IsFixed());

            for (IndexType j = 0; j < DomainSize; j++) {
                if (fix_displacements[j]) {
                    r_current_acceleration[j] = 0.0;
                }
            } // for DomainSize


             ////// ROTATION DEGRESS OF FREEDOM
            if (has_dof_for_rot_z) {
                const Double3DArray& nodal_inertia = it_node->GetValue(NODAL_INERTIA);
                const Double3DArray& r_current_residual_moment = it_node->FastGetSolutionStepValue(MOMENT_RESIDUAL);
                Double3DArray& r_current_angular_acceleration = it_node->FastGetSolutionStepValue(ANGULAR_ACCELERATION);

                const IndexType initial_k = DomainSize == 3 ? 0 : 2; // We do this because in 2D only the rotation Z is needed, then we start with 2, instead of 0
                for (IndexType kk = initial_k; kk < 3; ++kk) {
                    if (nodal_inertia[kk] > numerical_limit) {
                        r_current_angular_acceleration[kk] = r_current_residual_moment[kk] / nodal_inertia[kk];
                    } else {
                        r_current_angular_acceleration[kk] = 0.0;
                    }
                }

                std::array<bool, 3> fix_rotation = {false, false, false};
                if (DomainSize == 3) {
                    fix_rotation[0] = (it_node->GetDof(ROTATION_X, rotppos).IsFixed());
                    fix_rotation[1] = (it_node->GetDof(ROTATION_Y, rotppos + 1).IsFixed());
                }
                fix_rotation[2] = (it_node->GetDof(ROTATION_Z, rotppos + 2).IsFixed());

                for (IndexType j = initial_k; j < 3; j++) {
                    if (fix_rotation[j]) {
                        r_current_angular_acceleration[j] = 0.0;
                    }
                } // trans DoF
            }   // Rot DoF
        }     // for node parallel
        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param pElement The element to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void Calculate_RHS_Contribution(
        Element::Pointer pCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        this->TCalculate_RHS_Contribution(pCurrentElement, RHS_Contribution, rCurrentProcessInfo);
        KRATOS_CATCH("")
    }

    /**
     * @brief Functions that calculates the RHS of a "condition" object
     * @param pCondition The condition to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer pCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        this->TCalculate_RHS_Contribution(pCurrentCondition, RHS_Contribution, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }



    void CalculateAndAddRHS(ModelPart& rModelPart)
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
            Condition_Calculate_RHS_Contribution((*it_cond.base()), RHS_Contribution, equation_id_vector_dummy, r_current_process_info);
        }

        #pragma omp parallel for firstprivate(RHS_Contribution, equation_id_vector_dummy), schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
            auto it_elem = r_elements.begin() + i;
            Calculate_RHS_Contribution((*it_elem.base()), RHS_Contribution, equation_id_vector_dummy, r_current_process_info);
        }

        KRATOS_CATCH("")
    }



    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY;
        // The current process info
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // The array of nodes
        NodesArrayType& r_nodes = rModelPart.Nodes();

        /// Working in 2D/3D (the definition of DOMAIN_SIZE is check in the Check method)
        const SizeType dim = r_current_process_info[DOMAIN_SIZE];

        // Step Update
        // The first step is time =  initial_time ( 0.0) + delta time
        mTime.Delta = r_current_process_info[DELTA_TIME];
        mTime.MidStep = mTime.Delta * mDeltaTime.Fraction;

        // The iterator of the first node
        const auto it_node_begin = rModelPart.NodesBegin();
        const bool has_dof_for_rot_z = it_node_begin->HasDofFor(ROTATION_Z);


        // Getting dof position
        const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);
        const IndexType rotppos = has_dof_for_rot_z ? it_node_begin->GetDofPosition(ROTATION_X) : 0;

        // _____________________________________________________________________
        // ______________________________ STAGE 1 ______________________________
        // _____________________________________________________________________

        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            // Current step information "N+1" (before step update).
            this->UpdateDegreesOfFreedomStage1(it_node_begin + i,
                VELOCITY,DISPLACEMENT,ACCELERATION,dim);
        } // for Node parallel


        if (has_dof_for_rot_z){
            #pragma omp parallel for schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                // Current step information "N+1" (before step update).
                this->UpdateDegreesOfFreedomStage1(it_node_begin + i,
                    ANGULAR_VELOCITY,ROTATION,ANGULAR_ACCELERATION,dim);
            } // for Node parallel
        }


        InitializeResidual(rModelPart);
        CalculateAndAddRHS(rModelPart);

        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            // Current step information "N+1" (before step update).
            this->UpdateAccelerationStage(
                it_node_begin + i, disppos,FRACTIONAL_ACCELERATION,
                VELOCITY,NODAL_DISPLACEMENT_DAMPING,FORCE_RESIDUAL,NODAL_MASS,
                DISPLACEMENT_X,DISPLACEMENT_Y,DISPLACEMENT_Z,dim);
        } // for Node parallel

        if (has_dof_for_rot_z){
            #pragma omp parallel for schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                // Current step information "N+1" (before step update).
                this->UpdateAccelerationStage(
                    it_node_begin + i, rotppos,FRACTIONAL_ANGULAR_ACCELERATION,
                    ANGULAR_VELOCITY,NODAL_ROTATION_DAMPING,MOMENT_RESIDUAL,
                    NODAL_INERTIA,ROTATION_X,ROTATION_Y,ROTATION_Z,dim);
            } // for Node parallel
        }

        // _____________________________________________________________________
        // ______________________________ STAGE 2 ______________________________
        // _____________________________________________________________________

        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            // Current step information "N+1" (before step update).
            this->UpdateDegreesOfFreedomStage2(it_node_begin + i,
                VELOCITY,DISPLACEMENT,ACCELERATION,FRACTIONAL_ACCELERATION,dim);
        } // for Node parallel

        if (has_dof_for_rot_z){
            #pragma omp parallel for schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                // Current step information "N+1" (before step update).
                this->UpdateDegreesOfFreedomStage2(it_node_begin + i,
                    ANGULAR_VELOCITY,ROTATION,ANGULAR_ACCELERATION,
                    FRACTIONAL_ANGULAR_ACCELERATION,dim);
            } // for Node parallel
        }

        InitializeResidual(rModelPart);
        CalculateAndAddRHS(rModelPart);

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Operations
    ///@{

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

    ///@}
    ///@name Protected Structs
    ///@{

    /**
     * @brief This struct contains the information related with the increment od time step
     */
    struct DeltaTimeParameters {
        double Fraction;        // Fraction of the delta time
    };

    /**
     * @brief This struct contains the details of the time variables
     */
    struct TimeVariables {
        double MidStep;         // Time step*Fraction
        double Delta;          // Time step
    };

    ///@name Protected static Member Variables
    ///@{

    TimeVariables mTime;            /// This struct contains the details of the time variables
    DeltaTimeParameters mDeltaTime; /// This struct contains the information related with the increment od time step

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
    * @brief Functions that calculates the RHS of a "TObjectType" object
    * @param pCurrentEntity The TObjectType to compute
    * @param RHS_Contribution The RHS vector contribution
    * @param rCurrentProcessInfo The current process info instance
    */
    template <typename TObjectType>
    void TCalculate_RHS_Contribution(
        TObjectType pCurrentEntity,
        LocalSystemVectorType& RHS_Contribution,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        pCurrentEntity->CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);

        pCurrentEntity->AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);
        pCurrentEntity->AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, MOMENT_RESIDUAL, rCurrentProcessInfo);
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

    ///@}

}; /* Class ExplicitMultiStageKimScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_MULTI_STAGE_KIM_SCHEME_HPP_INCLUDED  defined */
