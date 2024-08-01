// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Klaus B Sautter (based on the work of MSantasusana)
//

#pragma once

// System includes

// External includes

// Project includes
#include "solving_strategies/schemes/scheme.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/explicit_integration_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/parallel_utilities.h"

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
 * @class ExplicitCentralDifferencesScheme
 * @ingroup StructuralMechanicsApplciation
 * @brief An explicit central difference scheme
 * @details This scheme is quite known, and can be consulted from different sources
 * For example from Wikipedia: https://en.wikipedia.org/wiki/Central_differencing_scheme
 * @author Klaus B Sautter
 */
template <class TSparseSpace,
          class TDenseSpace //= DenseSpace<double>
          >
class ExplicitCentralDifferencesScheme
    : public Scheme<TSparseSpace, TDenseSpace> {

public:
    ///@name Type Definitions
    ///@{

    /// The definition of the base type
    typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

    /// Definition of the current scheme
    typedef ExplicitCentralDifferencesScheme<TSparseSpace, TDenseSpace> ClassType;

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

    /// Counted pointer of ExplicitCentralDifferencesScheme
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitCentralDifferencesScheme);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details The ExplicitCentralDifferencesScheme method
     * @param MaximumDeltaTime The maximum delta time to be considered
     * @param DeltaTimeFraction The delta ttime fraction
     * @param DeltaTimePredictionLevel The prediction level
     */
    ExplicitCentralDifferencesScheme(
        const double MaximumDeltaTime,
        const double DeltaTimeFraction,
        const double DeltaTimePredictionLevel
        )
        : BaseType()
    {
        mDeltaTime.PredictionLevel = DeltaTimePredictionLevel;
        mDeltaTime.Maximum = MaximumDeltaTime;
        mDeltaTime.Fraction = DeltaTimeFraction;
    }

    /**
     * @brief Constructor with parameters
     * @details The ExplicitCentralDifferencesScheme method
     * @param ThisParameters The parameters containing the configuration parameters
     * @warning time_step_prediction_level should be an integer
     */
    ExplicitCentralDifferencesScheme(Parameters ThisParameters = Parameters(R"({})"))
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters); 
    }

    /** Destructor.
    */
    virtual ~ExplicitCentralDifferencesScheme() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided.
     * @details Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model of the problem to solve
     * @return Zero means  all ok
     */
    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY;

        BaseType::Check(rModelPart);

        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2) << "Insufficient buffer size for Central Difference Scheme. It has to be > 2" << std::endl;

        KRATOS_ERROR_IF_NOT(rModelPart.GetProcessInfo().Has(DOMAIN_SIZE)) << "DOMAIN_SIZE not defined on ProcessInfo. Please define" << std::endl;

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

        if ((mDeltaTime.PredictionLevel > 0) && (!BaseType::SchemeIsInitialized())) {
            Parameters prediction_parameters = Parameters(R"(
            {
                "time_step_prediction_level" : 2.0,
                "max_delta_time"             : 1.0e0,
                "safety_factor"              : 0.8
            })" );
            prediction_parameters["time_step_prediction_level"].SetDouble(mDeltaTime.PredictionLevel);
            prediction_parameters["max_delta_time"].SetDouble(mDeltaTime.Maximum);
            ExplicitIntegrationUtilities::CalculateDeltaTime(rModelPart, prediction_parameters);
        }

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Preparing the time values for the first step (where time = initial_time +
        // dt)
        mTime.Current = r_current_process_info[TIME] + r_current_process_info[DELTA_TIME];
        mTime.Delta = r_current_process_info[DELTA_TIME];
        mTime.Middle = mTime.Current - 0.5 * mTime.Delta;
        mTime.Previous = mTime.Current - mTime.Delta;
        mTime.PreviousMiddle = mTime.Current - 1.5 * mTime.Delta;

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
        if (mDeltaTime.PredictionLevel > 1) {
            Parameters prediction_parameters = Parameters(R"(
            {
                "time_step_prediction_level" : 2.0,
                "max_delta_time"             : 1.0e0,
                "safety_factor"              : 0.8
            })" );
            prediction_parameters["time_step_prediction_level"].SetDouble(mDeltaTime.PredictionLevel); // WARNING This could be problematic if PredictionLevel is a double and not a integer
            prediction_parameters["max_delta_time"].SetDouble(mDeltaTime.Maximum);
            ExplicitIntegrationUtilities::CalculateDeltaTime(rModelPart, prediction_parameters);
        }
        InitializeResidual(rModelPart);
        KRATOS_CATCH("")
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

        // Auxiliary values
        const array_1d<double, 3> zero_array = ZeroVector(3);
        // Initializing the variables
        VariableUtils().SetVariable(FORCE_RESIDUAL, zero_array,r_nodes);
        const bool has_dof_for_rot_z = !r_nodes.empty() && r_nodes.begin()->HasDofFor(ROTATION_Z);
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

        if (!r_nodes.empty()) {
            // The first iterator of the array of nodes
            const auto it_node_begin = rModelPart.NodesBegin();

            // Construct the node initializer lambda based on whether rot_z exists
            std::function<void(Node&)> initializer_base, initializer;
            const array_1d<double, 3> zero_array = ZeroVector(3);

            initializer_base = [&zero_array, DomainSize](Node& rNode){
                rNode.SetValue(NODAL_MASS, 0.0);
                rNode.FastGetSolutionStepValue(MIDDLE_VELOCITY) = zero_array;

                array_1d<double, 3>& r_middle_velocity = rNode.FastGetSolutionStepValue(MIDDLE_VELOCITY);
                const array_1d<double, 3>& r_current_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
                array_1d<double, 3>& r_current_residual = rNode.FastGetSolutionStepValue(FORCE_RESIDUAL);
                //array_1d<double,3>& r_current_displacement  = rNode.FastGetSolutionStepValue(DISPLACEMENT);

                for (IndexType j = 0; j < DomainSize; j++) {
                    r_middle_velocity[j] = r_current_velocity[j];
                    r_current_residual[j] = 0.0;
                    //r_current_displacement[j] = 0.0; // this might be wrong for presribed displacement // NOTE: then you should check if the dof is fixed
                }
            };

            const bool has_dof_for_rot_z = it_node_begin->HasDofFor(ROTATION_Z);
            if (has_dof_for_rot_z) {
                initializer = [&zero_array, DomainSize, &initializer_base](Node& rNode){
                    initializer_base(rNode);
                    rNode.SetValue(NODAL_INERTIA, zero_array);
                    rNode.FastGetSolutionStepValue(MIDDLE_ANGULAR_VELOCITY) = zero_array;

                    array_1d<double, 3>& r_middle_angular_velocity = rNode.FastGetSolutionStepValue(MIDDLE_ANGULAR_VELOCITY);
                    const array_1d<double, 3>& r_current_angular_velocity = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY);
                    array_1d<double, 3>& r_current_residual_moment = rNode.FastGetSolutionStepValue(MOMENT_RESIDUAL);
                    //array_1d<double,3>& current_rotation = rNode.FastGetSolutionStepValue(ROTATION);

                    const IndexType initial_j = DomainSize == 3 ? 0 : 2; // We do this because in 2D only the rotation Z is needed, then we start with 2, instead of 0
                    for (IndexType j = initial_j; j < 3; j++) {
                        r_middle_angular_velocity[j] = r_current_angular_velocity[j];
                        r_current_residual_moment[j] = 0.0;
                        // current_rotation[j] = 0.0; // this might be wrong for presribed rotations // NOTE: then you should check if the dof is fixed
                    }
                };
            } else {
                initializer = initializer_base;
            }

            /// Initialize all nodes
            block_for_each(r_nodes, initializer);
        } // if not nodes.empty()

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
        // The array of nodes
        NodesArrayType& r_nodes = rModelPart.Nodes();

        if (!r_nodes.empty()) {
            // The current process info
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            /// Working in 2D/3D (the definition of DOMAIN_SIZE is check in the Check method)
            const SizeType dim = r_current_process_info[DOMAIN_SIZE];

            // Step Update
            // The first step is time =  initial_time ( 0.0) + delta time
            mTime.Current = r_current_process_info[TIME];
            mTime.Delta = r_current_process_info[DELTA_TIME];

            mTime.Middle   = mTime.Current - 0.50*mTime.Delta;
            mTime.Previous = mTime.Current - 1.00*mTime.Delta;
            mTime.PreviousMiddle = mTime.Middle - 1.00*mTime.Delta;

            if (mTime.Previous<0.0) mTime.Previous=0.00;
            if (mTime.PreviousMiddle<0.0) mTime.PreviousMiddle=0.00;
            // The iterator of the first node
            const auto it_node_begin = rModelPart.NodesBegin();
            const bool has_dof_for_rot_z = it_node_begin->HasDofFor(ROTATION_Z);

            // Getting dof position
            const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);
            const IndexType rotppos = has_dof_for_rot_z ? it_node_begin->GetDofPosition(ROTATION_X) : 0;

            // Construct loop lambda
            std::function<void(SizeType)> updater_base, updater;
            updater_base = [it_node_begin, dim, &disppos, this](SizeType index){
                this->UpdateTranslationalDegreesOfFreedom(it_node_begin + index, disppos, dim);
            };

            if (has_dof_for_rot_z) {
                updater = [it_node_begin, dim, &rotppos, this, &updater_base](SizeType index){
                    updater_base(index);
                    this->UpdateRotationalDegreesOfFreedom(it_node_begin + index, rotppos, dim);
                };
            } else {
                updater = updater_base;
            }

            // Update nodes
            IndexPartition<SizeType>(r_nodes.size()).for_each(updater);
        } // if not nodes.empty()

        KRATOS_CATCH("")
    }

    /**
     * @brief This method updates the translation DoF
     * @param itCurrentNode The iterator of the current node
     * @param DisplacementPosition The position of the displacement dof on the database
     * @param DomainSize The current dimention of the problem
     */
    void UpdateTranslationalDegreesOfFreedom(
        NodeIterator itCurrentNode,
        const IndexType DisplacementPosition,
        const SizeType DomainSize = 3
        )
    {

        const double nodal_mass = itCurrentNode->GetValue(NODAL_MASS);
        const double nodal_displacement_damping = itCurrentNode->GetValue(NODAL_DISPLACEMENT_DAMPING);
        const array_1d<double, 3>& r_current_residual = itCurrentNode->FastGetSolutionStepValue(FORCE_RESIDUAL);


        array_1d<double, 3>& r_current_velocity = itCurrentNode->FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3>& r_current_displacement = itCurrentNode->FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3>& r_middle_velocity = itCurrentNode->FastGetSolutionStepValue(MIDDLE_VELOCITY);

        array_1d<double, 3>& r_current_acceleration = itCurrentNode->FastGetSolutionStepValue(ACCELERATION);

        const array_1d<double, 3>& r_previous_displacement = itCurrentNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
        const array_1d<double, 3>& r_previous_middle_velocity = itCurrentNode->FastGetSolutionStepValue(MIDDLE_VELOCITY, 1);
        // Solution of the explicit equation:
        if (nodal_mass > numerical_limit)
            noalias(r_current_acceleration) = (r_current_residual - nodal_displacement_damping * r_current_velocity) / nodal_mass;
        else
            noalias(r_current_acceleration) = ZeroVector(3);

        std::array<bool, 3> fix_displacements = {false, false, false};

        fix_displacements[0] = (itCurrentNode->GetDof(DISPLACEMENT_X, DisplacementPosition).IsFixed());
        fix_displacements[1] = (itCurrentNode->GetDof(DISPLACEMENT_Y, DisplacementPosition + 1).IsFixed());
        if (DomainSize == 3)
            fix_displacements[2] = (itCurrentNode->GetDof(DISPLACEMENT_Z, DisplacementPosition + 2).IsFixed());

        for (IndexType j = 0; j < DomainSize; j++) {
            if (fix_displacements[j]) {
                r_current_acceleration[j] = 0.0;
                r_middle_velocity[j] = 0.0;
            }

            r_current_velocity[j] =  r_previous_middle_velocity[j] + (mTime.Previous - mTime.PreviousMiddle) * r_current_acceleration[j]; //+ actual_velocity;
            r_middle_velocity[j] = r_current_velocity[j] + (mTime.Middle - mTime.Previous) * r_current_acceleration[j];
            r_current_displacement[j] = r_previous_displacement[j] + mTime.Delta * r_middle_velocity[j];
        } // for DomainSize
    }

    /**
     * @brief This method updates the rotation DoF
     * @param itCurrentNode The iterator of the current node
     * @param RotationPosition The position of the rotation dof on the database
     * @param DomainSize The current dimention of the problem
     */
    void UpdateRotationalDegreesOfFreedom(
        NodeIterator itCurrentNode,
        const IndexType RotationPosition,
        const SizeType DomainSize = 3
        )
    {
        ////// ROTATION DEGRESS OF FREEDOM
        const array_1d<double, 3>& nodal_inertia = itCurrentNode->GetValue(NODAL_INERTIA);
        const array_1d<double, 3>& nodal_rotational_damping = itCurrentNode->GetValue(NODAL_ROTATION_DAMPING);
        const array_1d<double, 3>& r_current_residual_moment = itCurrentNode->FastGetSolutionStepValue(MOMENT_RESIDUAL);
        array_1d<double, 3>& r_current_angular_velocity = itCurrentNode->FastGetSolutionStepValue(ANGULAR_VELOCITY);
        array_1d<double, 3>& r_current_rotation = itCurrentNode->FastGetSolutionStepValue(ROTATION);
        array_1d<double, 3>& r_middle_angular_velocity = itCurrentNode->FastGetSolutionStepValue(MIDDLE_ANGULAR_VELOCITY);
        array_1d<double, 3>& r_current_angular_acceleration = itCurrentNode->FastGetSolutionStepValue(ANGULAR_ACCELERATION);


        const array_1d<double, 3>& r_previous_rotation = itCurrentNode->FastGetSolutionStepValue(ROTATION, 1);
        const array_1d<double, 3>& r_previous_middle_angular_velocity = itCurrentNode->FastGetSolutionStepValue(MIDDLE_ANGULAR_VELOCITY, 1);

        const IndexType initial_k = DomainSize == 3 ? 0 : 2; // We do this because in 2D only the rotation Z is needed, then we start with 2, instead of 0
        for (IndexType kk = initial_k; kk < 3; ++kk) {
            if (nodal_inertia[kk] > numerical_limit)
                r_current_angular_acceleration[kk] = (r_current_residual_moment[kk] - nodal_rotational_damping[kk] * r_current_angular_velocity[kk]) / nodal_inertia[kk];
            else
                r_current_angular_acceleration[kk] = 0.0;
        }

        std::array<bool, 3> fix_rotation = {false, false, false};
        if (DomainSize == 3) {
            fix_rotation[0] = (itCurrentNode->GetDof(ROTATION_X, RotationPosition).IsFixed());
            fix_rotation[1] = (itCurrentNode->GetDof(ROTATION_Y, RotationPosition + 1).IsFixed());
        }
        fix_rotation[2] = (itCurrentNode->GetDof(ROTATION_Z, RotationPosition + 2).IsFixed());

        for (IndexType j = initial_k; j < 3; j++) {
            if (fix_rotation[j]) {
                r_current_angular_acceleration[j] = 0.0;
                r_middle_angular_velocity[j] = 0.0;
            }
            r_current_angular_velocity[j] = r_previous_middle_angular_velocity[j] + (mTime.Previous - mTime.PreviousMiddle) * r_current_angular_acceleration[j];
            r_middle_angular_velocity[j] = r_current_angular_velocity[j] + (mTime.Middle - mTime.Previous) * r_current_angular_acceleration[j];
            r_current_rotation[j] = r_previous_rotation[j] + mTime.Delta * r_middle_angular_velocity[j];
        }
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

        if (!r_nodes.empty()) {
            // The fisrt node interator
            const auto it_node_begin = rModelPart.NodesBegin();

            // If we consider the rotation DoF
            const bool has_dof_for_rot_z = it_node_begin->HasDofFor(ROTATION_Z);

            // Auxiliary zero array
            const array_1d<double, 3> zero_array = ZeroVector(3);

            // Getting dof position
            const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);
            const IndexType rotppos = has_dof_for_rot_z ? it_node_begin->GetDofPosition(ROTATION_X) : 0;

            // Construct initializer lambda
            std::function<void(Node&)> initializer_base, initializer;

            initializer_base = [&zero_array, &disppos, DomainSize, this](Node& rNode) {
                const double nodal_mass = rNode.GetValue(NODAL_MASS);
                const array_1d<double, 3>& r_current_residual = rNode.FastGetSolutionStepValue(FORCE_RESIDUAL);

                array_1d<double, 3>& r_current_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
                //array_1d<double,3>& r_current_displacement = rNode.FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double, 3>& r_middle_velocity = rNode.FastGetSolutionStepValue(MIDDLE_VELOCITY);

                array_1d<double, 3>& r_current_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION);

                // Solution of the explicit equation:
                if (nodal_mass > numerical_limit) {
                    r_current_acceleration = r_current_residual / nodal_mass;
                } else {
                    r_current_acceleration = zero_array;
                }

                std::array<bool, 3> fix_displacements = {false, false, false};

                fix_displacements[0] = (rNode.GetDof(DISPLACEMENT_X, disppos).IsFixed());
                fix_displacements[1] = (rNode.GetDof(DISPLACEMENT_Y, disppos + 1).IsFixed());
                if (DomainSize == 3)
                    fix_displacements[2] = (rNode.GetDof(DISPLACEMENT_Z, disppos + 2).IsFixed());

                for (IndexType j = 0; j < DomainSize; j++) {
                    if (fix_displacements[j]) {
                        r_current_acceleration[j] = 0.0;
                        r_middle_velocity[j] = 0.0;
                    }

                    r_middle_velocity[j] = 0.0 + (mTime.Middle - mTime.Previous) * r_current_acceleration[j];
                    r_current_velocity[j] = r_middle_velocity[j] + (mTime.Previous - mTime.PreviousMiddle) * r_current_acceleration[j]; //+ actual_velocity;
                    // r_current_displacement[j]  = 0.0;

                } // for DomainSize
            };

            if (has_dof_for_rot_z) {
                initializer = [&rotppos, &initializer_base, DomainSize, this](Node& rNode) {
                    initializer_base(rNode);
                    const array_1d<double, 3>& nodal_inertia = rNode.GetValue(NODAL_INERTIA);
                    const array_1d<double, 3>& r_current_residual_moment = rNode.FastGetSolutionStepValue(MOMENT_RESIDUAL);
                    array_1d<double, 3>& r_current_angular_velocity = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY);
                    // array_1d<double,3>& current_rotation = rNode.FastGetSolutionStepValue(ROTATION);
                    array_1d<double, 3>& r_middle_angular_velocity = rNode.FastGetSolutionStepValue(MIDDLE_ANGULAR_VELOCITY);
                    array_1d<double, 3>& r_current_angular_acceleration = rNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION);

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
                        fix_rotation[0] = (rNode.GetDof(ROTATION_X, rotppos).IsFixed());
                        fix_rotation[1] = (rNode.GetDof(ROTATION_Y, rotppos + 1).IsFixed());
                    }
                    fix_rotation[2] = (rNode.GetDof(ROTATION_Z, rotppos + 2).IsFixed());

                    for (IndexType j = initial_k; j < 3; j++) {
                        if (fix_rotation[j]) {
                            r_current_angular_acceleration[j] = 0.0;
                            r_middle_angular_velocity[j] = 0.0;
                        }

                        r_middle_angular_velocity[j] = 0.0 +  (mTime.Middle - mTime.Previous) * r_current_angular_acceleration[j];
                        r_current_angular_velocity[j] = r_middle_angular_velocity[j] +  (mTime.Previous - mTime.PreviousMiddle) *  r_current_angular_acceleration[j];
                        //current_rotation[j] = 0.0;
                    }
                };
            } else {
                initializer = initializer_base;
            }

            // Initialize all nodes
            block_for_each(r_nodes, initializer);
        } // if not nodes.empty()

        mTime.Previous = mTime.Current;
        mTime.PreviousMiddle = mTime.Middle;
        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rCurrentElement The element to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        this->TCalculateRHSContribution(rCurrentElement, RHS_Contribution, rCurrentProcessInfo);
        KRATOS_CATCH("")
    }

    /**
     * @brief Functions that calculates the RHS of a "condition" object
     * @param rCurrentCondition The condition to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        this->TCalculateRHSContribution(rCurrentCondition, RHS_Contribution, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }


    void CalculateAndAddRHS(ModelPart& rModelPart)
    {
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        auto& r_conditions = rModelPart.Conditions();
        auto& r_elements = rModelPart.Elements();

        using TLS = std::tuple<LocalSystemVectorType,Element::EquationIdVectorType>;
        TLS thread_local_storage; // {RHS_contribution, equation_id_vector_dummy}

        // Compute on conditions
        block_for_each(r_conditions, thread_local_storage, [&r_current_process_info, this](Condition& r_condition, TLS& r_rhs_contrib_and_equation_ids){
            CalculateRHSContribution(r_condition,
                                     std::get<0>(r_rhs_contrib_and_equation_ids),
                                     std::get<1>(r_rhs_contrib_and_equation_ids),
                                     r_current_process_info);
        });

        // Compute on elements
        block_for_each(r_elements, thread_local_storage, [&r_current_process_info, this](Element& r_element, TLS& r_rhs_contrib_and_equation_ids){
            CalculateRHSContribution(r_element,
                                     std::get<0>(r_rhs_contrib_and_equation_ids),
                                     std::get<1>(r_rhs_contrib_and_equation_ids),
                                     r_current_process_info);
        });

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
        CalculateAndAddRHS(rModelPart);
        KRATOS_CATCH("")
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                       : "central_differences",
            "time_step_prediction_level" : 0.0,
            "fraction_delta_time"        : 0.9,
            "max_delta_time"             : 1.0e0
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "central_differences";
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

    ///@}
    ///@name Protected Structs
    ///@{

    /**
     * @brief This struct contains the information related with the increment od time step
     */
    struct DeltaTimeParameters {
        double PredictionLevel; // 0, 1, 2 // NOTE: Should be a integer?
        double Maximum;         // Maximum delta time
        double Fraction;        // Fraction of the delta time
    };

    /**
     * @brief This struct contains the details of the time variables
     */
    struct TimeVariables {
        double PreviousMiddle; // n-1/2
        double Previous;       // n
        double Middle;         // n+1/2
        double Current;        // n+1

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
    
    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);
        mDeltaTime.PredictionLevel = ThisParameters["time_step_prediction_level"].GetDouble();
        mDeltaTime.Maximum = ThisParameters["max_delta_time"].GetDouble();
        mDeltaTime.Fraction = ThisParameters["fraction_delta_time"].GetDouble();
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
    * @brief Functions that calculates the RHS of a "TObjectType" object
    * @param rCurrentEntity The TObjectType to compute
    * @param RHS_Contribution The RHS vector contribution
    * @param rCurrentProcessInfo The current process info instance
    */
    template <typename TObjectType>
    void TCalculateRHSContribution(
        TObjectType& rCurrentEntity,
        LocalSystemVectorType& RHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        rCurrentEntity.CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);

        rCurrentEntity.AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);
        rCurrentEntity.AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, MOMENT_RESIDUAL, rCurrentProcessInfo);
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

}; /* Class ExplicitCentralDifferencesScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/
