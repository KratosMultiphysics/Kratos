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
//  Main authors:    Klaus B Sautter (based on the work of MSantasusana)
//
//

#if !defined(KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME_HPP_INCLUDED)
#define KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME_HPP_INCLUDED

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
        : Scheme<TSparseSpace, TDenseSpace>()
    {
        mDeltaTime.PredictionLevel = DeltaTimePredictionLevel;
        mDeltaTime.Maximum = MaximumDeltaTime;
        mDeltaTime.Fraction = DeltaTimeFraction;
    }

    /**
     * @brief Constructor with parameters
     * @details The ExplicitCentralDifferencesScheme method
     * @param rParameters The parameters containing the configuration parameters
     * @warning time_step_prediction_level should be an integer
     */
    ExplicitCentralDifferencesScheme(Parameters rParameters =  Parameters(R"({})"))
        : Scheme<TSparseSpace, TDenseSpace>()
    {
        Parameters default_parameters = Parameters(R"(
        {
            "time_step_prediction_level" : 0.0,
            "fraction_delta_time"        : 0.9,
            "max_delta_time"             : 1.0e0
        })" );

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mDeltaTime.PredictionLevel = rParameters["time_step_prediction_level"].GetDouble();
        mDeltaTime.Maximum = rParameters["max_delta_time"].GetDouble();
        mDeltaTime.Fraction = rParameters["fraction_delta_time"].GetDouble();
    }

    /** Destructor.
    */
    virtual ~ExplicitCentralDifferencesScheme() {}

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

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

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
        const array_1d<double, 3> zero_array = ZeroVector(3);
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
        const array_1d<double, 3> zero_array = ZeroVector(3);
        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            auto it_node = (it_node_begin + i);
            it_node->SetValue(NODAL_MASS, 0.0);
            array_1d<double, 3>& r_middle_velocity = it_node->FastGetSolutionStepValue(MIDDLE_VELOCITY);
            r_middle_velocity  = ZeroVector(3);
        }
        const bool has_dof_for_rot_z = it_node_begin->HasDofFor(ROTATION_Z);
        if (has_dof_for_rot_z) {
            #pragma omp parallel for schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                auto it_node = (it_node_begin + i);
                it_node->SetValue(NODAL_INERTIA, zero_array);
                array_1d<double, 3>& r_middle_angular_velocity = it_node->FastGetSolutionStepValue(MIDDLE_ANGULAR_VELOCITY);
                r_middle_angular_velocity  = ZeroVector(3);
            }
        }

        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            auto it_node = (it_node_begin + i);

            array_1d<double, 3>& r_middle_velocity = it_node->FastGetSolutionStepValue(MIDDLE_VELOCITY);
            const array_1d<double, 3>& r_current_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3>& r_current_residual = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);
//             array_1d<double,3>& r_current_displacement  = it_node->FastGetSolutionStepValue(DISPLACEMENT);

            for (IndexType j = 0; j < DomainSize; j++) {
                r_middle_velocity[j] = r_current_velocity[j];
                r_current_residual[j] = 0.0;
//                 r_current_displacement[j] = 0.0; // this might be wrong for presribed displacement // NOTE: then you should check if the dof is fixed
            }

            if (has_dof_for_rot_z) {
                array_1d<double, 3>& r_middle_angular_velocity = it_node->FastGetSolutionStepValue(MIDDLE_ANGULAR_VELOCITY);
                const array_1d<double, 3>& r_current_angular_velocity = it_node->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                array_1d<double, 3>& r_current_residual_moment = it_node->FastGetSolutionStepValue(MOMENT_RESIDUAL);
//                 array_1d<double,3>& current_rotation = it_node->FastGetSolutionStepValue(ROTATION);

                const IndexType initial_j = DomainSize == 3 ? 0 : 2; // We do this because in 2D only the rotation Z is needed, then we start with 2, instead of 0
                for (IndexType j = initial_j; j < 3; j++) {
                    r_middle_angular_velocity[j] = r_current_angular_velocity[j];
                    r_current_residual_moment[j] = 0.0;
                    // current_rotation[j] = 0.0; // this might be wrong for presribed rotations // NOTE: then you should check if the dof is fixed
                }
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

        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            // Current step information "N+1" (before step update).
            this->UpdateTranslationalDegreesOfFreedom(it_node_begin + i, disppos, dim);
        } // for Node parallel

        if (has_dof_for_rot_z){
            #pragma omp parallel for schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                this->UpdateRotationalDegreesOfFreedom(it_node_begin + i, rotppos, dim);
            } // for Node parallel
        }


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

        // The fisrt node interator
        const auto it_node_begin = rModelPart.NodesBegin();

        // If we consider the rotation DoF
        const bool has_dof_for_rot_z = it_node_begin->HasDofFor(ROTATION_Z);

        // Auxiliar zero array
        const array_1d<double, 3> zero_array = ZeroVector(3);

        // Getting dof position
        const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);
        const IndexType rotppos = has_dof_for_rot_z ? it_node_begin->GetDofPosition(ROTATION_X) : 0;

        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            // Current step information "N+1" (before step update).
            auto it_node = it_node_begin + i;

            const double nodal_mass = it_node->GetValue(NODAL_MASS);
            const array_1d<double, 3>& r_current_residual = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);

            array_1d<double, 3>& r_current_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
//             array_1d<double,3>& r_current_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3>& r_middle_velocity = it_node->FastGetSolutionStepValue(MIDDLE_VELOCITY);

            array_1d<double, 3>& r_current_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION);

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
                    r_middle_velocity[j] = 0.0;
                }

                r_middle_velocity[j] = 0.0 + (mTime.Middle - mTime.Previous) * r_current_acceleration[j];
                r_current_velocity[j] = r_middle_velocity[j] + (mTime.Previous - mTime.PreviousMiddle) * r_current_acceleration[j]; //+ actual_velocity;
                // r_current_displacement[j]  = 0.0;

            } // for DomainSize

            ////// ROTATION DEGRESS OF FREEDOM
            if (has_dof_for_rot_z) {
                const array_1d<double, 3>& nodal_inertia = it_node->GetValue(NODAL_INERTIA);
                const array_1d<double, 3>& r_current_residual_moment = it_node->FastGetSolutionStepValue(MOMENT_RESIDUAL);
                array_1d<double, 3>& r_current_angular_velocity = it_node->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                // array_1d<double,3>& current_rotation = it_node->FastGetSolutionStepValue(ROTATION);
                array_1d<double, 3>& r_middle_angular_velocity = it_node->FastGetSolutionStepValue(MIDDLE_ANGULAR_VELOCITY);
                array_1d<double, 3>& r_current_angular_acceleration = it_node->FastGetSolutionStepValue(ANGULAR_ACCELERATION);

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
                        r_middle_angular_velocity[j] = 0.0;
                    }

                    r_middle_angular_velocity[j] = 0.0 +  (mTime.Middle - mTime.Previous) * r_current_angular_acceleration[j];
                    r_current_angular_velocity[j] = r_middle_angular_velocity[j] +  (mTime.Previous - mTime.PreviousMiddle) *  r_current_angular_acceleration[j];
//                     current_rotation[j] = 0.0;
                } // trans DoF
            }   // Rot DoF
        }     // for node parallel

        mTime.Previous = mTime.Current;
        mTime.PreviousMiddle = mTime.Middle;
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

}; /* Class ExplicitCentralDifferencesScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_CENTRAL_DIFFERENCES_SCHEME  defined */
