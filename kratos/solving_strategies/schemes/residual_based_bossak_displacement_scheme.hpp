//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Main authors:  Josep Maria Carbonell
//                 Vicente Mataix Ferrandiz
//                 Andreas Winterstein (refactoring)
//

#if !defined(KRATOS_RESIDUAL_BASED_BOSSAK_DISPLACEMENT_SCHEME )
#define  KRATOS_RESIDUAL_BASED_BOSSAK_DISPLACEMENT_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/residual_based_implicit_time_scheme.h"
#include "includes/variables.h"
#include "includes/checks.h"

namespace Kratos
{
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
 * @class ResidualBasedBossakDisplacementScheme
 * @ingroup KratosCore
 * @brief Bossak integration scheme (for linear and nonlinear dynamic problems) for displacements
 * @details This is a dynamic implicit scheme based of the Bossak algorithm for displacements.
 * The parameter Alpha of Bossak introduces damping, the value of Bossak is from 0 to -0.5 (negative)
 * Implementation according to: "An alpha modification of Newmark's method; W.L. Wood, M. Bossak, O.C. Zienkiewicz;
 * Numerical Methods in Engineering; 1980"
 * @author Josep Maria Carbonell
 * @author Vicente Mataix Ferrandiz
 * @author Andreas Winterstein (refactoring)
 */
template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedBossakDisplacementScheme
    : public ResidualBasedImplicitTimeScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedBossakDisplacementScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                                  BaseType;

    typedef ResidualBasedImplicitTimeScheme<TSparseSpace,TDenseSpace> ImplicitBaseType;

    typedef typename ImplicitBaseType::TDataType                             TDataType;

    typedef typename ImplicitBaseType::DofsArrayType                     DofsArrayType;

    typedef typename Element::DofsVectorType                            DofsVectorType;

    typedef typename ImplicitBaseType::TSystemMatrixType             TSystemMatrixType;

    typedef typename ImplicitBaseType::TSystemVectorType             TSystemVectorType;

    typedef typename ImplicitBaseType::LocalSystemVectorType     LocalSystemVectorType;

    typedef typename ImplicitBaseType::LocalSystemMatrixType     LocalSystemMatrixType;

    typedef ModelPart::NodeIterator                                       NodeIterator;

    typedef ModelPart::NodesContainerType                               NodesArrayType;

    typedef ModelPart::ElementsContainerType                         ElementsArrayType;

    typedef ModelPart::ConditionsContainerType                     ConditionsArrayType;

    typedef typename BaseType::Pointer                                 BaseTypePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor. (with parameters)
     * @detail The bossak method
     * @param ThisParameters The parameters containing the configuration
     */
    explicit ResidualBasedBossakDisplacementScheme(Parameters ThisParameters)
        : ResidualBasedBossakDisplacementScheme(ThisParameters.Has("damp_factor_m") ? ThisParameters["damp_factor_m"].GetDouble() : -0.3,
                                                ThisParameters.Has("newmark_beta") ? ThisParameters["newmark_beta"].GetDouble() : 0.25)
    {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "name"          : "ResidualBasedBossakDisplacementScheme",
            "damp_factor_m" : -0.3,
            "newmark_beta"  : 0.25
        })" );
        ThisParameters.ValidateAndAssignDefaults(default_parameters);
    }

    /**
     * @brief Constructor.
     * @detail The bossak method
     * @param Alpha is the Bossak parameter. Default value is 0, which is the Newmark method
     * @param NewarkBeta the Newmark parameter. Default value is 0.25, for mean constant acceleration.
     */
    explicit ResidualBasedBossakDisplacementScheme(const double Alpha = 0.0)
        : ResidualBasedBossakDisplacementScheme(Alpha, 0.25)
    {
    }

    /**
     * @brief Constructor.
     * @detail The bossak method
     * @param Alpha is the Bossak parameter. Default value is 0, which is the Newmark method
     * @param NewarkBeta the Newmark parameter. Default value is 0.25, for mean constant acceleration.
     */
    explicit ResidualBasedBossakDisplacementScheme(const double Alpha, const double NewmarkBeta)
        :ImplicitBaseType()
    {
        // For pure Newmark Scheme
        mBossak.alpha = Alpha;
        mNewmark.beta = NewmarkBeta;
        mNewmark.gamma = 0.5;

        CalculateBossakCoefficients();

        // Allocate auxiliary memory
        const std::size_t num_threads = OpenMPUtils::GetNumThreads();

        mVector.v.resize(num_threads);
        mVector.a.resize(num_threads);
        mVector.ap.resize(num_threads);

        KRATOS_DETAIL("MECHANICAL SCHEME: The Bossak Time Integration Scheme ") << "[alpha_m= " << mBossak.alpha << " beta= " << mNewmark.beta << " gamma= " << mNewmark.gamma << "]" <<std::endl;
    }

    /**
     * @brief Copy Constructor.
     */
    explicit ResidualBasedBossakDisplacementScheme(ResidualBasedBossakDisplacementScheme& rOther)
        :ImplicitBaseType(rOther)
        ,mBossak(rOther.mBossak)
        ,mNewmark(rOther.mNewmark)
        ,mVector(rOther.mVector)
    {
    }

    /**
     * @brief Clone method
     */
    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new ResidualBasedBossakDisplacementScheme(*this) );
    }

    /** Destructor.
     */
    ~ResidualBasedBossakDisplacementScheme
    () override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Recalculates the Newmark coefficients, taking into account the alpha parameters
     * @param beta The Bossak beta coefficient
     * @param gamma The Bossak gamma coefficient
     */
    void CalculateBossakCoefficients()
    {
        mBossak.beta  = (1.0 - mBossak.alpha) * (1.0 - mBossak.alpha) * mNewmark.beta;
        mBossak.gamma = mNewmark.gamma  - mBossak.alpha;
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
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY;

        mpDofUpdater->UpdateDofs(rDofSet, rDx);

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );
        const auto it_node_begin = rModelPart.Nodes().begin();

        #pragma omp parallel for
        for(int i = 0;  i < num_nodes; ++i) {
            auto it_node = it_node_begin + i;

            array_1d<double, 3 > delta_displacement;

            noalias(delta_displacement) = it_node->FastGetSolutionStepValue(DISPLACEMENT) - it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);

            array_1d<double, 3>& r_current_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& r_previous_velocity = it_node->FastGetSolutionStepValue(VELOCITY, 1);

            array_1d<double, 3>& r_current_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION);
            const array_1d<double, 3>& r_previous_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION, 1);

            UpdateVelocity(it_node, r_current_velocity, delta_displacement, r_previous_velocity, r_previous_acceleration);
            UpdateAcceleration(it_node, r_current_acceleration, delta_displacement, r_previous_velocity, r_previous_acceleration);
        }

        KRATOS_CATCH( "" );
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
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY;

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const double delta_time = r_current_process_info[DELTA_TIME];

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );
        const auto it_node_begin = rModelPart.Nodes().begin();

        // Getting position
        const int disppos_x = it_node_begin->HasDofFor(DISPLACEMENT_X) ? it_node_begin->GetDofPosition(DISPLACEMENT_X) : -1;
        const int velpos_x = it_node_begin->HasDofFor(VELOCITY_X) ? it_node_begin->GetDofPosition(VELOCITY_X) : -1;
        const int accelpos_x = it_node_begin->HasDofFor(ACCELERATION_X) ? it_node_begin->GetDofPosition(ACCELERATION_X) : -1;
        const int disppos_y = it_node_begin->HasDofFor(DISPLACEMENT_Y) ? it_node_begin->GetDofPosition(DISPLACEMENT_Y) : -1;
        const int velpos_y = it_node_begin->HasDofFor(VELOCITY_Y) ? it_node_begin->GetDofPosition(VELOCITY_Y) : -1;
        const int accelpos_y = it_node_begin->HasDofFor(ACCELERATION_Y) ? it_node_begin->GetDofPosition(ACCELERATION_Y) : -1;
        const int disppos_z = it_node_begin->HasDofFor(DISPLACEMENT_Z) ? it_node_begin->GetDofPosition(DISPLACEMENT_Z) : -1;
        const int velpos_z = it_node_begin->HasDofFor(VELOCITY_Z) ? it_node_begin->GetDofPosition(VELOCITY_Z) : -1;
        const int accelpos_z = it_node_begin->HasDofFor(ACCELERATION_Z) ? it_node_begin->GetDofPosition(ACCELERATION_Z) : -1;

        // Auxiliar variables
        array_1d<double, 3 > delta_displacement;
        bool predicted_x, predicted_y, predicted_z;

        #pragma omp parallel for private(delta_displacement, predicted_x, predicted_y, predicted_z)
        for(int i = 0;  i < num_nodes; ++i) {
            auto it_node = it_node_begin + i;

            predicted_x = false;
            predicted_y = false;
            predicted_z = false;

            //Predicting: NewDisplacement = r_previous_displacement + r_previous_velocity * delta_time;
            //ATTENTION::: the prediction is performed only on free nodes

            const array_1d<double, 3>& r_previous_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double, 3>& r_previous_velocity     = it_node->FastGetSolutionStepValue(VELOCITY,     1);
            const array_1d<double, 3>& r_previous_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double, 3>& r_current_acceleration        = it_node->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3>& r_current_velocity            = it_node->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3>& r_current_displacement        = it_node->FastGetSolutionStepValue(DISPLACEMENT);

            if (accelpos_x > -1) {
                if (it_node->GetDof(ACCELERATION_X, accelpos_x).IsFixed()) {
                    delta_displacement[0] = (r_current_acceleration[0] + mBossak.c3 * r_previous_acceleration[0] +  mBossak.c2 * r_previous_velocity[0])/mBossak.c0;
                    r_current_displacement[0] =  r_previous_displacement[0] + delta_displacement[0];
                    predicted_x = true;
                }
            }
            if (velpos_x > -1 && !predicted_x) {
                if (it_node->GetDof(VELOCITY_X, velpos_x).IsFixed() && !predicted_x) {
                    delta_displacement[0] = (r_current_velocity[0] + mBossak.c4 * r_previous_velocity[0] + mBossak.c5 * r_previous_acceleration[0])/mBossak.c1;
                    r_current_displacement[0] =  r_previous_displacement[0] + delta_displacement[0];
                    predicted_x = true;
                }
            }
            if (disppos_x > -1 && !predicted_x) {
                if (!it_node->GetDof(DISPLACEMENT_X, disppos_x).IsFixed() && !predicted_x) {
                    r_current_displacement[0] = r_previous_displacement[0] + delta_time * r_previous_velocity[0] + 0.5 * std::pow(delta_time, 2) * r_previous_acceleration[0];
                }
            }

            if (accelpos_y > -1) {
                if (it_node->GetDof(ACCELERATION_Y, accelpos_y).IsFixed()) {
                    delta_displacement[1] = (r_current_acceleration[1] + mBossak.c3 * r_previous_acceleration[1] +  mBossak.c2 * r_previous_velocity[1])/mBossak.c0;
                    r_current_displacement[1] =  r_previous_displacement[1] + delta_displacement[1];
                    predicted_y = true;
                }
            }
            if (velpos_y > -1 && !predicted_y) {
                if (it_node->GetDof(VELOCITY_Y, velpos_y).IsFixed() && !predicted_y) {
                    delta_displacement[1] = (r_current_velocity[1] + mBossak.c4 * r_previous_velocity[1] + mBossak.c5 * r_previous_acceleration[1])/mBossak.c1;
                    r_current_displacement[1] =  r_previous_displacement[1] + delta_displacement[1];
                    predicted_y = true;
                }
            }
            if (disppos_y > -1 && !predicted_y) {
                if (!it_node->GetDof(DISPLACEMENT_Y, disppos_y).IsFixed() && !predicted_y) {
                    r_current_displacement[1] = r_previous_displacement[1] + delta_time * r_previous_velocity[1] + 0.5 * std::pow(delta_time, 2) * r_previous_acceleration[1];
                }
            }

            // For 3D cases
            if (accelpos_z > -1) {
                if (it_node->GetDof(ACCELERATION_Z, accelpos_z).IsFixed()) {
                    delta_displacement[2] = (r_current_acceleration[2] + mBossak.c3 * r_previous_acceleration[2] +  mBossak.c2 * r_previous_velocity[2])/mBossak.c0;
                    r_current_displacement[2] =  r_previous_displacement[2] + delta_displacement[2];
                    predicted_z = true;
                }
            }
            if (velpos_z > -1 && !predicted_z) {
                if (it_node->GetDof(VELOCITY_Z, velpos_z).IsFixed() && !predicted_z) {
                    delta_displacement[2] = (r_current_velocity[2] + mBossak.c4 * r_previous_velocity[2] + mBossak.c5 * r_previous_acceleration[2])/mBossak.c1;
                    r_current_displacement[2] =  r_previous_displacement[2] + delta_displacement[2];
                    predicted_z = true;
                }
            }
            if (disppos_z > -1 && !predicted_z) {
                if (!it_node->GetDof(DISPLACEMENT_Z, disppos_z).IsFixed() && !predicted_z) {
                    r_current_displacement[2] = r_previous_displacement[2] + delta_time * r_previous_velocity[2] + 0.5 * std::pow(delta_time, 2) * r_previous_acceleration[2];
                }
            }

            // Updating time derivatives ::: Please note that displacements and its time derivatives can not be consistently fixed separately
            noalias(delta_displacement) = r_current_displacement - r_previous_displacement;

            UpdateVelocity(it_node, r_current_velocity, delta_displacement, r_previous_velocity, r_previous_acceleration);

            UpdateAcceleration (it_node, r_current_acceleration, delta_displacement, r_previous_velocity, r_previous_acceleration);
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief It initializes time step solution. Only for reasons if the time step solution is restarted
     * @param rModelPart The model of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        ImplicitBaseType::InitializeSolutionStep(rModelPart, A, Dx, b);

        const double delta_time = r_current_process_info[DELTA_TIME];

        // Initializing Bossak constants
        mBossak.c0 = ( 1.0 / (mBossak.beta * delta_time * delta_time) );
        mBossak.c1 = ( mBossak.gamma / (mBossak.beta * delta_time) );
        mBossak.c2 = ( 1.0 / (mBossak.beta * delta_time) );
        mBossak.c3 = ( 0.5 / (mBossak.beta) - 1.0 );
        mBossak.c4 = ( (mBossak.gamma / mBossak.beta) - 1.0  );
        mBossak.c5 = ( delta_time * 0.5 * ( ( mBossak.gamma / mBossak.beta ) - 2.0 ) );

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );
        const auto it_node_begin = rModelPart.Nodes().begin();

        // Getting position
        const int disppos_x = it_node_begin->HasDofFor(DISPLACEMENT_X) ? it_node_begin->GetDofPosition(DISPLACEMENT_X) : -1;
        const int velpos_x = it_node_begin->HasDofFor(VELOCITY_X) ? it_node_begin->GetDofPosition(VELOCITY_X) : -1;
        const int accelpos_x = it_node_begin->HasDofFor(ACCELERATION_X) ? it_node_begin->GetDofPosition(ACCELERATION_X) : -1;
        const int disppos_y = it_node_begin->HasDofFor(DISPLACEMENT_Y) ? it_node_begin->GetDofPosition(DISPLACEMENT_Y) : -1;
        const int velpos_y = it_node_begin->HasDofFor(VELOCITY_Y) ? it_node_begin->GetDofPosition(VELOCITY_Y) : -1;
        const int accelpos_y = it_node_begin->HasDofFor(ACCELERATION_Y) ? it_node_begin->GetDofPosition(ACCELERATION_Y) : -1;
        const int disppos_z = it_node_begin->HasDofFor(DISPLACEMENT_Z) ? it_node_begin->GetDofPosition(DISPLACEMENT_Z) : -1;
        const int velpos_z = it_node_begin->HasDofFor(VELOCITY_Z) ? it_node_begin->GetDofPosition(VELOCITY_Z) : -1;
        const int accelpos_z = it_node_begin->HasDofFor(ACCELERATION_Z) ? it_node_begin->GetDofPosition(ACCELERATION_Z) : -1;

        bool fixed_x, fixed_y, fixed_z;

        #pragma omp parallel for private(fixed_x, fixed_y, fixed_z)
        for(int i = 0;  i < num_nodes; ++i) {
            auto it_node = it_node_begin + i;

            fixed_x = false;
            fixed_y = false;
            fixed_z = false;

            if (accelpos_x > -1) {
                if (it_node->GetDof(ACCELERATION_X, accelpos_x).IsFixed()) {
                    if (disppos_x > -1) it_node->Fix(DISPLACEMENT_X);
                    fixed_x = true;
                }
            }
            if (velpos_x > -1 && !fixed_x) {
                if (it_node->GetDof(VELOCITY_X, velpos_x).IsFixed() && !fixed_x) {
                    if (disppos_x > -1) it_node->Fix(DISPLACEMENT_X);
                }
            }
            if (accelpos_y > -1) {
                if (it_node->GetDof(ACCELERATION_Y, accelpos_y).IsFixed()) {
                    if (disppos_y > -1) it_node->Fix(DISPLACEMENT_Y);
                    fixed_y = true;
                }
            }
            if (velpos_y > -1 && !fixed_y) {
                if (it_node->GetDof(VELOCITY_Y, velpos_y).IsFixed() && !fixed_y) {
                    if (disppos_y > -1) it_node->Fix(DISPLACEMENT_Y);
                }
            }
            if (accelpos_z > -1) {
                if (it_node->GetDof(ACCELERATION_Z, accelpos_z).IsFixed()) {
                    if (disppos_z > -1) it_node->Fix(DISPLACEMENT_Z);
                    fixed_z = true;
                }
            }
            if (velpos_z > -1 && !fixed_z) {
                if (it_node->GetDof(VELOCITY_Z, velpos_z).IsFixed() && !fixed_z) {
                    if (disppos_z > -1) it_node->Fix(DISPLACEMENT_Z);
                }
            }
        }

        KRATOS_CATCH( "" );
    }

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

        const int err = ImplicitBaseType::Check(rModelPart);
        if(err != 0) return err;

        // Check for variables keys
        // Verify that the variables are correctly initialized
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
        KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)

        // Check that variables are correctly allocated
        for (const auto& rnode : rModelPart.Nodes()) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
        }

        // Check for minimum value of the buffer index
        // Verify buffer size
        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2)
            << "Insufficient buffer size. Buffer size should be greater than 2. Current size is: "
            << rModelPart.GetBufferSize() << std::endl;

        // Check for admissible value of the AlphaBossak
        KRATOS_ERROR_IF(mBossak.alpha > 0.0 || mBossak.alpha < -0.5) << "Value not admissible for "
            << "AlphaBossak. Admissible values are between 0.0 and -0.5\nCurrent value is: "
            << mBossak.alpha << std::endl;

        static const double epsilon = 1e-12;
        KRATOS_ERROR_IF_NOT(std::abs(mNewmark.beta - 0.0)   < epsilon ||
                            std::abs(mNewmark.beta - 0.167) < epsilon ||
                            std::abs(mNewmark.beta - 0.25)  < epsilon)
            << "Value not admissible for NewmarkBeta. Admissible values are:\n"
            << "0.0 for central-differencing\n"
            << "0.25 for mean-constant-acceleration\n"
            << "0.167 for linear-acceleration\n"
            << "Current value is: " << mNewmark.beta << std::endl;

        return 0;
        KRATOS_CATCH( "" );
    }

    /// Free memory allocated by this class.
    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedBossakDisplacementScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    /**
     * @brief The Bossak Alpha components
     */
    struct BossakAlphaMethod
    {
        double alpha; /// Alpha Bossak
        double beta;  /// Beta Bossak
        double gamma; /// Gamma Bossak

        // System constants
        double c0, c1, c2, c3, c4, c5;
    };

    /**
     * @brief The Newmark parameters used during integration
     */
    struct NewmarkMethod
    {
        // Newmark constants
        double beta;  ///Beta Newmark
        double gamma; //Gamma Newmark
    };

    /**
     * @brief Vector containing the velocity and acceleration used on integration
     */
    struct GeneralVectors
    {
        std::vector< Vector > v;  /// Velocity
        std::vector< Vector > a;  /// Acceleration
        std::vector< Vector > ap; /// Previous acceleration
    };

    BossakAlphaMethod mBossak;     /// The structure containing the Bossak components
    NewmarkMethod mNewmark;        /// The structure containing the Newmark parameters
    GeneralVectors mVector;        /// The structure containing the velocities and accelerations

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Updating first time Derivative
     * @param CurrentVelocity The current velocity
     * @param DeltaDisplacement The increment of displacement
     * @param PreviousVelocity The previous velocity
     * @param PreviousAcceleration The previous acceleration
     */
    inline void UpdateVelocity(
        NodeIterator itNode,
        array_1d<double, 3>& rCurrentVelocity,
        const array_1d<double, 3>& rDeltaDisplacement,
        const array_1d<double, 3>& rPreviousVelocity,
        const array_1d<double, 3>& rPreviousAcceleration
        )
    {
        noalias(rCurrentVelocity) = (mBossak.c1 * rDeltaDisplacement - mBossak.c4 * rPreviousVelocity - mBossak.c5 * rPreviousAcceleration);
    }

    /**
     * @brief Updating second time Derivative
     * @param CurrentAcceleration The current velocity
     * @param DeltaDisplacement The increment of displacement
     * @param PreviousVelocity The previous velocity
     * @param PreviousAcceleration The previous acceleration
     */
    inline void UpdateAcceleration(
        NodeIterator itNode,
        array_1d<double, 3>& rCurrentAcceleration,
        const array_1d<double, 3>& rDeltaDisplacement,
        const array_1d<double, 3>& rPreviousVelocity,
        const array_1d<double, 3>& rPreviousAcceleration
        )
    {
        noalias(rCurrentAcceleration) = (mBossak.c0 * rDeltaDisplacement - mBossak.c2 * rPreviousVelocity -  mBossak.c3 * rPreviousAcceleration);
    }

    /**
     * @brief It adds the dynamic LHS contribution of the elements M*c0 + D*c1 + K
     * @param LHS_Contribution The dynamic contribution for the LHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // Adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
            noalias(LHS_Contribution) += M * (1.0 - mBossak.alpha) * mBossak.c0;

        // Adding  damping contribution
        if (D.size1() != 0) // if D matrix declared
            noalias(LHS_Contribution) += D * mBossak.c1;
    }

    /**
     * @brief It adds the dynamic RHS contribution of the elements b - (1-alpha)*M*a_n+1 - alpha*M*a_n - D*v_n
     * @param pElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Element::Pointer pElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (M.size1() != 0) {
            pElement->GetSecondDerivativesVector(mVector.a[this_thread], 0);
            mVector.a[this_thread] *= (1.00 - mBossak.alpha);

            pElement->GetSecondDerivativesVector(mVector.ap[this_thread], 1);
            noalias(mVector.a[this_thread]) += mBossak.alpha * mVector.ap[this_thread];

            noalias(RHS_Contribution) -= prod(M, mVector.a[this_thread]);
        }

        // Adding damping contribution
        if (D.size1() != 0) {
            pElement->GetFirstDerivativesVector(mVector.v[this_thread], 0);

            noalias(RHS_Contribution) -= prod(D, mVector.v[this_thread]);
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the condition b - (1-alpha)*M*a_n+1 - alpha*M*a_n - D*v_n
     * @param pCondition The condition to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Condition::Pointer pCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (M.size1() != 0) {
            pCondition->GetSecondDerivativesVector(mVector.a[this_thread], 0);
            mVector.a[this_thread] *= (1.00 - mBossak.alpha);

            pCondition->GetSecondDerivativesVector(mVector.ap[this_thread], 1);
            noalias(mVector.a[this_thread]) += mBossak.alpha * mVector.ap[this_thread];

            noalias(RHS_Contribution) -= prod(M, mVector.a[this_thread]);
        }

        // Adding damping contribution
        // Damping contribution
        if (D.size1() != 0) {
            pCondition->GetFirstDerivativesVector(mVector.v[this_thread], 0);

            noalias(RHS_Contribution) -= prod(D, mVector.v[this_thread]);
        }
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
    ///@{

private:

    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}
}; /* Class ResidualBasedBossakDisplacementScheme */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BOSSAK_DISPLACEMENT_SCHEME defined */
