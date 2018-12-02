//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Main authors:  Josep Maria Carbonell
//                 Vicente Mataix Ferrandiz
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
 * @brief Bossak integration scheme (for dynamic problems) for displacements
 * @details This is a dynamic implicit scheme based of the Bossak algorithm for displacements.
 * The parameter Alpha of Bossak introduces damping, the value of Bossak is from 0 to -0.3 (negative)
 * @author Josep Maria Carbonell
 * @author Vicente Mataix Ferrandiz
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
        : ResidualBasedBossakDisplacementScheme(ThisParameters.Has("damp_factor_m") ? ThisParameters["damp_factor_m"].GetDouble() : -0.3)
    {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "damp_factor_m" : -0.3
        })" );
        ThisParameters.ValidateAndAssignDefaults(default_parameters);
    }

    /**
     * @brief Constructor.
     * @detail The bossak method
     * @param rAlpham The Bossak parameter. Default value is 0, which is the Newmark method
     */
    explicit ResidualBasedBossakDisplacementScheme(const double rAlpham = 0.0)
        :ImplicitBaseType()
    {
        // For pure Newmark Scheme
        mAlpha.f = 0.0;
        mAlpha.m = rAlpham;

        // Default values of the Newmark coefficients
        double beta  = 0.25;
        double gamma = 0.5;

        CalculateNewmarkCoefficients(beta, gamma);

        // Allocate auxiliary memory
        const std::size_t num_threads = OpenMPUtils::GetNumThreads();

        mVector.v.resize(num_threads);
        mVector.a.resize(num_threads);
        mVector.ap.resize(num_threads);

        KRATOS_DETAIL("MECHANICAL SCHEME: The Bossak Time Integration Scheme ") << "[alpha_m= " << mAlpha.m << " beta= " << mNewmark.beta << " gamma= " << mNewmark.gamma << "]" <<std::endl;
    }

    /**
     * @brief Copy Constructor.
     */
    explicit ResidualBasedBossakDisplacementScheme(ResidualBasedBossakDisplacementScheme& rOther)
        :ImplicitBaseType(rOther)
        ,mAlpha(rOther.mAlpha)
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
     * @param beta The Newmark beta coefficient
     * @param gamma The Newmark gamma coefficient
     */
    void CalculateNewmarkCoefficients(
            double beta,
            double gamma
            )
    {
        mNewmark.beta  = (1.0 + mAlpha.f - mAlpha.m) * (1.0 + mAlpha.f - mAlpha.m) * beta;
        mNewmark.gamma = gamma + mAlpha.f - mAlpha.m;
    }

    /**
     * @brief Performing the update of the solution
     * @details Incremental update within newton iteration. It updates the state variables at the end of the time step u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param rModelPart The model of the problem to solve
     * @param rDofSet Set of all primary variables
     * @param A LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;

        mpDofUpdater->UpdateDofs(rDofSet,Dx);

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>(rModelPart.NumberOfNodes());

        #pragma omp parallel for
        for(int i = 0;  i < num_nodes; ++i) {
            auto it_node = rModelPart.Nodes().begin() + i;

            array_1d<double, 3 > delta_displacement;

            noalias(delta_displacement) = it_node->FastGetSolutionStepValue(DISPLACEMENT) - it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);

            array_1d<double, 3>& current_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& previous_velocity = it_node->FastGetSolutionStepValue(VELOCITY, 1);

            array_1d<double, 3>& current_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION);
            const array_1d<double, 3>& previous_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION, 1);

            UpdateVelocity(current_velocity, delta_displacement, previous_velocity, previous_acceleration);
            UpdateAcceleration(current_acceleration, delta_displacement, previous_velocity, previous_acceleration);
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Performing the prediction of the solution
     * @details It predicts the solution for the current step x = xold + vold * Dt
     * @param rModelPart The model of the problem to solve
     * @param rDofSet set of all primary variables
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;

        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>(rModelPart.NumberOfNodes());

        array_1d<double, 3 > delta_displacement;

        #pragma omp parallel for private(delta_displacement)
        for(int i = 0;  i < num_nodes; ++i) {
            auto it_node = rModelPart.Nodes().begin() + i;

            //Predicting: NewDisplacement = previous_displacement + previous_velocity * delta_time;
            //ATTENTION::: the prediction is performed only on free nodes

            const array_1d<double, 3 > & previous_acceleration = (it_node)->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double, 3 > & previous_velocity     = (it_node)->FastGetSolutionStepValue(VELOCITY,     1);
            const array_1d<double, 3 > & previous_displacement = (it_node)->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double, 3 > & current_acceleration        = (it_node)->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3 > & current_velocity            = (it_node)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & current_displacement        = (it_node)->FastGetSolutionStepValue(DISPLACEMENT);

            if (it_node -> IsFixed(ACCELERATION_X)) {
                current_displacement[0] = previous_displacement[0] + delta_time * previous_velocity[0] + std::pow(delta_time, 2) * ( 0.5 * (1.0 -  2.0 * mNewmark.beta) * previous_acceleration[0] + mNewmark.beta * current_acceleration[0]);
            } else if (it_node -> IsFixed(VELOCITY_X)) {
                current_displacement[0] = previous_displacement[0] + 0.5 * delta_time * (previous_velocity[0] + current_velocity[0]) + 0.5 * std::pow(delta_time, 2) * previous_acceleration[0];
            } else if (it_node -> IsFixed(DISPLACEMENT_X) == false) {
                current_displacement[0] = previous_displacement[0] + delta_time * previous_velocity[0] + 0.5 * std::pow(delta_time, 2) * previous_acceleration[0];
            }

            if (it_node -> IsFixed(ACCELERATION_Y)) {
                current_displacement[1] = previous_displacement[1] + delta_time * previous_velocity[1] + std::pow(delta_time, 2) * ( 0.5 * (1.0 -  2.0 * mNewmark.beta) * previous_acceleration[1] + mNewmark.beta * current_acceleration[1]);
            } else if (it_node -> IsFixed(VELOCITY_Y)) {
                current_displacement[1] = previous_displacement[1] + 0.5 * delta_time * (previous_velocity[1] + current_velocity[1]) + 0.5 * std::pow(delta_time, 2) * previous_acceleration[1] ;
            } else if (it_node -> IsFixed(DISPLACEMENT_Y) == false) {
                current_displacement[1] = previous_displacement[1] + delta_time * previous_velocity[1] + 0.5 * std::pow(delta_time, 2) * previous_acceleration[1];
            }

            // For 3D cases
            if (it_node -> HasDofFor(DISPLACEMENT_Z)) {
                if (it_node -> IsFixed(ACCELERATION_Z)) {
                    current_displacement[2] = previous_displacement[2] + delta_time * previous_velocity[2] + std::pow(delta_time, 2) * ( 0.5 * (1.0 -  2.0 * mNewmark.beta) * previous_acceleration[2] + mNewmark.beta * current_acceleration[2]);
                } else if (it_node -> IsFixed(VELOCITY_Z)) {
                    current_displacement[2] = previous_displacement[2] + 0.5 * delta_time * (previous_velocity[2] + current_velocity[2]) + 0.5 * std::pow(delta_time, 2) * previous_acceleration[2] ;
                } else if (it_node -> IsFixed(DISPLACEMENT_Z) == false) {
                    current_displacement[2] = previous_displacement[2] + delta_time * previous_velocity[2] + 0.5 * std::pow(delta_time, 2) * previous_acceleration[2];
                }
            }


            // Updating time derivatives ::: Please note that displacements and its time derivatives can not be consistently fixed separately
            noalias(delta_displacement) = current_displacement - previous_displacement;

            UpdateVelocity     (current_velocity,     delta_displacement, previous_velocity, previous_acceleration);

            UpdateAcceleration (current_acceleration, delta_displacement, previous_velocity, previous_acceleration);
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

        ProcessInfo& current_process_info= rModelPart.GetProcessInfo();

        ImplicitBaseType::InitializeSolutionStep(rModelPart, A, Dx, b);

        const double delta_time = current_process_info[DELTA_TIME];

        double beta = 0.25;
        if (current_process_info.Has(NEWMARK_BETA))
            beta = current_process_info[NEWMARK_BETA];
        double gamma = 0.5;
        if (current_process_info.Has(NEWMARK_GAMMA))
            gamma = current_process_info[NEWMARK_GAMMA];

        CalculateNewmarkCoefficients(beta, gamma);

        // Initializing Newmark constants
        mNewmark.c0 = ( 1.0 / (mNewmark.beta * std::pow(delta_time, 2)) );
        mNewmark.c1 = ( mNewmark.gamma / (mNewmark.beta * delta_time) );
        mNewmark.c2 = ( 1.0 / (mNewmark.beta * delta_time) );
        mNewmark.c3 = ( 0.5 / (mNewmark.beta) - 1.0 );
        mNewmark.c4 = ( (mNewmark.gamma / mNewmark.beta) - 1.0  );
        mNewmark.c5 = ( delta_time * 0.5 * ( ( mNewmark.gamma / mNewmark.beta ) - 2.0 ) );

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
        for(auto& rnode : rModelPart.Nodes()) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
        }

        // Check for minimum value of the buffer index
        // Verify buffer size
        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2) << "Insufficient buffer size. Buffer size should be greater than 2. Current size is" << rModelPart.GetBufferSize() << std::endl;

        // Check for admissible value of the AlphaBossak
        KRATOS_ERROR_IF(mAlpha.m > 0.0 || mAlpha.m < -0.3) << "Value not admissible for AlphaBossak. Admissible values should be between 0.0 and -0.3. Current value is " << mAlpha.m << std::endl;

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
     * @brief The Generalized Alpha components
     * @detail For more about it:
     * J. Chung, G.M.Hubert. "A Time Integration Algorithm for Structural Dynamics with Improved Numerical Dissipation: The Generalized-α Method" ASME Journal of Applied Mechanics, 60, 371:375, 1993.
     */
    struct GeneralizedAlphaMethod
    {
        double f; /// Alpha Hilbert
        double m; /// Alpha Bosssak
    };

    /**
     * @brief The Newmark parameters used during integration
     */
    struct NewmarkMethod
    {
        // Newmark constants
        double beta, gamma;

        // System constants
        double c0, c1, c2, c3, c4, c5;
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

    GeneralizedAlphaMethod mAlpha; /// The structure containing the Generalized alpha components
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
        array_1d<double, 3>& CurrentVelocity,
        const array_1d<double, 3>& DeltaDisplacement,
        const array_1d<double, 3>& PreviousVelocity,
        const array_1d<double, 3>& PreviousAcceleration
        )
    {
        noalias(CurrentVelocity) = (mNewmark.c1 * DeltaDisplacement - mNewmark.c4 * PreviousVelocity
                                     - mNewmark.c5 * PreviousAcceleration);
    }

    /**
     * @brief Updating second time Derivative
     * @param CurrentAcceleration The current velocity
     * @param DeltaDisplacement The increment of displacement
     * @param PreviousVelocity The previous velocity
     * @param PreviousAcceleration The previous acceleration
     */
    inline void UpdateAcceleration(
        array_1d<double, 3>& CurrentAcceleration,
        const array_1d<double, 3>& DeltaDisplacement,
        const array_1d<double, 3>& PreviousVelocity,
        const array_1d<double, 3>& PreviousAcceleration
        )
    {
        noalias(CurrentAcceleration) = (mNewmark.c0 * DeltaDisplacement - mNewmark.c2 * PreviousVelocity
                                         -  mNewmark.c3 * PreviousAcceleration);
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
            noalias(LHS_Contribution) += M * (1.0 - mAlpha.m) * mNewmark.c0;

        // Adding  damping contribution
        if (D.size1() != 0) // if D matrix declared
            noalias(LHS_Contribution) += D * (1.0 - mAlpha.f) * mNewmark.c1;
    }

    /**
     * @brief It adds the dynamic RHS contribution of the elements b - M*a - D*v
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
            mVector.a[this_thread] *= (1.00 - mAlpha.m);

            pElement->GetSecondDerivativesVector(mVector.ap[this_thread], 1);
            noalias(mVector.a[this_thread]) += mAlpha.m * mVector.ap[this_thread];

            noalias(RHS_Contribution) -= prod(M, mVector.a[this_thread]);
        }

        // Adding damping contribution
        if (D.size1() != 0) {
            pElement->GetFirstDerivativesVector(mVector.v[this_thread], 0);

            noalias(RHS_Contribution) -= prod(D, mVector.v[this_thread]);
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the condition b - M*a - D*v
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
            mVector.a[this_thread] *= (1.00 - mAlpha.m);

            pCondition->GetSecondDerivativesVector(mVector.ap[this_thread], 1);
            noalias(mVector.a[this_thread]) += mAlpha.m * mVector.ap[this_thread];

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
