//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Josep Maria Carbonell
//                   Vicente Mataix Ferrandiz
//                   Andreas Winterstein (refactoring)
//

#pragma once

// Project includes
#include "solving_strategies/schemes/residual_based_implicit_time_scheme.h"
#include "includes/variables.h"
#include "includes/checks.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/** @brief Bossak integration scheme (for linear and nonlinear dynamic problems) for displacements
 *  @class ResidualBasedBossakDisplacementScheme
 *  @ingroup KratosCore
 *
 *  @details This is a dynamic implicit scheme based on the Bossak algorithm for displacements.
 *           The Bossak Alpha parameter ranges from 0 to -0.5, and introduces damping.
 *           Implementation according to: "An alpha modification of Newmark's method; W.L. Wood, M. Bossak, O.C. Zienkiewicz;
 *           Numerical Methods in Engineering; 1980"
 *  @author Josep Maria Carbonell
 *  @author Vicente Mataix Ferrandiz
 *  @author Andreas Winterstein (refactoring)
 */
template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedBossakDisplacementScheme
    : public ResidualBasedImplicitTimeScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedBossakDisplacementScheme );

    /// Base type for the scheme
    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    /// Implicit base type for the scheme
    using ImplicitBaseType = ResidualBasedImplicitTimeScheme<TSparseSpace, TDenseSpace>;

    /// Class type for the scheme
    using ClassType = ResidualBasedBossakDisplacementScheme<TSparseSpace, TDenseSpace>;

    /// Data type used within the ImplicitBaseType
    using TDataType = typename ImplicitBaseType::TDataType;

    /// Array type for degrees of freedom within ImplicitBaseType
    using DofsArrayType = typename ImplicitBaseType::DofsArrayType;

    /// Vector type for degrees of freedom within an Element
    using DofsVectorType = typename Element::DofsVectorType;

    /// Type for the system matrix within ImplicitBaseType
    using TSystemMatrixType = typename ImplicitBaseType::TSystemMatrixType;

    /// Type for the system vector within ImplicitBaseType
    using TSystemVectorType = typename ImplicitBaseType::TSystemVectorType;

    /// Type for local system vectors within ImplicitBaseType
    using LocalSystemVectorType = typename ImplicitBaseType::LocalSystemVectorType;

    /// Type for local system matrices within ImplicitBaseType
    using LocalSystemMatrixType = typename ImplicitBaseType::LocalSystemMatrixType;

    /// Iterator for nodes in a ModelPart
    using NodeIterator = ModelPart::NodeIterator;

    /// Container type for nodes in a ModelPart
    using NodesArrayType = ModelPart::NodesContainerType;

    /// Container type for elements in a ModelPart
    using ElementsArrayType = ModelPart::ElementsContainerType;

    /// Container type for conditions in a ModelPart
    using ConditionsArrayType = ModelPart::ConditionsContainerType;

    /// Pointer type for the BaseType
    using BaseTypePointer = typename BaseType::Pointer;

    /// Component type as 'double'
    using ComponentType = double;

    ///@}
    ///@name Life Cycle
    ///@{

    /** @brief Construct from a @ref Parameters object.
     *  @param ThisParameters The parameters containing the configuration.
     */
    explicit ResidualBasedBossakDisplacementScheme(Parameters ThisParameters)
        : ImplicitBaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        // For pure Newmark Scheme
        mNewmark.gamma = 0.5;

        AuxiliarInitializeBossak();
    }

    /** @brief Constructor from a Bossak parameter.
     *  @param Alpha is the Bossak parameter. Default value is 0, which is the Newmark method.
     *  @note The Newmark beta parameter is set to 0.25, for mean constant acceleration.
     */
    explicit ResidualBasedBossakDisplacementScheme(const double Alpha = 0.0)
        : ResidualBasedBossakDisplacementScheme(Alpha, 0.25)
    {
    }

    /** @brief Constructor.
     *  @param Alpha The Bossak parameter.
     *  @param NewarkBeta The Newmark parameter.
     *  @note A Bossak parameter (@a Alpha) of 0 results in the Newmark method.
     *  @note A Newmark parameter (@a NewmarkBeta) of 0.25 results in mean constant acceleration.
     */
    explicit ResidualBasedBossakDisplacementScheme(const double Alpha, const double NewmarkBeta)
        :ImplicitBaseType()
    {
        // For pure Newmark Scheme
        mBossak.alpha = Alpha;
        mNewmark.beta = NewmarkBeta;
        mNewmark.gamma = 0.5;

        AuxiliarInitializeBossak();
    }

    explicit ResidualBasedBossakDisplacementScheme(ResidualBasedBossakDisplacementScheme& rOther)
        :ImplicitBaseType(rOther)
        ,mBossak(rOther.mBossak)
        ,mNewmark(rOther.mNewmark)
        ,mVector(rOther.mVector)
    {
    }

    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new ResidualBasedBossakDisplacementScheme(*this) );
    }

    ~ResidualBasedBossakDisplacementScheme
    () override {}

    ///@}
    ///@name Operations
    ///@{

    /** @brief Construct a dynamically allocated new scheme from a @ref Parameters object.
     *  @param ThisParameters The configuration parameters.
     */
    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /// @brief Recalculate the Newmark coefficients, taking the alpha parameters into account.
    void CalculateBossakCoefficients()
    {
        mBossak.beta  = (1.0 - mBossak.alpha) * (1.0 - mBossak.alpha) * mNewmark.beta;
        mBossak.gamma = mNewmark.gamma  - mBossak.alpha;
    }

    /** @brief Update state variables within a newton iteration at the end of the time step.
     *  @details @f[ u_{n+1}^{k+1} = u_{n+1}^k+ \Delta u @f]
     *  @param rModelPart @ref ModelPart to update.
     *  @param rDofSet Set of all primary variables.
     *  @param rA Left hand side matrix.
     *  @param rDx Primary variable updates.
     *  @param rb Right hand side vector.
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
        block_for_each(rModelPart.Nodes(), array_1d<double,3>(), [&](Node& rNode, array_1d<double,3>& rDeltaDisplacementTLS){
            noalias(rDeltaDisplacementTLS) = rNode.FastGetSolutionStepValue(DISPLACEMENT) - rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);

            array_1d<double, 3>& r_current_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& r_previous_velocity = rNode.FastGetSolutionStepValue(VELOCITY, 1);

            array_1d<double, 3>& r_current_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION);
            const array_1d<double, 3>& r_previous_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

            UpdateVelocity(r_current_velocity, rDeltaDisplacementTLS, r_previous_velocity, r_previous_acceleration);
            UpdateAcceleration(r_current_acceleration, rDeltaDisplacementTLS, r_previous_velocity, r_previous_acceleration);
        });

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Apply the predictor.
     * @details @f[ x_{k+1} = x_{k} + v_{k} \cdot \Delta t @f]
     * @param rModelPart @ref ModelPart to update.
     * @param rDofSet Set of all primary variables.
     * @param rA Left hand side matrix.
     * @param rDx Primary variable updates.
     * @param rb Right hand side vector.
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
        if (rModelPart.Nodes().size() > 0) {
            const auto it_node_begin = rModelPart.Nodes().begin();

            // Getting position
            KRATOS_ERROR_IF_NOT(it_node_begin->HasDofFor(DISPLACEMENT_X)) << "ResidualBasedBossakDisplacementScheme:: DISPLACEMENT is not added" << std::endl;
            const int disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);
            const int velpos = it_node_begin->HasDofFor(VELOCITY_X) ? static_cast<int>(it_node_begin->GetDofPosition(VELOCITY_X)) : -1;
            const int accelpos = it_node_begin->HasDofFor(ACCELERATION_X) ? static_cast<int>(it_node_begin->GetDofPosition(ACCELERATION_X)) : -1;

            // Getting dimension
            KRATOS_WARNING_IF("ResidualBasedBossakDisplacementScheme", !r_current_process_info.Has(DOMAIN_SIZE)) << "DOMAIN_SIZE not defined. Please define DOMAIN_SIZE. 3D case will be assumed" << std::endl;
            const std::size_t dimension = r_current_process_info.Has(DOMAIN_SIZE) ? r_current_process_info.GetValue(DOMAIN_SIZE) : 3;

            // Auxiliar variables
            array_1d<double, 3 > delta_displacement;
            std::array<bool, 3> predicted = {false, false, false};
            const std::array<const Variable<ComponentType>*, 3> disp_components = {&DISPLACEMENT_X, &DISPLACEMENT_Y, &DISPLACEMENT_Z};
            const std::array<const Variable<ComponentType>*, 3> vel_components = {&VELOCITY_X, &VELOCITY_Y, &VELOCITY_Z};
            const std::array<const Variable<ComponentType>*, 3> accel_components = {&ACCELERATION_X, &ACCELERATION_Y, &ACCELERATION_Z};

            typedef std::tuple<array_1d<double,3>, std::array<bool,3>> TLSContainerType;
            TLSContainerType aux_TLS = std::make_tuple(delta_displacement, predicted);

            block_for_each(rModelPart.Nodes(), aux_TLS, [&](Node& rNode, TLSContainerType& rAuxTLS){
                auto& r_delta_displacement = std::get<0>(rAuxTLS);
                auto& r_predicted = std::get<1>(rAuxTLS);

                for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                    r_predicted[i_dim] = false;
                }

                // Predicting: NewDisplacement = r_previous_displacement + r_previous_velocity * delta_time;
                const array_1d<double, 3>& r_previous_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);
                const array_1d<double, 3>& r_previous_velocity     = rNode.FastGetSolutionStepValue(VELOCITY,     1);
                const array_1d<double, 3>& r_previous_displacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
                array_1d<double, 3>& r_current_acceleration        = rNode.FastGetSolutionStepValue(ACCELERATION);
                array_1d<double, 3>& r_current_velocity            = rNode.FastGetSolutionStepValue(VELOCITY);
                array_1d<double, 3>& r_current_displacement        = rNode.FastGetSolutionStepValue(DISPLACEMENT);

                if (accelpos > -1) {
                    for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                        if (rNode.GetDof(*accel_components[i_dim], accelpos + i_dim).IsFixed()) {
                            r_delta_displacement[i_dim] = (r_current_acceleration[i_dim] + mBossak.c3 * r_previous_acceleration[i_dim] +  mBossak.c2 * r_previous_velocity[i_dim])/mBossak.c0;
                            r_current_displacement[i_dim] =  r_previous_displacement[i_dim] + r_delta_displacement[i_dim];
                            r_predicted[i_dim] = true;
                        }
                    }
                }
                if (velpos > -1) {
                    for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                        if (rNode.GetDof(*vel_components[i_dim], velpos + i_dim).IsFixed() && !r_predicted[i_dim]) {
                            r_delta_displacement[i_dim] = (r_current_velocity[i_dim] + mBossak.c4 * r_previous_velocity[i_dim] + mBossak.c5 * r_previous_acceleration[i_dim])/mBossak.c1;
                            r_current_displacement[i_dim] =  r_previous_displacement[i_dim] + r_delta_displacement[i_dim];
                            r_predicted[i_dim] = true;
                        }
                    }
                }
                for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                    if (!rNode.GetDof(*disp_components[i_dim], disppos + i_dim).IsFixed() && !r_predicted[i_dim]) {
                        r_current_displacement[i_dim] = r_previous_displacement[i_dim] + delta_time * r_previous_velocity[i_dim] + 0.5 * std::pow(delta_time, 2) * r_previous_acceleration[i_dim];
                    }
                }

                // Updating time derivatives ::: Please note that displacements and its time derivatives can not be consistently fixed separately
                noalias(r_delta_displacement) = r_current_displacement - r_previous_displacement;
                UpdateVelocity(r_current_velocity, r_delta_displacement, r_previous_velocity, r_previous_acceleration);
                UpdateAcceleration(r_current_acceleration, r_delta_displacement, r_previous_velocity, r_previous_acceleration);
            });
        }

        KRATOS_CATCH( "" );
    }

    /** @brief Prepare state variables for a new solution step.
     *  @note This function should be called before solving a solution step,
     *        or restarting a failed one.
     *  @param rModelPart @ref ModelPart to update.
     *  @param A Left hand side matrix.
     *  @param Dx Primary variable updates.
     *  @param b Right hand side vector.
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
        if (rModelPart.Nodes().size() > 0) {
            const auto it_node_begin = rModelPart.Nodes().begin();

            // Getting dimension
            KRATOS_WARNING_IF("ResidualBasedBossakDisplacementScheme", !r_current_process_info.Has(DOMAIN_SIZE)) << "DOMAIN_SIZE not defined. Please define DOMAIN_SIZE. 3D case will be assumed" << std::endl;
            const std::size_t dimension = r_current_process_info.Has(DOMAIN_SIZE) ? r_current_process_info.GetValue(DOMAIN_SIZE) : 3;

            // Getting position
            const int velpos = it_node_begin->HasDofFor(VELOCITY_X) ? static_cast<int>(it_node_begin->GetDofPosition(VELOCITY_X)) : -1;
            const int accelpos = it_node_begin->HasDofFor(ACCELERATION_X) ? static_cast<int>(it_node_begin->GetDofPosition(ACCELERATION_X)) : -1;

            std::array<bool, 3> fixed = {false, false, false};
            const std::array<const Variable<ComponentType>*, 3> disp_components = {&DISPLACEMENT_X, &DISPLACEMENT_Y, &DISPLACEMENT_Z};
            const std::array<const Variable<ComponentType>*, 3> vel_components = {&VELOCITY_X, &VELOCITY_Y, &VELOCITY_Z};
            const std::array<const Variable<ComponentType>*, 3> accel_components = {&ACCELERATION_X, &ACCELERATION_Y, &ACCELERATION_Z};

            block_for_each(rModelPart.Nodes(), fixed, [&](Node& rNode, std::array<bool,3>& rFixedTLS){
                for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                    rFixedTLS[i_dim] = false;
                }

                if (accelpos > -1) {
                    for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                        if (rNode.GetDof(*accel_components[i_dim], accelpos + i_dim).IsFixed()) {
                            rNode.Fix(*disp_components[i_dim]);
                            rFixedTLS[i_dim] = true;
                        }
                    }
                }
                if (velpos > -1) {
                    for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                        if (rNode.GetDof(*vel_components[i_dim], velpos + i_dim).IsFixed() && !rFixedTLS[i_dim]) {
                            rNode.Fix(*disp_components[i_dim]);
                        }
                    }
                }
            });
        }

        KRATOS_CATCH( "" );
    }

    /** @brief Check whether the scheme and the provided @ref ModelPart are configured correctly.
     *  @details Checks may be expensive as the function is designed to catch user errors.
     *  @param rModelPart @ref ModelPart to check.
     *  @return 0 if ok, nonzero otherwise.
     */
    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY;

        const int err = ImplicitBaseType::Check(rModelPart);
        if(err != 0) return err;

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

    /// @brief Release dynamic memory allocated by this instance.
    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    /// @brief This function returns the default @ref Parameters to help avoiding conflicts between different constructors.
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"          : "bossak_scheme",
            "damp_factor_m" : -0.3,
            "newmark_beta"  : 0.25
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = ImplicitBaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /// @brief Return the name of the class as used in the settings (@a snake_case).
    static std::string Name()
    {
        return "bossak_scheme";
    }

    ///@}
    ///@name Input and output
    ///@{

    /// @brief Return information as a string.
    std::string Info() const override
    {
        return "ResidualBasedBossakDisplacementScheme";
    }

    /// @brief Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// @brief Print the instance's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}

protected:
    /// @todo Move to @ref ImplicitBaseType
    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

    /// @brief Bossak Alpha parameters.
    struct BossakAlphaMethod
    {
        /// @brief Bossak Alpha.
        double alpha;

        /// @brief Bossak Beta.
        double beta;

        /// @brief Bossak Gamma.
        double gamma;

        /// @brief System constants
        double c0, c1, c2, c3, c4, c5;
    };

    /// @brief Newmark parameters used for integration.
    struct NewmarkMethod
    {
        /// @brief Newmark Beta.
        double beta;

        /// @brief Newmark Gamma.
        double gamma;
    };

    /// @brief Velocities and accelerations used for integration.
    struct GeneralVectors
    {
        /// @brief Velocity.
        std::vector< Vector > v;

        /// @brief Acceleration.
        std::vector< Vector > a;

        /// @brief Previous acceleration.
        std::vector< Vector > ap;
    };

    /// @brief Bossak Alpha parameters.
    BossakAlphaMethod mBossak;

    /// @brief Newmark Beta parameters.
    NewmarkMethod mNewmark;

    /// @brief Aggregate struct for velocities and accelerations.
    GeneralVectors mVector;

    ///@name Protected Operations
    ///@{

    /** @brief Update the first time derivative.
     *  @param rCurrentVelocity Current velocity.
     *  @param rDeltaDisplacement Displacement increment.
     *  @param rPreviousVelocity Previous velocity.
     *  @param rPreviousAcceleration Previous acceleration.
     */
    inline void UpdateVelocity(
        array_1d<double, 3>& rCurrentVelocity,
        const array_1d<double, 3>& rDeltaDisplacement,
        const array_1d<double, 3>& rPreviousVelocity,
        const array_1d<double, 3>& rPreviousAcceleration
        )
    {
        noalias(rCurrentVelocity) = (mBossak.c1 * rDeltaDisplacement - mBossak.c4 * rPreviousVelocity - mBossak.c5 * rPreviousAcceleration);
    }

    /** @brief Update the second time derivative.
     *  @param rCurrentAcceleration Current velocity.
     *  @param rDeltaDisplacement Displacement increment.
     *  @param rPreviousVelocity Previous velocity.
     *  @param rPreviousAcceleration Previous acceleration.
     */
    inline void UpdateAcceleration(
        array_1d<double, 3>& rCurrentAcceleration,
        const array_1d<double, 3>& rDeltaDisplacement,
        const array_1d<double, 3>& rPreviousVelocity,
        const array_1d<double, 3>& rPreviousAcceleration
        )
    {
        noalias(rCurrentAcceleration) = (mBossak.c0 * rDeltaDisplacement - mBossak.c2 * rPreviousVelocity -  mBossak.c3 * rPreviousAcceleration);
    }

    /** @brief Add dynamic left hand side contribution from @ref Element s.
     *  @details @f[ M \cdot c_0 + D \cdot c_1 + K @f]
     *  @param LHS_Contribution Dynamic contribution for the left hand side.
     *  @param D Damping matrix.
     *  @param M Mass matrix.
     *  @param rCurrentProcessInfo Current @ref ProcessInfo.
     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // Adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
            noalias(LHS_Contribution) += M * (1.0 - mBossak.alpha) * mBossak.c0;

        // Adding  damping contribution
        if (D.size1() != 0) // if D matrix declared
            noalias(LHS_Contribution) += D * mBossak.c1;
    }

    /** @brief Add dynamic right hand side contribution from an @ref Element.
     *  @details @f[ b - (1 - \alpha) \cdot M \cdot a_{n+1} - \alpha \cdot M \cdot a_n - D \cdot v_n @f]
     *  @param rElement @ref Element to compute the contribution from.
     *  @param RHS_Contribution Dynamic contribution for the right hand side.
     *  @param D Damping matrix.
     *  @param M Mass matrix.
     *  @param rCurrentProcessInfo Current @ref ProcessInfo.
     */
    void AddDynamicsToRHS(
        Element& rElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        const auto& r_const_elem_ref = rElement;
        // Adding inertia contribution
        if (M.size1() != 0) {

            r_const_elem_ref.GetSecondDerivativesVector(mVector.a[this_thread], 0);
            mVector.a[this_thread] *= (1.00 - mBossak.alpha);

            r_const_elem_ref.GetSecondDerivativesVector(mVector.ap[this_thread], 1);
            noalias(mVector.a[this_thread]) += mBossak.alpha * mVector.ap[this_thread];

            noalias(RHS_Contribution) -= prod(M, mVector.a[this_thread]);
        }

        // Adding damping contribution
        if (D.size1() != 0) {
            r_const_elem_ref.GetFirstDerivativesVector(mVector.v[this_thread], 0);
            noalias(RHS_Contribution) -= prod(D, mVector.v[this_thread]);
        }
    }

    /** @brief Add dynamic right hand side contribution of a @ref Condition.
     *  @details @f[ b - (1-alpha)*M*a_n+1 - alpha*M*a_n - D*v_n @f]
     *  @param rCondition @ref Condition to compute the contribution from.
     *  @param RHS_Contribution Dynamic contribution for the right hand side.
     *  @param D Damping matrix.
     *  @param M Mass matrix.
     *  @param rCurrentProcessInfo Current @ref ProcessInfo.
     */
    void AddDynamicsToRHS(
        Condition& rCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();
        const auto& r_const_cond_ref = rCondition;

        // Adding inertia contribution
        if (M.size1() != 0) {
            r_const_cond_ref.GetSecondDerivativesVector(mVector.a[this_thread], 0);
            mVector.a[this_thread] *= (1.00 - mBossak.alpha);

            r_const_cond_ref.GetSecondDerivativesVector(mVector.ap[this_thread], 1);
            noalias(mVector.a[this_thread]) += mBossak.alpha * mVector.ap[this_thread];

            noalias(RHS_Contribution) -= prod(M, mVector.a[this_thread]);
        }

        // Adding damping contribution
        // Damping contribution
        if (D.size1() != 0) {
            r_const_cond_ref.GetFirstDerivativesVector(mVector.v[this_thread], 0);

            noalias(RHS_Contribution) -= prod(D, mVector.v[this_thread]);
        }
    }

    /// @brief Assign member variables from @ref Parameters.
    void AssignSettings(const Parameters ThisParameters) override
    {
        ImplicitBaseType::AssignSettings(ThisParameters);
        mBossak.alpha = ThisParameters["damp_factor_m"].GetDouble();
        mNewmark.beta = ThisParameters["newmark_beta"].GetDouble();
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    /// @brief Utility function for initializing some members.
    void AuxiliarInitializeBossak()
    {
        // Initialize Bossak coefficients
        CalculateBossakCoefficients();

        // Allocate auxiliary memory
        const std::size_t num_threads = ParallelUtilities::GetNumThreads();

        mVector.v.resize(num_threads);
        mVector.a.resize(num_threads);
        mVector.ap.resize(num_threads);

        KRATOS_DETAIL("MECHANICAL SCHEME: The Bossak Time Integration Scheme ") << "[alpha_m= " << mBossak.alpha << " beta= " << mNewmark.beta << " gamma= " << mNewmark.gamma << "]" <<std::endl;
    }

    ///@}
}; // Class ResidualBasedBossakDisplacementScheme

///@}

} // namespace Kratos
