//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Main authors:  Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_RESIDUAL_PSEUDO_STATIC_DISPLACEMENT_SCHEME )
#define  KRATOS_RESIDUAL_PSEUDO_STATIC_DISPLACEMENT_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"

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
 * @class ResidualBasedPseudoStaticDisplacementScheme
 * @ingroup KratosCore
 * @brief This is a pseudo-static scheme
 * @details For pseudo–static strategy: calculate the constant matrices D = Beta * M, "set" M = 0 after initializing the damping matrix
 * @note Based on Riccardo Rossi PhD Thesis: Light weight Structures: Structural Analysis and Coupling Issues
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedPseudoStaticDisplacementScheme
    : public ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedPseudoStaticDisplacementScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                                        BaseType;

    typedef typename BaseType::TDataType                                           TDataType;

    typedef typename BaseType::DofsArrayType                                   DofsArrayType;

    typedef typename Element::DofsVectorType                                  DofsVectorType;

    typedef typename BaseType::TSystemMatrixType                           TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                           TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType                   LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType                   LocalSystemMatrixType;

    typedef ModelPart::ElementsContainerType                               ElementsArrayType;

    typedef ModelPart::ConditionsContainerType                           ConditionsArrayType;

    typedef typename BaseType::Pointer                                       BaseTypePointer;

    typedef ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>  DerivedBaseType;

    typedef typename BaseType::LocalSystemComponents               LocalSystemComponentsType;

    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor. The pseudo static scheme (parameters)
     * @param ThisParameters Parameters with the Rayleigh variable
     */
    explicit ResidualBasedPseudoStaticDisplacementScheme(Parameters ThisParameters)
        : DerivedBaseType(0.0),
          mRayleighBeta(NODAL_MAUX)
    {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "name"                   : "ResidualBasedPseudoStaticDisplacementScheme",
            "rayleigh_beta_variable" : "RAYLEIGH_BETA"
        })" );
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        mRayleighBeta = KratosComponents<Variable<double>>::Get(ThisParameters["rayleigh_beta_variable"].GetString());
    }

    /**
     * @brief Default constructor. The pseudo static scheme
     */
    explicit ResidualBasedPseudoStaticDisplacementScheme(const Variable<double> RayleighBetaVariable)
        :DerivedBaseType(0.0),
        mRayleighBeta(RayleighBetaVariable)
    {
    }

    /** Copy Constructor.
     */
    explicit ResidualBasedPseudoStaticDisplacementScheme(ResidualBasedPseudoStaticDisplacementScheme& rOther)
        :DerivedBaseType(rOther),
        mRayleighBeta(rOther.mRayleighBeta)
    {
    }

    /**
     * Clone
     */
    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new ResidualBasedPseudoStaticDisplacementScheme(*this) );
    }

    /** Destructor.
     */
    ~ResidualBasedPseudoStaticDisplacementScheme() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY;

        // Process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Delta time
        const double delta_time = r_current_process_info[DELTA_TIME];

        // Updating time derivatives (nodally for efficiency)
        const auto it_node_begin = rModelPart.Nodes().begin();
        const int num_nodes = static_cast<int>(rModelPart.NumberOfNodes());

        array_1d<double, 3> zero_array = ZeroVector(3);
        array_1d<double, 3> delta_displacement = zero_array;

        #pragma omp parallel for private(zero_array, delta_displacement)
        for(int i = 0;  i < num_nodes; ++i) {
            auto it_node = it_node_begin + i;

            //Predicting: r_current_displacement = r_previous_displacement + r_previous_velocity * delta_time;
            //ATTENTION::: the prediction is performed only on free nodes

            const array_1d<double, 3>& r_previous_velocity     = (it_node)->FastGetSolutionStepValue(VELOCITY,     1);
            const array_1d<double, 3>& r_previous_displacement = (it_node)->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double, 3>& r_current_acceleration        = (it_node)->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3>& r_current_velocity            = (it_node)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3>& r_current_displacement        = (it_node)->FastGetSolutionStepValue(DISPLACEMENT);

            if (it_node->IsFixed(VELOCITY_X)) {
                delta_displacement[0] = (r_current_velocity[0] + DerivedBaseType::mBossak.c4 * r_previous_velocity[0])/DerivedBaseType::mBossak.c1;
                r_current_displacement[0] =  r_previous_displacement[0] + delta_displacement[0];
            } else if (!it_node->IsFixed(DISPLACEMENT_X)) {
                r_current_displacement[0] = r_previous_displacement[0] + delta_time * r_previous_velocity[0];
            }

            if (it_node->IsFixed(VELOCITY_Y)) {
                delta_displacement[1] = (r_current_velocity[1] + DerivedBaseType::mBossak.c4 * r_previous_velocity[1])/DerivedBaseType::mBossak.c1;
                r_current_displacement[1] =  r_previous_displacement[1] + delta_displacement[1];
            } else if (!it_node->IsFixed(DISPLACEMENT_Y)) {
                r_current_displacement[1] = r_previous_displacement[1] + delta_time * r_previous_velocity[1];
            }

            // For 3D cases
            if (it_node->HasDofFor(DISPLACEMENT_Z)) {
                if (it_node->IsFixed(VELOCITY_Z)) {
                    delta_displacement[2] = (r_current_velocity[2] + DerivedBaseType::mBossak.c4 * r_previous_velocity[2])/DerivedBaseType::mBossak.c1;
                    r_current_displacement[2] =  r_previous_displacement[2] + delta_displacement[2];
                } else if (!it_node->IsFixed(DISPLACEMENT_Z)) {
                    r_current_displacement[2] = r_previous_displacement[2] + delta_time * r_previous_velocity[2];
                }
            }

            // Updating time derivatives
            noalias(r_current_acceleration) = zero_array;

            // We check if the dofs are fixed
            if (!it_node->IsFixed(VELOCITY_X)) {
                r_current_velocity[0] = r_previous_velocity[0];
            }
            if (!it_node->IsFixed(VELOCITY_Y)) {
                r_current_velocity[1] = r_previous_velocity[1];
            }
            if (!it_node->IsFixed(VELOCITY_Z)) {
                r_current_velocity[2] = r_previous_velocity[2];
            }
        }

        KRATOS_CATCH( "" );
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
        return "ResidualBasedPseudoStaticDisplacementScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info() << ". Considering the following damping variable " << mRayleighBeta;
    }

    ///@}
    ///@name Friends
    ///@{

protected:
    ///@}
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Protected  Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief It adds the dynamic LHS contribution of the elements D*c1 + K
     * @param rLHSContribution The dynamic contribution for the LHS
     * @param rD The damping matrix
     * @param rM The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& rLHSContribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // Adding  damping contribution
        if (rD.size1() != 0 && TDenseSpace::TwoNorm(rD) > ZeroTolerance) // if D matrix declared
            noalias(rLHSContribution) += rD * DerivedBaseType::mBossak.c1;
        else if (rM.size1() != 0) {
            const double beta = rCurrentProcessInfo[mRayleighBeta];
            noalias(rLHSContribution) += rM * beta * DerivedBaseType::mBossak.c1;
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the elements b - D*v
     * @param pElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Element::Pointer pElement,
        LocalSystemVectorType& rRHSContribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding damping contribution
        if (rD.size1() != 0 && TDenseSpace::TwoNorm(rD) > ZeroTolerance) {
            pElement->GetFirstDerivativesVector(DerivedBaseType::mVector.v[this_thread], 0);
            noalias(rRHSContribution) -= prod(rD, DerivedBaseType::mVector.v[this_thread]);
        } else if (rM.size1() != 0) {
            const double beta = rCurrentProcessInfo[mRayleighBeta];
            pElement->GetFirstDerivativesVector(DerivedBaseType::mVector.v[this_thread], 0);
            noalias(rRHSContribution) -= beta * prod(rM, DerivedBaseType::mVector.v[this_thread]);
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the condition b - M*a - D*v
     * @param pCondition The condition to compute
     * @param rRHSContribution The dynamic contribution for the RHS
     * @param rD The damping matrix
     * @param rM The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Condition::Pointer pCondition,
        LocalSystemVectorType& rRHSContribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding damping contribution
        // Damping contribution
        if (rD.size1() != 0 && TDenseSpace::TwoNorm(rD) > ZeroTolerance) {
            pCondition->GetFirstDerivativesVector(DerivedBaseType::mVector.v[this_thread], 0);
            noalias(rRHSContribution) -= prod(rD, DerivedBaseType::mVector.v[this_thread]);
        } else if (rM.size1() != 0) {
            const double beta = rCurrentProcessInfo[mRayleighBeta];
            pCondition->GetFirstDerivativesVector(DerivedBaseType::mVector.v[this_thread], 0);
            noalias(rRHSContribution) -= beta * prod(rM, DerivedBaseType::mVector.v[this_thread]);
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

    Variable<double> mRayleighBeta; /// The Rayleigh Beta variable

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
}; /* Class ResidualBasedPseudoStaticDisplacementScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_PSEUDO_STATIC_DISPLACEMENT_SCHEME E defined */
