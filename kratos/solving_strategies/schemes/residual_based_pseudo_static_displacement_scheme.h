//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
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

    typedef ResidualBasedPseudoStaticDisplacementScheme<TSparseSpace, TDenseSpace> ClassType;

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

    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit ResidualBasedPseudoStaticDisplacementScheme()
        : DerivedBaseType(0.0),
          mpRayleighBeta(&NODAL_MAUX)
    {
    }

    /**
     * @brief Constructor. The pseudo static scheme (parameters)
     * @param ThisParameters Parameters with the Rayleigh variable
     */
    explicit ResidualBasedPseudoStaticDisplacementScheme(Parameters ThisParameters)
        : DerivedBaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Default constructor. The pseudo static scheme
     */
    explicit ResidualBasedPseudoStaticDisplacementScheme(const Variable<double>& RayleighBetaVariable)
        :DerivedBaseType(0.0),
        mpRayleighBeta(&RayleighBetaVariable)
    {
    }

    /** Copy Constructor.
     */
    explicit ResidualBasedPseudoStaticDisplacementScheme(ResidualBasedPseudoStaticDisplacementScheme& rOther)
        :DerivedBaseType(rOther),
        mpRayleighBeta(rOther.mpRayleighBeta)
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
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
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

        DerivedBaseType::mpDofUpdater->UpdateDofs(rDofSet, rDx);

        // Updating time derivatives (nodally for efficiency)
        array_1d<double, 3 > delta_displacement;
        block_for_each(rModelPart.Nodes(), delta_displacement, [&](Node& rNode, array_1d<double,3>& rDeltaDisplacementTLS){

            noalias(rDeltaDisplacementTLS) = rNode.FastGetSolutionStepValue(DISPLACEMENT) - rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);

            array_1d<double, 3>& r_current_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& r_previous_velocity = rNode.FastGetSolutionStepValue(VELOCITY, 1);

            noalias(r_current_velocity) = (this->mBossak.c1 * rDeltaDisplacementTLS - this->mBossak.c4 * r_previous_velocity);
        });

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

        // Process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Delta time
        const double delta_time = r_current_process_info[DELTA_TIME];

        // Updating time derivatives (nodally for efficiency)
        const auto it_node_begin = rModelPart.Nodes().begin();

        // Auxiliar variables
        const array_1d<double, 3> zero_array = ZeroVector(3);
        array_1d<double, 3 > delta_displacement = zero_array;

        // Getting position
        const int disppos_x = it_node_begin->HasDofFor(DISPLACEMENT_X) ? static_cast<int>(it_node_begin->GetDofPosition(DISPLACEMENT_X)) : -1;
        const int velpos_x = it_node_begin->HasDofFor(VELOCITY_X) ? static_cast<int>(it_node_begin->GetDofPosition(VELOCITY_X)) : -1;
        const int disppos_y = it_node_begin->HasDofFor(DISPLACEMENT_Y) ? static_cast<int>(it_node_begin->GetDofPosition(DISPLACEMENT_Y)) : -1;
        const int velpos_y = it_node_begin->HasDofFor(VELOCITY_Y) ? static_cast<int>(it_node_begin->GetDofPosition(VELOCITY_Y)) : -1;
        const int disppos_z = it_node_begin->HasDofFor(DISPLACEMENT_Z) ? static_cast<int>(it_node_begin->GetDofPosition(DISPLACEMENT_Z)) : -1;
        const int velpos_z = it_node_begin->HasDofFor(VELOCITY_Z) ? static_cast<int>(it_node_begin->GetDofPosition(VELOCITY_Z)) : -1;

        block_for_each(rModelPart.Nodes(), delta_displacement, [&](Node& rNode, array_1d<double,3>& rDeltaDisplacementTLS){

            bool predicted_x = false;
            bool predicted_y = false;
            bool predicted_z = false;

            //Predicting: r_current_displacement = r_previous_displacement + r_previous_velocity * delta_time;
            //ATTENTION::: the prediction is performed only on free nodes

            const array_1d<double, 3>& r_previous_velocity     = rNode.FastGetSolutionStepValue(VELOCITY,     1);
            const array_1d<double, 3>& r_previous_displacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double, 3>& r_current_acceleration        = rNode.FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3>& r_current_velocity            = rNode.FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3>& r_current_displacement        = rNode.FastGetSolutionStepValue(DISPLACEMENT);

            if (velpos_x > -1) {
                if (rNode.GetDof(VELOCITY_X, velpos_x).IsFixed()) {
                    rDeltaDisplacementTLS[0] = (r_current_velocity[0] + this->mBossak.c4 * r_previous_velocity[0])/this->mBossak.c1;
                    r_current_displacement[0] =  r_previous_displacement[0] + rDeltaDisplacementTLS[0];
                    predicted_x = true;
                }
            }
            if (disppos_x > -1 && !predicted_x) {
                if (!rNode.GetDof(DISPLACEMENT_X, disppos_x).IsFixed() && !predicted_x) {
                    r_current_displacement[0] = r_previous_displacement[0] + delta_time * r_previous_velocity[0];
                }
            }

            if (velpos_y > -1) {
                if (rNode.GetDof(VELOCITY_Y, velpos_y).IsFixed()) {
                    rDeltaDisplacementTLS[1] = (r_current_velocity[1] + this->mBossak.c4 * r_previous_velocity[1])/this->mBossak.c1;
                    r_current_displacement[1] =  r_previous_displacement[1] + rDeltaDisplacementTLS[1];
                    predicted_y = true;
                }
            }
            if (disppos_y > -1 && !predicted_y) {
                if (!rNode.GetDof(DISPLACEMENT_Y, disppos_y).IsFixed() && !predicted_y) {
                    r_current_displacement[1] = r_previous_displacement[1] + delta_time * r_previous_velocity[1];
                }
            }

            if (velpos_z > -1) {
                if (rNode.GetDof(VELOCITY_Z, velpos_z).IsFixed()) {
                    rDeltaDisplacementTLS[2] = (r_current_velocity[2] + this->mBossak.c4 * r_previous_velocity[2])/this->mBossak.c1;
                    r_current_displacement[2] =  r_previous_displacement[2] + rDeltaDisplacementTLS[2];
                    predicted_z = true;
                }
            }
            if (disppos_z > -1 && !predicted_z) {
                if (!rNode.GetDof(DISPLACEMENT_Z, disppos_z).IsFixed() && !predicted_z) {
                    r_current_displacement[2] = r_previous_displacement[2] + delta_time * r_previous_velocity[2];
                }
            }

            // Updating time derivatives
            noalias(r_current_acceleration) = zero_array;
            noalias(r_current_velocity) = r_previous_velocity;
        });

        KRATOS_CATCH( "" );
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                   : "pseudo_static_scheme",
            "rayleigh_beta_variable" : "RAYLEIGH_BETA"
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = DerivedBaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "pseudo_static_scheme";
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
        rOStream << Info() << ". Considering the following damping variable " << *mpRayleighBeta;
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
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // Adding  damping contribution
        if (rD.size1() != 0 && TDenseSpace::TwoNorm(rD) > ZeroTolerance) // if D matrix declared
            noalias(rLHSContribution) += rD * this->mBossak.c1;
        else if (rM.size1() != 0) {
            const double beta = rCurrentProcessInfo[*mpRayleighBeta];
            noalias(rLHSContribution) += rM * beta * this->mBossak.c1;
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the elements b - D*v
     * @param rElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Element& rElement,
        LocalSystemVectorType& rRHSContribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();
        const auto& r_const_elem_ref = rElement;
        // Adding damping contribution
        if (rD.size1() != 0 && TDenseSpace::TwoNorm(rD) > ZeroTolerance) {
            r_const_elem_ref.GetFirstDerivativesVector(this->mVector.v[this_thread], 0);
            noalias(rRHSContribution) -= prod(rD, this->mVector.v[this_thread]);
        } else if (rM.size1() != 0) {
            const double beta = rCurrentProcessInfo[*mpRayleighBeta];
            r_const_elem_ref.GetFirstDerivativesVector(this->mVector.v[this_thread], 0);
            noalias(rRHSContribution) -= beta * prod(rM, this->mVector.v[this_thread]);
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the condition b - M*a - D*v
     * @param rCondition The condition to compute
     * @param rRHSContribution The dynamic contribution for the RHS
     * @param rD The damping matrix
     * @param rM The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Condition& rCondition,
        LocalSystemVectorType& rRHSContribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();
        const auto& r_const_cond_ref = rCondition;
        // Adding damping contribution
        // Damping contribution
        if (rD.size1() != 0 && TDenseSpace::TwoNorm(rD) > ZeroTolerance) {
            r_const_cond_ref.GetFirstDerivativesVector(this->mVector.v[this_thread], 0);
            noalias(rRHSContribution) -= prod(rD, this->mVector.v[this_thread]);
        } else if (rM.size1() != 0) {
            const double beta = rCurrentProcessInfo[*mpRayleighBeta];
            r_const_cond_ref.GetFirstDerivativesVector(this->mVector.v[this_thread], 0);
            noalias(rRHSContribution) -= beta * prod(rM, this->mVector.v[this_thread]);
        }
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        DerivedBaseType::AssignSettings(ThisParameters);
        mpRayleighBeta = &KratosComponents<Variable<double>>::Get(ThisParameters["rayleigh_beta_variable"].GetString());
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

    const Variable<double>* mpRayleighBeta = nullptr; /// The Rayleigh Beta variable

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