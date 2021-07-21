//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_TRAPEZOIDAL_SHALLOW_SCHEME_H_INCLUDED
#define KRATOS_TRAPEZOIDAL_SHALLOW_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "shallow_water_application_variables.h"
#include "solving_strategies/schemes/residual_based_implicit_time_scheme.h"

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
 * @class TrapezoidalShallowScheme
 * @ingroup KratosShallowWaterApplication
 * @brief Trapezoidal integration scheme (for dynamic problems)
 * @details For \f$\theta=0\f$ Forward Euler, \f$\theta=0.5\f$ Crank Nicolson, \f$\theta=1\f$ Backward Euler.
 * @author Miguel Maso Sotomayor
 */
template<class TSparseSpace, class TDenseSpace>
class TrapezoidalShallowScheme : public ResidualBasedImplicitTimeScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( TrapezoidalShallowScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                           SchemeType;

    typedef typename SchemeType::Pointer                        SchemeTypePointer;

    typedef ResidualBasedImplicitTimeScheme<TSparseSpace,TDenseSpace>    BaseType;

    typedef typename BaseType::DofsArrayType                     DofsArrayType;

    typedef typename BaseType::TSystemMatrixType             TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType             TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType     LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType     LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                          NodesArrayType;

    typedef typename ModelPart::NodeType                                 NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor
    explicit TrapezoidalShallowScheme(const std::size_t Theta = 2, bool UpdateVelocities = false)
        : BaseType()
        , mTheta(Theta)
        , mUpdateVelocities(UpdateVelocities)
    {}

    // Constructor with parameters
    explicit TrapezoidalShallowScheme(Parameters ThisParameters)
        : BaseType(ThisParameters)
    {}

    // Copy Constructor
    explicit TrapezoidalShallowScheme(TrapezoidalShallowScheme& rOther)
        : BaseType(rOther)
        , mTheta(rOther.mTheta)
        , mUpdateVelocities(rOther.mUpdateVelocities)
    {}

    /**
     * Clone
     */
    SchemeTypePointer Clone() override
    {
        return Kratos::make_shared<TrapezoidalShallowScheme>(*this);
    }

    // Destructor
    ~TrapezoidalShallowScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /**
     * @brief This function is designed to be called once to perform all the checks needed on the input provided.
     * @details Checks can be "expensive" as the function is designed to catch user's errors.
     * @param rModelPart The model of the problem to solve
     * @return Zero means  all ok
     */
    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY;

        const int err = ImplicitBaseType::Check(rModelPart);
        if (err != 0) return err;

        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2) << "Insufficient buffer size. Buffer size should be greater than 2. Current size is " << rModelPart.GetBufferSize() << std::endl;

        KRATOS_CATCH("TrapezoidalShallowScheme.Check");

        return 0;
    }

    /**
     * @brief Performing the update of the solution
     * @details Incremental update within newton iteration. It updates the state variables at the end of the time step
     * \f[ u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u\f]
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

        // Update of displacement (by DOF)
        mpDofUpdater->UpdateDofs(rDofSet, rDx);

        UpdateDerivatives(rModelPart, rDofSet, rA, rDx, rb);

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Performing the prediction of the solution
     * @details It predicts the solution for the current step
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

        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];

        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );
        const auto it_node_begin = rModelPart.Nodes().begin();

        const std::array<const Variable<double>*, 3> var_components = {&MOMENTUM_X, &MOMENTUM_Y, &HEIGHT};
        const std::array<const Variable<double>*, 3> accel_components = {&ACCELERATION_X, &ACCELERATION_Y, &VERTICAL_VELOCITY};

        IndexPartition<std::size_t>(num_nodes).for_each([&](std::size_t i){
            auto it_node = it_node_begin + i;

            for (std::size_t j = 0; j < 3; ++j)
            {
                if (!it_node->IsFixed(*var_components[j])) {
                    double& un0 = it_node->FastGetSolutionStepValue(*var_components[j]);
                    double un1 = it_node->FastGetSolutionStepValue(*var_components[j], 1);
                    double dot_un1 = it_node->FastGetSolutionStepValue(*accel_components[j], 1);
                    un0 = un1 + delta_time * dot_un1;
                }
            }

            UpdateFirstDerivative(it_node);
        });

        KRATOS_CATCH("TrapezoidalShallowScheme.Predict");
    }

    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @param rCurrentElement The element to compute
     * @param rLHS_Contribution The LHS matrix contribution
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        BaseType::CalculateSystemContributions(
            rCurrentElement,
            rLHS_Contribution,
            rRHS_Contribution,
            rEquationId,
            rCurrentProcessInfo);

        mRotationTool.Rotate(rLHS_Contribution,rRHS_Contribution,rCurrentElement.GetGeometry());
        mRotationTool.ApplySlipCondition(rLHS_Contribution,rRHS_Contribution,rCurrentElement.GetGeometry());
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rCurrentElement The element to compute
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        BaseType::CalculateRHSContribution(
            rCurrentElement,
            rRHS_Contribution,
            rEquationId,
            rCurrentProcessInfo);

        mRotationTool.Rotate(rRHS_Contribution,rCurrentElement.GetGeometry());
        mRotationTool.ApplySlipCondition(rRHS_Contribution,rCurrentElement.GetGeometry());
    }

    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @param rCurrentCondition The condition to compute
     * @param rLHS_Contribution The LHS matrix contribution
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        BaseType::CalculateSystemContributions(
            rCurrentCondition,
            rLHS_Contribution,
            rRHS_Contribution,
            rEquationId,
            rCurrentProcessInfo);

        mRotationTool.Rotate(rLHS_Contribution,rRHS_Contribution,rCurrentCondition.GetGeometry());
        mRotationTool.ApplySlipCondition(rLHS_Contribution,rRHS_Contribution,rCurrentCondition.GetGeometry());
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rCurrentCondition The condition to compute
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        BaseType::CalculateRHSContribution(
            rCurrentCondition,
            rRHS_Contribution,
            rEquationId,
            rCurrentProcessInfo);

        mRotationTool.Rotate(rRHS_Contribution,rCurrentCondition.GetGeometry());
        mRotationTool.ApplySlipCondition(rRHS_Contribution,rCurrentCondition.GetGeometry());
    }

    /*
     * @brief Free memory allocated by this class.
     */
    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    /**
     * @brief This method provides the defaults parameters
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"              : "trapezoidal_shallow_scheme",
            "theta"             : 0.5,
            "update_velocities" : true
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.AddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        mTheta = ThisParameters["theta"].GetDouble();
        mUpdateVelocities = ThisParameters["updateVelocities"].GetBool();
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
        return "TrapezoidalShallowScheme";
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

    std::size_t mTheta;

    bool mUpdateVelocities;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Updating first time derivative
     * @param itNode the node interator
     */
    void UpdateFirstDerivative(NodesArrayType::iterator itNode) override
    {
        array_1d<double, 3>& dot_un0 = itNode->FastGetSolutionStepValue(ACCELERATION);
        double& dot_hn0 = itNode->FastGetSolutionStepValue(VERTICAL_VELOCITY);
        noalias(dot_un0) = BaseType::mBDF[0] * itNode->FastGetSolutionStepValue(MOMENTUM);
        dot_hn0 = BaseType::mBDF[0] * itNode->FastGetSolutionStepValue(HEIGHT);
        for (std::size_t i_order = 1; i_order < BaseType::mOrder + 1; ++i_order)
        {
            noalias(dot_un0) += BaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(MOMENTUM, i_order);
            dot_hn0 += BaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(HEIGHT, i_order);
        }
    }

    /**
     * @brief Updating second time derivative
     * @param itNode the node interator
     */
    void UpdateSecondDerivative(NodesArrayType::iterator itNode) override {}

    /**
     * @brief Updating the velocities
     * @param rModelPart The model part to compute
     */
    void UpdateVelocities(ModelPart& rModelPart)
    {
        block_for_each(rModelPart.Nodes(), [&](NodeType& r_node){
            auto& vel = r_node.FastGetSolutionStepValue(VELOCITY);
            const auto& q = r_node.FastGetSolutionStepValue(MOMENTUM);
            const auto& h = r_node.FastGetSolutionStepValue(HEIGHT);
            vel = q / h;
        });
    }

    /**
     * @brief It adds the dynamic LHS contribution of the elements
     * @param rLHS_Contribution The dynamic contribution for the LHS
     * @param rD The damping matrix
     * @param rM The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // Adding mass contribution to the dynamic stiffness
        if (rM.size1() != 0) { // if M matrix declared
            noalias(rLHS_Contribution) += rM * BaseType::mBDF[0];
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the elements
     * @param rElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Element& rElement,
        LocalSystemVectorType& rRHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const auto& r_const_element = rElement;
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (rM.size1() != 0) {
            r_const_element.GetFirstDerivativesVector(BaseType::mVector.dotun0[this_thread], 0);
            noalias(rRHS_Contribution) -= prod(rM, BaseType::mVector.dotun0[this_thread]);
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the condition
     * @param rCondition The condition to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Condition& rCondition,
        LocalSystemVectorType& rRHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const auto& r_const_condition = rCondition;
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (rM.size1() != 0) {
            r_const_condition.GetFirstDerivativesVector(BaseType::mVector.dotun0[this_thread], 0);
            noalias(rRHS_Contribution) -= prod(rM, BaseType::mVector.dotun0[this_thread]);
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

}; // Class TrapezoidalShallowScheme

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // Namespace Kratos

#endif // KRATOS_TRAPEZOIDAL_SHALLOW_SCHEME_H_INCLUDED defined
