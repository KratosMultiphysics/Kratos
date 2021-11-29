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


#ifndef KRATOS_RESIDUAL_BASED_ADAMS_MOULTON_SCHEME_H_INCLUDED
#define KRATOS_RESIDUAL_BASED_ADAMS_MOULTON_SCHEME_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/scheme.h"
#include "utilities/entities_utilities.h"

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
 * @class ResidualBasedAdamsMoultonScheme
 * @ingroup ShallowWaterApplication
 * @brief Predictor-corrector semi imlicit scheme for the Boussinesq element.
 * @details The element should assemble the RHS according to the scheme.
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 * @author Miguel Maso Sotomayor
 */
template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedAdamsMoultonScheme
    : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;

    typedef typename BaseType::DofsArrayType                 DofsArrayType;

    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType         TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ModelPart::NodeType                                   NodeType;

    ///@}
    ///@name Pointer definition
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedAdamsMoultonScheme );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     */
    explicit ResidualBasedAdamsMoultonScheme() : BaseType()
    {
        // Allocate auxiliary memory
        const std::size_t num_threads = ParallelUtilities::GetNumThreads();
        mM.resize(num_threads);
        mU0.resize(num_threads);
        mU1.resize(num_threads);
        mDU.resize(num_threads);
    }

    /**
     * @brief Constructor.
     * @param ThisParameters The configuration parameters
     */
    explicit ResidualBasedAdamsMoultonScheme(Parameters ThisParameters)
        : BaseType(ThisParameters)
    {
        // Allocate auxiliary memory
        const std::size_t num_threads = ParallelUtilities::GetNumThreads();
        mM.resize(num_threads);
        mU0.resize(num_threads);
        mU1.resize(num_threads);
        mDU.resize(num_threads);
    }

    /** 
     * @brief Copy Constructor.
     */
    explicit ResidualBasedAdamsMoultonScheme(ResidualBasedAdamsMoultonScheme& rOther)
        : BaseType(rOther)
        , mM(rOther.mM)
        , mU0(rOther.mU0)
        , mU1(rOther.mU1)
        , mDU(rOther.mDU)
    {
    }

    /**
     * @brief Clone
     */
    typename BaseType::Pointer Clone() override
    {
        return Kratos::make_shared<ResidualBasedAdamsMoultonScheme>(*this);
    }

    /** 
     * @brief Destructor.
     */
    ~ResidualBasedAdamsMoultonScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Perform the prediction using the explicit Adams-Bashforth scheme
     * @param rModelPart The model of the problem to solve
     * @param rDofSet set of all primary variables
     * @param rA LHS matrix
     * @param rDx Incremental prediction of primary variables
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
        KRATOS_ERROR_IF(delta_time < 1.0e-24) << "ERROR:: Detected delta_time near to zero" << std::endl;

        PredictDerivatives(rModelPart, rDofSet, rA, rDx, rb);

        // const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        // block_for_each(rModelPart.Elements(), [&](Element& rElement){
        //     rElement.AddExplicitContribution(r_process_info);
        // });

        KRATOS_CATCH( "" );
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
        KRATOS_TRY;

        // Update of displacement (by DOF)
        mpDofUpdater->UpdateDofs(rDofSet, rDx);

        UpdateDerivatives(rModelPart, rDofSet, rA, rDx, rb);

        KRATOS_CATCH( "" );
    }

    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @details It "asks" the matrix needed to the element and performs the operations needed to introduce the selected time integration scheme. This function calculates at the same time the contribution to the LHS and to the RHS of the system
     * @param rCurrentElement The element to compute
     * @param rLHSContribution The LHS matrix contribution
     * @param rRHSContribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHSContribution,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        TCalculateSystemContributions(rCurrentElement, rLHSContribution, rRHSContribution, rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedAdamsMoultonScheme.CalculateSystemContributions");
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
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        TCalculateRHSContribution(rCurrentElement, rRHSContribution, rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedAdamsMoultonScheme.CalculateRHSContribution");
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
        LocalSystemMatrixType& rLHSContribution,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        TCalculateSystemContributions(rCurrentCondition, rLHSContribution, rRHSContribution, rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedAdamsMoultonScheme.CalculateSystemContributions");
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
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        TCalculateRHSContribution(rCurrentCondition, rRHSContribution, rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedAdamsMoultonScheme.CalculateRHSContribution");
    }

    /**
     * @brief Free memory allocated by this class.
     */
    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "residual_based_adams_moulton_scheme"
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedAdamsMoultonScheme";
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    std::vector< Matrix > mM; /// First derivative matrix (usually mass matrix)
    std::vector< Vector > mU0; /// Values vector at the current time step
    std::vector< Vector > mU1; /// Values vector at the previous time step
    std::vector< Vector > mDU; /// Increment of the values vector

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Performing the prediction of the derivatives
     * @param rModelPart The model of the problem to solve
     * @param rDofSet Set of all primary variables
     * @param rA LHS matrix
     * @param rDx incremental update of primary variables
     * @param rb RHS Vector
     */
    inline void PredictDerivatives(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
        const double factor = 0.5 / delta_time;

        block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
            PredictDerivative(rNode, VELOCITY, ACCELERATION, factor);
            PredictDerivative(rNode, FREE_SURFACE_ELEVATION, VERTICAL_VELOCITY, factor);
        });
    }

    /**
     * @brief Performing the prediction of the derivative
     * @param rNode The node
     * @param rPrimitiveVariable The primitive variable
     * @param rDerivativeVariable The variable to predict
     * @param Coefficient The time integration coefficient
     */
    template<class TDataType>
    inline void PredictDerivative(
        NodeType& rNode,
        const Variable<TDataType>& rPrimitiveVariable,
        const Variable<TDataType>& rDerivativeVariable,
        const double Coefficient)
    {
        const TDataType& f1 = rNode.FastGetSolutionStepValue(rPrimitiveVariable, 1);
        const TDataType& f2 = rNode.FastGetSolutionStepValue(rPrimitiveVariable, 2);
        const TDataType& f3 = rNode.FastGetSolutionStepValue(rPrimitiveVariable, 3);

        TDataType& d1 = rNode.FastGetSolutionStepValue(rDerivativeVariable, 1);
        TDataType& d2 = rNode.FastGetSolutionStepValue(rDerivativeVariable, 2);
        TDataType& d3 = rNode.FastGetSolutionStepValue(rDerivativeVariable, 3);

        d1 = Coefficient * (3*f1 - 4*f2   + f3);
        d2 = Coefficient * (  f1          - f3);
        d3 = Coefficient * ( -f1 + 4*f2  -3*f3);
    }

    /**
     * @brief Performing the update of the derivatives
     * @param rModelPart The model of the problem to solve
     * @param rDofSet Set of all primary variables
     * @param rA LHS matrix
     * @param rDx incremental update of primary variables
     * @param rb RHS Vector
     */
    inline void UpdateDerivatives(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
        const double factor = 1.0 / (6.0 * delta_time);

        block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
            UpdateDerivative(rNode, VELOCITY, ACCELERATION, factor);
            UpdateDerivative(rNode, FREE_SURFACE_ELEVATION, VERTICAL_VELOCITY, factor);
        });
    }

    /**
     * @brief Performing the prediction of the derivative
     * @param rNode The node
     * @param rPrimitiveVariable The primitive variable
     * @param rDerivativeVariable The variable to predict
     * @param Coefficient The time integration coefficient
     */
    template<class TDataType>
    inline void UpdateDerivative(
        NodeType& rNode,
        const Variable<TDataType>& rPrimitiveVariable,
        const Variable<TDataType>& rDerivativeVariable,
        const double Coefficient)
    {
        const TDataType& f0 = rNode.FastGetSolutionStepValue(rPrimitiveVariable, 0);
        const TDataType& f1 = rNode.FastGetSolutionStepValue(rPrimitiveVariable, 1);
        const TDataType& f2 = rNode.FastGetSolutionStepValue(rPrimitiveVariable, 2);
        const TDataType& f3 = rNode.FastGetSolutionStepValue(rPrimitiveVariable, 3);

        TDataType& d0 = rNode.FastGetSolutionStepValue(rDerivativeVariable, 0);
        TDataType& d1 = rNode.FastGetSolutionStepValue(rDerivativeVariable, 1);
        TDataType& d2 = rNode.FastGetSolutionStepValue(rDerivativeVariable, 2);
        TDataType& d3 = rNode.FastGetSolutionStepValue(rDerivativeVariable, 3);

        d0 = Coefficient * (11*f0 - 18*f1  + 9*f2  - 2*f3);
        d1 = Coefficient * ( 2*f0  - 3*f1  - 6*f2    + f3);
        d2 = Coefficient * (  -f0  + 6*f1  - 3*f2  - 2*f3);
        d3 = Coefficient * ( 2*f0  - 9*f1 + 18*f2 - 11*f3);
    }

    /**
     * @brief It adds the dynamic LHS contribution of the elements LHS = 24/dt*M
     * @param rLHSContribution The dynamic contribution for the LHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& rLHSContribution,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        const double delta_time = rCurrentProcessInfo[DELTA_TIME];

        // Adding inertia contribution
        if (rM.size1() != 0) {
            rLHSContribution = 24.0 / delta_time * rM;
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the elements b - M*a - D*v
     * @param rCurrentElement The element to compute
     * @param rRHSContribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Element& rCurrentElement,
        LocalSystemVectorType& rRHSContribution,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();
        const double delta_time = rCurrentProcessInfo[DELTA_TIME];

        // Adding inertia contribution
        if (rM.size1() != 0) {
            rCurrentElement.GetValuesVector(mU0[this_thread], 0);
            rCurrentElement.GetValuesVector(mU1[this_thread], 1);
            mDU[this_thread] = mU0[this_thread] - mU1[this_thread];
            noalias(rRHSContribution) -= 24 / delta_time * prod(rM, mDU[this_thread]);
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the condition RHS = fext - M*an0 - D*vn0 - K*dn0
     * @param rCurrentCondition The condition to compute
     * @param rRHSContribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Condition& rCurrentCondition,
        LocalSystemVectorType& rRHSContribution,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();
        const double delta_time = rCurrentProcessInfo[DELTA_TIME];

        // Adding inertia contribution
        if (rM.size1() != 0) {
            rCurrentCondition.GetValuesVector(mU0[this_thread], 0);
            rCurrentCondition.GetValuesVector(mU1[this_thread], 1);
            mDU[this_thread] = mU0[this_thread] - mU1[this_thread];
            noalias(rRHSContribution) -= 24 / delta_time * prod(rM, mDU[this_thread]);
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

    /// Utility class to perform the update after solving the system, will be different in MPI runs.
    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @param rObject The object to compute
     * @param rLHSContribution The LHS matrix contribution
     * @param rRHSContribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    template <class TObjectType>
    void TCalculateSystemContributions(
        TObjectType& rObject,
        LocalSystemMatrixType& rLHSContribution,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY;

        const auto this_thread = OpenMPUtils::ThisThread();

        rObject.CalculateRightHandSide(rRHSContribution, rCurrentProcessInfo);

        rObject.CalculateMassMatrix(mM[this_thread], rCurrentProcessInfo);

        rObject.EquationIdVector(EquationId, rCurrentProcessInfo);

        AddDynamicsToLHS(rLHSContribution, mM[this_thread], rCurrentProcessInfo);

        AddDynamicsToRHS(rObject, rRHSContribution, mM[this_thread], rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedAdamsMoultonScheme.TCalculateSystemContributions");
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rObject The object to compute
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    template <class TObjectType>
    void TCalculateRHSContribution(
        TObjectType& rObject,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY;

        const auto this_thread = OpenMPUtils::ThisThread();

        rObject.CalculateRightHandSide(rRHSContribution, rCurrentProcessInfo);

        rObject.CalculateMassMatrix(mM[this_thread], rCurrentProcessInfo);

        rObject.EquationIdVector(rEquationId,rCurrentProcessInfo);

        AddDynamicsToRHS(rObject, rRHSContribution, mM[this_thread], rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedAdamsMoultonScheme.TCalculateRHSContribution");
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    ///@}

}; // Class ResidualBasedAdamsMoultonScheme

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_RESIDUAL_BASED_ADAMS_MOULTON_SCHEME_H_INCLUDED defined
