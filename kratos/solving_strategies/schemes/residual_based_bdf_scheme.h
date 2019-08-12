//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//


#if !defined(KRATOS_RESIDUAL_BASED_BDF_SCHEME )
#define  KRATOS_RESIDUAL_BASED_BDF_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/checks.h"
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
 * @class ResidualBasedBDFScheme
 * @ingroup KratosCore
 * @brief BDF integration scheme (for dynamic problems)
 * @details The \f$ n \f$ order Backward Differentiation Formula (BDF) method is a two step \f$ n \f$ order accurate method.
 * This scheme is designed to solve a system of the type:
 *\f[
 *   \mathbf{M} \frac{d^2(u_{n0})}{dt^2} + \mathbf{D} \frac{d(un0)}{dt} + \mathbf{K} u_{n0} = \mathbf{f}_{ext}
 * \f]
 *
 * If we call:
 *
 * - Second derivative:
 *      -# \f$ \ddot{u}_{ni} \f$ the second derivative at the step i
 * - First derivative:
 *      -# \f$ \dot{u}_{ni} \f$ the first derivative at the step i
 * - Third derivative:
 *      -# \f$ u_{ni} \f$ the variable  at the step i
 *
 * Then we assume:
 *  \f[ \frac{d^2(u_{n0})}{dt^2} \|t_{n0} = \sum_i c_i \dot{u}_{ni} \f]
 *  \f[ \frac{d(u_{n0})}{dt} \|t_{n0} = \sum_i c_i u_{n0} \f]
 * with for order 2 (BDF2):
 *  -# \f$ c_0 = \frac{1.5}{dt} \f$
 *  -# \f$ c_1 = \frac{-2.0}{dt} \f$
 *  -# \f$ c_2 = \frac{0.5}{dt} \f$
 *
 * The LHS and RHS can be defined as:
 *      \f[ RHS = \mathbf{f}_{ext} - \mathbf{M} \frac{d(\dot{u}_{n0})}{dt} - \mathbf{D} \frac{d(u_{n0})}{dt} - \mathbf{K} u_{n0} \f]
 * and
 *      \f[ LHS = \frac{d(-RHS)}{d(u_{n0})} = c_0^2 \mathbf{M} + c_0 \mathbf{D} + K \f]
 * @note This implies that elements are expected to be written in terms
 * of a variable with two time derivatives
 * <a href="https://mediatum.ub.tum.de/doc/1223319/80942.pdf">Main reference</a>
 * @todo Create a BibTeX file https://www.stack.nl/~dimitri/doxygen/manual/commands.html#cmdcite
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace,  class TDenseSpace>
class ResidualBasedBDFScheme
    : public ResidualBasedImplicitTimeScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedBDFScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                                  BaseType;

    typedef typename BaseType::Pointer                                 BaseTypePointer;

    typedef ResidualBasedImplicitTimeScheme<TSparseSpace,TDenseSpace> ImplicitBaseType;

    typedef typename ImplicitBaseType::TDataType                             TDataType;

    typedef typename ImplicitBaseType::DofsArrayType                     DofsArrayType;

    typedef typename Element::DofsVectorType                            DofsVectorType;

    typedef typename ImplicitBaseType::TSystemMatrixType             TSystemMatrixType;

    typedef typename ImplicitBaseType::TSystemVectorType             TSystemVectorType;

    typedef typename ImplicitBaseType::LocalSystemVectorType     LocalSystemVectorType;

    typedef typename ImplicitBaseType::LocalSystemMatrixType     LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                               NodesArrayType;

    /// Definition of epsilon
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor. The BDF method
     * @param Order The integration order
     * @todo The ideal would be to use directly the dof or the variable itself to identify the type of variable and is derivatives
     */
    explicit ResidualBasedBDFScheme(const std::size_t Order = 2)
        :ImplicitBaseType(),
         mOrder(Order)
    {
        // Allocate auxiliary memory
        const std::size_t num_threads = OpenMPUtils::GetNumThreads();

        mVector.dotun0.resize(num_threads);
        mVector.dot2un0.resize(num_threads);

        // Doing a minimal check
        KRATOS_ERROR_IF(mOrder < 1) << "ERROR:: Not possible to compute a BDF of order less than 1" << std::endl;

        // We resize the BDF coefficients
        if (mBDF.size() != (mOrder + 1))
            mBDF.resize(mOrder + 1);
    }

    /** Copy Constructor.
     */
    explicit ResidualBasedBDFScheme(ResidualBasedBDFScheme& rOther)
        :ImplicitBaseType(rOther)
        ,mOrder(rOther.mOrder)
        ,mBDF(rOther.mBDF)
        ,mVector(rOther.mVector)
    {
    }

    /**
     * Clone
     */
    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new ResidualBasedBDFScheme(*this) );
    }

    /** Destructor.
     */
    ~ResidualBasedBDFScheme
    () override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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

        KRATOS_ERROR << "Calling base BDF class" << std::endl;

        KRATOS_CATCH( "" );
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
        KRATOS_TRY;

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        ProcessInfo& r_previous_process_info = r_current_process_info.GetPreviousTimeStepInfo(1);

        ImplicitBaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        const double delta_time = r_current_process_info[DELTA_TIME];
        double previous_delta_time = r_previous_process_info[DELTA_TIME];

        KRATOS_ERROR_IF(delta_time < ZeroTolerance) << "Detected delta_time equal to zero or negative in the Solution Scheme DELTA_TIME: " << delta_time << ". PLEASE : check if the time step is created correctly for the current time step" << std::endl;
        KRATOS_WARNING_IF("ResidualBasedBDFScheme", previous_delta_time < ZeroTolerance) << "Detected previous_delta_time equal to zero or negative in the Solution Scheme DELTA_TIME: " << previous_delta_time << ". PLEASE : check if the time step is created correctly for the previous time step" << std::endl;
        previous_delta_time = std::abs(previous_delta_time) > ZeroTolerance ? previous_delta_time : delta_time;
        
        // Calculate the BDF coefficients
        const double rho = previous_delta_time / delta_time;
        double time_coeff = 0.0;
        for (std::size_t i_rho = 0; i_rho < mOrder; ++i_rho)
            time_coeff += delta_time * std::pow(rho, i_rho);
        time_coeff = 1.0/time_coeff;

        // We compute the BDF coefficients
        switch(mOrder) {
            case 1 : mBDF[0] =  time_coeff * rho; //coefficient for step n+1 (1Dt if Dt is constant)
                     mBDF[1] = -time_coeff * rho; //coefficient for step n (-1Dt if Dt is constant)
                     break;
            case 2 : mBDF[0] =  time_coeff * (std::pow(rho, 2) + 2.0 * rho); //coefficient for step n+1 (3/2Dt if Dt is constant)
                     mBDF[1] = -time_coeff * (std::pow(rho, 2) + 2.0 * rho + 1.0); //coefficient for step n (-4/2Dt if Dt is constant)
                     mBDF[2] =  time_coeff; //coefficient for step n-1 (1/2Dt if Dt is constant)
                     break;
            case 3 : mBDF[0] =  11.0/(6.0 * delta_time); //coefficient for step n+1 (11/6Dt if Dt is constant)
                     mBDF[1] = -18.0/(6.0 * delta_time);; //coefficient for step n (-18/6Dt if Dt is constant)
                     mBDF[2] =  9.0/(6.0 * delta_time);; //coefficient for step n-1 (9/6Dt if Dt is constant)
                     mBDF[3] = -2.0/(6.0 * delta_time);; //coefficient for step n-2 (2/6Dt if Dt is constant)
                     break;
            case 4 : mBDF[0] =  25.0/(12.0 * delta_time); //coefficient for step n+1 (25/12Dt if Dt is constant)
                     mBDF[1] = -48.0/(12.0 * delta_time); //coefficient for step n (-48/12Dt if Dt is constant)
                     mBDF[2] =  36.0/(12.0 * delta_time); //coefficient for step n-1 (36/12Dt if Dt is constant)
                     mBDF[3] = -16.0/(12.0 * delta_time); //coefficient for step n-2 (16/12Dt if Dt is constant)
                     mBDF[4] =  3.0/(12.0 * delta_time); //coefficient for step n-3 (3/12Dt if Dt is constant)
                     break;
            case 5 : mBDF[0] =  137.0/(60.0 * delta_time); //coefficient for step n+1 (137/60Dt if Dt is constant)
                     mBDF[1] = -300.0/(60.0 * delta_time); //coefficient for step n (-300/60Dt if Dt is constant)
                     mBDF[2] =  300.0/(60.0 * delta_time); //coefficient for step n-1 (300/60Dt if Dt is constant)
                     mBDF[3] = -200.0/(60.0 * delta_time); //coefficient for step n-2 (-200/60Dt if Dt is constant)
                     mBDF[4] =  75.0/(60.0 * delta_time); //coefficient for step n-3 (75/60Dt if Dt is constant)
                     mBDF[5] =  -12.0/(60.0 * delta_time); //coefficient for step n-4 (-12/60Dt if Dt is constant)
                     break;
            case 6 : mBDF[0] =  147.0/(60.0 * delta_time); //coefficient for step n+1 (147/60Dt if Dt is constant)
                     mBDF[1] = -360.0/(60.0 * delta_time); //coefficient for step n (-360/60Dt if Dt is constant)
                     mBDF[2] =  450.0/(60.0 * delta_time); //coefficient for step n-1 (450/60Dt if Dt is constant)
                     mBDF[3] = -400.0/(60.0 * delta_time); //coefficient for step n-2 (-400/60Dt if Dt is constant)
                     mBDF[4] =  225.0/(60.0 * delta_time); //coefficient for step n-3 (225/60Dt if Dt is constant)
                     mBDF[5] = -72.0/(60.0 * delta_time); //coefficient for step n-4 (-72/60Dt if Dt is constant)
                     mBDF[6] =  10.0/(60.0 * delta_time); //coefficient for step n-5 (10/60Dt if Dt is constant)
                     break;
            default : KRATOS_ERROR << "Methods with order > 6 are not zero-stable so they cannot be used" << std::endl;
        }

        const double tolerance = 1.0e-24;
        KRATOS_WARNING_IF("ResidualBasedBDFScheme", mOrder > 2 && std::abs(delta_time - previous_delta_time) > tolerance) << "For higher orders than 2 the time step is assumed to be constant. Sorry for the inconveniences" << std::endl;

        // Adding to the process info
        Vector bdf_vector(mOrder + 1);
        for (std::size_t i_order = 0; i_order < mOrder + 1; ++i_order)
            bdf_vector[i_order] = mBDF[i_order];
        r_current_process_info.SetValue(BDF_COEFFICIENTS, bdf_vector);

        KRATOS_CATCH( "" );
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed on the input provided.
     * @details Checks can be "expensive" as the function is designed to catch user's errors.
     * @param rModelPart The model of the problem to solve
     * @return Zero means  all ok
     */
    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        const int err = ImplicitBaseType::Check(rModelPart);
        if(err!=0) return err;

        // Check for minimum value of the buffer index
        // Verify buffer size
        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < mOrder + 1) << "Insufficient buffer size. Buffer size should be greater than" << mOrder + 1 << ". Current size is" << rModelPart.GetBufferSize() << std::endl;

        KRATOS_CATCH( "" );

        return 0;
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
        return "ResidualBasedBDFScheme";
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

    struct GeneralVectors
    {
        std::vector< Vector > dotun0; /// First derivative
        std::vector< Vector > dot2un0; /// Second derivative
    };

    const std::size_t mOrder; /// The integration order
    Vector mBDF;              /// The BDF coefficients
    GeneralVectors mVector;   /// The structure containing the  derivatives

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

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
        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );

        // Getting first node iterator
        const auto it_node_begin = rModelPart.Nodes().begin();

        #pragma omp parallel for
        for(int i = 0;  i< num_nodes; ++i) {
            auto it_node = it_node_begin + i;

            UpdateFirstDerivative(it_node);
            UpdateSecondDerivative(it_node);
        }
    }

    /**
     * @brief Updating first time derivative (velocity)
     * @param itNode the node interator
     */
    virtual inline void UpdateFirstDerivative(NodesArrayType::iterator itNode)
    {
        KRATOS_ERROR << "Calling base BDF class" << std::endl;
    }

    /**
     * @brief Updating second time derivative (acceleration)
     * @param itNode the node interator
     */
    virtual inline void UpdateSecondDerivative(NodesArrayType::iterator itNode)
    {
        KRATOS_ERROR << "Calling base BDF class" << std::endl;
    }

    /**
     * @brief It adds the dynamic LHS contribution of the elements
     * \f[ LHS = \frac{d(-RHS)}{d(u_{n0})} = c_0^2\mathbf{M} + c_0 \mathbf{D} + \mathbf{K} \f]
     * @param rLHS_Contribution The dynamic contribution for the LHS
     * @param rD The damping matrix
     * @param rM The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // Adding mass contribution to the dynamic stiffness
        if (rM.size1() != 0) { // if M matrix declared
            noalias(rLHS_Contribution) += rM * std::pow(mBDF[0], 2);
        }

        // Adding  damping contribution
        if (rD.size1() != 0) { // if D matrix declared
            noalias(rLHS_Contribution) += rD * mBDF[0];
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the objects
     * \f[ \mathbf{b} - \mathbf{M} a - \mathbf{D} v \f]
     * @param rObject The object to compute
     * @param rRHS_Contribution The dynamic contribution for the RHS
     * @param rD The damping matrix
     * @param rM The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    template <typename TObjectType>
    void TemplateAddDynamicsToRHS(
        TObjectType rObject,
        LocalSystemVectorType& rRHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (rM.size1() != 0) {
            rObject->GetSecondDerivativesVector(mVector.dot2un0[this_thread], 0);
            noalias(rRHS_Contribution) -= prod(rM, mVector.dot2un0[this_thread]);
        }

        // Adding damping contribution
        if (rD.size1() != 0) {
            rObject->GetFirstDerivativesVector(mVector.dotun0[this_thread], 0);
            noalias(rRHS_Contribution) -= prod(rD, mVector.dotun0[this_thread]);
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the elements
     * \f[ \mathbf{b} - \mathbf{M} a - \mathbf{D} v \f]
     * @param pElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Element::Pointer pElement,
        LocalSystemVectorType& rRHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        TemplateAddDynamicsToRHS<Element::Pointer>(pElement, rRHS_Contribution, rD, rM, rCurrentProcessInfo);
    }

    /**
     * @brief It adds the dynamic RHS contribution of the condition
     *  \f[ RHS = f_{ext} - \ddot{u}_{n0} \mathbf{M} + \dot{u}_{n0} \mathbf{D} + u_{n0} \mathbf{K} \f]
     * @param pCondition The condition to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Condition::Pointer pCondition,
        LocalSystemVectorType& rRHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        TemplateAddDynamicsToRHS<Condition::Pointer>(pCondition, rRHS_Contribution, rD, rM, rCurrentProcessInfo);
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
}; /* Class ResidualBasedBDFScheme */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BDF_SCHEME defined */
