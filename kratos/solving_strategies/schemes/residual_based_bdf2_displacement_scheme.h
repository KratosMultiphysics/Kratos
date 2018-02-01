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


#if !defined(KRATOS_RESIDUAL_BASED_BDF2_DISPLACEMENT_SCHEME )
#define  KRATOS_RESIDUAL_BASED_BDF2_DISPLACEMENT_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/residual_implicit_time_scheme.h"
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

/** @brief BDF2 integration scheme (for dynamic problems, displacement based)
 * @details The second order Backward Differentiation Formula (BDF2) method is a two step second order accurate method. 
 * This scheme is designed to solve a system of the type:
 *\f[
 *   \mathbf{M} \frac{d^2(u_{n0})}{dt^2} + \mathbf{D} \frac{d(un0)}{dt} + \mathbf{K} u_{n0} = \mathbf{f}_{ext}
 * \f]
 * 
 * If we call:
 * 
 * - Acelerations: 
 *      -# \f$ a_{n0} \f$ the acceleration at the current step
 *      -# \f$ a_{n1} \f$ the acceleration one step in the past
 *      -# \f$ a_{n2} \f$ the acceleration two steps in the past
 * - Velocities:
 *     -# \f$ v_{n0} \f$ the velocity at the current step
 *     -# \f$ v_{n1} \f$ the velocity one step in the past
 *     -# \f$ v_{n2} \f$ the velocity two steps in the past
 * - Displacements:
 *     -# \f$ u_{n0} \f$ the displacement at the current step
 *     -# \f$ u_{n1} \f$ the displacement one step in the past
 *     -# \f$ u_{n2} \f$ the displacement two steps in the past
 * 
 * Then we assume:
 *  \f[ \frac{d(vn0)}{dt} \|t_{n0} = c_0 v_{n0} + c_1 v_{n1} + c_2 v_{n2} \f]
 *  \f[ \frac{d(un0)}{dt} \|t_{n0} = c_0 u_{n0} + c_1 u_{n1} + c_2 u_{n2} \f]
 * with:
 *  -# \f$ c_0 = \frac{1.5}{dt} \f$ 
 *  -# \f$ c_1 = \frac{-2.0}{dt} \f$ 
 *  -# \f$ c_2 = \frac{0.5}{dt} \f$  
 * 
 * The LHS and RHS can be defined as:
 *      \f[ RHS = \mathbf{f}_{ext} - \mathbf{M} \frac{d(v_{n0})}{dt} - \mathbf{D} \frac{d(u_{n0})}{dt} - \mathbf{K} u_{n0} \f]
 * and 
 *      \f[ LHS = \frac{d(-RHS)}{d(u_{n0})} = c_0^2 \mathbf{M} + c_0 \mathbf{D} + K \f]
 * @note This implies that elements are expected to be written in terms 
 * of unknown DISPLACEMENTS
 * <a href="https://mediatum.ub.tum.de/doc/1223319/80942.pdf">Main reference</a> 
 * @todo Create a BibTeX file https://www.stack.nl/~dimitri/doxygen/manual/commands.html#cmdcite
 */
template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedBDF2DisplacementScheme
    : public ResidualBasedImplicitTimeScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedBDF2DisplacementScheme );

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
     * Constructor.
     * The DBF2 method
     */
    ResidualBasedBDF2DisplacementScheme()
        :ImplicitBaseType()
    {
    }

    /** Copy Constructor.
     */
    ResidualBasedBDF2DisplacementScheme(ResidualBasedBDF2DisplacementScheme& rOther)
        :ImplicitBaseType(rOther)
        ,mBDF2(rOther.mBDF2)
    {
    }

    /**
     * Clone
     */
    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new ResidualBasedBDF2DisplacementScheme(*this) );
    }

    /** Destructor.
     */
    ~ResidualBasedBDF2DisplacementScheme
    () override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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

        // Update of displacement (by DOF)
        const int num_dof = static_cast<int>(rDofSet.size());

        #pragma omp parallel for
        for(int i = 0;  i < num_dof; ++i) {
            auto it_dof = rDofSet.begin() + i;

            if (it_dof->IsFree())
                it_dof->GetSolutionStepValue() += TSparseSpace::GetValue(Dx,it_dof->EquationId());
        }

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>(rModelPart.Nodes().size());

        #pragma omp parallel for
        for(int i = 0;  i < num_nodes; ++i) {
            auto it_node = rModelPart.Nodes().begin() + i;
                        
            const array_1d<double, 3>& un0 = it_node->FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3>& un1 = it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double, 3>& un2 = it_node->FastGetSolutionStepValue(DISPLACEMENT, 2);
            
            array_1d<double, 3>& vn0 = it_node->FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& vn1 = it_node->FastGetSolutionStepValue(VELOCITY, 1);
            const array_1d<double, 3>& vn2 = it_node->FastGetSolutionStepValue(VELOCITY, 2);

            array_1d<double, 3>& an0 = it_node->FastGetSolutionStepValue(ACCELERATION);

            UpdateVelocity(vn0, un0, un1, un2);
            UpdateAcceleration(an0, vn0, vn1, vn2);
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

        ProcessInfo& current_process_info = rModelPart.GetProcessInfo();
        
        const double delta_time = current_process_info[DELTA_TIME];

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );

        #pragma omp parallel for
        for(int i = 0;  i< num_nodes; ++i) {
            auto it_node = rModelPart.Nodes().begin() + i;

            //Predicting: NewDisplacement = previous_displacement + previous_velocity * delta_time;
            //ATTENTION::: the prediction is performed only on free nodes

            const array_1d<double, 3>& vn2 = it_node->FastGetSolutionStepValue(VELOCITY,     2);
            const array_1d<double, 3>& un2 = it_node->FastGetSolutionStepValue(DISPLACEMENT, 2);
            const array_1d<double, 3>& an1 = it_node->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double, 3>& vn1 = it_node->FastGetSolutionStepValue(VELOCITY,     1);
            const array_1d<double, 3>& un1 = it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double, 3>& an0 = it_node->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3>& vn0 = it_node->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3>& un0 = it_node->FastGetSolutionStepValue(DISPLACEMENT);

            if (it_node -> IsFixed(ACCELERATION_X)) {
                vn0[0] = (an0[0] - mBDF2.c1 * vn1[0] - mBDF2.c2 * vn2[0])/mBDF2.c0;
                un0[0] = (vn0[0] - mBDF2.c1 * un1[0] - mBDF2.c2 * un2[0])/mBDF2.c0;
            } else if (it_node -> IsFixed(VELOCITY_X)) {
                un0[0] = (vn1[0] - mBDF2.c1 * un1[0] - mBDF2.c2 * un2[0])/mBDF2.c0;
            } else if (it_node -> IsFixed(DISPLACEMENT_X) == false) {
                un0[0] = un1[0] + delta_time * vn1[0] + 0.5 * std::pow(delta_time, 2) * an1[0];
            }

            if (it_node -> IsFixed(ACCELERATION_Y)) {
                vn0[1] = (an0[1] - mBDF2.c1 * vn1[1] - mBDF2.c2 * vn2[1])/mBDF2.c0;
                un0[1] = (vn0[1] - mBDF2.c1 * un1[1] - mBDF2.c2 * un2[1])/mBDF2.c0;
            } else if (it_node -> IsFixed(VELOCITY_Y)) {
                un0[1] = (vn1[1] - mBDF2.c1 * un1[1] - mBDF2.c2 * un2[1])/mBDF2.c0;
            } else if (it_node -> IsFixed(DISPLACEMENT_Y) == false) {
                un0[1] = un1[1] + delta_time * vn1[1] + 0.5 * std::pow(delta_time, 2) * an1[1];
            }

            // For 3D cases
            if (it_node -> HasDofFor(DISPLACEMENT_Z)) {
                if (it_node -> IsFixed(ACCELERATION_Z)) {
                    vn0[2] = (an0[2] - mBDF2.c1 * vn1[2] - mBDF2.c2 * vn2[2])/mBDF2.c0;
                    un0[2] = (vn0[2] - mBDF2.c1 * un1[2] - mBDF2.c2 * un2[2])/mBDF2.c0;
                } else if (it_node -> IsFixed(VELOCITY_Z)) {
                    un0[2] = (vn1[2] - mBDF2.c1 * un1[2] - mBDF2.c2 * un2[2])/mBDF2.c0;
                } else if (it_node -> IsFixed(DISPLACEMENT_Z) == false) {
                    un0[2] = un1[2] + delta_time * vn1[2] + 0.5 * std::pow(delta_time, 2) * an1[2];
                }
            }

            // Updating time derivatives ::: Please note that displacements and its time derivatives can not be consistently fixed separately
            UpdateVelocity(vn0, un0, un1, un2);
            UpdateAcceleration(an0, vn0, vn1, vn2);
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

        ProcessInfo& current_process_info = rModelPart.GetProcessInfo();

        ImplicitBaseType::InitializeSolutionStep(rModelPart, A, Dx, b);

        const double delta_time = current_process_info[DELTA_TIME];
        const double previous_delta_time = current_process_info.GetPreviousTimeStepInfo(1)[DELTA_TIME];
        
        // Calculate the BDF coefficients
        const double rho = previous_delta_time / delta_time;
        const double time_coeff = 1.0 / (delta_time * std::pow(rho, 2) + delta_time * rho);
        
        mBDF2.c0 =  time_coeff * (std::pow(rho, 2) + 2.0 * rho); //coefficient for step n+1 (3/2Dt if Dt is constant)
        mBDF2.c1 = -time_coeff * (std::pow(rho, 2) + 2.0 * rho + 1.0); //coefficient for step n (-4/2Dt if Dt is constant)
        mBDF2.c2 =  time_coeff; //coefficient for step n-1 (1/2Dt if Dt is constant)
        
        // Adding to the process info
        Vector bdf_vector(3);
        bdf_vector[0] = mBDF2.c0;
        bdf_vector[1] = mBDF2.c1;
        bdf_vector[2] = mBDF2.c2;
        current_process_info(BDF_COEFFICIENTS) = bdf_vector;
        
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
        if(err!=0) return err;

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
        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 3) << "Insufficient buffer size. Buffer size should be greater than 3. Current size is" << rModelPart.GetBufferSize() << std::endl;

        KRATOS_CATCH( "" );
        
        return 0;
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

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    struct BDF2Method
    {
        double c0, c1, c2;
    };

    BDF2Method mBDF2; // The BDF2 coefficients

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    
    /**
     * @brief Updating first time Derivative
     * @param vn0 the velocity at the current step
     * @param un0 the displacement at the current step
     * @param un1 the displacement one step in the past
     * @param un2 the displacement two steps in the past
     */
    
    inline void UpdateVelocity(
        array_1d<double, 3> & vn0,
        const array_1d<double, 3>& un0,
        const array_1d<double, 3>& un1,
        const array_1d<double, 3>& un2
        )
    {
        noalias(vn0) = (mBDF2.c0 * un0 + mBDF2.c1 * un1  + mBDF2.c2 * un2);
    }

    /**
     * @brief Updating second time Derivative
     * @param an0 the acceleration at the current step
     * @param vn0 the velocity at the current step
     * @param vn1 the velocity one step in the past
     * @param vn2 the velocity two steps in the past
     */

    inline void UpdateAcceleration(
        array_1d<double, 3> & an0,
        const array_1d<double, 3>& vn0,
        const array_1d<double, 3>& vn1,
        const array_1d<double, 3>& vn2
    )
    {
        noalias(an0) = (mBDF2.c0 * vn0 + mBDF2.c1 * vn1  + mBDF2.c2 * vn2);
    }

    /**
     * @brief It adds the dynamic LHS contribution of the elements LHS = d(-RHS)/d(un0) = c0*c0*M + c0*D + K
     * @param LHS_Contribution The dynamic contribution for the LHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param CurrentProcessInfo The current process info instance
     */

    void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo
        ) override
    {
        // Adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) { // if M matrix declared
            noalias(LHS_Contribution) += M * std::pow(mBDF2.c0, 2);
        }

        // Adding  damping contribution
        if (D.size1() != 0) { // if D matrix declared
            noalias(LHS_Contribution) += D * mBDF2.c0;
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the elements b - M*a - D*v
     * @param pElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param CurrentProcessInfo The current process info instance
     */

    void AddDynamicsToRHS(
        Element::Pointer pElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo
        ) override
    {
        // Adding inertia contribution
        if (M.size1() != 0) {
            Vector an0;
            pElement->GetSecondDerivativesVector(an0, 0);
            noalias(RHS_Contribution) -= prod(M, an0);
        }

        // Adding damping contribution
        if (D.size1() != 0) {
            Vector vn0;
            pElement->GetFirstDerivativesVector(vn0, 0);
            noalias(RHS_Contribution) -= prod(D, vn0);
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the condition RHS = fext - M*an0 - D*vn0 - K*dn0
     * @param pCondition The condition to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param CurrentProcessInfo The current process info instance
     */

    void AddDynamicsToRHS(
        Condition::Pointer pCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo
        ) override
    {
        // Adding inertia contribution
        if (M.size1() != 0) {
            Vector an0;
            pCondition->GetSecondDerivativesVector(an0, 0);
            noalias(RHS_Contribution) -= prod(M, an0);
        }

        // Adding damping contribution
        // Damping contribution
        if (D.size1() != 0) {
            Vector vn0;
            pCondition->GetFirstDerivativesVector(vn0, 0);
            noalias(RHS_Contribution) -= prod(D, vn0);
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
}; /* Class ResidualBasedBDF2DisplacementScheme */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BDF2_DISPLACEMENT_SCHEME defined */
