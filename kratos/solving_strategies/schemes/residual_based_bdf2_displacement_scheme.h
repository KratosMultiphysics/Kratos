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
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"

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
 * This scheme is designed to solve a system of the type M*d2(un0)/dt2 + D*d(un0)/dt + K*un0 = fext
 * If we call:
 * 
 *      an0 the acceleration at the current step
 *      an1 the acceleration one step in the past
 *      an2 the acceleration two steps in the past
 * 
 *      vn0 the velocity at the current step
 *      vn1 the velocity one step in the past
 *      vn2 the velocity two steps in the past
 * 
 *      un0 the displacemnt at the current step
 *      un1 the displacemnt one step in the past
 *      un2 the displacemnt two steps in the past
 * 
 * Then we assume:
 *      d(vn0)/dt |tn0 = c0*vn0 + c1*vn1 + c2*vn2 
 *      d(un0)/dt |tn0 = c0*un0 + c1*un1 + c2*un2
 * with:
 *      c0 = 1.5/dt
 *      c1 = -2.0/dt
 *      c2 = 0.5/dt 
 * 
 * The LHS and RHS can be defined as:
 *      RHS = fext - M*d(vn0)/dt - D*d(un0)/dt - K*dn0
 * and 
 *      LHS = d(-RHS)/d(un0) = c0*c0*M + c0*D + K
 * NOTE: this implies that elements are expected to be written in terms 
 * of unknown DISPLACEMENTS
 */
template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedBDF2DisplacementScheme
    : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedBDF2DisplacementScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;

    typedef typename BaseType::TDataType                         TDataType;

    typedef typename BaseType::DofsArrayType                 DofsArrayType;

    typedef typename Element::DofsVectorType                DofsVectorType;

    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType         TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                   NodesArrayType;

    typedef ModelPart::ElementsContainerType             ElementsArrayType;

    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;

    typedef typename BaseType::Pointer                     BaseTypePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     * The DBF2 method
     */
    ResidualBasedBDF2DisplacementScheme()
        :BaseType()
    {
        // Allocate auxiliary memory
        const unsigned int num_threads = OpenMPUtils::GetNumThreads();

        mMatrix.M.resize(num_threads);
        mMatrix.D.resize(num_threads);

        mVector.v.resize(num_threads);
        mVector.a.resize(num_threads);
    }

    /** Copy Constructor.
     */
    ResidualBasedBDF2DisplacementScheme(ResidualBasedBDF2DisplacementScheme& rOther)
        :BaseType(rOther)
        ,mBDF2(rOther.mBDF2)
        ,mMatrix(rOther.mMatrix)
        ,mVector(rOther.mVector)
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
     * Performing the update of the solution
     * Incremental update within newton iteration. It updates the state variables at the end of the time step u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
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

        const unsigned int num_threads = OpenMPUtils::GetNumThreads();

        // Update of displacement (by DOF)
        const int num_dof = static_cast<int>(rDofSet.size());
        OpenMPUtils::PartitionVector dof_partition;
        OpenMPUtils::DivideInPartitions(num_dof, num_threads, dof_partition);

        #pragma omp parallel for
        for(int i = 0;  i < num_dof; ++i)
        {
            auto it_dof = rDofSet.begin() + i;

            if (it_dof->IsFree())
            {
                it_dof->GetSolutionStepValue() += TSparseSpace::GetValue(Dx,it_dof->EquationId());
            }
        }

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>(rModelPart.Nodes().size());
        OpenMPUtils::PartitionVector node_partition;
        OpenMPUtils::DivideInPartitions(num_nodes, num_threads, node_partition);

        #pragma omp parallel for
        for(int i = 0;  i < num_nodes; ++i)
        {
            auto it_node = rModelPart.Nodes().begin() + i;
                        
            const array_1d<double, 3 > & un0 = it_node->FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3 > & un1 = it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double, 3 > & un2 = it_node->FastGetSolutionStepValue(DISPLACEMENT, 2);
            
            array_1d<double, 3 > & vn0 = it_node->FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3 > & vn1 = it_node->FastGetSolutionStepValue(VELOCITY, 1);
            const array_1d<double, 3 > & vn2 = it_node->FastGetSolutionStepValue(VELOCITY, 2);

            array_1d<double, 3 > & an0 = it_node->FastGetSolutionStepValue(ACCELERATION);

            UpdateVelocity(vn0, un0, un1, un2);
            UpdateAcceleration(an0, vn0, vn1, vn2);
        }

        KRATOS_CATCH( "" );
    }

    /**
     * Performing the prediction of the solution
     * It predicts the solution for the current step x = xold + vold * Dt
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
        const unsigned int num_threads = OpenMPUtils::GetNumThreads();
        
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );
        OpenMPUtils::PartitionVector node_partition;
        OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), num_threads, node_partition);

        #pragma omp parallel for
        for(int i = 0;  i< num_nodes; ++i)
        {
            auto it_node = rModelPart.Nodes().begin() + i;

            //Predicting: NewDisplacement = previous_displacement + previous_velocity * delta_time;
            //ATTENTION::: the prediction is performed only on free nodes

            const array_1d<double, 3 > & vn2 = it_node->FastGetSolutionStepValue(VELOCITY,     2);
            const array_1d<double, 3 > & un2 = it_node->FastGetSolutionStepValue(DISPLACEMENT, 2);
            const array_1d<double, 3 > & an1 = it_node->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double, 3 > & vn1 = it_node->FastGetSolutionStepValue(VELOCITY,     1);
            const array_1d<double, 3 > & un1 = it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double, 3 > & an0 = it_node->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3 > & vn0 = it_node->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & un0 = it_node->FastGetSolutionStepValue(DISPLACEMENT);

            if (it_node -> IsFixed(ACCELERATION_X))
            {
                vn0[0] = (an0[0] - mBDF2.c1 * vn1[0] - mBDF2.c2 * vn2[0])/mBDF2.c0;
                un0[0] = (vn0[0] - mBDF2.c1 * un1[0] - mBDF2.c2 * un2[0])/mBDF2.c0;
            }
            else if (it_node -> IsFixed(VELOCITY_X))
            {
                un0[0] = (vn1[0] - mBDF2.c1 * un1[0] - mBDF2.c2 * un2[0])/mBDF2.c0;
            }
            else if (it_node -> IsFixed(DISPLACEMENT_X) == false)
            {
                un0[0] = un1[0] + delta_time * vn1[0] + 0.5 * std::pow(delta_time, 2) * an1[0];
            }

            if (it_node -> IsFixed(ACCELERATION_Y))
            {
                vn0[1] = (an0[1] - mBDF2.c1 * vn1[1] - mBDF2.c2 * vn2[1])/mBDF2.c0;
                un0[1] = (vn0[1] - mBDF2.c1 * un1[1] - mBDF2.c2 * un2[1])/mBDF2.c0;
            }
            else if (it_node -> IsFixed(VELOCITY_Y))
            {
                un0[1] = (vn1[1] - mBDF2.c1 * un1[1] - mBDF2.c2 * un2[1])/mBDF2.c0;
            }
            else if (it_node -> IsFixed(DISPLACEMENT_Y) == false)
            {
                un0[1] = un1[1] + delta_time * vn1[1] + 0.5 * std::pow(delta_time, 2) * an1[1];
            }

            // For 3D cases
            if (it_node -> HasDofFor(DISPLACEMENT_Z))
            {
                if (it_node -> IsFixed(ACCELERATION_Z))
                {
                    vn0[2] = (an0[2] - mBDF2.c1 * vn1[2] - mBDF2.c2 * vn2[2])/mBDF2.c0;
                    un0[2] = (vn0[2] - mBDF2.c1 * un1[2] - mBDF2.c2 * un2[2])/mBDF2.c0;
                }
                else if (it_node -> IsFixed(VELOCITY_Z))
                {
                    un0[2] = (vn1[2] - mBDF2.c1 * un1[2] - mBDF2.c2 * un2[2])/mBDF2.c0;
                }
                else if (it_node -> IsFixed(DISPLACEMENT_Z) == false)
                {
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
     * It initializes time step solution. Only for reasons if the time step solution is restarted
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

        ProcessInfo current_process_info= rModelPart.GetProcessInfo();

        BaseType::InitializeSolutionStep(rModelPart, A, Dx, b);

        const double delta_time = current_process_info[DELTA_TIME];

        if (delta_time < 1.0e-24)
        {
            KRATOS_ERROR << " ERROR: detected delta_time = 0 in the Solution Scheme DELTA_TIME. PLEASE : check if the time step is created correctly for the current model part ";
        }
        
        // The BDF2 constants
        mBDF2.c0 =  1.5/delta_time;
        mBDF2.c1 = -2.0/delta_time;
        mBDF2.c2 =  0.5/delta_time;
        
        KRATOS_CATCH( "" );
    }

    /**
     * It initializes a non-linear iteration (for the element)
     * @param rModelPart The model of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */

    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;

        ProcessInfo& current_process_info = rModelPart.GetProcessInfo();
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++)
        {
            auto itElem = rModelPart.ElementsBegin() + i;
            itElem->InitializeNonLinearIteration(current_process_info);
        }
        
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++)
        {
            auto itElem = rModelPart.ConditionsBegin() + i;
            itElem->InitializeNonLinearIteration(current_process_info);
        }     
        
        KRATOS_CATCH( "" );
    }

    /**
     * It initializes a non-linear iteration (for an individual condition)
     * @param rCurrentConditiont The condition to compute
     * @param CurrentProcessInfo The current process info instance
     */

    void InitializeNonLinearIteration(
        Condition::Pointer rCurrentCondition,
        ProcessInfo& CurrentProcessInfo
    ) override
    {
        (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }

    /**
     * It initializes a non-linear iteration (for an individual element)
     * @param rCurrentElement The element to compute
     * @param CurrentProcessInfo The current process info instance
     */

    void InitializeNonLinearIteration(
        Element::Pointer rCurrentElement,
        ProcessInfo& CurrentProcessInfo
    ) override
    {
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }

    /**
     * This function is designed to be called in the builder and solver to introduce
     * @param rCurrentElement The element to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param CurrentProcessInfo The current process info instance
     */

    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const int thread = OpenMPUtils::ThisThread();

        //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        (rCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread],CurrentProcessInfo);

        (rCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

        AddDynamicsToLHS (LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    /**
     * This function is designed to calculate just the RHS contribution
     * @param rCurrentElemen The element to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param CurrentProcessInfo The current process info instance
     */

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {

        KRATOS_TRY;

        const int thread = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current element
        // (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        // Basic operations for the element considered
        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

        (rCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    /**
     * Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCurrentCondition The condition to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param CurrentProcessInfo The current process info instance
     */

    void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        const int thread = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current condition
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

        // Basic operations for the condition considered
        (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);

        (rCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

        (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

        AddDynamicsToLHS  (LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        // AssembleTimeSpaceLHS_Condition(rCurrentCondition, LHS_Contribution,DampMatrix, MassMatrix,CurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    /**
     * Functions that calculates the RHS of a "condition" object
     * @param rCurrentCondition The condition to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the condition degrees of freedom
     * @param CurrentProcessInfo The current process info instance
     */

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        const int thread = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current condition
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

        // Basic operations for the condition considered
        (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);

        (rCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

        (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread], CurrentProcessInfo);

        // Adding the dynamic contributions (static is already included)
        AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model of the problem to solve
     * @return Zero means  all ok
     */

    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        const int err = Scheme<TSparseSpace, TDenseSpace>::Check(rModelPart);
        if(err!=0) return err;

        // Check for variables keys
        // Verify that the variables are correctly initialized
        if(DISPLACEMENT.Key() == 0)
        {
            KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        }
        if(VELOCITY.Key() == 0)
        {
            KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is correctly registered" << std::endl;
        }
        if(ACCELERATION.Key() == 0)
        {
            KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;
        }

        // Check that variables are correctly allocated
        for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin();
                it!=rModelPart.NodesEnd(); it++)
        {
            if (it->SolutionStepsDataHas(DISPLACEMENT) == false)
            {
                KRATOS_ERROR << "DISPLACEMENT variable is not allocated for node " << it->Id() << std::endl;
            }
            if (it->SolutionStepsDataHas(VELOCITY) == false)
            {
                KRATOS_ERROR << "VELOCITY variable is not allocated for node " << it->Id() << std::endl;
            }
            if (it->SolutionStepsDataHas(ACCELERATION) == false)
            {
                KRATOS_ERROR << "ACCELERATION variable is not allocated for node " << it->Id() << std::endl;
            }
        }

        // Check that dofs exist
        for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin();
                it!=rModelPart.NodesEnd(); it++)
        {
            if(it->HasDofFor(DISPLACEMENT_X) == false)
            {
                KRATOS_ERROR << "missing DISPLACEMENT_X dof on node " << it->Id() << std::endl;
            }
            if(it->HasDofFor(DISPLACEMENT_Y) == false)
            {
                KRATOS_ERROR << "missing DISPLACEMENT_Y dof on node " << it->Id() << std::endl;
            }
            if(it->HasDofFor(DISPLACEMENT_Z) == false)
            {
                KRATOS_ERROR << "missing DISPLACEMENT_Z dof on node " << it->Id() << std::endl;
            }
        }

        // Check for minimum value of the buffer index
        // Verify buffer size
        if (rModelPart.GetBufferSize() < 2)
        {
            KRATOS_ERROR << "insufficient buffer size. Buffer size should be greater than 2. Current size is" << rModelPart.GetBufferSize() << std::endl;
        }

        return 0;
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
        double c0;
        double c1;
        double c2;
    };

    struct  GeneralMatrices
    {
        std::vector< Matrix > M;     // First derivative matrix  (usually mass matrix)
        std::vector< Matrix > D;     // Second derivative matrix (usually damping matrix)
    };

    struct GeneralVectors
    {
        std::vector< Vector > v;    // Velocity
        std::vector< Vector > a;    // Acceleration
    };

    BDF2Method        mBDF2; // The BDF2 coefficients
    GeneralMatrices mMatrix; // This contains the auxiliar matrices
    GeneralVectors  mVector; // This contains the auxiliar derivatives

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Updating first time Derivative
     * @param vn0 the velocity at the current step
     * @param un0 the displacemnt at the current step
     * @param un1 the displacemnt one step in the past
     * @param un2 the displacemnt two steps in the past
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
     * Updating second time Derivative
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
     * It adds the dynamic LHS contribution of the elements LHS = d(-RHS)/d(un0) = c0*c0*M + c0*D + K
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
        )
    {
        // Adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += M * std::pow(mBDF2.c0, 2);
        }

        // Adding  damping contribution
        if (D.size1() != 0) // if D matrix declared
        {
            noalias(LHS_Contribution) += D * mBDF2.c0;
        }
    }

    /**
     * It adds the dynamic RHS contribution of the elements b - M*a - D*v
     * @param rCurrentElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param CurrentProcessInfo The current process info instance
     */

    void AddDynamicsToRHS(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo
        )
    {
        const int thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentElement->GetSecondDerivativesVector(mVector.a[thread], 0);

            noalias(RHS_Contribution) -= prod(M, mVector.a[thread]);
        }

        // Adding damping contribution
        if (D.size1() != 0)
        {
            rCurrentElement->GetFirstDerivativesVector(mVector.v[thread], 0);

            noalias(RHS_Contribution) -= prod(D, mVector.v[thread]);
        }
    }

    /**
     * It adds the dynamic RHS contribution of the condition RHS = fext - M*an0 - D*vn0 - K*dn0
     * @param rCurrentCondition The condition to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param CurrentProcessInfo The current process info instance
     */

    void AddDynamicsToRHS(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo
        )
    {
        const int thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentCondition->GetSecondDerivativesVector(mVector.a[thread], 0);

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
        }

        // Adding damping contribution
        // Damping contribution
        if (D.size1() != 0)
        {
            rCurrentCondition->GetFirstDerivativesVector(mVector.v[thread], 0);

            noalias(RHS_Contribution) -= prod(D, mVector.v [thread]);
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
