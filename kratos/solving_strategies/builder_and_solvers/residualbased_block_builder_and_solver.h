//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborators:   Vicente Mataix
//
//
#if !defined(KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER )
#define  KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER

/* System includes */
#include <unordered_set>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"
#include "includes/kratos_flags.h"
#include "includes/lock_object.h"
#include "utilities/sparse_matrix_multiplication_utility.h"

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
 * @class ResidualBasedEliminationBuilderAndSolver
 * @ingroup KratosCore
 * @brief Current class provides an implementation for standard builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information.
 * Calculation of the reactions involves a cost very similiar to the calculation of the total residual
 * @tparam TSparseSpace The sparse system considered
 * @tparam TDenseSpace The dense system considered
 * @tparam TLinearSolver The linear solver considered
 * @author Riccardo Rossi
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedBlockBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the flags
    KRATOS_DEFINE_LOCAL_FLAG( SILENT_WARNINGS );

    // Scaling enum
    enum class SCALING_DIAGONAL {NO_SCALING = 0, CONSIDER_NORM_DIAGONAL = 1, CONSIDER_MAX_DIAGONAL = 2, CONSIDER_PRESCRIBED_DIAGONAL = 3};

    /// Definition of the pointer
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolver);

    /// Definition of the base class
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// The definition of the current class
    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

    // The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Definition of the classes from the base class
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    /// Additional definitions
    typedef PointerVectorSet<Element, IndexedObject> ElementsContainerType;
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;
    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrixType;

    /// DoF types definition
    typedef Node<3> NodeType;
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit ResidualBasedBlockBuilderAndSolver() : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ResidualBasedBlockBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) : BaseType(pNewLinearSystemSolver)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Default constructor.
     */
    explicit ResidualBasedBlockBuilderAndSolver(typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver)
    {
        mScalingDiagonal = SCALING_DIAGONAL::NO_SCALING;
    }

    /** Destructor.
     */
    ~ResidualBasedBlockBuilderAndSolver() override
    {
    }

    /**
     * @brief Create method
     * @param pNewLinearSystemSolver The linear solver for the system of equations
     * @param ThisParameters The configuration parameters
     */
    typename BaseType::Pointer Create(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(pNewLinearSystemSolver,ThisParameters);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to perform the build of the RHS. The vector could be sized as the total number
     * of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param b The RHS vector
     */
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
        {
            # pragma omp for  schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    element_is_active = (it)->Is(ACTIVE);

                if (element_is_active)
                {
                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
                }

            }

            #pragma omp for  schedule(guided, 512)
            for (int k = 0; k < nconditions; k++)
            {
                ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    condition_is_active = (it)->Is(ACTIVE);

                if (condition_is_active)
                {
                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
                }
            }
        }

        const double stop_build = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Build time: " << stop_build - start_build << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building of the LHS
     * @details Depending on the implementation choosen the size of the matrix could
     * be equal to the total number of Dofs or to the number of unrestrained dofs
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     */
    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA
        ) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const auto it_elem_begin = rModelPart.ElementsBegin();
        const auto it_cond_begin = rModelPart.ConditionsBegin();

        // Contributions to the system
        LocalSystemMatrixType lhs_contribution(0, 0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType equation_id;

        // Assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel firstprivate(nelements, nconditions, lhs_contribution, equation_id )
        {
            # pragma omp for  schedule(guided, 512) nowait
            for (int k = 0; k < nelements; ++k) {
                auto it_elem = it_elem_begin + k;

                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool element_is_active = true;
                if (it_elem->IsDefined(ACTIVE))
                    element_is_active = it_elem->Is(ACTIVE);

                if (element_is_active) {
                    // Calculate elemental contribution
                    pScheme->CalculateLHSContribution(*it_elem, lhs_contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleLHS(rA, lhs_contribution, equation_id);
                }

            }

            #pragma omp for  schedule(guided, 512)
            for (int k = 0; k < nconditions; ++k) {
                auto it_cond = it_cond_begin + k;

                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool condition_is_active = true;
                if (it_cond->IsDefined(ACTIVE))
                    condition_is_active = it_cond->Is(ACTIVE);

                if (condition_is_active)
                {
                    // Calculate elemental contribution
                    pScheme->CalculateLHSContribution(*it_cond, lhs_contribution, equation_id, r_current_process_info);

                    // Assemble the elemental contribution
                    AssembleLHS(rA, lhs_contribution, equation_id);
                }
            }
        }

        const double stop_build = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >= 1) << "Build time LHS: " << stop_build - start_build << std::endl;


        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() > 2) << "Finished parallel building LHS" << std::endl;


        KRATOS_CATCH("")
    }

    /**
     * @brief Build a rectangular matrix of size n*N where "n" is the number of unrestrained degrees of freedom
     * and "N" is the total number of degrees of freedom involved.
     * @details This matrix is obtained by building the total matrix without the lines corresponding to the fixed
     * degrees of freedom (but keeping the columns!!)
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     */
    void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A) override
    {
        KRATOS_TRY

        TSystemVectorType tmp(A.size1(), 0.0);
        this->Build(pScheme, rModelPart, A, tmp);

        KRATOS_CATCH("")
    }

    /**
     * @brief This is a call to the linear system solver
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void SystemSolve(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY
        double norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            //do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
        }
        else
            TSparseSpace::SetToZero(Dx);

        if(mT.size1() != 0) //if there are master-slave constraints
        {
            //recover solution of the original problem
            TSystemVectorType Dxmodified = Dx;

            TSparseSpace::Mult(mT, Dxmodified, Dx);
        }

        //prints informations about the current time
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() > 1) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    virtual void SystemSolveWithPhysics(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        ModelPart& rModelPart
        )
    {
        if(rModelPart.MasterSlaveConstraints().size() != 0) {
            TSystemVectorType Dxmodified(rb.size());

            InternalSystemSolveWithPhysics(rA, Dxmodified, rb, rModelPart);

            //recover solution of the original problem
            TSparseSpace::Mult(mT, Dxmodified, rDx);
        } else {
            InternalSystemSolveWithPhysics(rA, rDx, rb, rModelPart);
        }
    }

    /**
      *@brief This is a call to the linear system solver (taking into account some physical particularities of the problem)
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     * @param rModelPart The model part of the problem to solve
     */
    void InternalSystemSolveWithPhysics(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        ModelPart& rModelPart
    )
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if (norm_b != 0.00) {
            //provide physical data as needed
            if(BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded() )
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(A, Dx, b, BaseType::mDofSet, rModelPart);

            //do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
        } else {
            TSparseSpace::SetToZero(Dx);
            KRATOS_WARNING_IF("ResidualBasedBlockBuilderAndSolver", mOptions.IsNot(SILENT_WARNINGS)) << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

        // Prints informations about the current time
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() > 1) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building and solving phase at the same time.
     * @details It is ideally the fastest and safer function to use when it is possible to solve
     * just after building
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        Timer::Start("Build");

        Build(pScheme, rModelPart, A, b);

        Timer::Stop("Build");

        if(rModelPart.MasterSlaveConstraints().size() != 0) {
            Timer::Start("ApplyConstraints");
            ApplyConstraints(pScheme, rModelPart, A, b);
            Timer::Stop("ApplyConstraints");
        }

        ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        const double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        SystemSolveWithPhysics(A, Dx, b, rModelPart);

        Timer::Stop("Solve");
        const double stop_solve = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() >=1 && rModelPart.GetCommunicator().MyPID() == 0)) << "System solve time: " << stop_solve - start_solve << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building and solving phase at the same time Linearizing with the database at the old iteration
     * @details It is ideally the fastest and safer function to use when it is possible to solve just after building
     * @param pScheme The pointer to the integration scheme
     * @param rModelPart The model part to compute
     * @param rA The LHS matrix of the system of equations
     * @param rDx The vector of unkowns
     * @param rb The RHS vector of the system of equations
     * @param MoveMesh tells if the update of the scheme needs  to be performed when calling the Update of the scheme
     */
    void BuildAndSolveLinearizedOnPreviousIteration(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        const bool MoveMesh
        ) override
    {
        KRATOS_INFO_IF("BlockBuilderAndSolver", this->GetEchoLevel() > 0)
            << "Linearizing on Old iteration" << std::endl;

        KRATOS_ERROR_IF(rModelPart.GetBufferSize() == 1) << "BlockBuilderAndSolver: \n"
                << "The buffer size needs to be at least 2 in order to use \n"
                << "BuildAndSolveLinearizedOnPreviousIteration \n"
                << "current buffer size for modelpart: " << rModelPart.Name() << std::endl
                << "is :" << rModelPart.GetBufferSize()
                << " Please set IN THE STRATEGY SETTINGS "
                << " UseOldStiffnessInFirstIteration=false " << std::endl;

        DofsArrayType fixed_dofs;
        for(auto& r_dof : BaseType::mDofSet){
            if(r_dof.IsFixed()){
                fixed_dofs.push_back(&r_dof);
                r_dof.FreeDof();
            }
        }

        //TODO: Here we need to take the vector from other ones because
        // We cannot create a trilinos vector without a communicator. To be improved!
        TSystemVectorType dx_prediction(rDx);
        TSystemVectorType rhs_addition(rb); //we know it is zero here, so we do not need to set it

        // Here we bring back the database to before the prediction,
        // but we store the prediction increment in dx_prediction.
        // The goal is that the stiffness is computed with the
        // converged configuration at the end of the previous step.
        const auto it_dof_begin = BaseType::mDofSet.begin();
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(BaseType::mDofSet.size()); ++i) {
            auto it_dof = it_dof_begin + i;
            //NOTE: this is initialzed to - the value of dx prediction
            dx_prediction[it_dof->EquationId()] = -(it_dof->GetSolutionStepValue() - it_dof->GetSolutionStepValue(1));
        }

        // Use UpdateDatabase to bring back the solution to how it was at the end of the previous step
        pScheme->Update(rModelPart, BaseType::mDofSet, rA, dx_prediction, rb);
        if (MoveMesh) {
            VariableUtils().UpdateCurrentPosition(rModelPart.Nodes(),DISPLACEMENT,0);
        }

        this->Build(pScheme, rModelPart, rA, rb);

        // Put back the prediction into the database
        TSparseSpace::InplaceMult(dx_prediction, -1.0); //change sign to dx_prediction
        TSparseSpace::UnaliasedAdd(rDx, 1.0, dx_prediction);

        // Use UpdateDatabase to bring back the solution
        // to where it was taking into account BCs
        // it is done here so that constraints are correctly taken into account right after
        pScheme->Update(rModelPart, BaseType::mDofSet, rA, dx_prediction, rb);
        if (MoveMesh) {
            VariableUtils().UpdateCurrentPosition(rModelPart.Nodes(),DISPLACEMENT,0);
        }


        // Apply rb -= A*dx_prediction
        TSparseSpace::Mult(rA, dx_prediction, rhs_addition);
        TSparseSpace::UnaliasedAdd(rb, -1.0, rhs_addition);

        for(auto& dof : fixed_dofs)
            dof.FixDof();

        if (!rModelPart.MasterSlaveConstraints().empty()) {
            this->ApplyConstraints(pScheme, rModelPart, rA, rb);
        }
        this->ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);
        this->SystemSolveWithPhysics(rA, rDx, rb, rModelPart);
    }

    /**
     * @brief Corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        BuildRHS(pScheme, rModelPart, rb);

        if(rModelPart.MasterSlaveConstraints().size() != 0) {
            Timer::Start("ApplyRHSConstraints");
            ApplyRHSConstraints(pScheme, rModelPart, rb);
            Timer::Stop("ApplyRHSConstraints");
        }

        ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        const double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        Timer::Stop("Solve");
        const double stop_solve = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() >=1 && rModelPart.GetCommunicator().MyPID() == 0)) << "System solve time: " << stop_solve - start_solve << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the build of the RHS.
     * @details The vector could be sized as the total number of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        Timer::Start("BuildRHS");

        BuildRHSNoDirichlet(pScheme,rModelPart,b);

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for (int k = 0; k<ndofs; k++)
        {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;
            const std::size_t i = dof_iterator->EquationId();

            if (dof_iterator->IsFixed())
                b[i] = 0.0;
        }

        Timer::Stop("BuildRHS");

        KRATOS_CATCH("")
    }

    /**
     * @brief Builds the list of the DofSets involved in the problem by "asking" to each element
     * and condition its Dofs.
     * @details The list of dofs is stores insde the BuilderAndSolver as it is closely connected to the
     * way the matrix and RHS are built
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
    ) override
    {
        KRATOS_TRY;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Setting up the dofs" << std::endl;

        //Gets the array of elements from the modeler
        ElementsArrayType& r_elements_array = rModelPart.Elements();
        const int number_of_elements = static_cast<int>(r_elements_array.size());

        DofsVectorType dof_list, second_dof_list; // NOTE: The second dof list is only used on constraints to include master/slave relations

        unsigned int nthreads = OpenMPUtils::GetNumThreads();

        typedef std::unordered_set < NodeType::DofType::Pointer, DofPointerHasher>  set_type;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() > 2)) << "Number of threads" << nthreads << "\n" << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() > 2)) << "Initializing element loop" << std::endl;

        /**
         * Here we declare three sets.
         * - The global set: Contains all the DoF of the system
         * - The slave set: The DoF that are not going to be solved, due to MPC formulation
         */
        set_type dof_global_set;
        dof_global_set.reserve(number_of_elements*20);

        #pragma omp parallel firstprivate(dof_list, second_dof_list)
        {
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // We cleate the temporal set and we reserve some space on them
            set_type dofs_tmp_set;
            dofs_tmp_set.reserve(20000);

            // Gets the array of elements from the modeler
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < number_of_elements; ++i) {
                auto it_elem = r_elements_array.begin() + i;

                // Gets list of Dof involved on every element
                pScheme->GetDofList(*it_elem, dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }

            // Gets the array of conditions from the modeler
            ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
            const int number_of_conditions = static_cast<int>(r_conditions_array.size());
            #pragma omp for  schedule(guided, 512) nowait
            for (int i = 0; i < number_of_conditions; ++i) {
                auto it_cond = r_conditions_array.begin() + i;

                // Gets list of Dof involved on every element
                pScheme->GetDofList(*it_cond, dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }

            // Gets the array of constraints from the modeler
            auto& r_constraints_array = rModelPart.MasterSlaveConstraints();
            const int number_of_constraints = static_cast<int>(r_constraints_array.size());
            #pragma omp for  schedule(guided, 512) nowait
            for (int i = 0; i < number_of_constraints; ++i) {
                auto it_const = r_constraints_array.begin() + i;

                // Gets list of Dof involved on every element
                it_const->GetDofList(dof_list, second_dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                dofs_tmp_set.insert(second_dof_list.begin(), second_dof_list.end());
            }

            // We merge all the sets in one thread
            #pragma omp critical
            {
                dof_global_set.insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
            }
        }

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        Doftemp.reserve(dof_global_set.size());
        for (auto it= dof_global_set.begin(); it!= dof_global_set.end(); it++)
        {
            Doftemp.push_back( *it );
        }
        Doftemp.Sort();

        BaseType::mDofSet = Doftemp;

        //Throws an exception if there are no Degrees Of Freedom involved in the analysis
        KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseType::mDofSet.size() << std::endl;

        BaseType::mDofSetIsInitialized = true;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished setting up the dofs" << std::endl;

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        // This is tobe done only in debug mode
        if (BaseType::GetCalculateReactionsFlag()) {
            for (auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator) {
                    KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction()) << "Reaction variable not set for the following : " <<std::endl
                        << "Node : "<<dof_iterator->Id()<< std::endl
                        << "Dof : "<<(*dof_iterator)<<std::endl<<"Not possible to calculate reactions."<<std::endl;
            }
        }
#endif

        KRATOS_CATCH("");
    }

    /**
     * @brief Organises the dofset in order to speed up the building phase
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpSystem(
        ModelPart& rModelPart
    ) override
    {
        //int free_id = 0;
        BaseType::mEquationSystemSize = BaseType::mDofSet.size();
        int ndofs = static_cast<int>(BaseType::mDofSet.size());

        #pragma omp parallel for firstprivate(ndofs)
        for (int i = 0; i < static_cast<int>(ndofs); i++) {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + i;
            dof_iterator->SetEquationId(i);

        }
    }

    //**************************************************************************
    //**************************************************************************

    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ModelPart& rModelPart
    ) override
    {
        KRATOS_TRY
        if (pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
            pA.swap(pNewA);
        }
        if (pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
            pDx.swap(pNewDx);
        }
        if (pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
            pb.swap(pNewb);
        }

        TSystemMatrixType& A = *pA;
        TSystemVectorType& Dx = *pDx;
        TSystemVectorType& b = *pb;

        //resizing the system vectors and matrix
        if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            ConstructMatrixStructure(pScheme, A, rModelPart);
        }
        else
        {
            if (A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
            {
                KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permited."<<std::endl;
                A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
                ConstructMatrixStructure(pScheme, A, rModelPart);
            }
        }
        if (Dx.size() != BaseType::mEquationSystemSize)
            Dx.resize(BaseType::mEquationSystemSize, false);
        TSparseSpace::SetToZero(Dx);
        if (b.size() != BaseType::mEquationSystemSize) {
            b.resize(BaseType::mEquationSystemSize, false);
        }
        TSparseSpace::SetToZero(b);

        ConstructMasterSlaveConstraintsStructure(rModelPart);

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        KRATOS_TRY

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // Getting process info
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // Computing constraints
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        #pragma omp parallel for schedule(guided, 512) firstprivate(n_constraints, constraints_begin)
        for (int k = 0; k < n_constraints; ++k) {
            auto it = constraints_begin + k;
            it->InitializeSolutionStep(r_process_info); // Here each constraint constructs and stores its T and C matrices. Also its equation slave_ids.
        }

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        // Getting process info
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // Computing constraints
        const int n_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
        const auto constraints_begin = rModelPart.MasterSlaveConstraintsBegin();
        #pragma omp parallel for schedule(guided, 512) firstprivate(n_constraints, constraints_begin)
        for (int k = 0; k < n_constraints; ++k) {
            auto it = constraints_begin + k;
            it->FinalizeSolutionStep(r_process_info);
        }
    }

    //**************************************************************************
    //**************************************************************************

    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        TSparseSpace::SetToZero(b);

        //refresh RHS to have the correct reactions
        BuildRHSNoDirichlet(pScheme, rModelPart, b);

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for (int k = 0; k<ndofs; k++) {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;

            const int i = (dof_iterator)->EquationId();
            (dof_iterator)->GetSolutionStepReactionValue() = -b[i];
        }
    }

    /**
     * @brief Applies the dirichlet conditions. This operation may be very heavy or completely
     * unexpensive depending on the implementation choosen and on how the System Matrix is built.
     * @details For explanation of how it works for a particular implementation the user
     * should refer to the particular Builder And Solver choosen
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        const std::size_t system_size = rA.size1();
        Vector scaling_factors (system_size);

        const auto it_dof_iterator_begin = BaseType::mDofSet.begin();
        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        // NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for (int k = 0; k<ndofs; k++) {
            auto it_dof_iterator = it_dof_iterator_begin + k;
            if (it_dof_iterator->IsFixed()) {
                scaling_factors[k] = 0.0;
            } else {
                scaling_factors[k] = 1.0;
            }
        }

        double* Avalues = rA.value_data().begin();
        std::size_t* Arow_indices = rA.index1_data().begin();
        std::size_t* Acol_indices = rA.index2_data().begin();

        // The diagonal considered
        mScaleFactor = GetScaleNorm(rModelPart, rA);

        // Detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        #pragma omp parallel firstprivate(system_size)
        {
            std::size_t col_begin = 0, col_end  = 0;
            bool empty = true;

            #pragma omp for
            for (int k = 0; k < static_cast<int>(system_size); ++k) {
                col_begin = Arow_indices[k];
                col_end = Arow_indices[k + 1];
                empty = true;
                for (std::size_t j = col_begin; j < col_end; ++j) {
                    if(Avalues[j] != 0.0) {
                        empty = false;
                        break;
                    }
                }

                if(empty) {
                    rA(k, k) = mScaleFactor;
                    rb[k] = 0.0;
                }
            }
        }

        #pragma omp parallel for firstprivate(system_size)
        for (int k = 0; k < static_cast<int>(system_size); ++k) {
            std::size_t col_begin = Arow_indices[k];
            std::size_t col_end = Arow_indices[k+1];
            const double k_factor = scaling_factors[k];
            if (k_factor == 0.0) {
                // Zero out the whole row, except the diagonal
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if (static_cast<int>(Acol_indices[j]) != k )
                        Avalues[j] = 0.0;

                // Zero out the RHS
                rb[k] = 0.0;
            } else {
                // Zero out the column which is associated with the zero'ed row
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if(scaling_factors[ Acol_indices[j] ] == 0 )
                        Avalues[j] = 0.0;
            }
        }
    }

    /**
     * @brief Applies the constraints with master-slave relation matrix (RHS only)
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rb The RHS vector
     */
    void ApplyRHSConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            BuildMasterSlaveConstraints(rModelPart);

            // We compute the transposed matrix of the global relation matrix
            TSystemMatrixType T_transpose_matrix(mT.size2(), mT.size1());
            SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, mT, 1.0);

            TSystemVectorType b_modified(rb.size());
            TSparseSpace::Mult(T_transpose_matrix, rb, b_modified);
            TSparseSpace::Copy(b_modified, rb);

            // Apply diagonal values on slaves
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(mSlaveIds.size()); ++i) {
                const IndexType slave_equation_id = mSlaveIds[i];
                if (mInactiveSlaveDofs.find(slave_equation_id) == mInactiveSlaveDofs.end()) {
                    rb[slave_equation_id] = 0.0;
                }
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Applies the constraints with master-slave relation matrix
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     */
    void ApplyConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            BuildMasterSlaveConstraints(rModelPart);

            // We compute the transposed matrix of the global relation matrix
            TSystemMatrixType T_transpose_matrix(mT.size2(), mT.size1());
            SparseMatrixMultiplicationUtility::TransposeMatrix<TSystemMatrixType, TSystemMatrixType>(T_transpose_matrix, mT, 1.0);

            TSystemVectorType b_modified(rb.size());
            TSparseSpace::Mult(T_transpose_matrix, rb, b_modified);
            TSparseSpace::Copy(b_modified, rb);

            TSystemMatrixType auxiliar_A_matrix(mT.size2(), rA.size2());
            SparseMatrixMultiplicationUtility::MatrixMultiplication(T_transpose_matrix, rA, auxiliar_A_matrix); //auxiliar = T_transpose * rA
            T_transpose_matrix.resize(0, 0, false);                                                             //free memory

            SparseMatrixMultiplicationUtility::MatrixMultiplication(auxiliar_A_matrix, mT, rA); //A = auxilar * T   NOTE: here we are overwriting the old A matrix!
            auxiliar_A_matrix.resize(0, 0, false);                                              //free memory

            const double max_diag = GetMaxDiagonal(rA);

            // Apply diagonal values on slaves
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(mSlaveIds.size()); ++i) {
                const IndexType slave_equation_id = mSlaveIds[i];
                if (mInactiveSlaveDofs.find(slave_equation_id) == mInactiveSlaveDofs.end()) {
                    rA(slave_equation_id, slave_equation_id) = max_diag;
                    rb[slave_equation_id] = 0.0;
                }
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        BaseType::Clear();

        mSlaveIds.clear();
        mMasterIds.clear();
        mInactiveSlaveDofs.clear();
        mT.resize(0,0,false);
        mConstantVector.resize(0,false);
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model part of the problem to solve
     * @return 0 all ok
     */
    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        return 0;
        KRATOS_CATCH("");
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                 : "block_builder_and_solver",
            "block_builder"                        : true,
            "diagonal_values_for_dirichlet_dofs"   : "use_max_diagonal",
            "silent_warnings"                      : false
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "block_builder_and_solver";
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
        return "ResidualBasedBlockBuilderAndSolver";
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

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    TSystemMatrixType mT;              /// This is matrix containing the global relation for the constraints
    TSystemVectorType mConstantVector; /// This is vector containing the rigid movement of the constraint
    std::vector<IndexType> mSlaveIds;  /// The equation ids of the slaves
    std::vector<IndexType> mMasterIds; /// The equation ids of the master
    std::unordered_set<IndexType> mInactiveSlaveDofs; /// The set containing the inactive slave dofs
    double mScaleFactor = 1.0;         /// The manuallyset scale factor

    SCALING_DIAGONAL mScalingDiagonal; /// We identify the scaling considered for the dirichlet dofs
    Flags mOptions;                    /// Some flags used internally

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void BuildRHSNoDirichlet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        //Getting the Elements
        ElementsArrayType& pElements = rModelPart.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = rModelPart.Conditions();

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        //for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)

        const int nelements = static_cast<int>(pElements.size());
        #pragma omp parallel firstprivate(nelements, RHS_Contribution, EquationId)
        {
            #pragma omp for schedule(guided, 512) nowait
            for (int i=0; i<nelements; i++) {
                typename ElementsArrayType::iterator it = pElements.begin() + i;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if( (it)->IsDefined(ACTIVE) ) {
                    element_is_active = (it)->Is(ACTIVE);
                }

                if(element_is_active) {
                    //calculate elemental Right Hand Side Contribution
                    pScheme->CalculateRHSContribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }

            LHS_Contribution.resize(0, 0, false);
            RHS_Contribution.resize(0, false);

            // assemble all conditions
            const int nconditions = static_cast<int>(ConditionsArray.size());
            #pragma omp for schedule(guided, 512)
            for (int i = 0; i<nconditions; i++) {
                auto it = ConditionsArray.begin() + i;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if( (it)->IsDefined(ACTIVE) ) {
                    condition_is_active = (it)->Is(ACTIVE);
                }

                if(condition_is_active) {
                    //calculate elemental contribution
                    pScheme->CalculateRHSContribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }
        }

        KRATOS_CATCH("")

    }

    virtual void ConstructMasterSlaveConstraintsStructure(ModelPart& rModelPart)
    {
        if (rModelPart.MasterSlaveConstraints().size() > 0) {
            Timer::Start("ConstraintsRelationMatrixStructure");
            const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

            // Vector containing the localization in the system of the different terms
            DofsVectorType slave_dof_list, master_dof_list;

            // Constraint initial iterator
            const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();
            std::vector<std::unordered_set<IndexType>> indices(BaseType::mDofSet.size());

            std::vector<LockObject> lock_array(indices.size());

            #pragma omp parallel firstprivate(slave_dof_list, master_dof_list)
            {
                Element::EquationIdVectorType slave_ids(3);
                Element::EquationIdVectorType master_ids(3);
                std::unordered_map<IndexType, std::unordered_set<IndexType>> temp_indices;

                #pragma omp for schedule(guided, 512) nowait
                for (int i_const = 0; i_const < static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i_const) {
                    auto it_const = it_const_begin + i_const;

                    // Detect if the constraint is active or not. If the user did not make any choice the constraint
                    // It is active by default
                    bool constraint_is_active = true;
                    if( it_const->IsDefined(ACTIVE) ) {
                        constraint_is_active = it_const->Is(ACTIVE);
                    }

                    if(constraint_is_active) {
                        it_const->EquationIdVector(slave_ids, master_ids, r_current_process_info);

                        // Slave DoFs
                        for (auto &id_i : slave_ids) {
                            temp_indices[id_i].insert(master_ids.begin(), master_ids.end());
                        }
                    }
                }

                // Merging all the temporal indexes
                for (int i = 0; i < static_cast<int>(temp_indices.size()); ++i) {
                    lock_array[i].SetLock();
                    indices[i].insert(temp_indices[i].begin(), temp_indices[i].end());
                    lock_array[i].UnSetLock();
                }
            }

            mSlaveIds.clear();
            mMasterIds.clear();
            for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
                if (indices[i].size() == 0) // Master dof!
                    mMasterIds.push_back(i);
                else // Slave dof
                    mSlaveIds.push_back(i);
                indices[i].insert(i); // Ensure that the diagonal is there in T
            }

            // Count the row sizes
            std::size_t nnz = 0;
            for (IndexType i = 0; i < indices.size(); ++i)
                nnz += indices[i].size();

            mT = TSystemMatrixType(indices.size(), indices.size(), nnz);
            mConstantVector.resize(indices.size(), false);

            double *Tvalues = mT.value_data().begin();
            IndexType *Trow_indices = mT.index1_data().begin();
            IndexType *Tcol_indices = mT.index2_data().begin();

            // Filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
            Trow_indices[0] = 0;
            for (int i = 0; i < static_cast<int>(mT.size1()); i++)
                Trow_indices[i + 1] = Trow_indices[i] + indices[i].size();

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(mT.size1()); ++i) {
                const IndexType row_begin = Trow_indices[i];
                const IndexType row_end = Trow_indices[i + 1];
                IndexType k = row_begin;
                for (auto it = indices[i].begin(); it != indices[i].end(); ++it) {
                    Tcol_indices[k] = *it;
                    Tvalues[k] = 0.0;
                    k++;
                }

                indices[i].clear(); //deallocating the memory

                std::sort(&Tcol_indices[row_begin], &Tcol_indices[row_end]);
            }

            mT.set_filled(indices.size() + 1, nnz);

            Timer::Stop("ConstraintsRelationMatrixStructure");
        }
    }

    virtual void BuildMasterSlaveConstraints(ModelPart& rModelPart)
    {
        KRATOS_TRY

        TSparseSpace::SetToZero(mT);
        TSparseSpace::SetToZero(mConstantVector);

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Vector containing the localization in the system of the different terms
        DofsVectorType slave_dof_list, master_dof_list;

        // Contributions to the system
        Matrix transformation_matrix = LocalSystemMatrixType(0, 0);
        Vector constant_vector = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType slave_equation_ids, master_equation_ids;

        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

        // We clear the set
        mInactiveSlaveDofs.clear();

        #pragma omp parallel firstprivate(transformation_matrix, constant_vector, slave_equation_ids, master_equation_ids)
        {
            std::unordered_set<IndexType> auxiliar_inactive_slave_dofs;

            #pragma omp for schedule(guided, 512)
            for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
                auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;

                // Detect if the constraint is active or not. If the user did not make any choice the constraint
                // It is active by default
                bool constraint_is_active = true;
                if (it_const->IsDefined(ACTIVE))
                    constraint_is_active = it_const->Is(ACTIVE);

                if (constraint_is_active) {
                    it_const->CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);
                    it_const->EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);

                    for (IndexType i = 0; i < slave_equation_ids.size(); ++i) {
                        const IndexType i_global = slave_equation_ids[i];

                        // Assemble matrix row
                        AssembleRowContribution(mT, transformation_matrix, i_global, i, master_equation_ids);

                        // Assemble constant vector
                        const double constant_value = constant_vector[i];
                        double& r_value = mConstantVector[i_global];
                        #pragma omp atomic
                        r_value += constant_value;
                    }
                } else { // Taking into account inactive constraints
                    it_const->EquationIdVector(slave_equation_ids, master_equation_ids, r_current_process_info);
                    auxiliar_inactive_slave_dofs.insert(slave_equation_ids.begin(), slave_equation_ids.end());
                }
            }

            // We merge all the sets in one thread
            #pragma omp critical
            {
                mInactiveSlaveDofs.insert(auxiliar_inactive_slave_dofs.begin(), auxiliar_inactive_slave_dofs.end());
            }
        }

        // Setting the master dofs into the T and C system
        for (auto eq_id : mMasterIds) {
            mConstantVector[eq_id] = 0.0;
            mT(eq_id, eq_id) = 1.0;
        }

        // Setting inactive slave dofs in the T and C system
        for (auto eq_id : mInactiveSlaveDofs) {
            mConstantVector[eq_id] = 0.0;
            mT(eq_id, eq_id) = 1.0;
        }

        KRATOS_CATCH("")
    }

    virtual void ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType& A,
        ModelPart& rModelPart)
    {
        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        const std::size_t equation_size = BaseType::mEquationSystemSize;

        std::vector< LockObject > lock_array(equation_size);

        std::vector<std::unordered_set<std::size_t> > indices(equation_size);
        #pragma omp parallel for firstprivate(equation_size)
        for (int iii = 0; iii < static_cast<int>(equation_size); iii++) {
            indices[iii].reserve(40);
        }

        Element::EquationIdVectorType ids(3, 0);

        #pragma omp parallel for firstprivate(nelements, ids)
        for (int iii=0; iii<nelements; iii++) {
            typename ElementsContainerType::iterator i_element = el_begin + iii;
            pScheme->EquationId(*i_element, ids, CurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++) {
                lock_array[ids[i]].SetLock();
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
                lock_array[ids[i]].UnSetLock();
            }
        }

        #pragma omp parallel for firstprivate(nconditions, ids)
        for (int iii = 0; iii<nconditions; iii++) {
            typename ConditionsArrayType::iterator i_condition = cond_begin + iii;
            pScheme->EquationId(*i_condition, ids, CurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++) {
                lock_array[ids[i]].SetLock();
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
                lock_array[ids[i]].UnSetLock();
            }
        }

        //destroy locks
        lock_array = std::vector< LockObject >();

        //count the row sizes
        unsigned int nnz = 0;
        for (unsigned int i = 0; i < indices.size(); i++) {
            nnz += indices[i].size();
        }

        A = CompressedMatrixType(indices.size(), indices.size(), nnz);

        double* Avalues = A.value_data().begin();
        std::size_t* Arow_indices = A.index1_data().begin();
        std::size_t* Acol_indices = A.index2_data().begin();

        //filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
        Arow_indices[0] = 0;
        for (int i = 0; i < static_cast<int>(A.size1()); i++) {
            Arow_indices[i+1] = Arow_indices[i] + indices[i].size();
        }

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(A.size1()); i++) {
            const unsigned int row_begin = Arow_indices[i];
            const unsigned int row_end = Arow_indices[i+1];
            unsigned int k = row_begin;
            for (auto it = indices[i].begin(); it != indices[i].end(); it++) {
                Acol_indices[k] = *it;
                Avalues[k] = 0.0;
                k++;
            }

            indices[i].clear(); //deallocating the memory

            std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

        }

        A.set_filled(indices.size()+1, nnz);

        Timer::Stop("MatrixStructure");
    }

    void Assemble(
        TSystemMatrixType& A,
        TSystemVectorType& b,
        const LocalSystemMatrixType& LHS_Contribution,
        const LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++) {
            unsigned int i_global = EquationId[i_local];

            double& r_a = b[i_global];
            const double& v_a = RHS_Contribution(i_local);
            #pragma omp atomic
            r_a += v_a;

            AssembleRowContribution(A, LHS_Contribution, i_global, i_local, EquationId);
        }
    }

    //**************************************************************************

    void AssembleLHS(
        TSystemMatrixType& rA,
        const LocalSystemMatrixType& rLHSContribution,
        Element::EquationIdVectorType& rEquationId
        )
    {
        const SizeType local_size = rLHSContribution.size1();

        for (IndexType i_local = 0; i_local < local_size; i_local++) {
            const IndexType i_global = rEquationId[i_local];
            AssembleRowContribution(rA, rLHSContribution, i_global, i_local, rEquationId);
        }
    }

    //**************************************************************************

    void AssembleRHS(
        TSystemVectorType& b,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = RHS_Contribution.size();

        for (unsigned int i_local = 0; i_local < local_size; i_local++) {
            unsigned int i_global = EquationId[i_local];

            // ASSEMBLING THE SYSTEM VECTOR
            double& b_value = b[i_global];
            const double& rhs_value = RHS_Contribution[i_local];

            #pragma omp atomic
            b_value += rhs_value;
        }
    }

    inline void AssembleRowContribution(TSystemMatrixType& A, const Matrix& Alocal, const unsigned int i, const unsigned int i_local, Element::EquationIdVectorType& EquationId)
    {
        double* values_vector = A.value_data().begin();
        std::size_t* index1_vector = A.index1_data().begin();
        std::size_t* index2_vector = A.index2_data().begin();

        size_t left_limit = index1_vector[i];
//    size_t right_limit = index1_vector[i+1];

        //find the first entry
        size_t last_pos = ForwardFind(EquationId[0],left_limit,index2_vector);
        size_t last_found = EquationId[0];

        double& r_a = values_vector[last_pos];
        const double& v_a = Alocal(i_local,0);
        #pragma omp atomic
        r_a +=  v_a;

        //now find all of the other entries
        size_t pos = 0;
        for (unsigned int j=1; j<EquationId.size(); j++) {
            unsigned int id_to_find = EquationId[j];
            if(id_to_find > last_found) {
                pos = ForwardFind(id_to_find,last_pos+1,index2_vector);
            } else if(id_to_find < last_found) {
                pos = BackwardFind(id_to_find,last_pos-1,index2_vector);
            } else {
                pos = last_pos;
            }

            double& r = values_vector[pos];
            const double& v = Alocal(i_local,j);
            #pragma omp atomic
            r +=  v;

            last_found = id_to_find;
            last_pos = pos;
        }
    }

    /**
     * @brief This method returns the scale norm considering for scaling the diagonal
     * @param rModelPart The problem model part
     * @param rA The LHS matrix
     * @return The scale norm
     */
    double GetScaleNorm(
        ModelPart& rModelPart,
        TSystemMatrixType& rA
        )
    {
        switch (mScalingDiagonal) {
            case SCALING_DIAGONAL::NO_SCALING:
                return 1.0;
            case SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL: {
                const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
                KRATOS_ERROR_IF_NOT(r_current_process_info.Has(BUILD_SCALE_FACTOR)) << "Scale factor not defined at process info" << std::endl;
                return r_current_process_info.GetValue(BUILD_SCALE_FACTOR);
            }
            case SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL:
                return GetDiagonalNorm(rA)/static_cast<double>(rA.size1());
            case SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL:
                return GetMaxDiagonal(rA);
//                 return TSparseSpace::TwoNorm(rA)/static_cast<double>(rA.size1());
            default:
                return GetMaxDiagonal(rA);
        }
    }

    /**
     * @brief This method returns the diagonal norm considering for scaling the diagonal
     * @param rA The LHS matrix
     * @return The diagonal norm
     */
    double GetDiagonalNorm(TSystemMatrixType& rA)
    {
        double diagonal_norm = 0.0;
        #pragma omp parallel for reduction(+:diagonal_norm)
        for(int i = 0; i < static_cast<int>(TSparseSpace::Size1(rA)); ++i) {
            diagonal_norm += std::pow(rA(i,i), 2);
        }
        return std::sqrt(diagonal_norm);
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @return The diagonal  max value
     */
    double GetAveragevalueDiagonal(TSystemMatrixType& rA)
    {
        return 0.5 * (GetMaxDiagonal(rA) + GetMinDiagonal(rA));
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @return The diagonal  max value
     */
    double GetMaxDiagonal(TSystemMatrixType& rA)
    {
//         // NOTE: Reduction failing in MSVC
//         double max_diag = 0.0;
//         #pragma omp parallel for reduction(max:max_diag)
//         for(int i = 0; i < static_cast<int>(TSparseSpace::Size1(rA)); ++i) {
//             max_diag = std::max(max_diag, std::abs(rA(i,i)));
//         }
//         return max_diag;

        // Creating a buffer for parallel vector fill
        const int num_threads = OpenMPUtils::GetNumThreads();
        Vector max_vector(num_threads, 0.0);
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(TSparseSpace::Size1(rA)); ++i) {
            const int id = OpenMPUtils::ThisThread();
            const double abs_value_ii = std::abs(rA(i,i));
            if (abs_value_ii > max_vector[id])
                max_vector[id] = abs_value_ii;
        }

        double max_diag = 0.0;
        for(int i = 0; i < num_threads; ++i) {
            max_diag = std::max(max_diag, max_vector[i]);
        }
        return max_diag;
    }

    /**
     * @brief This method returns the diagonal min value
     * @param rA The LHS matrix
     * @return The diagonal min value
     */
    double GetMinDiagonal(TSystemMatrixType& rA)
    {
//         // NOTE: Reduction failing in MSVC
//         double min_diag = std::numeric_limits<double>::max();
//         #pragma omp parallel for reduction(min:min_diag)
//         for(int i = 0; i < static_cast<int>(TSparseSpace::Size1(rA)); ++i) {
//             min_diag = std::min(min_diag, std::abs(rA(i,i)));
//         }
//         return min_diag;

        // Creating a buffer for parallel vector fill
        const int num_threads = OpenMPUtils::GetNumThreads();
        Vector min_vector(num_threads, std::numeric_limits<double>::max());
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(TSparseSpace::Size1(rA)); ++i) {
            const int id = OpenMPUtils::ThisThread();
            const double abs_value_ii = std::abs(rA(i,i));
            if (abs_value_ii < min_vector[id])
                min_vector[id] = abs_value_ii;
        }

        double min_diag = std::numeric_limits<double>::max();
        for(int i = 0; i < num_threads; ++i) {
            min_diag = std::min(min_diag, min_vector[i]);
        }
        return min_diag;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        // Setting flags<
        const std::string& r_diagonal_values_for_dirichlet_dofs = ThisParameters["diagonal_values_for_dirichlet_dofs"].GetString();

        std::set<std::string> available_options_for_diagonal = {"no_scaling","use_max_diagonal","use_diagonal_norm","defined_in_process_info"};

        if (available_options_for_diagonal.find(r_diagonal_values_for_dirichlet_dofs) == available_options_for_diagonal.end()) {
            std::stringstream msg;
            msg << "Currently prescribed diagonal values for dirichlet dofs : " << r_diagonal_values_for_dirichlet_dofs << "\n";
            msg << "Admissible values for the diagonal scaling are : no_scaling, use_max_diagonal, use_diagonal_norm, or defined_in_process_info" << "\n";
            KRATOS_ERROR << msg.str() << std::endl;
        }

        // The first option will not consider any scaling (the diagonal values will be replaced with 1)
        if (r_diagonal_values_for_dirichlet_dofs == "no_scaling") {
            mScalingDiagonal = SCALING_DIAGONAL::NO_SCALING;
        } else if (r_diagonal_values_for_dirichlet_dofs == "use_max_diagonal") {
            mScalingDiagonal = SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL;
        } else if (r_diagonal_values_for_dirichlet_dofs == "use_diagonal_norm") { // On this case the norm of the diagonal will be considered
            mScalingDiagonal = SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL;
        } else { // Otherwise we will assume we impose a numerical value
            mScalingDiagonal = SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL;
        }
        mOptions.Set(SILENT_WARNINGS, ThisParameters["silent_warnings"].GetBool());
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

    ///@}

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

    inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate)
    {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();
        while (i != endit && (*i) != candidate) {
            i++;
        }
        if (i == endit) {
            v.push_back(candidate);
        }
    }

    //******************************************************************************************
    //******************************************************************************************

    inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, DenseVector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (unsigned int i = 1; i < number_of_threads; i++) {
            partitions[i] = partitions[i - 1] + partition_size;
        }
    }

    inline unsigned int ForwardFind(const unsigned int id_to_find,
                                    const unsigned int start,
                                    const size_t* index_vector)
    {
        unsigned int pos = start;
        while(id_to_find != index_vector[pos]) pos++;
        return pos;
    }

    inline unsigned int BackwardFind(const unsigned int id_to_find,
                                     const unsigned int start,
                                     const size_t* index_vector)
    {
        unsigned int pos = start;
        while(id_to_find != index_vector[pos]) pos--;
        return pos;
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class ResidualBasedBlockBuilderAndSolver */

///@}

///@name Type Definitions
///@{

// Here one should use the KRATOS_CREATE_LOCAL_FLAG, but it does not play nice with template parameters
template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
const Kratos::Flags ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::SILENT_WARNINGS(Kratos::Flags::Create(0));

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER  defined */
