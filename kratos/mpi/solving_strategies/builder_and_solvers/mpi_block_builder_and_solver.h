//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef KRATOS_MPI_BLOCK_BUILDER_AND_SOLVER_H_INCLUDED
#define KRATOS_MPI_BLOCK_BUILDER_AND_SOLVER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class MPIBlockBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MPIBlockBuilderAndSolver
    KRATOS_CLASS_POINTER_DEFINITION(MPIBlockBuilderAndSolver);

    /// Definition of the base class
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    // The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Definition of the classes from the base class
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit MPIBlockBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) : BaseType(pNewLinearSystemSolver)
    {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "MPIBlockBuilderAndSolver"
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);
    }

    /**
     * @brief Default constructor.
     */
    explicit MPIBlockBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver)
    {
    }

    /// Destructor.
    ~MPIBlockBuilderAndSolver() override = default;

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
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
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
        TSystemMatrixType& A) override
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
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
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
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
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    /**
     * @brief Corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
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
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
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
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    /**
     * @brief Organises the dofset in order to speed up the building phase
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpSystem(
        ModelPart& rModelPart
    ) override
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
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
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
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
        // TODO only local constraints
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
        // TODO only local constraints
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
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    /**
     * @brief Applies the dirichlet conditions. This operation may be very heavy or completely
     * unexpensive depending on the implementation choosen and on how the System Matrix is built.
     * @details For explanation of how it works for a particular implementation the user
     * should refer to the particular Builder And Solver choosen
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    void ApplyConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb) override
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
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

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MPIBlockBuilderAndSolver" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MPIBlockBuilderAndSolver";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

    ///@}

private:
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MPIBlockBuilderAndSolver& operator=(MPIBlockBuilderAndSolver const& rOther){}

    /// Copy constructor.
    MPIBlockBuilderAndSolver(MPIBlockBuilderAndSolver const& rOther){}

    ///@}

}; // Class MPIBlockBuilderAndSolver

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MPI_BLOCK_BUILDER_AND_SOLVER_H_INCLUDED defined
