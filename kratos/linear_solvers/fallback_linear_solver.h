//
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
#include <vector>

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "factories/linear_solver_factory.h"

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
 * @class FallbackLinearSolver
 * @ingroup KratosCore
 * @brief This is a simple linear solver that tries to solve the system using a list of solvers in order, until one of them succeeds.
 * @details Stores a vector of solvers and tries to solve the system using the first one. If the solver fails, it tries the next one, and so on.
 * @tparam TSparseSpaceType which specify type of the unknowns, coefficients, sparse matrix, vector of unknowns, right hand side vector and their respective operators.
 * @tparam TDenseMatrixType which specify type of the matrices used as temporary matrices or multi solve unknowns and right hand sides and their operators.
 * @tparam TReordererType which specify type of the Orderer that performs the reordering of matrix to optimize the solution.
 * @see SparseSpace
 * @see DenseSpace
 * @see Reorderer
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class FallbackLinearSolver
    : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearSolver
    KRATOS_CLASS_POINTER_DEFINITION(FallbackLinearSolver);

    /// Definition of the base class
    using BaseType = LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>;

    /// Definition of the base type pointer
    using LinearSolverPointer = typename BaseType::Pointer;

    /// The definition of the linear solver factory
    using LinearSolverFactoryType = LinearSolverFactory<TSparseSpaceType, TDenseSpaceType>;

    /// Type definition for sparse matrix
    using SparseMatrixType = typename BaseType::SparseMatrixType;

    /// Type definition for pointer to sparse matrix
    using SparseMatrixPointerType = typename BaseType::SparseMatrixPointerType;

    /// Type definition for vector
    using VectorType = typename BaseType::VectorType;

    /// Type definition for pointer to vector
    using VectorPointerType = typename BaseType::VectorPointerType;

    /// Type definition for dense matrix
    using DenseMatrixType = typename BaseType::DenseMatrixType;

    /// Type definition for dense vector
    using DenseVectorType = typename BaseType::DenseVectorType;

    /// Type definition for size
    using SizeType = std::size_t;

    /// Type definition for index
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief This is the default constructor
     */
    explicit FallbackLinearSolver() = default;

    /**
     * @brief This is the default constructor
     * @param ThisParameters The configuration parameters
     */
    explicit FallbackLinearSolver(Parameters ThisParameters)
        : mParameters(ThisParameters)
    {
        // Set the default parameters
        mParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Define the solvers
        const SizeType number_of_solver_parameters = mParameters["solvers"].size();
        KRATOS_ERROR_IF(number_of_solver_parameters == 0) << "No solvers defined in the parameters" << std::endl;
        mSolvers.reserve(number_of_solver_parameters);
        for (IndexType i = 0; i < number_of_solver_parameters; ++i) {
            mSolvers.push_back(ConstructLinearSolverFromSettings(mParameters["solvers"][i]));
            mRequiresAdditionalData.push_back(mSolvers.back()->AdditionalPhysicalDataIsNeeded());
        }

        // Set some member variables from the parameters
        CommonSettingsFromParameters();
    }

    /**
     * @brief Constructor with a list of solvers
     * @param rSolvers A vector of LinearSolverPointer to set.
     * @param ThisParameters The configuration parameters
     */
    FallbackLinearSolver(
        const std::vector<LinearSolverPointer>& rSolvers,
        Parameters ThisParameters = Parameters(R"({})")
        ) : mSolvers(rSolvers),
            mParameters(ThisParameters)
    {
        // Verify that linear solvers are not defined in the parameters
        KRATOS_ERROR_IF(mParameters.Has("solvers")) << "The solvers are already defined in the input parameters" << std::endl;

        // Set the default parameters
        mParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Fill the parameters with the solvers
        for (auto& p_solver : mSolvers) {
            FillParametersFromSolver(p_solver);
        }

        // Set some member variables from the parameters
        CommonSettingsFromParameters();
    }

    /// Copy constructor.
    FallbackLinearSolver(const FallbackLinearSolver& rOther) = delete;

    /// Destructor.
    ~FallbackLinearSolver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    FallbackLinearSolver& operator=(const FallbackLinearSolver& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called as few times as possible.
     * @details It creates the data structures that only depend on the connectivity of the matrix (and not on its coefficients) so that the memory can be allocated once and expensive operations can be done only when strictly  needed
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    void Initialize(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        // NOTE: We assume that initializes internal data of solvers, not system of equations
        for (auto& p_solver : mSolvers) {
            p_solver->Initialize(rA, rX, rB);
        }
    }

    /**
     * @brief This function is designed to be called every time the coefficients change in the system that is, normally at the beginning of each solve.
     * @details For example if we are implementing a direct solver, this is the place to do the factorization
     * so that then the backward substitution can be performed effectively more than once
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    void InitializeSolutionStep(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        GetCurrentSolver()->InitializeSolutionStep(rA, rX, rB);
    }

    /**
     * @brief This function is designed to be called at the end of the solve step.
     * @details for example this is the place to remove any data that we do not want to save for later
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    void FinalizeSolutionStep(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        GetCurrentSolver()->FinalizeSolutionStep(rA, rX, rB);
    }

    /**
     * @brief This function actually performs the solution work, eventually taking advantage of what was done before in the Initialize and InitializeSolutionStep functions.
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    void PerformSolutionStep(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        GetCurrentSolver()->PerformSolutionStep(rA, rX, rB);
    }

    /**
     * @brief This function is designed to clean up all internal data in the solver.
     * @details Clear is designed to leave the solver object as if newly created. After a clear a new Initialize is needed
     */
    void Clear() override
    {
        // Call all clears
        for (auto& p_solver : mSolvers) {
            p_solver->Clear();
        }

        // Clear the flags
        mRequiresAdditionalData.clear();

        // Clear the data
        mCurrentSolverIndex = 0;
    }

    /**
     * @brief Normal solve method.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    bool Solve(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        if (mResetSolverEachTry) mCurrentSolverIndex = 0;
        bool success = false;
        while (!success && mCurrentSolverIndex < mSolvers.size()) {
            success = GetCurrentSolver()->Solve(rA, rX, rB);
            // In case of failure
            if (!success) {
                // First, update the counter
                UpdateSolverIndex();

                // Call initialize methods
                InitializeSolutionStep(rA, rX, rB);

                // Provide additional data if needed
                if (mRequiresAdditionalData[mCurrentSolverIndex] && mpDoFSet != nullptr && mpModelPart != nullptr) {
                    ProvideAdditionalData(rA, rX, rB, *mpDoFSet, *mpModelPart);
                }
            }
        }
        // Reset pointers (to ensure is called only when needed)
        mpDoFSet = nullptr;
        mpModelPart = nullptr;
        return success;
    }

    /**
     * @brief Multi solve method for solving a set of linear systems with same coefficient matrix.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    bool Solve(
        SparseMatrixType& rA,
        DenseMatrixType& rX,
        DenseMatrixType& rB
        ) override
    {
        if (mResetSolverEachTry) mCurrentSolverIndex = 0;
        bool success = false;
        while (!success && mCurrentSolverIndex < mSolvers.size()) {
            success = GetCurrentSolver()->Solve(rA, rX, rB);
            // In case of failure
            if (!success) {
                // First, update the counter
                UpdateSolverIndex();

                // // Call initialize methods (NOTE: does not exist the method for dense matrices)
                // InitializeSolutionStep(rA, rX, rB);

                // // Provide additional data if needed (NOTE: does not exist the method for dense matrices)
                // if (mRequiresAdditionalData[mCurrentSolverIndex] && mpDoFSet != nullptr && mpModelPart != nullptr) {
                //     ProvideAdditionalData(rA, rX, rB, *mpDoFSet, *mpModelPart);
                // }
            }
        }
        // // Reset pointers (to ensure is called only when needed)
        // mpDoFSet = nullptr;
        // mpModelPart = nullptr;
        return success;
    }

    /**
     * @brief Eigenvalue and eigenvector solve method for derived eigensolvers
     * @param K The stiffness matrix
     * @param M The mass matrix
     * @param Eigenvalues The vector containing the eigen values
     * @param Eigenvectors The matrix containing the eigen vectors
     */
    void Solve(
        SparseMatrixType& K,
        SparseMatrixType& M,
        DenseVectorType& Eigenvalues,
        DenseMatrixType& Eigenvectors
        ) override
    {
        KRATOS_ERROR << "This is designed for linear system only, not for eigen values calculation" << std::endl;
    }

    /**
     * @brief Checks if additional physical data is needed by the solver.
     * @details Some solvers may require a minimum degree of knowledge of the structure of the matrix. 
     * For instance, when solving a mixed u-p problem, it is important to identify the row associated with v and p.
     * Another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers,
     * which require knowledge of the spatial position of the nodes associated with a given degree of freedom (DOF).
     * @return True if additional physical data is needed, false otherwise.
     */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return GetCurrentSolver()->AdditionalPhysicalDataIsNeeded();
    }

    /**
     * @brief Provides additional physical data required by the solver.
     * @details Some solvers may require a minimum degree of knowledge of the structure of the matrix.
     * For example, when solving a mixed u-p problem, it is important to identify the row associated with v and p.
     * Another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers,
     * which require knowledge of the spatial position of the nodes associated with a given degree of freedom (DOF).
     * This function provides the opportunity to provide such data if needed.
     * @param rA The sparse matrix.
     * @param rX The solution vector.
     * @param rB The right-hand side vector.
     * @param rDoFSet The set of degrees of freedom.
     * @param rModelPart The model part.
     */
    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rDoFSet,
        ModelPart& rModelPart
        ) override
    {
        GetCurrentSolver()->ProvideAdditionalData(rA, rX, rB, rDoFSet, rModelPart);

        // Prepare to provide additional data if needed in the next solver iteration
        if (mpDoFSet == nullptr && mpModelPart == nullptr) {
            mpDoFSet = &rDoFSet;
            mpModelPart = &rModelPart;
        }
    }

    /**
     * @brief Add a linear solver to the collection.
     * @details This function adds a linear solver to the collection without extending parameters.
     * @param pSolver Pointer to the linear solver to be added.
     */
    void AddSolver(LinearSolverPointer pSolver)
    {
        // Increase the solvers vector
        mSolvers.push_back(pSolver);

        // Extend the parameters
        FillParametersFromSolver(pSolver);
    }

    /**
     * @brief Add a linear solver to the collection with additional parameters.
     * @details This function adds a linear solver to the collection and extends the parameters.
     * @param ThisParameters Parameters associated with the linear solver.
     */
    void AddSolver(const Parameters ThisParameters)
    {
        // Create the solver
        auto p_solver = ConstructLinearSolverFromSettings(ThisParameters);
        mSolvers.push_back(p_solver);

        // Add the new solver parameters to the collection
        mParameters["solvers"].Append(ThisParameters);
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the Solvers list.
     * @return const std::vector<LinearSolverPointer>& Reference to the vector of solvers.
     */
    const std::vector<LinearSolverPointer>& GetSolvers() const
    {
        return mSolvers;
    }

    /**
     * @brief Set the Solvers list.
     * @param rSolvers A vector of LinearSolverPointer to set.
     */
    void SetSolvers(const std::vector<LinearSolverPointer>& rSolvers)
    {
        // Assign the solvers
        mSolvers = rSolvers;

        // Remove solvers and add again
        mParameters.RemoveValue("solvers");
        mParameters.AddEmptyArray("solvers");

        // Fill the parameters with the solvers
        for (auto& p_solver : mSolvers) {
            FillParametersFromSolver(p_solver);
        }
    }

    /**
     * @brief Check if the solver index should be reset for each try.
     * @return true If the solver index is reset for each try.
     * @return false Otherwise.
     */
    bool GetResetSolverEachTry() const
    {
        return mResetSolverEachTry;
    }

    /**
     * @brief Set whether the solver index should be reset for each try.
     * @param Reset Flag indicating whether to reset the solver index each try.
     */
    void SetResetSolverIndexEachTry(const bool Reset)
    {
        mResetSolverEachTry = Reset;
        mParameters["reset_solver_each_try"].SetBool(Reset);
    }

    /**
     * @brief Get the Parameters.
     * @return Parameters copy to the Parameters object.
     */
    Parameters GetParameters() const
    {
        return mParameters;
    }

    /**
     * @brief Get the Current Solver Index. (not mutable)
     * @return mCurrentSolverIndex The current solver index.
     */
    IndexType GetCurrentSolverIndex() const
    {
        return mCurrentSolverIndex;
    }

    /**
     * @brief Set the Current Solver Index to 0.
     */
    void ClearCurrentSolverIndex()
    {
        mCurrentSolverIndex = 0;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Get the number of solvers.
     * @details This function returns the number of solvers currently stored.
     * @return std::size_t The number of solvers.
     */
    std::size_t NumberOfSolvers() const
    {
        return mSolvers.size();
    }

    /**
     * @brief This method checks if the dimensions of the system of equations are consistent
     * @param rA The LHS of the system of equations
     * @param rX The vector containing the unknowns
     * @param rB The RHS of the system of equations
     * @return True if consistent, false otherwise
     */
    bool IsConsistent(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        return GetCurrentSolver()->IsConsistent(rA, rX, rB);
    }

    /**
     * @brief This method checks if the dimensions of the system of equations are consistent (dense matrix for RHS and unknowns version)
     * @param rA The LHS of the system of equations
     * @param rX The matrix containing the unknowns
     * @param rB The matrix containing the RHSs of the system of equations
     * @return True if consistent, false otherwise
     */
    bool IsConsistent(
        SparseMatrixType& rA,
        DenseMatrixType& rX,
        DenseMatrixType& rB
        ) override
    {
        return GetCurrentSolver()->IsConsistent(rA, rX, rB);
    }

    /**
     * @brief This method checks if the dimensions of the system of equations are not consistent
     * @param rA The LHS of the system of equations
     * @param rX The vector containing the unknowns
     * @param rB The RHS of the system of equations
     * @return False if consistent, true otherwise
     */
    bool IsNotConsistent(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        return GetCurrentSolver()->IsNotConsistent(rA, rX, rB);
    }

    /**
     * @brief This method checks if the dimensions of the system of equations are not consistent
     * @param rA The LHS of the system of equations
     * @param rX The matrix containing the unknowns
     * @param rB The matrix containing the RHSs of the system of equations
     * @return False if consistent, true otherwise
     */
    bool IsNotConsistent(
        SparseMatrixType& rA,
        DenseMatrixType& rX,
        DenseMatrixType& rB
        ) override
    {
        return GetCurrentSolver()->IsNotConsistent(rA, rX, rB);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Simple linear solver fallback";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Simple linear solver fallback";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Simple linear solver fallback data: ";
        for (auto& p_solver : mSolvers) {
            rOStream << "\nSolver: " << p_solver->Info() << "\n:";
            p_solver->PrintData(rOStream);
        }
        rOStream << "\nReset solver index each try: " << mResetSolverEachTry;
        rOStream << "\nGlobal parameters: " << mParameters;
        rOStream << "\nCurrent solver index: " << mCurrentSolverIndex << std::endl;
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

    /// The list of solvers to try
    std::vector<LinearSolverPointer> mSolvers;

    /// The list of flags to know if additional data is required
    std::vector<bool> mRequiresAdditionalData;

    /// The DoF set (for providing additional data if needed)
    ModelPart::DofsArrayType* mpDoFSet = nullptr;

    /// The model part (for providing additional data if needed)
    ModelPart* mpModelPart = nullptr;

    /// Flag to reset the solver index each try
    bool mResetSolverEachTry = false;

    /// The parameters
    Parameters mParameters;

    /// The current solver index
    IndexType mCurrentSolverIndex = 0;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Constructs a linear solver based on the provided settings.
     * @details This function is responsible for instantiating a linear solver object using the configuration specified in the `Settings` parameter. It leverages a factory pattern to select and construct the appropriate solver type based on the parameters provided. The selection process considers the solver type, tolerance, maximum iterations, and other solver-specific parameters contained within `Settings`.
     * The factory supports constructing various types of solvers, such as Conjugate Gradient, LU Decomposition, and others, depending on the `solver_type` parameter within `Settings`. It throws a `std::invalid_argument` exception if the settings specify an unsupported solver type or if required parameters for the selected solver type are missing or invalid.
     * @param Settings A `Parameters` object containing the configuration settings for the linear solver to be constructed. This includes type, tolerance, maximum iterations, and other solver-specific parameters.
     * @return Returns a `LinearSolverPointer` pointing to the newly constructed linear solver instance. If the factory implementation cannot find a matching solver type or if essential parameters are missing or invalid, it throws an exception.
     * @exception std::invalid_argument Thrown if the settings specify an unsupported solver type or if required parameters are missing or invalid.
     */
    LinearSolverPointer ConstructLinearSolverFromSettings(const Parameters Settings)
    {
        // Implementation that delegates to a factory method.
        // The factory method 'Create' is responsible for parsing 'Settings'
        // and constructing the appropriate solver.
        const auto linear_solver_factory = LinearSolverFactoryType();
        Parameters linear_solver_settings(Settings);
        const std::string linear_solver_type = linear_solver_settings["solver_type"].GetString();
        // If the solver type is a faster direct solver, try to use it
        if (linear_solver_type == "faster_direct_solver") {
            const std::vector<std::string> faster_direct_solvers = TSparseSpaceType::FastestDirectSolverList();
            for (const std::string& r_solver_name : faster_direct_solvers) {
                if (linear_solver_factory.Has(r_solver_name)) {
                    linear_solver_settings["solver_type"].SetString(r_solver_name);
                }
            }
        }
        return linear_solver_factory.Create(linear_solver_settings);
    }

    /**
     * @brief Adds solver parameters to the mParameters object.
     * @details This function generates a new set of parameters for a solver and adds it to the 'solvers' array within the mParameters object. Each new solver parameter set is given a unique identifier based on the current count of existing solvers. The function updates the 'solver_type' field with information from the provided solver, marking the settings as unknown.
     * @param pSolver A pointer to the solver from which information is extracted to fill the new parameters. This solver provides its own description through the Info() method, which is included in the 'solver_type'.
     * @throw std::invalid_argument If pSolver is nullptr, indicating an invalid solver pointer was provided.
     * @note The function assumes mParameters is properly initialized and contains a 'solvers' key that supports adding new values.
     */
    void FillParametersFromSolver(LinearSolverPointer pSolver)
    {
        // Ensure the solver pointer is not null
        KRATOS_ERROR_IF(!pSolver) << "Solver pointer is null." << std::endl;

        // Initialize dummy parameters for the new solver
        Parameters dummy_parameters = Parameters(R"({
            "solver_type": "a_name_for_the_solver"
        })");

        // Set the solver type to include the solver's information
        dummy_parameters["solver_type"].SetString(pSolver->Info() + " (settings unknown)");

        // Add the new solver parameters to the collection
        mParameters["solvers"].Append(dummy_parameters);

        // Append the boolean flag to the mRequiresAdditionalData vector
        mRequiresAdditionalData.push_back(pSolver->AdditionalPhysicalDataIsNeeded());
    }

    /**
     * @brief Updates common settings from the parameter set.
     * @details This method reads and applies common configuration settings from the `mParameters` member variable. It is designed to be called to update the object's configuration state based on external parameters, typically at the initialization stage or when parameters are updated.
     * @note This function currently updates settings related to resetting the solver index for each try. It can
     *       be extended to include more common settings as needed.
     * @see `mParameters` for the structure and expected parameters.
     */
    void CommonSettingsFromParameters()
    {
        // Set the member variables
        mResetSolverEachTry = mParameters["reset_solver_each_try"].GetBool();
    }

    /**
     * @brief Get the default parameters for this solver.
     * @details This function returns the default parameters for configuring this solver.
     * Empty in defaults. Should be filled with the solvers to try. For example:
     * {
     *     "solver_type": "amgcl"
     * },
     * {
     *     "solver_type": "skyline_lu_factorization"
     * }
     * Label of the solver is solver_x, where x is the index in the list
     * @return Default parameters for the solver.
     */
    const Parameters GetDefaultParameters() const
    {
        return Parameters(R"({
            "solver_type"           : "fallback_linear_solver",
            "solvers"               : [
                // As default is empty, if not defined it will throw an error
            ],
            "reset_solver_each_try" : false
        })");
    }

    /**
     * @brief Updates the solver index counter and logs transitions between solvers.
     * @details This method increments the current solver index to switch to the next solver in a sequence of fallback solvers.
     * It logs both the settings of the solver that failed (if any) and the settings of the next solver to be used.
     * The method is designed to facilitate easy tracking and debugging of solver transitions within a sequence.
     * @note This method now includes safety checks to prevent out-of-bounds access to the solver array. It also
     *       outlines a placeholder for enhanced future functionality, such as more comprehensive logging or
     *       additional transition actions.
     */
    void UpdateSolverIndex()
    {
        // Safety check to ensure we have solvers to work with
        KRATOS_ERROR_IF(mSolvers.empty()) << "No solvers are configured." << std::endl;

        // Log the settings of the current (failing) solver, if applicable
        if (mCurrentSolverIndex < mSolvers.size()) {
            KRATOS_INFO("FallbackLinearSolver") << "Current solver " << GetCurrentSolver()->Info() << " failed with the following settings: " << mParameters["solvers"][mCurrentSolverIndex].PrettyPrintJsonString() << std::endl;
        } else {
            KRATOS_WARNING("FallbackLinearSolver") << "Current solver index is out of bounds." << std::endl;
            return;
        }

        // Increment the counter
        mCurrentSolverIndex++;

        // Print the new solver settings
        if (mCurrentSolverIndex < mSolvers.size()) {
            KRATOS_INFO("FallbackLinearSolver") << "Switching to new solver " << GetCurrentSolver()->Info() << " with the following settings: " << mParameters["solvers"][mCurrentSolverIndex].PrettyPrintJsonString() << std::endl;
        } else {
            KRATOS_WARNING("FallbackLinearSolver") << "New solver index is out of bounds." << std::endl;
            return;
        }
    }

    /**
     * @brief Get the current linear solver.
     * @details This function returns the current linear solver being used.
     * @return A pointer to the current linear solver.
     * @note If the current solver index is out of range, an error is thrown.
     */
    LinearSolverPointer GetCurrentSolver()
    {
        // Safety check to ensure we have solvers to work with
        if (mCurrentSolverIndex < mSolvers.size()) {
            return mSolvers[mCurrentSolverIndex];
        } else { // Throw an error if the index is out of bounds
            KRATOS_ERROR << "Invalid solver index: " << mCurrentSolverIndex << std::endl;
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
}; // Class FallbackLinearSolver

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >> (std::istream& IStream,
                                  FallbackLinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FallbackLinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
