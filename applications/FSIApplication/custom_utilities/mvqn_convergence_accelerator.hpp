//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//

#if !defined(KRATOS_MVQN_CONVERGENCE_ACCELERATOR)
#define  KRATOS_MVQN_CONVERGENCE_ACCELERATOR

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"
#include "utilities/svd_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "ibqn_mvqn_convergence_accelerator.h"
#include "ibqn_mvqn_randomized_svd_convergence_accelerator.h"

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

/** @brief MVQN acceleration scheme
 * MultiVectorQuasiNewton convergence accelerator from Bogaers et al. 2016
 * @tparam TSparseSpace Linear algebra sparse space
 * @tparam TDenseSpace Linear algebra dense space
 */
template<class TSparseSpace, class TDenseSpace>
class MVQNFullJacobianConvergenceAccelerator: public ConvergenceAccelerator<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( MVQNFullJacobianConvergenceAccelerator );

    typedef ConvergenceAccelerator<TSparseSpace, TDenseSpace>              BaseType;
    typedef typename BaseType::Pointer                              BaseTypePointer;

    typedef typename BaseType::DenseVectorType                           VectorType;
    typedef typename BaseType::DenseVectorPointerType             VectorPointerType;

    typedef typename BaseType::DenseMatrixType                           MatrixType;
    typedef typename BaseType::DenseMatrixPointerType             MatrixPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new MVQNFullJacobianConvergenceAccelerator object
     * Default constructor of MVQNFullJacobianConvergenceAccelerator
     */
    explicit MVQNFullJacobianConvergenceAccelerator() = default;

    /**
     * @brief Construct a new MVQNFullJacobianConvergenceAccelerator object
     * MultiVector Quasi-Newton convergence accelerator with full Jacobian json settings constructor
     * This constructor also makes possible the MVQN usage for the interface block Newton equations
     * @param rConvAcceleratorParameters Json string encapsulating the settings
     */
    explicit MVQNFullJacobianConvergenceAccelerator(Parameters rConvAcceleratorParameters)
    {
        rConvAcceleratorParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        mOmega_0 = rConvAcceleratorParameters["w_0"].GetDouble();
        mAbsCutOff = rConvAcceleratorParameters["abs_cut_off_tol"].GetDouble();
        mUsedInBlockNewtonEquations = rConvAcceleratorParameters["interface_block_newton"].GetBool();
    }

    /**
     * @brief Construct a new MVQNFullJacobianConvergenceAccelerator object
     * This constructor intended to be used with the interface block Newton equations (@see IBQNMVQNConvergenceAccelerator)
     * @param AbsCutOff absolute tolerance for the observation cut off
     */
    explicit MVQNFullJacobianConvergenceAccelerator(const double AbsCutOff)
        : mOmega_0(0.0)
        , mAbsCutOff(AbsCutOff)
        , mUsedInBlockNewtonEquations(true)
    {
    }

    /**
     * Copy Constructor.
     */

    // Required to build the IBN-MVQN
    MVQNFullJacobianConvergenceAccelerator(const MVQNFullJacobianConvergenceAccelerator& rOther);

    /**
     * Destructor.
     */
    virtual ~MVQNFullJacobianConvergenceAccelerator(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initialize the internal iteration counter
     * This method initializes the convergence acceleratior iteration counter at the begining of the step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        mConvergenceAcceleratorIteration = 0;

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Performs the solution update
     * This method computes the solution update using a Jacobian approximation.
     * Such Jacobian approximation is computed with the MVQN (MultiVector Quasi-Newton method) algorithm.
     * @param rResidualVector Residual vector from the residual evaluation
     * @param rIterationGuess Current iteration guess to be corrected. Should be initialized to zero outside the convergence accelerator.
     */
    void UpdateSolution(
        const VectorType& rResidualVector,
        VectorType& rIterationGuess) override
    {
        KRATOS_TRY;

        // Update observation matrices and inverse Jacobian approximation
        UpdateInverseJacobianApproximation(rResidualVector, rIterationGuess);

        // Update the iteration guess values
        UpdateIterationGuess(rIterationGuess);

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Do the MVQN variables update
     * Updates the MVQN iteration variables for the next non-linear iteration
     */
    void FinalizeNonLinearIteration() override
    {
        KRATOS_TRY;

        //TODO: I THINK THAT IN HERE WE CAN CLEAR THE MATRIX M

        // Variables update
        mpIterationValue_0 = mpIterationValue_1;
        mpResidualVector_0 = mpResidualVector_1;
        mConvergenceAcceleratorIteration += 1;

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Save the current step Jacobian
     * This method saves the current step Jacobian as previous step Jacobian for the next time step iteration
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        // Update previous time step Jacobian as the last iteration Jacobian.
        // Note that it is required to check if the last iteration Jacobian exists. It exist a corner case (if the
        // very fist time step only requires the preliminary fixed point relaxation iteration to converge the
        // previous iteration Jacobian is not needed) that might set a nullptr as previous step Jacobian.
        if (mpJac_k1) {
            mpJac_n = mpJac_k1;
        }

        // Clear the observation matrices for the next step
        mpObsMatrixV = nullptr;
        mpObsMatrixW = nullptr;

        KRATOS_CATCH( "" );
    }

    virtual Parameters GetDefaultParameters() const
    {
        Parameters mvqn_default_parameters(R"({
            "solver_type"            : "MVQN",
            "w_0"                    : 0.825,
            "abs_cut_off_tol"        : 1e-8,
            "interface_block_newton" : false
        })");

        return mvqn_default_parameters;
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

    friend class IBQNMVQNConvergenceAccelerator<TSparseSpace, TDenseSpace>;
    friend class IBQNMVQNRandomizedSVDConvergenceAccelerator<TSparseSpace, TDenseSpace>;

    ///@}
protected:
    ///@name Protected  Operations
    ///@{

    /**
     * @brief Update the inverse Jacobian approximation
     * This method first appends the current iteration values to the observation matrix to then update the inverse Jacobian approximation
     * @param rResidualVector Current iteration residual vector (\tilde{u}_{k} - u_{k})
     * @param rIterationGuess Current iteration initial guess (previous iteration corrected value u_{k})
     */
    virtual void UpdateInverseJacobianApproximation(
        const VectorType& rResidualVector,
        const VectorType& rIterationGuess)
    {
        // Append and check the singularity of the current iteration information
        AppendCurrentIterationInformation(rResidualVector, rIterationGuess);

        // Update the inverse Jacobian approximation
        CalculateInverseJacobianApproximation();
    }

    /**
     * @brief Calculate the current iteration correction
     * This method calculates the current iteration correction and applies it to the provided current iteration guess
     * @param rIterationGuess Current iteration initial guess to be corrected (previous iteration corrected value u_{k})
     */
    void UpdateIterationGuess(VectorType& rIterationGuess)
    {
        if (mConvergenceAcceleratorFirstCorrectionPerformed == false) {
            // The very first correction of the problem is done with a fixed point iteration
            TSparseSpace::UnaliasedAdd(rIterationGuess, mOmega_0, *mpResidualVector_1);
            mConvergenceAcceleratorFirstCorrectionPerformed = true;
        } else {
            // Perform the correction
            // In the first iterations, the correction is performed with previous step Jacobian (stored in mpJac_k1)
            // In the subsequent iterations, the previous step Jacobian has been updated with the observation matrices
            VectorType AuxVec(mProblemSize);
            CalculateCorrectionWithJacobian(AuxVec);
            TSparseSpace::UnaliasedAdd(rIterationGuess, -1.0, AuxVec);
        }
    }

    /**
     * @brief Append the current iteration information to the observation matrices
     * This method appends the current iteration information to the observation matrices
     * First of all appends the provided information to then check the singularity of the new data
     * If the new columns are linear dependent to the existent ones, the new info is dropped
     * Nothing is appent in the first 0 iteration since the observation matrices collect the increments between iterations
     * @param rResidualVector Current iteration residual vector (\tilde{u}_{k} - u_{k})
     * @param rIterationGuess Current iteration initial guess (previous iteration corrected value u_{k})
     */
    void AppendCurrentIterationInformation(
        const VectorType& rResidualVector,
        const VectorType& rIterationGuess)
    {
        // Initialize the problem size with the first residual
        if (mProblemSize == 0 ) {
            mProblemSize = TSparseSpace::Size(rResidualVector);
        }

        // Update the current iteration and residual vector pointers
        VectorPointerType pAuxResidualVector(new VectorType(rResidualVector));
        VectorPointerType pAuxIterationGuess(new VectorType(rIterationGuess));
        std::swap(mpResidualVector_1, pAuxResidualVector);
        std::swap(mpIterationValue_1, pAuxIterationGuess);

        if (mConvergenceAcceleratorIteration != 0) {
            // Store current observation information
            StoreDataColumns();

            // Checks if the latest information is relevant. If not, the new columns are dropped
            // This also computes the trans(V)*V inverse matrix required for the Jacobian approximation
            CheckCurrentIterationInformationSingularity();
        }
    }

    /**
     * @brief Store the current iterations data columns
     * This method initializes the observation matrices for the current iteration and updates their value
     * If the number of iterations becomes larger thant the interface DOFs the oldest data is drop to prevent singularity
     */
    void StoreDataColumns()
    {
        if (mConvergenceAcceleratorIteration == 1) {
            // Resize and initalize the observation matrices
            InitializeDataColumns();
        } else {
            // Reshape the existent observation matrices
            const std::size_t n_old_cols = TDenseSpace::Size2(*mpObsMatrixV);
            if (n_old_cols < mProblemSize) {
                AppendDataColumns();
            } else {
                DropAndAppendDataColumns();
            }
        }
    }

    /**
     * @brief Check the new data singularity
     * This method checks that the new observation data is not linearly dependent to the existent one
     * To do that the ratio between the minimum and maximum eigenvalues is calculated with a small SVD
     * The matrices of the SVD are also employed to calculate the Moore-Penrose pseudoinverse (trans(V)*V)^-1
     */
    void CheckCurrentIterationInformationSingularity()
    {
        // Compute the matrix trans(V)*V
        const auto& r_V = *mpObsMatrixV;
        std::size_t data_cols = TDenseSpace::Size2(r_V);
        MatrixType transV_V(data_cols, data_cols);
        noalias(transV_V) = prod(trans(r_V), r_V);

        // Perform the auxiliary matrix trans(V)*V SVD decomposition
        MatrixType u_svd; // Orthogonal matrix (m x m)
        MatrixType w_svd; // Rectangular diagonal matrix (m x n)
        MatrixType v_svd; // Orthogonal matrix (n x n)
        std::string svd_type = "Jacobi"; // SVD decomposition type
        const double svd_rel_tol = 1.0e-6; // Relative tolerance of the SVD decomposition (it will be multiplied by the input matrix norm)
        const std::size_t max_iter = 500; // Maximum number of iterations of the Jacobi SVD decomposition
        SVDUtils<double>::SingularValueDecomposition(transV_V, u_svd, w_svd, v_svd, svd_type, svd_rel_tol, max_iter);

        // Get the maximum and minimum eigenvalues
        // Note that eigenvalues of trans(A)*A are equal to the eigenvalues of A^2
        double max_eig_V = 0.0;
        double min_eig_V = std::numeric_limits<double>::max();
        for (std::size_t i_col = 0; i_col < data_cols; ++i_col){
            const double i_eigval = std::sqrt(w_svd(i_col,i_col));
            if (max_eig_V < i_eigval) {
                max_eig_V = i_eigval;
            } else if (min_eig_V > i_eigval) {
                min_eig_V = i_eigval;
            }
        }

        // Check if the current information is relevant
        if (min_eig_V < mAbsCutOff * max_eig_V) {
            KRATOS_WARNING("MVQNFullJacobianConvergenceAccelerator")
                << "Dropping new observation columns information. Residual observation matrix min. eigval.: " << min_eig_V << " (tolerance " << mAbsCutOff * max_eig_V << ")" << std::endl;
            // Drop the observation matrices last column
            this->DropLastDataColumn();
            // Update the number of columns
            --data_cols;
            // Recompute trans(V)*V
            transV_V.resize(data_cols, data_cols, false);
            const auto& r_new_V = *(pGetResidualObservationMatrix());
            noalias(transV_V) = prod(trans(r_new_V),r_new_V);
            // Recompute the SVD for the matrix pseudo-inverse calculation
            SVDUtils<double>::SingularValueDecomposition(transV_V, u_svd, w_svd, v_svd, svd_type, svd_rel_tol);
        }

        // Calculate and save the residual observation matrices pseudo-inverse
        // This will be used later on for the inverse Jacobian approximation
        // Note that we take advantage of the fact that the matrix is always squared
        // MatrixType aux2_inv = ZeroMatrix(data_cols, data_cols);
        MatrixPointerType p_aux_inv(new MatrixType(data_cols, data_cols, 0.0));
        auto& r_aux_inv = *p_aux_inv;
        for (std::size_t i = 0; i < data_cols; ++i) {
            for (std::size_t j = 0; j < data_cols; ++j) {
                const double aux = v_svd(j,i) / w_svd(j,j);
                for (std::size_t k = 0; k < data_cols; ++k) {
                    r_aux_inv(i,k) += aux * u_svd(k,j);
                }
            }
        }
        std::swap(mpVtransVPseudoInv, p_aux_inv);
    }

    /**
     * @brief Calculates the inverse Jacobian approximation
     * This method calculates the MVQN inverse Jacobian approximation
     * If not initialized yet, the current and previous step Jacobians are initialized as required
     */
    void CalculateInverseJacobianApproximation()
    {
        if (!mJacobiansAreInitialized) {
            // Initialize the Jacobian matrices
            // This method already considers if the current MVQN accelerator is applied in a block Newton iteration
            InitializeJacobianMatrices();
            mJacobiansAreInitialized = true;
        } else {
            // Update the current Jacobian matrix if there are already observations (from iteration 1)
            // If there are no observations, the previous step Jacobian will be used to calculate the update
            UpdateCurrentJacobianMatrix();
        }
    }

    /**
     * @brief Calculates the current iteration correction with the Jacobian
     * This method calculates the current iteration correction with the inverse Jacobian approximation
     * @param rCorrection Auxiliary vector to store the correction
     */
    virtual void CalculateCorrectionWithJacobian(VectorType& rCorrection)
    {
        TDenseSpace::Mult(*mpJac_k1, *mpResidualVector_1, rCorrection);
    }

    /**
     * @brief Updates the inverse Jacobian approximation
     * This method updates the inverse Jacobian approximation with the current iteration information
     */
    virtual void UpdateCurrentJacobianMatrix()
    {
        // Compute the current inverse Jacobian approximation
        if (mpObsMatrixV != nullptr && mpObsMatrixW != nullptr) {
            // If there are already observations use these in the Jacobian update formula
            const auto& r_V = *mpObsMatrixV;
            const auto& r_W = *mpObsMatrixW;
            const std::size_t n_dofs = GetProblemSize();
            const std::size_t data_cols = TDenseSpace::Size2(r_V);
            MatrixType aux1(n_dofs, data_cols);
            MatrixType aux2(n_dofs, data_cols);
            noalias(aux1) = r_W - prod(*mpJac_n, r_V);
            noalias(aux2) = prod(aux1, *mpVtransVPseudoInv);
            MatrixPointerType p_aux_jac_k1 = MatrixPointerType(new MatrixType(*mpJac_n + prod(aux2,trans(r_V))));
            std::swap(mpJac_k1, p_aux_jac_k1);
        } else {
            // If there is no observations yet (iteration 0) use the previous step Jacobian as current one
            mpJac_k1 = mpJac_n;
        }
    }

    ///@}
    ///@name Protected  Access
    ///@{

    /**
     * @brief Get the Problem Size object
     * This method returns the problem size, which equals the number of interface DOFs
     * @return std::size_t The problem size
     */
    std::size_t GetProblemSize() const override
    {
        return mProblemSize;
    }

    std::size_t GetConvergenceAcceleratorIteration() const
    {
        return mConvergenceAcceleratorIteration;
    }

    std::size_t GetNumberOfObservations() const
    {
        return TDenseSpace::Size2(*mpObsMatrixV);
    }

    bool IsFirstCorrectionPerformed() const
    {
        return mConvergenceAcceleratorFirstCorrectionPerformed;
    }

    /**
     * @brief Returns the a bool indicating if IBQN are used
     * This method returns a bool flag indicating if the current MVQN instance is to be used in the IBQN equations
     * @return true If the MVQN instance is to be used in the IBQN equations
     * @return false If the MVQN instance is not to be used in the IBQN equations
     */
    bool IsUsedInBlockNewtonEquations() const
    {
        return mUsedInBlockNewtonEquations;
    }

    /**
     * @brief Get the residual vector
     * This method returns a pointer to the current iteration residual vector
     * @return VectorPointerType Pointer to the residual vector
     */
    VectorPointerType pGetCurrentIterationResidualVector()
    {
        return mpResidualVector_1;
    }

    /**
     * @brief Get the inverse Jacobian approximation
     * This method returns a pointer to the current inverse Jacobian approximation matrix
     * @return MatrixPointerType Pointer to the inverse Jacobian approximation matrix
     */
    MatrixPointerType pGetInverseJacobianApproximation() override
    {
        return mpJac_k1;
    }

    /**
     * @brief Get the observation matrix V
     * This method returns a pointer to the residual observation matrix
     * @return MatrixPointerType Pointer to the residual observation matrix
     */
    MatrixPointerType pGetResidualObservationMatrix()
    {
        return mpObsMatrixV;
    }

    /**
     * @brief Get the observation matrix W
     * This method returns a pointer to the solution observation matrix
     * @return MatrixPointerType Pointer to the solution observation matrix
     */
    MatrixPointerType pGetSolutionObservationMatrix()
    {
        return mpObsMatrixW;
    }

    /**
     * @brief Set the Initial Relaxation Omega
     * This method sets the relaxation parameter to be used in the very first correction
     * @param Omega Relaxation parameter to be set
     */
    void SetInitialRelaxationOmega(const double Omega)
    {
        mOmega_0 = Omega;
    }

    /**
     * @brief Set the Cut Off Tolerance
     * This method sets the cut off relative tolerance value to be used in the information singularity check
     * @param CutOffTolerance Cut off relative tolerance to be set
     */
    void SetCutOffTolerance(const double CutOffTolerance)
    {
        mAbsCutOff = CutOffTolerance;
    }

    virtual MatrixPointerType pGetJacobianDecompositionMatrixQU()
    {
        KRATOS_ERROR << "Jacobian decomposition not available for this convergence accelerator." << std::endl;
    }

    virtual MatrixPointerType pGetJacobianDecompositionMatrixSigmaV()
    {
        KRATOS_ERROR << "Jacobian decomposition not available for this convergence accelerator." << std::endl;
    }

    virtual MatrixPointerType pGetOldJacobianDecompositionMatrixQU()
    {
        KRATOS_ERROR << "Jacobian decomposition not available for this convergence accelerator." << std::endl;
    }

    virtual MatrixPointerType pGetOldJacobianDecompositionMatrixSigmaV()
    {
        KRATOS_ERROR << "Jacobian decomposition not available for this convergence accelerator." << std::endl;
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    double mOmega_0 = 0.825;                                        // Relaxation factor for the initial fixed point iteration
    //TODO: THIS IS NO LONGER ABSOLUTE --> CHANGE THE NAME
    double mAbsCutOff = 1.0e-8;                                     // Tolerance for the absolute cut-off criterion
    bool mUsedInBlockNewtonEquations = false;                       // Indicates if the current MVQN is to be used in the interface block Newton equations
    unsigned int mProblemSize = 0;                                  // Residual to minimize size
    unsigned int mConvergenceAcceleratorIteration = 0;              // Convergence accelerator iteration counter
    bool mJacobiansAreInitialized = false;                          // Indicates that the Jacobian matrices have been already initialized
    bool mConvergenceAcceleratorFirstCorrectionPerformed = false;   // Indicates that the initial fixed point iteration has been already performed

    VectorPointerType mpResidualVector_0;       // Previous iteration residual vector
    VectorPointerType mpResidualVector_1;       // Current iteration residual vector
    VectorPointerType mpIterationValue_0;       // Previous iteration guess
    VectorPointerType mpIterationValue_1;       // Current iteration guess

    MatrixPointerType mpJac_n;                  // Previous step Jacobian approximation
    MatrixPointerType mpJac_k1;                 // Current iteration Jacobian approximation
    MatrixPointerType mpObsMatrixV;             // Residual increment observation matrix
    MatrixPointerType mpObsMatrixW;             // Solution increment observation matrix
    MatrixPointerType mpVtransVPseudoInv;       // Auxiliary matrix trans(V)*V pseudo-inverse


    ///@}
    ///@name Private Operators
    ///@{

    void InitializeJacobianMatrices()
    {
        if (mUsedInBlockNewtonEquations) {
            // Initialize the previous step Jacobian to zero
            MatrixPointerType p_new_jac_n = Kratos::make_shared<MatrixType>(mProblemSize,mProblemSize);
            (*p_new_jac_n) = ZeroMatrix(mProblemSize,mProblemSize);
            std::swap(p_new_jac_n,mpJac_n);

            // Initialize the current Jacobian approximation to a zero matrix
            // Note that this is only required for the Interface Block Newton algorithm
            MatrixPointerType p_new_jac_k1 = Kratos::make_shared<MatrixType>(mProblemSize,mProblemSize);
            (*p_new_jac_k1) = ZeroMatrix(mProblemSize,mProblemSize);
            std::swap(p_new_jac_k1, mpJac_k1);
        } else {
            // Initialize the previous step Jacobian approximation to minus the diagonal matrix
            MatrixPointerType p_new_jac_n = Kratos::make_shared<MatrixType>(mProblemSize,mProblemSize);
            (*p_new_jac_n) = -1.0 * IdentityMatrix(mProblemSize,mProblemSize);
            std::swap(p_new_jac_n,mpJac_n);
        }
    }

    void InitializeDataColumns()
    {
        // Resize the observation matrices in accordance to the problem size
        MatrixPointerType p_aux_V = MatrixPointerType(new MatrixType(mProblemSize, 1));
        MatrixPointerType p_aux_W = MatrixPointerType(new MatrixType(mProblemSize, 1));
        std::swap(mpObsMatrixV, p_aux_V);
        std::swap(mpObsMatrixW, p_aux_W);

        // First observation matrices fill
        IndexPartition<unsigned int>(mProblemSize).for_each([&](unsigned int I){
            (*mpObsMatrixV)(I,0) = (*mpResidualVector_1)(I) - (*mpResidualVector_0)(I);
            (*mpObsMatrixW)(I,0) = (*mpIterationValue_1)(I) - (*mpIterationValue_0)(I);
        });
    }

    /**
     * @brief Append the new data to the observation matrices
     * This function appends the new data to both observation matrices
     */
    void AppendDataColumns()
    {
        // Create two auxilar observation matrices to reshape the existent ones
        const std::size_t n_old_cols = mpObsMatrixV->size2();
        MatrixPointerType p_new_V = Kratos::make_shared<MatrixType>(mProblemSize, n_old_cols + 1);
        MatrixPointerType p_new_W = Kratos::make_shared<MatrixType>(mProblemSize, n_old_cols + 1);

        // Recover the previous iterations information
        IndexPartition<unsigned int>(mProblemSize).for_each([&](unsigned int I){
            for (unsigned int j = 0; j < n_old_cols; j++){
                (*p_new_V)(I,j) = (*mpObsMatrixV)(I,j);
                (*p_new_W)(I,j) = (*mpObsMatrixW)(I,j);
            }
        });

        // Fill the attached column with the current iteration information
        IndexPartition<unsigned int>(mProblemSize).for_each([&](unsigned int I){
            (*p_new_V)(I, n_old_cols) = (*mpResidualVector_1)(I) - (*mpResidualVector_0)(I);
            (*p_new_W)(I, n_old_cols) = (*mpIterationValue_1)(I) - (*mpIterationValue_0)(I);
        });

        std::swap(mpObsMatrixV,p_new_V);
        std::swap(mpObsMatrixW,p_new_W);
    }

    /**
     * @brief Drop the oldest data columns to append the new data
     * This function drops the oldest information to append the new data to both observation matrices
     */
    void DropAndAppendDataColumns()
    {
        // Observation matrices size is equal to the interface DOFs number. Oldest columns are to be dropped.
        MatrixPointerType p_new_V = Kratos::make_shared<MatrixType>(mProblemSize, mProblemSize);
        MatrixPointerType p_new_W = Kratos::make_shared<MatrixType>(mProblemSize, mProblemSize);

        // Drop the oldest column and reorder data
        IndexPartition<unsigned int>(mProblemSize).for_each([&](unsigned int I){
            for (unsigned int j = 0; j < (mProblemSize-1); j++){
                (*p_new_V)(I,j) = (*mpObsMatrixV)(I,j+1);
                (*p_new_W)(I,j) = (*mpObsMatrixW)(I,j+1);
            }
        });

        // Fill the last observation matrices column
        IndexPartition<unsigned int>(mProblemSize).for_each([&](unsigned int I){
            (*p_new_V)(I, mProblemSize-1) = (*mpResidualVector_1)(I) - (*mpResidualVector_0)(I);
            (*p_new_W)(I, mProblemSize-1) = (*mpIterationValue_1)(I) - (*mpIterationValue_0)(I);
        });

        std::swap(mpObsMatrixV,p_new_V);
        std::swap(mpObsMatrixW,p_new_W);
    }

    /**
     * @brief Drop the last column of both observation matrices
     * This function drops the last column of both observation matrices
     */
    void DropLastDataColumn()
    {
        // Set two auxiliar observation matrices
        const auto n_cols = mpObsMatrixV->size2() - 1;
        MatrixPointerType p_aux_V = Kratos::make_shared<MatrixType>(mProblemSize, n_cols);
        MatrixPointerType p_aux_W = Kratos::make_shared<MatrixType>(mProblemSize, n_cols);

        // Drop the last column
        IndexPartition<std::size_t>(mProblemSize).for_each([&](unsigned int IRow){
            for (std::size_t i_col = 0; i_col < n_cols; ++i_col){
                (*p_aux_V)(IRow, i_col) = (*mpObsMatrixV)(IRow, i_col);
                (*p_aux_W)(IRow, i_col) = (*mpObsMatrixW)(IRow, i_col);
            }
        });

        // Set the member observation matrices pointers
        std::swap(mpObsMatrixV,p_aux_V);
        std::swap(mpObsMatrixW,p_aux_W);
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class MVQNFullJacobianConvergenceAccelerator */

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_MVQN_CONVERGENCE_ACCELERATOR defined */
