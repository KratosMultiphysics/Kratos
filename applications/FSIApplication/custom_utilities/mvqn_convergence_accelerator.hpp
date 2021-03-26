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
#include "utilities/qr_utility.h"
#include "utilities/svd_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "ibqn_mvqn_convergence_accelerator.h"

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

    typedef std::size_t                                                    SizeType;

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
     * MultiVector Quasi-Newton convergence accelerator with full Jacobian json settings constructor
     * This constructor also makes possible the MVQN usage for the interface block Newton equations
     * @param rConvAcceleratorParameters Json string encapsulating the settings
     */
    explicit MVQNFullJacobianConvergenceAccelerator(Parameters rConvAcceleratorParameters)
    {
        Parameters mvqn_default_parameters(R"({
            "solver_type"            : "MVQN",
            "w_0"                    : 0.825,
            "abs_cut_off_tol"        : 1e-8,
            "interface_block_newton" : false
        })");
        rConvAcceleratorParameters.ValidateAndAssignDefaults(mvqn_default_parameters);

        mOmega_0 = rConvAcceleratorParameters["w_0"].GetDouble();
        mAbsCutOff = rConvAcceleratorParameters["abs_cut_off_tol"].GetDouble();
        mUsedInBlockNewtonEquations = rConvAcceleratorParameters["interface_block_newton"].GetBool();
    }

    /**
     * @brief Construct a new MVQNFullJacobianConvergenceAccelerator object
     * This constructor also makes possible the MVQN usage for the interface block Newton equations
     * @param OmegaInitial initial relaxationo parameter
     * @param AbsCutOff absolute tolerance for the observation cut off
     * @param UsedInBlockNewtonEquations bool value indicating if the MVQN instance is to be used in the interface block Newton equations
     */
    explicit MVQNFullJacobianConvergenceAccelerator(
        const double OmegaInitial = 0.825,
        const double AbsCutOff = 1e-8,
        const bool UsedInBlockNewtonEquations = false)
        : mOmega_0(OmegaInitial)
        , mAbsCutOff(AbsCutOff)
        , mUsedInBlockNewtonEquations(UsedInBlockNewtonEquations)
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

        // Compute ((trans(V)*V)^-1)*trans(V)
        MatrixType aux_M;
        CalculateAuxiliaryOperationsV(aux_M);

        KRATOS_WATCH(aux_M.size1())
        KRATOS_WATCH(aux_M.size2())

        // Do the truncated SVD decomposition of the current Jacobian approximation to save it for the next iteration
        MatrixType y;
        MultiplyRight(aux_M, *mpOmega, y);
        KRATOS_WATCH(y.size1())
        KRATOS_WATCH(y.size2())

        //TODO: Substitute y.size2() by mModesNumber for the sake of clarity
        QR<double, row_major> qr_util;
        qr_util.compute(y.size1(), y.size2(), &(y(0,0)));
        qr_util.compute_q();
        MatrixType Q(mProblemSize, y.size2());
        for (SizeType i = 0; i < mProblemSize; ++i) { //TODO: This can be parallel
            for (SizeType j = 0; j < y.size2(); ++j) {
                Q(i,j) = qr_util.Q(i,j);
            }
        }

        KRATOS_WATCH(Q.size1())
        KRATOS_WATCH(Q.size2())

        MatrixType phi;
        MultiplyTransposeLeft(aux_M, Q, phi);

        KRATOS_WATCH(phi.size1())
        KRATOS_WATCH(phi.size2())
        KRATOS_WATCH("Before SVD")

        MatrixType u_svd; // Orthogonal matrix (m x m)
        MatrixType w_svd; // Rectangular diagonal matrix (m x n)
        MatrixType v_svd; // Orthogonal matrix (n x n)
        std::string svd_type = "Jacobi"; // SVD decomposition type
        double svd_rel_tol = 1.0e-4; // Relative tolerance of the SVD decomposition (it will be multiplied by the input matrix norm)
        const MatrixType trans_phi = trans(phi);
        SVDUtils<double>::SingularValueDecomposition(trans_phi, u_svd, w_svd, v_svd, svd_type, svd_rel_tol, 10);

        KRATOS_WATCH("After SVD")

        KRATOS_WATCH(u_svd.size1())
        KRATOS_WATCH(u_svd.size2())
        KRATOS_WATCH(w_svd.size1())
        KRATOS_WATCH(w_svd.size2())
        KRATOS_WATCH(v_svd.size1())
        KRATOS_WATCH(v_svd.size2())

        auto p_aux_QU = Kratos::make_unique<MatrixType>(prod(Q, trans(v_svd)));
        auto p_aux_SigmaV = Kratos::make_unique<MatrixType>(prod(trans(w_svd), trans(u_svd)));
        std::swap(mpOldJacQU, p_aux_QU);
        std::swap(mpOldJacSigmaV, p_aux_SigmaV);

        KRATOS_WATCH(mpOldJacQU->size1())
        KRATOS_WATCH(mpOldJacQU->size2())
        KRATOS_WATCH(mpOldJacSigmaV->size1())
        KRATOS_WATCH(mpOldJacSigmaV->size2())

        Matrix reconstructed_jacobian = prod(*mpOldJacQU, *mpOldJacSigmaV);
        double error_norm = 0.0;
        const auto& r_current_jac = *mpJac_k1;
        for (SizeType i = 0; i < reconstructed_jacobian.size1(); ++i) {
            for (SizeType j = 0; j < reconstructed_jacobian.size2(); ++j) {
                const double aux = i == j ? 1.0 + reconstructed_jacobian(i,j) : reconstructed_jacobian(i,j);
                error_norm += std::pow(aux - r_current_jac(i,j), 2);
            }
        }
        std::cout << "Jacobian error norm: " << std::sqrt(error_norm) << std::endl;

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

    friend class IBQNMVQNConvergenceAccelerator<TSparseSpace, TDenseSpace>;

    ///@}
protected:
    ///@name Protected  Operations
    ///@{

    void UpdateInverseJacobianApproximation(
        const VectorType& rResidualVector,
        const VectorType& rIterationGuess)
    {
        VectorPointerType pAuxResidualVector(new VectorType(rResidualVector));
        VectorPointerType pAuxIterationGuess(new VectorType(rIterationGuess));
        std::swap(mpResidualVector_1, pAuxResidualVector);
        std::swap(mpIterationValue_1, pAuxIterationGuess);

        if (mConvergenceAcceleratorIteration == 0) {
            if (!mJacobiansAreInitialized) {
                // Initialize the problem size with the first residual
                mProblemSize = TSparseSpace::Size(rResidualVector);

                // Initialize the Jacobian matrices
                // This method already considers if the current MVQN accelerator is applied in a block Newton iteration
                InitializeJacobianMatrices();
                mJacobiansAreInitialized = true;
            } else {
                // Set as current iteration Jacobian the previous step one
                // This is required since the first iteration needs to be performed with the previous step Jacobian
                MatrixPointerType p_aux_jac_n = Kratos::make_shared<MatrixType>(*mpJac_n);
                std::swap(p_aux_jac_n, mpJac_k1);
            }

            if (!mRandomValuesAreInitialized) {
                // If not initialized yet, create the random values matrix for the truncated SVD
                InitializeRandomValuesMatrix();
            }
        } else {
            // Store current observation information
            if (mConvergenceAcceleratorIteration == 1) {
                // Resize and initalize the observation matrices
                InitializeDataColumns();
            } else {
                // Reshape the existent observation matrices
                const std::size_t n_old_cols = mpObsMatrixV->size2();
                if (n_old_cols < mProblemSize) {
                    AppendDataColumns();
                } else {
                    DropAndAppendDataColumns();
                }
            }

            // Compute the jacobian approximation
            std::size_t data_cols = mpObsMatrixV->size2();
            MatrixType aux2(data_cols, data_cols);
            noalias(aux2) = prod(trans(*mpObsMatrixV),*mpObsMatrixV);

            // Perform the observation matrix V Singular Value Decomposition (SVD) such that
            // matrix V (m x n) is equal to the SVD matrices product V = u_svd * w_svd * v_svd
            MatrixType u_svd; // Orthogonal matrix (m x m)
            MatrixType w_svd; // Rectangular diagonal matrix (m x n)
            MatrixType v_svd; // Orthogonal matrix (n x n)
            std::string svd_type = "Jacobi"; // SVD decomposition type
            double svd_rel_tol = 1.0e-6; // Relative tolerance of the SVD decomposition (it will be multiplied by the input matrix norm)
            SVDUtils<double>::SingularValueDecomposition(aux2, u_svd, w_svd, v_svd, svd_type, svd_rel_tol);

            // Get the eigenvalues vector. Remember that eigenvalues
            // of trans(A)*A are equal to the eigenvalues of A^2
            std::vector<double> eig_vector(data_cols);
            for (std::size_t i_col = 0; i_col < data_cols; ++i_col){
                const double i_col_eig = std::sqrt(w_svd(i_col,i_col));
                eig_vector[i_col] = i_col_eig;
            }

            // Get the maximum and minimum eigenvalues
            double max_eig_V = 0.0;
            double min_eig_V = std::numeric_limits<double>::max();
            for (std::size_t i_col = 0; i_col < data_cols; ++i_col){
                if (max_eig_V < eig_vector[i_col]) {
                    max_eig_V = eig_vector[i_col];
                } else if (min_eig_V > eig_vector[i_col]) {
                    min_eig_V = eig_vector[i_col];
                }
            }

            if (min_eig_V < mAbsCutOff * max_eig_V){
                KRATOS_WARNING("MVQNFullJacobianConvergenceAccelerator")
                    << "Dropping new observation columns information. Residual observation matrix min. eigval.: " << min_eig_V << " (tolerance " << mAbsCutOff * max_eig_V << ")" << std::endl;
                // Drop the observation matrices last column
                this->DropLastDataColumn();
                // Update the number of columns
                --data_cols;
                // Recompute trans(V)*V
                aux2.resize(data_cols, data_cols);
                noalias(aux2) = prod(trans(*mpObsMatrixV),*mpObsMatrixV);
                // Recompute the SVD for the matrix pseudo-inverse calculation
                SVDUtils<double>::SingularValueDecomposition(aux2, u_svd, w_svd, v_svd, svd_type, svd_rel_tol);
            }

            // Compute the matrix pseudo-inverse
            // Note that we take advantage of the fact that the matrix is always squared
            MatrixType aux2_inv = ZeroMatrix(data_cols, data_cols);
            for (std::size_t i = 0; i < data_cols; ++i) {
                for (std::size_t j = 0; j < data_cols; ++j) {
                    const double aux = v_svd(j,i) / w_svd(j,j);
                    for (std::size_t k = 0; k < data_cols; ++k) {
                        aux2_inv(i,k) += aux * u_svd(k,j);
                    }
                }
            }

            // Compute the current inverse Jacobian approximation
            MatrixType aux1(mProblemSize, data_cols);
            MatrixType aux3(mProblemSize, data_cols);
            noalias(aux1) = *mpObsMatrixW - prod(*mpJac_n,*mpObsMatrixV);
            noalias(aux3) = prod(aux1,aux2_inv);

            MatrixPointerType p_aux_jac_k1 = MatrixPointerType(new MatrixType(*mpJac_n + prod(aux3,trans(*mpObsMatrixV))));
            std::swap(mpJac_k1, p_aux_jac_k1);
        }
    }

    void UpdateIterationGuess(VectorType& rIterationGuess)
    {
        if (mConvergenceAcceleratorFirstCorrectionPerformed == false) {
            // The very first correction of the problem is done with a fixed point iteration
            TSparseSpace::UnaliasedAdd(rIterationGuess, mOmega_0, *mpResidualVector_1);
            mConvergenceAcceleratorFirstCorrectionPerformed = true;
        } else {
            // Perform the correction
            // In the first iterations, the correction is performed with previous step Jacobian (stored in mpJac_k1)
            // In the subsequent iterations, the previous step Jacobian is updated with the observation matrices
            VectorType AuxVec(mProblemSize);
            TSparseSpace::Mult(*mpJac_k1, *mpResidualVector_1, AuxVec);
            TSparseSpace::UnaliasedAdd(rIterationGuess, -1.0, AuxVec);
        }
    }

    ///@}
    ///@name Protected  Access
    ///@{

    std::size_t GetProblemSize() const override
    {
        return mProblemSize;
    }

    MatrixPointerType pGetInverseJacobianApproximation() override
    {
        return mpJac_k1;
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    double mOmega_0;                                                // Relaxation factor for the initial fixed point iteration
    double mAbsCutOff;                                              // Tolerance for the absolute cut-off criterion
    bool mUsedInBlockNewtonEquations;                               // Indicates if the current MVQN is to be used in the interface block Newton equations
    unsigned int mProblemSize = 0;                                  // Residual to minimize size
    unsigned int mConvergenceAcceleratorIteration = 0;              // Convergence accelerator iteration counter
    bool mJacobiansAreInitialized = false;                          // Indicates that the Jacobian matrices have been already initialized
    bool mRandomValuesAreInitialized = false;                       // Indicates if the random values for the truncated SVD have been already set
    bool mConvergenceAcceleratorFirstCorrectionPerformed = false;   // Indicates that the initial fixed point iteration has been already performed

    VectorPointerType mpResidualVector_0;       // Previous iteration residual vector
    VectorPointerType mpResidualVector_1;       // Current iteration residual vector
    VectorPointerType mpIterationValue_0;       // Previous iteration guess
    VectorPointerType mpIterationValue_1;       // Current iteration guess

    MatrixPointerType mpJac_n;                  // Previous step Jacobian approximation
    MatrixPointerType mpJac_k1;                 // Current iteration Jacobian approximation
    MatrixPointerType mpObsMatrixV;             // Residual increment observation matrix
    MatrixPointerType mpObsMatrixW;             // Solution increment observation matrix

    Kratos::unique_ptr<MatrixType> mpOmega;     // Matrix with random values for truncated SVD
    Kratos::unique_ptr<MatrixType> mpOldJacQU = nullptr;
    Kratos::unique_ptr<MatrixType> mpOldJacSigmaV = nullptr;


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

    void InitializeRandomValuesMatrix()
    {
        // Set the random values matrix pointer
        const SizeType n_modes = 10; //TODO: MAKE THIS USER-DEFINABLE
        auto p_aux_omega = Kratos::make_unique<MatrixType>(mProblemSize, n_modes);
        std::swap(p_aux_omega, mpOmega);

        // Create the random values generator
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> distribution(0.0, 1.0);

        // Fill the random values matrix
        auto& r_omega_matrix = *mpOmega;
        for (SizeType i = 0; i < mProblemSize; ++i) {
            for (SizeType j = 0; j < n_modes; ++j) {
                // Use distribution to transform the random unsigned int generated by generator into a
                // double in [0.0, 1.0). Each call to distribution(generator) generates a new random double
                r_omega_matrix(i,j) = distribution(generator);
            }
        }
    }

    void MultiplyRight(
        const Matrix& rAuxM,
        const Matrix& rRightMatrix,
        Matrix& rSolution)
    {
        KRATOS_ERROR_IF(rRightMatrix.size1() != mProblemSize) << "Obtained right multiplication matrix size " << rRightMatrix.size1() << " does not match the problem size " << mProblemSize << " expected one." << std::endl;

        // Add to the solution matrix
        MatrixType aux_prod_omega = prod(rAuxM, rRightMatrix);
        rSolution = prod(*mpObsMatrixW, aux_prod_omega);
        if (mpOldJacQU == nullptr && mpOldJacSigmaV == nullptr) {
            for (SizeType i = 0; i < rSolution.size1(); ++i) { //TODO: DO THIS PARALLEL
                for (SizeType j = 0; j < rSolution.size2(); ++j) {
                    rSolution(i,j) -= rRightMatrix(i,j);
                }
            }
        } else {
            KRATOS_WATCH(mpOldJacQU->size1())
            KRATOS_WATCH(mpOldJacQU->size2())
            KRATOS_WATCH(mpOldJacSigmaV->size1())
            KRATOS_WATCH(mpOldJacSigmaV->size2())
            const auto& r_V = *mpObsMatrixV;
            MatrixType V_aux_prod_omega = prod(r_V, aux_prod_omega);
            MatrixType aux_old_jac_1 = prod(*mpOldJacSigmaV, rRightMatrix);
            MatrixType aux_old_jac_2 = prod(*mpOldJacSigmaV, V_aux_prod_omega);
            MatrixType aux_1 = prod(*mpOldJacQU, aux_old_jac_1);
            MatrixType aux_2 = prod(*mpOldJacQU, aux_old_jac_2);
            for (SizeType i = 0; i < rSolution.size1(); ++i) { //TODO: DO THIS PARALLEL
                for (SizeType j = 0; j < rSolution.size2(); ++j) {
                    rSolution(i,j) += aux_1(i,j) - aux_2(i,j) - V_aux_prod_omega(i,j);
                }
            }
        }
    }

    void MultiplyTransposeLeft(
        const Matrix& rAuxM,
        const Matrix& rLeftMatrix,
        Matrix& rSolution)
    {
        KRATOS_ERROR_IF(rLeftMatrix.size1() != mProblemSize) << "Obtained left multiplication matrix size " << rLeftMatrix.size1() << " does not match the problem size " << mProblemSize << " expected one." << std::endl;

        // Compute first (common) term
        //TODO: DO THIS PARALLEL
        const auto& r_W = *mpObsMatrixW;
        const SizeType n_data_cols = r_W.size2();
        const SizeType n_modes = rLeftMatrix.size2();
        rSolution = ZeroMatrix(n_modes, mProblemSize);
        for (SizeType i = 0; i < n_modes; ++i) {
            for (SizeType j = 0; j < mProblemSize; ++j) {
                double& r_solution_ij = rSolution(i,j);
                for (SizeType k = 0; k < mProblemSize; ++k) {
                    const double aux_trans_Q = rLeftMatrix(k,i);
                    for (SizeType m = 0; m < n_data_cols; ++m) {
                        r_solution_ij += aux_trans_Q * r_W(k,m) * rAuxM(m,j);
                    }
                }
            }
        }

        // Add to the solution matrix
        if (mpOldJacQU == nullptr && mpOldJacSigmaV == nullptr) {
            for (SizeType i = 0; i < rSolution.size1(); ++i) { //TODO: DO THIS PARALLEL
                for (SizeType j = 0; j < rSolution.size2(); ++j) {
                    rSolution(i,j) -= rLeftMatrix(j,i);
                }
            }
        } else {
            const auto& r_V = *mpObsMatrixV;
            MatrixType aux_left_V = prod(trans(rLeftMatrix), r_V);
            MatrixType aux_left_V_M = prod(aux_left_V, rAuxM);
            MatrixType aux_left_jac_1 = prod(trans(rLeftMatrix), *mpOldJacQU);
            MatrixType aux_left_jac_2 = prod(aux_left_jac_1, *mpOldJacSigmaV);
            //TODO: 3 AND 4 CAN BE DONE IN A UNIQUE NESTED LOOP
            MatrixType aux_left_jac_3 = prod(aux_left_jac_2, r_V);
            MatrixType aux_left_jac_4 = prod(aux_left_jac_3, rAuxM);
            for (SizeType i = 0; i < rSolution.size1(); ++i) { //TODO: DO THIS PARALLEL
                for (SizeType j = 0; j < rSolution.size2(); ++j) {
                    rSolution(i,j) += aux_left_jac_2(i,j) - aux_left_V_M(i,j) - aux_left_jac_4(i,j);
                }
            }
        }
    }

    void CalculateAuxiliaryOperationsV(MatrixType& rAuxM)
    {
        // Compute (trans(V)*V)^-1
        const auto& r_V = *mpObsMatrixV;
        const SizeType n_data_cols = r_V.size2();
        MatrixType transV_V(n_data_cols, n_data_cols);
        noalias(transV_V) = prod(trans(r_V), r_V);

        // Perform the observation matrix V Singular Value Decomposition (SVD) such that
        // matrix V (m x n) is equal to the SVD matrices product V = u_svd * w_svd * v_svd
        MatrixType u_svd; // Orthogonal matrix (m x m)
        MatrixType w_svd; // Rectangular diagonal matrix (m x n)
        MatrixType v_svd; // Orthogonal matrix (n x n)
        std::string svd_type = "Jacobi"; // SVD decomposition type
        double svd_rel_tol = 1.0e-6; // Relative tolerance of the SVD decomposition (it will be multiplied by the input matrix norm)
        SVDUtils<double>::SingularValueDecomposition(transV_V, u_svd, w_svd, v_svd, svd_type, svd_rel_tol);

        // Compute the matrix pseudo-inverse
        // Note that we take advantage of the fact that the matrix is always squared
        MatrixType transV_V_inv = ZeroMatrix(n_data_cols, n_data_cols);
        for (std::size_t i = 0; i < n_data_cols; ++i) {
            for (std::size_t j = 0; j < n_data_cols; ++j) {
                const double aux = v_svd(j,i) / w_svd(j,j);
                for (std::size_t k = 0; k < n_data_cols; ++k) {
                    transV_V_inv(i,k) += aux * u_svd(k,j);
                }
            }
        }

        //TODO: PARALLELIZE THIS
        // Compute ((trans(V)*V)^-1)*trans(V)
        rAuxM = ZeroMatrix(n_data_cols, mProblemSize);
        for (SizeType i = 0; i < n_data_cols; ++i) {
            for (SizeType j = 0; j < mProblemSize; ++j) {
                for (SizeType k = 0; k < n_data_cols; ++k) {
                    rAuxM(i,j) += transV_V_inv(i,k) * r_V(j,k);
                }
            }
        }
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
