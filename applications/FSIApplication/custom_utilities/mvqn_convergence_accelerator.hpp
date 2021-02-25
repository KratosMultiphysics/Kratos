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
#include "utilities/math_utils.h"

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

        // VectorPointerType pAuxResidualVector(new VectorType(rResidualVector));
        // VectorPointerType pAuxIterationGuess(new VectorType(rIterationGuess));
        // std::swap(mpResidualVector_1, pAuxResidualVector);
        // std::swap(mpIterationValue_1, pAuxIterationGuess);

        // if (mConvergenceAcceleratorIteration == 0) {
        //     if (mConvergenceAcceleratorFirstCorrectionPerformed == false) {
        //         // The very first correction of the problem is done with a fixed point iteration
        //         TSparseSpace::UnaliasedAdd(rIterationGuess, mOmega_0, *mpResidualVector_1);

        //         // Initialize the problem size with the first residual 
        //         mProblemSize = TSparseSpace::Size(rResidualVector);

        //         // Initialize the Jacobian approximation matrix as minus the diagonal matrix
        //         // Note that this is exclusively done in the very fist iteration
        //         MatrixPointerType p_new_jac_n = Kratos::make_shared<MatrixType>(mProblemSize,mProblemSize);
        //         (*p_new_jac_n) = -1.0 * IdentityMatrix(mProblemSize,mProblemSize);
        //         std::swap(p_new_jac_n,mpJac_n);

        //         mConvergenceAcceleratorFirstCorrectionPerformed = true;
        //     } else {
        //         // Fist step correction is done with the previous step Jacobian
        //         VectorType AuxVec(mProblemSize);
        //         TSparseSpace::Mult(*mpJac_n, *mpResidualVector_1, AuxVec);
        //         TSparseSpace::UnaliasedAdd(rIterationGuess, -1.0, AuxVec);
        //     }
        // } else {
        //     if (mConvergenceAcceleratorIteration == 1) {
        //         // Resize the observation matrices in accordance to the problem size
        //         MatrixPointerType pNewObsMatrixV = MatrixPointerType(new MatrixType(mProblemSize, 1));
        //         MatrixPointerType pNewObsMatrixW = MatrixPointerType(new MatrixType(mProblemSize, 1));

        //         std::swap(mpObsMatrixV,pNewObsMatrixV);
        //         std::swap(mpObsMatrixW,pNewObsMatrixW);

        //         // First observation matrices fill
        //         for (unsigned int i = 0; i < mProblemSize; i++)
        //         {
        //             (*mpObsMatrixV)(i,0) = (*mpResidualVector_1)(i) - (*mpResidualVector_0)(i);
        //             (*mpObsMatrixW)(i,0) = (*mpIterationValue_1)(i) - (*mpIterationValue_0)(i);
        //         }
        //     } else {
        //         // Reshape the existent observation matrices
        //         const std::size_t n_old_cols = mpObsMatrixV->size2();
        //         if (n_old_cols < mProblemSize){
        //             this->AppendDataColumns();
        //         } else {
        //             this->DropAndAppendDataColumns();
        //         }

        //     }

        //     // Compute the jacobian approximation
        //     std::size_t data_cols = mpObsMatrixV->size2();
        //     DenseMatrixType aux2(data_cols, data_cols);
        //     noalias(aux2) = prod(trans(*mpObsMatrixV),*mpObsMatrixV);

        //     // Perform the observation matrix V Singular Value Decomposition (SVD) such that
        //     // matrix V (m x n) is equal to the SVD matrices product V = u_svd * w_svd * v_svd
        //     DenseMatrixType u_svd; // Orthogonal matrix (m x m)
        //     DenseMatrixType w_svd; // Rectangular diagonal matrix (m x n)
        //     DenseMatrixType v_svd; // Orthogonal matrix (n x n)
        //     std::string svd_type = "Jacobi"; // SVD decomposition type
        //     double svd_rel_tol = 1.0e-6; // Relative tolerance of the SVD decomposition (it will be multiplied by the input matrix norm)
        //     SVDUtils<double>::SingularValueDecomposition(aux2, u_svd, w_svd, v_svd, svd_type, svd_rel_tol);

        //     // Get the eigenvalues vector. Remember that eigenvalues
        //     // of trans(A)*A are equal to the eigenvalues of A^2
        //     std::vector<double> eig_vector(data_cols);
        //     for (std::size_t i_col = 0; i_col < data_cols; ++i_col){
        //         const double i_col_eig = std::sqrt(w_svd(i_col,i_col));
        //         eig_vector[i_col] = i_col_eig;
        //     }

        //     // Get the maximum and minimum eigenvalues
        //     double max_eig_V = 0.0;
        //     double min_eig_V = std::numeric_limits<double>::max();
        //     for (std::size_t i_col = 0; i_col < data_cols; ++i_col){
        //         if (max_eig_V < eig_vector[i_col]) {
        //             max_eig_V = eig_vector[i_col];
        //         } else if (min_eig_V > eig_vector[i_col]) {
        //             min_eig_V = eig_vector[i_col];
        //         }
        //     }

        //     if (min_eig_V < mAbsCutOff * max_eig_V){
        //         KRATOS_WARNING("MVQNFullJacobianConvergenceAccelerator")
        //             << "Dropping info for new observation columns. Residual observation matrix min. eigval.: " << min_eig_V << " (tolerance " << mAbsCutOff * max_eig_V << " )" << std::endl;
        //         // Drop the observation matrices last column
        //         this->DropLastDataColumn();
        //         // Update the number of columns
        //         --data_cols;
        //         // Recompute trans(V)*V
        //         aux2.resize(data_cols, data_cols);
        //         noalias(aux2) = prod(trans(*mpObsMatrixV),*mpObsMatrixV);
        //     }

        //     // Perform the matrix inversion
        //     double det_aux2;
        //     MatrixType aux2inv(data_cols, data_cols);
        //     MathUtils<double>::InvertMatrix(aux2, aux2inv, det_aux2, -1.0e-15);
        //     // KRATOS_WARNING_IF("MVQNFullJacobianConvergenceAccelerator", std::pow(max_eig_V,2) / std::pow(min_eig_V,2) < std::pow(mAbsCutOff,2))
        //     //     << "Inverted matrix determinant is close to be singular" << std::endl;

        //     // Compute the current inverse Jacobian approximation
        //     MatrixType aux1(mProblemSize, data_cols);
        //     MatrixType aux3(mProblemSize, data_cols);
        //     noalias(aux1) = *mpObsMatrixW - prod(*mpJac_n,*mpObsMatrixV);
        //     noalias(aux3) = prod(aux1,aux2inv);

        //     MatrixPointerType pJac_k1 = MatrixPointerType(new MatrixType(mProblemSize, mProblemSize));
        //     std::swap(mpJac_k1,pJac_k1);
        //     noalias(*mpJac_k1) = *mpJac_n + prod(aux3,trans(*mpObsMatrixV));

        //     // Perform the correction
        //     VectorType AuxVec(mProblemSize);
        //     TSparseSpace::Mult(*mpJac_k1, *mpResidualVector_1, AuxVec);
        //     TSparseSpace::UnaliasedAdd(rIterationGuess, -1.0, AuxVec);
        // }

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

                // // Initialize the previous step Jacobian approximation to minus the diagonal matrix
                // // Note that this is exclusively done in the very fist iteration
                // MatrixPointerType p_new_jac_n = Kratos::make_shared<MatrixType>(mProblemSize,mProblemSize);
                // (*p_new_jac_n) = -1.0 * IdentityMatrix(mProblemSize,mProblemSize);
                // std::swap(p_new_jac_n,mpJac_n);

                // // Initialize the current Jacobian approximation to a zero matrix
                // // This is required in for the Interface Block Newton algorithm
                // MatrixPointerType p_new_jac_k1 = Kratos::make_shared<MatrixType>(mProblemSize,mProblemSize);
                // (*p_new_jac_k1) = ZeroMatrix(mProblemSize,mProblemSize);
                // std::swap(p_new_jac_k1, mpJac_k1);
            } else {
                // Set as current iteration Jacobian the previous step one
                // This is required since the first iteration needs to be performed with the previous step Jacobian
                MatrixPointerType p_aux_jac_n = Kratos::make_shared<MatrixType>(*mpJac_n);
                std::swap(p_aux_jac_n, mpJac_k1);
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
                    << "Dropping info for new observation columns. Residual observation matrix min. eigval.: " << min_eig_V << " (tolerance " << mAbsCutOff * max_eig_V << " )" << std::endl;
                // Drop the observation matrices last column
                this->DropLastDataColumn();
                // Update the number of columns
                --data_cols;
                // Recompute trans(V)*V
                aux2.resize(data_cols, data_cols);
                noalias(aux2) = prod(trans(*mpObsMatrixV),*mpObsMatrixV);
            }

            // Perform the matrix inversion
            double det_aux2;
            MatrixType aux2_inv(data_cols, data_cols);
            MathUtils<double>::InvertMatrix(aux2, aux2_inv, det_aux2, -1.0e-15);
            // KRATOS_WARNING_IF("MVQNFullJacobianConvergenceAccelerator", std::pow(max_eig_V,2) / std::pow(min_eig_V,2) < std::pow(mAbsCutOff,2))
            //     << "Inverted matrix determinant is close to be singular" << std::endl;

            // Compute the current inverse Jacobian approximation
            MatrixType aux1(mProblemSize, data_cols);
            MatrixType aux3(mProblemSize, data_cols);
            noalias(aux1) = *mpObsMatrixW - prod(*mpJac_n,*mpObsMatrixV);
            noalias(aux3) = prod(aux1,aux2_inv);

            MatrixPointerType p_aux_jac_k1 = MatrixPointerType(new MatrixType(*mpJac_n + prod(aux3,trans(*mpObsMatrixV))));
            std::swap(mpJac_k1, p_aux_jac_k1);
            // MatrixPointerType pJac_k1 = MatrixPointerType(new MatrixType(mProblemSize, mProblemSize));
            // std::swap(mpJac_k1,pJac_k1);
            // noalias(*mpJac_k1) = *mpJac_n + prod(aux3,trans(*mpObsMatrixV));
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
    
        // if (mConvergenceAcceleratorIteration == 0) {
        //     if (mConvergenceAcceleratorFirstCorrectionPerformed == false) {
        //         // The very first correction of the problem is done with a fixed point iteration
        //         TSparseSpace::UnaliasedAdd(rIterationGuess, mOmega_0, *mpResidualVector_1);
        //         mConvergenceAcceleratorFirstCorrectionPerformed = true;
        //     } else {
        //         // Fist step correction is done with the previous step Jacobian
        //         VectorType AuxVec(mProblemSize);
        //         TSparseSpace::Mult(*mpJac_n, *mpResidualVector_1, AuxVec);
        //         TSparseSpace::UnaliasedAdd(rIterationGuess, -1.0, AuxVec);
        //     }
        // } else {
        //     // Perform the correction
        //     VectorType AuxVec(mProblemSize);
        //     TSparseSpace::Mult(*mpJac_k1, *mpResidualVector_1, AuxVec);
        //     TSparseSpace::UnaliasedAdd(rIterationGuess, -1.0, AuxVec);
        // }
    }

    ///@}
    ///@name Protected  Access
    ///@{

    unsigned int GetProblemSize() override
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
    bool mConvergenceAcceleratorFirstCorrectionPerformed = false;   // Indicates that the initial fixed point iteration has been already performed

    VectorPointerType mpResidualVector_0;       // Previous iteration residual vector
    VectorPointerType mpResidualVector_1;       // Current iteration residual vector
    VectorPointerType mpIterationValue_0;       // Previous iteration guess
    VectorPointerType mpIterationValue_1;       // Current iteration guess

    MatrixPointerType mpJac_n;                  // Previous step Jacobian approximation
    MatrixPointerType mpJac_k1;                 // Current iteration Jacobian approximation
    MatrixPointerType mpObsMatrixV;             // Residual increment observation matrix
    MatrixPointerType mpObsMatrixW;             // Solution increment observation matrix


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
        //TODO: I think this can be parallel
        for (unsigned int i = 0; i < mProblemSize; i++) {
            (*mpObsMatrixV)(i,0) = (*mpResidualVector_1)(i) - (*mpResidualVector_0)(i);
            (*mpObsMatrixW)(i,0) = (*mpIterationValue_1)(i) - (*mpIterationValue_0)(i);
        }
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
        for (unsigned int i = 0; i < mProblemSize; i++){
            for (unsigned int j = 0; j < n_old_cols; j++){
                (*p_new_V)(i,j) = (*mpObsMatrixV)(i,j);
                (*p_new_W)(i,j) = (*mpObsMatrixW)(i,j);
            }
        }

        // Fill the attached column with the current iteration information
        for (unsigned int i = 0; i < mProblemSize; i++){
            (*p_new_V)(i, n_old_cols) = (*mpResidualVector_1)(i) - (*mpResidualVector_0)(i);
            (*p_new_W)(i, n_old_cols) = (*mpIterationValue_1)(i) - (*mpIterationValue_0)(i);
        }

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
        for (unsigned int i = 0; i < mProblemSize; i++){
            for (unsigned int j = 0; j < (mProblemSize-1); j++){
                (*p_new_V)(i,j) = (*mpObsMatrixV)(i,j+1);
                (*p_new_W)(i,j) = (*mpObsMatrixW)(i,j+1);
            }
        }

        // Fill the last observation matrices column
        for (unsigned int i = 0; i < mProblemSize; i++){
            (*p_new_V)(i, mProblemSize-1) = (*mpResidualVector_1)(i) - (*mpResidualVector_0)(i);
            (*p_new_W)(i, mProblemSize-1) = (*mpIterationValue_1)(i) - (*mpIterationValue_0)(i);
        }

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
        for (std::size_t i_row = 0; i_row < mProblemSize; ++i_row){
            for (std::size_t i_col = 0; i_col < n_cols; ++i_col){
                (*p_aux_V)(i_row, i_col) = (*mpObsMatrixV)(i_row, i_col);
                (*p_aux_W)(i_row, i_col) = (*mpObsMatrixW)(i_row, i_col);
            }
        }

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
