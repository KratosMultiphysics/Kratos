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

#if !defined(KRATOS_MVQN_RANDOMIZED_SVD_CONVERGENCE_ACCELERATOR)
#define  KRATOS_MVQN_RANDOMIZED_SVD_CONVERGENCE_ACCELERATOR

// System includes
#include <random>

// External includes

// Project includes
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"
#include "utilities/dense_qr_decomposition.h"
#include "utilities/dense_svd_decomposition.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "mvqn_convergence_accelerator.hpp"
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
class MVQNRandomizedSVDConvergenceAccelerator: public MVQNFullJacobianConvergenceAccelerator<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( MVQNRandomizedSVDConvergenceAccelerator );

    typedef std::size_t SizeType;

    typedef MVQNFullJacobianConvergenceAccelerator<TSparseSpace, TDenseSpace> BaseType;
    typedef typename BaseType::Pointer BaseTypePointer;

    typedef typename BaseType::DenseVectorType VectorType;
    typedef typename BaseType::DenseVectorPointerType VectorPointerType;

    typedef typename BaseType::DenseMatrixType MatrixType;
    typedef typename BaseType::DenseMatrixPointerType MatrixPointerType;

    typedef typename DenseQRDecomposition<TDenseSpace>::Pointer DenseQRPointerType;
    typedef typename DenseSingularValueDecomposition<TDenseSpace>::Pointer DenseSVDPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new MVQNRandomizedSVDConvergenceAccelerator object
     * MultiVector Quasi-Newton convergence accelerator with randomized SVD Jacobian json settings constructor
     * @param pDenseSVD Pointer to the dense SVD utility instance
     * @param ConvAcceleratorParameters Json string encapsulating the settings
     */
    explicit MVQNRandomizedSVDConvergenceAccelerator(
        DenseQRPointerType pDenseQR,
        DenseSVDPointerType pDenseSVD,
        Parameters ConvAcceleratorParameters)
    : BaseType()
    , mpDenseQR(pDenseQR)
    , mpDenseSVD(pDenseSVD)
    {
        ConvAcceleratorParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        mUserNumberOfModes = ConvAcceleratorParameters["jacobian_modes"].GetInt();
        mAutomaticJacobianModes = ConvAcceleratorParameters["automatic_jacobian_modes"].GetBool();
        mLimitModesToIterations = ConvAcceleratorParameters["limit_modes_to_iterations"].GetBool(); //TODO: THIS FIELD MUST BE RENAMED TO SOMETHING LIKE "restart_random_matrix"
        mMinRandSVDExtraModes = ConvAcceleratorParameters["min_rand_svd_extra_modes"].GetInt();
        BaseType::SetInitialRelaxationOmega(ConvAcceleratorParameters["w_0"].GetDouble());
        BaseType::SetCutOffTolerance(ConvAcceleratorParameters["cut_off_tol"].GetDouble());

        KRATOS_WARNING_IF("MVQNRandomizedSVDConvergenceAccelerator", mAutomaticJacobianModes && mUserNumberOfModes != 0)
            << "Custom and automatic number of modes have been selected. Automatic will be used." << std::endl;
    }

    /**
     * Copy Constructor.
     */

    // Required to build the IBQN-MVQN
    MVQNRandomizedSVDConvergenceAccelerator(const MVQNRandomizedSVDConvergenceAccelerator& rOther) = default;

    /**
     * Destructor.
     */
    virtual ~MVQNRandomizedSVDConvergenceAccelerator(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        // If required, calculate the SVD of the inverse Jacobian approximation
        // It is important to check if the residual minimization has been done to avoid throwing errors in the first do nothing steps
        // Note that in the IBQN case we avoid doing the decomposition as we reuse the last iteration one
        if (!BaseType::IsUsedInBlockNewtonEquations() && BaseType::IsFirstCorrectionPerformed()) {
            CalculateInverseJacobianRandomizedSVD();
        }

        // Save the current (last) inverse Jacobian decomposition data for the next step
        // Note that the current decomposition is checked in case the current step converges in one iteration with the old Jacobian
        if (mpJacQU != nullptr && mpJacSigmaV != nullptr) {
            mpOldJacQU = mpJacQU;
            mpOldJacSigmaV = mpJacSigmaV;
            std::size_t n_obs = BaseType::GetNumberOfObservations();
            if (n_obs > mOldMaxRank) {
                mOldMaxRank = n_obs;
            }
        }

        // Call the base class FinalizeSolutionStep
        // Note that this clears the observation matrices so it needs to be called after the inverse Jacobian SVD
        BaseType::FinalizeSolutionStep();

        // Clear the current step Jacobian decomposition pointers
        mpJacQU = nullptr;
        mpJacSigmaV = nullptr;

        // Update the seed for the next time step
        mSeed++;

        KRATOS_CATCH( "" );
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters mvqn_randomized_svd_default_parameters(R"({
            "solver_type"               : "MVQN_randomized_SVD",
            "automatic_jacobian_modes"  : true,
            "jacobian_modes"            : 0,
            "w_0"                       : 0.825,
            "cut_off_tol"               : 1e-8,
            "interface_block_newton"    : false,
            "limit_modes_to_iterations" : false,
            "min_rand_svd_extra_modes"  : 10
        })");

        return mvqn_randomized_svd_default_parameters;
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

    friend class IBQNMVQNRandomizedSVDConvergenceAccelerator<TSparseSpace, TDenseSpace>;

    ///@}
protected:
    ///@name Protected  Life Cycle
    ///@{

    /**
     * @brief Construct a new MVQNRandomizedSVDConvergenceAccelerator object
     * Multivector Quasi-Newton convergence accelerator with randomized SVD jacobian constructor
     * This constructor is intended to be used to set the subdomain MVQN pointers with the IBQN equations
     * @param pDenseQR Pointer to the dense QR utility instance
     * @param pDenseSVD Pointer to the dense SVD utility instance
     * @param JacobianModes Number of modes to be kept in the randomization of the Jacobian
     * @param CutOffTolerance Tolerance for the new observation data cut off (relative to the maximum residual observation matrix eigenvalue)
     * @param LimitModesToIterations If true limits the randomized SVD decomposition to the current number of iterations
     * @param MinRandSVDExtraModes Minimum number of extra modes to be used in the randomized SVD decomposition
     */
    explicit MVQNRandomizedSVDConvergenceAccelerator(
        DenseQRPointerType pDenseQR,
        DenseSVDPointerType pDenseSVD,
        const bool AutomaticJacobianModes,
        const unsigned int JacobianModes,
        const double CutOffTolerance,
        const bool LimitModesToIterations,
        const unsigned int MinRandSVDExtraModes)
    : BaseType(CutOffTolerance)
    , mUserNumberOfModes(JacobianModes)
    , mMinRandSVDExtraModes(MinRandSVDExtraModes)
    , mAutomaticJacobianModes(AutomaticJacobianModes)
    , mLimitModesToIterations(LimitModesToIterations)
    , mpDenseQR(pDenseQR)
    , mpDenseSVD(pDenseSVD)
    {
        KRATOS_WARNING_IF("MVQNRandomizedSVDConvergenceAccelerator", mAutomaticJacobianModes && mUserNumberOfModes != 0)
            << "Custom and automatic number of modes have been selected. Automatic will be used." << std::endl;
    }

    ///@}
    ///@name Protected  Operations
    ///@{

    void UpdateInverseJacobianApproximation(
        const VectorType& rResidualVector,
        const VectorType& rIterationGuess) override
    {
        // Append the current iteration information to the observation matrices
        // This method also checks the linear dependency of the new data to dismiss it
        BaseType::AppendCurrentIterationInformation(rResidualVector, rIterationGuess);

        // If the complete Jacobian is required for the IBQN equations calculate it
        // Note that this is only done after the first iteration since the observation matrices are empty at the first one
        if (BaseType::IsUsedInBlockNewtonEquations() && BaseType::GetConvergenceAcceleratorIteration() != 0) {
            CalculateInverseJacobianRandomizedSVD();
        }
    }

    void CalculateCorrectionWithJacobian(VectorType& rCorrection) override
    {
        const SizeType n_dofs = BaseType::GetProblemSize();
        const auto p_V = BaseType::pGetResidualObservationMatrix();
        const auto p_W = BaseType::pGetSolutionObservationMatrix();
        auto& r_res_vect = *(BaseType::pGetCurrentIterationResidualVector());

        if (p_V != nullptr && p_W != nullptr) {
            // Get the required arrays
            const auto& r_V = *p_V;
            const auto& r_W = *p_W;

            // Calculate the correction as the inverse Jacobian approximation times the current iteration residual
            MatrixType aux_M;
            CalculateAuxiliaryMatrixM(aux_M);
            VectorType M_res = prod(aux_M, r_res_vect);
            VectorType V_M_res = prod(r_V, M_res);
            VectorType W_M_res = prod(r_W, M_res);

            noalias(rCorrection) = ZeroVector(n_dofs);
            TDenseSpace::UnaliasedAdd(rCorrection, 1.0, V_M_res);
            TDenseSpace::UnaliasedAdd(rCorrection, 1.0, W_M_res);
            TDenseSpace::UnaliasedAdd(rCorrection, -1.0, r_res_vect);

            if (mpOldJacQU != nullptr && mpOldJacSigmaV != nullptr) {
                const auto& r_A = *mpOldJacQU;
                const auto& r_B = *mpOldJacSigmaV;
                VectorType B_res = prod(r_B, r_res_vect);
                VectorType A_B_res = prod(r_A, B_res);
                VectorType B_V_M_res = prod(r_B, V_M_res);
                VectorType A_B_V_M_res = prod(r_A, B_V_M_res);
                TDenseSpace::UnaliasedAdd(rCorrection, 1.0, A_B_res);
                TDenseSpace::UnaliasedAdd(rCorrection, -1.0, A_B_V_M_res);
            }
        } else {
            if (mpOldJacQU != nullptr && mpOldJacSigmaV != nullptr) {
                const auto& r_A = *mpOldJacQU;
                const auto& r_B = *mpOldJacSigmaV;
                const SizeType n_modes = TDenseSpace::Size2(r_A);
                VectorType B_res(n_modes);
                TDenseSpace::Mult(r_B, r_res_vect, B_res);
                TDenseSpace::Mult(r_A, B_res, rCorrection);
                TDenseSpace::UnaliasedAdd(rCorrection, -1.0, r_res_vect);
            } else {
                KRATOS_ERROR << "There is neither observation nor old Jacobian decomposition. Correction cannot be computed." << std::endl;
            }
        }
    }

    void CalculateInverseJacobianRandomizedSVD()
    {
        const auto p_W = BaseType::pGetSolutionObservationMatrix();
        const auto p_V = BaseType::pGetResidualObservationMatrix();

        // If there exists new observations, use these to calculate the inverse Jacobian SVD
        // If not, the old Jacobian decomposition will be used
        if (p_V != nullptr && p_W != nullptr) {
            // Compute ((trans(V)*V)^-1)*trans(V)
            MatrixType aux_M;
            CalculateAuxiliaryMatrixM(aux_M);

            // If not initialized yet, create the random values matrix for the truncated SVD
            if (!mRandomValuesAreInitialized || mpOmega == nullptr) {
                InitializeRandomValuesMatrix();
            }

            // Do the randomized SVD decomposition of the current Jacobian approximation
            // Note that the identity part of the Jacobian will not be considered in the randomized SVD as this is full rank
            MatrixType y;
            MultiplyRight(aux_M, *mpOmega, y);

            mpDenseQR->Compute(y);
            const SizeType n_dofs = BaseType::GetProblemSize();
            const SizeType n_modes = TDenseSpace::Size2(*mpOmega);
            MatrixType Q(n_dofs, n_modes);
            mpDenseQR->MatrixQ(Q);

            MatrixType phi;
            MultiplyTransposeLeft(aux_M, Q, phi);

            VectorType s_svd; // Singular values vector
            MatrixType u_svd; // Left orthogonal matrix
            MatrixType v_svd; // Right orthogonal matrix
            Parameters svd_settings(R"({
                "compute_thin_u" : true,
                "compute_thin_v" : true
            })");
            mpDenseSVD->Compute(phi, s_svd, u_svd, v_svd, svd_settings);

            // Save the decomposition matrices
            // Note that the added extra modes are discarded in here
            SizeType n_modes_final = n_modes - mNumberOfExtraModes;
            auto p_aux_Q_U = Kratos::make_shared<MatrixType>(ZeroMatrix(n_dofs, n_modes_final));
            auto p_aux_sigma_V = Kratos::make_shared<MatrixType>(n_modes_final, n_dofs);
            auto& r_aux_Q_U = *p_aux_Q_U;
            auto& r_aux_sigma_V = *p_aux_sigma_V;
            IndexPartition<SizeType>(n_dofs).for_each([&](SizeType I){
                for (SizeType j = 0; j < n_modes_final; ++j) {
                    for (SizeType k = 0; k < n_modes_final; ++k) {
                        r_aux_Q_U(I,j) += Q(I,k) * u_svd(k,j);
                    }
                    r_aux_sigma_V(j,I) = s_svd(j) * v_svd(I,j);
                }
            });

            std::swap(mpJacQU, p_aux_Q_U);
            std::swap(mpJacSigmaV, p_aux_sigma_V);

            // Clear the random values matrix
            if (mLimitModesToIterations) {
                mpOmega = nullptr;
            }

        } else {
            mpJacQU = nullptr;
            mpJacSigmaV = nullptr;
        }
    }

    ///@}
    ///@name Protected  Access
    ///@{

    MatrixPointerType pGetJacobianDecompositionMatrixQU() override
    {
        return mpJacQU;
    }

    MatrixPointerType pGetJacobianDecompositionMatrixSigmaV() override
    {
        return mpJacSigmaV;
    }

    MatrixPointerType pGetOldJacobianDecompositionMatrixQU() override
    {
        return mpOldJacQU;
    }

    MatrixPointerType pGetOldJacobianDecompositionMatrixSigmaV() override
    {
        return mpOldJacSigmaV;
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    SizeType mSeed = 0; // Seed to be used in the random matrix generation (to be updated at each step)
    SizeType mUserNumberOfModes; // User-defined number of modes to be kept in the Jacobian randomized SVD
    SizeType mNumberOfExtraModes; // Number of extra modes used in the randomization
    SizeType mMinRandSVDExtraModes; // Minimum number of extra modes used in the randomization
    SizeType mCurrentNumberOfModes = 0; // Current number of modes to be kept in the Jacobian randomized SVD //TODO: I THINK WE CAN REMOVE THIS MEMBER VAR
    bool mAutomaticJacobianModes = true; // Automatic selection of randomized SVD number of modes
    bool mLimitModesToIterations = true; // Limits the number of modes to the current iterations
    bool mRandomValuesAreInitialized = false; // Indicates if the random values for the truncated SVD have been already set

    DenseQRPointerType mpDenseQR; // Pointer to the dense QR utility
    DenseSVDPointerType mpDenseSVD; // Pointer to the dense SVD utility

    MatrixPointerType mpOmega = nullptr; // Matrix with random values for truncated SVD
    MatrixPointerType mpJacQU = nullptr; // Left DOFs x modes matrix from the current Jacobian decomposition
    MatrixPointerType mpJacSigmaV = nullptr; // Right modes x DOFs matrix from the current Jacobian decomposition
    MatrixPointerType mpOldJacQU = nullptr; // Left DOFs x modes matrix from the previous Jacobian decomposition
    MatrixPointerType mpOldJacSigmaV = nullptr; // Right modes x DOFs matrix from the previous Jacobian decomposition

    SizeType mOldMaxRank = 0; // Previous Jacobian rank (computed from the SVD decomposition)

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void InitializeRandomValuesMatrix()
    {
        // Initialize the number of modes
        // Note that if automatic modes are selected we set the previous step rank plus the current number of observations
        // This comes from the additivity property rank(A+B) \leq rank(A) + rank(B). Hence, the current Jacobian rank is at
        // most the summation of the old Jacobian rank plus the number of current of observations
        const SizeType n_obs = BaseType::GetNumberOfObservations();
        const SizeType full_rank_modes = mOldMaxRank + n_obs;
        SizeType n_modes;
        if (mAutomaticJacobianModes) {
            n_modes = full_rank_modes;
        } else {
            // If the number of modes is not automatic we need to check if it needs to be limited to the number of iterations (e.g. IBQN)
            // Note that we check it with the current number of observations rather than using the iteration counter
            // This is important in case some iteration information has been dropped because its linear dependency
            // Also note that we consider the previous Jacobian (via the full_rank_modes variable) as this contributes to the rank
            //TODO: HERE WE SHOULD SET A CHECK IF MAXIMUM MODES IS NOT ZERO DO THIS. THE AUTOMATIC IS SAME AS LIMIT MODES TO ITERATIONS
            n_modes = mLimitModesToIterations ? std::min(full_rank_modes, mUserNumberOfModes) : mUserNumberOfModes;
        }

        // Initialize auxiliary variables
        const SizeType aux_extra_modes = std::ceil(0.2 * n_modes);
        mNumberOfExtraModes = mMinRandSVDExtraModes > aux_extra_modes ? mMinRandSVDExtraModes : aux_extra_modes;
        mCurrentNumberOfModes = n_modes + mNumberOfExtraModes;

        // Set the random values matrix pointer
        const SizeType n_dofs = BaseType::GetProblemSize();
        auto p_aux_omega = Kratos::make_shared<MatrixType>(n_dofs, mCurrentNumberOfModes);
        std::swap(p_aux_omega, mpOmega);

        // Create the random values generator
        std::mt19937 generator(mSeed); //Standard mersenne_twister_engine seeded with the step number as seed
        std::uniform_real_distribution<> distribution(0.0, 1.0);

        // Fill the random values matrix
        auto& r_omega_matrix = *mpOmega;
        for (SizeType i = 0; i < n_dofs; ++i) {
            for (SizeType j = 0; j < mCurrentNumberOfModes; ++j) {
                // Use distribution to transform the random unsigned int generated by generator into a
                // double in [0.0, 1.0). Each call to distribution(generator) generates a new random double
                r_omega_matrix(i,j) = distribution(generator);
            }
        }

        // Set the flag to avoid performing this operation again (the random values are stored)
        mRandomValuesAreInitialized = !mLimitModesToIterations;
    }

    void MultiplyRight(
        const Matrix& rAuxM,
        const Matrix& rRightMatrix,
        Matrix& rSolution)
    {
        const SizeType n_dofs = BaseType::GetProblemSize();
        KRATOS_ERROR_IF(TDenseSpace::Size1(rRightMatrix) != n_dofs) << "Obtained right multiplication matrix size " << TDenseSpace::Size1(rRightMatrix) << " does not match the problem size " << n_dofs << " expected one." << std::endl;
        if (TDenseSpace::Size1(rSolution) != n_dofs || TDenseSpace::Size2(rSolution) != mCurrentNumberOfModes) {
            rSolution.resize(n_dofs, mCurrentNumberOfModes);
        }

        // Get observation matrices from base MVQN convergence accelerator
        const auto& r_W = *(BaseType::pGetSolutionObservationMatrix());
        const auto& r_V = *(BaseType::pGetResidualObservationMatrix());

        // Add to the solution matrix
        MatrixType M_omega = prod(rAuxM, rRightMatrix);
        noalias(rSolution) = prod(r_W, M_omega);

        if (mpOldJacQU == nullptr && mpOldJacSigmaV == nullptr) {
            if (!BaseType::IsUsedInBlockNewtonEquations()) {
                // Old Jacobian initialized to minus the identity matrix
                MatrixType V_M_omega = prod(r_V, M_omega);
                IndexPartition<SizeType>(n_dofs).for_each([&rSolution,&V_M_omega,this](SizeType I){
                    for (SizeType j = 0; j < mCurrentNumberOfModes; ++j) {
                        rSolution(I,j) += V_M_omega(I,j);
                    }
                });
            }
        } else {
            MatrixType V_M_omega = prod(r_V, M_omega);
            MatrixType B_omega = prod(*mpOldJacSigmaV, rRightMatrix);
            MatrixType A_B_omega = prod(*mpOldJacQU, B_omega);
            MatrixType B_V_M_omega = prod(*mpOldJacSigmaV, V_M_omega);
            MatrixType A_B_V_M_omega = prod(*mpOldJacQU, B_V_M_omega);
            if (!BaseType::IsUsedInBlockNewtonEquations()) {
                IndexPartition<SizeType>(n_dofs).for_each([&rSolution,&A_B_omega,&V_M_omega,&A_B_V_M_omega,this](SizeType I){
                    for (SizeType j = 0; j < mCurrentNumberOfModes; ++j) {
                        rSolution(I,j) += A_B_omega(I,j) + V_M_omega(I,j) - A_B_V_M_omega(I,j);
                    }
                });
            } else {
                IndexPartition<SizeType>(n_dofs).for_each([&rSolution,&A_B_omega,&A_B_V_M_omega,this](SizeType I){
                    for (SizeType j = 0; j < mCurrentNumberOfModes; ++j) {
                        rSolution(I,j) += A_B_omega(I,j) - A_B_V_M_omega(I,j);
                    }
                });
            }
        }
    }

    void MultiplyTransposeLeft(
        const Matrix& rAuxM,
        const Matrix& rLeftMatrix,
        Matrix& rSolution)
    {
        const SizeType n_dofs = BaseType::GetProblemSize();
        KRATOS_ERROR_IF(TDenseSpace::Size1(rLeftMatrix) != n_dofs) << "Obtained left multiplication matrix size " << TDenseSpace::Size1(rLeftMatrix) << " does not match the problem size " << n_dofs << " expected one." << std::endl;
        if (TDenseSpace::Size1(rSolution) != mCurrentNumberOfModes|| TDenseSpace::Size2(rSolution) != n_dofs) {
            rSolution.resize(mCurrentNumberOfModes, n_dofs);
        }

        // Get observation matrices from base MVQN convergence accelerator
        const auto& r_W = *(BaseType::pGetSolutionObservationMatrix());
        const auto& r_V = *(BaseType::pGetResidualObservationMatrix());

        // Add to the solution matrix
        const auto Qtrans = trans(rLeftMatrix);
        MatrixType Qtrans_W = prod(Qtrans, r_W);
        noalias(rSolution) = prod(Qtrans_W, rAuxM);

        if (mpOldJacQU == nullptr && mpOldJacSigmaV == nullptr) {
            if (!BaseType::IsUsedInBlockNewtonEquations()) {
                // Old Jacobian initialized to minus the identity matrix
                MatrixType Qtrans_V = prod(Qtrans, r_V);
                MatrixType Qtrans_V_M = prod(Qtrans_V, rAuxM);
                IndexPartition<SizeType>(n_dofs).for_each([&rSolution,&Qtrans_V_M,this](SizeType J){
                    for (SizeType i = 0; i < mCurrentNumberOfModes; ++i) {
                        rSolution(i,J) += Qtrans_V_M(i,J);
                    }
                });
            }
        } else {
            MatrixType Qtrans_V = prod(Qtrans, r_V);
            MatrixType Qtrans_V_M = prod(Qtrans_V, rAuxM);
            MatrixType Qtrans_A = prod(Qtrans, *mpOldJacQU);
            MatrixType Qtrans_A_B = prod(Qtrans_A, *mpOldJacSigmaV);
            MatrixType Qtrans_A_B_V = prod(Qtrans_A_B, r_V);
            MatrixType Qtrans_A_B_V_M = prod(Qtrans_A_B_V, rAuxM);
            if (!BaseType::IsUsedInBlockNewtonEquations()) {
                IndexPartition<SizeType>(n_dofs).for_each([&rSolution,&Qtrans_A_B,&Qtrans_V_M,&Qtrans_A_B_V_M,this](SizeType J){
                    for (SizeType i = 0; i < mCurrentNumberOfModes; ++i) {
                        rSolution(i,J) += Qtrans_A_B(i,J) + Qtrans_V_M(i,J) - Qtrans_A_B_V_M(i,J);
                    }
                });
            } else {
                IndexPartition<SizeType>(n_dofs).for_each([&rSolution,&Qtrans_A_B,&Qtrans_A_B_V_M,this](SizeType J){
                    for (SizeType i = 0; i < mCurrentNumberOfModes; ++i) {
                        rSolution(i,J) += Qtrans_A_B(i,J) - Qtrans_A_B_V_M(i,J);
                    }
                });
            }
        }
    }

    void CalculateAuxiliaryMatrixM(MatrixType& rAuxM)
    {
        // Do the QR decomposition of V
        auto& r_V = *(BaseType::pGetResidualObservationMatrix());
        mpDenseQR->Compute(r_V);

        // Get the QR decomposition matrices
        // Note that in here we are assuming that the pivoting QR is used
        const std::size_t n_dofs = TDenseSpace::Size1(r_V);
        const std::size_t n_data_cols = TDenseSpace::Size2(r_V);
        MatrixType Q(n_dofs, n_data_cols);
        MatrixType R(n_data_cols, n_data_cols);
        MatrixType P(n_data_cols, n_data_cols);
        mpDenseQR->MatrixQ(Q);
        mpDenseQR->MatrixR(R);
        mpDenseQR->MatrixP(P);

        // Calculate the SVD decomposition of R
        // Note that we firstly undo the pivoting by doing R*inv(P)
        // Also note that since P is orthogonal we actually do R*trans(P)
        MatrixType R_transP(n_data_cols, n_data_cols);
        noalias(R_transP) = prod(R, trans(P));

        // Compute the Moore-Penrose pseudo-inverse of R
        MatrixType R_transP_inv;
        CalculateMoorePenroseInverse(R_transP, R_transP_inv);

        // Calculate the solution of the problem inv(trans(V)*V)*trans(V)
        // The problem is written to avoid the transpose multiplication as V*M=I
        // Hence, from the previous QR decomposition we can do Q*R*M=I that is M=pseudoinv(R)*trans(Q)
        const std::size_t m = TDenseSpace::Size1(rAuxM);
        const std::size_t n = TDenseSpace::Size2(rAuxM);
        if (m != n_data_cols || n != n_dofs) {
            TDenseSpace::Resize(rAuxM, n_data_cols, n_dofs);
        }
        const MatrixType trans_Q = trans(Q);
        noalias(rAuxM) = prod(R_transP_inv, trans_Q);
    }

    //TODO: Whe should add this to the MathUtils
    void CalculateMoorePenroseInverse(
        const MatrixType& rInputMatrix,
        MatrixType& rMoorePenroseInverse)
    {
        IndexType aux_size_1 = TDenseSpace::Size1(rInputMatrix);
        IndexType aux_size_2 = TDenseSpace::Size2(rInputMatrix);
        KRATOS_ERROR_IF_NOT(aux_size_1 == aux_size_2) << "Input matrix is not squared." << std::endl;

        VectorType s_svd; // Singular values vector
        MatrixType u_svd; // Left orthogonal matrix
        MatrixType v_svd; // Right orthogonal matrix
        Parameters svd_settings(R"({
            "compute_thin_u" : true,
            "compute_thin_v" : true
        })");
        mpDenseSVD->Compute(const_cast<MatrixType&>(rInputMatrix), s_svd, u_svd, v_svd, svd_settings);
        const std::size_t n_sing_val = s_svd.size();
        //TODO: This allocation can be avoided
        MatrixType s_inv = ZeroMatrix(n_sing_val, n_sing_val);
        for (std::size_t i = 0; i < n_sing_val; ++i) {
            s_inv(i,i) = 1.0 / s_svd(i);
        }

        // Calculate and save the input matrix pseudo-inverse
        rMoorePenroseInverse = ZeroMatrix(aux_size_2, aux_size_1);

        // Note that we take advantage of the fact that the input matrix is always square
        for (std::size_t i = 0; i < aux_size_2; ++i) {
            for (std::size_t j = 0; j < aux_size_1; ++j) {
                double& r_value = rMoorePenroseInverse(i,j);
                for (std::size_t k = 0; k < n_sing_val; ++k) {
                    const double v_ik = v_svd(i,k);
                    for (std::size_t m = 0; m < n_sing_val; ++m) {
                        const double ut_mj = u_svd(j,m);
                        const double s_inv_km = s_inv(k,m);
                        r_value += v_ik * s_inv_km * ut_mj;
                    }
                }
            }
        }
    }

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
}; /* Class MVQNRandomizedSVDConvergenceAccelerator */

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_MVQN_RANDOMIZED_SVD_CONVERGENCE_ACCELERATOR defined */
