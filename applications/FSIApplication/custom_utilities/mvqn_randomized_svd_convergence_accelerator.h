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
#include "includes/ublas_interface.h"
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"
#include "utilities/dense_svd_decomposition.h"
#include "utilities/parallel_utilities.h"
#include "utilities/qr_utility.h"
#include "utilities/svd_utils.h"

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

    typedef typename DenseSingularValueDecomposition<TDenseSpace>::Pointer DenseSVDPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new MVQNRandomizedSVDConvergenceAccelerator object
     * MultiVector Quasi-Newton convergence accelerator with randomized SVD Jacobian json settings constructor
     * This constructor also makes possible the MVQN usage for the interface block Newton equations
     * @param pDenseSVD Pointer to the dense SVD instance
     * @param ConvAcceleratorParameters Json string encapsulating the settings
     */
    explicit MVQNRandomizedSVDConvergenceAccelerator(
        DenseSVDPointerType pDenseSVD,
        Parameters ConvAcceleratorParameters)
    : BaseType()
    , mpDenseSVD(pDenseSVD)
    {
        ConvAcceleratorParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        mNumberOfModes = ConvAcceleratorParameters["jacobian_modes"].GetInt();
        BaseType::SetInitialRelaxationOmega(ConvAcceleratorParameters["w_0"].GetDouble());
        BaseType::SetCutOffTolerance(ConvAcceleratorParameters["cut_off_tol"].GetDouble());
    }

    explicit MVQNRandomizedSVDConvergenceAccelerator(
        DenseSVDPointerType pDenseSVD,
        const unsigned int JacobianModes = 10,
        const double AbsCutOff = 1e-8)
    : BaseType(AbsCutOff)
    , mNumberOfModes(JacobianModes)
    , mpDenseSVD(pDenseSVD)
    {
    }

    /**
     * Copy Constructor.
     */

    // Required to build the IBN-MVQN
    MVQNRandomizedSVDConvergenceAccelerator(const MVQNRandomizedSVDConvergenceAccelerator& rOther);

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

    /**
     * @brief Save the current step Jacobian
     * This method saves the current step Jacobian as previous step Jacobian for the next time step iteration
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        //TODO: REMOVE THIS AFTER DEVELOPMENT (IT ONLY SETS THE CURRENT AS OLD JACOBIAN)
        // BaseType::FinalizeSolutionStep();
        //TODO: REMOVE THIS AFTER DEVELOPMENT (IT ONLY SETS THE CURRENT AS OLD JACOBIAN)

        // Compute ((trans(V)*V)^-1)*trans(V)
        MatrixType aux_M;
        CalculateAuxiliaryMatrixM(aux_M);

        KRATOS_WATCH(aux_M.size1())
        KRATOS_WATCH(aux_M.size2())

        // If not initialized yet, create the random values matrix for the truncated SVD
        if (!mRandomValuesAreInitialized) {
            InitializeRandomValuesMatrix();
        }

        // Do the truncated SVD decomposition of the current Jacobian approximation to save it for the next iteration
        MatrixType y;
        MultiplyRight(aux_M, *mpOmega, y);
        KRATOS_WATCH(y.size1())
        KRATOS_WATCH(y.size2())

        QR<double, row_major> qr_util;
        const SizeType n_dofs = BaseType::GetProblemSize();
        qr_util.compute(n_dofs, mNumberOfModes, &(y(0,0)));
        qr_util.compute_q();
        MatrixType Q(n_dofs, mNumberOfModes);
        for (SizeType i = 0; i < n_dofs; ++i) { //TODO: This can be parallel
            for (SizeType j = 0; j < mNumberOfModes; ++j) {
                Q(i,j) = qr_util.Q(i,j);
            }
        }

        KRATOS_WATCH(Q.size1())
        KRATOS_WATCH(Q.size2())
        auto q_qtrans_check = prod(Q, trans(Q));
        // auto q_qtrans_check = prod(trans(Q), Q);
        // KRATOS_WATCH(q_qtrans_check)
        double identity_error = 0.0;
        for (SizeType i = 0; i < q_qtrans_check.size1(); ++i) {
            for (SizeType j = 0; j < q_qtrans_check.size2(); ++j) {
                const double ref = i == j ? 1.0 : 0.0;
                identity_error += std::pow(q_qtrans_check(i,j)-ref, 2);
            }
        }
        KRATOS_WATCH(std::sqrt(identity_error))


        MatrixType phi;
        MultiplyTransposeLeft(aux_M, Q, phi);

        KRATOS_WATCH(phi.size1())
        KRATOS_WATCH(phi.size2())
        KRATOS_WATCH("Before SVD")

        VectorType s_svd; // Eigenvalues vector
        MatrixType u_svd; // Left orthogonal matrix
        MatrixType v_svd; // Right orthogonal matrix
        Parameters svd_settings(R"({
            "compute_thin_u" : true,
            "compute_thin_v" : true
        })");
        mpDenseSVD->Compute(phi, s_svd, u_svd, v_svd, svd_settings);

        KRATOS_WATCH("After SVD")

        KRATOS_WATCH(u_svd.size1())
        KRATOS_WATCH(u_svd.size2())
        KRATOS_WATCH(s_svd.size())
        KRATOS_WATCH(v_svd.size1())
        KRATOS_WATCH(v_svd.size2())

        auto p_aux_Q_U = Kratos::make_unique<MatrixType>(prod(Q, u_svd));
        auto p_aux_sigma_V = Kratos::make_unique<MatrixType>(v_svd.size2(), v_svd.size1());
        auto& r_aux_sigma_V = *p_aux_sigma_V;
        for (SizeType i = 0; i < r_aux_sigma_V.size1(); ++i) {
            const double aux_s = s_svd(i);
            for (SizeType j = 0; j < r_aux_sigma_V.size2(); ++j) {
                r_aux_sigma_V(i,j) = aux_s * v_svd(j,i);
            }
        }
        std::swap(mpOldJacQU, p_aux_Q_U);
        std::swap(mpOldJacSigmaV, p_aux_sigma_V);

        KRATOS_WATCH(mpOldJacQU->size1())
        KRATOS_WATCH(mpOldJacQU->size2())
        KRATOS_WATCH(mpOldJacSigmaV->size1())
        KRATOS_WATCH(mpOldJacSigmaV->size2())

        // MatrixType reconstructed_jacobian = prod(*mpOldJacQU, *mpOldJacSigmaV);
        // for (SizeType i = 0; i < reconstructed_jacobian.size1(); ++i) {
        //     reconstructed_jacobian(i,i) -= 1.0;
        // }

        // double error_norm = 0.0;
        // const auto& r_current_jac = *(BaseType::pGetInverseJacobianApproximation());
        // for (SizeType i = 0; i < reconstructed_jacobian.size1(); ++i) {
        //     for (SizeType j = 0; j < reconstructed_jacobian.size2(); ++j) {
        //         error_norm += std::pow(reconstructed_jacobian(i,j) - r_current_jac(i,j), 2);
        //     }
        // }
        // std::cout << "Jacobian error norm: " << std::sqrt(error_norm) << std::endl;

        // KRATOS_WATCH(reconstructed_jacobian(0,0))
        // KRATOS_WATCH(r_current_jac(0,0))
        // KRATOS_WATCH(reconstructed_jacobian(50,50))
        // KRATOS_WATCH(r_current_jac(50,50))
        // KRATOS_WATCH(reconstructed_jacobian(200,200))
        // KRATOS_WATCH(r_current_jac(200,200))
        // KRATOS_WATCH(reconstructed_jacobian(10,50))
        // KRATOS_WATCH(r_current_jac(10,50))
        // KRATOS_WATCH(reconstructed_jacobian(100,50))
        // KRATOS_WATCH(r_current_jac(100,50))

        // // Output values to a file to do Python testing
        // std::ofstream out_jacobian("out_jacobian.txt");
        // std::cout.rdbuf(out_jacobian.rdbuf());
        // for (SizeType i = 0; i < r_current_jac.size1(); ++i) {
        //     for (SizeType j = 0; j < r_current_jac.size2(); ++j) {
        //         std::cout << r_current_jac(i,j);
        //         if (j != r_current_jac.size2() - 1) {
        //             std::cout << " ";
        //         } else {
        //             std::cout << std::endl;
        //         }
        //     }
        // }

        // const auto& r_V = *(BaseType::pGetResidualObservationMatrix());
        // std::ofstream out_V("out_V.txt");
        // std::cout.rdbuf(out_V.rdbuf());
        // for (SizeType i = 0; i < r_V.size1(); ++i) {
        //     for (SizeType j = 0; j < r_V.size2(); ++j) {
        //         std::cout << r_V(i,j);
        //         if (j != r_V.size2() - 1) {
        //             std::cout << " ";
        //         } else {
        //             std::cout << std::endl;
        //         }
        //     }
        // }

        // const auto& r_W = *(BaseType::pGetSolutionObservationMatrix());
        // std::ofstream out_W("out_W.txt");
        // std::cout.rdbuf(out_W.rdbuf());
        // for (SizeType i = 0; i < r_W.size1(); ++i) {
        //     for (SizeType j = 0; j < r_W.size2(); ++j) {
        //         std::cout << r_W(i,j);
        //         if (j != r_W.size2() - 1) {
        //             std::cout << " ";
        //         } else {
        //             std::cout << std::endl;
        //         }
        //     }
        // }

        KRATOS_CATCH( "" );
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters mvqn_randomized_svd_default_parameters(R"({
            "solver_type"            : "MVQN_randomized_SVD",
            "jacobian_modes"         : 10,
            "w_0"                    : 0.825,
            "cut_off_tol"            : 1e-8,
            "interface_block_newton" : false
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
        //FIXME: THIS IS A TEMPORARY SOLUTION UNTIL WE IMPLEMENT THE MATRIX FREE RESOLUTION METHOD
        // if (BaseType::IsUsedInBlockNewtonEquations()) {
        //     CalculateInverseJacobianApproximation()
        // }
    }

    void CalculateCorrectionWithJacobian(VectorType& rCorrection) override
    {
        KRATOS_WATCH("Inside CalculateCorrectionWithJacobian")

        // Get the required arrays
        const auto& r_V = *(BaseType::pGetResidualObservationMatrix());
        const auto& r_W = *(BaseType::pGetSolutionObservationMatrix());
        const auto& r_res_vect = *(BaseType::pGetCurrentIterationResidualVector());

        // Calculate the correction as the inverse Jacobian approximation times the current iteration residual
        MatrixType aux_M;
        CalculateAuxiliaryMatrixM(aux_M);
        VectorType M_res = prod(aux_M, r_res_vect);
        VectorType V_M_res = prod(r_V, M_res);
        VectorType W_M_res = prod(r_W, M_res);

        const SizeType n_dofs = BaseType::GetProblemSize();
        noalias(rCorrection) = ZeroVector(n_dofs);
        TDenseSpace::UnaliasedAdd(rCorrection, 1.0, V_M_res);
        TDenseSpace::UnaliasedAdd(rCorrection, 1.0, W_M_res);
        TDenseSpace::UnaliasedAdd(rCorrection, -1.0, r_res_vect);

        // KRATOS_WATCH(*mpOldJacQU)
        // KRATOS_WATCH(*mpOldJacSigmaV)
        if (mpOldJacQU != nullptr && mpOldJacSigmaV != nullptr) {
            KRATOS_WATCH("Adding AB info")
            const auto& r_A = *mpOldJacQU;
            const auto& r_B = *mpOldJacSigmaV;
            VectorType B_res = prod(r_B, r_res_vect);
            VectorType A_B_res = prod(r_A, B_res);
            VectorType B_V_M_res = prod(r_B, V_M_res);
            VectorType A_B_V_M_res = prod(r_A, B_V_M_res);
            TDenseSpace::UnaliasedAdd(rCorrection, 1.0, A_B_res);
            TDenseSpace::UnaliasedAdd(rCorrection, -1.0, A_B_V_M_res);
        }
    }

    void UpdateCurrentJacobianMatrix() override
    {
        // Compute the current inverse Jacobian approximation
        auto&& rp_jac_k1 = BaseType::pGetInverseJacobianApproximation();
        MatrixPointerType p_aux_jac_k1 = MatrixPointerType(new MatrixType(prod(*mpOldJacQU,*mpOldJacSigmaV)));
        std::swap(rp_jac_k1, p_aux_jac_k1);
    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    // double mOmega_0;                                                // Relaxation factor for the initial fixed point iteration
    // double mAbsCutOff;                                              // Tolerance for the absolute cut-off criterion
    // bool mUsedInBlockNewtonEquations;                               // Indicates if the current MVQN is to be used in the interface block Newton equations
    // unsigned int mProblemSize = 0;                                  // Residual to minimize size
    // unsigned int mConvergenceAcceleratorIteration = 0;              // Convergence accelerator iteration counter
    // bool mJacobiansAreInitialized = false;                          // Indicates that the Jacobian matrices have been already initialized
    SizeType mNumberOfModes;
    bool mRandomValuesAreInitialized = false;                       // Indicates if the random values for the truncated SVD have been already set
    // bool mConvergenceAcceleratorFirstCorrectionPerformed = false;   // Indicates that the initial fixed point iteration has been already performed

    DenseSVDPointerType mpDenseSVD;

    Kratos::unique_ptr<MatrixType> mpOmega;     // Matrix with random values for truncated SVD
    Kratos::unique_ptr<MatrixType> mpOldJacQU = nullptr;
    Kratos::unique_ptr<MatrixType> mpOldJacSigmaV = nullptr;

    ///@}
    ///@name Private Operators
    ///@{

    // void InitializeJacobianMatrices()
    // {
    //     if (mUsedInBlockNewtonEquations) {
    //         // Initialize the previous step Jacobian to zero
    //         MatrixPointerType p_new_jac_n = Kratos::make_shared<MatrixType>(mProblemSize,mProblemSize);
    //         (*p_new_jac_n) = ZeroMatrix(mProblemSize,mProblemSize);
    //         std::swap(p_new_jac_n,mpJac_n);

    //         // Initialize the current Jacobian approximation to a zero matrix
    //         // Note that this is only required for the Interface Block Newton algorithm
    //         MatrixPointerType p_new_jac_k1 = Kratos::make_shared<MatrixType>(mProblemSize,mProblemSize);
    //         (*p_new_jac_k1) = ZeroMatrix(mProblemSize,mProblemSize);
    //         std::swap(p_new_jac_k1, mpJac_k1);
    //     } else {
    //         // Initialize the previous step Jacobian approximation to minus the diagonal matrix
    //         MatrixPointerType p_new_jac_n = Kratos::make_shared<MatrixType>(mProblemSize,mProblemSize);
    //         (*p_new_jac_n) = -1.0 * IdentityMatrix(mProblemSize,mProblemSize);
    //         std::swap(p_new_jac_n,mpJac_n);
    //     }
    // }

    // /**
    //  * @brief Drop the last column of both observation matrices
    //  * This function drops the last column of both observation matrices
    //  */
    // void DropLastDataColumn()
    // {
    //     // Set two auxiliar observation matrices
    //     const auto n_cols = mpObsMatrixV->size2() - 1;
    //     MatrixPointerType p_aux_V = Kratos::make_shared<MatrixType>(mProblemSize, n_cols);
    //     MatrixPointerType p_aux_W = Kratos::make_shared<MatrixType>(mProblemSize, n_cols);

    //     // Drop the last column
    //     IndexPartition<std::size_t>(mProblemSize).for_each([&](unsigned int IRow){
    //         for (std::size_t i_col = 0; i_col < n_cols; ++i_col){
    //             (*p_aux_V)(IRow, i_col) = (*mpObsMatrixV)(IRow, i_col);
    //             (*p_aux_W)(IRow, i_col) = (*mpObsMatrixW)(IRow, i_col);
    //         }
    //     });

    //     // Set the member observation matrices pointers
    //     std::swap(mpObsMatrixV,p_aux_V);
    //     std::swap(mpObsMatrixW,p_aux_W);
    // }

    void InitializeRandomValuesMatrix()
    {
        // Set the random values matrix pointer
        const SizeType n_dofs = BaseType::GetProblemSize();
        auto p_aux_omega = Kratos::make_unique<MatrixType>(n_dofs, mNumberOfModes);
        std::swap(p_aux_omega, mpOmega);

        // Create the random values generator
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> distribution(0.0, 1.0);

        // Fill the random values matrix
        auto& r_omega_matrix = *mpOmega;
        for (SizeType i = 0; i < n_dofs; ++i) {
            for (SizeType j = 0; j < mNumberOfModes; ++j) {
                // Use distribution to transform the random unsigned int generated by generator into a
                // double in [0.0, 1.0). Each call to distribution(generator) generates a new random double
                r_omega_matrix(i,j) = distribution(generator);
            }
        }

        // Set the flag to avoid performing this operation again (the random values are stored)
        mRandomValuesAreInitialized = true;
    }

    void MultiplyRight(
        const Matrix& rAuxM,
        const Matrix& rRightMatrix,
        Matrix& rSolution)
    {
        const SizeType n_dofs = BaseType::GetProblemSize();
        KRATOS_ERROR_IF(rRightMatrix.size1() != n_dofs) << "Obtained right multiplication matrix size " << rRightMatrix.size1() << " does not match the problem size " << n_dofs << " expected one." << std::endl;
        if (rSolution.size1() != n_dofs|| rSolution.size2() != mNumberOfModes) {
            rSolution.resize(n_dofs, mNumberOfModes);
        }

        // Get observation matrices from base MVQN convergence accelerator
        const auto& r_W = *(BaseType::pGetSolutionObservationMatrix());
        const auto& r_V = *(BaseType::pGetResidualObservationMatrix());

        // Add to the solution matrix
        MatrixType M_omega = prod(rAuxM, rRightMatrix);
        noalias(rSolution) = prod(r_W, M_omega);

        if (mpOldJacQU == nullptr && mpOldJacSigmaV == nullptr) {
            if (!BaseType::IsUsedInBlockNewtonEquations()) {
                KRATOS_WATCH("Old Jacobian initialized as minus identity!")
                // Old Jacobian initialized to minus the identity matrix
                MatrixType V_M_omega = prod(r_V, M_omega);
                for (SizeType i = 0; i < rSolution.size1(); ++i) { //TODO: DO THIS PARALLEL
                    for (SizeType j = 0; j < rSolution.size2(); ++j) {
                        rSolution(i,j) += V_M_omega(i,j);
                    }
                }
            } else {
                // Old Jacobian initalized to zero
                // FIXME: THIS NEEDS TO BE STUDIED FURTHER
                for (SizeType i = 0; i < rSolution.size1(); ++i) { //TODO: DO THIS PARALLEL
                    for (SizeType j = 0; j < rSolution.size2(); ++j) {
                        rSolution(i,j) -= rRightMatrix(i,j);
                    }
                }
            }
        } else {
            KRATOS_WATCH(mpOldJacQU->size1())
            KRATOS_WATCH(mpOldJacQU->size2())
            KRATOS_WATCH(mpOldJacSigmaV->size1())
            KRATOS_WATCH(mpOldJacSigmaV->size2())
            MatrixType V_M_omega = prod(r_V, M_omega);
            MatrixType B_omega = prod(*mpOldJacSigmaV, rRightMatrix);
            MatrixType A_B_omega = prod(*mpOldJacQU, B_omega);
            MatrixType B_V_M_omega = prod(*mpOldJacSigmaV, V_M_omega);
            MatrixType A_B_V_M_omega = prod(*mpOldJacQU, B_V_M_omega);
            for (SizeType i = 0; i < rSolution.size1(); ++i) { //TODO: DO THIS PARALLEL
                for (SizeType j = 0; j < rSolution.size2(); ++j) {
                    rSolution(i,j) += A_B_omega(i,j) + V_M_omega(i,j) - A_B_V_M_omega(i,j);
                }
            }
        }
    }

    void MultiplyTransposeLeft(
        const Matrix& rAuxM,
        const Matrix& rLeftMatrix,
        Matrix& rSolution)
    {
        const SizeType n_dofs = BaseType::GetProblemSize();
        KRATOS_ERROR_IF(rLeftMatrix.size1() != n_dofs) << "Obtained left multiplication matrix size " << rLeftMatrix.size1() << " does not match the problem size " << n_dofs << " expected one." << std::endl;
        if (rSolution.size1() != mNumberOfModes|| rSolution.size2() != n_dofs) {
            rSolution.resize(mNumberOfModes, n_dofs);
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
                for (SizeType i = 0; i < rSolution.size1(); ++i) { //TODO: DO THIS PARALLEL
                    for (SizeType j = 0; j < rSolution.size2(); ++j) {
                        rSolution(i,j) += Qtrans_V_M(i,j);
                    }
                }
            } else {
                // Old Jacobian initialized to zero
                // FIXME: THIS NEEDS TO BE STUDIED FURTHER
                for (SizeType i = 0; i < rSolution.size1(); ++i) { //TODO: DO THIS PARALLEL
                    for (SizeType j = 0; j < rSolution.size2(); ++j) {
                        rSolution(i,j) -= Qtrans(i,j);
                    }
                }
            }
        } else {
            KRATOS_WATCH(mpOldJacQU->size1())
            KRATOS_WATCH(mpOldJacQU->size2())
            KRATOS_WATCH(mpOldJacSigmaV->size1())
            KRATOS_WATCH(mpOldJacSigmaV->size2())
            MatrixType Qtrans_V = prod(Qtrans, r_V);
            MatrixType Qtrans_V_M = prod(Qtrans_V, rAuxM);
            MatrixType Qtrans_A = prod(Qtrans, *mpOldJacQU);
            MatrixType Qtrans_A_B = prod(Qtrans_A, *mpOldJacSigmaV);
            MatrixType Qtrans_A_B_V = prod(Qtrans_A_B, r_V);
            MatrixType Qtrans_A_B_V_M = prod(Qtrans_A_B_V, rAuxM);
            for (SizeType i = 0; i < rSolution.size1(); ++i) { //TODO: DO THIS PARALLEL
                for (SizeType j = 0; j < rSolution.size2(); ++j) {
                    rSolution(i,j) += Qtrans_A_B(i,j) + Qtrans_V_M(i,j) - Qtrans_A_B_V_M(i,j);
                }
            }
        }
    }

    void CalculateAuxiliaryMatrixM(MatrixType& rAuxM)
    {
        // Compute (trans(V)*V)^-1
        const auto& r_V = *(BaseType::pGetResidualObservationMatrix());
        const auto Vtrans = trans(r_V);
        const SizeType n_data_cols = r_V.size2();
        MatrixType transV_V(n_data_cols, n_data_cols);
        noalias(transV_V) = prod(Vtrans, r_V);

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

        // Compute ((trans(V)*V)^-1)*trans(V)
        const SizeType n_dofs = BaseType::GetProblemSize();
        if (rAuxM.size1() != n_data_cols || rAuxM.size2() != n_dofs) {
            rAuxM.resize(n_data_cols, n_dofs);
        }
        noalias(rAuxM) = prod(transV_V_inv, Vtrans);
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
