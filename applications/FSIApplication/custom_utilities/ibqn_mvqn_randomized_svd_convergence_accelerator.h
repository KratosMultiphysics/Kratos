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

#if !defined(KRATOS_IBQN_MVQN_RANDOMIZED_SVD_CONVERGENCE_ACCELERATOR)
#define  KRATOS_IBQN_MVQN_RANDOMIZED_SVD_CONVERGENCE_ACCELERATOR

// System includes

// External includes

// Project includes
#include "utilities/dense_svd_decomposition.h"

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

/// Forward declaration of MVQN with ranzomized SVD Jacobian convergence accelerator
/// This is required since the include ward avoids the inclusion of the MVQNRandomizedSVDConvergenceAccelerator
template<class TSparseSpace, class TDenseSpace>
class MVQNRandomizedSVDConvergenceAccelerator;

/** @brief Interface Block Newton MVQN with randomized Jacobian convergence accelerator
 * Interface Block Newton equations convergence accelerator with MVQN randomized SVD Jacobian for the subdomains
 * @tparam TSparseSpace Linear algebra sparse space
 * @tparam TDenseSpace Linear algebra dense space
 */
template<class TSparseSpace, class TDenseSpace>
class IBQNMVQNRandomizedSVDConvergenceAccelerator: public IBQNMVQNConvergenceAccelerator<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( IBQNMVQNRandomizedSVDConvergenceAccelerator );

    typedef IBQNMVQNConvergenceAccelerator<TSparseSpace, TDenseSpace>      BaseType;
    typedef typename BaseType::Pointer                              BaseTypePointer;

    typedef typename BaseType::DenseVectorType                           VectorType;
    typedef typename BaseType::DenseVectorPointerType             VectorPointerType;

    typedef typename BaseType::DenseMatrixType                           MatrixType;
    typedef typename BaseType::DenseMatrixPointerType             MatrixPointerType;

    typedef MVQNRandomizedSVDConvergenceAccelerator<TSparseSpace, TDenseSpace>        MVQNType;
    typedef typename BaseType::MVQNPointerType                                 MVQNPointerType;

    typedef typename DenseSingularValueDecomposition<TDenseSpace>::Pointer DenseSVDPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new IBQNMVQNRandomizedSVDConvergenceAccelerator object
     * IBQN-MVQN with randomized SVD Jacobian convergence accelerator Json settings constructor
     * @param pDenseSVD Pointer to the dense SVD utility instance
     * @param rConvAcceleratorParameters Json string encapsulating the settings
     */
    explicit IBQNMVQNRandomizedSVDConvergenceAccelerator(
        DenseSVDPointerType pDenseSVD,
        Parameters ConvAcceleratorParameters)
    : BaseType()
    , mpSmallMatrixDenseSVD(pDenseSVD)
    {
        ConvAcceleratorParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        BaseType::SetInitialRelaxationOmega(ConvAcceleratorParameters["w_0"].GetDouble());

        // Set the subdomains MVQN randomized SVD convergence accelerator pointers
        // Note that we call the simplified constructor with default zero relaxation omega and IBQN switch
        const double cut_off_tol = ConvAcceleratorParameters["cut_off_tol"].GetDouble();
        const unsigned int jacobian_modes = ConvAcceleratorParameters["jacobian_modes"].GetInt();
        const bool limit_modes_to_iterations = ConvAcceleratorParameters["limit_modes_to_iterations"].GetBool();
        MVQNPointerType p_convergence_accelerator_left = Kratos::make_unique<MVQNType>(pDenseSVD, jacobian_modes, cut_off_tol, limit_modes_to_iterations);
        MVQNPointerType p_convergence_accelerator_right = Kratos::make_unique<MVQNType>(pDenseSVD, jacobian_modes, cut_off_tol, limit_modes_to_iterations);
        BaseType::SetLeftConvergenceAcceleratorPointer(p_convergence_accelerator_left);
        BaseType::SetRightConvergenceAcceleratorPointer(p_convergence_accelerator_right);
    }

    /**
     * Copy Constructor.
     */
    IBQNMVQNRandomizedSVDConvergenceAccelerator(const IBQNMVQNRandomizedSVDConvergenceAccelerator& rOther) = delete;

    /**
     * Destructor.
     */
    virtual ~IBQNMVQNRandomizedSVDConvergenceAccelerator(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Parameters GetDefaultParameters() const override
    {
        Parameters ibqn_mvqn_default_parameters(R"({
            "solver_type"               : "IBQN_MVQN_randomized_SVD",
            "jacobian_modes"            : 10,
            "w_0"                       : 0.825,
            "cut_off_tol"               : 1e-8,
            "limit_modes_to_iterations" : true
        })");

        return ibqn_mvqn_default_parameters;
    }

    void SolveInterfaceBlockQuasiNewtonLeftUpdate(
        const VectorType& rDisplacementInputVector,
        const VectorType& rForceOutputVector,
        const VectorType& rIterationGuess,
        VectorType& rLeftCorrection) override
    {
        // Get the subdomains convergence accelerators
        auto p_conv_acc_left = BaseType::pGetConvergenceAcceleratorLeft();
        auto p_conv_acc_right = BaseType::pGetConvergenceAcceleratorRight();

        // Get both left and inverse Jacobian
        auto p_inv_jac_QU_left = p_conv_acc_left->pGetJacobianDecompositionMatixQU();
        auto p_inv_jac_SigmaV_left = p_conv_acc_left->pGetJacobianDecompositionMatixSigmaV();
        auto p_inv_jac_QU_right = p_conv_acc_right->pGetJacobianDecompositionMatixQU();
        auto p_inv_jac_SigmaV_right = p_conv_acc_right->pGetJacobianDecompositionMatixSigmaV();
        auto p_uncorrected_displacement_vector = BaseType::pGetUncorrectedDisplacementVector();

        // Set the residual of the update problem to be solved
        VectorType aux_RHS = rForceOutputVector - rIterationGuess;
        if (p_inv_jac_QU_left != nullptr && p_inv_jac_SigmaV_left != nullptr) {
            const std::size_t n_dofs = TDenseSpace::Size(*p_uncorrected_displacement_vector);
            const std::size_t n_modes = TDenseSpace::Size1(*p_inv_jac_SigmaV_left);
            VectorType aux_vect(n_modes);
            VectorType aux_left_onto_right(n_dofs);
            VectorType displacement_iteration_update = *p_uncorrected_displacement_vector - rDisplacementInputVector;
            TSparseSpace::Mult(*p_inv_jac_SigmaV_left, displacement_iteration_update, aux_vect);
            TSparseSpace::Mult(*p_inv_jac_QU_left, aux_vect, aux_left_onto_right);
            TSparseSpace::UnaliasedAdd(aux_RHS, 1.0, aux_left_onto_right);

            if (p_inv_jac_QU_right != nullptr && p_inv_jac_SigmaV_right != nullptr) {
                // Calculate auxiliary intermediate small matrices
                MatrixType aux_C = ZeroMatrix(n_modes, n_modes);
                IndexPartition<std::size_t>(n_modes).for_each(
                    VectorType(n_modes),
                    [&aux_C, &p_inv_jac_SigmaV_left, &p_inv_jac_QU_right, &n_modes](std::size_t Col, VectorType& rAuxColumnTLS)
                {
                    TDenseSpace::GetColumn(Col, *p_inv_jac_QU_right, rAuxColumnTLS);
                    for (std::size_t row = 0; row < n_modes; ++row) {
                        aux_C(row,Col) = TDenseSpace::RowDot(row, *p_inv_jac_SigmaV_left, rAuxColumnTLS);
                    }
                });

                // Calculate the correction with the Woodbury matrix identity
                ApplyWoodburyMatrixIdentity(*p_inv_jac_QU_left, *p_inv_jac_SigmaV_right, aux_C, aux_RHS, rLeftCorrection);

            } else {
                rLeftCorrection = aux_RHS;
            }
        } else {
            KRATOS_WATCH("NO POINTERS LEFT")
            rLeftCorrection = aux_RHS;
        }
    }

    void SolveInterfaceBlockQuasiNewtonRightUpdate(
        const VectorType& rForceInputVector,
        const VectorType& rDisplacementOutputVector,
        const VectorType& rIterationGuess,
        VectorType& rRightCorrection) override
    {
        // Get the subdomains convergence accelerators
        auto p_conv_acc_left = BaseType::pGetConvergenceAcceleratorLeft();
        auto p_conv_acc_right = BaseType::pGetConvergenceAcceleratorRight();

        // Get both left and inverse Jacobian
        auto p_inv_jac_QU_left = p_conv_acc_left->pGetJacobianDecompositionMatixQU();
        auto p_inv_jac_SigmaV_left = p_conv_acc_left->pGetJacobianDecompositionMatixSigmaV();
        auto p_inv_jac_QU_right = p_conv_acc_right->pGetJacobianDecompositionMatixQU();
        auto p_inv_jac_SigmaV_right = p_conv_acc_right->pGetJacobianDecompositionMatixSigmaV();
        auto p_uncorrected_force_vector = BaseType::pGetUncorrectedForceVector();

        // Set the residual of the update problem to be solved
        VectorType aux_RHS = rDisplacementOutputVector - rIterationGuess;
        if (p_inv_jac_QU_right != nullptr && p_inv_jac_SigmaV_right != nullptr) {
            const std::size_t n_dofs = TDenseSpace::Size(*p_uncorrected_force_vector);
            const std::size_t n_modes_right = TDenseSpace::Size1(*p_inv_jac_SigmaV_right);
            VectorType aux_vect(n_modes_right);
            VectorType aux_right_onto_left(n_dofs);
            VectorType force_iteration_update = *p_uncorrected_force_vector - rForceInputVector;
            TSparseSpace::Mult(*p_inv_jac_SigmaV_right, force_iteration_update, aux_vect);
            TSparseSpace::Mult(*p_inv_jac_QU_right, aux_vect, aux_right_onto_left);
            TSparseSpace::UnaliasedAdd(aux_RHS, 1.0, aux_right_onto_left);

            KRATOS_WATCH(aux_RHS(0))
            KRATOS_WATCH(aux_RHS(10))
            KRATOS_WATCH(aux_RHS(13))
            KRATOS_WATCH(aux_RHS(23))

            // Calculate auxiliary intermediate small matrices
            if (p_inv_jac_QU_left != nullptr && p_inv_jac_SigmaV_left != nullptr) {

                auto jac_left = prod(*p_inv_jac_QU_left, *p_inv_jac_SigmaV_left);
                auto jac_right = prod(*p_inv_jac_QU_right, *p_inv_jac_SigmaV_right);
                KRATOS_WATCH(jac_left(100,100))
                KRATOS_WATCH(jac_right(100,100))
                KRATOS_WATCH(jac_left(100,200))
                KRATOS_WATCH(jac_right(100,200))

                const std::size_t n_modes_left = TDenseSpace::Size1(*p_inv_jac_SigmaV_left);
                MatrixType aux_C = ZeroMatrix(n_modes_right, n_modes_left);
                IndexPartition<std::size_t>(n_modes_left).for_each(
                    VectorType(n_dofs),
                    [&aux_C, &p_inv_jac_SigmaV_right, &p_inv_jac_QU_left, &n_modes_right](std::size_t Col, VectorType& rAuxColumnTLS)
                {
                    TDenseSpace::GetColumn(Col, *p_inv_jac_QU_left, rAuxColumnTLS);
                    for (std::size_t row = 0; row < n_modes_right; ++row) {
                        aux_C(row,Col) = TDenseSpace::RowDot(row, *p_inv_jac_SigmaV_right, rAuxColumnTLS);
                    }
                });

                // // Calculate the correction with the Woodbury matrix identity
                // ApplyWoodburyMatrixIdentity(*p_inv_jac_QU_right, *p_inv_jac_SigmaV_left, aux_C, aux_RHS, rRightCorrection);

                //TODO:REMOVE THIS AFTER TESTING
                MatrixType aux_LHS = IdentityMatrix(n_dofs, n_dofs);
                for (std::size_t i = 0; i < n_dofs; ++i) {
                    for (std::size_t j = 0; j < n_dofs; ++j) {
                        for (std::size_t k = 0; k < p_inv_jac_QU_right->size2(); ++k) {
                            for (std::size_t m = 0; m < p_inv_jac_SigmaV_left->size1(); ++m) {
                                aux_LHS(i,j) -= (*p_inv_jac_QU_right)(i,k) * aux_C(k,m) * (*p_inv_jac_SigmaV_left)(m,j);
                            }
                        }
                    }
                }
                QR<double, row_major> qr_util;
                qr_util.compute(n_dofs, n_dofs, &(aux_LHS)(0,0));
                qr_util.solve(&(aux_RHS)(0), &(rRightCorrection)(0));



            } else {
                rRightCorrection = aux_RHS;
            }
        } else {
            KRATOS_WATCH("NO POINTERS RIGHT")
            rRightCorrection = aux_RHS;
        }
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


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    DenseSVDPointerType mpSmallMatrixDenseSVD = nullptr;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Apply the Woodbury matrix identity to solve the update
     * This method applies the Woodbury matrix identity to solve the problem (I - ACB) * x = y
     */
    void ApplyWoodburyMatrixIdentity(
        MatrixType& rA,
        MatrixType& rB,
        const MatrixType& rC,
        VectorType& rY,
        VectorType& rX)
    {
        TDenseSpace::Assign(rX, 1.0, rY);

        MatrixType aux_C_inv;
        CalculateMoorePenroseInverse(rC, aux_C_inv);

        KRATOS_WATCH(rC)
        KRATOS_WATCH(aux_C_inv)

        KRATOS_WATCH(prod(MatrixType(prod(rC, aux_C_inv)), rC))

        const std::size_t n_dofs = TDenseSpace::Size(rY);
        const std::size_t n_col_A = TDenseSpace::Size2(rA);
        const std::size_t n_row_B = TDenseSpace::Size1(rB);

        MatrixType aux_D = ZeroMatrix(n_row_B, n_col_A);
        IndexPartition<std::size_t>(n_col_A).for_each(
            VectorType(n_dofs),
            [&aux_D, &rA, &rB, &aux_C_inv, &n_row_B](std::size_t Col, VectorType& rAuxColumnTLS)
        {
            TDenseSpace::GetColumn(Col, rA, rAuxColumnTLS);
            for (std::size_t row = 0; row < n_row_B; ++row) {
                aux_D(row,Col) = aux_C_inv(row,Col) - TDenseSpace::RowDot(row, rB, rAuxColumnTLS);
            }
        });

        MatrixType aux_D_inv;
        CalculateMoorePenroseInverse(aux_D, aux_D_inv);

        VectorType aux_vect(n_row_B);
        VectorType aux_vect_2(n_col_A);
        VectorType aux_vect_3(n_dofs);
        TDenseSpace::Mult(rB, rY, aux_vect);
        TDenseSpace::Mult(aux_D_inv, aux_vect, aux_vect_2);
        TDenseSpace::Mult(rA, aux_vect_2, aux_vect_3);

        TDenseSpace::UnaliasedAdd(rX, 1.0, aux_vect_3);
    }

    void CalculateMoorePenroseInverse(
        const MatrixType& rInputMatrix,
        MatrixType& rMoorePenroseInverse)
    {
        IndexType aux_size_1 = TDenseSpace::Size1(rInputMatrix);
        IndexType aux_size_2 = TDenseSpace::Size2(rInputMatrix);

        // MatrixType u_svd; // Orthogonal matrix (m x m)
        // MatrixType w_svd; // Rectangular diagonal matrix (m x n)
        // MatrixType v_svd; // Orthogonal matrix (n x n)
        // std::string svd_type = "Jacobi"; // SVD decomposition type
        // const double svd_rel_tol = 1.0e-6; // Relative tolerance of the SVD decomposition (it will be multiplied by the input matrix norm)
        // SVDUtils<double>::SingularValueDecomposition(rInputMatrix, u_svd, w_svd, v_svd, svd_type, svd_rel_tol);

        VectorType s_svd; // Eigenvalues vector
        MatrixType u_svd; // Left orthogonal matrix
        MatrixType v_svd; // Right orthogonal matrix
        Parameters svd_settings(R"({
            "compute_thin_u" : true,
            "compute_thin_v" : true
        })");
        // MatrixType input_test;
        // mpSmallMatrixDenseSVD->Compute(input_test, s_svd, u_svd, v_svd, svd_settings);
        mpSmallMatrixDenseSVD->Compute(const_cast<MatrixType&>(rInputMatrix), s_svd, u_svd, v_svd, svd_settings);
        const std::size_t n_eigs = s_svd.size();
        MatrixType s_inv = ZeroMatrix(n_eigs, n_eigs);
        for (std::size_t i = 0; i < n_eigs; ++i) {
            s_inv(i,i) = 1.0 / s_svd(i);
        }

        // Calculate and save the input matrix pseudo-inverse
        rMoorePenroseInverse = ZeroMatrix(aux_size_2, aux_size_1);

        // // Note that we take advantage of the fact that the input matrix is always square
        // for (std::size_t i = 0; i < aux_size; ++i) {
        //     for (std::size_t j = 0; j < aux_size; ++j) {
        //         const double aux = v_svd(j,i) / w_svd(j,j);
        //         for (std::size_t k = 0; k < aux_size; ++k) {
        //             rMoorePenroseInverse(i,k) += aux * u_svd(k,j);
        //         }
        //     }
        // }
        for (std::size_t i = 0; i < aux_size_2; ++i) {
            for (std::size_t j = 0; j < aux_size_1; ++j) {
                double& r_value = rMoorePenroseInverse(i,j);
                for (std::size_t k = 0; k < n_eigs; ++k) {
                    const double v_ik = v_svd(i,k);
                    for (std::size_t m = 0; m < n_eigs; ++m) {
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
}; /* Class IBQNMVQNRandomizedSVDConvergenceAccelerator */

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
}  /* namespace Kratos.*/

#endif /* KRATOS_IBQN_MVQN_RANDOMIZED_SVD_CONVERGENCE_ACCELERATOR defined */
