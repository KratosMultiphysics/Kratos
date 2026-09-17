//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolas Sibuet
//

#if !defined(ROM_RBF_UTILITY_H_INCLUDED)
#define ROM_RBF_UTILITY_H_INCLUDED

// System includes
#include <unordered_map>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <limits>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/builtin_timer.h"
#include "custom_utilities/ublas_wrapper.h"
#include "includes/ublas_interface.h"

// Application includes
#include "rom_application_variables.h"


namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>
class RomRBFUtility
{
public:

    using SizeType = std::size_t;

    using EigenDynamicMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using EigenDynamicVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

    /**
     * @brief Run the decoder on the input latent vector and also get its derivative over the input.
     * @param rRomUnknowns The latent vector to decode.
     * @param rx The resulting snapshot vector.
     * @param rPhiGlobal The PhiEffective matrix to be updated
     * @param rSVDPhiMatrices The list of projection matrices for the snapshot
     * @param r_W_mat The weights matrix for the RBF
     * @param r_centers_mat The collection of center vectors in matrix form
     * @param r_kernel_type Index indicating the kernel type
     * @param r_kernel_eps Epsilon value for the RBF's kernel
     * @param rReferenceSnapshot Reference snapshot to sum within the decoder
     */
    static void GetXAndDecoderGradient(
        Vector rRomUnknowns,
        Vector& rx,
        Matrix& rPhiGlobal,
        vector<Matrix>& rSVDPhiMatrices,
        Matrix& r_W_mat,
        Matrix& r_centers_mat,
        IndexType& r_kernel_type,
        double& r_kernel_eps,
        Vector& rReferenceSnapshot
    )
    {
        Eigen::Map<EigenDynamicMatrix> eigen_phi_inf(rSVDPhiMatrices[0].data().begin(), rSVDPhiMatrices[0].size1(), rSVDPhiMatrices[0].size2());
        Eigen::Map<EigenDynamicMatrix> eigen_phisig_sup(rSVDPhiMatrices[1].data().begin(), rSVDPhiMatrices[1].size1(), rSVDPhiMatrices[1].size2());  // Shape is num_dofs x n_sup
        Eigen::Map<EigenDynamicMatrix> eigen_sig_inv_inf(rSVDPhiMatrices[2].data().begin(), rSVDPhiMatrices[2].size1(), rSVDPhiMatrices[2].size2()); // Shape is n_inf x n_inf

        Eigen::Map<EigenDynamicVector> eigen_rom_unknowns(rRomUnknowns.data().begin(), rRomUnknowns.size());
        Eigen::Map<EigenDynamicMatrix> eigen_phi_global(rPhiGlobal.data().begin(), rPhiGlobal.size1(), rPhiGlobal.size2());
        Eigen::Map<EigenDynamicVector> eigen_rx(rx.data().begin(), rx.size());
        
        EigenDynamicVector q_inf_pred = eigen_sig_inv_inf*eigen_rom_unknowns;
        Eigen::Map<EigenDynamicMatrix> distance_mat(r_centers_mat.data().begin(), r_centers_mat.size1(), r_centers_mat.size2());

        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "(" << q_inf_pred.rows() << "," << q_inf_pred.cols() << ")" << std::endl;
        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "(" << distance_mat.rows() << "," << distance_mat.cols() << ")" << std::endl;

        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << q_inf_pred << std::endl;
        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << distance_mat(Eigen::seq(0,3),Eigen::placeholders::all) << std::endl;
        
        distance_mat.rowwise() -= q_inf_pred.transpose();
        distance_mat *= -1;  // Shape is n_centers x n_inf
        
        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "(" << distance_mat.rows() << "," << distance_mat.cols() << ")" << std::endl;
        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << distance_mat(Eigen::seq(0,3),Eigen::placeholders::all) << std::endl;

        EigenDynamicVector norms_vec = distance_mat.rowwise().norm();

        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "(" << norms_vec.rows() << "," << norms_vec.cols() << ")" << std::endl;
        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << norms_vec(Eigen::seq(0,3),Eigen::placeholders::all) << std::endl;

        EigenDynamicVector kernel_results_vec;
        EigenDynamicVector aux_derivatives_vec;
        
        if (r_kernel_type == 0) {
            kernel_results_vec = (-(r_kernel_eps * norms_vec.array()).square()).exp();
            aux_derivatives_vec = -2*std::pow(r_kernel_eps,2) * kernel_results_vec; // The derivative of the kernel would be missing multiplication by the norm,
                            // but then we would divide by the norm in the next step anyways as part of the derivative of the norm.
        } else if (r_kernel_type == 1) {
            EigenDynamicVector aux_term = 1.0 + (r_kernel_eps * norms_vec.array()).square();

            // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "Epsilon" << r_kernel_eps << std::endl;
            // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "(" << aux_term.rows() << "," << aux_term.cols() << ")" << std::endl;
            // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << aux_term(Eigen::seq(0,3),Eigen::placeholders::all) << std::endl;

            kernel_results_vec = 1.0 / (aux_term.array()).sqrt();
            aux_derivatives_vec = -aux_term.array().pow(-3.0/2.0)*std::pow(r_kernel_eps,2); // Also skips the multiplication by the norm

            // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "(" << kernel_results_vec.rows() << "," << kernel_results_vec.cols() << ")" << std::endl;
            // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << kernel_results_vec(Eigen::seq(0,3),Eigen::placeholders::all) << std::endl;
            // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "(" << aux_derivatives_vec.rows() << "," << aux_derivatives_vec.cols() << ")" << std::endl;
            // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << aux_derivatives_vec(Eigen::seq(0,3),Eigen::placeholders::all) << std::endl;

            // Tested until here
        }
        
        Eigen::Map<EigenDynamicMatrix> eigen_W_mat(r_W_mat.data().begin(), r_W_mat.size1(), r_W_mat.size2());  // Shape is num_centers x n_sup

        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "eigen_W_mat: " << "(" << eigen_W_mat.rows() << "," << eigen_W_mat.cols() << ")" << std::endl;

        EigenDynamicVector q_sup_aprox = eigen_W_mat.transpose() * kernel_results_vec;

        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "q_sup_aprox: " << "(" << q_sup_aprox.rows() << "," << q_sup_aprox.cols() << ")" << std::endl;
        KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "q_sup_aprox: " << q_sup_aprox(Eigen::seq(0,3),Eigen::placeholders::all) << std::endl;

        EigenDynamicMatrix rbf_gradient_aux =  aux_derivatives_vec.asDiagonal() * distance_mat;

        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "rbf_gradient_aux: " << "(" << rbf_gradient_aux.rows() << "," << rbf_gradient_aux.cols() << ")" << std::endl;

        EigenDynamicMatrix rbf_gradient = eigen_W_mat.transpose() * rbf_gradient_aux;  // Shape is n_sup x n_inf

        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "rbf_gradient: " << "(" << rbf_gradient.rows() << "," << rbf_gradient.cols() << ")" << std::endl;
        
        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "eigen_phisig_sup: " << "(" << eigen_phisig_sup.rows() << "," << eigen_phisig_sup.cols() << ")" << std::endl;
        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "eigen_sig_inv_inf: " << "(" << eigen_sig_inv_inf.rows() << "," << eigen_sig_inv_inf.cols() << ")" << std::endl;
        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "eigen_phi_inf: " << "(" << eigen_phi_inf.rows() << "," << eigen_phi_inf.cols() << ")" << std::endl;

        eigen_phi_global=eigen_phi_inf+eigen_phisig_sup*rbf_gradient*eigen_sig_inv_inf;  // Shape is num_dofs x n_inf

        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "eigen_phi_global: " << "(" << eigen_phi_global.rows() << "," << eigen_phi_global.cols() << ")" << std::endl;

        Eigen::Map<EigenDynamicVector> eigen_ref_snapshot(rReferenceSnapshot.data().begin(), rReferenceSnapshot.size());
        eigen_rx = eigen_ref_snapshot + eigen_phi_inf*eigen_rom_unknowns + eigen_phisig_sup*q_sup_aprox;

        // KRATOS_INFO("TESTING WITHIN RBF UTILITY:") << "eigen_rx: " << "(" << eigen_rx.rows() << "," << eigen_rx.cols() << ")" << std::endl;
    }

    /**
     * @brief Run the decoder on the input latent vector.
     * @param rRomUnknowns The latent vector to decode.
     * @param rx The resulting snapshot vector.
     * @param r_W_mat The weights matrix for the RBF
     * @param r_centers_mat The collection of center vectors in matrix form
     * @param r_kernel_type Index indicating the kernel type
     * @param r_kernel_eps Epsilon value for the RBF's kernel
     * @param rNNLayers The neural network's dense layers weights (no bias)
     * @param rReferenceSnapshot Reference snapshot to sum within the decoder
     */
    static void GetXFromDecoder(
        Vector rRomUnknowns,
        Vector& rx,
        vector<Matrix>& rSVDPhiMatrices,
        Matrix& r_W_mat,
        Matrix& r_centers_mat,
        IndexType& r_kernel_type,
        double& r_kernel_eps,
        Vector& rReferenceSnapshot
        )
    {
        Eigen::Map<EigenDynamicMatrix> eigen_phi_inf(rSVDPhiMatrices[0].data().begin(), rSVDPhiMatrices[0].size1(), rSVDPhiMatrices[0].size2());
        Eigen::Map<EigenDynamicMatrix> eigen_phisig_sup(rSVDPhiMatrices[1].data().begin(), rSVDPhiMatrices[1].size1(), rSVDPhiMatrices[1].size2());
        Eigen::Map<EigenDynamicMatrix> eigen_sig_inv_inf(rSVDPhiMatrices[2].data().begin(), rSVDPhiMatrices[2].size1(), rSVDPhiMatrices[2].size2());

        Eigen::Map<EigenDynamicVector> eigen_rom_unknowns(rRomUnknowns.data().begin(), rRomUnknowns.size());
        Eigen::Map<EigenDynamicVector> eigen_rx(rx.data().begin(), rx.size());
        
        EigenDynamicVector q_inf_pred = eigen_sig_inv_inf*eigen_rom_unknowns;
        Eigen::Map<EigenDynamicMatrix> distance_mat(r_centers_mat.data().begin(), r_centers_mat.size1(), r_centers_mat.size2());
        
        distance_mat.rowwise() -= q_inf_pred.transpose();
        distance_mat *= -1;

        EigenDynamicVector norms_vec = distance_mat.rowwise().norm();
        EigenDynamicVector kernel_results_vec;
        
        if (r_kernel_type == 0) {
            kernel_results_vec = (-(r_kernel_eps * norms_vec.array()).square()).exp();
        } else if (r_kernel_type == 1) {
            EigenDynamicVector aux_term = 1.0 + (r_kernel_eps * norms_vec.array()).square();
            kernel_results_vec = 1.0 / (aux_term.array()).sqrt();
        }
        
        Eigen::Map<EigenDynamicMatrix> eigen_W_mat(r_W_mat.data().begin(), r_W_mat.size1(), r_W_mat.size2());
        EigenDynamicVector q_sup_aprox = eigen_W_mat.transpose() * kernel_results_vec;

        Eigen::Map<EigenDynamicVector> eigen_ref_snapshot(rReferenceSnapshot.data().begin(), rReferenceSnapshot.size());
        eigen_rx = eigen_ref_snapshot + eigen_phi_inf*eigen_rom_unknowns + eigen_phisig_sup*q_sup_aprox;
    }


    /**
     * @brief Get the derivative of the decoder over the input.
     * @param rRomUnknowns The latent vector to decode.
     * @param rPhiGlobal The PhiEffective matrix to be updated
     * @param rSVDPhiMatrices The list of projection matrices for the snapshot
     * @param r_W_mat The weights matrix for the RBF
     * @param r_centers_mat The collection of center vectors in matrix form
     * @param r_kernel_type Index indicating the kernel type
     * @param r_kernel_eps Epsilon value for the RBF's kernel
     */
    static void GetDecoderGradient(
        Vector rRomUnknowns,
        Matrix& rPhiGlobal,
        vector<Matrix>& rSVDPhiMatrices,
        Matrix& r_W_mat,
        Matrix& r_centers_mat,
        IndexType& r_kernel_type,
        double& r_kernel_eps
    )
    {
        Eigen::Map<EigenDynamicMatrix> eigen_phi_inf(rSVDPhiMatrices[0].data().begin(), rSVDPhiMatrices[0].size1(), rSVDPhiMatrices[0].size2());
        Eigen::Map<EigenDynamicMatrix> eigen_phisig_sup(rSVDPhiMatrices[1].data().begin(), rSVDPhiMatrices[1].size1(), rSVDPhiMatrices[1].size2());
        Eigen::Map<EigenDynamicMatrix> eigen_sig_inv_inf(rSVDPhiMatrices[2].data().begin(), rSVDPhiMatrices[2].size1(), rSVDPhiMatrices[2].size2());

        Eigen::Map<EigenDynamicVector> eigen_rom_unknowns(rRomUnknowns.data().begin(), rRomUnknowns.size());
        Eigen::Map<EigenDynamicMatrix> eigen_phi_global(rPhiGlobal.data().begin(), rPhiGlobal.size1(), rPhiGlobal.size2());

        EigenDynamicVector q_inf_pred = eigen_sig_inv_inf*eigen_rom_unknowns;
        Eigen::Map<EigenDynamicMatrix> distance_mat(r_centers_mat.data().begin(), r_centers_mat.size1(), r_centers_mat.size2());
        
        distance_mat.rowwise() -= q_inf_pred.transpose();
        distance_mat *= -1;

        EigenDynamicVector norms_vec = distance_mat.rowwise().norm();
        EigenDynamicVector aux_derivatives_vec;
        
        if (r_kernel_type == 0) {
            EigenDynamicVector kernel_results_vec = (-(r_kernel_eps * norms_vec.array()).square()).exp();
            aux_derivatives_vec = -2*std::pow(r_kernel_eps,2) * kernel_results_vec; // The derivative of the kernel would be missing multiplication by the norm,
                            // but then we would divide by the norm in the next step anyways as part of the derivative of the norm.
        } else if (r_kernel_type == 1) {
            EigenDynamicVector aux_term = 1.0 + (r_kernel_eps * norms_vec.array()).square();
            aux_derivatives_vec = -aux_term.array().pow(-3.0/2.0)*std::pow(r_kernel_eps,2); // Also skips the multiplication by the norm
        }
        
        Eigen::Map<EigenDynamicMatrix> eigen_W_mat(r_W_mat.data().begin(), r_W_mat.size1(), r_W_mat.size2());
        EigenDynamicMatrix rbf_gradient_aux =  aux_derivatives_vec.asDiagonal() * distance_mat;
        EigenDynamicMatrix rbf_gradient = eigen_W_mat.transpose() * rbf_gradient_aux;  // Shape is n_sup x n_inf

        eigen_phi_global=eigen_phi_inf+eigen_phisig_sup*rbf_gradient*eigen_sig_inv_inf;  // Shape is num_dofs x n_inf
    }

};

///@}

} // namespace Kratos

#endif // ROM_RBF_UTILITY_H_INCLUDED
