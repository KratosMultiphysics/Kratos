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
#include <Eigen/Core>
#include <Eigen/Dense>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

// Application includes
#include "rom_application_variables.h"


namespace Kratos
{

/**
 * @class RomRBFUtility
 * @ingroup RomApplication
 * @brief Decoder of the RBF-enhanced ROMs: x = x_ref + Phi_inf*q + PhiSig_sup*W^T*k(SigInv_inf*q)
 * @details k contains the kernel evaluated at the distances from the (scaled) reduced coordinates to each of the centers.
 * The kernel types are 0: gaussian, exp(-(eps*r)^2), and 1: inverse multiquadric, 1/sqrt(1+(eps*r)^2)
 */
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
     * @param rWeights The weights matrix of the RBF (number of centers x number of superior modes)
     * @param rCenters The centers of the RBF, one per row
     * @param KernelType Index indicating the kernel type
     * @param KernelEps Epsilon value of the kernel
     * @param rReferenceSnapshot Reference snapshot to sum within the decoder
     */
    static void GetXAndDecoderGradient(
        const Vector& rRomUnknowns,
        Vector& rx,
        Matrix& rPhiGlobal,
        const vector<Matrix>& rSVDPhiMatrices,
        const Matrix& rWeights,
        const Matrix& rCenters,
        const IndexType KernelType,
        const double KernelEps,
        const Vector& rReferenceSnapshot
    )
    {
        GetXFromDecoder(rRomUnknowns, rx, rSVDPhiMatrices, rWeights, rCenters, KernelType, KernelEps, rReferenceSnapshot);
        GetDecoderGradient(rRomUnknowns, rPhiGlobal, rSVDPhiMatrices, rWeights, rCenters, KernelType, KernelEps);
    }

    /**
     * @brief Run the decoder on the input latent vector.
     * @param rRomUnknowns The latent vector to decode.
     * @param rx The resulting snapshot vector.
     * @param rSVDPhiMatrices The list of projection matrices for the snapshot
     * @param rWeights The weights matrix of the RBF (number of centers x number of superior modes)
     * @param rCenters The centers of the RBF, one per row
     * @param KernelType Index indicating the kernel type
     * @param KernelEps Epsilon value of the kernel
     * @param rReferenceSnapshot Reference snapshot to sum within the decoder
     */
    static void GetXFromDecoder(
        const Vector& rRomUnknowns,
        Vector& rx,
        const vector<Matrix>& rSVDPhiMatrices,
        const Matrix& rWeights,
        const Matrix& rCenters,
        const IndexType KernelType,
        const double KernelEps,
        const Vector& rReferenceSnapshot
        )
    {
        Eigen::Map<const EigenDynamicMatrix> eigen_phi_inf(rSVDPhiMatrices[0].data().begin(), rSVDPhiMatrices[0].size1(), rSVDPhiMatrices[0].size2());
        Eigen::Map<const EigenDynamicMatrix> eigen_phisig_sup(rSVDPhiMatrices[1].data().begin(), rSVDPhiMatrices[1].size1(), rSVDPhiMatrices[1].size2());
        Eigen::Map<const EigenDynamicMatrix> eigen_weights(rWeights.data().begin(), rWeights.size1(), rWeights.size2());
        Eigen::Map<const EigenDynamicVector> eigen_rom_unknowns(rRomUnknowns.data().begin(), rRomUnknowns.size());
        Eigen::Map<const EigenDynamicVector> eigen_ref_snapshot(rReferenceSnapshot.data().begin(), rReferenceSnapshot.size());
        Eigen::Map<EigenDynamicVector> eigen_rx(rx.data().begin(), rx.size());

        EigenDynamicMatrix differences;
        EigenDynamicVector kernels;
        EigenDynamicVector kernel_derivative_factors;
        EvaluateKernels(rRomUnknowns, rSVDPhiMatrices, rCenters, KernelType, KernelEps, differences, kernels, kernel_derivative_factors);

        const EigenDynamicVector q_sup = eigen_weights.transpose() * kernels;
        eigen_rx = eigen_ref_snapshot + eigen_phi_inf*eigen_rom_unknowns + eigen_phisig_sup*q_sup;
    }

    /**
     * @brief Get the derivative of the decoder over the input.
     * @param rRomUnknowns The latent vector to decode.
     * @param rPhiGlobal The PhiEffective matrix to be updated
     * @param rSVDPhiMatrices The list of projection matrices for the snapshot
     * @param rWeights The weights matrix of the RBF (number of centers x number of superior modes)
     * @param rCenters The centers of the RBF, one per row
     * @param KernelType Index indicating the kernel type
     * @param KernelEps Epsilon value of the kernel
     */
    static void GetDecoderGradient(
        const Vector& rRomUnknowns,
        Matrix& rPhiGlobal,
        const vector<Matrix>& rSVDPhiMatrices,
        const Matrix& rWeights,
        const Matrix& rCenters,
        const IndexType KernelType,
        const double KernelEps
    )
    {
        Eigen::Map<const EigenDynamicMatrix> eigen_phi_inf(rSVDPhiMatrices[0].data().begin(), rSVDPhiMatrices[0].size1(), rSVDPhiMatrices[0].size2());
        Eigen::Map<const EigenDynamicMatrix> eigen_phisig_sup(rSVDPhiMatrices[1].data().begin(), rSVDPhiMatrices[1].size1(), rSVDPhiMatrices[1].size2());
        Eigen::Map<const EigenDynamicMatrix> eigen_sig_inv_inf(rSVDPhiMatrices[2].data().begin(), rSVDPhiMatrices[2].size1(), rSVDPhiMatrices[2].size2());
        Eigen::Map<const EigenDynamicMatrix> eigen_weights(rWeights.data().begin(), rWeights.size1(), rWeights.size2());
        Eigen::Map<EigenDynamicMatrix> eigen_phi_global(rPhiGlobal.data().begin(), rPhiGlobal.size1(), rPhiGlobal.size2());

        EigenDynamicMatrix differences;
        EigenDynamicVector kernels;
        EigenDynamicVector kernel_derivative_factors;
        EvaluateKernels(rRomUnknowns, rSVDPhiMatrices, rCenters, KernelType, KernelEps, differences, kernels, kernel_derivative_factors);

        // Gradient of the RBF with respect to the scaled reduced coordinates (number of superior modes x number of inferior modes)
        const EigenDynamicMatrix rbf_gradient = eigen_weights.transpose() * (kernel_derivative_factors.asDiagonal() * differences);

        eigen_phi_global = eigen_phi_inf + eigen_phisig_sup*rbf_gradient*eigen_sig_inv_inf;
    }

private:

    /**
     * @brief Evaluates the kernel of each center at the (scaled) reduced coordinates.
     * @param rDifferences The scaled reduced coordinates minus each center, one per row
     * @param rKernels The value of the kernel of each center
     * @param rKernelDerivativeFactors Factors such that the gradient of the kernel of the i-th center is rKernelDerivativeFactors[i]*rDifferences.row(i)
     */
    static void EvaluateKernels(
        const Vector& rRomUnknowns,
        const vector<Matrix>& rSVDPhiMatrices,
        const Matrix& rCenters,
        const IndexType KernelType,
        const double KernelEps,
        EigenDynamicMatrix& rDifferences,
        EigenDynamicVector& rKernels,
        EigenDynamicVector& rKernelDerivativeFactors
    )
    {
        Eigen::Map<const EigenDynamicMatrix> eigen_sig_inv_inf(rSVDPhiMatrices[2].data().begin(), rSVDPhiMatrices[2].size1(), rSVDPhiMatrices[2].size2());
        Eigen::Map<const EigenDynamicVector> eigen_rom_unknowns(rRomUnknowns.data().begin(), rRomUnknowns.size());
        Eigen::Map<const EigenDynamicMatrix> eigen_centers(rCenters.data().begin(), rCenters.size1(), rCenters.size2());

        const EigenDynamicVector scaled_rom_unknowns = eigen_sig_inv_inf*eigen_rom_unknowns;

        // Note that the centers are not modified
        rDifferences = (-eigen_centers).rowwise() + scaled_rom_unknowns.transpose();
        const EigenDynamicVector squared_scaled_distances = std::pow(KernelEps, 2) * rDifferences.rowwise().squaredNorm();

        if (KernelType == 0) { // Gaussian
            rKernels = (-squared_scaled_distances.array()).exp();
            rKernelDerivativeFactors = -2.0*std::pow(KernelEps, 2) * rKernels;
        } else if (KernelType == 1) { // Inverse multiquadric
            rKernels = (1.0 + squared_scaled_distances.array()).rsqrt();
            rKernelDerivativeFactors = -std::pow(KernelEps, 2) * rKernels.array().cube();
        } else {
            KRATOS_ERROR << "Unknown RBF kernel type " << KernelType << ". Available options are 0 (gaussian) and 1 (inverse multiquadric)." << std::endl;
        }
    }

};

///@}

} // namespace Kratos

#endif // ROM_RBF_UTILITY_H_INCLUDED
