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

#if !defined(ROM_NN_UTILITY_H_INCLUDED)
#define ROM_NN_UTILITY_H_INCLUDED

// System includes
#include <unordered_map>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

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
class RomNNUtility
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
     * @param rNNLayers The neural network's dense layers weights (no bias)
     * @param rReferenceSnapshot Reference snapshot to sum within the decoder
     */
    static void GetXAndDecoderGradient(
        Vector rRomUnknowns,
        Vector& rx,
        Matrix& rPhiGlobal,
        vector<Matrix>& rSVDPhiMatrices,
        vector<Matrix>& rNNLayers,
        Vector& rReferenceSnapshot
    )
    {
        Eigen::Map<EigenDynamicMatrix> eigen_phi_inf(rSVDPhiMatrices[0].data().begin(), rSVDPhiMatrices[0].size1(), rSVDPhiMatrices[0].size2());
        Eigen::Map<EigenDynamicMatrix> eigen_phisig_sup(rSVDPhiMatrices[1].data().begin(), rSVDPhiMatrices[1].size1(), rSVDPhiMatrices[1].size2());
        Eigen::Map<EigenDynamicMatrix> eigen_sig_inv_inf(rSVDPhiMatrices[2].data().begin(), rSVDPhiMatrices[2].size1(), rSVDPhiMatrices[2].size2());

        Eigen::Map<EigenDynamicVector> eigen_rom_unknowns(rRomUnknowns.data().begin(), rRomUnknowns.size());
        Eigen::Map<EigenDynamicMatrix> eigen_phi_global(rPhiGlobal.data().begin(), rPhiGlobal.size1(), rPhiGlobal.size2());
        Eigen::Map<EigenDynamicVector> eigen_rx(rx.data().begin(), rx.size());

        vector<EigenDynamicMatrix> mGradientLayers(rNNLayers.size());
        
        EigenDynamicVector layerOutTemp = eigen_sig_inv_inf*eigen_rom_unknowns;
        EigenDynamicVector layerOut = layerOutTemp;

        SizeType i=0;
        for(auto& layers_item : rNNLayers)
        {   
            Eigen::Map<EigenDynamicMatrix> eigen_layer_item(layers_item.data().begin(), layers_item.size1(), layers_item.size2());
            layerOutTemp = eigen_layer_item.transpose()*layerOut;
            layerOut = layerOutTemp;

            mGradientLayers[i] = eigen_layer_item;

            if(i<rNNLayers.size()-1){
                IndexPartition<IndexType>(layerOutTemp.size()).for_each([&](IndexType Index)
                {   
                    if (layerOutTemp[Index]<0.0){
                        mGradientLayers[i].col(Index)*=exp(layerOutTemp[Index]);
                        layerOut[Index] = exp(layerOutTemp[Index])-1.0;
                    }
                });
            }
            i++;
        }

        EigenDynamicMatrix intermediateGradientsTemp = mGradientLayers[0];
        EigenDynamicMatrix intermediateGradients = intermediateGradientsTemp;

        for (SizeType i = 0; i < mGradientLayers.size()-1; i++)
        {
            intermediateGradientsTemp = intermediateGradients*mGradientLayers[i+1];
            intermediateGradients = intermediateGradientsTemp;
        }

        eigen_phi_global=eigen_phi_inf+eigen_phisig_sup*intermediateGradients.transpose()*eigen_sig_inv_inf;

        Eigen::Map<EigenDynamicVector> eigen_ref_snapshot(rReferenceSnapshot.data().begin(), rReferenceSnapshot.size());
        eigen_rx = eigen_ref_snapshot + eigen_phi_inf*eigen_rom_unknowns + eigen_phisig_sup*layerOut;
    }

    /**
     * @brief Run the decoder on the input latent vector.
     * @param rRomUnknowns The latent vector to decode.
     * @param rx The resulting snapshot vector.
     * @param rSVDPhiMatrices The list of projection matrices for the snapshot
     * @param rNNLayers The neural network's dense layers weights (no bias)
     * @param rReferenceSnapshot Reference snapshot to sum within the decoder
     */
    static void GetXFromDecoder(
        Vector rRomUnknowns,
        Vector& rx,
        vector<Matrix>& rSVDPhiMatrices,
        vector<Matrix>& rNNLayers,
        Vector& rReferenceSnapshot
        )
    {
        Eigen::Map<EigenDynamicMatrix> eigen_phi_inf(rSVDPhiMatrices[0].data().begin(), rSVDPhiMatrices[0].size1(), rSVDPhiMatrices[0].size2());
        Eigen::Map<EigenDynamicMatrix> eigen_phisig_sup(rSVDPhiMatrices[1].data().begin(), rSVDPhiMatrices[1].size1(), rSVDPhiMatrices[1].size2());
        Eigen::Map<EigenDynamicMatrix> eigen_sig_inv_inf(rSVDPhiMatrices[2].data().begin(), rSVDPhiMatrices[2].size1(), rSVDPhiMatrices[2].size2());

        Eigen::Map<EigenDynamicVector> eigen_rom_unknowns(rRomUnknowns.data().begin(), rRomUnknowns.size());
        Eigen::Map<EigenDynamicVector> eigen_rx(rx.data().begin(), rx.size());
        
        EigenDynamicVector layerOutTemp = eigen_sig_inv_inf*eigen_rom_unknowns;
        EigenDynamicVector layerOut = layerOutTemp;

        SizeType i=0;
        for(auto& layers_item : rNNLayers)
        {   
            Eigen::Map<EigenDynamicMatrix> eigen_layer_item(layers_item.data().begin(), layers_item.size1(), layers_item.size2());
            layerOutTemp = eigen_layer_item.transpose()*layerOut;
            layerOut = layerOutTemp;

            if(i<rNNLayers.size()-1){
                IndexPartition<IndexType>(layerOutTemp.size()).for_each([&](IndexType Index)
                {   
                    if (layerOutTemp[Index]<0.0){
                        layerOut[Index] = exp(layerOutTemp[Index])-1.0;
                    }
                });
            }
            i++;
        }

        Eigen::Map<EigenDynamicVector> eigen_ref_snapshot(rReferenceSnapshot.data().begin(), rReferenceSnapshot.size());
        eigen_rx = eigen_ref_snapshot + eigen_phi_inf*eigen_rom_unknowns + eigen_phisig_sup*layerOut;
    }


    /**
     * @brief Get the derivative of the decoder over the input.
     * @param rRomUnknowns The latent vector to decode.
     * @param rPhiGlobal The PhiEffective matrix to be updated
     * @param rSVDPhiMatrices The list of projection matrices for the snapshot
     * @param rNNLayers The neural network's dense layers weights (no bias)
     */
    static void GetDecoderGradient(
        Vector rRomUnknowns,
        Matrix& rPhiGlobal,
        vector<Matrix>& rSVDPhiMatrices,
        vector<Matrix>& rNNLayers
    )
    {
        Eigen::Map<EigenDynamicMatrix> eigen_phi_inf(rSVDPhiMatrices[0].data().begin(), rSVDPhiMatrices[0].size1(), rSVDPhiMatrices[0].size2());
        Eigen::Map<EigenDynamicMatrix> eigen_phisig_sup(rSVDPhiMatrices[1].data().begin(), rSVDPhiMatrices[1].size1(), rSVDPhiMatrices[1].size2());
        Eigen::Map<EigenDynamicMatrix> eigen_sig_inv_inf(rSVDPhiMatrices[2].data().begin(), rSVDPhiMatrices[2].size1(), rSVDPhiMatrices[2].size2());

        Eigen::Map<EigenDynamicVector> eigen_rom_unknowns(rRomUnknowns.data().begin(), rRomUnknowns.size());
        Eigen::Map<EigenDynamicMatrix> eigen_phi_global(rPhiGlobal.data().begin(), rPhiGlobal.size1(), rPhiGlobal.size2());

        vector<EigenDynamicMatrix> mGradientLayers(rNNLayers.size());
        
        EigenDynamicVector layerOutTemp = eigen_sig_inv_inf*eigen_rom_unknowns;
        EigenDynamicVector layerOut = layerOutTemp;

        SizeType i=0;
        for(auto& layers_item : rNNLayers)
        {   
            Eigen::Map<EigenDynamicMatrix> eigen_layer_item(layers_item.data().begin(), layers_item.size1(), layers_item.size2());
            layerOutTemp = eigen_layer_item.transpose()*layerOut;
            layerOut = layerOutTemp;

            mGradientLayers[i] = eigen_layer_item;

            if(i<rNNLayers.size()-1){
                IndexPartition<IndexType>(layerOutTemp.size()).for_each([&](IndexType Index)
                {   
                    if (layerOutTemp[Index]<0.0){
                        mGradientLayers[i].col(Index)*=exp(layerOutTemp[Index]);
                        layerOut[Index] = exp(layerOutTemp[Index])-1.0;
                    }
                });
            }
            i++;
        }

        EigenDynamicMatrix intermediateGradientsTemp = mGradientLayers[0];
        EigenDynamicMatrix intermediateGradients = intermediateGradientsTemp;

        for (SizeType i = 0; i < mGradientLayers.size()-1; i++)
        {
            intermediateGradientsTemp = intermediateGradients*mGradientLayers[i+1];
            intermediateGradients = intermediateGradientsTemp;
        }

        eigen_phi_global=eigen_phi_inf+eigen_phisig_sup*intermediateGradients.transpose()*eigen_sig_inv_inf;
    }

};

///@}

} // namespace Kratos

#endif // ROM_NN_UTILITY_H_INCLUDED
