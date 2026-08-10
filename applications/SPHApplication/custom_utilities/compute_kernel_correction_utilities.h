//  ____  ____  _   _                   _ _           _   _             
// / ___||  _ \| | | | __ _ _ __  _ __ | (_) ___ __ _| |_(_) ___  _ __  
// \___ \| |_) | |_| |/ _` | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_ \ 
//  ___) |  __/|  _  | (_| | |_) | |_) | | | (_| (_| | |_| | (_) | | | |
// |____/|_|   |_| |_|\__,_| .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                         |_|   |_|                                    

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once 

#include <cmath>
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "utilities/parallel_utilities.h"
#include "sph_application_variables.h"
#include "custom_utilities/custom_kernels/kernel_factory.h"

/**
 * @class ComputeKernelCorrectionUtilities
 * @brief This class provides utility functions to compute the kernel correction 
 * and respective checks for the SPH method
 * @details The methods are static, so it can be called without constructing the class
 */

namespace Kratos
{
class ComputeKernelCorrectionUtilities
{
public:

    using SizeType = std::size_t;
    using VectorType = Vector;
    using MatrixType = Matrix;

    /**
     * @brief This function computes the integration correction which ensure fisrt-order consistency in the boundaries of the domain.
     * @details J. Bonet and T.S.L. Lok "Variational and momentum preservation aspects of Smooth Particle Hydrodynamic formulations"
     */
    static void ComputeWeightedSums(ModelPart& rThisModelPart);

    static void ComputeGradientCorrection(ModelPart& rThisModelPart);
    
    /**
     * @brief This function computes the integration correction which ensure fisrt-order consistency in the boundaries of the domain.
     * @details J. Bonet and S. Kulasegaram "Correction and stabilization of smooth particle 
     * hydrodynamics methods with applications  in metal forming simulations"
     */
    static void ComputeIntegrationCorrection(ModelPart& rThisModelPart, Parameters& rThisParameters, unsigned int& rIter);

    /**
     * @brief This functions applies the gradient and the kernel corrections
     */
    static void ApplyKernelCorrection(Element& IP, double& kernel_target);
    
    static void ApplyKernelGradientCorrection(Element& rIntegrationParticle, double& rKernel, Vector& rKernelGradient);

    static void ApplyKernelGradientCorrectionInverted(Element& rNeighbouringParticle, double& kernel_target, Vector& dkernel_target);

    static void ApplyIntegrationCorrection(Element& IP, double& kernel_target, Vector& dkernel_target, bool IsParticleItself);

    /**
     * @brief These functions check the effectiveness of the kernel corrections
     */
    static bool VerifyKernelCorrection(ModelPart& rThisModelPart, Parameters& rThisParameters);

    static bool VerifyIntegrationCorrection(ModelPart& rThisModelPart, Parameters& rThisParameters);

};

}
