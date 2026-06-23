
#pragma once 

#include <cmath>
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "sph_application_variables.h"
#include "custom_utilities/custom_kernels/kernel_factory.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"

/**
 * @class ComputeKernelCorrectionUtilities
 * @brief 
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

    static void ComputeWeightedSums(ModelPart& rThisModelPart);

    /**
     * @brief This function computes the integration correction which ensure fisrt-order consistency in the boundaries of the domain.
     * @details J. Bonet and T.S.L. Lok "Variational and momentum preservation aspects of Smooth Particle Hydrodynamic formulations"
     */
    static void ComputeGradientCorrection(ModelPart& rThisModelPart);
    
    /**
     * @brief This function computes the integration correction which ensure fisrt-order consistency in the boundaries of the domain.
     * @details J. Bonet and S. Kulasegaram "Correction and stabilization of smooth particle 
     * hydrodynamics methods with applications  in metal forming simulations"
     */
    static void ComputeIntegrationCorrection(ModelPart& rThisModelPart, Parameters& rThisParameters, unsigned int& iter);

    /**
     * @brief This function applies the kernel correction
     */
    static void ApplyKernelCorrection(Element& IP, double& kernel_target);
    
    /**
     * @brief This function applies the gradient and the kernel corrections
     */
    static void ApplyKernelGradientCorrection(Element& IP, double& kernel_target, Vector& dkernel_target);

    static void ApplyKernelGradientCorrectionInverted(Element& JP, double& kernel_target, Vector& dkernel_target);

    /**
     * @brief This function applies the integration correction
     */
    static void ApplyIntegrationCorrection(Element& IP, double& kernel_target, Vector& dkernel_target, bool IsParticleItself);

    /**
     * @brief These functions check the effectiveness of the kernel corrections
     */
    static bool VerifyKernelCorrection(ModelPart& rThisModelPart, Parameters& rThisParameters);

    static bool VerifyIntegrationCorrection(ModelPart& rThisModelPart, Parameters& rThisParameters);

};

}