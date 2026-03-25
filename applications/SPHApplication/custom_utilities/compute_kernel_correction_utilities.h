
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

    static void ComputeGradientCorrection(ModelPart& rThisModelPart);

    static void ComputeIntegrationCorrectionVector(ModelPart& rThisModelPart);

    static void CalculateLinearSystemForCorrection(
        ModelPart& rThisModelPart, 
        MatrixType& rLHS,
        VectorType& rRHS
    );

    static bool VerifyKernelCorrection(ModelPart& rThisModelPart, Parameters& rThisParameters);

    static void ApplyKernelCorrection(Element& IP, double& kernel_target);
    
    static void ApplyKernelGradientCorrection(Element& IP, double& kernel_target, Vector& dkernel_target);

};

}