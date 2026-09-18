#pragma once 

#include "includes/ublas_interface.h" // to include the boost interface
#include "includes/global_variables.h"

/**
 * @class CubicKernel3D
 * @brief This class implements a kernel function 2D used in SPH simulations
 */

namespace Kratos

{
class CubicKernel3D
{

public:

    static void ComputeKernelValue(double& kernel, const double h, const Vector& X_AB_target);

    static void ComputeKernelGradientValue(Vector& dkernel, const double h, const Vector& X_AB_target);

    static void ComputeKernelLaplacianValue();
};
    

}