#pragma once 

#include "custom_utilities/custom_kernels/cubic_kernel_3d.h"

namespace Kratos
{

void CubicKernel3D::ComputeKernelValue(double& kernel, const double h, const Vector& X_AB_target)
{ 
    kernel = 0.0;

    const double dist = norm_2(X_AB_target);
    double q = dist / h;
    double ConstFactor = 3.0 / (2.0 * Globals::Pi * std::pow(h, 3.0));

    if (q <= 1.0){
        kernel = ConstFactor * (2.0/3.0 - q * q + 0.5 * q * q * q);
    } else if (q > 1.0 && q <= 2.0){
        double t = 2.0 - q;
        kernel = ConstFactor * (1.0/6.0) * (t * t * t);
    }
}

void CubicKernel3D::ComputeKernelGradientValue(Vector& dkernel, const double h, const Vector& X_AB_target)
{
    noalias(dkernel) = ZeroVector(dkernel.size());
    
    const double dist = norm_2(X_AB_target);
    double q = dist / h;
    double ConstFactor = 3.0 / (2.0 * Globals::Pi * std::pow(h, 3.0));
    Vector VectorFactor = X_AB_target / (h * dist);

    if (q > 0.0 && q <= 1.0){
        dkernel = ConstFactor * (- 2.0 * q + 1.5 * q * q) * VectorFactor;
    } else if (q > 1.0 && q <= 2.0) {
        double t = 2 - q; 
        dkernel = - ConstFactor * 0.5 * t * t * VectorFactor;
    }
}

void CubicKernel3D::ComputeKernelLaplacianValue()
{
}

}