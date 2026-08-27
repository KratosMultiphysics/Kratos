#pragma once 

#include "custom_utilities/custom_kernels/cubic_kernel_2d.h"


namespace Kratos
{

void CubicKernel2D::ComputeKernelValue(double& kernel, const double h, const Vector& X_AB_target)
{
    kernel = 0.0;
    
    const double dist = norm_2(X_AB_target);
    double q = dist / h;
    double ConstFactor = 15.0 / (7.0 * Globals::Pi * h * h);

    if (q <= 1.0){
        kernel = ConstFactor * (2.0/3.0 - q * q + 0.5 * q * q * q);
    } else if (q > 1.0 && q <= 2.0){
        double t = 2.0 - q;
        kernel = ConstFactor * (1.0/6.0) * (t * t * t);
    }
}

void CubicKernel2D::ComputeKernelGradientValue(Vector& dkernel, const double h, const Vector& X_AB_target)
{
    noalias(dkernel) = ZeroVector(dkernel.size());
    
    const double dist = norm_2(X_AB_target);
    double q = dist / h;
    double ConstFactor = 15.0 / (7.0 * Globals::Pi * h * h);
    Vector VectorFactor = X_AB_target / (h * dist);

    if (q > 0.0 && q <= 1.0){
        dkernel = ConstFactor * (- 2.0 * q + 1.5 * q * q) * VectorFactor;
    } else if (q > 1.0 && q <= 2.0) {
        double t = 2 - q; 
        dkernel = - ConstFactor * 0.5 * t * t * VectorFactor;
    }
}

void CubicKernel2D::ComputeKernelLaplacianValue()
{
}

}