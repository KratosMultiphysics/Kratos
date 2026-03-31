#if !defined(KRATOS_BOUSSINESQ_COUPLING_UTILITIES_H_INCLUDED)
#define KRATOS_BOUSSINESQ_COUPLING_UTILITIES_H_INCLUDED

#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{

class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) BoussinesqCouplingUtilities
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BoussinesqCouplingUtilities);

    /// Computes the relative L2 residual between a historical variable and a non-historical variable
    /// integrated over the domains elements.
    static double ComputeRelativeResidual(
        ModelPart& rModelPart,
        const Variable<double>& rVarNew,
        const Variable<double>& rVarOld);

    /// Applies relaxation to the nodal values using OpenMP
    static void ApplyRelaxation(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const Variable<double>& rOldVariable,
        const double Omega);

    /// Compute the x and residual vectors for the convergence accelerator (OpenMP)
    static void ComputeQuasiNewtonUpdateVectors(
        ModelPart& rModelPart,
        const Variable<double>& rVarNew,
        const Variable<double>& rVarOld,
        Vector& rX,
        Vector& rR);

    /// Write back the updated values from the quasi-newton step into the modelpart (OpenMP)
    static void UpdateConvergenceVariables(
        ModelPart& rModelPart,
        const Variable<double>& rVarNew,
        const Variable<double>& rVarOld,
        const Vector& rX);
};

} // namespace Kratos

#endif // KRATOS_BOUSSINESQ_COUPLING_UTILITIES_H_INCLUDED
