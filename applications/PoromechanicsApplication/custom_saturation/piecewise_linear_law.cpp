//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Danilo Cavalcanti
//                   Ignasi de Pouplana
//

// Application includes
#include "custom_saturation/piecewise_linear_law.hpp"

namespace Kratos
{

void PiecewiseLinearLaw::CalculateLiquidSaturationDegree (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rSl = rValues.GetSl();
    double& rdSldPc = rValues.GetdSldPc();

    const double Slmin = (rVariables.Slr > 0.0) ? rVariables.Slr : 0.2;
    const double Slmax = (rVariables.Sgr > 0.0) ? (1.0 - rVariables.Sgr) : 0.8;
    const double pc_min = 0.0;
    const double pc_max = (rVariables.pb > 0.0) ? rVariables.pb : 1.0;

    const double slope = -(Slmax - Slmin) / (pc_max - pc_min);

    if (rVariables.pc <= pc_min) {
        rSl = Slmax;
    } else if (rVariables.pc >= pc_max) {
        rSl = Slmin;
    } else {
        rSl = Slmax + slope * (rVariables.pc - pc_min);
    }

    // Derivative: PorousLab UMAT GetCurveDerivative preserves slope everywhere
    rdSldPc = slope;
}

//------------------------------------------------------------------------------------------------

void PiecewiseLinearLaw::CalculateLiquidRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrl = rValues.Getkrl();

    this->CalculateEffectiveSaturation(rVariables, rValues);

    const double krmin = (rVariables.krmin > 0.0) ? rVariables.krmin : 1.0e-4;

    if (rVariables.Se >= 1.0) {
        rkrl = 1.0;
    } else if (rVariables.Se <= 0.0) {
        rkrl = krmin;
    } else {
        const double lambda = (rVariables.lambda > 0.0) ? rVariables.lambda : 2.0;
        const double nl = (2.0 + 3.0 * lambda) / lambda;
        rkrl = std::pow(rVariables.Se, nl);
        rkrl = std::max(rkrl, krmin);
    }
}

//------------------------------------------------------------------------------------------------

void PiecewiseLinearLaw::CalculateGasRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrg = rValues.Getkrg();

    this->CalculateEffectiveSaturation(rVariables, rValues);

    const double krmin = (rVariables.krmin > 0.0) ? rVariables.krmin : 1.0e-4;

    if (rVariables.Se >= 1.0) {
        rkrg = krmin;
    } else if (rVariables.Se <= 0.0) {
        rkrg = 1.0;
    } else {
        const double lambda = (rVariables.lambda > 0.0) ? rVariables.lambda : 2.0;
        const double ng = (2.0 + lambda) / lambda;
        rkrg = (1.0 - rVariables.Se) * (1.0 - rVariables.Se) * (1.0 - std::pow(rVariables.Se, ng));
        rkrg = std::max(rkrg, krmin);
    }
}

} // namespace Kratos

