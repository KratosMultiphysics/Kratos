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
#include "custom_saturation/liakopoulos_law.hpp"

namespace Kratos
{

void LiakopoulosLaw::CalculateLiquidSaturationDegree (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rSl = rValues.GetSl();
    double& rdSldPc = rValues.GetdSldPc();

    const double a_ret = 1.9722e-11;
    const double b_ret = 2.4279;
    const double Slmin = (rVariables.Slr > 0.0) ? rVariables.Slr : 0.2;

    if (rVariables.pc <= 0.0) {
        rSl = 1.0;
        rdSldPc = 0.0;
    } else {
        rSl = 1.0 - a_ret * std::pow(rVariables.pc, b_ret);
        rSl = std::max(rSl, Slmin);

        const double pc_max = std::pow((1.0 - Slmin) / a_ret, 1.0 / b_ret);
        if (rVariables.pc >= pc_max) {
            rdSldPc = 0.0;
        } else {
            rdSldPc = -a_ret * b_ret * std::pow(rVariables.pc, b_ret - 1.0);
        }
    }
}

//------------------------------------------------------------------------------------------------

void LiakopoulosLaw::CalculateLiquidRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrl = rValues.Getkrl();
    const double& rSl = rValues.GetSl();

    const double a_perm = 2.207;
    const double b_perm = 1.0121;
    const double Slmin = (rVariables.Slr > 0.0) ? rVariables.Slr : 0.2;

    if (rSl <= Slmin) {
        rkrl = rVariables.krmin;
    } else if (rSl >= 1.0) {
        rkrl = 1.0;
    } else {
        rkrl = 1.0 - a_perm * std::pow(1.0 - rSl, b_perm);
        rkrl = std::max(rkrl, rVariables.krmin);
    }
}

//------------------------------------------------------------------------------------------------

void LiakopoulosLaw::CalculateGasRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrg = rValues.Getkrg();

    this->CalculateEffectiveSaturation(rVariables, rValues);

    const double krmin = (rVariables.krmin > 0.0) ? rVariables.krmin : 1.0e-4;

    if (rVariables.Se >= 1.0) {
        rkrg = krmin;
    } else if (rVariables.Se <= 0.0) {
        rkrg = 1.0;
    } else {
        const double lambda = (rVariables.lambda > 0.0) ? rVariables.lambda : 3.0;
        const double ng = (2.0 + lambda) / lambda;
        rkrg = (1.0 - rVariables.Se) * (1.0 - rVariables.Se) * (1.0 - std::pow(rVariables.Se, ng));
        rkrg = std::max(rkrg, krmin);
    }
}

} // namespace Kratos
