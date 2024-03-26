//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Danilo Cavalcanti
//

/* Project includes */
#include "custom_saturation/brooksandcorey_law.hpp"

namespace Kratos
{

void BrooksAndCoreyLaw::CalculateLiquidSaturationDegree (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rSl = rValues.GetSl();
    double& rdSldPc = rValues.GetdSldPc();

    //TODO. Ignasi
    // This is only used in the Liakopoulos test
    // rSl = 1.0 - 0.10152*std::pow(rVariables.pc/(9806.0),2.4279);
    //

    // If the capillar pressure is lower than the gas-entry pressure, the porous media is fully saturated with the wetting phase.
    rSl = 1.0 - rVariables.Sgr;
    rdSldPc = 0.0;

    if(rVariables.pc > rVariables.pb)
    {
        // Liquid saturation degree
        rSl = (1.0 - rVariables.Sgr - rVariables.Slr)*std::pow(rVariables.pb/rVariables.pc,rVariables.lambda)
                + rVariables.Slr;

        // Derivative of the liquid saturation degree with respect to the capillary pressure
        rdSldPc = -rVariables.lambda * (1.0 - rVariables.Sgr - rVariables.Slr) * 
                    std::pow(rVariables.pb,rVariables.lambda) / std::pow(rVariables.pc,rVariables.lambda+1.0);
    }
}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::CalculateLiquidRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrl = rValues.Getkrl();

    //TODO. Ignasi
    // This is only used in the Liakopoulos test
    // double& rSl = rValues.GetSl();
    // rkrl = 1.0 - 2.207*std::pow(1.0-rSl,1.0121);
    //

    if (rVariables.Se >= 1.0) {
        // Fully saturated medium
        rkrl = 1.0;
    } else if (rVariables.Se <= 0.0) {
        // Dry medium
        rkrl = rVariables.krmin;
    } else {
        const double nl = (2.0 + 3.0*rVariables.lambda)/rVariables.lambda;
        rkrl = std::pow(rVariables.Se,nl);
        rkrl = std::max(rkrl,rVariables.krmin);
    }
}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::CalculateGasRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrg = rValues.Getkrg();

    if (rVariables.Se >= 1.0) {
        // Fully saturated medium
        rkrg = rVariables.krmin;
    } else if (rVariables.Se <= 0.0) {
        // Dry medium
        rkrg = 1.0;
    } else {
        const double ng = (2.0 + rVariables.lambda)/rVariables.lambda;
        rkrg = (1.0-rVariables.Se)*(1.0-rVariables.Se)*(1.0 - std::pow(rVariables.Se,ng));
        rkrg = std::max(rkrg,rVariables.krmin);
    }
}

}