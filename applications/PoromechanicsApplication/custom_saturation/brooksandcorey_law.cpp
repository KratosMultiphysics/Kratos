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

    // If the capillar pressure is lower than the gas-entry pressure, the porous media is fully saturated with the wetting phase.
    rSl = 1.0;
    rdSldPc = 0.0;

    if(rVariables.pc > rVariables.pb)
    {
        // Liquid saturation degree
        rSl = (1.0 - rVariables.Sgr - rVariables.Slr)*std::pow(rVariables.pb/rVariables.pc,rVariables.lambda)
                + rVariables.Slr;

        // Derivative of the liquid saturation degree with respect to the capillary pressure
        rdSldPc = (1.0 - rVariables.Sgr - rVariables.Slr) * 
                    rVariables.lambda * rVariables.pb * std::pow(rVariables.pb/rVariables.pc,rVariables.lambda-1.0) /
                    (rVariables.pc * rVariables.pc);
    }
}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::CalculateLiquidRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrl = rValues.Getkrl();

    const double nw = (2.0 + 3.0*rVariables.lambda)/rVariables.lambda;

    rkrl = std::pow(rVariables.Se,nw);
    rkrl = std::max(rkrl,rVariables.krmin);
}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::CalculateGasRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrg = rValues.Getkrg();
    
    const double ng = (2.0 + rVariables.lambda)/rVariables.lambda;

    rkrg = std::pow(1.0-rVariables.Se,2.0)*(1.0 - std::pow(rVariables.Se,ng));
    rkrg = std::max(rkrg,rVariables.krmin);
}

}