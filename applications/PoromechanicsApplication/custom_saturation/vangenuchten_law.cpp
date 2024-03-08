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
#include "custom_saturation/vangenuchten_law.hpp"

namespace Kratos
{

void VanGenuchtenLaw::CalculateLiquidSaturationDegree (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rSl = rValues.GetSl();
    double& rdSldPc = rValues.GetdSldPc();

    // If the capillar pressure is lower than 0.0, the porous media is fully saturated with the wetting phase.
    rSl = 1.0 - rVariables.Sgr;
    rdSldPc = 0.0;

    if(rVariables.pc > 0.0)
    {
        // Liquid saturation degree
        rSl = (1.0 - rVariables.Sgr - rVariables.Slr)*std::pow(1.0 + std::pow(rVariables.pc/rVariables.pb,1.0/(1.0-rVariables.lambda)),-rVariables.lambda) 
                + rVariables.Slr;

        if (rSl <= rVariables.Slr) {
            rSl = rVariables.Slr;
        } else if (rSl >= 1.0 - rVariables.Sgr) {
            rSl = 1.0 - rVariables.Sgr;
        } else {
            // Derivative of the liquid saturation degree with respect to the capillary pressure
            rdSldPc = -rVariables.lambda * (1.0 - rVariables.Sgr - rVariables.Slr) *
                        std::pow(1.0 + std::pow(rVariables.pc/rVariables.pb,1.0/(1.0-rVariables.lambda)),-(rVariables.lambda+1.0)) *
                        std::pow(std::pow(rVariables.pc,rVariables.lambda)/rVariables.pb,1.0/(1.0-rVariables.lambda))/(1.0-rVariables.lambda);
        }
    }
}

//------------------------------------------------------------------------------------------------

void VanGenuchtenLaw::CalculateLiquidRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrl = rValues.Getkrl();

    const double nl = 1.0/rVariables.lambda;

    rkrl = std::sqrt(rVariables.Se)*std::pow(1.0-std::pow(1.0-std::pow(rVariables.Se,nl),rVariables.lambda),2.0);
    rkrl = std::max(rkrl,rVariables.krmin);
}

//------------------------------------------------------------------------------------------------

void VanGenuchtenLaw::CalculateGasRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrg = rValues.Getkrg();

    const double ng = 1.0/rVariables.lambda;
    
    rkrg = std::sqrt(1.0-rVariables.Se)*std::pow(1.0-std::pow(rVariables.Se,ng),2.0*rVariables.lambda);
    rkrg = std::max(rkrg,rVariables.krmin);
}

}