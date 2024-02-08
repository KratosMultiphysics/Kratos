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

void BrooksAndCoreyLaw::CalculateWaterSaturationDegree (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rSw = rValues.GetSw();
    double& rdSwdPc = rValues.GetdSwdPc();

    // If the capillar pressure is lower than the gas-entry pressure, the porous media is fully saturated with the wetting phase.
    rSw = 1.0;
    rdSwdPc = 0.0;

    if(rVariables.pc > rVariables.pb)
    {
        // Water saturation degree
        rSw = (1.0 - rVariables.Swr)*std::pow(rVariables.pb/rVariables.pc,rVariables.lambda)
                + rVariables.Swr;

        // Derivative of the water saturation degree with respect to the capilar pressure
        rdSwdPc = (1.0 - rVariables.Swr) * 
                    rVariables.lambda * rVariables.pb * std::pow(rVariables.pb/rVariables.pc,rVariables.lambda-1.0) /
                    (rVariables.pc * rVariables.pc);
    }
}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::CalculateWaterRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrw = rValues.Getkrw();

    const double nw = (2.0 + 3.0*rVariables.lambda)/rVariables.lambda;

    rkrw = std::pow(rVariables.Se,nw);
    // TODO. Review this number
    rkrw = std::min(rkrw,0.0001);
}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::CalculateGasRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrg = rValues.Getkrg();
    
    const double ng = (2.0 + rVariables.lambda)/rVariables.lambda;

    rkrg = std::pow(1.0-rVariables.Se,2.0)*(1.0 - std::pow(rVariables.Se,ng));
    // TODO. Review this number
    rkrg = std::min(rkrg,0.0001);
}

}