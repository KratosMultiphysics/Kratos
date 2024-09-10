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
    // rSl = 1.0 - rVariables.Sgr;
    // rdSldPc = 0.0;
    
    // if(rVariables.pc > rVariables.pb)
    // {
    //     Liquid saturation degree
    //     rSl = (1.0 - rVariables.Sgr - rVariables.Slr)*std::pow(rVariables.pb/rVariables.pc,rVariables.lambda)
    //             + rVariables.Slr;

    //     Derivative of the liquid saturation degree with respect to the capillary pressure
    //     rdSldPc = -rVariables.lambda * (1.0 - rVariables.Sgr - rVariables.Slr) * 
    //                 std::pow(rVariables.pb,rVariables.lambda) / std::pow(rVariables.pc,rVariables.lambda+1.0);
    // }

    //TODO. Ignasi
    // Provisional OGS implementation. This is only used in the Liakopoulos test
    // if (rVariables.pc < 0.0) {
    //     rSl = 1.0;
    // } else {
    //     rSl = std::max(rVariables.Slr,1.0 - 1.9722e-11*std::pow(rVariables.pc,2.4279));
    // }
    // double pc_max = std::pow((1.0-rVariables.Slr)/1.9722e-11,(1.0/2.4279));
    // double pc_restr = std::min(rVariables.pc,pc_max);
    // if (rVariables.pc <= 0.0) {
    //     rdSldPc = 0.0;
    // } else {
    //     rdSldPc = -1.9722e-11*2.4279*std::pow(pc_restr,1.4279);
    // }
    // rSl = 1.0 - 0.10152*std::pow(rVariables.pc/(9806.0),2.4279);
    // rdSldPc = (-2.4279*0.10152/std::pow(9806.0,2.4279))*std::pow(rVariables.pc,1.4279);



    // NOTE. This implementation is just done to validate the Khoei example
    // B = 101325
    // v = 5; 
    
    if (rVariables.pc < 0.0) {
        rSl = 1.0;
    } else {
        rSl = std::exp(-rVariables.pc/101325);
    }
    
    if (rVariables.pc < 0.0) {
        rdSldPc = 0.0;
    } else {
        rdSldPc = (-1.0/101325)*std::exp(-rVariables.pc/101325);
    }


}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::CalculateLiquidRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrl = rValues.Getkrl();

    // if (rVariables.Se >= 1.0) {
    //     // Fully saturated medium
    //     rkrl = 1.0;
    // } else if (rVariables.Se <= 0.0) {
    //     // Dry medium
    //     rkrl = rVariables.krmin;
    // } else {
    //     const double nl = (2.0 + 3.0*rVariables.lambda)/rVariables.lambda;
    //     rkrl = std::pow(rVariables.Se,nl);
    // }

    //TODO. Ignasi
    // Provisional OGS implementation. This is only used in the Liakopoulos test
    // double& rSl = rValues.GetSl();
    // if (rSl <= rVariables.Slr) {
    //     rkrl = 0.0;
    // } else if (rSl >= 1.0) {
    //     rkrl = 1.0;
    // } else {
    //     rkrl = 1.0 - 2.207*std::pow(1.0-rSl,1.0121);
    //     rkrl = std::max(rkrl,0.0);
    // }
    // rkrl = 1.0 - 2.207*std::pow(1.0-rSl,1.0121);
    // rkrl = std::max(rkrl,rVariables.krmin);



     // NOTE. This implementation is just done to validate the Khoei example
    // B = 101325
    // v = 5;
    double& rSl = rValues.GetSl();
    if (rSl <= rVariables.Slr) { 
        rkrl = 0.0;
    } else if (rSl >= 1.0) {
        rkrl = 1.0;
    } else {
        rkrl = std::pow(rSl, 5);
    }
    

}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::CalculateGasRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double& rkrg = rValues.Getkrg();

    // if (rVariables.Se >= 1.0) {
    //     // Fully saturated medium
    //     rkrg = rVariables.krmin;
    // } else if (rVariables.Se <= 0.0) {
    //     // Dry medium
    //     rkrg = 1.0;
    // } else {
    //     const double ng = (2.0 + rVariables.lambda)/rVariables.lambda;
    //     rkrg = (1.0-rVariables.Se)*(1.0-rVariables.Se)*(1.0 - std::pow(rVariables.Se,ng));
    //     rkrg = std::max(rkrg,rVariables.krmin);
    // }


 
    // NOTE. This implementation is just done to validate the Khoei example
    // B = 101325
    // v = 5;
    const double& rSl = rValues.GetSl();
    if (rVariables.Se >= 1.0) {
        rkrg = rVariables.krmin;
    } else if (rVariables.Se <= 0.0) {
        rkrg = 1.0;
    } else {
        rkrg = std::pow((1.0 - rSl), 5);
    }

}

}