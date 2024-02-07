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

int BrooksAndCoreyLaw::Check(const Properties& rMaterialProperties,
                           const GeometryType& rElementGeometry,
                           const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int ierr = SaturationLaw::Check(rCurrentProcessInfo);
    if(ierr != 0)
        return ierr;

    return ierr;

    KRATOS_CATCH("");
}

//------------------------------------------------------------------------------------------------

// void BrooksAndCoreyLaw::InitializeMaterial(const Properties& rMaterialProperties,
//         const GeometryType& rElementGeometry,
//         const Vector& rShapeFunctionsValues)
// {
//     SaturationLaw::InitializeMaterial(rMaterialProperties,rElementGeometry,rShapeFunctionsValues);
// }

//------------------------------------------------------------------------------------------------

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
        rSw = (1.0 - rVariables.Swr)*pow(rVariables.pb/rVariables.pc,rVariables.lambda) + rVariables.Swr;

        // Derivative of the water saturation degree with respect to the capilar pressure
        rdSwdPc = (1.0 - rVariables.Swr) * 
                    rVariables.lambda * rVariables.pb * pow(rVariables.pb/rVariables.pc,rVariables.lambda-1.0) /
                    (rVariables.pc * rVariables.pc);
    }
}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::WaterRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double krw, nw, lambda;

    // Compute the water relative permeability according with the consitutive model
        // -- Brooks and Corey
        //           (see pg. 479 from Khoei's 2015 book: Extended Finite Element: theory and applications, ISBN 978-1-118-45768-9) ----
           
            lambda = Variables.PoreSizeFactor;
            nw     = (2.0 + 3.0*lambda)/lambda;
            krw    = pow(Se,nw);

    krw = std::min(krw,0.0001);
    return krw;    
}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::GasRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double krg, ng, lambda;

    // Compute the gas relative permeability according with the consitutive model
        // -- Brooks and Corey
        //           (see pg. 479 from Khoei's 2015 book: Extended Finite Element: theory and applications, ISBN 978-1-118-45768-9) -----
            
            lambda = Variables.PoreSizeFactor;
            ng     = (2.0 + lambda)/lambda;
            krg    = pow(1.0-Se,2.0)*(1.0 - pow(Se,ng));

    krg = std::min(krg,0.0001);
    return krg;    
}

}