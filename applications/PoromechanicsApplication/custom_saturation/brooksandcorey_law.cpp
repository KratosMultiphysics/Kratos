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

void BrooksAndCoreyLaw::CalculateMaterialResponse (Parameters& rValues)
{
    //Check
    rValues.CheckAllParameters();

    //Initialize main variables
    SaturationLawVariables Variables;
    this->InitializeConstitutiveLawVariables(Variables,rValues);

        // //Compute the capilar pressure at the integration point
        // Variables.ipCapilarPressure = inner_prod(Variables.Np,Variables.CapilarPressureVector);

    this->CalculateWaterSaturationDegree(Variables,rValues);

    this->EffectiveSaturation(Variables,rValues);
    this->WaterRelativePermeability(Variables,rValues);
    this->GasRelativePermeability(Variables,rValues);
    



}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::CalculateSaturation (Parameters& rValues)
{
    //Check
    rValues.CheckAllParameters();

    //Initialize main variables
    SaturationLawVariables Variables;
    this->InitializeConstitutiveLawVariables(Variables,rValues);

        // //Compute the capilar pressure at the integration point
        // Variables.ipCapilarPressure = inner_prod(Variables.Np,Variables.CapilarPressureVector);

    this->CalculateWaterSaturationDegree(Variables,rValues);
    



}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::InitializeSaturationLawVariables (SaturationLawVariables& rVariables, Parameters& rValues)
{
    
}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::CalculateWaterSaturationDegree (SaturationLawVariables& rVariables, Parameters& rValues)
{
    //Get material parameters
    double Swr    = rVariables.ResidualWaterSaturation;
    double lambda = rVariables.PoreSizeFactor;
    double pb     = rVariables.GasEntryPressure;
    double pc     = rVariables.ipCapilarPressure;

    // If the capillar pressure is lower than the gas-entry pressure, the porous media is fully saturated with the wetting phase.
    rVariables.Sw     = 1.0;
    rVariables.dSwdPc = 0.0;

    if(pc > pb)
    {
            // -- Brooks and Corey
            //           (see pg. 479 from Khoei's 2015 book: Extended Finite Element: theory and applications, ISBN 978-1-118-45768-9) -----
                
                // Water saturation degree
                rVariables.Sw = (1.0 - Swr)*pow(pb/pc,lambda) + Swr;

                // Derivative of the water saturation degree with respect to the capilar pressure
                rVariables.dSwdPc = (1.0 - Swr) * lambda * pb * pow(pb/pc,lambda-1) / (pc * pc);
    }
}

//------------------------------------------------------------------------------------------------

void BrooksAndCoreyLaw::EffectiveSaturation (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double Se = (Sw - Swr)/(1.0 - Swr);
    return Se;    
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