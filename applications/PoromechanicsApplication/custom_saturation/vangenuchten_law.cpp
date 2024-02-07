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

int VanGenuchtenLaw::Check(const Properties& rMaterialProperties,
                           const GeometryType& rElementGeometry,
                           const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int ierr = BrooksAndCoreyLaw::Check(rCurrentProcessInfo);
    if(ierr != 0)
        return ierr;

    return ierr;

    KRATOS_CATCH("");
}

//------------------------------------------------------------------------------------------------

// void VanGenuchtenLaw::InitializeMaterial(const Properties& rMaterialProperties,
//         const GeometryType& rElementGeometry,
//         const Vector& rShapeFunctionsValues)
// {
//     BrooksAndCoreyLaw::InitializeMaterial(rMaterialProperties,rElementGeometry,rShapeFunctionsValues);
// }

//------------------------------------------------------------------------------------------------

void VanGenuchtenLaw::CalculateMaterialResponse (Parameters& rValues)
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

void VanGenuchtenLaw::CalculateSaturation (Parameters& rValues)
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

void VanGenuchtenLaw::InitializeSaturationLawVariables (SaturationLawVariables& rVariables, Parameters& rValues)
{
    
}

//------------------------------------------------------------------------------------------------

void VanGenuchtenLaw::CalculateWaterSaturationDegree (SaturationLawVariables& rVariables, Parameters& rValues)
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
            // -- van Genuchten (https://www.sciencedirect.com/science/article/pii/S0266352X22004657) -----

            // Water saturation degree
            rVariables.Sw = (1.0 - Swr)*pow(1.0 + pow(pc/pb,1.0/(1.0-lambda)),-lambda) + Swr;

            // Derivative of the water saturation degree with respect to the capilar pressure
            rVariables.dSwdPc = (1.0 - Swr)*lambda/(pb * pow(pc/pb,1.0+1.0/(lambda-1.0)) * (lambda - 1.0) * pow(1.0+1.0/(pow(pc/pb,1/(lambda-1.0))),lambda+1.0));
                
    }
}

//------------------------------------------------------------------------------------------------

void VanGenuchtenLaw::WaterRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double krw, nw, lambda;

    // Compute the water relative permeability according with the consitutive model
        // -- van Genuchten (https://www.sciencedirect.com/science/article/pii/S0266352X22004657) -----

            nw  = 1.5;
            krw = pow(Se,nw);

    krw = std::min(krw,0.0001);
    return krw;    
}

//------------------------------------------------------------------------------------------------

void VanGenuchtenLaw::GasRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    double krg, ng, lambda;

    // Compute the gas relative permeability according with the consitutive model
        // -- van Genuchten (https://www.sciencedirect.com/science/article/pii/S0266352X22004657) -----

            ng  = 3.0;
            krg = pow(1.0-Se,ng);

    krg = std::min(krg,0.0001);
    return krg;    
}
}