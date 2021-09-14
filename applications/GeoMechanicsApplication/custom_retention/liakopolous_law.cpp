// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Hoang-Giang Bui
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_retention/liakopolous_law.h"

namespace Kratos
{
//-------------------------------------------------------------------------------------------------
LiakopolousLaw::LiakopolousLaw()
    : RetentionLaw()
{
}

//-------------------------------------------------------------------------------------------------
LiakopolousLaw::LiakopolousLaw(const LiakopolousLaw& rOther)
    : RetentionLaw(rOther)
{
}

//-------------------------------------------------------------------------------------------------
RetentionLaw::Pointer LiakopolousLaw::Clone() const
{
    return Kratos::make_shared<LiakopolousLaw>(*this);
}

//-------------------------------------------------------------------------------------------------
LiakopolousLaw::~LiakopolousLaw()
{
}

//-------------------------------------------------------------------------------------------------
double& LiakopolousLaw::
    CalculateValue(RetentionLaw::Parameters& rParameters,
                   const Variable<double>& rThisVariable,
                   double& rValue)
{
    if (rThisVariable == DENSITY_AIR)
    {
        const double initial_density = rParameters.GetMaterialProperties()[DENSITY_AIR];
        const double bulk_air = rParameters.GetMaterialProperties()[BULK_AIR];
        rValue = initial_density + bulk_air*rParameters.GetAirPressure();
        return rValue;
    }
    else if (rThisVariable == SATURATION)
    {
        double capillaryPressure;
        capillaryPressure = rParameters.GetCapillaryPressure(capillaryPressure);
        if (capillaryPressure <= 0.0) // no matric suction
        {
            rValue = 1.0;
        }
        else
        {
            rValue = 1.0-1.9722*1e-11*pow(capillaryPressure, 2.4279);
        }
        return rValue;
    }
    else if (rThisVariable == PERMEABILITY_WATER)
    {
        double saturation;
        this->CalculateValue(rParameters, SATURATION, saturation);
        rValue = 1.0 - 2.207*pow(1.0-saturation, 1.0121);
        return rValue;
    }
    else if (rThisVariable == PERMEABILITY_AIR)
    {
        double saturation;
        this->CalculateValue(rParameters, SATURATION, saturation);
        double relSat = (saturation-0.2)/0.8;
        rValue = 0.0001+pow((1.0-relSat),2)*(1-pow(relSat,5.0/3.0));
        return rValue;
    }

    return( rValue );
}

//-------------------------------------------------------------------------------------------------
double& LiakopolousLaw::
    CalculateDerivative(RetentionLaw::Parameters& rParameters,
                        const Variable<double>& rThisVariable,
                        double &rValue)
{
    if (rThisVariable == DENSITY_AIR)
    {
        const double bulk_air = rParameters.GetMaterialProperties()[BULK_AIR];
        rValue = bulk_air;
        return rValue;
    }
    else if (rThisVariable == SATURATION)
    {
        double capillaryPressure;
        capillaryPressure = rParameters.GetCapillaryPressure(capillaryPressure);
        if (capillaryPressure <= 0.0) // no matric suction
        {
            rValue = 0.0;
        }
        else
        {
            rValue = -1.9722*2.4279*1e-11*pow(capillaryPressure, 1.4279);
        }
        return rValue;
    }
    else if (rThisVariable == PERMEABILITY_WATER)
    {
        double saturation;
        this->CalculateValue(rParameters, SATURATION, saturation);
        rValue = 2.207*1.0121*pow(1.0-saturation, 0.0121);
        return rValue;
    }
    else if (rThisVariable == PERMEABILITY_AIR)
    {
        double saturation;
        this->CalculateValue(rParameters, SATURATION, saturation);
        double relSat = (saturation-0.2)/0.8;
        rValue = -2*(1.0-relSat)*(-1)/0.8*(1-pow(relSat,5.0/3.0))
            + pow((1.0-relSat),2)*(-1)*5.0/3.0*pow(relSat,2.0/3.0)/0.8;
        return rValue;
    }

    return( rValue );
}

//-------------------------------------------------------------------------------------------------
double& LiakopolousLaw::
    CalculateSecondDerivative(RetentionLaw::Parameters& rParameters,
                              const Variable<double>& rThisVariable,
                              double &rValue)
{
    if (rThisVariable == SATURATION)
    {
        double capillaryPressure;
        capillaryPressure = rParameters.GetCapillaryPressure(capillaryPressure);
        if (capillaryPressure <= 0.0) // no matric suction
        {
            rValue = 0.0;
        }
        else
        {
            rValue = -1.9722*2.4279*1.4279*1e-11*pow(capillaryPressure, 0.4279);
        }
        return rValue;
    }

    return( rValue );
}

//------------------------- RETENSION LAW GENERAL FEATURES ----------------------------------------
//-------------------------------------------------------------------------------------------------
void LiakopolousLaw::
    InitializeMaterial(const Properties& rMaterialProperties,
                       const GeometryType& rElementGeometry,
                       const Vector& rShapeFunctionsValues)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void LiakopolousLaw::
    Initialize(RetentionLaw::Parameters& rParameters)
{
    // nothing is needed
}
//-------------------------------------------------------------------------------------------------
void LiakopolousLaw::
    InitializeSolutionStep(RetentionLaw::Parameters& rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void LiakopolousLaw::
    Finalize(RetentionLaw::Parameters& rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void LiakopolousLaw::
    FinalizeSolutionStep(RetentionLaw::Parameters& rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
int LiakopolousLaw::Check(const Properties& rMaterialProperties,
                           const ProcessInfo& rCurrentProcessInfo)
{
    if (rMaterialProperties.Has(SATURATED_SATURATION))
    {
        KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] < 0.0)
                        << "SATURATED_SATURATION cannot be less than 0 " << std::endl;

        KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] > 1.0)
                        << "SATURATED_SATURATION cannot be greater than 1.0 " << std::endl;
    }

    return 0;
}

} // Namespace Kratos
