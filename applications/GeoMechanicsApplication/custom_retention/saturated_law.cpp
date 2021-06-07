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
//  Main authors:    Vahid Galavi
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_retention/saturated_law.h"

namespace Kratos
{
//-------------------------------------------------------------------------------------------------
SaturatedLaw::SaturatedLaw()
    : RetentionLaw()
{
}

//-------------------------------------------------------------------------------------------------
SaturatedLaw::SaturatedLaw(const SaturatedLaw& rOther)
    : RetentionLaw(rOther)
{
}

//-------------------------------------------------------------------------------------------------
RetentionLaw::Pointer SaturatedLaw::Clone() const
{
    return Kratos::make_shared<SaturatedLaw>(*this);
}

//-------------------------------------------------------------------------------------------------
SaturatedLaw::~SaturatedLaw()
{
}

//-------------------------------------------------------------------------------------------------
double SaturatedLaw::
    CalculateSaturation(Parameters &rParameters)
{
    KRATOS_TRY;

    const Properties &rMaterialProperties = rParameters.GetMaterialProperties();

    if (rMaterialProperties.Has(SATURATED_SATURATION))
    {
        return rMaterialProperties[SATURATED_SATURATION];
    }
    else
    {
        return 1.0;
    }

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double SaturatedLaw::
    CalculateEffectiveSaturation(Parameters &rParameters)
{
    KRATOS_TRY;

    return 1.0;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double SaturatedLaw::
    CalculateDerivativeOfSaturation(Parameters &rParameters)
{
    KRATOS_TRY;

    return 0.0;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double SaturatedLaw::
    CalculateRelativePermeability(Parameters &rParameters)
{
    KRATOS_TRY;

    return 1.0;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double SaturatedLaw::
    CalculateBishopCoefficient(Parameters &rParameters)
{
    KRATOS_TRY;

    return CalculateEffectiveSaturation(rParameters);

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double& SaturatedLaw::CalculateValue(RetentionLaw::Parameters& rParameterValues,
                                        const Variable<double>& rThisVariable,
                                        double& rValue)
{
    if (rThisVariable == DEGREE_OF_SATURATION)
    {
        rValue = this->CalculateSaturation(rParameterValues);
        return rValue;
    }
    else if (rThisVariable == EFFECTIVE_SATURATION)
    {
        rValue = this->CalculateEffectiveSaturation(rParameterValues);
        return rValue;
    }
    else if (rThisVariable == BISHOP_COEFICIENT)
    {
        rValue = this->CalculateBishopCoefficient(rParameterValues);
        return rValue;
    }
    else if (rThisVariable == DERIVATIVE_OF_SATURATION)
    {
        rValue = this->CalculateDerivativeOfSaturation(rParameterValues);
        return rValue;
    }
    else if (rThisVariable == RELATIVE_PERMEABILITY)
    {
        rValue = this->CalculateRelativePermeability(rParameterValues);
        return rValue;
    }

    return( rValue );
}

//------------------------- RETENSION LAW GENERAL FEATURES ----------------------------------------
//-------------------------------------------------------------------------------------------------
void SaturatedLaw::
    InitializeMaterial(const Properties& rMaterialProperties,
                       const GeometryType& rElementGeometry,
                       const Vector& rShapeFunctionsValues)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void SaturatedLaw::
    Initialize(Parameters &rParameters)
{
    // nothing is needed
}
//-------------------------------------------------------------------------------------------------
void SaturatedLaw::
    InitializeSolutionStep(Parameters &rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void SaturatedLaw::
    Finalize(Parameters &rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void SaturatedLaw::
    FinalizeSolutionStep(Parameters &rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
int SaturatedLaw::Check(const Properties& rMaterialProperties,
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
