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
#include "custom_retention/van_genuchten_law.h"

namespace Kratos
{
//-------------------------------------------------------------------------------------------------
VanGenuchtenLaw::VanGenuchtenLaw()
    : RetentionLaw()
{
}

//-------------------------------------------------------------------------------------------------
VanGenuchtenLaw::VanGenuchtenLaw(const VanGenuchtenLaw& rOther)
    : RetentionLaw(rOther)
{
}

//-------------------------------------------------------------------------------------------------
RetentionLaw::Pointer VanGenuchtenLaw::Clone() const
{
    return Kratos::make_shared<VanGenuchtenLaw>(*this);
}

//-------------------------------------------------------------------------------------------------
VanGenuchtenLaw::~VanGenuchtenLaw()
{
}

//-------------------------------------------------------------------------------------------------
double VanGenuchtenLaw::
    CalculateSaturation(Parameters &rParameters)
{
    KRATOS_TRY;

    const double &p = rParameters.GetFluidPressure();
    const Properties &rMaterialProperties = rParameters.GetMaterialProperties();

    if (p > 0.0)
    {
        const double &satMax = rMaterialProperties[SATURATED_SATURATION];
        const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
        const double &pb     = rMaterialProperties[VAN_GENUCHTEN_AIR_ENTRY_PRESSURE];
        const double &gn     = rMaterialProperties[VAN_GENUCHTEN_GN];
        const double gc      = (1.0 - gn) / gn;

        double sat = satMin + (satMax - satMin) * pow(1.0 + pow(p/pb, gn), gc);
        return sat;
    }
    else
    {
        return rMaterialProperties[SATURATED_SATURATION];
    }

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double VanGenuchtenLaw::
    CalculateEffectiveSaturation(Parameters &rParameters)
{
    KRATOS_TRY;

    const double sat = CalculateSaturation(rParameters);

    const auto &rMaterialProperties = rParameters.GetMaterialProperties();
    const double &satMax = rMaterialProperties[SATURATED_SATURATION];
    const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];

    double effectiveSat = (sat - satMin) / (satMax - satMin);

    return effectiveSat;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double VanGenuchtenLaw::
    CalculateDerivativeOfSaturation(Parameters &rParameters)
{
    KRATOS_TRY;
    const double &p = rParameters.GetFluidPressure();

    if (p > 0.0)
    {
        const auto &rMaterialProperties = rParameters.GetMaterialProperties();
        const double &satMax = rMaterialProperties[SATURATED_SATURATION];
        const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
        const double &pb     = rMaterialProperties[VAN_GENUCHTEN_AIR_ENTRY_PRESSURE];
        const double &gn     = rMaterialProperties[VAN_GENUCHTEN_GN];
        const double gc      = (1.0 - gn) / gn;

        double dSdp = (satMax - satMin) * gc * pow((1.0 + pow(p/pb, gn)), gc-1.0)
                                        * gn * pow(pb,-gn) * pow(p, gn-1.0);
        return dSdp;
    }
    else
    {
        return 0.0;
    }

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double VanGenuchtenLaw::
    CalculateRelativePermeability(Parameters &rParameters)
{
    KRATOS_TRY;

    const double effSat = CalculateEffectiveSaturation(rParameters);

    const auto &rMaterialProperties = rParameters.GetMaterialProperties();
    const double &gl = rMaterialProperties[VAN_GENUCHTEN_GL];
    const double &gn = rMaterialProperties[VAN_GENUCHTEN_GN];

    double relPerm = pow(effSat, gl) * pow(1.0 - pow(1.0 - pow(effSat, gn/(gn-1.0)), (gn-1.0)/gn), 2);

    const double &minRelPerm = rMaterialProperties[MINIMUM_RELATIVE_PERMEABILITY];

    relPerm = std::max(relPerm, minRelPerm);

    return relPerm;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double VanGenuchtenLaw::
    CalculateBishopCoefficient(Parameters &rParameters)
{
    KRATOS_TRY;

    const double effSat = CalculateEffectiveSaturation(rParameters);
    return effSat;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double& VanGenuchtenLaw::CalculateValue(RetentionLaw::Parameters& rParameterValues,
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
void VanGenuchtenLaw::
    InitializeMaterial(const Properties& rMaterialProperties,
                       const GeometryType& rElementGeometry,
                       const Vector& rShapeFunctionsValues)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void VanGenuchtenLaw::
    Initialize(Parameters &rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void VanGenuchtenLaw::
    InitializeSolutionStep(Parameters &rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void VanGenuchtenLaw::
    Finalize(Parameters &rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void VanGenuchtenLaw::
    FinalizeSolutionStep(Parameters &rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
int VanGenuchtenLaw::Check(const Properties& rMaterialProperties,
                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] < 0.0)
                    << "SATURATED_SATURATION cannot be less than 0 " << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[RESIDUAL_SATURATION] < 0.0)
                    << "RESIDUAL_SATURATION cannot be less than 0 " << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] > 1.0)
                    << "SATURATED_SATURATION cannot be greater than 1.0 " << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[RESIDUAL_SATURATION] > 1.0)
                    << "RESIDUAL_SATURATION cannot be greater than 1.0 " << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] < rMaterialProperties[RESIDUAL_SATURATION])
                    << "RESIDUAL_SATURATION cannot be greater than SATURATED_SATURATION " << std::endl;

    KRATOS_ERROR_IF(!(rMaterialProperties[VAN_GENUCHTEN_AIR_ENTRY_PRESSURE] > 0.0))
                    << "VAN_GENUCHTEN_AIR_ENTRY_PRESSURE must be greater than 0 " << std::endl;

    KRATOS_ERROR_IF(!(rMaterialProperties[MINIMUM_RELATIVE_PERMEABILITY] > 0.0))
                    << "MINIMUM_RELATIVE_PERMEABILITY must be greater than 0 " << std::endl;

    return 0;
}

} // Namespace Kratos
