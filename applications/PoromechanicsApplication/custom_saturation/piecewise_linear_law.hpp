//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Danilo Cavalcanti
//                   Ignasi de Pouplana
//

#pragma once

/* System includes */
#include <cmath>

/* Project includes */
#include "includes/serializer.h"
#include "includes/checks.h"

// Application includes
#include "custom_saturation/saturation_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos {

/// PiecewiseLinearLaw / UMAT Saturation Law
/// Piecewise linear retention model for Buckley-Leverett and user-defined capillary curves
/// Ref: PorousLab v1.1.1 (Cavalcanti et al.)

class KRATOS_API(POROMECHANICS_APPLICATION) PiecewiseLinearLaw : public SaturationLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(PiecewiseLinearLaw);

    PiecewiseLinearLaw() = default;

    PiecewiseLinearLaw (const PiecewiseLinearLaw& rOther) : SaturationLaw(rOther)
    {
    }

    ~PiecewiseLinearLaw() override
    {
    }

    SaturationLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<PiecewiseLinearLaw>(PiecewiseLinearLaw(*this));
    }

protected:

    void CalculateLiquidSaturationDegree(SaturationLawVariables& rVariables, Parameters& rValues) override;

    void CalculateLiquidRelativePermeability(SaturationLawVariables& rVariables, Parameters& rValues) override;

    void CalculateGasRelativePermeability(SaturationLawVariables& rVariables, Parameters& rValues) override;

private:

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SaturationLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SaturationLaw )
    }

}; /* Class PiecewiseLinearLaw */

} /* namespace Kratos.*/

