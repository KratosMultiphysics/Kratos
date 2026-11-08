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

/* External includes */

/* Project includes */
#include "includes/serializer.h"
#include "includes/checks.h"

// Application includes
#include "custom_saturation/saturation_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos {

/// Liakopoulos Saturation Law
/// Empirical retention and relative permeability models for Liakopoulos drainage benchmark
/// Ref: Asadi, R., Ataie-Ashtiani, B. (2015), PorousLab v1.1.0

class KRATOS_API(POROMECHANICS_APPLICATION) LiakopoulosLaw : public SaturationLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(LiakopoulosLaw);

    LiakopoulosLaw() = default;

    LiakopoulosLaw (const LiakopoulosLaw& rOther) : SaturationLaw(rOther)
    {
    }

    ~LiakopoulosLaw() override
    {
    }

    SaturationLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<LiakopoulosLaw>(LiakopoulosLaw(*this));
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

}; /* Class LiakopoulosLaw */

} /* namespace Kratos.*/
