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

/// Van Genuchten Saturation Law. TODO(U_Pl_pg): change reference according to the new expression used for kr
/// see https://www.sciencedirect.com/science/article/pii/S0266352X22004657

class KRATOS_API(POROMECHANICS_APPLICATION) VanGenuchtenLaw : public SaturationLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(VanGenuchtenLaw);

    ///------------------------------------------------------------------------------------------------

    VanGenuchtenLaw() = default;

    VanGenuchtenLaw (const VanGenuchtenLaw& rOther) : SaturationLaw(rOther)
    {
    }

    ~VanGenuchtenLaw() override
    {
    }

    SaturationLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<VanGenuchtenLaw>(VanGenuchtenLaw(*this));
    }

    ///------------------------------------------------------------------------------------------------

protected:

    void CalculateLiquidSaturationDegree(SaturationLawVariables& rVariables, Parameters& rValues) override;

    void CalculateLiquidRelativePermeability(SaturationLawVariables& rVariables, Parameters& rValues) override;

    void CalculateGasRelativePermeability(SaturationLawVariables& rVariables, Parameters& rValues) override;

    ///------------------------------------------------------------------------------------------------

private:

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SaturationLaw )
    }

    void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SaturationLaw )
    }

}; /* Class VanGenuchtenLaw */

} /* namespace Kratos.*/
