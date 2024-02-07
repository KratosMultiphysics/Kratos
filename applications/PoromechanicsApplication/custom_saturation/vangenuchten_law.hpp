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
#include "custom_saturation/brooksandcorey_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos {

class KRATOS_API(POROMECHANICS_APPLICATION) VanGenuchtenLaw : public BrooksAndCoreyLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(VanGenuchtenLaw);

    ///------------------------------------------------------------------------------------------------

    VanGenuchtenLaw() = default;

    VanGenuchtenLaw (const VanGenuchtenLaw& rOther) : BrooksAndCoreyLaw(rOther)
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

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;

    // void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues ) override;

    ///------------------------------------------------------------------------------------------------

    void CalculateMaterialResponse (Parameters & rValues) override;

    void CalculateSaturation (Parameters & rValues) override;

    ///------------------------------------------------------------------------------------------------

protected:

    void InitializeSaturationLawVariables(SaturationLawVariables& rVariables, Parameters& rValues) override;

    void CalculateWaterSaturationDegree(SaturationLawVariables& rVariables, Parameters& rValues) override;

    void WaterRelativePermeability(SaturationLawVariables& rVariables, Parameters& rValues) override;

    void GasRelativePermeability(SaturationLawVariables& rVariables, Parameters& rValues) override;

    ///------------------------------------------------------------------------------------------------

private:

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BrooksAndCoreyLaw )
    }

    void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BrooksAndCoreyLaw )
    }

}; /* Class VanGenuchtenLaw */

} /* namespace Kratos.*/
