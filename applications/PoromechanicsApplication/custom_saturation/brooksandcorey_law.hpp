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

class KRATOS_API(POROMECHANICS_APPLICATION) BrooksAndCoreyLaw : public SaturationLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(BrooksAndCoreyLaw);

    ///------------------------------------------------------------------------------------------------

    BrooksAndCoreyLaw() = default;

    BrooksAndCoreyLaw (const BrooksAndCoreyLaw& rOther) : SaturationLaw(rOther)
    {
    }

    ~BrooksAndCoreyLaw() override
    {
    }

    SaturationLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<BrooksAndCoreyLaw>(BrooksAndCoreyLaw(*this));
    }

    ///------------------------------------------------------------------------------------------------

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;

    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues ) override;

    ///------------------------------------------------------------------------------------------------

    void CalculateMaterialResponse (Parameters & rValues) override;

    ///------------------------------------------------------------------------------------------------

protected:

    virtual void InitializeSaturationLawVariables(SaturationLawVariables& rVariables, Parameters& rValues);

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

}; /* Class BrooksAndCoreyLaw */

} /* namespace Kratos.*/
