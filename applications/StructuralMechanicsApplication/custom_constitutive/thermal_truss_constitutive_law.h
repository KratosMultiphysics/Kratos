// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Talhah Ansari
//

#pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/truss_constitutive_law.h"

namespace Kratos
{

/**
 * @namespace ThermalTrussConstitutiveLaw
 *
 * @brief This constitutive law represents a thermo dependent linear elastic 1D law, including the addition of thermal expansion strains
 *
 * @author Talhah Ansari
 */


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ThermalTrussConstitutiveLaw : public TrussConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( ThermalTrussConstitutiveLaw );

    ThermalTrussConstitutiveLaw();
    ConstitutiveLaw::Pointer Clone() const override;
    ThermalTrussConstitutiveLaw(const ThermalTrussConstitutiveLaw& rOther);
    ~ThermalTrussConstitutiveLaw() override;

    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
    ) override;

    double CalculateStressElastic(ConstitutiveLaw::Parameters& rParameterValues) override;

private:
    double mReferenceTemperature = 0.0;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw);
        rSerializer.save("ReferenceTemperature", mReferenceTemperature);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw);
        rSerializer.load("ReferenceTemperature", mReferenceTemperature);
    }
};

} // namespace Kratos
