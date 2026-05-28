// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: IgaApplication/license.txt
//

#pragma once

#include "includes/constitutive_law.h"

#include "custom_constitutive/structural_elements_constitutive_laws/thickness_integrated_composite_constitutive_law.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) IgaThicknessIntegratedCompositeLaw
    : public ThicknessIntegratedCompositeConstitutiveLaw
{
public:

    using BaseType = ThicknessIntegratedCompositeConstitutiveLaw;
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(IgaThicknessIntegratedCompositeLaw);

    IgaThicknessIntegratedCompositeLaw();

    IgaThicknessIntegratedCompositeLaw(
        const std::vector<double>& rZCoordinates,
        const std::vector<double>& rEulerAngles,
        const std::vector<double>& rThicknesses);

    IgaThicknessIntegratedCompositeLaw(const IgaThicknessIntegratedCompositeLaw& rOther);

    ~IgaThicknessIntegratedCompositeLaw() override;

    ConstitutiveLaw::Pointer Clone() const override;

    ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override;

    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
    void CalculateMaterialResponsePK2   (ConstitutiveLaw::Parameters& rValues) override;

    void InitializeShearReductionFactors(const Properties& rMaterialProperties) override;

    SizeType GetNumberOfLayers() const { return mConstitutiveLaws.size(); }

    void CalculateLayerGeneralizedStiffness(
        IndexType LayerIndex,
        ConstitutiveLaw::Parameters& rValues,
        Matrix& rQbar,
        Matrix& rShear,
        double& rZLo,
        double& rZHi);

    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo) const override;

private:

    enum class PlyStressMeasure { Cauchy, PK2 };
    void CalculateLayeredResponse(
        ConstitutiveLaw::Parameters& rValues,
        PlyStressMeasure PlyMeasure);

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    }

};

}
