// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#pragma once

#include "includes/constitutive_law.h"

namespace Kratos
{

/**
 * @brief Empty shell for adapting the reference 3D Mohr-Coulomb law to GeoMechanics plane strain.
 *
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) Reference3DMohrCoulombPlaneStrainLaw
    : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(Reference3DMohrCoulombPlaneStrainLaw);

    Reference3DMohrCoulombPlaneStrainLaw() = default;

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    SizeType WorkingSpaceDimension() override;
    [[nodiscard]] SizeType GetStrainSize() const override;
    StrainMeasure GetStrainMeasure() override;
    StressMeasure GetStressMeasure() override;
    void GetLawFeatures(Features& rFeatures) override;

    void InitializeMaterial(const Properties& rMaterialProperties,
                            const Geometry<Node>& rElementGeometry,
                            const Vector& rShapeFunctionsValues) override;
    bool RequiresInitializeMaterialResponse() override;
    bool RequiresFinalizeMaterialResponse() override;
    void InitializeMaterialResponseCauchy(Parameters& rValues) override;
    void CalculateMaterialResponseCauchy(Parameters& rValues) override;
    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;

    bool Has(const Variable<double>& rVariable) override;
    bool Has(const Variable<Vector>& rVariable) override;
    using ConstitutiveLaw::Has;

    double& GetValue(const Variable<double>& rVariable, double& rValue) override;
    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override;
    using ConstitutiveLaw::GetValue;

    std::string Info() const override;

private:
    static constexpr SizeType GeoPlaneStrainSize = 4;
    static constexpr SizeType Reference3DSize     = 6;

    ConstitutiveLaw::Pointer mpReferenceLaw;
    bool   mHasReferenceStrengthParameters = false;
    double mLastReferenceCohesion           = 0.0;
    double mLastReferenceFrictionAngle      = 0.0;

    void EnsureReferenceLawIsCreated();
    void SynchronizeReferenceStrengthParameters(
        const Properties& rReferenceProperties,
        const Geometry<Node>& rElementGeometry,
        const Vector& rShapeFunctionsValues);

    static Vector MapPlaneStrainVectorTo3D(const Vector& rPlaneStrainVector);
    static Vector Map3DVectorToPlaneStrain(const Vector& rThreeDimensionalVector);
    static Matrix Map3DMatrixToPlaneStrain(const Matrix& rThreeDimensionalMatrix);
};

} // namespace Kratos
