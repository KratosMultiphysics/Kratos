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
//  Main authors:    Wijtze Pieter Kikstra
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "includes/table.h"

namespace Kratos
{

/**
 * @namespace TrussBackboneConstitutiveLaw
 * @brief This constitutive law represents a 1D backbone law with linear elastic un- and reloading
 * @author Wijtze Pieter Kikstra
 */

class KRATOS_API(GEO_MECHANICS_APPLICATION) TrussBackboneConstitutiveLaw : public ConstitutiveLaw
{
public:
    using BaseType = ConstitutiveLaw;

    KRATOS_CLASS_POINTER_DEFINITION(TrussBackboneConstitutiveLaw);

    void InitializeMaterial(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector&       rShapeFunctionsValues) override;

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    void GetLawFeatures(Features& rFeatures) override;

    void SetValue(const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    using BaseType::SetValue;

    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;
    using BaseType::GetValue;

    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<double>&      rThisVariable,
                           double&                      rValue) override;
    using BaseType::CalculateValue;

    void FinalizeMaterialResponsePK2(Parameters& rValues) override;

    void CalculateMaterialResponsePK2(Parameters& rValues) override;

    bool RequiresInitializeMaterialResponse() override;

    [[nodiscard]] SizeType GetStrainSize() const override;

    [[nodiscard]] int Check(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const ProcessInfo&  rCurrentProcessInfo) const override;

    std::string Info() const override;

private:
    double                mAccumulatedStrain   = 0.;
    double                mPreviousAxialStrain = 0.;
    double                mUnReLoadCenter      = 0.;
    Table<double, double> mStressStrainTable;

    [[nodiscard]] double BackboneStress(double Strain) const;
    [[nodiscard]] double BackboneStiffness(double Strain) const;
    [[nodiscard]] double CalculateUnReLoadAmplitude(double YoungsModulus) const;
    [[nodiscard]] bool   IsWithinUnReLoading(double Strain, double YoungsModulus) const;

    static void CheckStressStrainDiagram(const Properties& rMaterialProperties);
    static void CheckStrainValuesAreAscending(const Vector& rStrains);
    static void CheckBackboneStiffnessesDontExceedYoungsModulus(const Vector& rStrains,
                                                                const Vector& rStresses,
                                                                double        YoungsModulus);

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class TrussBackboneConstitutiveLaw
} // namespace Kratos.
