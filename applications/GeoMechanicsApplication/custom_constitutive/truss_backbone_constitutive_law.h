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
#include "includes/checks.h"
#include "includes/constitutive_law.h"

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

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    void GetLawFeatures(Features& rFeatures) override;

    void SetValue(const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    using ConstitutiveLaw::SetValue;

    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;
    using ConstitutiveLaw::GetValue;

    array_1d<double, 3>& GetValue(const Variable<array_1d<double, 3>>& rThisVariable,
                                  array_1d<double, 3>&                 rValue) override;

    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<double>&      rThisVariable,
                           double&                      rValue) override;
    using ConstitutiveLaw::CalculateValue;

    void FinalizeMaterialResponsePK2(Parameters& rValues) override;

    void CalculateMaterialResponsePK2(Parameters& rValues) override;

    bool RequiresInitializeMaterialResponse() override;

    [[nodiscard]] SizeType GetStrainSize() const override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rMaterialProperties: The properties of the material
     * @param rElementGeometry: The geometry of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    [[nodiscard]] int Check(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const ProcessInfo&  rCurrentProcessInfo) const override;

private:
    double mAccumulatedStrain   = 0.;
    double mPreviousAxialStrain = 0.;
    double mUnReLoadCenter      = 0.;

    [[nodiscard]] double BackboneStress(double Strain) const;
    [[nodiscard]] double BackboneStiffness(double Strain) const;
    [[nodiscard]] double CalculateUnReLoadAmplitude(double YoungsModulus) const;
    [[nodiscard]] bool IsWithinUnReLoading(double Strain, double YoungsModulus) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class TrussBackboneConstitutiveLaw
} // namespace Kratos.
