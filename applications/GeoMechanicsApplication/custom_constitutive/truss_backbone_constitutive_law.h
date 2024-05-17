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
#include "includes/checks.h"

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
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw  BaseType;
    typedef std::size_t      SizeType;

    KRATOS_CLASS_POINTER_DEFINITION( TrussBackboneConstitutiveLaw );

    TrussBackboneConstitutiveLaw();

    ConstitutiveLaw::Pointer Clone() const override;

    TrussBackboneConstitutiveLaw(const TrussBackboneConstitutiveLaw& rOther);

    ~TrussBackboneConstitutiveLaw();
    void GetLawFeatures(Features& rFeatures) override;

    void SetValue(const Variable<double>& rThisVariable,
                  const double&           rValue,
                  const ProcessInfo&      rCurrentProcessInfo) override;

    double& GetValue(const Variable<double>& rThisVariable,
                     double&                 rValue) override;

    array_1d<double, 3>& GetValue(const Variable<array_1d<double, 3>>& rThisVariable,
                                  array_1d<double, 3>&                 rValue) override;

    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<double>&      rThisVariable,
                           double&                      rValue) override;

    void FinalizeMaterialResponsePK2(Parameters& rValues) override;

    void CalculateMaterialResponsePK2(Parameters& rValues) override;

    bool RequiresInitializeMaterialResponse() override
    {
        return false;
    }

    SizeType GetStrainSize() const override
    {
        return 1;
    }

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rMaterialProperties: The properties of the material
     * @param rElementGeometry: The geometry of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    int Check(const Properties&   rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo&  rCurrentProcessInfo) const override;

private:
    double mAccumulatedStrain = 0.;

    double BackboneStress(const double Strain) const;
    double BackboneStiffness(const double Strain) const;

    friend class Serializer;

    void save(Serializer& rSerializer) const override {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw);
        rSerializer.save("AccumulatedStrain", mAccumulatedStrain);
    }

    void load(Serializer& rSerializer) override {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw);
        rSerializer.load("AccumulatedStrain", mAccumulatedStrain);
    }

}; // Class TrussBackboneConstitutiveLaw
}  // namespace Kratos.
