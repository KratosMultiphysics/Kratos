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
//  Main authors:    Mohamed Nabi
//                   Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#pragma once

#include "coulomb_yield_surface.h"
#include "tension_cutoff.h"

namespace Kratos
{

class Properties;
class Serializer;

namespace Geo
{
struct PrincipalStresses;
struct SigmaTau;
} // namespace Geo

class CoulombWithTensionCutOffImpl
{
public:
    CoulombWithTensionCutOffImpl() = default;
    explicit CoulombWithTensionCutOffImpl(const Properties& rMaterialProperties);

    [[nodiscard]] bool IsAdmissibleStressState(const Geo::SigmaTau& rTrialSigmaTau) const;
    [[nodiscard]] bool IsAdmissibleStressState(const Geo::PrincipalStresses& rTrialPrincipalStresses) const;

    [[nodiscard]] Vector DoReturnMapping(const Vector&                             rTrialSigmaTau,
                                         CoulombYieldSurface::CoulombAveragingType AveragingType);
    [[nodiscard]] Geo::PrincipalStresses DoReturnMapping(const Geo::PrincipalStresses& rTrialSigmaTau,
                                                         CoulombYieldSurface::CoulombAveragingType AveragingType) const;

    void SaveKappaOfCoulombYieldSurface();
    void RestoreKappaOfCoulombYieldSurface();

private:
    CoulombYieldSurface mCoulombYieldSurface;
    TensionCutoff       mTensionCutOff;
    double              mSavedKappaOfCoulombYieldSurface{0.0};
    double              mAbsoluteYieldFunctionValueTolerance{1.0e-8};
    std::size_t         mMaxNumberOfPlasticIterations{100};

    template <typename StressStateType>
    [[nodiscard]] bool IsAdmissibleStressState(const StressStateType& rTrialStressState) const;

    [[nodiscard]] Vector CalculateCornerPoint() const;
    [[nodiscard]] bool   IsStressAtTensionApexReturnZone(const Vector& rTrialSigmaTau) const;
    [[nodiscard]] bool   IsStressAtTensionCutoffReturnZone(const Vector& rTrialSigmaTau) const;
    [[nodiscard]] bool   IsStressAtCornerReturnZone(const Vector& rTrialSigmaTau,
                                                    CoulombYieldSurface::CoulombAveragingType AveragingType) const;
    [[nodiscard]] Vector ReturnStressAtTensionApexReturnZone() const;
    [[nodiscard]] Vector ReturnStressAtTensionCutoffReturnZone(const Vector& rSigmaTau) const;
    [[nodiscard]] Vector ReturnStressAtRegularFailureZone(const Vector& rSigmaTau,
                                                          CoulombYieldSurface::CoulombAveragingType AveragingType) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const;
    void load(Serializer& rSerializer);
};

} // namespace Kratos
