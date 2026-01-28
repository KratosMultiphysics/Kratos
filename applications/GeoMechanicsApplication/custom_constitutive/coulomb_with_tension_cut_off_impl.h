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

    [[nodiscard]] Geo::SigmaTau DoReturnMapping(const Geo::SigmaTau& rTrialSigmaTau,
                                                CoulombYieldSurface::CoulombAveragingType AveragingType);
    [[nodiscard]] Geo::PrincipalStresses DoReturnMapping(const Geo::PrincipalStresses& rTrialSigmaTau,
                                                         CoulombYieldSurface::CoulombAveragingType AveragingType);

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
    template <typename StressStateType>
    [[nodiscard]] StressStateType DoReturnMapping(const StressStateType& rTrialStressState,
                                                  CoulombYieldSurface::CoulombAveragingType AveragingType);

    [[nodiscard]] Geo::SigmaTau CalculateCornerPoint(const Geo::SigmaTau&) const;
    [[nodiscard]] Geo::PrincipalStresses CalculateCornerPoint(const Geo::PrincipalStresses& rPrincipalStresses) const;
    [[nodiscard]] bool IsStressAtTensionApexReturnZone(const Vector& rTrialSigmaTau) const;
    [[nodiscard]] bool IsStressAtTensionApexReturnZone(const Geo::SigmaTau& rTrialSigmaTau) const;
    [[nodiscard]] bool IsStressAtTensionApexReturnZone(const Geo::PrincipalStresses& rTrialPrincipalStresses) const;
    [[nodiscard]] bool IsStressAtTensionCutoffReturnZone(const Vector& rTrialSigmaTau) const;
    [[nodiscard]] bool IsStressAtTensionCutoffReturnZone(const Geo::SigmaTau& rTrialSigmaTau) const;
    [[nodiscard]] bool IsStressAtTensionCutoffReturnZone(const Geo::PrincipalStresses& rTrialPrincipalStresses) const;
    [[nodiscard]] bool IsStressAtCornerReturnZone(const Vector& rTrialSigmaTau,
                                                  CoulombYieldSurface::CoulombAveragingType AveragingType) const;
    [[nodiscard]] bool IsStressAtCornerReturnZone(const Geo::SigmaTau& rTrialSigmaTau,
                                                  CoulombYieldSurface::CoulombAveragingType AveragingType) const;
    [[nodiscard]] bool IsStressAtCornerReturnZone(const Geo::PrincipalStresses& rTrialPrincipalStresses,
                                                  CoulombYieldSurface::CoulombAveragingType AveragingType) const;

    [[nodiscard]] Vector        ReturnStressAtTensionApexReturnZone() const;
    [[nodiscard]] Geo::SigmaTau ReturnStressAtTensionApexReturnZone(const Geo::SigmaTau&) const;
    [[nodiscard]] Geo::PrincipalStresses ReturnStressAtTensionApexReturnZone(const Geo::PrincipalStresses& rPrincipalStresses) const;
    [[nodiscard]] Vector ReturnStressAtTensionCutoffReturnZone(const Vector& rSigmaTau) const;
    [[nodiscard]] Geo::SigmaTau ReturnStressAtTensionCutoffReturnZone(const Geo::SigmaTau& rSigmaTau) const;
    [[nodiscard]] Geo::PrincipalStresses ReturnStressAtTensionCutoffReturnZone(const Geo::PrincipalStresses& rPrincipalStresses) const;
    [[nodiscard]] Vector        ReturnStressAtRegularFailureZone(const Vector& rSigmaTau,
                                                                 CoulombYieldSurface::CoulombAveragingType AveragingType) const;
    [[nodiscard]] Geo::SigmaTau ReturnStressAtRegularFailureZone(const Geo::SigmaTau& rSigmaTau,
                                                                 CoulombYieldSurface::CoulombAveragingType AveragingType) const;
    [[nodiscard]] Geo::PrincipalStresses ReturnStressAtRegularFailureZone(
        const Geo::PrincipalStresses&             rPrincipalStresses,
        CoulombYieldSurface::CoulombAveragingType AveragingType) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const;
    void load(Serializer& rSerializer);
};

} // namespace Kratos
