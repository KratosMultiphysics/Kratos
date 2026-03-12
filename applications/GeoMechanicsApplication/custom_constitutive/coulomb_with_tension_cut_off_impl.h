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
#include "custom_constitutive/principal_stresses.hpp"
#include "geo_mechanics_application_constants.h"
#include "tension_cutoff.h"

namespace Kratos
{

class Properties;
class Serializer;

namespace Geo
{
class SigmaTau;
} // namespace Geo

class CoulombWithTensionCutOffImpl
{
public:
    CoulombWithTensionCutOffImpl() = default;
    explicit CoulombWithTensionCutOffImpl(const Properties& rMaterialProperties);

    [[nodiscard]] bool IsAdmissibleStressState(const Geo::SigmaTau& rTrialTraction) const;
    [[nodiscard]] bool IsAdmissibleStressState(const Geo::PrincipalStresses& rTrialPrincipalStresses) const;

    [[nodiscard]] Geo::SigmaTau DoReturnMapping(const Geo::SigmaTau& rTrialTraction,
                                                const Matrix&        rElasticConstitutiveTensor,
                                                Geo::PrincipalStresses::AveragingType AveragingType);
    [[nodiscard]] Geo::PrincipalStresses DoReturnMapping(const Geo::PrincipalStresses& rTrialSigmaTau,
                                                         const Matrix& rElasticConstitutiveTensor,
                                                         Geo::PrincipalStresses::AveragingType AveragingType);

    void                           SaveKappaOfCoulombYieldSurface();
    void                           RestoreKappaOfCoulombYieldSurface();
    [[nodiscard]] PlasticityStatus GetPlasticityStatus() const;

private:
    CoulombYieldSurface mCoulombYieldSurface;
    TensionCutoff       mTensionCutOff;
    double              mSavedKappaOfCoulombYieldSurface{0.0};
    double              mAbsoluteYieldFunctionValueTolerance{1.0e-8};
    std::size_t         mMaxNumberOfPlasticIterations{100};
    PlasticityStatus    mPlasticityStatus{PlasticityStatus::ELASTIC};

    template <typename StressStateType>
    [[nodiscard]] bool IsAdmissibleStressState(const StressStateType& rTrialStressState) const;
    template <typename StressStateType, typename StressStateToSigmaTauFunctionType>
    [[nodiscard]] StressStateType DoReturnMapping(const StressStateType& rTrialStressState,
                                                  const StressStateToSigmaTauFunctionType& rStressStateToSigmaTau,
                                                  const Matrix& rElasticConstitutiveTensor,
                                                  Geo::PrincipalStresses::AveragingType AveragingType);

    [[nodiscard]] Geo::SigmaTau CalculateCornerPoint() const;
    [[nodiscard]] bool IsStressAtTensionApexReturnZone(const Geo::SigmaTau& rTrialTraction) const;
    [[nodiscard]] bool IsStressAtTensionCutoffReturnZone(const Geo::SigmaTau& rTrialTraction) const;
    [[nodiscard]] bool IsStressAtCornerReturnZone(const Geo::SigmaTau& rTrialTraction,
                                                  Geo::PrincipalStresses::AveragingType AveragingType) const;

    [[nodiscard]] Geo::SigmaTau ReturnStressAtTensionApexReturnZone(const Geo::SigmaTau&) const;
    [[nodiscard]] Geo::PrincipalStresses ReturnStressAtTensionApexReturnZone(const Geo::PrincipalStresses& rTrialPrincipalStresses) const;
    [[nodiscard]] Geo::SigmaTau ReturnStressAtTensionCutoffReturnZone(const Geo::SigmaTau& rTrialTraction,
                                                                      const Matrix& rElasticConstitutiveTensor,
                                                                      Geo::PrincipalStresses::AveragingType AveragingType) const;
    [[nodiscard]] Geo::PrincipalStresses ReturnStressAtTensionCutoffReturnZone(
        const Geo::PrincipalStresses& rTrialPrincipalStresses,
        const Matrix&                 rElasticConstitutiveTensor,
        Geo::PrincipalStresses::AveragingType) const;
    [[nodiscard]] Geo::SigmaTau ReturnStressAtRegularFailureZone(const Geo::SigmaTau& rTrialTraction,
                                                                 const Matrix& rElasticConstitutiveTensor,
                                                                 Geo::PrincipalStresses::AveragingType AveragingType) const;
    [[nodiscard]] Geo::PrincipalStresses ReturnStressAtRegularFailureZone(
        const Geo::PrincipalStresses&         rTrialPrincipalStresses,
        const Matrix&                         rElasticConstitutiveTensor,
        Geo::PrincipalStresses::AveragingType AveragingType) const;
    [[nodiscard]] Geo::PrincipalStresses ReturnStressAtCornerPoint(const Geo::PrincipalStresses& rTrialPrincipalStresses,
                                                                   const Matrix& rElasticConstitutiveTensor,
                                                                   Geo::PrincipalStresses::AveragingType AveragingType) const;
    [[nodiscard]] Geo::SigmaTau ReturnStressAtCornerPoint(const Geo::SigmaTau&,
                                                          const Matrix&,
                                                          Geo::PrincipalStresses::AveragingType AveragingType) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const;
    void load(Serializer& rSerializer);
};

} // namespace Kratos
