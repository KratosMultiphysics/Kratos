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
//  Main authors:    Mohamed Nabi,
//                   Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#pragma once

#include "custom_constitutive/p_q.hpp"
#include "geo_aliases.h"
#include "includes/properties.h"

namespace Kratos
{

class CheckProperties;

class KRATOS_API(GEO_MECHANICS_APPLICATION) CompressionCapYieldSurface
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CompressionCapYieldSurface);

    CompressionCapYieldSurface();
    explicit CompressionCapYieldSurface(const Properties& rMaterialProperties);

    [[nodiscard]] double GetCapSize() const;
    [[nodiscard]] double GetPreconsolidationStress() const;

    [[nodiscard]] double YieldFunctionValue(const Geo::PQ& rPQ) const;
    [[nodiscard]] Vector DerivativeOfFlowFunction(const Geo::PQ& rPQ) const;

private:
    void InitializeKappaDependentFunctions();
    void CheckMaterialProperties() const;

    friend class Serializer;
    void save(Serializer& rSerializer) const;
    void load(Serializer& rSerializer);

    double                      mKappa = 0.0;
    Properties                  mMaterialProperties;
    Geo::KappaDependentFunction mCapSizeCalculator;
    Geo::KappaDependentFunction mPreconsolidationStressCalculator;

}; // Class CompressionCapYieldSurface

} // namespace Kratos