// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#pragma once

#include "includes/properties.h"

#include <string>
#include <variant>

namespace Kratos
{
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoYoungsModulusFormulations
{
public:
    struct Constant {
        static constexpr const char* Name = "Constant";
        double                       operator()(const Properties&, double, const Vector&) const;
    };

    struct Eur {
        static constexpr const char* Name = "Eur";
        double                       operator()(const Properties&, double, const Vector&) const;
    };

    using YoungsModulusVariant = std::variant<Constant, Eur>;

    static void                 CheckInputData(const Properties& rMaterialProperties);
    static YoungsModulusVariant InitializeFormulation(const std::string& rFormulation);
    static std::string          GetYoungsModulusFormulation(const Properties& rProperties);
    static std::string GetYoungsModulusFormulation(const YoungsModulusVariant& rFormulation);
    static double      GetYoungsModulus(YoungsModulusVariant& rFormulation,
                                        const Properties&     rProperties,
                                        double                YoungsModulus,
                                        const Vector&         rStressVectorFinalized);
    static double      CalculateYoungsModulusForEur(const Properties& rProperties,
                                                    double            YoungsModulus,
                                                    const Vector&     rStressVectorFinalized);
    static double      CalculateMinorPrincipalEffectiveStress(const Vector& rStressVectorFinalized);
};
} // namespace Kratos
