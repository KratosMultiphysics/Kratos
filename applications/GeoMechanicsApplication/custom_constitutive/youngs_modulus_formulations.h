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
        static constexpr const char* Name = "CONSTANT";
        double                       operator()(const Properties&, double, const Vector&) const;
    };

    struct SchanzVermeer {
        static constexpr const char* Name = "SCHANZ_VERMEER";
        double                       operator()(const Properties&, double, const Vector&) const;
    };

    using YoungsModulusVariant = std::variant<Constant, SchanzVermeer>;

    static void                 CheckInputData(const Properties& rMaterialProperties);
    static YoungsModulusVariant InitializeFormulation(const std::string& rFormulation);
    static std::string          GetYoungsModulusFormulation(const Properties& rProperties);
    static std::string GetYoungsModulusFormulation(const YoungsModulusVariant& rFormulation);
    static double      GetYoungsModulus(YoungsModulusVariant& rFormulation,
                                        const Properties&     rProperties,
                                        double                ReferenceYoungsModulus,
                                        const Vector&         rStressVectorFinalized);
    static double      CalculateYoungsModulusSchanzVermeer(const Properties& rProperties,
                                                           double            ReferenceYoungsModulus,
                                                           const Vector&     rStressVectorFinalized);
    static double      CalculateMinorPrincipalEffectiveStress(const Vector& rStressVectorFinalized);
};
} // namespace Kratos
