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

namespace Formulations
{
struct Constant {
    static constexpr const char* Name = "Constant";
    double                       operator()(const Properties&, double, Vector&) const;
};

struct Eur {
    static constexpr const char* Name = "Eur";
    double                       operator()(const Properties&, double, Vector&) const;
};

using YoungsModulusVariant = std::variant<Constant, Eur>;

YoungsModulusVariant InitializeFormulation(const std::string& rFormulation);
std::string          GetYoungsModulusFormulation(const Properties& rProperties);
std::string          GetYoungsModulusFormulation(const YoungsModulusVariant& rFormulation);
double               GetYoungsModulus(YoungsModulusVariant& rFormulation,
                                      const Properties&     rProperties,
                                      double                YoungsModulus,
                                      Vector&               rStressVectorFinalized);
double CalculateYoungsModulusForEur(const Properties& rProperties, double YoungsModulus, Vector& rStressVectorFinalized);
double CalculateMinorPrincipalEffectiveStress(Vector& rStressVectorFinalized);
} // namespace Formulations
} // namespace Kratos
