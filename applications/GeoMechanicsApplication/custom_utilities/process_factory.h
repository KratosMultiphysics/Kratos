// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#pragma once

#include "processes/process.h"

#include <functional>
#include <map>
#include <memory>
#include <string>

namespace Kratos
{

class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ProcessFactory
{
public:
    using ProductType = std::unique_ptr<Process>;

    [[nodiscard]] ProductType Create(const std::string& rProcessClassName, const Parameters& rProcessSettings) const;
    void AddCreator(const std::string& rProcessClassName, std::function<ProductType(const Parameters&)> Creator);

    void SetCallBackWhenProcessIsUnknown(const std::function<void(const std::string&)>& function);

private:
    std::map<std::string, std::function<ProductType(const Parameters&)>, std::less<>> mCreatorMap;
    std::function<void(const std::string&)> mCallBackIfProcessIsUnknown;
};

} // namespace Kratos
