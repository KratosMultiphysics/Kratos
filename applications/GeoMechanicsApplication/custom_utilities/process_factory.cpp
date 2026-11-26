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

#include "custom_utilities/process_factory.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

ProcessFactory::ProductType ProcessFactory::Create(const std::string& rProcessClassName,
                                                   const Parameters&  rProcessSettings) const
{
    auto pos = mCreatorMap.find(rProcessClassName);
    if (pos == mCreatorMap.end()) {
        if (mCallBackIfProcessIsUnknown != nullptr) {
            mCallBackIfProcessIsUnknown(rProcessClassName);
        }

        return nullptr;
    }

    return pos->second ? pos->second(rProcessSettings) : nullptr;
}

void ProcessFactory::AddCreator(const std::string&                            rProcessClassName,
                                std::function<ProductType(const Parameters&)> Creator)
{
    mCreatorMap[rProcessClassName] = std::move(Creator);
}

void ProcessFactory::SetCallBackWhenProcessIsUnknown(const std::function<void(const std::string&)>& function)
{
    mCallBackIfProcessIsUnknown = function;
}

} // namespace Kratos
