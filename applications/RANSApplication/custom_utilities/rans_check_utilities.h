//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_APPLICATION_CHECK_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_CHECK_UTILITIES_H_INCLUDED

// System includes
#include <string>

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
namespace RansCheckUtilities
{
bool CheckIfModelPartExists(const Model& rModel, const std::string& rModelPartName);

template <typename TVariableType>
bool CheckIfVariableExistsInModelPart(const ModelPart& rModelPart,
                                      const TVariableType& rVariable);
} // namespace RansCheckUtilities

} // namespace Kratos

#endif // KRATOS_RANS_APPLICATION_CHECK_UTILITIES_H_INCLUDED defined