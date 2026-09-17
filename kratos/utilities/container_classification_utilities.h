#pragma once

#include "includes/model_part.h"
#include "containers/pointer_vector.h"
#include "geometries/geometry_data.h"
#include <vector>

namespace Kratos {

class KRATOS_API(KRATOS_CORE) ContainerClassificationUtilities {
public:
    KRATOS_CLASS_POINTER_DEFINITION(ContainerClassificationUtilities);

    ContainerClassificationUtilities() = default;
    ~ContainerClassificationUtilities() = default;

    /// Classifies Elements/Conditions/Geometries using both Entity dynamic type (typeid) and GeometryType enum
    template <typename TContainerType>
    std::vector<TContainerType> Classify(TContainerType& rContainer);

    /// Classifies Elements/Conditions/Geometries strictly by GeometryType enum
    /// Returns a std::vector<TContainerType> sized exactly to GeometryData::KratosGeometryType::NumberOfGeometryTypes
    template <typename TContainerType>
    std::vector<TContainerType> ClassifyByGeometryType(TContainerType& rContainer);
};

} // namespace Kratos