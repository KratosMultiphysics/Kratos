//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/condition.h"
#include "utilities/input_output_utilities.h"

namespace Kratos
{

template<typename TEntityType>
bool InputOutputUtilities::SkippableEntity(
    const TEntityType& rEntity,
    const std::string& rWarningLabel
    )
{
    constexpr bool is_element = std::is_same_v<TEntityType, Element>;
    constexpr bool is_condition = std::is_same_v<TEntityType, Condition>;
    if constexpr (is_element || is_condition) {
        const auto& r_geometry = rEntity.GetGeometry();
        switch (r_geometry.GetGeometryType()) {
            case GeometryData::KratosGeometryType::Kratos_Nurbs_Curve:
                KRATOS_WARNING(rWarningLabel) << "Skipping geometry type: Kratos_Nurbs_Curve" << std::endl;
                return true;
            case GeometryData::KratosGeometryType::Kratos_Nurbs_Surface:
                KRATOS_WARNING(rWarningLabel) << "Skipping geometry type: Kratos_Nurbs_Surface" << std::endl;
                return true;
            case GeometryData::KratosGeometryType::Kratos_Nurbs_Volume:
                KRATOS_WARNING(rWarningLabel) << "Skipping geometry type: Kratos_Nurbs_Volume" << std::endl;
                return true;
            case GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface:
                KRATOS_WARNING(rWarningLabel) << "Skipping geometry type: Kratos_Nurbs_Curve_On_Surface" << std::endl;
                return true;
            case GeometryData::KratosGeometryType::Kratos_Surface_In_Nurbs_Volume:
                KRATOS_WARNING(rWarningLabel) << "Skipping geometry type: Kratos_Surface_In_Nurbs_Volume" << std::endl;
                return true;
            case GeometryData::KratosGeometryType::Kratos_Brep_Curve:
                KRATOS_WARNING(rWarningLabel) << "Skipping geometry type: Kratos_Brep_Curve" << std::endl;
                return true;
            case GeometryData::KratosGeometryType::Kratos_Brep_Surface:
                KRATOS_WARNING(rWarningLabel) << "Skipping geometry type: Kratos_Brep_Surface" << std::endl;
                return true;
            case GeometryData::KratosGeometryType::Kratos_Brep_Curve_On_Surface:
                KRATOS_WARNING(rWarningLabel) << "Skipping geometry type: Kratos_Brep_Curve_On_Surface" << std::endl;
                return true;
            case GeometryData::KratosGeometryType::Kratos_Coupling_Geometry:
                KRATOS_WARNING(rWarningLabel) << "Skipping geometry type: Kratos_Coupling_Geometry" << std::endl;
                return true;
            case GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry:
                KRATOS_WARNING(rWarningLabel) << "Skipping geometry type: Kratos_Quadrature_Point_Curve_On_Surface_Geometry" << std::endl;
                return true;
            case GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Surface_In_Volume_Geometry:
                KRATOS_WARNING(rWarningLabel) << "Skipping geometry type: Kratos_Quadrature_Point_Surface_In_Volume_Geometry" << std::endl;
                return true;
            default:
                // For any other geometry type, do nothing and continue execution
                break;
        }
    }

    return false;
}

/// Explicit instantation for Element and Condition (also Node as is used in VtkOutput)
template bool InputOutputUtilities::SkippableEntity<Node>(const Node&, const std::string&);
template bool InputOutputUtilities::SkippableEntity<Element>(const Element&, const std::string&);
template bool InputOutputUtilities::SkippableEntity<Condition>(const Condition&, const std::string&);

} // namespace Kratos