//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// Include base h
#include "vtk_definitions.h"

namespace Kratos {

const std::map<GeometryData::KratosGeometryType, char> VtkDefinitions::KratosVtkGeometryTypes = {
        {GeometryData::KratosGeometryType::Kratos_Point2D, 1},
        {GeometryData::KratosGeometryType::Kratos_Point3D, 1},
        {GeometryData::KratosGeometryType::Kratos_Line2D2, 3},
        {GeometryData::KratosGeometryType::Kratos_Line3D2, 3},
        {GeometryData::KratosGeometryType::Kratos_Triangle2D3, 5},
        {GeometryData::KratosGeometryType::Kratos_Triangle3D3, 5},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, 9},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, 9},
        {GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4, 10},
        {GeometryData::KratosGeometryType::Kratos_Hexahedra3D8, 12},
        {GeometryData::KratosGeometryType::Kratos_Prism3D6, 13},
        {GeometryData::KratosGeometryType::Kratos_Pyramid3D5, 14},
        {GeometryData::KratosGeometryType::Kratos_Line2D3, 21},
        {GeometryData::KratosGeometryType::Kratos_Line3D3, 21},
        {GeometryData::KratosGeometryType::Kratos_Triangle2D6, 22},
        {GeometryData::KratosGeometryType::Kratos_Triangle3D6, 22},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8, 23},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, 23},
        {GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10, 24},
        {GeometryData::KratosGeometryType::Kratos_Hexahedra3D20, 25},
        {GeometryData::KratosGeometryType::Kratos_Prism3D15, 26},
        {GeometryData::KratosGeometryType::Kratos_Pyramid3D13, 27},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9, 28},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9, 28},
        {GeometryData::KratosGeometryType::Kratos_Hexahedra3D27, 29}};

} // namespace Kratos