//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// System includes

// External includes

// Project includes
#include "geometries/geometry_data.h"
#include "add_geometry_data_to_python.h"

namespace Kratos::Python
{

void AddGeometryDataToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto integration_method = py::enum_<GeometryData::IntegrationMethod>(m, "GeometryData_IntegrationMethod")
        .value("GI_GAUSS_1", GeometryData::IntegrationMethod::GI_GAUSS_1)
        .value("GI_GAUSS_2", GeometryData::IntegrationMethod::GI_GAUSS_2)
        .value("GI_GAUSS_3", GeometryData::IntegrationMethod::GI_GAUSS_3)
        .value("GI_GAUSS_4", GeometryData::IntegrationMethod::GI_GAUSS_4)
        .value("GI_GAUSS_5", GeometryData::IntegrationMethod::GI_GAUSS_5)
        .value("GI_EXTENDED_GAUSS_1", GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1)
        .value("GI_EXTENDED_GAUSS_2", GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2)
        .value("GI_EXTENDED_GAUSS_3", GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_3)
        .value("GI_EXTENDED_GAUSS_4", GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_4)
        .value("GI_EXTENDED_GAUSS_5", GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_5);

    auto geometry_family = py::enum_<GeometryData::KratosGeometryFamily>(m, "GeometryData_KratosGeometryFamily")
        .value("Kratos_NoElement", GeometryData::KratosGeometryFamily::Kratos_NoElement)
        .value("Kratos_Point", GeometryData::KratosGeometryFamily::Kratos_Point)
        .value("Kratos_Linear", GeometryData::KratosGeometryFamily::Kratos_Linear)
        .value("Kratos_Triangle", GeometryData::KratosGeometryFamily::Kratos_Triangle)
        .value("Kratos_Quadrilateral", GeometryData::KratosGeometryFamily::Kratos_Quadrilateral)
        .value("Kratos_Tetrahedra", GeometryData::KratosGeometryFamily::Kratos_Tetrahedra)
        .value("Kratos_Hexahedra", GeometryData::KratosGeometryFamily::Kratos_Hexahedra)
        .value("Kratos_Prism", GeometryData::KratosGeometryFamily::Kratos_Prism)
        .value("Kratos_Pyramid", GeometryData::KratosGeometryFamily::Kratos_Pyramid)
        .value("Kratos_Nurbs", GeometryData::KratosGeometryFamily::Kratos_Nurbs)
        .value("Kratos_Brep", GeometryData::KratosGeometryFamily::Kratos_Brep)
        .value("Kratos_Quadrature_Geometry", GeometryData::KratosGeometryFamily::Kratos_Quadrature_Geometry)
        .value("Kratos_Composite", GeometryData::KratosGeometryFamily::Kratos_Composite)
        .value("Kratos_generic_family", GeometryData::KratosGeometryFamily::Kratos_generic_family);

    auto geometry_type = py::enum_<GeometryData::KratosGeometryType>(m, "GeometryData_KratosGeometryType")
        .value("Kratos_generic_type", GeometryData::KratosGeometryType::Kratos_generic_type)
        .value("Kratos_Hexahedra3D20", GeometryData::KratosGeometryType::Kratos_Hexahedra3D20)
        .value("Kratos_Hexahedra3D27", GeometryData::KratosGeometryType::Kratos_Hexahedra3D27)
        .value("Kratos_Hexahedra3D8", GeometryData::KratosGeometryType::Kratos_Hexahedra3D8)
        .value("Kratos_Prism3D15", GeometryData::KratosGeometryType::Kratos_Prism3D15)
        .value("Kratos_Prism3D6", GeometryData::KratosGeometryType::Kratos_Prism3D6)
        .value("Kratos_Pyramid3D13", GeometryData::KratosGeometryType::Kratos_Pyramid3D13)
        .value("Kratos_Pyramid3D5", GeometryData::KratosGeometryType::Kratos_Pyramid3D5)
        .value("Kratos_Quadrilateral2D4", GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
        .value("Kratos_Quadrilateral2D8", GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8)
        .value("Kratos_Quadrilateral2D9", GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9)
        .value("Kratos_Quadrilateral3D4", GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)
        .value("Kratos_Quadrilateral3D8", GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8)
        .value("Kratos_Quadrilateral3D9", GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9)
        .value("Kratos_Tetrahedra3D10", GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10)
        .value("Kratos_Tetrahedra3D4", GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
        .value("Kratos_Triangle2D3", GeometryData::KratosGeometryType::Kratos_Triangle2D3)
        .value("Kratos_Triangle2D6", GeometryData::KratosGeometryType::Kratos_Triangle2D6)
        .value("Kratos_Triangle2D10", GeometryData::KratosGeometryType::Kratos_Triangle2D10)
        .value("Kratos_Triangle2D15", GeometryData::KratosGeometryType::Kratos_Triangle2D15)
        .value("Kratos_Triangle3D3", GeometryData::KratosGeometryType::Kratos_Triangle3D3)
        .value("Kratos_Triangle3D6", GeometryData::KratosGeometryType::Kratos_Triangle3D6)
        .value("Kratos_Line2D2", GeometryData::KratosGeometryType::Kratos_Line2D2)
        .value("Kratos_Line2D3", GeometryData::KratosGeometryType::Kratos_Line2D3)
        .value("Kratos_Line2D4", GeometryData::KratosGeometryType::Kratos_Line2D4)
        .value("Kratos_Line2D5", GeometryData::KratosGeometryType::Kratos_Line2D5)
        .value("Kratos_Line3D2", GeometryData::KratosGeometryType::Kratos_Line3D2)
        .value("Kratos_Line3D3", GeometryData::KratosGeometryType::Kratos_Line3D3)
        .value("Kratos_Point2D", GeometryData::KratosGeometryType::Kratos_Point2D)
        .value("Kratos_Point3D", GeometryData::KratosGeometryType::Kratos_Point3D)
        .value("Kratos_Sphere3D1", GeometryData::KratosGeometryType::Kratos_Sphere3D1)
        .value("Kratos_Nurbs_Curve", GeometryData::KratosGeometryType::Kratos_Nurbs_Curve)
        .value("Kratos_Nurbs_Surface", GeometryData::KratosGeometryType::Kratos_Nurbs_Surface)
        .value("Kratos_Nurbs_Volume", GeometryData::KratosGeometryType::Kratos_Nurbs_Volume)
        .value("Kratos_Nurbs_Curve_On_Surface", GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface)
        .value("Kratos_Surface_In_Nurbs_Volume", GeometryData::KratosGeometryType::Kratos_Surface_In_Nurbs_Volume)
        .value("Kratos_Brep_Curve", GeometryData::KratosGeometryType::Kratos_Brep_Curve)
        .value("Kratos_Brep_Surface", GeometryData::KratosGeometryType::Kratos_Brep_Surface)
        .value("Kratos_Brep_Curve_On_Surface", GeometryData::KratosGeometryType::Kratos_Brep_Curve_On_Surface)
        .value("Kratos_Quadrature_Point_Geometry", GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Geometry)
        .value("Kratos_Coupling_Geometry", GeometryData::KratosGeometryType::Kratos_Coupling_Geometry)
        .value("Kratos_Quadrature_Point_Curve_On_Surface_Geometry", GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry)
        .value("Kratos_Quadrature_Point_Surface_In_Volume_Geometry", GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Surface_In_Volume_Geometry);

    py::class_<GeometryData, GeometryData::Pointer>(m, "GeometryData")
        .def_property_readonly_static("IntegrationMethod", [integration_method](py::object) { return integration_method; })
        .def_property_readonly_static("KratosGeometryFamily", [geometry_family](py::object) { return geometry_family; })
        .def_property_readonly_static("KratosGeometryType", [geometry_type](py::object) { return geometry_type; });

}

}