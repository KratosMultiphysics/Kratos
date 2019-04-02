//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "geometries/point.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_3d_6.h"
// #include "geometries/prism_3d_15.h"
#include "geometries/hexahedra_3d_8.h"
// #include "geometries/hexahedra_3d_20.h"
// #include "geometries/hexahedra_3d_27.h"
#include "python/add_geometries_to_python.h"

namespace Kratos
{

namespace Python
{
    typedef Geometry<Node<3> > GeometryType;
    typedef GeometryType::PointsArrayType NodesArrayType;
    typedef GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef Point::CoordinatesArrayType CoordinatesArrayType;

    const PointerVector< Node<3> >& ConstGetPoints( GeometryType& geom ) { return geom.Points(); }
    PointerVector< Node<3> >& GetPoints( GeometryType& geom ) { return geom.Points(); }

    array_1d<double,3> GetNormal(
        GeometryType& dummy,
        CoordinatesArrayType& LocalCoords
        )
    {
        return( dummy.Normal(LocalCoords) );
    }

    array_1d<double,3> FastGetNormal(GeometryType& dummy)
    {
        CoordinatesArrayType LocalCoords;
        LocalCoords.clear();
        return( dummy.Normal(LocalCoords) );
    }

    array_1d<double,3> GetUnitNormal(
        GeometryType& dummy,
        CoordinatesArrayType& LocalCoords
        )
    {
        return( dummy.UnitNormal(LocalCoords) );
    }

    array_1d<double,3> FastGetUnitNormal(GeometryType& dummy)
    {
        CoordinatesArrayType LocalCoords;
        LocalCoords.clear();
        return( dummy.UnitNormal(LocalCoords) );
    }

void  AddGeometriesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef Node<3> NodeType;
    typedef NodeType::Pointer pNodeType;
    typedef Geometry<NodeType > GeometryType;

    py::class_<GeometryType, GeometryType::Pointer >(m,"Geometry")
    .def(py::init<>())
    .def(py::init< GeometryType::PointsArrayType& >())
    .def("WorkingSpaceDimension",&GeometryType::WorkingSpaceDimension)
    .def("LocalSpaceDimension",&GeometryType::LocalSpaceDimension)
    .def("DomainSize",&GeometryType::DomainSize)
    .def("PointsNumber",&GeometryType::PointsNumber)
    .def("Normal",GetNormal)
    .def("Normal",FastGetNormal)
    .def("UnitNormal",GetUnitNormal)
    .def("UnitNormal",FastGetUnitNormal)
    .def("Center",&GeometryType::Center)
    .def("Length",&GeometryType::Length)
    .def("Area",&GeometryType::Area)
    .def("Volume",&GeometryType::Volume)
    .def("__str__", PrintObject<GeometryType>)
//     .def("Points", &GeometryType::ConstGetPoints)
//     .def("Points", &GeometryType::GetPoints)
    ;

    // 2D
    py::class_<Line2D2<NodeType>, Line2D2<NodeType>::Pointer,  GeometryType  >(m,"Line2D2").def(py::init<pNodeType, pNodeType>())
    ;
    py::class_<Line2D3<NodeType>, Line2D3<NodeType>::Pointer,  GeometryType  >(m,"Line2D3").def(py::init<pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Triangle2D3<NodeType>, Triangle2D3<NodeType>::Pointer,  GeometryType  >(m,"Triangle2D3").def(py::init<pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Triangle2D6<NodeType>, Triangle2D6<NodeType>::Pointer,  GeometryType  >(m,"Triangle2D6").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Quadrilateral2D4<NodeType>, Quadrilateral2D4<NodeType>::Pointer,  GeometryType  >(m,"Quadrilateral2D4").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Quadrilateral2D8<NodeType>, Quadrilateral2D8<NodeType>::Pointer,  GeometryType  >(m,"Quadrilateral2D8").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Quadrilateral2D9<NodeType>, Quadrilateral2D9<NodeType>::Pointer,  GeometryType  >(m,"Quadrilateral2D9").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;

    // 3D
    py::class_<Line3D2<NodeType>, Line3D2<NodeType>::Pointer,  GeometryType  >(m,"Line3D2").def(py::init<pNodeType, pNodeType>())
    ;
    py::class_<Line3D3<NodeType>, Line3D3<NodeType>::Pointer,  GeometryType  >(m,"Line3D3").def(py::init<pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Triangle3D3<NodeType>, Triangle3D3<NodeType>::Pointer,  GeometryType  >(m,"Triangle3D3").def(py::init<pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Triangle3D6<NodeType>, Triangle3D6<NodeType>::Pointer,  GeometryType  >(m,"Triangle3D6").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Quadrilateral3D4<NodeType>, Quadrilateral3D4<NodeType>::Pointer,  GeometryType  >(m,"Quadrilateral3D4").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Quadrilateral3D8<NodeType>, Quadrilateral3D8<NodeType>::Pointer,  GeometryType  >(m,"Quadrilateral3D8").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Quadrilateral3D9<NodeType>, Quadrilateral3D9<NodeType>::Pointer,  GeometryType  >(m,"Quadrilateral3D9").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Tetrahedra3D4<NodeType>, Tetrahedra3D4<NodeType>::Pointer,  GeometryType  >(m,"Tetrahedra3D4").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Tetrahedra3D10<NodeType>, Tetrahedra3D10<NodeType>::Pointer,  GeometryType  >(m,"Tetrahedra3D10").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    py::class_<Prism3D6<NodeType>, Prism3D6<NodeType>::Pointer,  GeometryType  >(m,"Prism3D6").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
//     py::class_<Prism3D15<NodeType>, Prism3D15<NodeType>::Pointer,  GeometryType  >(m,"Prism3D15").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
//     ;
    py::class_<Hexahedra3D8<NodeType>, Hexahedra3D8<NodeType>::Pointer,  GeometryType  >(m,"Hexahedra3D8").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
//     py::class_<Hexahedra3D20<NodeType>, Hexahedra3D20<NodeType>::Pointer,  GeometryType  >(m,"Hexahedra3D20").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
//     ;
//     py::class_<Hexahedra3D27<NodeType>, Hexahedra3D27<NodeType>::Pointer,  GeometryType  >(m,"Hexahedra3D27").def(py::init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
//     ;

}

}  // namespace Python.

} // Namespace Kratos

