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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
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
#include "python/bounded_vector_python_interface.h"
#include "python/vector_scalar_operator_python.h"
#include "python/vector_vector_operator_python.h"

namespace Kratos
{

namespace Python
{
    const PointerVector< Node<3> >& ConstGetPoints( Geometry<Node<3> >& geom ) { return geom.Points(); }
    PointerVector< Node<3> >& GetPoints( Geometry<Node<3> >& geom ) { return geom.Points(); }
void  AddGeometriesToPython()
{

    typedef Node<3> NodeType;
    typedef NodeType::Pointer pNodeType;
    typedef Geometry<NodeType > GeometryType;

    class_<GeometryType, GeometryType::Pointer >("Geometry", init<>())
    .def(init< GeometryType::PointsArrayType& >())
    .def("DomainSize",&GeometryType::DomainSize)
    .def("PointsNumber",&GeometryType::PointsNumber)
    .def("AreaNormal",&GeometryType::AreaNormal)
    .def("UnitNormal",&GeometryType::UnitNormal)
    .def("Center",&GeometryType::Center)
    .def("Length",&GeometryType::Length)
    .def("Area",&GeometryType::Area)
    .def("Volume",&GeometryType::Volume)
    .def("Info",&GeometryType::Info)
//     .def("Points", &GeometryType::ConstGetPoints)
//     .def("Points", &GeometryType::GetPoints)
    ;
    
    // 2D
    class_<Line2D2<NodeType>, Line2D2<NodeType>::Pointer, bases< GeometryType > >("Line2D2", init<pNodeType, pNodeType>())
    ;
    class_<Line2D3<NodeType>, Line2D3<NodeType>::Pointer, bases< GeometryType > >("Line2D3", init<pNodeType, pNodeType, pNodeType>())
    ;
    class_<Triangle2D3<NodeType>, Triangle2D3<NodeType>::Pointer, bases< GeometryType > >("Triangle2D3", init<pNodeType, pNodeType, pNodeType>())
    ;
    class_<Triangle2D6<NodeType>, Triangle2D6<NodeType>::Pointer, bases< GeometryType > >("Triangle2D6", init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    class_<Quadrilateral2D4<NodeType>, Quadrilateral2D4<NodeType>::Pointer, bases< GeometryType > >("Quadrilateral2D4", init<pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    class_<Quadrilateral2D8<NodeType>, Quadrilateral2D8<NodeType>::Pointer, bases< GeometryType > >("Quadrilateral2D8", init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    class_<Quadrilateral2D9<NodeType>, Quadrilateral2D9<NodeType>::Pointer, bases< GeometryType > >("Quadrilateral2D9", init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;

    // 3D
    class_<Line3D2<NodeType>, Line3D2<NodeType>::Pointer, bases< GeometryType > >("Line3D2", init<pNodeType, pNodeType>())
    ;
    class_<Line3D3<NodeType>, Line3D3<NodeType>::Pointer, bases< GeometryType > >("Line3D3", init<pNodeType, pNodeType, pNodeType>())
    ;
    class_<Triangle3D3<NodeType>, Triangle3D3<NodeType>::Pointer, bases< GeometryType > >("Triangle3D3", init<pNodeType, pNodeType, pNodeType>())
    ;
    class_<Triangle3D6<NodeType>, Triangle3D6<NodeType>::Pointer, bases< GeometryType > >("Triangle3D6", init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    class_<Quadrilateral3D4<NodeType>, Quadrilateral3D4<NodeType>::Pointer, bases< GeometryType > >("Quadrilateral3D4", init<pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    class_<Quadrilateral3D8<NodeType>, Quadrilateral3D8<NodeType>::Pointer, bases< GeometryType > >("Quadrilateral3D8", init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    class_<Quadrilateral3D9<NodeType>, Quadrilateral3D9<NodeType>::Pointer, bases< GeometryType > >("Quadrilateral3D9", init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    class_<Tetrahedra3D4<NodeType>, Tetrahedra3D4<NodeType>::Pointer, bases< GeometryType > >("Tetrahedra3D4", init<pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    class_<Tetrahedra3D10<NodeType>, Tetrahedra3D10<NodeType>::Pointer, bases< GeometryType > >("Tetrahedra3D10", init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
    class_<Prism3D6<NodeType>, Prism3D6<NodeType>::Pointer, bases< GeometryType > >("Prism3D6", init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
//     class_<Prism3D15<NodeType>, Prism3D15<NodeType>::Pointer, bases< GeometryType > >("Prism3D15", init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
//     ;
    class_<Hexahedra3D8<NodeType>, Hexahedra3D8<NodeType>::Pointer, bases< GeometryType > >("Hexahedra3D8", init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
    ;
//     class_<Hexahedra3D20<NodeType>, Hexahedra3D20<NodeType>::Pointer, bases< GeometryType > >("Hexahedra3D20", init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
//     ;
//     class_<Hexahedra3D27<NodeType>, Hexahedra3D27<NodeType>::Pointer, bases< GeometryType > >("Hexahedra3D27", init<pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType, pNodeType>())
//     ;
     
}

}  // namespace Python.

} // Namespace Kratos

