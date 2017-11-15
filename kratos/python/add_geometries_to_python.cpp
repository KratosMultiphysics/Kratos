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
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
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

    typedef Node<3> Node3D;
    typedef Node3D::Pointer pNode3D;
    typedef Geometry<Node3D > GeometryType;

    class_<GeometryType, GeometryType::Pointer >("Geometry", init<>())
    .def(init< GeometryType::PointsArrayType& >())
	.def("DomainSize",&GeometryType::DomainSize)
//     .def("Points", &GeometryType::ConstGetPoints)
//     .def("Points", &GeometryType::GetPoints)
    ;
    
    class_<Triangle2D3<Node3D>, Triangle2D3<Node3D>::Pointer, bases< GeometryType > >("Triangle2D3", init<pNode3D, pNode3D, pNode3D>())
    ;    
    
    class_<Tetrahedra3D4<Node3D>, Tetrahedra3D4<Node3D>::Pointer, bases< GeometryType > >("Tetrahedra3D4", init<pNode3D, pNode3D, pNode3D, pNode3D>())
    ;

    class_<Hexahedra3D8<Node3D>, Hexahedra3D8<Node3D>::Pointer, bases< GeometryType > >("Hexahedra3D8", init<pNode3D, pNode3D, pNode3D, pNode3D, pNode3D, pNode3D, pNode3D, pNode3D>())
    ;

//     class_<GeometryType, GeometryType::Pointer, bases<PointerVector< Node<3> > > >("Geometry", init<>())
//      .def(init< GeometryType::PointsArrayType& >())
//      ;
     
}

}  // namespace Python.

} // Namespace Kratos

