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
#include "python/add_points_to_python.h"


namespace Kratos
{

namespace Python
{
template<class TPointType> inline
void PointSetX(TPointType& ThisPoint, double Value)
{
    ThisPoint.X() = Value;
}

template<class TPointType> inline
void PointSetY(TPointType& ThisPoint, double Value)
{
    ThisPoint.Y() = Value;
}

template<class TPointType> inline
void PointSetZ(TPointType& ThisPoint, double Value)
{
    ThisPoint.Z() = Value;
}

template<class TPointType> inline
double PointGetX(TPointType& ThisPoint)
{
    return ThisPoint.X();
}

template<class TPointType> inline
double PointGetY(TPointType& ThisPoint)
{
    return ThisPoint.Y();
}

template<class TPointType> inline
double PointGetZ(TPointType& ThisPoint)
{
    return ThisPoint.Z();
}

void  AddPointsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<Point, Point::Pointer, array_1d<double,3> >(m,"Point") //WARNING: this was previously called Point3D
    .def(py::init<double, double, double>())
    .def(py::init<vector_expression<Vector> >())
    .def_property("X", PointGetX<Point >, PointSetX<Point >)
    .def_property("Y", PointGetY<Point >, PointSetY<Point >)
    .def_property("Z", PointGetZ<Point >, PointSetZ<Point >)
    ;
}

}  // namespace Python.

} // Namespace Kratos

