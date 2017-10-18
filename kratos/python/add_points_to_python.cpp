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
#include "python/add_points_to_python.h"
#include "python/bounded_vector_python_interface.h"
#include "python/vector_scalar_operator_python.h"
#include "python/vector_vector_operator_python.h"

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

void  AddPointsToPython()
{
    BoundedVectorPythonInterface<Point, 3>::CreateInterface("Point3D")
    .def(init<double>())
    .def(init<double, double>())
    .def(init<double, double, double>())
    .def(init<vector_expression<vector<double> > >())
    .add_property("X", PointGetX<Point >, PointSetX<Point >)
    .add_property("Y", PointGetY<Point >, PointSetY<Point >)
    .add_property("Z", PointGetZ<Point >, PointSetZ<Point >)
    .def(VectorScalarOperatorPython<Point, double, Point >())
    .def(VectorVectorOperatorPython<Point, zero_vector<double>, Point >())
    .def(VectorVectorOperatorPython<Point, unit_vector<double>, Point >())
    .def(VectorVectorOperatorPython<Point, scalar_vector<double>, Point >())
    .def(VectorVectorOperatorPython<Point, vector<double>, Point >())
    ;
}

}  // namespace Python.

} // Namespace Kratos

