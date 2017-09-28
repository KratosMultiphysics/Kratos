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
    /*
    BoundedVectorPythonInterface<Point<1>, 1>::CreateInterface("Point1D")
      .def(init<double>())
      .def(init<vector_expression<vector<double> > >())
      .add_property("X", PointGetX<Point<1> >, PointSetX<Point<1> >)
      .def(VectorScalarOperatorPython<Point<1>, double, Point<1> >())
      .def(VectorVectorOperatorPython<Point<1>, zero_vector<double>, Point<1> >())
      .def(VectorVectorOperatorPython<Point<1>, unit_vector<double>, Point<1> >())
      .def(VectorVectorOperatorPython<Point<1>, scalar_vector<double>, Point<1> >())
      .def(VectorVectorOperatorPython<Point<1>, vector<double>, Point<1> >())
      .def(VectorVectorOperatorPython<Point<1>, mapped_vector<double>, Point<1> >())
      .def(VectorVectorOperatorPython<Point<1>, compressed_vector<double>, Point<1> >())
      .def(VectorVectorOperatorPython<Point<1>, coordinate_vector<double>, Point<1> >())
       ;

     BoundedVectorPythonInterface<Point<2>, 2>::CreateInterface("Point2D")
       .def(init<double>())
       .def(init<double, double>())
       .def(init<vector_expression<vector<double> > >())
       .add_property("X", PointGetX<Point<2> >, PointSetX<Point<2> >)
       .add_property("Y", PointGetY<Point<2> >, PointSetY<Point<2> >)
       .def(VectorScalarOperatorPython<Point<2>, double, Point<2> >())
       .def(VectorVectorOperatorPython<Point<2>, zero_vector<double>, Point<2> >())
       .def(VectorVectorOperatorPython<Point<2>, unit_vector<double>, Point<2> >())
       .def(VectorVectorOperatorPython<Point<2>, scalar_vector<double>, Point<2> >())
       .def(VectorVectorOperatorPython<Point<2>, vector<double>, Point<2> >())
       .def(VectorVectorOperatorPython<Point<2>, mapped_vector<double>, Point<2> >())
       .def(VectorVectorOperatorPython<Point<2>, compressed_vector<double>, Point<2> >())
       .def(VectorVectorOperatorPython<Point<2>, coordinate_vector<double>, Point<2> >())
        ;
    */
    BoundedVectorPythonInterface<Point<3>, 3>::CreateInterface("Point3D")
    .def(init<double>())
    .def(init<double, double>())
    .def(init<double, double, double>())
    .def(init<vector_expression<vector<double> > >())
    .add_property("X", PointGetX<Point<3> >, PointSetX<Point<3> >)
    .add_property("Y", PointGetY<Point<3> >, PointSetY<Point<3> >)
    .add_property("Z", PointGetZ<Point<3> >, PointSetZ<Point<3> >)
    .def(VectorScalarOperatorPython<Point<3>, double, Point<3> >())
    .def(VectorVectorOperatorPython<Point<3>, zero_vector<double>, Point<3> >())
    .def(VectorVectorOperatorPython<Point<3>, unit_vector<double>, Point<3> >())
    .def(VectorVectorOperatorPython<Point<3>, scalar_vector<double>, Point<3> >())
    .def(VectorVectorOperatorPython<Point<3>, vector<double>, Point<3> >())
    //.def(VectorVectorOperatorPython<Point<3>, mapped_vector<double>, Point<3> >())
    //.def(VectorVectorOperatorPython<Point<3>, compressed_vector<double>, Point<3> >())
    //       .def(VectorVectorOperatorPython<Point<3>, coordinate_vector<double>, Point<3> >())
    ;
}

}  // namespace Python.

} // Namespace Kratos

