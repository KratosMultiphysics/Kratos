// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.2 $
//
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

