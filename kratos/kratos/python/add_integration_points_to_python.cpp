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



// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "integration/integration_point.h"
//#include "quadratures/integration_point.h"
#include "python/add_geometries_to_python.h"
#include "python/bounded_vector_python_interface.h"
#include "python/vector_scalar_operator_python.h"
#include "python/vector_vector_operator_python.h"

namespace Kratos
{

namespace Python
{
template<class TIntegrationPointType> inline
void IntegrationPointSetWeight(TIntegrationPointType& ThisIntegrationPoint, double Value)
{
    ThisIntegrationPoint.SetWeight(Value);
}

template<class TIntegrationPointType> inline
double IntegrationPointGetWeight(TIntegrationPointType& ThisIntegrationPoint)
{
    return ThisIntegrationPoint.Weight();
}

void  AddIntegrationPointsToPython()
{
    BoundedVectorPythonInterface<IntegrationPoint<1>, 1>::CreateInterfaceWithBase("IntegrationPoint1D", Point<1>())
    .def(init<double>())
    .def(init<Point<1>, double>())
    .def(init<vector_expression<vector<double> > >())
    .def(init<vector_expression<vector<double> >, double >())
    .add_property("Weight", IntegrationPointGetWeight<IntegrationPoint<1> >, IntegrationPointSetWeight<IntegrationPoint<1> >)
    .def(VectorScalarOperatorPython<IntegrationPoint<1>, double, IntegrationPoint<1> >())
    .def(VectorVectorOperatorPython<IntegrationPoint<1>, zero_vector<double>, IntegrationPoint<1> >())
    .def(VectorVectorOperatorPython<IntegrationPoint<1>, unit_vector<double>, IntegrationPoint<1> >())
    .def(VectorVectorOperatorPython<IntegrationPoint<1>, scalar_vector<double>, IntegrationPoint<1> >())
    .def(VectorVectorOperatorPython<IntegrationPoint<1>, vector<double>, IntegrationPoint<1> >())
    //.def(VectorVectorOperatorPython<IntegrationPoint<1>, mapped_vector<double>, IntegrationPoint<1> >())
    //.def(VectorVectorOperatorPython<IntegrationPoint<1>, compressed_vector<double>, IntegrationPoint<1> >())
    //.def(VectorVectorOperatorPython<IntegrationPoint<1>, coordinate_vector<double>, IntegrationPoint<1> >())
    ;

    BoundedVectorPythonInterface<IntegrationPoint<2>, 2>::CreateInterfaceWithBase("IntegrationPoint2D", Point<2>())
    .def(init<double>())
    .def(init<double, double>())
    .def(init<Point<2>, double>())
    .def(init<vector_expression<vector<double> > >())
    .def(init<vector_expression<vector<double> >, double >())
    .add_property("Weight", IntegrationPointGetWeight<IntegrationPoint<2> >, IntegrationPointSetWeight<IntegrationPoint<2> >)
    .def(VectorScalarOperatorPython<IntegrationPoint<2>, double, IntegrationPoint<2> >())
    .def(VectorVectorOperatorPython<IntegrationPoint<2>, zero_vector<double>, IntegrationPoint<2> >())
    .def(VectorVectorOperatorPython<IntegrationPoint<2>, unit_vector<double>, IntegrationPoint<2> >())
    .def(VectorVectorOperatorPython<IntegrationPoint<2>, scalar_vector<double>, IntegrationPoint<2> >())
    .def(VectorVectorOperatorPython<IntegrationPoint<2>, vector<double>, IntegrationPoint<2> >())
    //.def(VectorVectorOperatorPython<IntegrationPoint<2>, mapped_vector<double>, IntegrationPoint<2> >())
    //.def(VectorVectorOperatorPython<IntegrationPoint<2>, compressed_vector<double>, IntegrationPoint<2> >())
    //.def(VectorVectorOperatorPython<IntegrationPoint<2>, coordinate_vector<double>, IntegrationPoint<2> >())
    ;

    BoundedVectorPythonInterface<IntegrationPoint<3>, 3>::CreateInterfaceWithBase("IntegrationPoint3D", Point<3>())
    .def(init<double>())
    .def(init<double, double>())
    .def(init<double, double, double>())
    .def(init<Point<3>, double>())
    .def(init<vector_expression<vector<double> > >())
    .def(init<vector_expression<vector<double> >, double >())
    .add_property("Weight", IntegrationPointGetWeight<IntegrationPoint<3> >, IntegrationPointSetWeight<IntegrationPoint<3> >)
    .def(VectorScalarOperatorPython<IntegrationPoint<3>, double, IntegrationPoint<3> >())
    .def(VectorVectorOperatorPython<IntegrationPoint<3>, zero_vector<double>, IntegrationPoint<3> >())
    .def(VectorVectorOperatorPython<IntegrationPoint<3>, unit_vector<double>, IntegrationPoint<3> >())
    .def(VectorVectorOperatorPython<IntegrationPoint<3>, scalar_vector<double>, IntegrationPoint<3> >())
    .def(VectorVectorOperatorPython<IntegrationPoint<3>, vector<double>, IntegrationPoint<3> >())
    //.def(VectorVectorOperatorPython<IntegrationPoint<3>, mapped_vector<double>, IntegrationPoint<3> >())
    //.def(VectorVectorOperatorPython<IntegrationPoint<3>, compressed_vector<double>, IntegrationPoint<3> >())
    //.def(VectorVectorOperatorPython<IntegrationPoint<3>, coordinate_vector<double>, IntegrationPoint<3> >())
    ;
}

}  // namespace Python.

} // Namespace Kratos

