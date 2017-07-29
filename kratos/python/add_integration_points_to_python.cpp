//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//


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

