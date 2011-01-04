/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
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

