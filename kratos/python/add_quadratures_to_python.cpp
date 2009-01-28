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
//   Last modified by:    $Author: janosch $
//   Date:                $Date: 2007-12-13 14:50:05 $
//   Revision:            $Revision: 1.4 $
//
//


// System includes 
#include <cstddef>

// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


// Project includes
#include "includes/define.h"
#include "integration/quadrature.h" // to bo changed to the new version below.
#include "integration/gauss_legendre_integration_points.h"
//#include "quadratures/gauss_legendre_integration_points.h"
//#include "quadratures/quadrature.h"
#include "python/container_from_python.h"

namespace Kratos
{

	template<std::size_t TDimension>
  std::stringstream& operator << (std::stringstream& rOStream,
				  std::vector<IntegrationPoint<TDimension> > const & rThis)
    { 
 	typename std::vector<IntegrationPoint<TDimension> >::size_type i;
 	for(i = 0 ; i < rThis.size() - 1 ; i++)
	  rOStream << rThis[i] << " , " << std::endl;
 	rOStream << rThis[i];
      return rOStream;
    }

namespace Python
{

  using namespace boost::python;

  template<class TArrayType>
  void AddIntegrationPointsArray(std::string const& rArrayName, TArrayType const & Dummy)
  {
    ContainerFromPython< TArrayType >();

     class_<TArrayType>(rArrayName.c_str(), init<int>())
	.def(init<const TArrayType&>())
	.def(vector_indexing_suite<TArrayType>())
 	.def(self_ns::str(self))
	;

  }

  void  AddQuadraturesToPython()
  {
    AddIntegrationPointsArray("IntegrationPointsArray", Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints1,1,Kratos::IntegrationPoint<3> >::IntegrationPointsArrayType());
    //AddIntegrationPointsArray("IntegrationPoint2DsArray", Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints<1>,2,Kratos::IntegrationPoint<3> >::IntegrationPointsArrayType());
    //AddIntegrationPointsArray("IntegrationPoint3DsArray", Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints<1>,3,Kratos::IntegrationPoint<3> >::IntegrationPointsArrayType());
    //AddIntegrationPointsArray("IntegrationPoint4DsArray", Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints<1>,4,Kratos::IntegrationPoint<3> >::IntegrationPointsArrayType());

    scope().attr("GaussLegendreQuadrature1D1") = Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints1,1,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature1D2") = Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints2,1,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature1D3") = Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints3,1,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature2D1") = Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints1,2,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature2D2") = Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints2,2,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature2D3") = Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints3,2,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature3D1") = Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints1,3,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature3D2") = Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints2,3,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature3D3") = Kratos::Quadrature<Kratos::GaussLegendreIntegrationPoints3,3,Kratos::IntegrationPoint<3> >::IntegrationPoints();
  }
	
}  // namespace Python.

} // Namespace Kratos

