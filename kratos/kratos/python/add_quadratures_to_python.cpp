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


// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// System includes

// Project includes
#include "integration/line_gauss_legendre_integration_points.h"

#include "python/container_from_python.h"

namespace Kratos
{

template<class TStream, std::size_t TDimension>
TStream& operator << (TStream& rOStream,
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
    AddIntegrationPointsArray("IntegrationPointsArray", Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints1,1,Kratos::IntegrationPoint<3> >::IntegrationPointsArrayType());
    //AddIntegrationPointsArray("IntegrationPoint2DsArray", Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints<1>,2,Kratos::IntegrationPoint<3> >::IntegrationPointsArrayType());
    //AddIntegrationPointsArray("IntegrationPoint3DsArray", Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints<1>,3,Kratos::IntegrationPoint<3> >::IntegrationPointsArrayType());
    //AddIntegrationPointsArray("IntegrationPoint4DsArray", Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints<1>,4,Kratos::IntegrationPoint<3> >::IntegrationPointsArrayType());

    scope().attr("GaussLegendreQuadrature1D1") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints1,1,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature1D2") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints2,1,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature1D3") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints3,1,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature2D1") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints1,2,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature2D2") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints2,2,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature2D3") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints3,2,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature3D1") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints1,3,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature3D2") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints2,3,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    scope().attr("GaussLegendreQuadrature3D3") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints3,3,Kratos::IntegrationPoint<3> >::IntegrationPoints();
}

}  // namespace Python.

} // Namespace Kratos

