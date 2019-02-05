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


// External includes

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

namespace py = pybind11;

template<class TArrayType>
void AddIntegrationPointsArray(std::string ArrayName, TArrayType const & Dummy)
{
    ContainerInterface< TArrayType >::CreateInterfac(m, ArrayName)
    ContainerFromPython< TArrayType >();

    py::class_<TArrayType>(rArrayName.c_str(), init<int>())
    .def(py::init<const TArrayType&>())
    .def(vector_indexing_suite<TArrayType>())
    .def(self_ns::str(self))
    ;

}

void  AddQuadraturesToPython(pybind11::module& m)
{
    ContainerInterface<  Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints1,1,Kratos::IntegrationPoint<3> > >::CreateInterface(m, "IntegrationPointsArray");
    ContainerInterface<  Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints1,2,Kratos::IntegrationPoint<3> > >::CreateInterface(m, "IntegrationPoints2DArray");
    ContainerInterface<  Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints1,3,Kratos::IntegrationPoint<3> > >::CreateInterface(m, "IntegrationPoints3DArray");

//     AddIntegrationPointsArray("IntegrationPointsArray", Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints1,1,Kratos::IntegrationPoint<3> >::IntegrationPointsArrayType());
//     AddIntegrationPointsArray("IntegrationPoint2DsArray", Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints<1>,2,Kratos::IntegrationPoint<3> >::IntegrationPointsArrayType());
//     AddIntegrationPointsArray("IntegrationPoint3DsArray", Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints<1>,3,Kratos::IntegrationPoint<3> >::IntegrationPointsArrayType());
    //AddIntegrationPointsArray("IntegrationPoint4DsArray", Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints<1>,4,Kratos::IntegrationPoint<3> >::IntegrationPointsArrayType());

    m.attr("GaussLegendreQuadrature1D1") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints1,1,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    m.attr("GaussLegendreQuadrature1D2") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints2,1,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    m.attr("GaussLegendreQuadrature1D3") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints3,1,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    m.attr("GaussLegendreQuadrature2D1") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints1,2,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    m.attr("GaussLegendreQuadrature2D2") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints2,2,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    m.attr("GaussLegendreQuadrature2D3") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints3,2,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    m.attr("GaussLegendreQuadrature3D1") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints1,3,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    m.attr("GaussLegendreQuadrature3D2") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints2,3,Kratos::IntegrationPoint<3> >::IntegrationPoints();
    m.attr("GaussLegendreQuadrature3D3") = Kratos::Quadrature<Kratos::LineGaussLegendreIntegrationPoints3,3,Kratos::IntegrationPoint<3> >::IntegrationPoints();
}

}  // namespace Python.

} // Namespace Kratos

