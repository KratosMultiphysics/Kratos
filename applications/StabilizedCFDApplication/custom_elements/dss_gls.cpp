//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//



#include "dss_gls.h"
#include "custom_utilities/turbulence_statistics_container.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
DSS_GLS<TDim>::DSS_GLS(IndexType NewId):
    DSS_FIC<TDim>(NewId)

{}

template< unsigned int TDim >
DSS_GLS<TDim>::DSS_GLS(IndexType NewId, const NodesArrayType& ThisNodes):
    DSS_FIC<TDim>(NewId,ThisNodes)
{}


template< unsigned int TDim >
DSS_GLS<TDim>::DSS_GLS(IndexType NewId, GeometryType::Pointer pGeometry):
    DSS_FIC<TDim>(NewId,pGeometry)
{}


template< unsigned int TDim >
DSS_GLS<TDim>::DSS_GLS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    DSS_FIC<TDim>(NewId,pGeometry,pProperties)
{}


template< unsigned int TDim >
DSS_GLS<TDim>::~DSS_GLS()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer DSS_GLS<TDim>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DSS_GLS(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS_GLS<TDim>::CalculateStabilizationParameters(double Density,
                                       double KinematicVisc,
                                       const array_1d<double,3> &Velocity,
                                       const ProcessInfo& rProcessInfo,
                                       double &TauIncompr,
                                       double &TauMomentum,
                                       array_1d<double, 3> &TauGrad)
{
    const double c1 = 8.0;
    const double c2 = 2.0;

    double Havg = this->AverageElementSize();

    double VelNorm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < TDim; d++)
        VelNorm += Velocity[d]*Velocity[d];
    VelNorm = std::sqrt(VelNorm);

    double InvTau = Density * ( c1 * KinematicVisc / (Havg*Havg) + c2 * VelNorm / Havg );
    TauIncompr = 1.0/InvTau;
    TauMomentum = TauIncompr;
}



///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< unsigned int TDim >
void DSS_GLS<TDim>::save(Serializer& rSerializer) const
{
    typedef DSS_FIC<TDim> _BaseT;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, _BaseT );
}


template< unsigned int TDim >
void DSS_GLS<TDim>::load(Serializer& rSerializer)
{
    typedef DSS_FIC<TDim> _BaseT;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, _BaseT );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class DSS_GLS<2>;
template class DSS_GLS<3>;

}
