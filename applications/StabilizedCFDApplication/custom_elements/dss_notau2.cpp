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

#include "dss_notau2.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
DSS_notau2<TDim>::DSS_notau2(IndexType NewId):
    DSS<TDim>(NewId)
{}

template< unsigned int TDim >
DSS_notau2<TDim>::DSS_notau2(IndexType NewId, const NodesArrayType& ThisNodes):
    DSS<TDim>(NewId,ThisNodes)
{}


template< unsigned int TDim >
DSS_notau2<TDim>::DSS_notau2(IndexType NewId, GeometryType::Pointer pGeometry):
    DSS<TDim>(NewId,pGeometry)
{}


template< unsigned int TDim >
DSS_notau2<TDim>::DSS_notau2(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    DSS<TDim>(NewId,pGeometry,pProperties)
{}


template< unsigned int TDim >
DSS_notau2<TDim>::~DSS_notau2()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer DSS_notau2<TDim>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DSS_notau2(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output


template< unsigned int TDim >
std::string DSS_notau2<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "DSS_notau2 #" << this->Id();
    return buffer.str();
}


template< unsigned int TDim >
void DSS_notau2<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "DSS_notau2" << TDim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS_notau2<TDim>::CalculateStaticTau(double Density,
                                          double KinematicVisc,
                                          const array_1d<double,3> &Velocity,
                                          double ElemSize,
                                          const ProcessInfo& rProcessInfo,
                                          double &TauOne,
                                          double &TauTwo)
{
    const double c1 = 8.0;
    const double c2 = 2.0;

    double VelNorm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < TDim; d++)
        VelNorm += Velocity[d]*Velocity[d];
    VelNorm = std::sqrt(VelNorm);

    double InvTau = Density * ( c1 * KinematicVisc / (ElemSize*ElemSize) + c2 * VelNorm / ElemSize );
    TauOne = 1.0/InvTau;
    TauTwo = 0.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< unsigned int TDim >
void DSS_notau2<TDim>::save(Serializer& rSerializer) const
{
    typedef DSS<TDim> _basetype;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, _basetype );
}


template< unsigned int TDim >
void DSS_notau2<TDim>::load(Serializer& rSerializer)
{
    typedef DSS<TDim> _basetype;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, _basetype);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class DSS_notau2<2>;
template class DSS_notau2<3>;

}
