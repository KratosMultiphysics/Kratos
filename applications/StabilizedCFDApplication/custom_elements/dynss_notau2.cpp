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

#include "dynss_notau2.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
DYNSS_NOTAU2<TDim>::DYNSS_NOTAU2(IndexType NewId):
    DynSS<TDim>(NewId)
{}

template< unsigned int TDim >
DYNSS_NOTAU2<TDim>::DYNSS_NOTAU2(IndexType NewId, const NodesArrayType& ThisNodes):
    DynSS<TDim>(NewId,ThisNodes)
{}


template< unsigned int TDim >
DYNSS_NOTAU2<TDim>::DYNSS_NOTAU2(IndexType NewId, GeometryType::Pointer pGeometry):
    DynSS<TDim>(NewId,pGeometry)
{}


template< unsigned int TDim >
DYNSS_NOTAU2<TDim>::DYNSS_NOTAU2(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    DynSS<TDim>(NewId,pGeometry,pProperties)
{}


template< unsigned int TDim >
DYNSS_NOTAU2<TDim>::~DYNSS_NOTAU2()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer DYNSS_NOTAU2<TDim>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DYNSS_NOTAU2(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output


template< unsigned int TDim >
std::string DYNSS_NOTAU2<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "DYNSS_NOTAU2 #" << this->Id();
    return buffer.str();
}


template< unsigned int TDim >
void DYNSS_NOTAU2<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "DYNSS_NOTAU2" << TDim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////



template< unsigned int TDim >
void DYNSS_NOTAU2<TDim>::CalculateTau(double Density,
                                      double KinematicVisc,
                                      const array_1d<double,3> &Velocity,
                                      const ProcessInfo& rProcessInfo,
                                      double ElemSize,
                                      double &TauOne,
                                      double &TauTwo,
                                      double &TauP)
{
    const double c1 = 8.0;
    const double c2 = 2.0;
    const double dt = rProcessInfo[DELTA_TIME];

    double VelNorm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < TDim; d++)
        VelNorm += Velocity[d]*Velocity[d];
    VelNorm = std::sqrt(VelNorm);

    double InvTau = Density * ( 1.0/dt + c1 * KinematicVisc / (ElemSize*ElemSize) + c2 * VelNorm / ElemSize );
    TauOne = 1.0/InvTau;
    //TauTwo = Density * (KinematicVisc + c2 * VelNorm * ElemSize / c1);
    TauTwo = 0.0;

    // Auxiliary coefficient StaticTauOne*TauTwo/Dt that appears on the pressure subscale model
    //TauP = Density * ElemSize*ElemSize / (c1*dt);
    TauP = 0.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< unsigned int TDim >
void DYNSS_NOTAU2<TDim>::save(Serializer& rSerializer) const
{
    typedef DynSS<TDim> _basetype;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, _basetype );
}


template< unsigned int TDim >
void DYNSS_NOTAU2<TDim>::load(Serializer& rSerializer)
{
    typedef DynSS<TDim> _basetype;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, _basetype);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class DYNSS_NOTAU2<2>;
template class DYNSS_NOTAU2<3>;

}
