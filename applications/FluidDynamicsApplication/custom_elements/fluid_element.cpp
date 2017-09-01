//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include "fluid_element.h"
#include "fluid_element_data.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
FluidElement<TElementData>::FluidElement(IndexType NewId):
    Element(NewId)

{}

template< class TElementData >
FluidElement<TElementData>::FluidElement(IndexType NewId, const NodesArrayType& ThisNodes):
    Element(NewId,ThisNodes)
{}


template< class TElementData >
FluidElement<TElementData>::FluidElement(IndexType NewId, GeometryType::Pointer pGeometry):
    Element(NewId,pGeometry)
{}


template< class TElementData >
FluidElement<TElementData>::FluidElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    Element(NewId,pGeometry,pProperties)
{}


template< class TElementData >
FluidElement<TElementData>::~FluidElement()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer FluidElement<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR << "Attempting to Create base FluidElement<" << TElementData::Dim << "> instances, but this is an abstract element." << std::endl;
}

template< class TElementData >
Element::Pointer FluidElement<TElementData>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR << "Attempting to Create base FluidElement<" << TElementData::Dim << "> instances, but this is an abstract element." << std::endl;
}

template <class TElementData>
void FluidElement<TElementData>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                      VectorType& rRightHandSideVector,
                                                      ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if( rLeftHandSideMatrix.size1() != TElementData::LocalSize )
        rLeftHandSideMatrix.resize(TElementData::LocalSize,TElementData::LocalSize,false);

    if( rRightHandSideVector.size() != TElementData::LocalSize )
        rRightHandSideVector.resize(TElementData::LocalSize,false);

    rLeftHandSideMatrix = ZeroMatrix(TElementData::LocalSize,TElementData::LocalSize);
    rRightHandSideVector = ZeroVector(TElementData::LocalSize);
}

template <class TElementData>
void FluidElement<TElementData>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                       ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if( rLeftHandSideMatrix.size1() != TElementData::LocalSize )
        rLeftHandSideMatrix.resize(TElementData::LocalSize,TElementData::LocalSize,false);

    rLeftHandSideMatrix = ZeroMatrix(TElementData::LocalSize,TElementData::LocalSize);
}

template <class TElementData>
void FluidElement<TElementData>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                        ProcessInfo& rCurrentProcessInfo)
{
    if( rRightHandSideVector.size() != TElementData::LocalSize )
        rRightHandSideVector.resize(TElementData::LocalSize,false);

    rRightHandSideVector = ZeroVector(TElementData::LocalSize);
}

template <class TElementData>
void FluidElement<TElementData>::CalculateLocalVelocityContribution(
    MatrixType& rDampMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if( rDampMatrix.size1() != TElementData::LocalSize )
        rDampMatrix.resize(TElementData::LocalSize,TElementData::LocalSize,false);

    if( rRightHandSideVector.size() != TElementData::LocalSize )
        rRightHandSideVector.resize(TElementData::LocalSize,false);

    rDampMatrix = ZeroMatrix(TElementData::LocalSize,TElementData::LocalSize);
    rRightHandSideVector = ZeroVector(TElementData::LocalSize);

    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    TElementData Data(this->GetGeometry());
    IntegrationPointData<TElementData> IPData;

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        IntegrationPointData<TElementData>::FillIntegrationPointData(
            IPData, Data, g, GaussWeights, ShapeFunctions, ShapeDerivatives);

        this->AddSystemTerms(Data,IPData,rCurrentProcessInfo,rDampMatrix,rRightHandSideVector);
    }

    // Rewrite local contribution into residual form (A*dx = b - A*x)
    VectorType U = ZeroVector(TElementData::LocalSize);
    int LocalIndex = 0;

    for (unsigned int i = 0; i < TElementData::NumNodes; ++i)
    {
        for (unsigned int d = 0; d < TElementData::Dim; ++d) // Velocity Dofs
            U[LocalIndex++] = Data.Velocity(i,d);
        U[LocalIndex++] = Data.Pressure[i]; // Pressure Dof
    }

    noalias(rRightHandSideVector) -= prod(rDampMatrix, U);
}

template <class TElementData>
void FluidElement<TElementData>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                     ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if( rMassMatrix.size1() != TElementData::LocalSize )
        rMassMatrix.resize(TElementData::LocalSize,TElementData::LocalSize,false);

    rMassMatrix = ZeroMatrix(TElementData::LocalSize,TElementData::LocalSize);

    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    TElementData Data(this->GetGeometry());
    IntegrationPointData<TElementData> IPData;

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        IntegrationPointData<TElementData>::FillIntegrationPointData(
            IPData, Data, g, GaussWeights, ShapeFunctions, ShapeDerivatives);

        this->AddMassTerms(Data,IPData,rCurrentProcessInfo, rMassMatrix);
    }
}


template< class TElementData >
void FluidElement<TElementData>::Calculate(const Variable<double> &rVariable,
                          double &rOutput,
                          const ProcessInfo &rCurrentProcessInfo)
{

}


template< class TElementData >
void FluidElement<TElementData>::Calculate(const Variable<array_1d<double,3> > &rVariable,
                          array_1d<double,3> &rOutput,
                          const ProcessInfo &rCurrentProcessInfo)
{
    // Lumped projection terms
    if (rVariable == ADVPROJ)
    {
        this->CalculateProjections();
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// For TElementData::Dim == 2
///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
void FluidElement< FluidElementData<2,3> >::EquationIdVector(EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    unsigned int LocalIndex = 0;

    if (rResult.size() != FluidElementData<2,3>::LocalSize)
        rResult.resize(FluidElementData<2,3>::LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (unsigned int i = 0; i < FluidElementData<2,3>::NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE,ppos).EquationId();
    }
}


template<>
void FluidElement< FluidElementData<2,3> >::GetDofList(DofsVectorType &rElementalDofList, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();

     if (rElementalDofList.size() != FluidElementData<2,3>::LocalSize)
         rElementalDofList.resize(FluidElementData<2,3>::LocalSize);

     unsigned int LocalIndex = 0;

     for (unsigned int i = 0; i < FluidElementData<2,3>::NumNodes; ++i)
     {
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
     }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// For TElementData::Dim == 3
///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
void FluidElement< FluidElementData<3,4> >::EquationIdVector(EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();

    unsigned int LocalIndex = 0;

    if (rResult.size() != FluidElementData<3,4>::LocalSize)
        rResult.resize(FluidElementData<3,4>::LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (unsigned int i = 0; i < FluidElementData<3,4>::NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE,ppos).EquationId();
    }
}


template<>
void FluidElement< FluidElementData<3,4> >::GetDofList(DofsVectorType &rElementalDofList, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();

     if (rElementalDofList.size() != FluidElementData<3,4>::LocalSize)
         rElementalDofList.resize(FluidElementData<3,4>::LocalSize);

     unsigned int LocalIndex = 0;

     for (unsigned int i = 0; i < FluidElementData<3,4>::NumNodes; ++i)
     {
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
     }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::GetFirstDerivativesVector(Vector &rValues, int Step)
{
    GeometryType& rGeom = this->GetGeometry();

    if (rValues.size() != TElementData::LocalSize)
        rValues.resize(TElementData::LocalSize,false);

    noalias(rValues) = ZeroVector(TElementData::LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        const array_1d<double,3>& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY,Step);
        for (unsigned int d = 0; d < TElementData::Dim; d++)
            rValues[Index++] = rVel[d];
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(PRESSURE,Step);
    }
}


template< class TElementData >
void FluidElement<TElementData>::GetSecondDerivativesVector(Vector &rValues, int Step)
{
    GeometryType& rGeom = this->GetGeometry();

    if (rValues.size() != TElementData::LocalSize)
        rValues.resize(TElementData::LocalSize,false);

    noalias(rValues) = ZeroVector(TElementData::LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        const array_1d<double,3>& rAcc = rGeom[i].FastGetSolutionStepValue(ACCELERATION,Step);
        for (unsigned int d = 0; d < TElementData::Dim; d++)
            rValues[Index++] = rAcc[d];
        rValues[Index++] = 0.0; // skip pressure Dof
    }
}


template< class TElementData >
GeometryData::IntegrationMethod FluidElement<TElementData>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int FluidElement<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);

    // Check variables used by TElementData
    out = TElementData::Check(*this);

    // Extra variables used in computing projections
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(NODAL_AREA);

    // Elemental data
    KRATOS_CHECK_VARIABLE_KEY(C_SMAGORINSKY);

    for(unsigned int i=0; i<TElementData::NumNodes; ++i)
    {
        Node<3>& rNode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_AREA,rNode);

        // Check that required dofs exist
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X,rNode);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y,rNode);
        if (TElementData::Dim == 3) KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z,rNode);
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE,rNode);
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if ( TElementData::Dim == 2)
    {
        for (unsigned int i=0; i<TElementData::NumNodes; ++i)
        {
            if (this->GetGeometry()[i].Z() != 0.0)
                KRATOS_ERROR << "Node " << this->GetGeometry()[i].Id() << "has non-zero Z coordinate." << std::endl;
        }
    }

    return out;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output


template< class TElementData >
std::string FluidElement<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "FluidElement #" << Id();
    return buffer.str();
}


template< class TElementData >
void FluidElement<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "FluidElement" << TElementData::Dim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::CalculateGeometryData(Vector &rGaussWeights,
                                      Matrix &rNContainer,
                                      ShapeFunctionDerivativesArrayType &rDN_DX)
{
    const GeometryData::IntegrationMethod IntMethod = this->GetIntegrationMethod();
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumGauss = rGeom.IntegrationPointsNumber(IntMethod);

    Vector DetJ;
    rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX,DetJ,IntMethod);

    rNContainer.resize(NumGauss,TElementData::NumNodes,false);
    rNContainer = rGeom.ShapeFunctionsValues(IntMethod);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(IntMethod);

    rGaussWeights.resize(NumGauss,false);

    for (unsigned int g = 0; g < NumGauss; g++)
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
}

/** Calculate characteristic element length.
 * @return Minimum element height
 */
template< class TElementData >
double FluidElement<TElementData>::ElementSize()
{
    KRATOS_TRY;
    GeometryType& rGeom = this->GetGeometry();
    GeometryData::KratosGeometryFamily GeoFamily = rGeom.GetGeometryFamily();

    switch (GeoFamily)
    {
    case GeometryData::Kratos_Triangle:
    {
        /* Calculate node-edge distances */
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();

        double x20 = rGeom[2].X() - rGeom[0].X();
        double y20 = rGeom[2].Y() - rGeom[0].Y();

        // node 0, edge 12
        double nx = -(y20-y10);
        double ny = x20-x10;
        double Hsq = x10*nx + y10*ny;
        Hsq *= Hsq / (nx*nx + ny*ny);

        // node 1, edge 20
        nx = -y20;
        ny = x20;
        double hsq = x10*nx + y10*ny;
        hsq *= hsq / (nx*nx + ny*ny);
        Hsq = ( hsq < Hsq ) ? hsq : Hsq;

        // node 2, edge 10
        nx = -y10;
        ny = x10;
        hsq = x20*nx + y20*ny;
        hsq *= hsq / (nx*nx + ny*ny);
        Hsq = ( hsq < Hsq ) ? hsq : Hsq;
        return std::sqrt(Hsq);
    }
    case GeometryData::Kratos_Quadrilateral:
    {
        /* Calculate node-edge distances, assuming parallel faces */
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();

        double x20 = rGeom[2].X() - rGeom[0].X();
        double y20 = rGeom[2].Y() - rGeom[0].Y();

        double x30 = rGeom[3].X() - rGeom[0].X();
        double y30 = rGeom[3].Y() - rGeom[0].Y();

        // node 0, edge 12
        double nx = -(y20-y10);
        double ny = x20-y10;
        double Hsq = x10*nx + y10*ny;
        Hsq *= Hsq / std::sqrt(nx*nx + ny*ny);

        // node 0, edge 23
        nx = -(y30-y20);
        ny = x30-y20;
        double hsq = x20*nx + y20*ny;
        hsq *= hsq / (nx*nx + ny*ny);
        Hsq = ( hsq < Hsq ) ? hsq : Hsq;
        return std::sqrt(Hsq);
    }
    case GeometryData::Kratos_Tetrahedra:
    {
        /* Calculate distances between each node and the opposite face */
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();
        double z10 = rGeom[1].Z() - rGeom[0].Z();

        double x20 = rGeom[2].X() - rGeom[0].X();
        double y20 = rGeom[2].Y() - rGeom[0].Y();
        double z20 = rGeom[2].Z() - rGeom[0].Z();

        double x30 = rGeom[3].X() - rGeom[0].X();
        double y30 = rGeom[3].Y() - rGeom[0].Y();
        double z30 = rGeom[3].Z() - rGeom[0].Z();

        // face 123
        double nx = (y30-y10)*(z20-z10) - (z30-z10)*(y20-y10);
        double ny = (z30-z10)*(x20-x10) - (x30-x10)*(z20-z10);
        double nz = (x30-x10)*(y20-y10) - (y30-y10)*(x20-x10);
        double Hsq = x10*nx + y10*ny + z10*nz; // scalar product x10*n
        Hsq *= Hsq / (nx*nx + ny*ny + nz*nz); // H^2 = (x10*n)^2 / ||n||^2

        // face 230
        nx = y30*z20 - z30*y20;
        ny = z30*x20 - x30*z20;
        nz = x30*y20 - y30*x20;
        double hsq = x10*nx + y10*ny + z10*nz;
        hsq *= hsq / (nx*nx + ny*ny + nz*nz);
        Hsq = (hsq < Hsq) ? hsq : Hsq;

        // face 301
        nx = y10*z30 - z10*y30;
        ny = z10*x30 - x10*z30;
        nz = x10*y30 - y10*x30;
        hsq = x20*nx + y20*ny + z20*nz;
        hsq *= hsq / (nx*nx + ny*ny + nz*nz);
        Hsq = (hsq < Hsq) ? hsq : Hsq;

        // face 012
        nx = y10*z20 - z10*y20;
        ny = z10*x20 - x10*z20;
        nz = x10*y20 - y10*x20;
        hsq = x30*nx + y30*ny + z30*nz;
        hsq *= hsq / (nx*nx + ny*ny + nz*nz);
        Hsq = (hsq < Hsq) ? hsq : Hsq;
        return std::sqrt(Hsq);
    }
    case GeometryData::Kratos_Hexahedra:
    {
        /* Numbering assumes bottom face nodes 0123, top nodes 4567
         * considering the distance between a few face--node pairs:
         * face node
         *  034  1
         *  014  3
         *  013  4
         * This assumes parallel faces!
         * (otherwise all nodes should be checked against opposite faces)
         */
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();
        double z10 = rGeom[1].Z() - rGeom[0].Z();

        double x30 = rGeom[3].X() - rGeom[0].X();
        double y30 = rGeom[3].Y() - rGeom[0].Y();
        double z30 = rGeom[3].Z() - rGeom[0].Z();

        double x40 = rGeom[4].X() - rGeom[0].X();
        double y40 = rGeom[4].Y() - rGeom[0].Y();
        double z40 = rGeom[4].Z() - rGeom[0].Z();


        // Face 034
        double nx = y30*z40 - z30*y40;
        double ny = z30*x40 - x30*z40;
        double nz = x30*y40 - y30*x40;
        double Hsq = x10*nx + y10*ny + z10*nz; // scalar product x10*n
        Hsq *= Hsq / (nx*nx + ny*ny + nz*nz); // H^2 = (x10*n)^2 / ||n||^2

        // face 014
        nx = y10*z40 - z10*y40;
        ny = z10*x40 - x10*z40;
        nz = x10*y40 - y10*x40;
        double hsq = x30*nx + y30*ny + z30*nz;
        hsq *= hsq / (nx*nx + ny*ny + nz*nz);
        Hsq = (hsq < Hsq) ? hsq : Hsq;

        // face 013
        nx = y10*z30 - z10*y30;
        ny = z10*x30 - x10*z30;
        nz = x10*y30 - y10*x30;
        hsq = x40*nx + y40*ny + z40*nz;
        hsq *= hsq / (nx*nx + ny*ny + nz*nz);
        Hsq = (hsq < Hsq) ? hsq : Hsq;
        return std::sqrt(Hsq);
    }
    default:
    {
        KRATOS_THROW_ERROR(std::invalid_argument,"FluidElement::ElementSize not implemented for this geometry type","");
        return 0;
    }
    }
    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::CalculateStaticTau(double Density,
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
    for (unsigned int d = 1; d < TElementData::Dim; d++)
        VelNorm += Velocity[d]*Velocity[d];
    VelNorm = std::sqrt(VelNorm);

    double InvTau = Density * ( c1 * KinematicVisc / (ElemSize*ElemSize) + c2 * VelNorm / ElemSize );
    TauOne = 1.0/InvTau;
    TauTwo = Density * (KinematicVisc + c2 * VelNorm * ElemSize / c1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////


template< class TElementData >
double FluidElement<TElementData>::EffectiveViscosity(
    const TElementData& rData,
    const IntegrationPointData<TElementData>& rIPData,
    double ElemSize,
    const ProcessInfo &rCurrentProcessInfo)
{
    const FluidElement* const_this = static_cast<const FluidElement*>(this);
    double Csmag = const_this->GetValue(C_SMAGORINSKY);

    double KinViscosity = rIPData.Viscosity;

    if (Csmag != 0.0 )
    {
        // Calculate Symetric gradient
        MatrixType S = ZeroMatrix(TElementData::Dim,TElementData::Dim);
        for (unsigned int n = 0; n < TElementData::NumNodes; ++n)
        {
            for (unsigned int i = 0; i < TElementData::Dim; ++i)
                for (unsigned int j = 0; j < TElementData::Dim; ++j)
                    S(i,j) += 0.5 * ( rIPData.DN_DX(n,j) * rData.Velocity(n,i) + rIPData.DN_DX(n,i) * rData.Velocity(n,j) );
        }

        // Norm of symetric gradient
        double NormS = 0.0;
        for (unsigned int i = 0; i < TElementData::Dim; ++i)
            for (unsigned int j = 0; j < TElementData::Dim; ++j)
                NormS += S(i,j) * S(i,j);
        NormS = sqrt(2.0*NormS);

        // Nu_sgs = (Csmag * Delta)^2 * (2*Sij*Sij)^(1/2)
        KinViscosity += Csmag * Csmag * ElemSize * ElemSize * NormS;
    }

    return KinViscosity;
}


template< class TElementData >
void FluidElement<TElementData>::ResolvedConvectiveVelocity(
    const TElementData& rData,
    const ShapeFunctionsType &rN,
    array_1d<double,3> &rConvVel)
{
    for (unsigned int d = 0; d < TElementData::Dim; d++)
        rConvVel[d] = rN[0] * ( rData.Velocity(0,d) - rData.MeshVelocity(0,d) );

    for (unsigned int i = 1; i < TElementData::NumNodes; i++)
    {
        for (unsigned int d = 0; d < TElementData::Dim; d++)
            rConvVel[d] += rN[i] * ( rData.Velocity(i,d) - rData.MeshVelocity(i,d) );
    }
}


template< class TElementData >
void FluidElement<TElementData>::FullConvectiveVelocity(array_1d<double,3> &rConvVel,
                                       const ShapeFunctionsType &rN,
                                       const array_1d<double,3> &rSubscaleVel)
{
    GeometryType& rGeom = this->GetGeometry();

    rConvVel = rSubscaleVel;

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
        rConvVel += rN[i]*(rGeom[i].FastGetSolutionStepValue(VELOCITY)-rGeom[i].FastGetSolutionStepValue(MESH_VELOCITY));
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::ConvectionOperator(Vector &rResult,
                                   const array_1d<double,3> &rConvVel,
                                   const ShapeFunctionDerivativesType &DN_DX)
{
    if(rResult.size() != TElementData::NumNodes) rResult.resize(TElementData::NumNodes,false);

    for (unsigned int i = 0; i < TElementData::NumNodes; i++)
    {
        rResult[i] = rConvVel[0]*DN_DX(i,0);
        for(unsigned int k = 1; k < TElementData::Dim; k++)
            rResult[i] += rConvVel[k]*DN_DX(i,k);
    }
}


template<>
void FluidElement< FluidElementData<2,3> >::IntegrationPointVorticity(const ShapeFunctionDerivativesType &rDN_DX,array_1d<double,3>& rVorticity) const
{
    rVorticity = array_1d<double,3>(3,0.0);

    for (unsigned int i = 0; i < FluidElementData<2,3>::NumNodes; ++i)
    {
        const array_1d<double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        rVorticity[2] += rDN_DX(i,0)*rVelocity[1] - rDN_DX(i,1)*rVelocity[0];
    }
}


template<>
void FluidElement< FluidElementData<3,4> >::IntegrationPointVorticity(const ShapeFunctionDerivativesType &rDN_DX,array_1d<double,3>& rVorticity) const
{
    rVorticity = array_1d<double,3>(3,0.0);

    for (unsigned int i = 0; i < FluidElementData<3,4>::NumNodes; ++i)
    {
        const array_1d<double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        rVorticity[0] += rDN_DX(i,1)*rVelocity[2] - rDN_DX(i,2)*rVelocity[1];
        rVorticity[1] += rDN_DX(i,2)*rVelocity[0] - rDN_DX(i,0)*rVelocity[2];
        rVorticity[2] += rDN_DX(i,0)*rVelocity[1] - rDN_DX(i,1)*rVelocity[0];
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void FluidElement<TElementData>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
}


template< class TElementData >
void FluidElement<TElementData>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class FluidElement< FluidElementData<2,3> >;
template class FluidElement< FluidElementData<3,4> >;

}
