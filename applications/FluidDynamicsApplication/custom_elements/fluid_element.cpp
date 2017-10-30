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
    KRATOS_ERROR << "Attempting to Create base FluidElement instances, but this is an abstract element." << std::endl;
}

template< class TElementData >
Element::Pointer FluidElement<TElementData>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR << "Attempting to Create base FluidElement instances, but this is an abstract element." << std::endl;
}

template <class TElementData>
void FluidElement<TElementData>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                      VectorType& rRightHandSideVector,
                                                      ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if( rLeftHandSideMatrix.size1() != LocalSize )
        rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);

    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize,false);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);
    rRightHandSideVector = ZeroVector(LocalSize);
}

template <class TElementData>
void FluidElement<TElementData>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                       ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if( rLeftHandSideMatrix.size1() != LocalSize )
        rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);
}

template <class TElementData>
void FluidElement<TElementData>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                        ProcessInfo& rCurrentProcessInfo)
{
    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize,false);

    rRightHandSideVector = ZeroVector(LocalSize);
}

template <class TElementData>
void FluidElement<TElementData>::CalculateLocalVelocityContribution(
    MatrixType& rDampMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if( rDampMatrix.size1() != LocalSize )
        rDampMatrix.resize(LocalSize,LocalSize,false);

    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize,false);

    rDampMatrix = ZeroMatrix(LocalSize,LocalSize);
    rRightHandSideVector = ZeroVector(LocalSize);

    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    TElementData data;
    data.Initialize(*this,rCurrentProcessInfo);

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        IntegrationPointGeometryData integration_point(
            GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g]);

        this->AddSystemTerms(data,integration_point,rCurrentProcessInfo,rDampMatrix,rRightHandSideVector);
    }

    // Rewrite local contribution into residual form (A*dx = b - A*x)
    VectorType U = ZeroVector(LocalSize);
    int LocalIndex = 0;

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        for (unsigned int d = 0; d < Dim; ++d) // Velocity Dofs
            U[LocalIndex++] = data.GetVELOCITY().Get()(i,d);
        U[LocalIndex++] = data.GetPRESSURE().Get()[i]; // Pressure Dof
    }

    noalias(rRightHandSideVector) -= prod(rDampMatrix, U);
}

template <class TElementData>
void FluidElement<TElementData>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                     ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if( rMassMatrix.size1() != LocalSize )
        rMassMatrix.resize(LocalSize,LocalSize,false);

    rMassMatrix = ZeroMatrix(LocalSize,LocalSize);

    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    TElementData data;
    data.Initialize(*this,rCurrentProcessInfo);

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        IntegrationPointGeometryData integration_point(
            GaussWeights[g], row(ShapeFunctions, g), ShapeDerivatives[g]);

        this->AddMassTerms(data,integration_point, rCurrentProcessInfo, rMassMatrix);
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
        this->CalculateProjections(rCurrentProcessInfo);
    }
}

template< class TElementData >
void FluidElement< TElementData >::EquationIdVector(EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        if (Dim == 3) rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE,ppos).EquationId();
    }
}


template< class TElementData >
void FluidElement< TElementData >::GetDofList(DofsVectorType &rElementalDofList, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();

     if (rElementalDofList.size() != LocalSize)
         rElementalDofList.resize(LocalSize);

     unsigned int LocalIndex = 0;

     for (unsigned int i = 0; i < NumNodes; ++i)
     {
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
         if (Dim == 3) rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
     }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::GetFirstDerivativesVector(Vector &rValues, int Step)
{
    GeometryType& rGeom = this->GetGeometry();

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);

    noalias(rValues) = ZeroVector(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY,Step);
        for (unsigned int d = 0; d < Dim; d++)
            rValues[Index++] = rVel[d];
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(PRESSURE,Step);
    }
}


template< class TElementData >
void FluidElement<TElementData>::GetSecondDerivativesVector(Vector &rValues, int Step)
{
    GeometryType& rGeom = this->GetGeometry();

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);

    noalias(rValues) = ZeroVector(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rAcc = rGeom[i].FastGetSolutionStepValue(ACCELERATION,Step);
        for (unsigned int d = 0; d < Dim; d++)
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
    TElementData data;
    out = data.Check(*this);

    // Extra variables used in computing projections
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(NODAL_AREA);

    // Elemental data
    KRATOS_CHECK_VARIABLE_KEY(C_SMAGORINSKY);

    for(unsigned int i=0; i<NumNodes; ++i)
    {
        Node<3>& rNode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_AREA,rNode);

        // Check that required dofs exist
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X,rNode);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y,rNode);
        if (Dim == 3) KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z,rNode);
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE,rNode);
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if ( Dim == 2)
    {
        for (unsigned int i=0; i<NumNodes; ++i)
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
    rOStream << "FluidElement" << Dim << "D";
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

    rNContainer.resize(NumGauss,NumNodes,false);
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
    for (unsigned int d = 1; d < Dim; d++)
        VelNorm += Velocity[d]*Velocity[d];
    VelNorm = std::sqrt(VelNorm);

    double InvTau = Density * ( c1 * KinematicVisc / (ElemSize*ElemSize) + c2 * VelNorm / ElemSize );
    TauOne = 1.0/InvTau;
    TauTwo = Density * (KinematicVisc + c2 * VelNorm * ElemSize / c1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////


template< class TElementData >
double FluidElement<TElementData>::EffectiveViscosity(
    TElementData& rData,
    const IntegrationPointGeometryData& rIntegrationPoint,
    double ElemSize,
    const ProcessInfo &rCurrentProcessInfo)
{
    const FluidElement* const_this = static_cast<const FluidElement*>(this);
    double Csmag = const_this->GetValue(C_SMAGORINSKY);

    double KinViscosity = rData.GetVISCOSITY().Interpolate(rIntegrationPoint.N,this);

    if (Csmag != 0.0 )
    {
        // Calculate Symetric gradient
        MatrixType S = ZeroMatrix(Dim,Dim);
        for (unsigned int n = 0; n < NumNodes; ++n)
        {
            for (unsigned int i = 0; i < Dim; ++i)
                for (unsigned int j = 0; j < Dim; ++j)
                    S(i,j) += 0.5 * ( rIntegrationPoint.DN_DX(n,j) * rData.GetVELOCITY().Get()(n,i) + rIntegrationPoint.DN_DX(n,i) * rData.GetVELOCITY().Get()(n,j) );
        }

        // Norm of symetric gradient
        double NormS = 0.0;
        for (unsigned int i = 0; i < Dim; ++i)
            for (unsigned int j = 0; j < Dim; ++j)
                NormS += S(i,j) * S(i,j);
        NormS = sqrt(2.0*NormS);

        // Nu_sgs = (Csmag * Delta)^2 * (2*Sij*Sij)^(1/2)
        KinViscosity += Csmag * Csmag * ElemSize * ElemSize * NormS;
    }

    return KinViscosity;
}


template< class TElementData >
void FluidElement<TElementData>::ResolvedConvectiveVelocity(
    TElementData& rData,
    const ShapeFunctionsType &rN,
    array_1d<double,3> &rConvVel)
{
    for (unsigned int d = 0; d < Dim; d++)
        rConvVel[d] = rN[0] * ( rData.GetVELOCITY().Get()(0,d) - rData.GetMESH_VELOCITY().Get()(0,d) );

    for (unsigned int i = 1; i < NumNodes; i++)
    {
        for (unsigned int d = 0; d < Dim; d++)
            rConvVel[d] += rN[i] * ( rData.GetVELOCITY().Get()(i,d) - rData.GetMESH_VELOCITY().Get()(i,d) );
    }
}


template< class TElementData >
void FluidElement<TElementData>::FullConvectiveVelocity(array_1d<double,3> &rConvVel,
                                       const ShapeFunctionsType &rN,
                                       const array_1d<double,3> &rSubscaleVel)
{
    GeometryType& rGeom = this->GetGeometry();

    rConvVel = rSubscaleVel;

    for (unsigned int i = 0; i < NumNodes; i++)
        rConvVel += rN[i]*(rGeom[i].FastGetSolutionStepValue(VELOCITY)-rGeom[i].FastGetSolutionStepValue(MESH_VELOCITY));
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidElement<TElementData>::ConvectionOperator(Vector &rResult,
                                   const array_1d<double,3> &rConvVel,
                                   const ShapeFunctionDerivativesType &DN_DX)
{
    if(rResult.size() != NumNodes) rResult.resize(NumNodes,false);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        rResult[i] = rConvVel[0]*DN_DX(i,0);
        for(unsigned int k = 1; k < Dim; k++)
            rResult[i] += rConvVel[k]*DN_DX(i,k);
    }
}

template <class TElementData>
void FluidElement<TElementData>::IntegrationPointVorticity(
    const ShapeFunctionDerivativesType& rDN_DX, array_1d<double, 3>& rVorticity) const
{
    rVorticity = array_1d<double, 3>(3, 0.0);

    if (Dim == 2)
    {
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            const array_1d<double, 3>& rVelocity =
                this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            rVorticity[2] += rDN_DX(i, 0) * rVelocity[1] - rDN_DX(i, 1) * rVelocity[0];
        }
    }
    else
    {
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            const array_1d<double, 3>& rVelocity =
                this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            rVorticity[0] += rDN_DX(i, 1) * rVelocity[2] - rDN_DX(i, 2) * rVelocity[1];
            rVorticity[1] += rDN_DX(i, 2) * rVelocity[0] - rDN_DX(i, 0) * rVelocity[2];
            rVorticity[2] += rDN_DX(i, 0) * rVelocity[1] - rDN_DX(i, 1) * rVelocity[0];
        }
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
template class FluidElement< DSSData2D3N >;
template class FluidElement< DSSData3D4N >;

}
