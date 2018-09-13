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

#include "dss_ps.h"

#include "utilities/math_utils.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
DSS_PS<TDim>::DSS_PS(IndexType NewId):
    DSS<TDim>(NewId)
{}

template< unsigned int TDim >
DSS_PS<TDim>::DSS_PS(IndexType NewId, const NodesArrayType& ThisNodes):
    DSS<TDim>(NewId,ThisNodes)
{}


template< unsigned int TDim >
DSS_PS<TDim>::DSS_PS(IndexType NewId, GeometryType::Pointer pGeometry):
    DSS<TDim>(NewId,pGeometry)
{}


template< unsigned int TDim >
DSS_PS<TDim>::DSS_PS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    DSS<TDim>(NewId,pGeometry,pProperties)
{}


template< unsigned int TDim >
DSS_PS<TDim>::~DSS_PS()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer DSS_PS<TDim>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DSS_PS(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}


template< unsigned int TDim >
void DSS_PS<TDim>::CalculateLocalVelocityContribution(MatrixType &rDampMatrix, VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes*(TDim+1);

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

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& rN = row(ShapeFunctions,g);
        const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

        this->AddSystemTerms(g,GaussWeight,rN,rDN_DX,rCurrentProcessInfo,rDampMatrix,rRightHandSideVector);
    }

    // Rewrite local contribution into residual form (A*dx = b - A*x)
    VectorType U = ZeroVector(LocalSize);
    int LocalIndex = 0;

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        const array_1d<double,3> &rVel = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int d = 0; d < TDim; ++d) // Velocity Dofs
            U[LocalIndex++] = rVel[d];
        U[LocalIndex++] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
    }

    noalias(rRightHandSideVector) -= prod(rDampMatrix, U);


    this->AddPressureSubscale(rDampMatrix,rRightHandSideVector,rCurrentProcessInfo);
}


//template< unsigned int TDim >
//void DSS_PS<TDim>::CalculateMassMatrix(MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo)
//{
//    const GeometryType& rGeom = this->GetGeometry();
//    const unsigned int NumNodes = rGeom.PointsNumber();
//    const unsigned int LocalSize = NumNodes*(TDim+1);

//    // Resize and intialize output
//    if( rMassMatrix.size1() != LocalSize )
//        rMassMatrix.resize(LocalSize,LocalSize,false);

//    rMassMatrix = ZeroMatrix(LocalSize,LocalSize);

//    // Get Shape function data
//    Vector GaussWeights;
//    Matrix ShapeFunctions;
//    ShapeFunctionDerivativesArrayType ShapeDerivatives;
//    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
//    const unsigned int NumGauss = GaussWeights.size();

//    // Iterate over integration points to evaluate local contribution
//    for (unsigned int g = 0; g < NumGauss; g++)
//    {
//        const double GaussWeight = GaussWeights[g];
//        const ShapeFunctionsType& rN = row(ShapeFunctions,g);
//        const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

//        this->AddMassTerms(GaussWeight,rN,rMassMatrix);

//        /* Note on OSS and full projection: Riccardo says that adding the terms provided by
//         * AddMassStabilization (and incluiding their corresponding terms in the projeciton)
//         * could help reduce the non-linearity of the coupling between projection and u,p
//         * However, leaving them on gives a lot of trouble whith the Bossak scheme:
//         * think that we solve F - (1-alpha)*M*u^(n+1) - alpha*M*u^(n) - K(u^(n+1)) = 0
//         * so the projection of the dynamic terms should be Pi( (1-alpha)*u^(n+1) - alpha*u^(n) )
//         */
//        if ( rCurrentProcessInfo[OSS_SWITCH] != 1.0 )
//            this->AddMassStabilization(g,GaussWeight,rN,rDN_DX,rCurrentProcessInfo,rMassMatrix);
//    }

//    this->AddPressureSubscaleMass(rMassMatrix,rCurrentProcessInfo);
//}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output


template< unsigned int TDim >
std::string DSS_PS<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "DSS_PS #" << this->Id();
    return buffer.str();
}


template< unsigned int TDim >
void DSS_PS<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "DSS_PS" << TDim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////


template< unsigned int TDim >
void DSS_PS<TDim>::CalculateStaticTau(double Density,
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


template< unsigned int TDim >
void DSS_PS<TDim>::AddPressureSubscale(MatrixType& rLHS, VectorType& rRHS, ProcessInfo& rProcessInfo)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = TDim+1;
    const unsigned int LocalSize = BlockSize*NumNodes;
    const unsigned int ExtraDofs = NumNodes;

    /*
     * Full system matrix is a block matrix including extra Dofs which represent the pressure subscale.
     * These extera dofs are local to each element and will be statically condesed in this funcion.
     * Notation of the block matrix
     *
     *  | rLHS  Right |  | u,p |   | rRHS |
     *  | Down    A   | *|  ps | = |   b  |
     */
//    MatrixType Right = ZeroMatrix(LocalSize,ExtraDofs);
//    MatrixType Down = ZeroMatrix(ExtraDofs,LocalSize);
//    MatrixType A = ZeroMatrix(ExtraDofs,ExtraDofs);
//    VectorType b = ZeroVector(ExtraDofs);

//    MatrixType DynAux = ZeroMatrix(LocalSize,ExtraDofs);

    MatrixType Right = ZeroMatrix(LocalSize,ExtraDofs+1);
    MatrixType Down = ZeroMatrix(ExtraDofs+1,LocalSize);
    MatrixType A = ZeroMatrix(ExtraDofs+1,ExtraDofs+1);
    VectorType b = ZeroVector(ExtraDofs+1);

//    MatrixType DynAux = ZeroMatrix(LocalSize,ExtraDofs+1);

    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& rN = row(ShapeFunctions,g);
        const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

        // Interpolate nodal data on the integration point
        double Density;
        this->EvaluateInPoint(Density,DENSITY,rN);

        array_1d<double,3> BodyForce(3,0.0);
        this->EvaluateInPoint(BodyForce,BODY_FORCE,rN);
        //this->BodyForceTest(rProcessInfo,rN,BodyForce);

        array_1d<double,3> ConvVel(3,0.0);
        this->ResolvedConvectiveVelocity(ConvVel,rN);

        double ElemSize = this->ElementSize();
        //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
        double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

        double TauOne;
        double TauTwo;
        this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

        Vector AGradN;
        this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

        // Multiplying some quantities by density to have correct units
        Viscosity *= Density; // Dynamic viscosity
        BodyForce *= Density; // Force per unit of volume
        AGradN *= Density; // Convective term is always multiplied by density

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            for (unsigned int j = 0; j < NumNodes; j++)
            {
                //A(i,j) -= GaussWeight * rN[i] * rN[j];
                for (unsigned int d = 0; d < TDim; d++)
                {
                    A(i,j) += GaussWeight * TauOne * rDN_DX(i,d) * rDN_DX(j,d);
                    //A(i,j) += GaussWeight * rDN_DX(i,d) * TauOne * rDN_DX(j,d);

                    Down(i,BlockSize*j+d) += GaussWeight * rN[i] * rDN_DX(j,d);
                    Right(BlockSize*i+d,j) -= GaussWeight * rDN_DX(i,d) * rN[j];
//                    Right(BlockSize*i+d,j) += GaussWeight * rN[i] * rDN_DX(j,d);
//                    Down(i,BlockSize*j+d) -= GaussWeight * rDN_DX(i,d) * rN[j];

//                    Down(i,BlockSize*j+d) -= GaussWeight * rDN_DX(i,d) * TauOne * AGradN[j];
//                    Down(i,BlockSize*j+TDim) -= GaussWeight * rDN_DX(i,d) * TauOne * rDN_DX(j,d);

//                    DynAux(i,BlockSize*j+d) -= GaussWeight * rDN_DX(i,d) * TauOne * rN[j];

//                    Right(BlockSize*i+d,j) -= GaussWeight * rDN_DX(j,d) * TauOne * AGradN[i];
//                    Right(BlockSize*i+TDim,j) -= GaussWeight * rDN_DX(i,d) * TauOne * rDN_DX(j,d);
//                    Right(BlockSize*i+d,j) += GaussWeight * AGradN[i] * rDN_DX(j,d);
                    //Right(BlockSize*i+TDim,j) += GaussWeight * rDN_DX(i,d) * rDN_DX(j,d);

                    //Right(BlockSize*i+d,j) -= GaussWeight * rDN_DX(i,d) * rN[j];

//                    b[i] -= GaussWeight * rDN_DX(i,d) * TauOne * BodyForce[d];
                }
            }
        }
    }

    // Impose a restriction so that the local laplacian can be inverted: average of the solution of Ax = b is zero
    for (unsigned int i = 0; i < ExtraDofs; i++)
    {
        A(ExtraDofs,i) = 1.0;
        A(i,ExtraDofs) = 1.0;
    }


    /* FOR BDF2
    // Dynamic contribution
    const Vector& rBDFcoeffs = rProcessInfo[BDF_COEFFICIENTS];
    Down += rBDFcoeffs[0] * DynAux;

    for (unsigned int step = 1; step < 3; step++ )
    {
        VectorType U = ZeroVector(LocalSize);
        int LocalIndex = 0;

        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            const array_1d<double,3> &rVel = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,step);
            for (unsigned int d = 0; d < TDim; ++d) // Velocity Dofs
                U[LocalIndex++] = rVel[d];
            //U[LocalIndex++] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE,step); // Pressure Dof
            LocalIndex++; // nothing dynamic about pressure
        }

        U *= rBDFcoeffs[step];
        noalias(b) -= prod(DynAux, U);
    }
*/

    // Rewrite local contribution into residual form (A*dx = b - A*x)
    VectorType U = ZeroVector(LocalSize);
    int LocalIndex = 0;

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        const array_1d<double,3> &rVel = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int d = 0; d < TDim; ++d) // Velocity Dofs
            U[LocalIndex++] = rVel[d];
        U[LocalIndex++] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
    }

    noalias(b) -= prod(Down, U);


    // Condense the extra block back into the system
    MatrixType InvA;
    double det;
    MathUtils<double>::InvertMatrix(A,InvA,det);

    rLHS -= prod( Right, MatrixType( prod(InvA,Down)) );
    rRHS -= prod( Right, VectorType( prod(InvA,b)) );
}



template< unsigned int TDim >
void DSS_PS<TDim>::AddPressureSubscaleMass(MatrixType& rLHS, ProcessInfo& rProcessInfo)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = TDim+1;
    const unsigned int LocalSize = BlockSize*NumNodes;
    const unsigned int ExtraDofs = NumNodes;

    /*
     * Full system matrix is a block matrix including extra Dofs which represent the pressure subscale.
     * These extra dofs are local to each element and will be statically condesed in this funcion.
     * Notation of the block matrix
     *
     *  | rLHS  Right |  | u,p |   | rRHS |
     *  | Down    A   | *|  ps | = |   b  |
     */
//    MatrixType Right = ZeroMatrix(LocalSize,ExtraDofs);
//    MatrixType Down = ZeroMatrix(ExtraDofs,LocalSize);
//    MatrixType A = ZeroMatrix(ExtraDofs,ExtraDofs);
//    VectorType b = ZeroVector(ExtraDofs);

//    MatrixType DynAux = ZeroMatrix(LocalSize,ExtraDofs);

    MatrixType Right = ZeroMatrix(LocalSize,ExtraDofs+1);
    MatrixType Down = ZeroMatrix(ExtraDofs+1,LocalSize);
    MatrixType A = ZeroMatrix(ExtraDofs+1,ExtraDofs+1);
    VectorType b = ZeroVector(ExtraDofs+1);

    // Get Shape function data
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& rN = row(ShapeFunctions,g);
        const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

        // Interpolate nodal data on the integration point
        double Density;
        this->EvaluateInPoint(Density,DENSITY,rN);

        array_1d<double,3> BodyForce(3,0.0);
        this->EvaluateInPoint(BodyForce,BODY_FORCE,rN);
        //this->BodyForceTest(rProcessInfo,rN,BodyForce);

        array_1d<double,3> ConvVel(3,0.0);
        this->ResolvedConvectiveVelocity(ConvVel,rN);

        double ElemSize = this->ElementSize();
        //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
        double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

        double TauOne;
        double TauTwo;
        this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

        Vector AGradN;
        this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

        // Multiplying some quantities by density to have correct units
        Viscosity *= Density; // Dynamic viscosity
        BodyForce *= Density; // Force per unit of volume
        AGradN *= Density; // Convective term is always multiplied by density

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            for (unsigned int j = 0; j < NumNodes; j++)
            {
                //A(i,j) -= GaussWeight * rN[i] * rN[j];
                for (unsigned int d = 0; d < TDim; d++)
                {
                    A(i,j) += GaussWeight * rDN_DX(i,d) * rDN_DX(j,d);
                    //A(i,j) += GaussWeight * rDN_DX(i,d) * TauOne * rDN_DX(j,d);

                    Down(i,BlockSize*j+d) -= GaussWeight * rDN_DX(i,d) * TauOne * rN[j];

                    //Right(BlockSize*j+d,i) -= GaussWeight * rDN_DX(i,d) * TauOne * AGradN[j];
                    //Right(BlockSize*j+d,i) -= GaussWeight * rDN_DX(i,d) * TauOne * rDN_DX(j,d);
                    Right(BlockSize*j+d,i) -= GaussWeight * rDN_DX(i,d) * AGradN[j];
                    Right(BlockSize*j+d,i) -= GaussWeight * rDN_DX(i,d) * rDN_DX(j,d);

                    b[i] -= GaussWeight * rDN_DX(i,d) * TauOne * BodyForce[d];
                }
            }
        }
    }

    // Impose a restriction so that the local laplacian can be inverted: average of the solution of Ax = b is zero
    for (unsigned int i = 0; i < ExtraDofs; i++)
    {
        A(ExtraDofs,i) = 1.0;
        A(i,ExtraDofs) = 1.0;
    }

    // Condense the extra block back into the system
    MatrixType InvA;
    double det;
    MathUtils<double>::InvertMatrix(A,InvA,det);

    rLHS -= prod( Right, MatrixType( prod(InvA,Down)) );
}



///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< unsigned int TDim >
void DSS_PS<TDim>::save(Serializer& rSerializer) const
{
    typedef DSS<TDim> _basetype;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, _basetype );
}


template< unsigned int TDim >
void DSS_PS<TDim>::load(Serializer& rSerializer)
{
    typedef DSS<TDim> _basetype;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, _basetype);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class DSS_PS<2>;
template class DSS_PS<3>;

}
