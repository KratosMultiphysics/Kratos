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

#include "utilities/math_utils.h"

#include "dss_fic.h"
#include "custom_utilities/turbulence_statistics_container.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
DSS_FIC<TDim>::DSS_FIC(IndexType NewId):
    DSS<TDim>(NewId)

{}

template< unsigned int TDim >
DSS_FIC<TDim>::DSS_FIC(IndexType NewId, const NodesArrayType& ThisNodes):
    DSS<TDim>(NewId,ThisNodes)
{}


template< unsigned int TDim >
DSS_FIC<TDim>::DSS_FIC(IndexType NewId, GeometryType::Pointer pGeometry):
    DSS<TDim>(NewId,pGeometry)
{}


template< unsigned int TDim >
DSS_FIC<TDim>::DSS_FIC(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    DSS<TDim>(NewId,pGeometry,pProperties)
{}


template< unsigned int TDim >
DSS_FIC<TDim>::~DSS_FIC()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer DSS_FIC<TDim>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DSS_FIC(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}




template< unsigned int TDim >
void DSS_FIC<TDim>::CalculateMassMatrix(MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes*(TDim+1);

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

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& rN = row(ShapeFunctions,g);
        const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

        this->AddMassTerms(GaussWeight,rN,rMassMatrix);

        /* NOTE: in FIC we have momentum stabilization terms on the mass matrix irregardless of OSS
         */
        this->AddMassStabilization(g,GaussWeight,rN,rDN_DX,rCurrentProcessInfo,rMassMatrix);
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS_FIC<TDim>::CalculateStabilizationParameters(double Density,
                                       double KinematicVisc,
                                       const array_1d<double,3> &Velocity,
                                       const ProcessInfo& rProcessInfo,
                                       double &TauIncompr,
                                       double &TauMomentum,
                                       array_1d<double,3> &TauGrad)
{
    const double c1 = 8.0;
    const double c2 = 2.0;

    const double Beta = rProcessInfo[FIC_BETA];
    const double Nobeta = 1.0-Beta;

    double Havg = this->AverageElementSize();

    double VelNorm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < TDim; d++)
        VelNorm += Velocity[d]*Velocity[d];
    VelNorm = std::sqrt(VelNorm);

    double Hvel = Havg;
    if (VelNorm > 1e-6)
    {
        Hvel = this->ProjectedElementSize(Velocity);
    }
/*
    // Velocity term in incompressibility tau is c2/t, with t = min{ h/u, dt }
    // NOW TRYING THE OPPOSITE: VelTerm = min{ c2 u/h, c2/dt }
    double VelTerm = VelNorm / Hmin;
    double TimeTerm = 1.0/rProcessInfo[DELTA_TIME];
    if (TimeTerm < VelTerm)
        VelTerm = TimeTerm;

    double InvTau = Density * ( c1 * KinematicVisc / (Hmin*Hmin) + c2 * VelTerm );
*/
    double InvTau = Density * ( c1 * KinematicVisc / (Havg*Havg) + c2 * VelNorm / Havg );
    TauIncompr = 1.0/InvTau;
    TauMomentum = (Hvel / (Density * c2 * VelNorm) );

    // TAU limiter for momentum equation: tau = min{ h/2u, dt }
    double TimeTerm = rProcessInfo[DELTA_TIME]/Density;
    if (TauMomentum > TimeTerm)
    {
        TauMomentum = TimeTerm;
    }

    TauMomentum *= Beta;

    // Coefficients for FIC shock-capturing term
    this->CalculateTauGrad(TauGrad);
    TauGrad /= Density;
    for (unsigned int d = 0; d < TDim; d++)
        if (TauGrad[d] > Havg*TimeTerm)
            TauGrad[d] = Havg*TimeTerm;

    TauGrad *= Nobeta;
}

template< unsigned int TDim >
void DSS_FIC<TDim>::CalculateTauGrad(array_1d<double,3> &TauGrad)
{
    // Small constant to prevent division by zero
    const double Small = 1e-12;

    // Evaluate velocity gradient
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    ShapeFunctionDerivativesArrayType DN_DX;
    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX,GeometryData::GI_GAUSS_1);
    ShapeFunctionDerivativesType& rDN_DX = DN_DX[0];

    BoundedMatrix<double,3,3> Gradient = ZeroMatrix(3,3);
    for (unsigned int n = 0; n < NumNodes; n++)
    {
        const array_1d<double,3>& rU = rGeom[n].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i = 0; i < TDim; i++)
        {
            for (unsigned int j = 0; j < TDim; j++)
            {
                Gradient(i,j) += rDN_DX(n,j)*rU[i];
            }
        }
    }

    // Calculate characteristic lenghts on the gradient directions and gradient norms
    array_1d<double,3> Hg(3,0.0);
    array_1d<double,3> GradNorm(3,0.0);
    for (unsigned int d = 0; d < TDim; d++)
    {
        array_1d<double,3> Gi(3,0.0);
        Gi[0] = Gradient(d,0);
        Gi[1] = Gradient(d,1);
        Gi[2] = Gradient(d,2);

        Hg[d] = this->ProjectedElementSize(Gi);

        GradNorm[d] = std::sqrt(Gi[0]*Gi[0]+Gi[1]*Gi[1]+Gi[2]*Gi[2]);

        // Coefficients for shock capturing term (remember to multiply by momentum residual!)
        TauGrad[d] = Hg[d] / (2.*GradNorm[d] + Small);
    }
}

template< unsigned int TDim >
double DSS_FIC<TDim>::AverageElementSize()
{
    KRATOS_TRY;
    GeometryType& rGeom = this->GetGeometry();
    GeometryData::KratosGeometryFamily GeoFamily = rGeom.GetGeometryFamily();

    double Havg;

    switch (GeoFamily)
    {
    case GeometryData::Kratos_Triangle:
    {
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();

        double x20 = rGeom[2].X() - rGeom[0].X();
        double y20 = rGeom[2].Y() - rGeom[0].Y();

        Havg = std::sqrt(0.5 * (x10*y20-x20*y10) );
        break;
    }
    case GeometryData::Kratos_Quadrilateral:
    {
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();

        double x30 = rGeom[3].X() - rGeom[0].X();
        double y30 = rGeom[3].Y() - rGeom[0].Y();

        Havg = std::sqrt(x10*y30-x30*y10);
        break;
    }
    case GeometryData::Kratos_Tetrahedra:
    {
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();
        double z10 = rGeom[1].Z() - rGeom[0].Z();

        double x20 = rGeom[2].X() - rGeom[0].X();
        double y20 = rGeom[2].Y() - rGeom[0].Y();
        double z20 = rGeom[2].Z() - rGeom[0].Z();

        double x30 = rGeom[3].X() - rGeom[0].X();
        double y30 = rGeom[3].Y() - rGeom[0].Y();
        double z30 = rGeom[3].Z() - rGeom[0].Z();

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
        //double Volume = detJ*0.1666666666666666666667;
        Havg = pow(detJ/6.0,1./3.);
        break;
    }
    case GeometryData::Kratos_Hexahedra:
    {
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();
        double z10 = rGeom[1].Z() - rGeom[0].Z();

        double x30 = rGeom[3].X() - rGeom[0].X();
        double y30 = rGeom[3].Y() - rGeom[0].Y();
        double z30 = rGeom[3].Z() - rGeom[0].Z();

        double x40 = rGeom[4].X() - rGeom[0].X();
        double y40 = rGeom[4].Y() - rGeom[0].Y();
        double z40 = rGeom[4].Z() - rGeom[0].Z();
/* REMOVE THIS **/
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
        Havg = std::sqrt(Hsq);
		break;
/* UP TO HERE **/

/*        double detJ = x10 * y30 * z40 - x10 * y40 * z30 + y10 * z30 * x40 - y10 * x30 * z40 + z10 * x30 * y40 - z10 * y30 * x40;
        Havg = pow(detJ,1./3.);
        break;*/

    }
    default:
    {
        KRATOS_THROW_ERROR(std::invalid_argument,"DSS_FIC::AverageElementSize not implemented for this geometry type","");
        break;
    }
    }

    return Havg;

    KRATOS_CATCH("");
}


template< unsigned int TDim >
double DSS_FIC<TDim>::ProjectedElementSize(const array_1d<double,3>& rVelocity)
{
    KRATOS_TRY;
    GeometryType& rGeom = this->GetGeometry();
    GeometryData::KratosGeometryFamily GeoFamily = rGeom.GetGeometryFamily();

    double Hvel = 0.0;

    switch (GeoFamily)
    {
    case GeometryData::Kratos_Tetrahedra:
    {
        const unsigned int NumNodes = 4;

        // Loop over edges looking for maximum 'projected' length
        array_1d<double,3> Edge(3,0.0);
        double lu = 0.0;
        for(unsigned int i = 0; i < NumNodes; ++i)
        {
            for(unsigned int j = i+1; j < NumNodes; ++j)
            {
                Edge = rGeom[j] - rGeom[i];
                lu = rVelocity[0] * Edge[0];
                for (unsigned int d = 1; d < TDim; ++d)
                    lu += rVelocity[d] * Edge[d];
                lu = fabs(lu);
                if(Hvel < lu) Hvel = lu;
            }
        }

        break;
    }
    case GeometryData::Kratos_Hexahedra:
    {
        Hvel = this->ProjectedSizeHexa(rVelocity);
        /*
        // Using edges THIS IS WHAT WAS USED FOR MY THESIS (JC)
        double x10 = rGeom[1].X() - rGeom[0].X();
        double y10 = rGeom[1].Y() - rGeom[0].Y();
        double z10 = rGeom[1].Z() - rGeom[0].Z();

        double x30 = rGeom[3].X() - rGeom[0].X();
        double y30 = rGeom[3].Y() - rGeom[0].Y();
        double z30 = rGeom[3].Z() - rGeom[0].Z();

        double x40 = rGeom[4].X() - rGeom[0].X();
        double y40 = rGeom[4].Y() - rGeom[0].Y();
        double z40 = rGeom[4].Z() - rGeom[0].Z();

        double h10 = fabs(rVelocity[0]*x10 + rVelocity[1]*y10 + rVelocity[2]*z10);
        double h30 = fabs(rVelocity[0]*x30 + rVelocity[1]*y30 + rVelocity[2]*z30);
        double h40 = fabs(rVelocity[0]*x40 + rVelocity[1]*y40 + rVelocity[2]*z40);

        Hvel = h10;
        if (Hvel < h30) Hvel = h30;
        if (Hvel < h40) Hvel = h40;
        */

        /*
        // Using diagonals
        double x60 = rGeom[6].X() - rGeom[0].X();
        double y60 = rGeom[6].Y() - rGeom[0].Y();
        double z60 = rGeom[6].Z() - rGeom[0].Z();

        double x71 = rGeom[7].X() - rGeom[1].X();
        double y71 = rGeom[7].Y() - rGeom[1].Y();
        double z71 = rGeom[7].Z() - rGeom[1].Z();

        double x24 = rGeom[2].X() - rGeom[4].X();
        double y24 = rGeom[2].Y() - rGeom[4].Y();
        double z24 = rGeom[2].Z() - rGeom[4].Z();

        double x35 = rGeom[3].X() - rGeom[5].X();
        double y35 = rGeom[3].Y() - rGeom[5].Y();
        double z35 = rGeom[3].Z() - rGeom[5].Z();

        double h60 = fabs(rVelocity[0]*x60 + rVelocity[1]*y60 + rVelocity[2]*z60);
        double h71 = fabs(rVelocity[0]*x71 + rVelocity[1]*y71 + rVelocity[2]*z71);
        double h24 = fabs(rVelocity[0]*x24 + rVelocity[1]*y24 + rVelocity[2]*z24);
        double h35 = fabs(rVelocity[0]*x35 + rVelocity[1]*y35 + rVelocity[2]*z35);

        Hvel = h60;
        if (Hvel < h71) Hvel = h71;
        if (Hvel < h24) Hvel = h24;
        if (Hvel < h35) Hvel = h35;
        */

        break;
    }
    case GeometryData::Kratos_Triangle:
    case GeometryData::Kratos_Quadrilateral:
    {
        const unsigned int NumNodes = 3;

        // Loop over edges looking for maximum 'projected' length
        array_1d<double,3> Edge(3,0.0);
        double lu = 0.0;
        for(unsigned int i = 0; i < NumNodes; ++i)
        {
            unsigned int j = (i+1) % NumNodes;
            Edge = rGeom[j] - rGeom[i];
            lu = rVelocity[0] * Edge[0];
            for (unsigned int d = 1; d < TDim; ++d)
                lu += rVelocity[d] * Edge[d];
            lu = fabs(lu);
            if(Hvel < lu) Hvel = lu;
        }

        break;
    }
    default:
    {
        KRATOS_THROW_ERROR(std::invalid_argument,"DSS_FIC::VelocityElementSize not implemented for this geometry type","");
        break;
    }
    }

    if (Hvel > 0.0)
    {
        double VelNorm = std::sqrt(rVelocity[0]*rVelocity[0] + rVelocity[1]*rVelocity[1] + rVelocity[2]*rVelocity[2]);
        Hvel /= VelNorm;
    }

    return Hvel;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points
///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS_FIC<TDim>::AddSystemTerms(unsigned int GaussIndex,
                                   double GaussWeight,
                                   const ShapeFunctionsType &rN,
                                   const ShapeFunctionDerivativesType &rDN_DX,
                                   const ProcessInfo &rProcessInfo,
                                   MatrixType &rLHS,
                                   VectorType &rRHS)
{
    // Interpolate nodal data on the integration point
    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> BodyForce(3,0.0);
    this->EvaluateInPoint(BodyForce,BODY_FORCE,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    double ElemSize = this->ElementSize();
    double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

    double TauIncompr;
    double TauMomentum;
    array_1d<double,3> TauGrad(3,0.0);
    this->CalculateStabilizationParameters(Density,Viscosity,ConvVel,rProcessInfo,TauIncompr,TauMomentum,TauGrad);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    // These two should be zero unless we are using OSS
    array_1d<double,3> MomentumProj(3,0.0);
    this->EvaluateInPoint(MomentumProj,ADVPROJ,rN);

    // Residual (used by FIC shock-capturing term)
    array_1d<double,3> MomRes(3,0.0);
    this->ASGSMomentumResidual(GaussIndex,rN,rDN_DX,MomRes);

    // Multiplying some quantities by density to have correct units
    Viscosity *= Density; // Dynamic viscosity
    BodyForce *= Density; // Force per unit of volume
    AGradN *= Density; // Convective term is always multiplied by density

    // Auxiliary variables for matrix looping
    const unsigned int NumNodes = rN.size();
    const unsigned int BlockSize = TDim+1;
    unsigned int Row = 0;
    unsigned int Col = 0;

    // Temporary containers
    double K,L,G,PDivV,qF;


    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        Row = i*BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            Col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j
            K = 0.5*(rN[i]*AGradN[j] - AGradN[i]*rN[j]); // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            K += AGradN[i]*TauMomentum*(AGradN[j]); // Stabilization: u*grad(v) * TauOne * u*grad(u)
            K *= GaussWeight;

            // q-p stabilization block (reset result)
            L = 0;

            // The following lines implement the viscous term as a Laplacian
            //for (unsigned int d = 0; d < TDim; d++)
            //    K += GaussWeight * Density * Viscosity * rDN_DX(i, d) * rDN_DX(j, d);

            for (unsigned int d = 0; d < TDim; d++)
            {
                //K += GaussWeight * Density * Viscosity * rDN_DX(i, d) * rDN_DX(j, d);
                rLHS(Row+d,Col+d) += K;

                // v * Grad(p) block
                G = AGradN[i] * rDN_DX(j,d); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                PDivV = rDN_DX(i,d) * rN[j]; // Div(v) * p

                // Write v * Grad(p) component
                rLHS(Row+d,Col+TDim) += GaussWeight * (TauMomentum*G - PDivV);
                // Use symmetry to write the q * Div(u) component
                rLHS(Col+TDim,Row+d) += GaussWeight * (TauIncompr*G + PDivV);

                // q-p stabilization block
                L += rDN_DX(i,d) * rDN_DX(j,d); // Stabilization: Grad(q) * TauOne * Grad(p)
            }

            // Write q-p term
            rLHS(Row+TDim,Col+TDim) += GaussWeight*TauIncompr*L;

            // FIC shock capturing term
            for (unsigned int d = 0; d < TDim; d++)
                rLHS(Row+d,Col+d) += GaussWeight * Density * std::fabs(TauGrad[d]*MomRes[d]) * L;
        }

        // RHS terms
        qF = 0.0;
        for (unsigned int d = 0; d < TDim; ++d)
        {
            rRHS[Row+d] += GaussWeight * rN[i] * BodyForce[d]; // v*BodyForce
            rRHS[Row+d] += GaussWeight * TauMomentum * AGradN[i] * BodyForce[d]; // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
            qF += rDN_DX(i, d) * (BodyForce[d] - MomentumProj[d]);
        }
        rRHS[Row + TDim] += GaussWeight * TauIncompr * qF; // Grad(q) * TauOne * (Density * BodyForce)
    }

    // Viscous contribution (symmetric gradient E(u) - 1/3 Tr(E) )
    // This could potentially be optimized, as it can be integrated exactly using one less integration order when compared to previous terms.
    this->AddViscousTerm(Viscosity,GaussWeight,rDN_DX,rLHS);
}


///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS_FIC<TDim>::AddMassStabilization(unsigned int GaussIndex,
                                         double GaussWeight,
                                         const ShapeFunctionsType &rN,
                                         const ShapeFunctionDerivativesType &rDN_DX,
                                         const ProcessInfo &rProcessInfo,
                                         MatrixType &rMassMatrix)
{
    // Interpolate nodal data on the integration point
    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    double ElemSize = this->ElementSize();
    double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

    double TauIncompr;
    double TauMomentum;
    array_1d<double,3> TauGrad(3,0.0);
    this->CalculateStabilizationParameters(Density,Viscosity,ConvVel,rProcessInfo,TauIncompr,TauMomentum,TauGrad);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    // Multiplying some quantities by density to have correct units
    AGradN *= Density; // Convective term is always multiplied by density

    // Auxiliary variables for matrix looping
    const unsigned int NumNodes = rN.size();
    const unsigned int BlockSize = TDim+1;
    unsigned int Row = 0;
    unsigned int Col = 0;

    // Temporary container
    double K;
    double W = GaussWeight * Density; // This density is for the dynamic term in the residual (rho*Du/Dt)
    const int oss_switch = rProcessInfo[OSS_SWITCH];

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        Row = i*BlockSize;

        for (unsigned int j = 0; j < NumNodes; j++)
        {
            Col = j*BlockSize;

            K = TauMomentum * W * AGradN[i] * rN[j];

            for (unsigned int d = 0; d < TDim; d++)
            {
                rMassMatrix(Row+d,Col+d) += K;
                if (oss_switch != 1)
                    rMassMatrix(Row+TDim,Col+d) += TauIncompr * W*rDN_DX(i,d)*rN[j];
            }
        }
    }
}


template< unsigned int TDim >
double DSS_FIC<TDim>::Module(const array_1d<double,3> &rVector)
{
    return std::sqrt(rVector[0]*rVector[0] + rVector[1]*rVector[1] + rVector[2]*rVector[2]);
}


template< unsigned int TDim >
void DSS_FIC<TDim>::Normalize(array_1d<double,3> &rVector)
{
    rVector /= this->Module(rVector);
}

template< unsigned int TDim >
double DSS_FIC<TDim>::ProjectedSizeHexa(const array_1d<double,3> &rDirection)
{
    // Logic: define a box given by hexahedra edges 10,30,40 (which I'm assuming to be orthogonal)
    // Transform the given direction to the coordinate system defined by the edges.
    // In this reference frame, place a vector aligned with rDirection on the origin and see
    // which side of the box intersects the vector (that is: which component of the vector is larger)
    // Then scale the vector to hit the limit of the box.
    // The lenght of the vector is what we are looking for.
    // NOTE: we use absolute values on all checks, this allows us to simplify the problem
    // to a single sector (otherwise we would start by determining in which of the eight sectors
    // given by +-x, +-y, +-z we have to look for the intersection).
    array_1d<double,3> U = rDirection;
    this->Normalize(U);

    Geometry< Node<3> > &rGeom = this->GetGeometry();

    array_1d<double,3> v10 = rGeom[1].Coordinates() - rGeom[0].Coordinates();
    array_1d<double,3> v30 = rGeom[3].Coordinates() - rGeom[0].Coordinates();
    array_1d<double,3> v40 = rGeom[4].Coordinates() - rGeom[0].Coordinates();

    // Express U in the coordinate system defined by {v10,v30,v40}
    Matrix Q = ZeroMatrix(3,3);
    for (unsigned int i = 0; i < 3; i++)
    {
        Q(i,0) = v10[i];
        Q(i,1) = v30[i];
        Q(i,2) = v40[i];
    }

    Matrix QInv;
    double det;
    MathUtils<double>::InvertMatrix(Q,QInv,det);

    array_1d<double,3> Uq(3,0.0);
    for (unsigned int i = 0; i < 3; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            Uq[i] += QInv(i,j)*U[j];
        }

        // Work in absolute values
        Uq[i] = std::fabs(Uq[i]);
    }

    double max_v = Uq[0];
    for (unsigned int d = 1; d < 3; d++)
        if (Uq[d] > max_v)
            max_v = Uq[d];

    double scale = 1.0/max_v;
    Uq *= scale;

    // Undo the transform
    for (unsigned int i = 0; i < 3; i++)
    {
        U[i] = 0.0;
        for (unsigned int j = 0; j < 3; j++)
        {
            U[i] += Q(i,j)*Uq[j];
        }
    }

    return this->Module(U);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< unsigned int TDim >
void DSS_FIC<TDim>::save(Serializer& rSerializer) const
{
    typedef DSS<TDim> _BaseT;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, _BaseT );
}


template< unsigned int TDim >
void DSS_FIC<TDim>::load(Serializer& rSerializer)
{
    typedef DSS<TDim> _BaseT;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, _BaseT );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class DSS_FIC<2>;
template class DSS_FIC<3>;

}
