#include "dss_fic_limited.h"
#include "custom_utilities/turbulence_statistics_container.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
DSS_FIC_LIMITED<TDim>::DSS_FIC_LIMITED(IndexType NewId):
    DSS_FIC<TDim>(NewId),
    mBeta(3,0.0)
{}

template< unsigned int TDim >
DSS_FIC_LIMITED<TDim>::DSS_FIC_LIMITED(IndexType NewId, const NodesArrayType& ThisNodes):
    DSS_FIC<TDim>(NewId,ThisNodes),
    mBeta(3,0.0)
{}


template< unsigned int TDim >
DSS_FIC_LIMITED<TDim>::DSS_FIC_LIMITED(IndexType NewId, GeometryType::Pointer pGeometry):
    DSS_FIC<TDim>(NewId,pGeometry),
    mBeta(3,0.0)
{}


template< unsigned int TDim >
DSS_FIC_LIMITED<TDim>::DSS_FIC_LIMITED(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    DSS_FIC<TDim>(NewId,pGeometry,pProperties),
    mBeta(3,0.0)
{}


template< unsigned int TDim >
DSS_FIC_LIMITED<TDim>::~DSS_FIC_LIMITED()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer DSS_FIC_LIMITED<TDim>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DSS_FIC_LIMITED(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

template< unsigned int TDim >
Element::Pointer DSS_FIC_LIMITED<TDim>::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DSS_FIC_LIMITED(NewId, pGeom, pProperties));
}

template< unsigned int TDim >
void DSS_FIC_LIMITED<TDim>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    const double BetaLimit = rCurrentProcessInfo[FIC_BETA];

    // Evaluate velocity gradient
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    Vector DetJ;
    ShapeFunctionDerivativesArrayType DN_DX;
    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX,DetJ,GeometryData::GI_GAUSS_1);
    ShapeFunctionDerivativesType& rDN_DX = DN_DX[0];

    Matrix NContainer;
    NContainer.resize(1,NumNodes,false);
    NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

    Matrix Gradient = ZeroMatrix(3,3);
    array_1d<double,3> Vel(3,0.0);
    for (unsigned int n = 0; n < NumNodes; n++)
    {
        const array_1d<double,3>& rU = rGeom[n].FastGetSolutionStepValue(VELOCITY);
        Vel += rU * NContainer(0,n);
        for (unsigned int i = 0; i < TDim; i++)
        {
            for (unsigned int j = 0; j < TDim; j++)
            {
                Gradient(i,j) += rDN_DX(n,j)*rU[i];
            }
        }
    }

    double VelNorm = std::sqrt(Vel[0]*Vel[0] + Vel[1]*Vel[1] + Vel[2]*Vel[2]);
    for (unsigned int d = 0; d < TDim; d++)
    {
        double Prod = std::fabs( Vel[0]*Gradient(d,0) + Vel[1]*Gradient(d,1) + Vel[2]*Gradient(d,2) );
        double Norm = std::sqrt( Gradient(d,0)*Gradient(d,0) + Gradient(d,1)*Gradient(d,1) + Gradient(d,2)*Gradient(d,2));

        mBeta[d] = ( VelNorm*Norm > 1e-6 ) ? 1.0 - Prod / (Norm*VelNorm) : 1.0;
        if (mBeta[d] < BetaLimit) mBeta[d] = BetaLimit;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////


template< unsigned int TDim >
void DSS_FIC_LIMITED<TDim>::CalculateAllTaus(double Density,
                                             double KinematicVisc,
                                             const array_1d<double,3>& Velocity,
                                             const ProcessInfo& rProcessInfo,
                                             double& TauIncompr,
                                             array_1d<double,3>& TauMomentum,
                                             array_1d<double,3>& TauGrad)
{
    const double c1 = 8.0;
    const double c2 = 2.0;

    // Velocity norm
    double VelNorm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < TDim; d++)
        VelNorm += Velocity[d]*Velocity[d];
    VelNorm = std::sqrt(VelNorm);

    // Average Element Size
    double Havg = this->AverageElementSize();

    // Element size in the direction of velocity
    double Hvel = Havg;
    if (VelNorm > 1e-6)
    {
        Hvel = this->ProjectedElementSize(Velocity);
    }

    // Base tau for the momentum equation
    // (We will multiply it by a directional coefficient later)
    double TauM = (Hvel / (Density * c2 * VelNorm) );

    // TAU limiter for momentum equation: tau = min{ h/2u, dt }
    double TimeTerm = rProcessInfo[DELTA_TIME]/Density;
    if (TauM > TimeTerm)
    {
        TauM = TimeTerm;
    }


    // Coefficients for FIC shock-capturing term
    this->CalculateTauGrad(TauGrad);
    TauGrad /= Density;
    for (unsigned int d = 0; d < TDim; d++)
        if (TauGrad[d] > Havg*TimeTerm)
            TauGrad[d] = Havg*TimeTerm;


    for( unsigned int d = 0; d < TDim; d++)
    {
        // Momentum tau for direction d
        TauMomentum[d] = mBeta[d] * TauM;

        // Gradient term for direction d
        TauGrad[d] *= 1.0-mBeta[d];
    }


    // GLS tau for the incompressibility equation
    double InvTau = Density * ( c1 * KinematicVisc / (Havg*Havg) + c2 * VelNorm / Havg );
    TauIncompr = 1.0/InvTau;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points
///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS_FIC_LIMITED<TDim>::AddSystemTerms(unsigned int GaussIndex,
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
    array_1d<double,3> TauMomentum(3,0.0);
    array_1d<double,3> TauGrad(3,0.0);
    this->CalculateAllTaus(Density,Viscosity,ConvVel,rProcessInfo,TauIncompr,TauMomentum,TauGrad);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    // These two should be zero unless we are using OSS
    array_1d<double,3> MomentumProj(3,0.0);
    double MassProj = 0;
    this->EvaluateInPoint(MomentumProj,ADVPROJ,rN);
    this->EvaluateInPoint(MassProj,DIVPROJ,rN);

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
            //K += AGradN[i]*TauMomentum*(AGradN[j]); // Stabilization: u*grad(v) * TauOne * u*grad(u)
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
                rLHS(Row+d,Col+d) += GaussWeight * AGradN[i]*TauMomentum[d]*(AGradN[j]); // Stabilization: u*grad(v) * TauOne * u*grad(u)

                // v * Grad(p) block
                G = AGradN[i] * rDN_DX(j,d); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                PDivV = rDN_DX(i,d) * rN[j]; // Div(v) * p

                // Write v * Grad(p) component
                rLHS(Row+d,Col+TDim) += GaussWeight * (TauMomentum[d]*G - PDivV);
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
            rRHS[Row+d] += GaussWeight * TauMomentum[d] * AGradN[i] * BodyForce[d]; // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
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
void DSS_FIC_LIMITED<TDim>::AddMassStabilization(unsigned int GaussIndex,
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
    array_1d<double,3> TauMomentum(3,0.0);
    array_1d<double,3> TauGrad(3,0.0);
    this->CalculateAllTaus(Density,Viscosity,ConvVel,rProcessInfo,TauIncompr,TauMomentum,TauGrad);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    // Multiplying some quantities by density to have correct units
    AGradN *= Density; // Convective term is always multiplied by density

    // Auxiliary variables for matrix looping
    const unsigned int NumNodes = rN.size();
    const unsigned int BlockSize = TDim+1;
    unsigned int Row = 0;
    unsigned int Col = 0;

    double W = GaussWeight * Density; // This density is for the dynamic term in the residual (rho*Du/Dt)

	const int oss_switch = rProcessInfo[OSS_SWITCH];

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        Row = i*BlockSize;

        for (unsigned int j = 0; j < NumNodes; j++)
        {
            Col = j*BlockSize;


            for (unsigned int d = 0; d < TDim; d++)
            {
                rMassMatrix(Row+d,Col+d) += TauMomentum[d] * W * AGradN[i] * rN[j];
				if (oss_switch != 1)
					rMassMatrix(Row+TDim,Col+d) += TauIncompr * W*rDN_DX(i,d)*rN[j];
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< unsigned int TDim >
void DSS_FIC_LIMITED<TDim>::save(Serializer& rSerializer) const
{
    typedef DSS_FIC<TDim> _BaseT;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, _BaseT );
}


template< unsigned int TDim >
void DSS_FIC_LIMITED<TDim>::load(Serializer& rSerializer)
{
    typedef DSS_FIC<TDim> _BaseT;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, _BaseT );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class DSS_FIC_LIMITED<2>;
template class DSS_FIC_LIMITED<3>;

}
