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


#include "dynss.h"
#include "custom_utilities/turbulence_statistics_container.h"
#include "utilities/math_utils.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
DynSS<TDim>::DynSS(IndexType NewId):
    DSS<TDim>(NewId),
    mPredSsVel(),
    mOldSsVel()
{}

template< unsigned int TDim >
DynSS<TDim>::DynSS(IndexType NewId, const NodesArrayType& ThisNodes):
    DSS<TDim>(NewId,ThisNodes),
    mPredSsVel(),
    mOldSsVel()
{}


template< unsigned int TDim >
DynSS<TDim>::DynSS(IndexType NewId, GeometryType::Pointer pGeometry):
    DSS<TDim>(NewId,pGeometry),
    mPredSsVel(),
    mOldSsVel()
{}


template< unsigned int TDim >
DynSS<TDim>::DynSS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    DSS<TDim>(NewId,pGeometry,pProperties),
    mPredSsVel(),
    mOldSsVel()
{}


template< unsigned int TDim >
DynSS<TDim>::~DynSS()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer DynSS<TDim>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DynSS(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

template< unsigned int TDim >
Element::Pointer DynSS<TDim>::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DynSS(NewId, pGeom, pProperties));
}

template< unsigned int TDim >
void DynSS<TDim>::Initialize()
{
    unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());

    // I'm checking the size because I could be doing this after a restart
    // and I want to keep loaded values in that case
    if (mOldSsVel.size() != NumGauss)
    {
        mOldSsVel.resize(0);
        mOldSsVel.reserve(NumGauss);
        for (unsigned int g = 0; g < NumGauss; g++)
            mOldSsVel.push_back( array_1d<double,3>(3,0.0) );
    }

    // The prediction is updated before each non-linear iteration, so it is not stored
    mPredSsVel.resize(0);
    mPredSsVel.reserve(NumGauss);
    for (unsigned int g = 0; g < NumGauss; g++)
        mPredSsVel.push_back( array_1d<double,3>(3,0.0) );
}


template< unsigned int TDim >
void DynSS<TDim>::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    // Update the old small scale value for next iteration
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const ShapeFunctionsType& N = row(ShapeFunctions,g);
        const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

        // Not doing the update "in place" because SubscaleVelocity uses mOldSsVel
        // (although I probably could do it safely)
        array_1d<double,3> UpdatedValue(3,0.0);
        this->SubscaleVelocity(g,N,rDN_DX,rCurrentProcessInfo,UpdatedValue);
        mOldSsVel[g] = UpdatedValue;
    }
}


template< unsigned int TDim >
void DynSS<TDim>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{
/*
    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    // Update the old small scale value for next iteration
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const ShapeFunctionsType& N = row(ShapeFunctions,g);
        const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

        // Not doing the update "in place" because SubscaleVelocity uses mOldSsVel
        // (although I probably could do it safely)
        array_1d<double,3> UpdatedValue(3,0.0);
        this->SubscaleVelocity(g,N,rDN_DX,rCurrentProcessInfo,UpdatedValue);
        mPredSsVel[g] = UpdatedValue;
    }
    */

    // Get Shape function data
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    Vector GaussWeights;
    Matrix ShapeFunctions;
    ShapeFunctionDerivativesArrayType ShapeDerivatives;
    this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
    const unsigned int NumGauss = GaussWeights.size();

    double Density = 0.0;
    array_1d<double,3> CoarseConvVel(3,0.0);

    const double Dt = rCurrentProcessInfo.GetValue(DELTA_TIME);

    // A non-linear equation must be solved for each integration point
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const ShapeFunctionsType& N = row(ShapeFunctions,g);
        const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

        // Elemental large-scale velocity gradient
        Matrix CoarseVelGradient = ZeroMatrix(TDim,TDim);
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            const array_1d<double,3>& rCoarseVelocity = rGeom[i].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int m = 0; m < TDim; m++)
                for (unsigned int n = 0; n < TDim; n++)
                    CoarseVelGradient(m,n) += rDN_DX(i,n) * rCoarseVelocity[m];
        }

        const array_1d<double,3>& OldSubscale = mOldSsVel[g];

        this->ResolvedConvectiveVelocity(CoarseConvVel,N);
        double h = this->ElementSize();
        //double h = this->ElementSize(CoarseConvVel+OldSubscale,rDN_DX);

        this->EvaluateInPoint(Density,DENSITY,N);
        const double Viscosity = this->EffectiveViscosity(N,rDN_DX,h,rCurrentProcessInfo);

        // Part of the residual that does not depend on the subscale
        array_1d<double,3> StaticResidual(3,0.0);
        // Note I'm only using large scale convection here, small-scale convection is re-evaluated at each iteration.
        array_1d<double,3> ResolvedConvVel(3,0.0);
        this->ResolvedConvectiveVelocity(ResolvedConvVel,N);

        if (rCurrentProcessInfo[OSS_SWITCH] != 1.0)
            this->DASGSMomentumResidual(N,rDN_DX,ResolvedConvVel,StaticResidual);
        else
            this->DOSSMomentumResidual(N,rDN_DX,ResolvedConvVel,StaticResidual);


        // Add the time discretization term to obtain the part of the residual that does not change during iteration
        const double MassTerm = Density / Dt;
        for (unsigned int d = 0; d < TDim; d++)
            StaticResidual[d] += MassTerm * OldSubscale[d];

        ///@todo: redefine
        const double SubscaleTol = 1e-14;
        const double SubscaleRHSTol = 1e-14;

        // Newton-Raphson iterations for the subscale
        unsigned int Iter = 0;
        double SubscaleError = 2.0 * SubscaleTol;
        double SubscaleNorm = 0.0;
        double ResidualNorm = 2.0 * SubscaleRHSTol;
        double InvTau = 0.0;
        Matrix J = ZeroMatrix(TDim,TDim);
        Vector RHS = ZeroVector(TDim);
        Vector U = ZeroVector(TDim);
        Vector dU = ZeroVector(TDim);
        // Use last result as initial guess
        const array_1d<double,3>& rSubscale = mPredSsVel[g];
        for(unsigned int d = 0; d < TDim; d++)
            U[d] = rSubscale[d];


        while (Iter < 10 && SubscaleError > SubscaleTol )
        {
            // Update iteration counter
            ++Iter;

            // Calculate new Tau
            double ConvVelNorm = 0.0;
            for (unsigned int d = 0; d < TDim; d++)
            {
                double Vd = CoarseConvVel[d] + U[d];
                ConvVelNorm += Vd*Vd;
            }
            ConvVelNorm = sqrt(ConvVelNorm);
            InvTau = Density * ( 1.0/Dt + 8.0*Viscosity/(h*h) + 2.0*ConvVelNorm/h );

            // Newton-Raphson LHS
            noalias(J) = Density * CoarseVelGradient;
            for (unsigned int d = 0; d < TDim; d++)
                J(d,d) += InvTau;

            // Newton-Raphson RHS
            //noalias(RHS) = StaticResidual;
            for (unsigned int d = 0; d < TDim; d++)
                RHS[d] = StaticResidual[d];
            noalias(RHS) -= prod(J,U);

            ResidualNorm = RHS[0]*RHS[0];
            for (unsigned int d = 1; d < TDim; d++)
                ResidualNorm += RHS[d]*RHS[d];

            if (ResidualNorm > SubscaleRHSTol)
            {
                this->DenseSystemSolve(J,RHS,dU);

                // Update
                noalias(U) += dU;

                // Convergence check
                SubscaleError = dU[0]*dU[0];
                SubscaleNorm = U[0]*U[0];
                for(unsigned int d = 1; d < TDim; ++d)
                {
                    SubscaleError += dU[d]*dU[d];
                    SubscaleNorm += U[d]*U[d];
                }
                SubscaleError /= SubscaleNorm;
            }
            else
            {
                break; // If RHS is zero, dU is zero too, converged.
            }
        }

        //KRATOS_WATCH(Iter)

//        std::cout << "Id: " << this->Id() << " g: " << g << " iter: " << Iter << " ss err: " << SubscaleError << " ss norm: " << SubscaleNorm << " residual: " << ResidualNorm << std::endl;

        // Store new subscale values
        array_1d<double,3>& rOut = mPredSsVel[g];
        for(unsigned int d = 0; d < TDim; d++)
            rOut[d] = U[d];
    }
}

template< unsigned int TDim >
void DynSS<TDim>::GetValueOnIntegrationPoints(const Variable<double> &rVariable,
                                              std::vector<double> &rValues,
                                              const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == SUBSCALE_PRESSURE)
    {
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        const unsigned int NumGauss = GaussWeights.size();

        rValues.resize(NumGauss);

        // Update the old small scale value for next iteration
        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& N = row(ShapeFunctions,g);
            const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

            this->SubscalePressure(g,N,rDN_DX,rCurrentProcessInfo,rValues[g]);
        }
    }
    else
    {
        DSS<TDim>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}

template< unsigned int TDim >
void DynSS<TDim>::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> > &rVariable,
                                              std::vector<array_1d<double,3> > &rValues,
                                              const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == SUBSCALE_VELOCITY)
    {
        unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
        rValues.resize(NumGauss);
        // check if initialized
        if (mOldSsVel.size() == NumGauss)
        {
            for (unsigned int g = 0; g < NumGauss; g++)
                rValues[g] = mOldSsVel[g];
        }
        else
        {
            for (unsigned int g = 0; g < NumGauss; g++)
                rValues[g] = array_1d<double,3>(3,0.0);
        }
    }
    else
    {
        DSS<TDim>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}

/**
 * TO COPY SUBSCALE VALUES AFTER SERIALIZATION
 */
template< unsigned int TDim >
void DynSS<TDim>::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                      std::vector<double>& rOutput,
                      const ProcessInfo& rCurrentProcessInfo)
{
}

template< unsigned int TDim >
void DynSS<TDim>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                               std::vector< array_1d<double, 3 > >& rOutput,
                                               const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == SUBSCALE_VELOCITY)
        this->GetValueOnIntegrationPoints(SUBSCALE_VELOCITY,rOutput,rCurrentProcessInfo);
}

template< unsigned int TDim >
void DynSS<TDim>::CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
                      std::vector< Vector >& rOutput,
                      const ProcessInfo& rCurrentProcessInfo)
{
}

template< unsigned int TDim >
void DynSS<TDim>::CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable,
                      std::vector< Matrix >& rOutput,
                      const ProcessInfo& rCurrentProcessInfo)
{
}



/**
 * TO COPY SUBSCALE VALUES AFTER SERIALIZATION
 */
template< unsigned int TDim >
void DynSS<TDim>::SetValueOnIntegrationPoints(const Variable<double>& rVariable,
                     std::vector<double>& rValues,
                     const ProcessInfo& rCurrentProcessInfo)
{
}

template< unsigned int TDim >
void DynSS<TDim>::SetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                     std::vector<array_1d<double, 3 > > rValues,
                     const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    if (rVariable == SUBSCALE_VELOCITY)
    {
        unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() == NumGauss)
        {
            mOldSsVel.resize(NumGauss);
            for (unsigned int g = 0; g < NumGauss; g++)
                mOldSsVel[g] = rValues[g];

        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error, "Getting something strange as old velocity subscale!!!!","");
        }
    }
    KRATOS_CATCH("")
}


template< unsigned int TDim >
void DynSS<TDim>::SetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                     std::vector<array_1d<double, 6 > > rValues,
                     const ProcessInfo& rCurrentProcessInfo)
{
}

template< unsigned int TDim >
void DynSS<TDim>::SetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                     std::vector<Vector>& rValues,
                     const ProcessInfo& rCurrentProcessInfo)
{
}

template< unsigned int TDim >
void DynSS<TDim>::SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                     std::vector<Matrix>& rValues,
                     const ProcessInfo& rCurrentProcessInfo)
{
}

template< unsigned int TDim >
void DynSS<TDim>::SetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                     std::vector<ConstitutiveLaw::Pointer>& rValues,
                     const ProcessInfo& rCurrentProcessInfo)
{
}



///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output


template< unsigned int TDim >
std::string DynSS<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "DynSS #" << this->Id();
    return buffer.str();
}


template< unsigned int TDim >
void DynSS<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "DynSS" << TDim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////



template< unsigned int TDim >
void DynSS<TDim>::CalculateTau(double Density,
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
    TauTwo = Density * (KinematicVisc + c2 * VelNorm * ElemSize / c1);
    //TauTwo = 0.0;

    // Auxiliary coefficient StaticTauOne*TauTwo/Dt that appears on the pressure subscale model
    //TauP = Density * ElemSize*ElemSize / (c1*dt);
    TauP = 0.0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DynSS<TDim>::DASGSMomentumResidual(const ShapeFunctionsType &rN,
                                       const ShapeFunctionDerivativesType &rDN_DX,
                                       const array_1d<double,3>& rConvVel,
                                       array_1d<double,3> &rMomentumRes)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    Vector AGradN;
    this->ConvectionOperator(AGradN,rConvVel,rDN_DX);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rBodyForce = rGeom[i].FastGetSolutionStepValue(BODY_FORCE);
        const array_1d<double,3>& rAcc = rGeom[i].FastGetSolutionStepValue(ACCELERATION);
        const array_1d<double,3>& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        const double Press = rGeom[i].FastGetSolutionStepValue(PRESSURE);

        for (unsigned int d = 0; d < TDim; d++)
        {
            rMomentumRes[d] += Density * ( rN[i]*(rBodyForce[d] - rAcc[d]) - AGradN[i]*rVel[d]) - rDN_DX(i,d)*Press;
        }
    }
}


template< unsigned int TDim >
void DynSS<TDim>::ASGSMassResidual(double GaussIndex,
                                   const ShapeFunctionsType &rN,
                                   const ShapeFunctionDerivativesType &rDN_DX,
                                   double &rMomentumRes)
{
    this->MassProjTerm(GaussIndex,rN,rDN_DX,rMomentumRes);
}


template< unsigned int TDim >
void DynSS<TDim>::DOSSMomentumResidual(const ShapeFunctionsType &rN,
                                      const ShapeFunctionDerivativesType &rDN_DX,
                                      const array_1d<double,3>& rConvVel,
                                      array_1d<double,3> &rMomentumRes)
{
    this->DynamicMomentumProjTerm(rN,rDN_DX,rConvVel,rMomentumRes);

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rProj = rGeom[i].FastGetSolutionStepValue(ADVPROJ);

        for (unsigned int d = 0; d < TDim; d++)
            rMomentumRes[d] -= rN[i]*rProj[d];
    }
}


template< unsigned int TDim >
void DynSS<TDim>::OSSMassResidual(double GaussIndex,
                                  const ShapeFunctionsType &rN,
                                  const ShapeFunctionDerivativesType &rDN_DX,
                                  double &rMassRes)
{
    this->MassProjTerm(GaussIndex,rN,rDN_DX,rMassRes);

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const double Proj = rGeom[i].FastGetSolutionStepValue(DIVPROJ);
        rMassRes -= rN[i]*Proj;
    }
}


template< unsigned int TDim >
void DynSS<TDim>::DynamicMomentumProjTerm(const ShapeFunctionsType &rN,
                                   const ShapeFunctionDerivativesType &rDN_DX,
                                   const array_1d<double,3> &rConvVel,
                                   array_1d<double,3> &rMomentumRHS)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    Vector AGradN;
    this->ConvectionOperator(AGradN,rConvVel,rDN_DX);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rBodyForce = rGeom[i].FastGetSolutionStepValue(BODY_FORCE);
        const array_1d<double,3>& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        const double Press = rGeom[i].FastGetSolutionStepValue(PRESSURE);

        for (unsigned int d = 0; d < TDim; d++)
        {
            rMomentumRHS[d] += Density * ( rN[i]*rBodyForce[d] - AGradN[i]*rVel[d]) - rDN_DX(i,d)*Press;
        }
    }
}


template< unsigned int TDim >
void DynSS<TDim>::MassProjTerm(double GaussIndex,
                               const ShapeFunctionsType &rN,
                               const ShapeFunctionDerivativesType &rDN_DX,
                               double &rMassRHS)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);

        for (unsigned int d = 0; d < TDim; d++)
            rMassRHS -= rDN_DX(i,d)*rVel[d];
    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////


template< unsigned int TDim >
void DynSS<TDim>::FullConvectiveVelocity(array_1d<double,3> &rConvVel,
                                         const ShapeFunctionsType &rN,
                                         const array_1d<double,3> &rSubscaleVel)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    rConvVel = rSubscaleVel;

    for (unsigned int i = 0; i < NumNodes; i++)
        rConvVel += rN[i]*(rGeom[i].FastGetSolutionStepValue(VELOCITY)-rGeom[i].FastGetSolutionStepValue(MESH_VELOCITY));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points
///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DynSS<TDim>::AddSystemTerms(unsigned int GaussIndex,
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
    //this->BodyForceTest(rProcessInfo,rN,BodyForce);

    array_1d<double,3> FullConvVel(3,0.0);
    this->FullConvectiveVelocity(FullConvVel,rN,mPredSsVel[GaussIndex]);

    double ElemSize = this->ElementSize();
    //double ElemSize = this->ElementSize(FullConvVel,rDN_DX);
    double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    double TauP;
    this->CalculateTau(Density,Viscosity,FullConvVel,rProcessInfo,ElemSize,TauOne,TauTwo,TauP);
    //this->CalculateTau(Density,Viscosity,FullConvVel,rProcessInfo,TauOne,TauTwo,TauP);
    double Dt = rProcessInfo[DELTA_TIME];

    // small scale velocity contributions (subscale tracking)
    array_1d<double,3> OldUssTerm = (Density/Dt) * mOldSsVel[GaussIndex]; // rho * u_ss^{n-1}/dt

    // Old mass residual for dynamic pressure subscale: -Div(u^n)
    double OldResidual = 0.0;
    for (unsigned int a = 0; a < rN.size(); a++)
    {
        const array_1d<double,3>& rOldVel = this->GetGeometry()[a].FastGetSolutionStepValue(VELOCITY,1);
        double OldDivProj = this->GetGeometry()[a].FastGetSolutionStepValue(DIVPROJ,1);
        for (unsigned int d = 0; d < TDim; d++)
            OldResidual -= rDN_DX(a,d)*rOldVel[d] + rN[a]*OldDivProj;
    }

    Vector AGradN;
    this->ConvectionOperator(AGradN,FullConvVel,rDN_DX);

    // These two should be zero unless we are using OSS
    array_1d<double,3> MomentumProj(3,0.0);
    double MassProj = 0;
    this->EvaluateInPoint(MomentumProj,ADVPROJ,rN);
    this->EvaluateInPoint(MassProj,DIVPROJ,rN);

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
    double K,L,G,PDivV,qF,Sg;


    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        Row = i*BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            Col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j

            // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            K = 0.5*(rN[i]*AGradN[j] - AGradN[i]*rN[j]);
            //K = rN[i]*AGradN[j];

            // Stabilization: u*grad(v) * TauOne * u*grad(u) - vh * TauOne/Dt u*grad(u)
            // The last term comes from vh*d(u_ss)/dt
            K += (AGradN[i] - Density*rN[i]/Dt)*TauOne*(AGradN[j]);
            K *= GaussWeight;

            // q-p stabilization block (reset result)
            L = 0;

            // The following lines implement the viscous term as a Laplacian
            //for (unsigned int d = 0; d < TDim; d++)
            //    K += GaussWeight * Density * Viscosity * rDN_DX(i, d) * rDN_DX(j, d);

            for (unsigned int d = 0; d < TDim; d++)
            {
                rLHS(Row+d,Col+d) += K;

                /* v * Grad(p) block */
                // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                G = TauOne * AGradN[i] * rDN_DX(j,d);

                // From vh*d(u_ss)/dt: vh * TauOne/Dt * Grad(p)
                Sg = TauOne*Density*rN[i]/Dt * rDN_DX(j,d);

                // Galerkin pressure term: Div(v) * p
                PDivV = rDN_DX(i,d) * rN[j];

                // Write v * Grad(p) component
                rLHS(Row+d,Col+TDim) += GaussWeight * (G - Sg - PDivV);
                // Use symmetry to write the q * Div(u) component
                rLHS(Col+TDim,Row+d) += GaussWeight * (G + PDivV);

                /* q-p stabilization block */
                // Stabilization: Grad(q) * TauOne * Grad(p)
                L += rDN_DX(i,d) * rDN_DX(j,d);

                /* v-u block */
                // Stabilization: Div(v) * TauTwo * Div(u)
                for (unsigned int e = 0; e < TDim; e++)
                    rLHS(Row+d,Col+e) += GaussWeight*(TauTwo +TauP )*rDN_DX(i,d)*rDN_DX(j,e);
            }

            // Write q-p term
            rLHS(Row+TDim,Col+TDim) += GaussWeight*TauOne*L;
        }

        // RHS terms
        qF = 0.0;
        for (unsigned int d = 0; d < TDim; ++d)
        {
            // v*BodyForce + v * du_ss/dt
            rRHS[Row+d] += GaussWeight * rN[i] * (BodyForce[d] + OldUssTerm[d]);

            // ( a * Grad(v) ) * TauOne * (Density * BodyForce - Projection)
            // vh * TauOne/Dt * f (from vh*d(uss)/dt
            rRHS[Row+d] += GaussWeight * TauOne * (AGradN[i] - Density*rN[i]/Dt ) * ( BodyForce[d] - MomentumProj[d] + OldUssTerm[d] );

            // OSS pressure subscale projection
            rRHS[Row+d] -= GaussWeight * rDN_DX(i,d) * (TauTwo +TauP ) * MassProj;

            // Dynamic term in pressure subscale div(vh) * h^2/(c1*dt) * (-div(uh^n) )
            rRHS[Row+d] -= GaussWeight * rDN_DX(i,d) * TauP * OldResidual;

            // Grad(q) * TauOne * (Density * BodyForce - Projection)
            qF += rDN_DX(i, d) * (BodyForce[d] - MomentumProj[d] + OldUssTerm[d]);
        }
        rRHS[Row + TDim] += GaussWeight * TauOne * qF; // Grad(q) * TauOne * (Density * BodyForce)
    }

    // Viscous contribution (with symmetric gradient 2*nu*{E(u) - 1/3 Tr(E)} )
    // This could potentially be optimized, as it can be integrated exactly using one less integration order when compared to previous terms.
    this->AddViscousTerm(Viscosity,GaussWeight,rDN_DX,rLHS);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DynSS<TDim>::AddMassStabilization(unsigned int GaussIndex,
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
    this->FullConvectiveVelocity(ConvVel,rN,mPredSsVel[GaussIndex]);

    double ElemSize = this->ElementSize();
    //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
    double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    double TauP;
    this->CalculateTau(Density,Viscosity,ConvVel,rProcessInfo,ElemSize,TauOne,TauTwo,TauP);
    const double Dt = rProcessInfo[DELTA_TIME];

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    // Multiplying some quantities by density to have correct units
    //Viscosity *= Density; // Dynamic viscosity
    AGradN *= Density; // Convective term is always multiplied by density

    // Auxiliary variables for matrix looping
    const unsigned int NumNodes = rN.size();
    const unsigned int BlockSize = TDim+1;
    unsigned int Row = 0;
    unsigned int Col = 0;

    // Temporary container
    double K;
    double W = GaussWeight * TauOne * Density; // This density is for the dynamic term in the residual (rho*Du/Dt)

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        Row = i*BlockSize;

        for (unsigned int j = 0; j < NumNodes; j++)
        {
            Col = j*BlockSize;

            // u*grad(v) * TauOne * du/dt
            // v * TauOne/dt * du/dt (from v*d(uss)/dt)
            K = W * (AGradN[i] - rN[i]/Dt) * rN[j];

            for (unsigned int d = 0; d < TDim; d++)
            {
                rMassMatrix(Row+d,Col+d) += K;
                // grad(q) * TauOne * du/dt
                rMassMatrix(Row+TDim,Col+d) += W*rDN_DX(i,d)*rN[j];
            }
        }
    }
}


template<>
void DynSS<2>::DenseSystemSolve(const Matrix &rA, const Vector &rB, Vector &rX) const
{
    Matrix Inv = ZeroMatrix(2,2);
    double Det = 0.0;
    MathUtils<double>::InvertMatrix2(rA,Inv,Det);
    noalias(rX) = prod(Inv,rB);
}


template<>
void DynSS<3>::DenseSystemSolve(const Matrix &rA, const Vector &rB, Vector &rX) const
{
    Matrix Inv = ZeroMatrix(3,3);
    double Det = 0.0;
    MathUtils<double>::InvertMatrix3(rA,Inv,Det);
    noalias(rX) = prod(Inv,rB);
}


template< unsigned int TDim >
void DynSS<TDim>::SubscaleVelocity(unsigned int GaussIndex,
                                   const ShapeFunctionsType &rN,
                                   const ShapeFunctionDerivativesType &rDN_DX,
                                   const ProcessInfo &rProcessInfo,
                                   array_1d<double,3> &rVelocitySubscale)
{
    if (mOldSsVel.size() != 0)
    {
        // Interpolate nodal data on the integration point
        double Density;
        this->EvaluateInPoint(Density,DENSITY,rN);

        array_1d<double,3> ConvVel(3,0.0);
        if (mPredSsVel.size() > 0 )
            this->FullConvectiveVelocity(ConvVel,rN,mPredSsVel[GaussIndex]);
        else
            this->ResolvedConvectiveVelocity(ConvVel,rN);

        double ElemSize = this->ElementSize();
        //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
        double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

        double TauOne;
        double TauTwo;
        double TauP;
        this->CalculateTau(Density,Viscosity,ConvVel,rProcessInfo,ElemSize,TauOne,TauTwo,TauP);

        double Dt = rProcessInfo[DELTA_TIME];

        array_1d<double,3> Residual(3,0.0);

        if (rProcessInfo[OSS_SWITCH] != 1.0)
            this->DASGSMomentumResidual(rN,rDN_DX,ConvVel,Residual);
        else
            this->DOSSMomentumResidual(rN,rDN_DX,ConvVel,Residual);

        rVelocitySubscale = TauOne*(Residual + (Density/Dt)*mOldSsVel[GaussIndex]);
    }
    else
    {
        rVelocitySubscale = array_1d<double,3>(3,0.0);
    }
}

template< unsigned int TDim >
void DynSS<TDim>::SubscalePressure(unsigned int GaussIndex,
                                   const ShapeFunctionsType &rN,
                                   const ShapeFunctionDerivativesType &rDN_DX,
                                   const ProcessInfo &rProcessInfo,
                                   double &rPressureSubscale)
{
    if (mPredSsVel.size() != 0)
    {
        // Interpolate nodal data on the integration point
        double Density;
        this->EvaluateInPoint(Density,DENSITY,rN);

        array_1d<double,3> ConvVel(3,0.0);
        if (mPredSsVel.size() > 0)
            this->FullConvectiveVelocity(ConvVel,rN,mPredSsVel[GaussIndex]);
        else
            this->ResolvedConvectiveVelocity(ConvVel,rN);

        double ElemSize = this->ElementSize();
        //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
        double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

        double TauOne;
        double TauTwo;
        double TauP;
        this->CalculateTau(Density,Viscosity,ConvVel,rProcessInfo,ElemSize,TauOne,TauTwo,TauP);

        // Old mass residual for dynamic pressure subscale
        // Note: Residual is defined as -Div(u)
        double OldResidual = 0.0;
        for (unsigned int a = 0; a < rN.size(); a++)
        {
            const array_1d<double,3>& rOldVel = this->GetGeometry()[a].FastGetSolutionStepValue(VELOCITY,1);
            double OldDivProj = this->GetGeometry()[a].FastGetSolutionStepValue(DIVPROJ,1);
            for (unsigned int d = 0; d < TDim; d++)
                OldResidual -= rDN_DX(a,d)*rOldVel[d] + rN[a]*OldDivProj;
        }

        double Residual = 0.0;

        if (rProcessInfo[OSS_SWITCH] != 1.0)
            this->ASGSMassResidual(GaussIndex,rN,rDN_DX,Residual);
        else
            this->OSSMassResidual(GaussIndex,rN,rDN_DX,Residual);

        rPressureSubscale = (TauTwo+TauP)*Residual - TauP*OldResidual;
        //rPressureSubscale = TauTwo*Residual;
    }
    else
    {
        rPressureSubscale = 0.0;
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< unsigned int TDim >
void DynSS<TDim>::save(Serializer& rSerializer) const
{
    typedef DSS<TDim> _basetype;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, _basetype );
    rSerializer.save("mOldSsVel",mOldSsVel);
}


template< unsigned int TDim >
void DynSS<TDim>::load(Serializer& rSerializer)
{
    typedef DSS<TDim> _basetype;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, _basetype);
    rSerializer.load("mOldSsVel",mOldSsVel);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class DynSS<2>;
template class DynSS<3>;

}
