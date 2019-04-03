//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#include "dss.h"
#include "custom_utilities/turbulence_statistics_container.h"
#include "includes/cfd_variables.h"
#include "../FluidDynamicsApplication/fluid_dynamics_application_variables.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
DSS<TDim>::DSS(IndexType NewId):
    Element(NewId)

{}

template< unsigned int TDim >
DSS<TDim>::DSS(IndexType NewId, const NodesArrayType& ThisNodes):
    Element(NewId,ThisNodes)
{}


template< unsigned int TDim >
DSS<TDim>::DSS(IndexType NewId, GeometryType::Pointer pGeometry):
    Element(NewId,pGeometry)
{}


template< unsigned int TDim >
DSS<TDim>::DSS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    Element(NewId,pGeometry,pProperties)
{}


template< unsigned int TDim >
DSS<TDim>::~DSS()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer DSS<TDim>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DSS(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}


template< unsigned int TDim >
void DSS<TDim>::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                                 VectorType &rRightHandSideVector,
                                                 ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int LocalSize = this->GetGeometry().PointsNumber()*(TDim+1);

    // Resize and intialize output
    if( rLeftHandSideMatrix.size1() != LocalSize )
        rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);

    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize,false);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);
    rRightHandSideVector = ZeroVector(LocalSize);
}


template< unsigned int TDim >
void DSS<TDim>::CalculateLeftHandSide(MatrixType &rLeftHandSideMatrix, ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int LocalSize = this->GetGeometry().PointsNumber()*(TDim+1);

    // Resize and intialize output
    if( rLeftHandSideMatrix.size1() != LocalSize )
        rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);
}


template< unsigned int TDim >
void DSS<TDim>::CalculateRightHandSide(VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int LocalSize = this->GetGeometry().PointsNumber()*(TDim+1);

    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize,false);

    rRightHandSideVector = ZeroVector(LocalSize);
}


template< unsigned int TDim >
void DSS<TDim>::CalculateLocalVelocityContribution(MatrixType &rDampMatrix, VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
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
}


template< unsigned int TDim >
void DSS<TDim>::CalculateMassMatrix(MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo)
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

        /* Note on OSS and full projection: Riccardo says that adding the terms provided by
         * AddMassStabilization (and incluiding their corresponding terms in the projeciton)
         * could help reduce the non-linearity of the coupling between projection and u,p
         * However, leaving them on gives a lot of trouble whith the Bossak scheme:
         * think that we solve F - (1-alpha)*M*u^(n+1) - alpha*M*u^(n) - K(u^(n+1)) = 0
         * so the projection of the dynamic terms should be Pi( (1-alpha)*u^(n+1) - alpha*u^(n) )
         */
        if ( rCurrentProcessInfo[OSS_SWITCH] != 1.0 )
            this->AddMassStabilization(g,GaussWeight,rN,rDN_DX,rCurrentProcessInfo,rMassMatrix);
    }
}


template< unsigned int TDim >
void DSS<TDim>::AddExplicitContribution(ProcessInfo &rCurrentProcessInfo)
{
    /*
     * WARNING: this funciton is assuming that the caller will have taken
     * the necessary steps to avoid race conditions when writing on nodes.
     */
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

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

        array_1d<double,3> Mg = rN[0]*rGeom[0].FastGetSolutionStepValue(MOMENTUM_PROJECTION);
        double Cg = rN[0]*rGeom[0].FastGetSolutionStepValue(MASS_PROJECTION);
        for (unsigned int i = 1; i < NumNodes; i++)
        {
            Mg += rN[i]*rGeom[i].FastGetSolutionStepValue(MOMENTUM_PROJECTION);
            Cg += rN[i]*rGeom[i].FastGetSolutionStepValue(MASS_PROJECTION);
        }

        for (unsigned int j = 0; j < NumNodes; j++)
        {
            rGeom[j].FastGetSolutionStepValue(MOMENTUM_PROJECTION_RHS) -= GaussWeight*rN[j]*Mg;
            rGeom[j].FastGetSolutionStepValue(MASS_PROJECTION_RHS) -= GaussWeight*rN[j]*Cg;
        }
    }
}


template< unsigned int TDim >
void DSS<TDim>::Calculate(const Variable<double> &rVariable,
                          double &rOutput,
                          const ProcessInfo &rCurrentProcessInfo)
{

}


template< unsigned int TDim >
void DSS<TDim>::Calculate(const Variable<array_1d<double,3> > &rVariable,
                          array_1d<double,3> &rOutput,
                          const ProcessInfo &rCurrentProcessInfo)
{
    // Lumped projection terms
    if (rVariable == ADVPROJ)
    {
        // Get Shape function data
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        const unsigned int NumGauss = GaussWeights.size();

        GeometryType& rGeom = this->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();
        VectorType MomentumRHS = ZeroVector(NumNodes*TDim);
        VectorType MassRHS = ZeroVector(NumNodes);
        VectorType NodalArea = ZeroVector(NumNodes);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const double GaussWeight = GaussWeights[g];
            const ShapeFunctionsType& rN = row(ShapeFunctions,g);
            const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

            array_1d<double,3> MomentumRes(3,0.0);
            double MassRes = 0.0;

            this->MomentumProjTerm(g,rN,rDN_DX,MomentumRes);
            this->MassProjTerm(g,rN,rDN_DX,MassRes);

            for (unsigned int i = 0; i < NumNodes; i++)
            {
                double W = GaussWeight*rN[i];
                unsigned int Row = i*TDim;
                for (unsigned int d = 0; d < TDim; d++)
                    MomentumRHS[Row+d] += W*MomentumRes[d];
                NodalArea[i] += W;
                MassRHS[i] += W*MassRes;
            }
        }

        // Add carefully to nodal variables to avoid OpenMP race condition
        unsigned int Row = 0;
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            rGeom[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double,3>& rMomValue = rGeom[i].FastGetSolutionStepValue(ADVPROJ);
            for (unsigned int d = 0; d < TDim; ++d)
                rMomValue[d] += MomentumRHS[Row++];
            rGeom[i].FastGetSolutionStepValue(DIVPROJ) += MassRHS[i];
            rGeom[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
            rGeom[i].UnSetLock(); // Free the node for other threads
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// For TDim == 2
///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
void DSS<2>::EquationIdVector(EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = 3*NumNodes;

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE,ppos).EquationId();
    }
}


template<>
void DSS<2>::GetDofList(DofsVectorType &rElementalDofList, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
     const unsigned int NumNodes = rGeom.PointsNumber();
     const unsigned int LocalSize = 3*NumNodes;

     if (rElementalDofList.size() != LocalSize)
         rElementalDofList.resize(LocalSize);

     unsigned int LocalIndex = 0;

     for (unsigned int i = 0; i < NumNodes; ++i)
     {
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
     }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// For TDim == 3
///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
void DSS<3>::EquationIdVector(EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = 4*NumNodes;

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE,ppos).EquationId();
    }
}


template<>
void DSS<3>::GetDofList(DofsVectorType &rElementalDofList, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
     const unsigned int NumNodes = rGeom.PointsNumber();
     const unsigned int LocalSize = 4*NumNodes;

     if (rElementalDofList.size() != LocalSize)
         rElementalDofList.resize(LocalSize);

     unsigned int LocalIndex = 0;

     for (unsigned int i = 0; i < NumNodes; ++i)
     {
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
     }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS<TDim>::GetFirstDerivativesVector(Vector &rValues, int Step)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes * (TDim+1);

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);

    noalias(rValues) = ZeroVector(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY,Step);
        for (unsigned int d = 0; d < TDim; d++)
            rValues[Index++] = rVel[d];
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(PRESSURE,Step);
    }
}


template< unsigned int TDim >
void DSS<TDim>::GetSecondDerivativesVector(Vector &rValues, int Step)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes * (TDim+1);

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);

    noalias(rValues) = ZeroVector(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rAcc = rGeom[i].FastGetSolutionStepValue(ACCELERATION,Step);
        for (unsigned int d = 0; d < TDim; d++)
            rValues[Index++] = rAcc[d];
        rValues[Index++] = 0.0; // skip pressure Dof
    }
}


template< unsigned int TDim >
GeometryData::IntegrationMethod DSS<TDim>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
//    return GeometryData::IntegrationMethod(this->pGetGeometry()->GetDefaultIntegrationMethod()+1);
}


template< unsigned int TDim >
void DSS<TDim>::GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                            std::vector<array_1d<double, 3 > >& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == SUBSCALE_VELOCITY)
    {
        // Get Shape function data
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        const unsigned int NumGauss = GaussWeights.size();

        rValues.resize(NumGauss);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& rN = row(ShapeFunctions,g);
            const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];
            this->SubscaleVelocity(g,rN,rDN_DX,rCurrentProcessInfo,rValues[g]);
        }
    }
    else if (rVariable == MEAN_VELOCITY)
    {
        unsigned int NumSteps = rCurrentProcessInfo.GetValue(RECORDED_STEPS);
        const DSS<TDim>* const_this = static_cast<const DSS<TDim>*>(this);

        if (const_this->GetValue(TURBULENCE_STATISTICS) != 0)
        {
            const unsigned int NumGauss = const_this->GetGeometry().IntegrationPointsNumber(const_this->GetIntegrationMethod());
            rValues.resize(NumGauss);

            typedef TurbulenceStatisticsContainer TSC;

            for (unsigned int g = 0; g < NumGauss; g++)
            {
                array_1d<double,3>& MeanU = rValues[g];
                MeanU = array_1d<double,3>(3,0.0);
                MeanU[0] = const_this->GetValue(TURBULENCE_STATISTICS)->GetValue(TSC::U,g,NumSteps);
                MeanU[1] = const_this->GetValue(TURBULENCE_STATISTICS)->GetValue(TSC::V,g,NumSteps);
                MeanU[2] = const_this->GetValue(TURBULENCE_STATISTICS)->GetValue(TSC::W,g,NumSteps);

            }
        }
    }
    else if (rVariable == VORTICITY)
    {
        // Get Shape function data
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        const unsigned int NumGauss = GaussWeights.size();

        rValues.resize(NumGauss);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            this->IntegrationPointVorticity(ShapeDerivatives[g],rValues[g]);
        }
    }
}


template< unsigned int TDim >
void DSS<TDim>::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                            std::vector<double>& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == SUBSCALE_PRESSURE)
    {
        // Get Shape function data
        Vector GaussWeights;
        Matrix ShapeFunctions;
        ShapeFunctionDerivativesArrayType ShapeDerivatives;
        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
        const unsigned int NumGauss = GaussWeights.size();

        rValues.resize(NumGauss);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& rN = row(ShapeFunctions,g);
            const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

            this->SubscalePressure(g,rN,rDN_DX,rCurrentProcessInfo,rValues[g]);
        }

    }
    else if (rVariable == TAU) // Data about stabilization coefficients
    {
//        // Get Shape function data
//        Vector GaussWeights;
//        Matrix ShapeFunctions;
//        ShapeFunctionDerivativesArrayType ShapeDerivatives;
//        this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
//        //const unsigned int NumNodes = this->GetGeometry().PointsNumber();
//        const unsigned int NumGauss = GaussWeights.size();

//        rValues.resize(NumGauss);

//        typedef TurbulenceStatisticsContainer TSC;
//        std::vector<double> Values(TSC::NUM_INST_DATA,0.0);

//        for (unsigned int g = 0; g < NumGauss; g++)
//        {
//            const ShapeFunctionsType& rN = row(ShapeFunctions,g);
//            const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

//            double ElemSize = this->ElementSize();
//            double Density;
//            this->EvaluateInPoint(Density,DENSITY,rN);
//            double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rCurrentProcessInfo);
//            array_1d<double,3> ConvVel(3,0.0);
//            this->ResolvedConvectiveVelocity(ConvVel,rN);

//            double TauOne;
//            double TauTwo;
//            this->CalculateStaticTau(Density,Viscosity,ConvVel,rCurrentProcessInfo,TauOne,TauTwo);

//            double MassResidual = 0.0;
//            this->MassProjTerm(g,rN,rDN_DX,MassResidual);

//            double VelNorm = std::sqrt(ConvVel[0]*ConvVel[0] + ConvVel[1]*ConvVel[1] + ConvVel[2]*ConvVel[2]);
//            double Hused = this->Size_MinY();
//            double Hvel = 0.0;

//            if (VelNorm > 1e-12)
//            {
//                array_1d<double,3> Direction = ConvVel/VelNorm;
//                Hvel = this->Size_Projected(Direction);
//            }

//            double Hmin = 0.0;
//            double Hmax = 0.0;
//            this->Size_BoundingBox(Hmin,Hmax);

//            Values[TSC::Time_Label] = rCurrentProcessInfo.GetValue(TIME);
//            Values[TSC::Tau_Momentum] = TauOne;
//            Values[TSC::Tau_Mass] = TauTwo;
//            Values[TSC::Tau_VelNorm] = VelNorm;
//            Values[TSC::H_min] = Hused;
//            Values[TSC::H_vel] = Hvel;
//            Values[TSC::H_min] = Hmin;
//            Values[TSC::H_max] = Hmax;
//            Values[TSC::Div_error] = MassResidual;

//            this->GetValue(TURBULENCE_STATISTICS)->StoreInstantData(Values,g);
//        }
    }
    else if (rVariable == MEAN_KINETIC_ENERGY) // Turbulence statistics
    {
        unsigned int NumSteps = rCurrentProcessInfo.GetValue(RECORDED_STEPS);
        if (NumSteps > 0)
        {
            // Get Shape function data
            Vector GaussWeights;
            Matrix ShapeFunctions;
            ShapeFunctionDerivativesArrayType ShapeDerivatives;
            this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
            const unsigned int NumNodes = this->GetGeometry().PointsNumber();
            const unsigned int NumGauss = GaussWeights.size();

            rValues.resize(NumGauss);

            typedef TurbulenceStatisticsContainer TSC;
            std::vector<double> Values(TSC::NUM_DATA,0.0);

            for (unsigned int g = 0; g < NumGauss; g++)
            {
                const ShapeFunctionsType& rN = row(ShapeFunctions,g);
                const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

                Matrix Gradient = ZeroMatrix(3,3);
                for (unsigned int n = 0; n < NumNodes; n++)
                {
                    const array_1d<double,3>& rU = this->GetGeometry()[n].FastGetSolutionStepValue(VELOCITY);
                    for (unsigned int i = 0; i < TDim; i++)
                    {
                        for (unsigned int j = 0; j < TDim; j++)
                        {
                            Gradient(i,j) += rDN_DX(n,j)*rU[i];
                        }
                    }
                }

                array_1d<double,3> Vel(3,0.0);
                this->EvaluateInPoint(Vel,VELOCITY,rN);
                double P = 0.0;
                this->EvaluateInPoint(P,PRESSURE,rN);

                Values[TSC::U] = Vel[0];
                Values[TSC::V] = Vel[1];
                Values[TSC::W] = Vel[2];

                Values[TSC::P] = P;

                Values[TSC::DU_DX] = Gradient(0,0);
                Values[TSC::DU_DY] = Gradient(0,1);
                Values[TSC::DU_DZ] = Gradient(0,2);

                Values[TSC::DV_DX] = Gradient(1,0);
                Values[TSC::DV_DY] = Gradient(1,1);
                Values[TSC::DV_DZ] = Gradient(1,2);

                Values[TSC::DW_DX] = Gradient(2,0);
                Values[TSC::DW_DY] = Gradient(2,1);
                Values[TSC::DW_DZ] = Gradient(2,2);

                array_1d<double,3> Vss(3,0.0);
                double Pss = 0.0;
                this->SubscaleVelocity(g,rN,rDN_DX,rCurrentProcessInfo,Vss);
                this->SubscalePressure(g,rN,rDN_DX,rCurrentProcessInfo,Pss);

                Values[TSC::Uss] = Vss[0];
                Values[TSC::Vss] = Vss[1];
                Values[TSC::Wss] = Vss[2];

                Values[TSC::Pss] = Pss;

                this->GetValue(TURBULENCE_STATISTICS)->AddStep(Values,g,NumSteps);
            }
        }
    }
    else if (rVariable == MEAN_PRESSURE)
    {
        unsigned int NumSteps = rCurrentProcessInfo.GetValue(RECORDED_STEPS);
        const DSS<TDim>* const_this = static_cast<const DSS<TDim>*>(this);

        if (const_this->GetValue(TURBULENCE_STATISTICS) != 0)
        {
            const unsigned int NumGauss = const_this->GetGeometry().IntegrationPointsNumber(const_this->GetIntegrationMethod());
            rValues.resize(NumGauss);

            typedef TurbulenceStatisticsContainer TSC;

            for (unsigned int g = 0; g < NumGauss; g++)
                rValues[g] = const_this->GetValue(TURBULENCE_STATISTICS)->GetValue(TSC::P,g,NumSteps);
        }
    }
    else if (rVariable == Q_VALUE)
    {
		Vector GaussWeights;
		Matrix ShapeFunctions;
		ShapeFunctionDerivativesArrayType ShapeDerivatives;
		this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
		const unsigned int NumNodes = this->GetGeometry().PointsNumber();
		const unsigned int NumGauss = GaussWeights.size();

		rValues.resize(NumGauss);
		Matrix GradVel;

		// Loop on integration points
		for (unsigned int g = 0; g < NumGauss; g++)
		{
			GradVel = ZeroMatrix(TDim,TDim);
			const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];

			// Compute velocity gradient
			for (unsigned int i=0; i < TDim; ++i)
				for (unsigned int j=0; j < TDim; ++j)
					for (unsigned int iNode=0; iNode < NumNodes; ++iNode)
					{
						array_1d<double,3>& Vel =
							this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
						GradVel(i,j) += Vel[i] * rDN_DX(iNode,j);
					}

			// Compute Q-value
			double qval = 0.0;
			for (unsigned int i=0; i < TDim; ++i)
				for (unsigned int j=0; j < TDim; ++j)
					qval += GradVel(i,j) * GradVel(j,i);

			qval *= -0.5;
			rValues[g] = qval;
		}
	}
	else if (rVariable == VORTICITY_MAGNITUDE)
	{
		Vector GaussWeights;
		Matrix ShapeFunctions;
		ShapeFunctionDerivativesArrayType ShapeDerivatives;
		this->CalculateGeometryData(GaussWeights,ShapeFunctions,ShapeDerivatives);
		const unsigned int NumNodes = this->GetGeometry().PointsNumber();
		const unsigned int NumGauss = GaussWeights.size();

		rValues.resize(NumGauss);

  		// Loop on integration points
		for (unsigned int g = 0; g < NumGauss; g++)
		{
			const ShapeFunctionDerivativesType& rDN_DX = ShapeDerivatives[g];
			array_1d<double,3> Vorticity(3,0.0);

			if(TDim == 2)
			{
				for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
				{
					array_1d<double,3>& Vel =
						this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
					Vorticity[2] += Vel[1] * rDN_DX(iNode,0) - Vel[0] * rDN_DX(iNode,1);
				}
			}
			else
			{
				for (unsigned int iNode = 0; iNode < this->GetGeometry().size(); iNode++)
				{
					array_1d<double,3>& Vel =
						this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
					Vorticity[0] += Vel[2] * rDN_DX(iNode,1) - Vel[1] * rDN_DX(iNode,2);
					Vorticity[1] += Vel[0] * rDN_DX(iNode,2) - Vel[2] * rDN_DX(iNode,0);
					Vorticity[2] += Vel[1] * rDN_DX(iNode,0) - Vel[0] * rDN_DX(iNode,1);
				}
			}

			rValues[g] = sqrt(Vorticity[0] * Vorticity[0] + Vorticity[1] * Vorticity[1]
					+ Vorticity[2] * Vorticity[2]);
		}
	}
}


template< unsigned int TDim >
void DSS<TDim>::GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                            std::vector<array_1d<double, 6 > >& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
{

}


template< unsigned int TDim >
void DSS<TDim>::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                            std::vector<Vector>& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
{

}


template< unsigned int TDim >
void DSS<TDim>::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                            std::vector<Matrix>& rValues,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == VELOCITY_COVARIANCES)
    {
        unsigned int NumSteps = rCurrentProcessInfo.GetValue(RECORDED_STEPS);
        const DSS<TDim>* const_this = static_cast<const DSS<TDim>*>(this);

        if (const_this->GetValue(TURBULENCE_STATISTICS) != 0)
        {
            const unsigned int NumGauss = const_this->GetGeometry().IntegrationPointsNumber(const_this->GetIntegrationMethod());
            rValues.resize(NumGauss);

            typedef TurbulenceStatisticsContainer TSC;

            for (unsigned int g = 0; g < NumGauss; g++)
            {
                Matrix& Covariances = rValues[g];
                Covariances = ZeroMatrix(3,3);
                Covariances(0,0) = const_this->GetValue(TURBULENCE_STATISTICS)->GetValue(TSC::UU,g,NumSteps);
                Covariances(0,1) = const_this->GetValue(TURBULENCE_STATISTICS)->GetValue(TSC::UV,g,NumSteps);
                Covariances(0,2) = const_this->GetValue(TURBULENCE_STATISTICS)->GetValue(TSC::UW,g,NumSteps);
                Covariances(1,0) = Covariances(1,0);
                Covariances(1,1) = const_this->GetValue(TURBULENCE_STATISTICS)->GetValue(TSC::VV,g,NumSteps);
                Covariances(1,2) = const_this->GetValue(TURBULENCE_STATISTICS)->GetValue(TSC::VW,g,NumSteps);
                Covariances(2,0) = Covariances(0,2);
                Covariances(2,1) = Covariances(1,2);
                Covariances(2,2) = const_this->GetValue(TURBULENCE_STATISTICS)->GetValue(TSC::WW,g,NumSteps);
            }
        }
    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< unsigned int TDim >
int DSS<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);

    // Check that required variables are registered

    // Nodal data
    if(VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check that required variables have been correctly registered.","");
    if(PRESSURE.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check that required variables have been correctly registered.","");
    if(BODY_FORCE.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"BODY_FORCE Key is 0. Check that the application was correctly registered.","");
    if(DENSITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check that the application was correctly registered.","");
    if(VISCOSITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"VISCOSITY Key is 0. Check that the application was correctly registered.","");
    if(MESH_VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check that the application was correctly registered.","");

    // Variables used in computing projections
    if(ACCELERATION.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION Key is 0. Check that the application was correctly registered.","");
    if(NODAL_AREA.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"NODAL_AREA Key is 0. Check that the application was correctly registered.","");
    if(ADVPROJ.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"ADVPROJ Key is 0. Check that the application was correctly registered.","");
    if(DIVPROJ.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"DIVPROJ Key is 0. Check that the application was correctly registered.","");

    // Output variables (for Calculate() functions)
    if(SUBSCALE_VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"SUBSCALE_VELOCITY Key is 0. Check that the application was correctly registered.","");
    if(SUBSCALE_PRESSURE.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"SUBSCALE_PRESSURE Key is 0. Check that the application was correctly registered.","");

    // Elemental data
    if(C_SMAGORINSKY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"C_SMAGORINSKY Key is 0. Check that the application was correctly registered.","");

    // Process Info data
    // none at the moment, consider OSS_SWITCH


    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
    {
        // Check that required nodal variables are included in SolutionStepData
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing BODY_FORCE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(VISCOSITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());

        if(this->GetGeometry()[i].SolutionStepsDataHas(ACCELERATION) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing ACCELERATION variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(NODAL_AREA) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing NODAL_AREA variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(ADVPROJ) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing ADVPROJ variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(DIVPROJ) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing DIVPROJ variable on solution step data for node ",this->GetGeometry()[i].Id());

        // Check that required dofs exist
        if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
           this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
           this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if (this->GetGeometry().WorkingSpaceDimension() == 2)
    {
        for (unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if (this->GetGeometry()[i].Z() != 0.0)
                KRATOS_THROW_ERROR(std::invalid_argument,"Node with non-zero Z coordinate found. Id: ",this->GetGeometry()[i].Id());
        }
    }

    return out;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output


template< unsigned int TDim >
std::string DSS<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "DSS #" << Id();
    return buffer.str();
}


template< unsigned int TDim >
void DSS<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "DSS" << TDim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS<TDim>::CalculateGeometryData(Vector &rGaussWeights,
                                      Matrix &rNContainer,
                                      ShapeFunctionDerivativesArrayType &rDN_DX)
{
    const GeometryData::IntegrationMethod IntMethod = this->GetIntegrationMethod();
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
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

/*
template< unsigned int TDim >
double DSS<TDim>::ElementSize() const
{
//    double Area = this->GetGeometry().Area();
//    return 1.128379167 * sqrt(Area); //Diameter of circumference of given Area

    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

//    // calculate minimum element length (used in stabilization Tau)
//    array_1d<double,3> Edge(3,0.0);
//    Edge = rGeom[1].Coordinates() - rGeom[0].Coordinates();
//    double ElemSize = Edge[0]*Edge[0];
//    for (SizeType d = 1; d < TDim; d++)
//        ElemSize += Edge[d]*Edge[d];

//    for (SizeType i = 2; i < NumNodes; i++)
//        for(SizeType j = 0; j < i; j++)
//        {
//            Edge = rGeom[i].Coordinates() - rGeom[j].Coordinates();
//            double Length = Edge[0]*Edge[0];
//            for (SizeType d = 1; d < TDim; d++)
//                Length += Edge[d]*Edge[d];
//            if (Length < ElemSize) ElemSize = Length;
//        }
//    return sqrt(ElemSize);

    double zmax = rGeom[0].Y();
    double zmin = zmax;
//    double h = sqrt(ElemSize);
    for (unsigned int i = 1; i < NumNodes; i++)
    {
        double z = rGeom[i].Y();
        if (z < zmin) zmin = z;
        else if (z > zmax) zmax = z;
    }

    return zmax-zmin;
}*/

/** Calculate characteristic element length.
 * @return Minimum element height
 */
template< unsigned int TDim >
double DSS<TDim>::ElementSize()
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
        KRATOS_THROW_ERROR(std::invalid_argument,"DSS::ElementSize not implemented for this geometry type","");
        return 0;
    }
    }
    KRATOS_CATCH("");
}

template< unsigned int TDim >
double DSS<TDim>::ElementSize(const array_1d<double,3> &rVel,
                              const ShapeFunctionDerivativesType &rDN_DX)
{
    double Num = 2.0*std::sqrt(rVel[0]*rVel[0]+rVel[1]*rVel[1]+rVel[2]*rVel[2]);

    if (Num > 1e-10)
    {
        double Den = 0.0;
        double Tmp = 0.0;
        for (unsigned int i = 0; i < rDN_DX.size1(); i++)
        {
            Tmp = rVel[0]*rDN_DX(i,0);
            for (unsigned int d = 1; d < TDim; d++)
                Tmp += rVel[d]*rDN_DX(i,d);
            Den += std::fabs(Tmp);
        }

        return Num/Den;
    }
    else
        return this->ElementSize();
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS<TDim>::CalculateStaticTau(double Density,
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
    TauTwo = Density * (KinematicVisc + c2 * VelNorm * ElemSize / c1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS<TDim>::ASGSMomentumResidual(double GaussIndex,
                                     const ShapeFunctionsType &rN,
                                     const ShapeFunctionDerivativesType &rDN_DX,
                                     array_1d<double,3> &rMomentumRes)
{
    const GeometryType rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

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
void DSS<TDim>::ASGSMassResidual(double GaussIndex,
                                 const ShapeFunctionsType &rN,
                                 const ShapeFunctionDerivativesType &rDN_DX,
                                 double &rMomentumRes)
{
    this->MassProjTerm(GaussIndex,rN,rDN_DX,rMomentumRes);
}


template< unsigned int TDim >
void DSS<TDim>::OSSMomentumResidual(double GaussIndex,
                                    const ShapeFunctionsType &rN,
                                    const ShapeFunctionDerivativesType &rDN_DX,
                                    array_1d<double,3> &rMomentumRes)
{
    this->MomentumProjTerm(GaussIndex,rN,rDN_DX,rMomentumRes);

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
void DSS<TDim>::OSSMassResidual(double GaussIndex,
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
void DSS<TDim>::MomentumProjTerm(double GaussIndex,
                                 const ShapeFunctionsType &rN,
                                 const ShapeFunctionDerivativesType &rDN_DX,
                                 array_1d<double,3> &rMomentumRHS)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    Vector AGradN;
    this->ConvectionOperator(AGradN,ConvVel,rDN_DX);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rBodyForce = rGeom[i].FastGetSolutionStepValue(BODY_FORCE);
        //const array_1d<double,3>& rAcc = rGeom[i].FastGetSolutionStepValue(ACCELERATION);
        //const array_1d<double,3>& rAccOld = rGeom[i].FastGetSolutionStepValue(ACCELERATION,1);
        //array_1d<double,3> BossakAcc = 1.3*rAcc - 0.3*rAccOld;
        const array_1d<double,3>& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        const double Press = rGeom[i].FastGetSolutionStepValue(PRESSURE);

        for (unsigned int d = 0; d < TDim; d++)
        {
            rMomentumRHS[d] += Density * ( rN[i]*(rBodyForce[d] /*- BossakAcc[d]*/ /*- rAcc[d]*/) - AGradN[i]*rVel[d]) - rDN_DX(i,d)*Press;
        }
    }
}


template< unsigned int TDim >
void DSS<TDim>::MassProjTerm(double GaussIndex,
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
double DSS<TDim>::EffectiveViscosity(const ShapeFunctionsType &rN,
                                     const ShapeFunctionDerivativesType &rDN_DX,
                                     double ElemSize,
                                     const ProcessInfo &rCurrentProcessInfo)
{
    const DSS* const_this = static_cast<const DSS*>(this);
    double Csmag = const_this->GetValue(C_SMAGORINSKY);

    double KinViscosity = 0.0;
    this->EvaluateInPoint(KinViscosity,VISCOSITY,rN);

    if (Csmag != 0.0 )
    {
        const unsigned int NumNodes = this->GetGeometry().PointsNumber();

        // Calculate Symetric gradient
        MatrixType S = ZeroMatrix(TDim,TDim);
        for (unsigned int n = 0; n < NumNodes; ++n)
        {
            const array_1d<double,3>& rVel = this->GetGeometry()[n].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int i = 0; i < TDim; ++i)
                for (unsigned int j = 0; j < TDim; ++j)
                    S(i,j) += 0.5 * ( rDN_DX(n,j) * rVel[i] + rDN_DX(n,i) * rVel[j] );
        }

        // Norm of symetric gradient
        double NormS = 0.0;
        for (unsigned int i = 0; i < TDim; ++i)
            for (unsigned int j = 0; j < TDim; ++j)
                NormS += S(i,j) * S(i,j);
        NormS = sqrt(2.0*NormS);

        // Nu_sgs = (Csmag * Delta)^2 * (2*Sij*Sij)^(1/2)
        KinViscosity += Csmag * Csmag * ElemSize * ElemSize * NormS;
    }

    return KinViscosity;
}


template< unsigned int TDim >
void DSS<TDim>::ResolvedConvectiveVelocity(array_1d<double,3> &rConvVel, const ShapeFunctionsType &rN)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    array_1d<double,3> NodeVel = rGeom[0].FastGetSolutionStepValue(VELOCITY);
    NodeVel -= rGeom[0].FastGetSolutionStepValue(MESH_VELOCITY);
    rConvVel = rN[0] * NodeVel;

    for (unsigned int i = 1; i < NumNodes; i++)
    {
        NodeVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        NodeVel -= rGeom[i].FastGetSolutionStepValue(MESH_VELOCITY);
        rConvVel += rN[i] * NodeVel;
    }
}


template< unsigned int TDim >
void DSS<TDim>::FullConvectiveVelocity(array_1d<double,3> &rConvVel,
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

template< unsigned int TDim >
void DSS<TDim>::ConvectionOperator(Vector &rResult,
                                   const array_1d<double,3> &rConvVel,
                                   const ShapeFunctionDerivativesType &DN_DX)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();

    if(rResult.size() != NumNodes) rResult.resize(NumNodes,false);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        rResult[i] = rConvVel[0]*DN_DX(i,0);
        for(unsigned int k = 1; k < TDim; k++)
            rResult[i] += rConvVel[k]*DN_DX(i,k);
    }
}



///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points
///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS<TDim>::AddSystemTerms(unsigned int GaussIndex,
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
            //K = 0.5*(rN[i]*AGradN[j] - AGradN[i]*rN[j]); // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
            K = rN[i]*AGradN[j];
            K += AGradN[i]*TauOne*(AGradN[j]); // Stabilization: u*grad(v) * TauOne * u*grad(u)
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
                G = TauOne * AGradN[i] * rDN_DX(j,d); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                PDivV = rDN_DX(i,d) * rN[j]; // Div(v) * p

                // Write v * Grad(p) component
                rLHS(Row+d,Col+TDim) += GaussWeight * (G - PDivV);
                // Use symmetry to write the q * Div(u) component
                rLHS(Col+TDim,Row+d) += GaussWeight * (G + PDivV);

                // q-p stabilization block
                L += rDN_DX(i,d) * rDN_DX(j,d); // Stabilization: Grad(q) * TauOne * Grad(p)

                for (unsigned int e = 0; e < TDim; e++) // Stabilization: Div(v) * TauTwo * Div(u)
                    rLHS(Row+d,Col+e) += GaussWeight*TauTwo*rDN_DX(i,d)*rDN_DX(j,e);
            }

            // Write q-p term
            rLHS(Row+TDim,Col+TDim) += GaussWeight*TauOne*L;
        }

        // RHS terms
        qF = 0.0;
        for (unsigned int d = 0; d < TDim; ++d)
        {
            rRHS[Row+d] += GaussWeight * rN[i] * BodyForce[d]; // v*BodyForce
            rRHS[Row+d] += GaussWeight * TauOne * AGradN[i] * ( BodyForce[d] - MomentumProj[d]); // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
            rRHS[Row+d] -= GaussWeight * TauTwo * rDN_DX(i,d) * MassProj;
            qF += rDN_DX(i, d) * (BodyForce[d] - MomentumProj[d]);
        }
        rRHS[Row + TDim] += GaussWeight * TauOne * qF; // Grad(q) * TauOne * (Density * BodyForce)
    }

    // Viscous contribution (with symmetric gradient 2*nu*{E(u) - 1/3 Tr(E)} )
    // This could potentially be optimized, as it can be integrated exactly using one less integration order when compared to previous terms.
    this->AddViscousTerm(Viscosity,GaussWeight,rDN_DX,rLHS);

    // Modulated gradient diffusion
    //this->ModulatedGradientDiffusion(rLHS,rDN_DX,GaussWeight);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS<TDim>::AddMassTerms(double GaussWeight, const ShapeFunctionsType &rN, MatrixType &rMassMatrix)
{
    const unsigned int NumNodes = rN.size();
    const unsigned int BlockSize = TDim+1;

    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    unsigned int Row = 0;
    unsigned int Col = 0;

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        Row = i*BlockSize;
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            Col = j*BlockSize;
            const double Mij = GaussWeight * Density * rN[i] * rN[j];
            for (unsigned int d = 0; d < TDim; d++)
                rMassMatrix(Row+d,Col+d) += Mij;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS<TDim>::AddMassStabilization(unsigned int GaussIndex,
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
    //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
    double Viscosity = EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

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

            K = W * AGradN[i] * rN[j];

            for (unsigned int d = 0; d < TDim; d++)
            {
                rMassMatrix(Row+d,Col+d) += K;
                rMassMatrix(Row+TDim,Col+d) += W*rDN_DX(i,d)*rN[j];
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// For TDim == 2
///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
void DSS<2>::AddViscousTerm(double DynamicViscosity,
                            double GaussWeight,
                            const ShapeFunctionDerivativesType &rDN_DX,
                            MatrixType &rLHS)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = 3;
    double Weight = GaussWeight * DynamicViscosity;

    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    // The implementation is equivalent to the following matrix product, made traceless
    /*
    MatrixType B = ZeroMatrix(3,NumNodes*BlockSize);

    for (int a = 0; a < NumNodes; a++)
    {
        int col = BlockSize*a;

        B(0,col)     = rDN_DX(a,0);
        B(0,col+1)   = 0.0;
        B(1,col)     = 0.0;
        B(1,col+1)   = rDN_DX(a,1);
        B(2,col)     = rDN_DX(a,1);
        B(2,col+1)   = rDN_DX(a,0);
    }

    MatrixType C = ZeroMatrix(3,3);
    C(0,0) = 2.0*DynamicViscosity;
    C(1,1) = 2.0*DynamicViscosity;
    C(2,2) = 1.0*DynamicViscosity;

    MatrixType Temp = prod(C,B);
    rLHS += GaussWeight * prod(trans(B),Temp);
    */

    unsigned int Row(0),Col(0);

    for (unsigned int a = 0; a < NumNodes; ++a)
    {
        Row = a*BlockSize;
        for (unsigned int b = 0; b < NumNodes; ++b)
        {
            Col = b*BlockSize;

            // First Row
            rLHS(Row,Col) += Weight * ( FourThirds * rDN_DX(a,0) * rDN_DX(b,0) + rDN_DX(a,1) * rDN_DX(b,1) );
            rLHS(Row,Col+1) += Weight * ( nTwoThirds * rDN_DX(a,0) * rDN_DX(b,1) + rDN_DX(a,1) * rDN_DX(b,0) );

            // Second Row
            rLHS(Row+1,Col) += Weight * ( nTwoThirds * rDN_DX(a,1) * rDN_DX(b,0) + rDN_DX(a,0) * rDN_DX(b,1) );
            rLHS(Row+1,Col+1) += Weight * ( FourThirds * rDN_DX(a,1) * rDN_DX(b,1) + rDN_DX(a,0) * rDN_DX(b,0) );
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// For TDim == 3
///////////////////////////////////////////////////////////////////////////////////////////////////

template<>
void DSS<3>::AddViscousTerm(double DynamicViscosity,
                            double GaussWeight,
                            const ShapeFunctionDerivativesType &rDN_DX,
                            MatrixType &rLHS)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = 4;
    double Weight = GaussWeight * DynamicViscosity;

    const double OneThird = 1.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    unsigned int Row(0),Col(0);

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        Row = i*BlockSize;
        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            Col = j*BlockSize;
            // (dN_i/dx_k dN_j/dx_k)
            const double Diag =  rDN_DX(i,0) * rDN_DX(j,0) + rDN_DX(i,1) * rDN_DX(j,1) + rDN_DX(i,2) * rDN_DX(j,2);

            // First Row
            rLHS(Row,Col) += Weight * ( OneThird * rDN_DX(i,0) * rDN_DX(j,0) + Diag );
            rLHS(Row,Col+1) += Weight * ( nTwoThirds * rDN_DX(i,0) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,0) );
            rLHS(Row,Col+2) += Weight * ( nTwoThirds * rDN_DX(i,0) * rDN_DX(j,2) + rDN_DX(i,2) * rDN_DX(j,0) );

            // Second Row
            rLHS(Row+1,Col) += Weight * ( nTwoThirds * rDN_DX(i,1) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,1) );
            rLHS(Row+1,Col+1) += Weight * ( OneThird * rDN_DX(i,1) * rDN_DX(j,1) + Diag );
            rLHS(Row+1,Col+2) += Weight * ( nTwoThirds * rDN_DX(i,1) * rDN_DX(j,2) + rDN_DX(i,2) * rDN_DX(j,1) );

            // Third Row
            rLHS(Row+2,Col) += Weight * ( nTwoThirds * rDN_DX(i,2) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,2) );
            rLHS(Row+2,Col+1) += Weight * ( nTwoThirds * rDN_DX(i,2) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,2) );
            rLHS(Row+2,Col+2) += Weight * ( OneThird * rDN_DX(i,2) * rDN_DX(j,2) + Diag );
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Modulated gradient diffusion
///////////////////////////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void DSS<TDim>::ModulatedGradientDiffusion(MatrixType &rDampMatrix, const ShapeFunctionDerivativesType &rDN_DX, const double Weight)
{

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Velocity gradient
    MatrixType GradU = ZeroMatrix(TDim,TDim);
    for (unsigned int n = 0; n < NumNodes; n++)
    {
        const array_1d<double,3>& rVel = this->GetGeometry()[n].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
                GradU(i,j) += rDN_DX(n,j)*rVel[i];
    }

    // Element lengths
    array_1d<double,3> Delta(3,0.0);
    Delta[0] = fabs(rGeom[NumNodes-1].X()-rGeom[0].X());
    Delta[1] = fabs(rGeom[NumNodes-1].Y()-rGeom[0].Y());
    Delta[2] = fabs(rGeom[NumNodes-1].Z()-rGeom[0].Z());

    for (unsigned int n = 1; n < NumNodes; n++)
    {
        double hx = fabs(rGeom[n].X()-rGeom[n-1].X());
        if (hx > Delta[0]) Delta[0] = hx;
        double hy = fabs(rGeom[n].Y()-rGeom[n-1].Y());
        if (hy > Delta[1]) Delta[1] = hy;
        double hz = fabs(rGeom[n].Z()-rGeom[n-1].Z());
        if (hz > Delta[2]) Delta[2] = hz;
    }

    double AvgDeltaSq = Delta[0];
    for (unsigned int d = 1; d < TDim; d++)
        AvgDeltaSq *= Delta[d];
    AvgDeltaSq = std::pow(AvgDeltaSq,2./TDim);

    Delta[0] = Delta[0]*Delta[0]/12.0;
    Delta[1] = Delta[1]*Delta[1]/12.0;
    Delta[2] = Delta[2]*Delta[2]/12.0;

    // Gij
    MatrixType G = ZeroMatrix(TDim,TDim);
    for (unsigned int i = 0; i < TDim; i++)
        for (unsigned int j = 0; j < TDim; j++)
            for (unsigned int d = 0; d < TDim; d++)
                G(i,j) += Delta[d]*GradU(i,d)*GradU(j,d);

    // Gij:Sij
    double GijSij = 0.0;
    for (unsigned int i = 0; i < TDim; i++)
        for (unsigned int j = 0; j < TDim; j++)
            GijSij += 0.5*G(i,j)*( GradU(i,j) + GradU(j,i) );

    if (GijSij < 0.0) // Otherwise model term is clipped
    {
        // Gkk
        double Gkk = G(0,0);
        for (unsigned int d = 1; d < TDim; d++)
            Gkk += G(d,d);

        // C_epsilon
        const double Ce = 1.0;

        // ksgs
        double ksgs = -4*AvgDeltaSq*GijSij/(Ce*Ce*Gkk);

        // Assembly of model term
        unsigned int RowIndex = 0;
        unsigned int ColIndex = 0;

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            for (unsigned int j = 0; j < NumNodes; j++)
            {
                for (unsigned int d = 0; d < TDim; d++)
                {
                    double Aux = Delta[0] * G(d,0)*rDN_DX(j,0);
                    for (unsigned int k = 1; k < TDim; k++)
                        Aux += Delta[k] * G(d,k)*rDN_DX(j,k);
                    rDampMatrix(RowIndex+d,ColIndex+d) += Weight * 2.0*ksgs * rDN_DX(i,d) * Aux;
                }

                ColIndex += TDim;
            }
            RowIndex += TDim;
            ColIndex = 0;
        }
    }
}

template< unsigned int TDim >
void DSS<TDim>::SubscaleVelocity(unsigned int GaussIndex,
                                 const ShapeFunctionsType &rN,
                                 const ShapeFunctionDerivativesType &rDN_DX,
                                 const ProcessInfo &rProcessInfo,
                                 array_1d<double,3> &rVelocitySubscale)
{
    // Interpolate nodal data on the integration point
    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    double ElemSize = this->ElementSize();
    //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
    double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

    array_1d<double,3> Residual(3,0.0);

    if (rProcessInfo[OSS_SWITCH] != 1.0)
        this->ASGSMomentumResidual(GaussIndex,rN,rDN_DX,Residual);
    else
        this->OSSMomentumResidual(GaussIndex,rN,rDN_DX,Residual);

    rVelocitySubscale = TauOne*Residual;
}

template< unsigned int TDim >
void DSS<TDim>::SubscalePressure(unsigned int GaussIndex,
                                 const ShapeFunctionsType &rN,
                                 const ShapeFunctionDerivativesType &rDN_DX,
                                 const ProcessInfo &rProcessInfo,
                                 double &rPressureSubscale)
{
    // Interpolate nodal data on the integration point
    double Density;
    this->EvaluateInPoint(Density,DENSITY,rN);

    array_1d<double,3> ConvVel(3,0.0);
    this->ResolvedConvectiveVelocity(ConvVel,rN);

    //double ElemSize = this->ElementSize(ConvVel,rDN_DX);
    double ElemSize = this->ElementSize();
    double Viscosity = this->EffectiveViscosity(rN,rDN_DX,ElemSize,rProcessInfo);

    double TauOne;
    double TauTwo;
    this->CalculateStaticTau(Density,Viscosity,ConvVel,ElemSize,rProcessInfo,TauOne,TauTwo);

    double Residual = 0.0;

    if (rProcessInfo[OSS_SWITCH] != 1.0)
        this->ASGSMassResidual(GaussIndex,rN,rDN_DX,Residual);
    else
        this->OSSMassResidual(GaussIndex,rN,rDN_DX,Residual);

    rPressureSubscale = TauTwo*Residual;
}


template<>
void DSS<2>::IntegrationPointVorticity(const ShapeFunctionDerivativesType &rDN_DX,array_1d<double,3>& rVorticity) const
{
    rVorticity = array_1d<double,3>(3,0.0);
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        const array_1d<double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        rVorticity[2] += rDN_DX(i,0)*rVelocity[1] - rDN_DX(i,1)*rVelocity[0];
    }
}


template<>
void DSS<3>::IntegrationPointVorticity(const ShapeFunctionDerivativesType &rDN_DX,array_1d<double,3>& rVorticity) const
{
    rVorticity = array_1d<double,3>(3,0.0);
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();

    for (unsigned int i = 0; i < NumNodes; ++i)
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

template< unsigned int TDim >
void DSS<TDim>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
}


template< unsigned int TDim >
void DSS<TDim>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class DSS<2>;
template class DSS<3>;

}
