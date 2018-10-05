#include "dynamic_vms.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

// public DynamicVMS methods **************************************************

template< unsigned int TDim >
DynamicVMS<TDim>::DynamicVMS(IndexType NewId, const NodesArrayType &ThisNodes):
    Element(NewId,ThisNodes),
    mIntegrationMethod(GeometryData::GI_GAUSS_1)
{
    unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->mIntegrationMethod);
    mSubscaleVel.resize(NumGauss,array_1d<double,3>(3,0.0));
    mOldSubscaleVel.resize(NumGauss,array_1d<double,3>(3,0.0));
    mIterCount.resize(NumGauss,0);
    this->CalculateGeometryData();
}


template< unsigned int TDim >
DynamicVMS<TDim>::DynamicVMS(IndexType NewId, GeometryType::Pointer pGeometry):
    Element(NewId,pGeometry),
    mIntegrationMethod(GeometryData::GI_GAUSS_1)
{
    unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->mIntegrationMethod);
    mSubscaleVel.resize(NumGauss,array_1d<double,3>(3,0.0));
    mOldSubscaleVel.resize(NumGauss,array_1d<double,3>(3,0.0));
    mIterCount.resize(NumGauss,0);
    this->CalculateGeometryData();
}

template< unsigned int TDim >
DynamicVMS<TDim>::DynamicVMS(IndexType NewId, GeometryType::Pointer pGeometry, const GeometryData::IntegrationMethod ThisIntegrationMethod):
    Element(NewId,pGeometry),
    mIntegrationMethod(ThisIntegrationMethod)
{
    unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->mIntegrationMethod);
    mSubscaleVel.resize(NumGauss,array_1d<double,3>(3,0.0));
    mOldSubscaleVel.resize(NumGauss,array_1d<double,3>(3,0.0));
    mIterCount.resize(NumGauss,0);
    this->CalculateGeometryData();
}

template< unsigned int TDim >
DynamicVMS<TDim>::DynamicVMS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    Element(NewId,pGeometry,pProperties),
    mIntegrationMethod(GeometryData::GI_GAUSS_1)
{
    unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->mIntegrationMethod);
    mSubscaleVel.resize(NumGauss,array_1d<double,3>(3,0.0));
    mOldSubscaleVel.resize(NumGauss,array_1d<double,3>(3,0.0));
    mIterCount.resize(NumGauss,0);
    this->CalculateGeometryData();
}


template< unsigned int TDim >
DynamicVMS<TDim>::DynamicVMS(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties, const GeometryData::IntegrationMethod ThisIntegrationMethod):
    Element(NewId,pGeometry,pProperties),
    mIntegrationMethod(ThisIntegrationMethod)
{
    unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->mIntegrationMethod);
    mSubscaleVel.resize(NumGauss,array_1d<double,3>(3,0.0));
    mOldSubscaleVel.resize(NumGauss,array_1d<double,3>(3,0.0));
    mIterCount.resize(NumGauss,0);
    this->CalculateGeometryData();
}

template< unsigned int TDim >
DynamicVMS<TDim>::~DynamicVMS()
{}


template< unsigned int TDim >
Element::Pointer DynamicVMS<TDim>::Create(IndexType NewId, const NodesArrayType &ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<DynamicVMS>(NewId,this->GetGeometry().Create(ThisNodes),pProperties,this->mIntegrationMethod);
}

template< unsigned int TDim >
Element::Pointer DynamicVMS<TDim>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<DynamicVMS>(NewId,pGeom,pProperties,this->mIntegrationMethod);
}


template< unsigned int TDim >
void DynamicVMS<TDim>::Initialize()
{}

template< unsigned int TDim >
void DynamicVMS<TDim>::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    this->CalculateGeometryData();

    // Advance the subscale values in time
    mOldSubscaleVel = mSubscaleVel;
}


template< unsigned int TDim >
void DynamicVMS<TDim>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{
//    this->LinearUpdateSubscale(rCurrentProcessInfo);
    this->UpdateSubscale(rCurrentProcessInfo);
//    for (unsigned int g = 0; g < this->GetGeometry().IntegrationPointsNumber(this->mIntegrationMethod); g++)
//        mSubscaleVel[g] = ZeroVector(3);
}


template< unsigned int TDim >
void DynamicVMS<TDim>::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                            VectorType &rRightHandSideVector,
                                            ProcessInfo &rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes * (TDim + 1);

    if(rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);

    if(rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize,false);

    noalias(rRightHandSideVector) = ZeroVector(LocalSize);
}


template< unsigned int TDim >
void DynamicVMS<TDim>::CalculateRightHandSide(VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes * (TDim + 1);

    if(rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize,false);

    noalias(rRightHandSideVector) = ZeroVector(LocalSize);
}


template< unsigned int TDim >
void DynamicVMS<TDim>::CalculateLocalVelocityContribution(MatrixType &rDampingMatrix,
                                                          VectorType &rRightHandSideVector,
                                                          ProcessInfo &rCurrentProcessInfo)
{
    if (rCurrentProcessInfo[OSS_SWITCH]==1.0)
        this->CalculateOSSVelocityContribution(rDampingMatrix,rRightHandSideVector,rCurrentProcessInfo);
    else
        this->CalculateASGSVelocityContribution(rDampingMatrix,rRightHandSideVector,rCurrentProcessInfo);
}


template< unsigned int TDim >
void DynamicVMS<TDim>::CalculateMassMatrix(MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int BlockSize = TDim + 1;
    const unsigned int LocalSize = NumNodes * BlockSize;

    if(rMassMatrix.size1() != LocalSize)
        rMassMatrix.resize(LocalSize,LocalSize,false);

    noalias(rMassMatrix) = ZeroMatrix(LocalSize,LocalSize);

    if (this->mIntegrationMethod == GeometryData::GI_GAUSS_1)
        this->LumpedMassMatrix(rMassMatrix);
    else
        this->ConsistentMassMatrix(rMassMatrix);

    // Additional ASGS terms
    if(rCurrentProcessInfo[OSS_SWITCH]!=1.0)
    {
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(this->mIntegrationMethod);
        const unsigned int NumGauss = IntegrationPoints.size();
        const Matrix& NContainer = rGeom.ShapeFunctionsValues(this->mIntegrationMethod);

        const double Dt = rCurrentProcessInfo[DELTA_TIME];

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& N = row(NContainer,g);
            const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();

            // Evaluate problem parameters at integration point
            double Density = 0.0;
            double Viscosity = 0.0;
            array_1d<double,3> ConvVel = ZeroVector(3);
            array_1d<double,3> BodyForce = ZeroVector(3);
            Vector Convection = ZeroVector(NumNodes);

            this->EvaluateInPoint(Density,DENSITY,N);
            this->EvaluateViscosity(Viscosity,N);
            this->EvaluateConvVelocity(ConvVel,mSubscaleVel[g],N);
            this->EvaluateInPoint(BodyForce,BODY_FORCE,N);
            this->ConvectionOperator(Convection,ConvVel);
            Convection *= Density;

            // Calculate stablitization parameters
            double ConvVelNorm = ConvVel[0] * ConvVel[0];
            for (unsigned int d = 1; d < TDim; d++)
                ConvVelNorm += ConvVel[d] * ConvVel[d];
            ConvVelNorm = sqrt(ConvVelNorm);

            double Tau_t = this->TauTime(Density,Viscosity,ConvVelNorm,Dt);

            unsigned int RowIndex = 0;
            unsigned int ColIndex = 0;

            const double GaussMass = GaussWeight * Density;

            for (unsigned int i = 0; i < NumNodes; i++)
            {
                double Conv_i = Convection[i] * Tau_t;
                for (unsigned int j = 0; j < NumNodes; j++)
                {
                    double Mij = GaussMass * Conv_i * N[j];

                    for (unsigned int d = 0; d < TDim; d++)
                    {
                        rMassMatrix(RowIndex+d,ColIndex+d) += Mij;
                        rMassMatrix(RowIndex+TDim,ColIndex+d) += GaussWeight * mDN_DX(i,d) * Tau_t * Density * N[j];
                    }

                    // Update iteration indices
                    ColIndex += BlockSize;
                }

                // Update iteration indices
                RowIndex += BlockSize;
                ColIndex = 0;
            }
        }
    }
}


template<>
void DynamicVMS<2>::EquationIdVector(EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes * 3;
    unsigned int Index = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        NodeType& rNode = rGeom[i];
        rResult[Index++] = rNode.GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = rNode.GetDof(VELOCITY_Y).EquationId();
        rResult[Index++] = rNode.GetDof(PRESSURE).EquationId();
    }
}


template<>
void DynamicVMS<3>::EquationIdVector(EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes * 4;
    unsigned int Index = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        NodeType& rNode = rGeom[i];
        rResult[Index++] = rNode.GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = rNode.GetDof(VELOCITY_Y).EquationId();
        rResult[Index++] = rNode.GetDof(VELOCITY_Z).EquationId();
        rResult[Index++] = rNode.GetDof(PRESSURE).EquationId();
    }
}


template<>
void DynamicVMS<2>::GetDofList(DofsVectorType &rElementalDofList, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes * 3;
    unsigned int Index = 0;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        NodeType& rNode = rGeom[i];
        rElementalDofList[Index++] = rNode.pGetDof(VELOCITY_X);
        rElementalDofList[Index++] = rNode.pGetDof(VELOCITY_Y);
        rElementalDofList[Index++] = rNode.pGetDof(PRESSURE);
    }
}


template<>
void DynamicVMS<3>::GetDofList(DofsVectorType &rElementalDofList, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes * 4;
    unsigned int Index = 0;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        NodeType& rNode = rGeom[i];
        rElementalDofList[Index++] = rNode.pGetDof(VELOCITY_X);
        rElementalDofList[Index++] = rNode.pGetDof(VELOCITY_Y);
        rElementalDofList[Index++] = rNode.pGetDof(VELOCITY_Z);
        rElementalDofList[Index++] = rNode.pGetDof(PRESSURE);
    }
}


template< unsigned int TDim >
void DynamicVMS<TDim>::GetFirstDerivativesVector(Vector &rValues, int Step)
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
void DynamicVMS<TDim>::GetSecondDerivativesVector(Vector &rValues, int Step)
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
void DynamicVMS<TDim>::Calculate(const Variable<array_1d<double,3> > &rVariable, array_1d<double,3> &rOutput, const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == ADVPROJ)
    {
        GeometryType& rGeom = this->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();

        // Shape functions and integration points
        const Matrix& NContainer = rGeom.ShapeFunctionsValues(this->mIntegrationMethod);
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(this->mIntegrationMethod);
        const unsigned int NumGauss = IntegrationPoints.size();

        VectorType MomentumTerm = ZeroVector(NumNodes*TDim);
        VectorType MassTerm = ZeroVector(NumNodes);
        VectorType NodalArea = ZeroVector(NumNodes);

        for ( unsigned int g = 0; g < NumGauss; g++ )
        {
            const ShapeFunctionsType& N = row(NContainer,g);
            const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();

            double Density = 0.0;
            array_1d<double,3> ConvVel = ZeroVector(3);

            this->EvaluateInPoint(Density, DENSITY, N);
            this->EvaluateConvVelocity(ConvVel,mSubscaleVel[g],N);

            // Output containers
            array_1d< double, 3 > ElementalMomRes = ZeroVector(3);
            double ElementalMassRes(0.0);

            this->OSSMomentumResidual(ElementalMomRes,Density,ConvVel,N);
            this->MassResidual(ElementalMassRes);

            unsigned int Index = 0;
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                for (unsigned int d = 0; d < TDim; d++)
                    MomentumTerm[Index++] += GaussWeight * N[i] * ElementalMomRes[d];
                MassTerm[i] += GaussWeight * N[i] * ElementalMassRes;
                NodalArea[i] += GaussWeight * N[i];
            }
        }

        // Carefully write results to nodal variables, to avoid parallelism problems
        unsigned int Index = 0;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            NodeType& rNode = rGeom[i];
            rNode.SetLock(); // So it is safe to write in the node in OpenMP
            array_1d< double, 3 > & rAdvProj = rNode.FastGetSolutionStepValue(ADVPROJ);
            for (unsigned int d = 0; d < TDim; d++)
                rAdvProj[d] += MomentumTerm[Index++];

            rNode.FastGetSolutionStepValue(DIVPROJ) += MassTerm[i];
            rNode.FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
            rNode.UnSetLock(); // Free the node for other threads
        }
    }
}

template< unsigned int TDim >
void DynamicVMS<TDim>::CalculateOnIntegrationPoints(const Variable<double> &rVariable, std::vector<double> &rValues, const ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(this->mIntegrationMethod);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(this->mIntegrationMethod);
    const unsigned int NumGauss = IntegrationPoints.size();

    if (rVariable == SUBSCALE_PRESSURE)
    {
        rValues.resize(NumGauss);

        double Density = 0.0;
        double Viscosity = 0.0;
        array_1d<double,3> ConvVel = ZeroVector(3);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& N = row(NContainer,g);

            this->EvaluateInPoint(Density,DENSITY,N);
            this->EvaluateViscosity(Viscosity,N);
            this->EvaluateConvVelocity(ConvVel,mSubscaleVel[g],N);
            double ConvVelNorm = ConvVel[0]*ConvVel[0];
            for (unsigned int d = 1; d < TDim; d++)
                ConvVelNorm += ConvVel[d]*ConvVel[d];
            ConvVelNorm = sqrt(ConvVelNorm);

            double Tau_2 = this->TauTwo(Density,Viscosity,ConvVelNorm);
            double ElemMassRes = 0.0;
            this->MassResidual(ElemMassRes);

            if (rCurrentProcessInfo[OSS_SWITCH]==1.0)
            {
                double MassProj = 0.0;
                this->EvaluateInPoint(MassProj,DIVPROJ,N);
                ElemMassRes -= MassProj;
            }
            rValues[g] = Tau_2 * ElemMassRes;
        }
    }
    else if(rVariable == FLAG_VARIABLE) /// @todo Use suitable variable as iteration counter
    {
        rValues.resize(NumGauss);

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            rValues[g] = double(mIterCount[g]);
            mIterCount[g] = 0;
        }
    }
}

template< unsigned int TDim >
void DynamicVMS<TDim>::CalculateOnIntegrationPoints(const Variable<array_1d<double,3> > &rVariable, std::vector<array_1d<double, 3 > > &rValues, const ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(this->mIntegrationMethod);
    const unsigned int NumGauss = IntegrationPoints.size();

    if (rVariable == SUBSCALE_VELOCITY)
    {
        rValues = mSubscaleVel;
    }
    else if (rVariable == VORTICITY)
    {
        rValues.resize(NumGauss);

        for (unsigned int g = 0; g < NumGauss; g++)
            this->EvaluateVorticity(rValues[g],mDN_DX);
    }
}


template< unsigned int TDim >
void DynamicVMS<TDim>::SetValuesOnIntegrationPoints(const Variable<double> &rVariable, std::vector<double> &rValues, const ProcessInfo &rCurrentProcessInfo)
{
    unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(this->mIntegrationMethod);

    KRATOS_TRY;
    if (rValues.size() != NumGauss)
    {
        KRATOS_THROW_ERROR(std::runtime_error,"Wrong number of input values for DynamicVMS::SetValueOnIntegrationPoints. Expected number of input values is ",NumGauss);
    }
    KRATOS_CATCH("");

//    mExtraTerm = std::vector< array_1d<double,3> >(NumGauss,array_1d<double,3>(3,0.0));

//    if (rVariable == TAUONE/*SUBSCALE_VELOCITY_X*/)
//    {
//        for (unsigned int g = 0; g < NumGauss; g++)
//            mExtraTerm[g][0] = rValues[g];
//    }
//    else if (rVariable == TAUTWO/*SUBSCALE_VELOCITY_Y*/)
//    {
//        for (unsigned int g = 0; g < NumGauss; g++)
//            mExtraTerm[g][1] = rValues[g];
//    }
//    else if (TDim == 2 && rVariable == SUBSCALE_VELOCITY_Z)
//    {
//        for (unsigned int g = 0; g < NumGauss; g++)
//            mExtraTerm[g][2] = rValues[g];
//    }
}

template< unsigned int TDim >
GeometryData::IntegrationMethod DynamicVMS<TDim>::GetIntegrationMethod() const
{
    return this->mIntegrationMethod;
}


template< unsigned int TDim >
std::string DynamicVMS<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "DynamicVMS" << TDim << "D #" << Id();
    return buffer.str();
}

template< unsigned int TDim >
void DynamicVMS<TDim>::PrintInfo(std::ostream &rOStream) const
{
    rOStream << "DynamicVMS" << TDim << "D #" << Id() << std::endl;
}

template< unsigned int TDim >
void DynamicVMS<TDim>::PrintData(std::ostream &rOStream) const
{
    rOStream << "DynamicVMS" << TDim << "D #" << Id();
    rOStream << "Geometry:" << std::endl;
    this->GetGeometry().PrintData(rOStream);
    rOStream << "Integration method: " << this->mIntegrationMethod << std::endl;
}

// protected DynamicVMS methods ***********************************************

template< unsigned int TDim >
void DynamicVMS<TDim>::CalculateGeometryData()
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients( GeometryData::GI_GAUSS_1 );

    // Temporary container for inverse of J
    Matrix InvJ;

    GeometryType::JacobiansType J;
    rGeom.Jacobian( J, GeometryData::GI_GAUSS_1 );

    // calculate inverse of the jacobian and its determinant
    MathUtils<double>::InvertMatrix( J[0], InvJ, mDetJ );

    // calculate the shape function derivatives in global coordinates
    mDN_DX.resize(NumNodes,TDim, false);
    noalias( mDN_DX ) = prod( DN_De[0], InvJ );

    // calculate minimum element length (used in stabilization Tau)
    /// @todo Use heights
    array_1d<double,3> Edge = ZeroVector(3);
    Edge = rGeom[1].Coordinates() - rGeom[0].Coordinates();
    mElemSize = Edge[0]*Edge[0];
    for (SizeType d = 1; d < TDim; d++)
        mElemSize += Edge[d]*Edge[d];

    for (SizeType i = 2; i < NumNodes; i++)
        for(SizeType j = 0; j < i; j++)
        {
            Edge = rGeom[i].Coordinates() - rGeom[j].Coordinates();
            double Length = Edge[0]*Edge[0];
            for (SizeType d = 1; d < TDim; d++)
                Length += Edge[d]*Edge[d];
            if (Length < mElemSize) mElemSize = Length;
        }
    mElemSize = sqrt(mElemSize);
}

template< unsigned int TDim >
void DynamicVMS<TDim>::UpdateSubscale(const ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int NumGauss = rGeom.IntegrationPointsNumber(this->mIntegrationMethod);
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(this->mIntegrationMethod);

    double Density = 0.0;
    double Viscosity = 0.0;
    array_1d<double,3> CoarseConvVel = ZeroVector(3);

    // Elemental large-scale velocity gradient
    Matrix CoarseVelGradient = ZeroMatrix(TDim,TDim);
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rCoarseVelocity = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int m = 0; m < TDim; m++)
            for (unsigned int n = 0; n < TDim; n++)
                CoarseVelGradient(m,n) += mDN_DX(i,n) * rCoarseVelocity[m];
    }

    const double Dt = rCurrentProcessInfo[DELTA_TIME];

    // A non-linear equation must be solved for each integration point
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const ShapeFunctionsType& N = row(NContainer,g);

        array_1d<double,3>& rSubscale = mSubscaleVel[g];
        const array_1d<double,3>& rOldSubscale = mOldSubscaleVel[g];
        this->EvaluateConvVelocity(CoarseConvVel,N);
        this->EvaluateInPoint(Density,DENSITY,N);
        this->EvaluateViscosity(Viscosity,N);

        // Part of the residual that does not depend on the subscale
        array_1d<double,3> StaticResidual = ZeroVector(3);
        if ( rCurrentProcessInfo[OSS_SWITCH] == 1.0 )
        {
            this->OSSMomentumResidual(StaticResidual,Density,CoarseConvVel,N);
            array_1d<double,3> MomentumProjection = ZeroVector(3);
            this->EvaluateInPoint(MomentumProjection,ADVPROJ,N);
            StaticResidual -= MomentumProjection;
        }
        else
            this->ASGSMomentumResidual(StaticResidual,Density,CoarseConvVel,N);

        // Add the time discretization term to obtain the part of the residual that does not change during iteration
        const double MassTerm = Density / Dt;
        for (unsigned int d = 0; d < TDim; d++)
            StaticResidual[d] += MassTerm * rOldSubscale[d];

//        if (mExtraTerm.size() == NumGauss)
//            StaticResidual += mExtraTerm[g];

        // Newton-Raphson iterations for the subscale
        unsigned int Iter = 0;
        double SubscaleError = 2.0 * mSubscaleTol;
        double SubscaleNorm = 0.0;
        double ResidualNorm = 2.0 * mSubscaleTol;
        double InvTau = 0.0;
        Matrix J = ZeroMatrix(TDim,TDim);
        Vector RHS = ZeroVector(TDim);
        Vector U = ZeroVector(TDim);
        Vector dU = ZeroVector(TDim);
        for(unsigned int d = 0; d < TDim; d++)
            U[d] = rSubscale[d];

        while (Iter < 10 && SubscaleError > mSubscaleTol )
        {
            // Update iteration counter
            ++Iter;

            // Calculate new Tau
            double TmpV = CoarseConvVel[0] + U[0];
            double ConvVelNorm = TmpV*TmpV;
            for (unsigned int d = 1; d < TDim; d++)
            {
                TmpV = CoarseConvVel[d] + U[d];
                ConvVelNorm += TmpV * TmpV;
            }

            ConvVelNorm = sqrt(ConvVelNorm);
            InvTau = InvTauTime(Density,Viscosity,ConvVelNorm,Dt);

            // Newton-Raphson LHS
            noalias(J) = Density * CoarseVelGradient;
            for (unsigned int d = 0; d < TDim; d++)
                J(d,d) += InvTau;

            // Newton-Raphson RHS
            noalias(RHS) = StaticResidual;
            noalias(RHS) -= prod(J,U);

            ResidualNorm = RHS[0]*RHS[0];
            for (unsigned int d = 1; d < TDim; d++)
                ResidualNorm += RHS[d]*RHS[d];

            if (ResidualNorm > mSubscaleRHSTol)
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

//        std::cout << "Id: " << this->Id() << " g: " << g << " iter: " << Iter << " ss err: " << SubscaleError << " ss norm: " << SubscaleNorm << " residual: " << ResidualNorm << std::endl;

        // Store new subscale values
        for(unsigned int d = 0; d < TDim; d++)
            rSubscale[d] = U[d];
        mIterCount[g] += Iter;
    }
}

template< unsigned int TDim >
void DynamicVMS<TDim>::LinearUpdateSubscale(const ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumGauss = rGeom.IntegrationPointsNumber(this->mIntegrationMethod);
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(this->mIntegrationMethod);

    const double Dt = rCurrentProcessInfo[DELTA_TIME];
    if (Dt > 0.0)
    {
        const double InvDt = 1.0 /Dt;

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& N = row(NContainer,g);

            // Evaluate problem parameters at integration point
            double Density = 0.0;
            double Viscosity = 0.0;
            array_1d<double,3> ConvVel= ZeroVector(3);
            array_1d<double,3> BodyForce = ZeroVector(3);

            this->EvaluateInPoint(Density,DENSITY,N);
            this->EvaluateViscosity(Viscosity,N);
            this->EvaluateConvVelocity(ConvVel,N);
            this->EvaluateInPoint(BodyForce,BODY_FORCE,N);

            // Calculate stablitization parameters
            double ConvVelNorm = ConvVel[0] * ConvVel[0];
            for (unsigned int d = 1; d < TDim; d++)
                ConvVelNorm += ConvVel[d] * ConvVel[d];
            ConvVelNorm = sqrt(ConvVelNorm);

            double Tau_t = this->TauTime(Density,Viscosity,ConvVelNorm,Dt);

            array_1d<double,3> ElementalMomRes = ZeroVector(3);
            if (rCurrentProcessInfo[OSS_SWITCH] == 1.0)
            {
                this->OSSMomentumResidual(ElementalMomRes,Density,ConvVel+mOldSubscaleVel[g],N);
                array_1d<double,3> ElementalMomProj = ZeroVector(3);
                this->EvaluateInPoint(ElementalMomProj,ADVPROJ,N);
                ElementalMomRes -= ElementalMomProj;
            }
            else
                this->ASGSMomentumResidual(ElementalMomRes,Density,ConvVel,N);
            array_1d<double,3> OldSubscale = mOldSubscaleVel[g];
            mSubscaleVel[g] = Tau_t * ( ElementalMomRes + Density * OldSubscale * InvDt);
        }
    }
}


template< unsigned int TDim >
void DynamicVMS<TDim>::CalculateASGSVelocityContribution(MatrixType &rDampingMatrix, VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int BlockSize = TDim + 1;
    const unsigned int LocalSize = NumNodes * BlockSize;

    if(rDampingMatrix.size1() != LocalSize)
        rDampingMatrix.resize(LocalSize,LocalSize,false);

    noalias(rDampingMatrix) = ZeroMatrix(LocalSize,LocalSize);

    if(rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize,false);

    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(this->mIntegrationMethod);
    const unsigned int NumGauss = IntegrationPoints.size();
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(this->mIntegrationMethod);

    const double Dt = rCurrentProcessInfo[DELTA_TIME];
    const double InvDt = 1.0 /Dt;

    double ViscCoeff = 0.0;
    Matrix ViscousTerm = ZeroMatrix(LocalSize,LocalSize);
    this->AddViscousTerm(ViscousTerm,1.0,mDN_DX);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const ShapeFunctionsType& N = row(NContainer,g);
        const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();

        // Evaluate problem parameters at integration point
        double Density = 0.0;
        double Viscosity = 0.0;
        const array_1d<double,3>& rSubscaleVel = mSubscaleVel[g];
        array_1d<double,3> ConvVel = ZeroVector(3);
        array_1d<double,3> BodyForce = ZeroVector(3);
        Vector Convection = ZeroVector(NumNodes);

        this->EvaluateInPoint(Density,DENSITY,N);
        this->EvaluateViscosity(Viscosity,N);
        this->EvaluateConvVelocity(ConvVel,rSubscaleVel,N);
        this->EvaluateInPoint(BodyForce,BODY_FORCE,N);
        this->ConvectionOperator(Convection,ConvVel);
        // These terms are always multiplied by density in the equations
        BodyForce *= Density;
        Convection *= Density;

        array_1d<double,3> OldSubscaleTerm = mOldSubscaleVel[g];
        OldSubscaleTerm *= Density * InvDt;

        // Calculate stablitization parameters
        double ConvVelNorm = ConvVel[0] * ConvVel[0];
        for (unsigned int d = 1; d < TDim; d++)
            ConvVelNorm += ConvVel[d] * ConvVel[d];
        ConvVelNorm = sqrt(ConvVelNorm);

        double Tau_t = this->TauTime(Density,Viscosity,ConvVelNorm,Dt);
        double Tau_2 = this->TauTwo(Density,Viscosity,ConvVelNorm);

        // Contribution to viscous term (added after the Gauss point loop, as viscous matrix is constant)
        ViscCoeff += Density * Viscosity * GaussWeight;

        // Add integration point contributions to local system
        unsigned int RowIndex = 0;
        unsigned int ColIndex = 0;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            double Ni_mod = N[i] + Convection[i] * Tau_t; // Perturbed velocity shape function

            for (unsigned int j = 0; j < NumNodes; j++)
            {
                // Convective term + stabilization a * grad(v) tau_t a * grad(u)
                double Cij = GaussWeight * Ni_mod * Convection[j];

                for (unsigned int d = 0; d < TDim; d++)
                    rDampingMatrix(RowIndex+d,ColIndex+d) += Cij;

                // Pressure subscale term div(v) tau_2 div(u)
                for (unsigned int m = 0; m < TDim; m++)
                    for (unsigned int n = 0; n < TDim; n++)
                        rDampingMatrix(RowIndex+m,ColIndex+n) += GaussWeight * mDN_DX(i,m) * Tau_2 * mDN_DX(j,n);


                // Pressure term and velocity divergence
                for (unsigned int d = 0; d < TDim; d++)
                {
                    // div(v) * p and q * div(u)
                    double Gij = GaussWeight * mDN_DX(i,d) * N[j];
                    // a * grad(v) Tau_t grad(p) and grad(q) Tau_t a * grad(u)
                    double StabGij = GaussWeight * Convection[i] * Tau_t * mDN_DX(j,d);

                    rDampingMatrix(RowIndex+d,ColIndex+TDim) += StabGij - Gij;
                    rDampingMatrix(ColIndex+TDim,RowIndex+d) += Gij + StabGij;
                }

                // Pressure Laplacian grad(q) Tau_t grad(p)
                double Lij =  mDN_DX(i,0) * mDN_DX(j,0);
                for (unsigned int d = 1; d < TDim; d++)
                    Lij += mDN_DX(i,d) * mDN_DX(j,d);
                rDampingMatrix(RowIndex+TDim,ColIndex+TDim) += GaussWeight * Tau_t * Lij;

                // Update iteration indices
                ColIndex += BlockSize;
            }

            // Body force
            double vF;
            double qF = 0.0;
            for (unsigned int d = 0; d < TDim; d++)
            {
                // v f and a*grad(v) Tau_t f
                vF = BodyForce[d] * Ni_mod;
                // Tracking term a*grad(v) * OldSubscaleU / Dt
                vF += Convection[i] * Tau_t * OldSubscaleTerm[d];
                rRightHandSideVector[RowIndex+d] += GaussWeight * vF;
                // Grad(q) Tau_t f + Tracking term grad(q) * OldSubscaleU / Dt
                qF += mDN_DX(i,d) * Tau_t * ( BodyForce[d] + OldSubscaleTerm[d] );
            }
            rRightHandSideVector[RowIndex+TDim] += GaussWeight * qF;

            // Update iteration indices
            RowIndex += BlockSize;
            ColIndex = 0;
        }

    }

    // Add viscous contribution
    noalias(rDampingMatrix) += ViscCoeff * ViscousTerm;

    // Write local contributions in residual form
    Vector U = ZeroVector(LocalSize);
    this->GetFirstDerivativesVector(U);
    noalias(rRightHandSideVector) -= prod(rDampingMatrix,U);
}


template< unsigned int TDim >
void DynamicVMS<TDim>::CalculateOSSVelocityContribution(MatrixType &rDampingMatrix, VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int BlockSize = TDim + 1;
    const unsigned int LocalSize = NumNodes * BlockSize;

    if(rDampingMatrix.size1() != LocalSize)
        rDampingMatrix.resize(LocalSize,LocalSize,false);

    noalias(rDampingMatrix) = ZeroMatrix(LocalSize,LocalSize);

    if(rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize,false);

    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(this->mIntegrationMethod);
    const unsigned int NumGauss = IntegrationPoints.size();
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(this->mIntegrationMethod);

    const double Dt = rCurrentProcessInfo[DELTA_TIME];
    const double InvDt = 1.0 /Dt;

    double ViscCoeff = 0.0;
    Matrix ViscousTerm = ZeroMatrix(LocalSize,LocalSize);
    this->AddViscousTerm(ViscousTerm,1.0,mDN_DX);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const ShapeFunctionsType& N = row(NContainer,g);
        const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();

        // Evaluate problem parameters at integration point
        double Density = 0.0;
        double Viscosity = 0.0;
        array_1d<double,3> ConvVel = ZeroVector(3);
        array_1d<double,3> BodyForce = ZeroVector(3);
        Vector Convection = ZeroVector(NumNodes);

        this->EvaluateInPoint(Density,DENSITY,N);
        this->EvaluateViscosity(Viscosity,N);
        this->EvaluateConvVelocity(ConvVel,mSubscaleVel[g],N);
        this->EvaluateInPoint(BodyForce,BODY_FORCE,N);
        this->ConvectionOperator(Convection,ConvVel);
        // These terms are always multiplied by density in the equations
        BodyForce *= Density;
        Convection *= Density;

        array_1d<double,3> MomentumProjection = ZeroVector(3);
        double MassProjection = 0.0;
        this->EvaluateInPoint(MomentumProjection,ADVPROJ,N);
        this->EvaluateInPoint(MassProjection,DIVPROJ,N);


        // Calculate stablitization parameters
        double ConvVelNorm = ConvVel[0] * ConvVel[0];
        for (unsigned int d = 1; d < TDim; d++)
            ConvVelNorm += ConvVel[d] * ConvVel[d];
        ConvVelNorm = sqrt(ConvVelNorm);

        double Tau_t = this->TauTime(Density,Viscosity,ConvVelNorm,Dt);
        double Tau_2 = this->TauTwo(Density,Viscosity,ConvVelNorm);

        // Constant part of the stabilization RHS
        array_1d<double,3> MomentumStabRHS = BodyForce + Density * mOldSubscaleVel[g] * InvDt - MomentumProjection;
        MomentumStabRHS *= Tau_t;

        // Contribution to viscous term (added after the Gauss point loop, as viscous matrix is constant)
        ViscCoeff += Density * Viscosity * GaussWeight;

        // Add integration point contributions to local system
        unsigned int RowIndex = 0;
        unsigned int ColIndex = 0;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            double Ni_mod = N[i] + Convection[i] * Tau_t; // Perturbed velocity shape function

            for (unsigned int j = 0; j < NumNodes; j++)
            {
                // Convective term + stabilization a * grad(v) tau_t a * grad(u)
                double Cij = GaussWeight * Ni_mod * Convection[j];

                for (unsigned int d = 0; d < TDim; d++)
                    rDampingMatrix(RowIndex+d,ColIndex+d) += Cij;

                // Pressure subscale term div(v) tau_2 div(u)
                for (unsigned int m = 0; m < TDim; m++)
                    for (unsigned int n = 0; n < TDim; n++)
                        rDampingMatrix(RowIndex+m,ColIndex+n) += GaussWeight * mDN_DX(i,m) * Tau_2 * mDN_DX(j,n);


                // Pressure term and velocity divergence
                for (unsigned int d = 0; d < TDim; d++)
                {
                    // div(v) * p and q * div(u)
                    double Gij = GaussWeight * mDN_DX(i,d) * N[j];
                    // a * grad(v) Tau_t grad(p) and grad(q) Tau_t a * grad(u)
                    double StabGij = GaussWeight * Convection[i] * Tau_t * mDN_DX(j,d);

                    rDampingMatrix(RowIndex+d,ColIndex+TDim) += StabGij - Gij;
                    rDampingMatrix(ColIndex+TDim,RowIndex+d) += Gij + StabGij;
                }


                // Pressure Laplacian grad(q) Tau_t grad(p)
                double Lij =  mDN_DX(i,0) * mDN_DX(j,0);
                for (unsigned int d = 1; d < TDim; d++)
                    Lij += mDN_DX(i,d) * mDN_DX(j,d);
                rDampingMatrix(RowIndex+TDim,ColIndex+TDim) += GaussWeight * Tau_t * Lij;

                // Update iteration indices
                ColIndex += BlockSize;
            }

            // Body force
            double vF;
            double qF = 0.0;
            for (unsigned int d = 0; d < TDim; d++)
            {
                // v f
                vF = BodyForce[d] * N[i];
                // a*grad(v) Tau_t f and Tracking term and projection: a*grad(v) * ( OldSubscaleU / Dt - MomentumProjection )
                vF += Convection[i] * MomentumStabRHS[d];
                // Mass Projection term  div(v) * Tau_2 * MassProj
                vF -= mDN_DX(i,d) * Tau_2 * MassProjection;
                rRightHandSideVector[RowIndex+d] += GaussWeight * vF;
                // Grad(q) Tau_t f + Tracking term grad(q) * OldSubscaleU / Dt - Projection term grad(q) * MomentumProjection
                qF += mDN_DX(i,d) * MomentumStabRHS[d];
            }
            rRightHandSideVector[RowIndex+TDim] += GaussWeight * qF;

            // Update iteration indices
            RowIndex += BlockSize;
            ColIndex = 0;
        }

    }

    // Add viscous contribution
    noalias(rDampingMatrix) += ViscCoeff * ViscousTerm;

    // Write local contributions in residual form
    Vector U = ZeroVector(LocalSize);
    this->GetFirstDerivativesVector(U);
    noalias(rRightHandSideVector) -= prod(rDampingMatrix,U);
}

template< unsigned int TDim >
void DynamicVMS<TDim>::LumpedMassMatrix(MatrixType &rMassMatrix)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int BlockSize = TDim + 1;

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(this->mIntegrationMethod);
    const unsigned int NumGauss = IntegrationPoints.size();
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(this->mIntegrationMethod);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const ShapeFunctionsType& N = row(NContainer,g);
        const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();

        // Evaluate problem parameters at integration point
        double Density = 0.0;
        this->EvaluateInPoint(Density,DENSITY,N);

        unsigned int Index = 0;

        const double GaussMass = GaussWeight * Density;

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            double Mij = GaussMass * N[i];

            for (unsigned int d = 0; d < TDim; d++)
                rMassMatrix(Index+d,Index+d) += Mij;

            // Update iteration indices
            Index += BlockSize;
        }
    }
}

template< unsigned int TDim >
void DynamicVMS<TDim>::ConsistentMassMatrix(MatrixType &rMassMatrix)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int BlockSize = TDim + 1;

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(this->mIntegrationMethod);
    const unsigned int NumGauss = IntegrationPoints.size();
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(this->mIntegrationMethod);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const ShapeFunctionsType& N = row(NContainer,g);
        const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();

        // Evaluate problem parameters at integration point
        double Density = 0.0;
        this->EvaluateInPoint(Density,DENSITY,N);

        unsigned int RowIndex = 0;
        unsigned int ColIndex = 0;

        const double GaussMass = GaussWeight * Density;

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            for (unsigned int j = 0; j < NumNodes; j++)
            {
                double Mij = GaussMass * N[i] * N[j];

                for (unsigned int d = 0; d < TDim; d++)
                {
                    rMassMatrix(RowIndex+d,ColIndex+d) += Mij;
                }

                // Update iteration indices
                ColIndex += BlockSize;
            }

            // Update iteration indices
            RowIndex += BlockSize;
            ColIndex = 0;
        }
    }
}

template< unsigned int TDim >
void DynamicVMS<TDim>::EvaluateConvVelocity(array_1d<double,3>& rConvection, const ShapeFunctionsType &rN)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    array_1d<double,3> NodeVel = rGeom[0].FastGetSolutionStepValue(VELOCITY);
    NodeVel -= rGeom[0].FastGetSolutionStepValue(MESH_VELOCITY);
    rConvection = rN[0] * NodeVel;

    for (unsigned int i = 1; i < NumNodes; i++)
    {
        NodeVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        NodeVel -= rGeom[i].FastGetSolutionStepValue(MESH_VELOCITY);
        rConvection += rN[i] * NodeVel;
    }
}

template< unsigned int TDim >
void DynamicVMS<TDim>::EvaluateConvVelocity(array_1d<double,3> &rConvection, const array_1d<double,3> &rSubscaleVel, const ShapeFunctionsType &rN)
{
    this->EvaluateConvVelocity(rConvection,rN);
    rConvection += rSubscaleVel;
}

template< unsigned int TDim >
void  DynamicVMS<TDim>::EvaluateViscosity(double &rViscosity, const ShapeFunctionsType &rN)
{
    this->EvaluateInPoint(rViscosity,VISCOSITY,rN);
}

template< unsigned int TDim >
void DynamicVMS<TDim>::ConvectionOperator(Vector &rResult, const array_1d<double,3> &rConvVel)
{
    const unsigned int NumNodes = rResult.size();
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        rResult[i] = rConvVel[0] * mDN_DX(i,0);
        for (unsigned int d = 1; d < TDim; d++)
            rResult[i] += rConvVel[d] * mDN_DX(i,d);
    }
}

template< unsigned int TDim >
double DynamicVMS<TDim>::TauOne(const double Density,const double Viscosity,const double ConvVel)
{
    const double c1 = 4.0;
    const double c2 = 2.0;
    double InvTau = c1 * Viscosity / (mElemSize * mElemSize) + c2 * ConvVel / mElemSize;
    return 1.0 / (Density * InvTau);
}

template< unsigned int TDim >
double DynamicVMS<TDim>::TauTwo(const double Density, const double Viscosity, const double ConvVel)
{
    return Density * (Viscosity + 0.5 * ConvVel * mElemSize);
}

template< unsigned int TDim >
double DynamicVMS<TDim>::TauTime(const double Density, const double Viscosity, const double ConvVel, const double Dt)
{
    const double c1 = 4.0;
    const double c2 = 2.0;
    double InvTau = 1.0 / Dt + c1 * Viscosity / (mElemSize * mElemSize) + c2 * ConvVel / mElemSize;
    return 1.0 / (Density * InvTau);
}

template< unsigned int TDim >
double DynamicVMS<TDim>::InvTauTime(const double Density, const double Viscosity, const double ConvVel, const double Dt)
{
    const double c1 = 4.0;
    const double c2 = 2.0;
    double InvTau = 1.0 / Dt + c1 * Viscosity / (mElemSize * mElemSize) + c2 * ConvVel / mElemSize;
    return Density * InvTau;
}

template< unsigned int TDim >
void DynamicVMS<TDim>::ASGSMomentumResidual(array_1d<double,3> &rResult,
                                              const double Density,
                                              const array_1d<double,3> &rConvVel,
                                              const ShapeFunctionsType &rN)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        NodeType& rNode = rGeom[i];
        const array_1d<double,3>& rVelocity = rNode.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double,3>& rAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION);
        const array_1d<double,3>& rBodyForce = rNode.FastGetSolutionStepValue(BODY_FORCE);
        double Pressure = rNode.FastGetSolutionStepValue(PRESSURE);

        double Conv = 0.0;
        for (unsigned int d = 0; d < TDim; d++)
            Conv += rConvVel[d] * mDN_DX(i,d);

        rResult += Density * ( rN[i] * (rBodyForce - rAcceleration) - Conv * rVelocity );

        for (unsigned int d = 0; d < TDim; d++)
            rResult[d] -= mDN_DX(i,d) * Pressure;
    }
}

template< unsigned int TDim >
void DynamicVMS<TDim>::OSSMomentumResidual(array_1d<double,3> &rResult,
                                           const double Density,
                                           const array_1d<double,3> &rConvVel,
                                           const ShapeFunctionsType &rN)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    rResult = ZeroVector(3);;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        NodeType& rNode = rGeom[i];
        const array_1d<double,3>& rVelocity = rNode.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double,3>& rBodyForce = rNode.FastGetSolutionStepValue(BODY_FORCE);
        double Pressure = rNode.FastGetSolutionStepValue(PRESSURE);

        double Conv = 0.0;
        for (unsigned int d = 0; d < TDim; d++)
            Conv += rConvVel[d] * mDN_DX(i,d);

        rResult += Density * ( rN[i] * rBodyForce - Conv * rVelocity );

        for (unsigned int d = 0; d < TDim; d++)
            rResult[d] -= mDN_DX(i,d) * Pressure;
    }
}

template< unsigned int TDim >
void DynamicVMS<TDim>::MassResidual(double &rResult)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rVelocity = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int d = 0; d < TDim; d++)
            rResult -= mDN_DX(i,d) * rVelocity[d];
    }
}

template<>
void DynamicVMS<2>::AddViscousTerm(MatrixType &rDampingMatrix, const double Weight, const ShapeDerivativesType& rDN_DX)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = 3; // TDim + 1

    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    unsigned int FirstRow(0),FirstCol(0);

    for (unsigned int j = 0; j < NumNodes; ++j)
    {
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            // First Row
            rDampingMatrix(FirstRow,FirstCol) += Weight * ( FourThirds * rDN_DX(i,0) * rDN_DX(j,0) + rDN_DX(i,1) * rDN_DX(j,1) );
            rDampingMatrix(FirstRow,FirstCol+1) += Weight * ( nTwoThirds * rDN_DX(i,0) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,0) );

            // Second Row
            rDampingMatrix(FirstRow+1,FirstCol) += Weight * ( nTwoThirds * rDN_DX(i,1) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,1) );
            rDampingMatrix(FirstRow+1,FirstCol+1) += Weight * ( FourThirds * rDN_DX(i,1) * rDN_DX(j,1) + rDN_DX(i,0) * rDN_DX(j,0) );

            // Update Counter
            FirstRow += BlockSize;
        }
        FirstRow = 0;
        FirstCol += BlockSize;
    }
}

template<>
void DynamicVMS<3>::AddViscousTerm(MatrixType &rDampingMatrix, const double Weight, const ShapeDerivativesType& rDN_DX)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = 4; // TDim + 1

    const double OneThird = 1.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    unsigned int FirstRow(0),FirstCol(0);

    for (unsigned int j = 0; j < NumNodes; ++j)
    {
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            // (dN_i/dx_k dN_j/dx_k)
            const double Diag =  rDN_DX(i,0) * rDN_DX(j,0) + rDN_DX(i,1) * rDN_DX(j,1) + rDN_DX(i,2) * rDN_DX(j,2);

            // First Row
            rDampingMatrix(FirstRow,FirstCol) += Weight * ( OneThird * rDN_DX(i,0) * rDN_DX(j,0) + Diag );
            rDampingMatrix(FirstRow,FirstCol+1) += Weight * ( nTwoThirds * rDN_DX(i,0) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,0) );
            rDampingMatrix(FirstRow,FirstCol+2) += Weight * ( nTwoThirds * rDN_DX(i,0) * rDN_DX(j,2) + rDN_DX(i,2) * rDN_DX(j,0) );

            // Second Row
            rDampingMatrix(FirstRow+1,FirstCol) += Weight * ( nTwoThirds * rDN_DX(i,1) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,1) );
            rDampingMatrix(FirstRow+1,FirstCol+1) += Weight * ( OneThird * rDN_DX(i,1) * rDN_DX(j,1) + Diag );
            rDampingMatrix(FirstRow+1,FirstCol+2) += Weight * ( nTwoThirds * rDN_DX(i,1) * rDN_DX(j,2) + rDN_DX(i,2) * rDN_DX(j,1) );

            // Third Row
            rDampingMatrix(FirstRow+2,FirstCol) += Weight * ( nTwoThirds * rDN_DX(i,2) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,2) );
            rDampingMatrix(FirstRow+2,FirstCol+1) += Weight * ( nTwoThirds * rDN_DX(i,2) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,2) );
            rDampingMatrix(FirstRow+2,FirstCol+2) += Weight * ( OneThird * rDN_DX(i,2) * rDN_DX(j,2) + Diag );

            // Update Counter
            FirstRow += BlockSize;
        }
        FirstRow = 0;
        FirstCol += BlockSize;
    }
}

template<>
void DynamicVMS<2>::DenseSystemSolve(const Matrix &rA, const Vector &rB, Vector &rX)
{
    Matrix Inv = ZeroMatrix(2,2);
    double Det = 0.0;
    MathUtils<double>::InvertMatrix2(rA,Inv,Det);
    noalias(rX) = prod(Inv,rB);
}

template<>
void DynamicVMS<3>::DenseSystemSolve(const Matrix &rA, const Vector &rB, Vector &rX)
{
    Matrix Inv = ZeroMatrix(3,3);
    double Det = 0.0;
    MathUtils<double>::InvertMatrix3(rA,Inv,Det);
    noalias(rX) = prod(Inv,rB);
}

template<>
void DynamicVMS<2>::EvaluateVorticity(array_1d<double,3>& rVorticity,
                                      const ShapeDerivativesType& rDN_DX)
{
    rVorticity[0] = 0.0;
    rVorticity[1] = 0.0;
    rVorticity[2] = 0.0;

    const unsigned int NumNodes = this->GetGeometry().PointsNumber();

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        rVorticity[2] += (rDN_DX(i,0)*rVelocity[1] - rDN_DX(i,1)*rVelocity[0]);
    }

}

template<>
void DynamicVMS<3>::EvaluateVorticity(array_1d<double,3>& rVorticity,
                                      const ShapeDerivativesType& rDN_DX)
{
    rVorticity[0] = 0.0;
    rVorticity[1] = 0.0;
    rVorticity[2] = 0.0;

    const unsigned int NumNodes = this->GetGeometry().PointsNumber();

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        rVorticity[0] += (rDN_DX(i,1)*rVelocity[2] - rDN_DX(i,2)*rVelocity[1]);
        rVorticity[1] += (rDN_DX(i,2)*rVelocity[0] - rDN_DX(i,0)*rVelocity[2]);
        rVorticity[2] += (rDN_DX(i,0)*rVelocity[1] - rDN_DX(i,1)*rVelocity[0]);
    }
}

// private DynamicVMS methods *************************************************

template< unsigned int TDim >
const double DynamicVMS<TDim>::mSubscaleTol = 1e-16;

template< unsigned int TDim >
const double DynamicVMS<TDim>::mSubscaleRHSTol = 1e-32;

template< unsigned int TDim >
DynamicVMS<TDim>::DynamicVMS():
    Element(),
    mIntegrationMethod(GeometryData::GI_GAUSS_1)
{}

template< unsigned int TDim >
void DynamicVMS<TDim>::save(Serializer &rSerializer) const
{
    KRATOS_TRY;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    // Save Integration method enum as integer
    int IntMethod = int(mIntegrationMethod);
    rSerializer.save("IntMethod",IntMethod);
    rSerializer.save("mSubscaleVel",mSubscaleVel);
    rSerializer.save("mOldSubscaleVel",mOldSubscaleVel);
    KRATOS_CATCH("");
}

template< unsigned int TDim >
void DynamicVMS<TDim>::load(Serializer &rSerializer)
{
    KRATOS_TRY;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    // Integration method enum is stored as int
    int IntMethod = 0;
    rSerializer.load("IntMethod",IntMethod);
    mIntegrationMethod = GeometryData::IntegrationMethod(IntMethod);
    rSerializer.load("mSubscaleVel",mSubscaleVel);
    rSerializer.load("mOldSubscaleVel",mOldSubscaleVel);
    this->CalculateGeometryData();
    mIterCount.resize(this->GetGeometry().IntegrationPointsNumber(this->mIntegrationMethod),0);
    KRATOS_CATCH("");
}

// template class declarations ************************************************

template class DynamicVMS<2>;
template class DynamicVMS<3>;

}
