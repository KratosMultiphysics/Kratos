#include "fractional_step.h"

namespace Kratos {

/*
 * public FractionalStep<TDim> functions
 */

template< unsigned int TDim >
void FractionalStep<TDim>::Initialize()
{
//    this->CalculateGeometryData();
}

template< unsigned int TDim >
void FractionalStep<TDim>::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    this->CalculateGeometryData();
}

template< unsigned int TDim >
void FractionalStep<TDim>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{
//    this->CalculateGeometryData();
    double Csmag = this->GetValue(C_SMAGORINSKY);

    // Total viscosity
    const ShapeFunctionsType& N = row( this->GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_1), 0);

    double Viscosity = 0.0;
    this->EvaluateInPoint(Viscosity,VISCOSITY,N);

    if (Csmag != 0.0 )
    {
        const unsigned int NumNodes = TDim + 1;

        // Calculate Symetric gradient
        MatrixType S = ZeroMatrix(TDim,TDim);
        for (unsigned int n = 0; n < NumNodes; ++n)
        {
            const array_1d<double,3>& rVel = this->GetGeometry()[n].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int i = 0; i < TDim; ++i)
                for (unsigned int j = 0; j < TDim; ++j)
                    S(i,j) += 0.5 * ( mDN_DX(n,j) * rVel[i] + mDN_DX(n,i) * rVel[j] );
        }

        // Norm of symetric gradient
        double NormS = 0.0;
        for (unsigned int i = 0; i < TDim; ++i)
            for (unsigned int j = 0; j < TDim; ++j)
                NormS += S(i,j) * S(i,j);
        NormS = sqrt(2.0*NormS);

        // Nu_sgs = (Csmag * Delta)^2 * (2*Sij*Sij)^(1/2)
        Viscosity += Csmag * Csmag * mElemSize * mElemSize * NormS;
    }

    this->SetValue(VISCOSITY,Viscosity);
}

template< unsigned int TDim >
void FractionalStep<TDim>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                VectorType& rRightHandSideVector,
                                                ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
    case 1:
    {
        this->CalculateLocalFractionalVelocitySystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        break;
    }
    case 5:
    {
        this->CalculateLocalPressureSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        break;
    }
    case 6:
    {
        KRATOS_ERROR(std::logic_error,"Full solution of end of step velocity is not implemented, see Calculate(VELOCITY)","");
        break;
    }
    default:
    {
        KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
    }
    }

    KRATOS_CATCH("");
}

template< unsigned int TDim >
void FractionalStep<TDim>::Calculate(const Variable<double> &rVariable,
                                     double &rOutput,
                                     const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == DIVPROJ)
    {
        const SizeType NumNodes = TDim + 1;
        const GeometryType& rGeom = this->GetGeometry();

        // Shape functions and integration points
        const Matrix& NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const SizeType NumGauss = IntegrationPoints.size();

        array_1d<double,3> MomentumRHS(3,0.0);
        double MassRHS = 0.0;

        // Loop on integration points
        for (SizeType g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& N = row(NContainer,g);
            const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();

            this->CalculateProjectionRHS(MomentumRHS,MassRHS,N,mDN_DX,GaussWeight);
        }

        // distribute residual to nodes, to compute projections
        const ShapeFunctionsType CenterN = row(rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1),0);

        // Carefully write results to nodal variables, to avoid parallelism problems
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
            this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ) += MassRHS * CenterN[i];
            this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ) += MomentumRHS * CenterN[i];
            this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
        }
    }
}

template< unsigned int TDim >
void FractionalStep<TDim>::Calculate(const Variable<array_1d<double,3> > &rVariable,
                                     array_1d<double,3> &rOutput,
                                     const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == ADVPROJ)
    {
        double Tmp = 0.0;
        this->Calculate(DIVPROJ,Tmp,rCurrentProcessInfo);
    }
    else if (rVariable == CONV_PROJ)
    {
        const unsigned int NumNodes = TDim + 1;
        const unsigned int LocalSize = TDim * NumNodes;
        const GeometryType& rGeom = this->GetGeometry();

        // Shape functions and integration points
        const Matrix& NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const SizeType NumGauss = IntegrationPoints.size();

        VectorType ConvTerm = ZeroVector(LocalSize);
        VectorType PresTerm = ZeroVector(LocalSize);
        VectorType DivTerm = ZeroVector(NumNodes);
        VectorType NodalArea = ZeroVector(NumNodes);

        // Loop on integration points
        for (unsigned int g = 0; g < NumGauss; g++)
        {
            const ShapeFunctionsType& N = row(NContainer,g);
            const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();

            for (unsigned int i = 0; i < NumNodes; i++)
                NodalArea[i] += N[i] * GaussWeight;

            this->CalculateProjectionRHS(ConvTerm,PresTerm,DivTerm,N,mDN_DX,GaussWeight);
        }

        // Carefully write results to nodal variables, to avoid parallelism problems
        unsigned int RowIndex = 0;
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double,3>& rConvVal = this->GetGeometry()[i].FastGetSolutionStepValue(CONV_PROJ);
            array_1d<double,3>& rPresVal = this->GetGeometry()[i].FastGetSolutionStepValue(PRESS_PROJ);
            for (unsigned int d = 0; d < TDim; ++d)
            {
                rConvVal[d] += ConvTerm[RowIndex];
                rPresVal[d] += PresTerm[RowIndex];
                ++RowIndex;
            }
            this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ) += DivTerm[i];
            this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
            this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
        }
    }
    else if (rVariable == VELOCITY)
    {
        const SizeType NumNodes = TDim + 1;
        const SizeType LocalSize = TDim * NumNodes;

        // Shape functions and integration points
        GeometryType& rGeom = this->GetGeometry();
        const Matrix& NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const SizeType NumGauss = IntegrationPoints.size();

        VectorType NodalVelCorrection = ZeroVector(LocalSize);

        // Loop on integration points
        for (SizeType g = 0; g < NumGauss; ++g)
        {
            const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();
            const ShapeFunctionsType& N = row(NContainer,g);

            double Density;

            this->EvaluateInPoint(Density,DENSITY,N);

            const double Coeff = GaussWeight / ( Density * rCurrentProcessInfo[BDF_COEFFICIENTS][0] );

            // Calculate contribution to the gradient term (RHS)
            double DeltaPressure;
            this->EvaluateInPoint(DeltaPressure,PRESSURE_OLD_IT,N);

            SizeType RowIndex = 0;

            for (SizeType i = 0; i < NumNodes; ++i)
            {
                for (SizeType d = 0; d < TDim; ++d)
                {
                    NodalVelCorrection[RowIndex++] += Coeff * mDN_DX(i,d) * DeltaPressure;
                }
            }
        }

        SizeType Index = 0;

        for (SizeType i = 0; i < NumNodes; ++i)
        {
            rGeom[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double,3>& rTemp = rGeom[i].FastGetSolutionStepValue(FRACT_VEL);
            for (SizeType d = 0; d < TDim; ++d)
            {
                rTemp[d] += NodalVelCorrection[Index++];
            }
            rGeom[i].UnSetLock(); // Free the node for other threads
        }

    }
}

template< unsigned int TDim >
void FractionalStep<TDim>::EquationIdVector(EquationIdVectorType& rResult,
                                            ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
    case 1:
    {
        this->VelocityEquationIdVector(rResult,rCurrentProcessInfo);
        break;
    }
    case 5:
    {
        this->PressureEquationIdVector(rResult,rCurrentProcessInfo);
        break;
    }
    case 6:
    {
        this->VelocityEquationIdVector(rResult,rCurrentProcessInfo);
        break;
    }
    default:
    {
        KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
    }
    }

    KRATOS_CATCH("");
}

template< unsigned int TDim >
void FractionalStep<TDim>::GetDofList(DofsVectorType& rElementalDofList,
                                      ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
    case 1:
    {
        this->GetVelocityDofList(rElementalDofList,rCurrentProcessInfo);
        break;
    }
    case 5:
    {
        this->GetPressureDofList(rElementalDofList,rCurrentProcessInfo);
        break;
    }
    case 6:
    {
        this->GetVelocityDofList(rElementalDofList,rCurrentProcessInfo);
        break;
    }
    default:
    {
        KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
    }
    }

    KRATOS_CATCH("");
}

template< unsigned int TDim >
void FractionalStep<TDim>::CalculateLocalFractionalVelocitySystem(MatrixType& rLeftHandSideMatrix,
                                                                  VectorType& rRightHandSideVector,
                                                                  ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = TDim + 1;
    const SizeType LocalSize = TDim * NumNodes;
    const GeometryType& rGeom = this->GetGeometry();

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != LocalSize )
        rLeftHandSideMatrix.resize(LocalSize,LocalSize);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);

    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize);

    rRightHandSideVector = ZeroVector(LocalSize);

    // Shape functions and integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const SizeType NumGauss = IntegrationPoints.size();

    MatrixType MassMatrix = ZeroMatrix(LocalSize,LocalSize);

    // Stabilization parameters
    double TauOne;
    double TauTwo;
    this->CalculateTau(TauOne,TauTwo,rCurrentProcessInfo);

    // Loop on integration points
    for (SizeType g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();
        const ShapeFunctionsType& N = row(NContainer,g);

        // Evaluate required variables at the integration point
        double Density;
        double MassProjection;
        array_1d<double,3> Velocity(3,0.0);
        array_1d<double,3> MeshVelocity(3,0.0);
        array_1d<double,3> BodyForce(3,0.0);
        array_1d<double,3> MomentumProjection(3,0.0);

        this->EvaluateInPoint(Density,DENSITY,N);
        this->EvaluateInPoint(MassProjection,DIVPROJ,N);
        this->EvaluateInPoint(Velocity,VELOCITY,N);
        this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);
        this->EvaluateInPoint(BodyForce,BODY_FORCE,N);
//        this->EvaluateInPoint(MomentumProjection,ADVPROJ,N);
        this->EvaluateInPoint(MomentumProjection,CONV_PROJ,N);

        // Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
        double OldPressure;
        this->EvaluateInPoint(OldPressure,PRESSURE,N,0);

        // For ALE: convective velocity
        array_1d<double,3> ConvVel = Velocity - MeshVelocity;

//        // Stabilization parameters
//        double TauOne;
//        double TauTwo;
//        this->CalculateTau(TauOne,TauTwo,ConvVel,Density,Viscosity,rCurrentProcessInfo);

        // Evaluate convection operator Velocity * Grad(N)
        Vector UGradN(NumNodes);
        this->EvaluateConvection(UGradN,ConvVel,mDN_DX);

        // Add integration point contribution to the local mass matrix
        this->AddMomentumMassTerm(MassMatrix,N,GaussWeight*Density);

        // Add convection, stabilization and RHS contributions to the local system equation
        this->AddMomentumSystemTerms(rLeftHandSideMatrix,rRightHandSideVector,Density,UGradN,BodyForce,OldPressure,
                                     TauOne,TauTwo,MomentumProjection,MassProjection,N,mDN_DX,GaussWeight);
    }

    // Add viscous term (using a single integration point)
    const ShapeFunctionsType CenterN = row(rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1),0);
    const double CenterWeight = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_1)[0].Weight();

    double Density;
    this->EvaluateInPoint(Density,DENSITY,CenterN);
    double Viscosity = this->GetValue(VISCOSITY);

    const double ViscousCoeff = Density * Viscosity * CenterWeight * mDetJ;
    this->AddViscousTerm(rLeftHandSideMatrix,mDN_DX,ViscousCoeff);

    // Add residual of previous iteration to RHS
    VectorType LastValues = ZeroVector(LocalSize);
    this->GetVelocityValues(LastValues,0);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,LastValues);

    // Add dynamic term
    const Vector& rBDFCoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
    noalias( rLeftHandSideMatrix ) += rBDFCoeffs[0] * MassMatrix;

    VectorType TimeTerm = rBDFCoeffs[0] * LastValues;
    for(SizeType i = 1; i < rBDFCoeffs.size(); i++)
    {
        this->GetVelocityValues(LastValues,i);
        noalias( TimeTerm ) += rBDFCoeffs[i] * LastValues;
    }

    noalias( rRightHandSideVector ) -= prod(MassMatrix,TimeTerm);
}

template< unsigned int TDim >
void FractionalStep<TDim>::CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix,
                                                        VectorType& rRightHandSideVector,
                                                        ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = TDim + 1;
    GeometryType& rGeom = this->GetGeometry();

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != NumNodes )
        rLeftHandSideMatrix.resize(NumNodes,NumNodes);

    rLeftHandSideMatrix = ZeroMatrix(NumNodes,NumNodes);

    if( rRightHandSideVector.size() != NumNodes )
        rRightHandSideVector.resize(NumNodes);

    rRightHandSideVector = ZeroVector(NumNodes);

    // Shape functions and integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const SizeType NumGauss = IntegrationPoints.size();

    // Stabilization parameters
    double TauOne;
    double TauTwo;
    this->CalculateTau(TauOne,TauTwo,rCurrentProcessInfo);

    // Loop on integration points
    for (SizeType g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();
        const ShapeFunctionsType& N = row(NContainer,g);

        // Evaluate required variables at the integration point
        double Density;
//        double Viscosity;
//        array_1d<double,3> Velocity(3,0.0);
//        array_1d<double,3> MeshVelocity(3,0.0);
        array_1d<double,3> BodyForce(3,0.0);
        array_1d<double,3> MomentumProjection(3,0.0);

        this->EvaluateInPoint(Density,DENSITY,N);
//        this->EvaluateInPoint(Viscosity,VISCOSITY,N);
//        this->EvaluateInPoint(Velocity,VELOCITY,N);
//        this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);
        this->EvaluateInPoint(BodyForce,BODY_FORCE,N);
//        this->EvaluateInPoint(MomentumProjection,ADVPROJ,N);
        this->EvaluateInPoint(MomentumProjection,PRESS_PROJ,N);

//        // Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
//        double OldPressure;
//        this->EvaluateInPoint(OldPressure,PRESSURE,N,0);

        array_1d<double,TDim> OldPressureGradient(TDim,0.0);
        this->EvaluateGradientInPoint(OldPressureGradient,PRESSURE,mDN_DX);

//        // For ALE: convective velocity
//        array_1d<double,3> ConvVel = Velocity - MeshVelocity;

//        // Stabilization parameters
//        double TauOne;
//        double TauTwo;
//        this->CalculateTau(TauOne,TauTwo,ConvVel,Density,Viscosity,rCurrentProcessInfo);

//        // Evaluate convection operator Velocity * Grad(N)
//        Vector UGradN(NumNodes);
//        this->EvaluateConvection(UGradN,ConvVel,mDN_DX);

        double DivU;
        this->EvaluateDivergenceInPoint(DivU,VELOCITY,mDN_DX);

        // constant coefficient multiplying the pressure Laplacian (See Codina, Badia 2006 paper for details in case of a BDF2 time scheme)
        const double LaplacianCoeff = 1.0 / (Density * rCurrentProcessInfo[BDF_COEFFICIENTS][0]) ;

        // Add convection, stabilization and RHS contributions to the local system equation
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            // LHS contribution
            for (SizeType j = 0; j < NumNodes; ++j)
            {
                double Lij = 0.0;
                for (SizeType d = 0; d < TDim; ++d)
                    Lij += mDN_DX(i,d) * mDN_DX(j,d);
                Lij *= (LaplacianCoeff + TauOne);

                rLeftHandSideMatrix(i,j) += GaussWeight * Lij;
            }

            // RHS contribution

            // Velocity divergence
            double RHSi = - N[i] * DivU;

            for (SizeType d = 0; d < TDim; ++d)
            {
//                double Conv = UGradN[0] * rGeom[0].FastGetSolutionStepValue(VELOCITY)[d];
//                for (SizeType j = 1; j < NumNodes; ++j) Conv += UGradN[i] * rGeom[i].FastGetSolutionStepValue(VELOCITY)[d];
                // Momentum stabilization
                RHSi += mDN_DX(i,d) * TauOne * ( Density  * ( BodyForce[d]/* - Conv*/ ) - OldPressureGradient[d] - MomentumProjection[d] );
            }

            rRightHandSideVector[i] += GaussWeight * RHSi;
        }
    }
}

template< unsigned int TDim >
int FractionalStep<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    if(VELOCITY.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check that the application was correctly registered.","");
    if(PRESSURE.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check that the application was correctly registered.","");
    if(BODY_FORCE.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"BODY_FORCE Key is 0. Check that the application was correctly registered.","");
    if(DENSITY.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"DENSITY Key is 0. Check that the application was correctly registered.","");
    if(VISCOSITY.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"VISCOSITY Key is 0. Check that the application was correctly registered.","");
    if(MESH_VELOCITY.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check that the application was correctly registered.","");
    if(FRACT_VEL.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"FRACT_VEL Key is 0. Check that the application was correctly registered.","");
    if(PRESSURE_OLD_IT.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"PRESSURE_OLD_IT Key is 0. Check that the application was correctly registered.","");
    if(NODAL_AREA.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"NODAL_AREA Key is 0. Check that the application was correctly registered.","");
    if(CONV_PROJ.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"CONV_PROJ Key is 0. Check that the application was correctly registered.","");
    if(PRESS_PROJ.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"PRESS_PROJ Key is 0. Check that the application was correctly registered.","");
    if(DIVPROJ.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"DIVPROJ Key is 0. Check that the application was correctly registered.","");
    if(BDF_COEFFICIENTS.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"BDF_COEFFICIENTS Key is 0. Check that the application was correctly registered.","");
    if(DELTA_TIME.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check that the application was correctly registered.","");
    if(DYNAMIC_TAU.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"DYNAMIC_TAU Key is 0. Check that the application was correctly registered.","");

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
            KRATOS_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
            KRATOS_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE) == false)
            KRATOS_ERROR(std::invalid_argument,"missing BODY_FORCE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
            KRATOS_ERROR(std::invalid_argument,"missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(VISCOSITY) == false)
            KRATOS_ERROR(std::invalid_argument,"missing VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
            KRATOS_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(FRACT_VEL) == false)
            KRATOS_ERROR(std::invalid_argument,"missing FRACT_VEL variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE_OLD_IT) == false)
            KRATOS_ERROR(std::invalid_argument,"missing PRESSURE_OLD_IT variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(NODAL_AREA) == false)
            KRATOS_ERROR(std::invalid_argument,"missing NODAL_AREA variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(CONV_PROJ) == false)
            KRATOS_ERROR(std::invalid_argument,"missing CONV_PROJ variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESS_PROJ) == false)
            KRATOS_ERROR(std::invalid_argument,"missing PRESS_PROJ variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(DIVPROJ) == false)
            KRATOS_ERROR(std::invalid_argument,"missing DIVPROJ variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
           this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
           this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
            KRATOS_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
            KRATOS_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if (this->GetGeometry().WorkingSpaceDimension() == 2)
    {
        for (unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if (this->GetGeometry()[i].Z() != 0.0)
                KRATOS_ERROR(std::invalid_argument,"Node with non-zero Z coordinate found. Id: ",this->GetGeometry()[i].Id());
        }
    }

    return ierr;

    KRATOS_CATCH("");
}

template<>
void FractionalStep<2>::VelocityEquationIdVector(EquationIdVectorType& rResult,
                                                 ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 6;
    GeometryType& rGeom = this->GetGeometry();

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
    }
}

template<>
void FractionalStep<3>::VelocityEquationIdVector(EquationIdVectorType& rResult,
                                                 ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 4;
    const SizeType LocalSize = 12;
    GeometryType& rGeom = this->GetGeometry();

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
    }
}

template< unsigned int TDim >
void FractionalStep<TDim>::PressureEquationIdVector(EquationIdVectorType& rResult,
                                                    ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = TDim + 1;
    GeometryType& rGeom = this->GetGeometry();

    if (rResult.size() != NumNodes)
        rResult.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; ++i)
        rResult[i] = rGeom[i].GetDof(PRESSURE).EquationId();
}

template<>
void FractionalStep<2>::GetVelocityDofList(DofsVectorType& rElementalDofList,
                                           ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 6;
    GeometryType& rGeom = this->GetGeometry();

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
    }
}

template<>
void FractionalStep<3>::GetVelocityDofList(DofsVectorType& rElementalDofList,
                                           ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 4;
    const SizeType LocalSize = 12;
    GeometryType& rGeom = this->GetGeometry();

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
    }
}

template< unsigned int TDim >
void FractionalStep<TDim>::GetPressureDofList(DofsVectorType& rElementalDofList,
                                              ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = TDim + 1;
    GeometryType& rGeom = this->GetGeometry();

    if (rElementalDofList.size() != NumNodes)
        rElementalDofList.resize(NumNodes);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
    }
}

template< unsigned int TDim >
void FractionalStep<TDim>::GetPressureValues(Vector& rValues,
                                            const int Step)
{
    const SizeType NumNodes = TDim + 1;
    GeometryType& rGeom = this->GetGeometry();

    if (rValues.size() != NumNodes) rValues.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; ++i)
        rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE,Step);
}

template<>
void FractionalStep<2>::GetVelocityValues(Vector& rValues,
                                          const int Step)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 6;

    GeometryType& rGeom = this->GetGeometry();

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
    }
}

template<>
void FractionalStep<3>::GetVelocityValues(Vector& rValues,
                                          const int Step)
{
    const SizeType NumNodes = 4;
    const SizeType LocalSize = 12;

    GeometryType& rGeom = this->GetGeometry();

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Z,Step);
    }
}


/*
 * protected FractionalStep<TDim> functions
 */

template< unsigned int TDim >
void FractionalStep<TDim>::CalculateGeometryData()
{
    const SizeType NumNodes = TDim + 1;
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rGeom.ShapeFunctionsLocalGradients( GeometryData::GI_GAUSS_1 );

    // Temporary container for inverse of J
    Matrix InvJ;

    GeometryType::JacobiansType J;
    rGeom.Jacobian( J, GeometryData::GI_GAUSS_1 );

    // calculate inverse of the jacobian and its determinant
    MathUtils<double>::InvertMatrix( J[0], InvJ, mDetJ );

    // calculate the shape function derivatives in global coordinates
    mDN_DX.resize(NumNodes,TDim);
    noalias( mDN_DX ) = prod( DN_De[0], InvJ );

//    // calculate nodal area, will be used to (approximately) compute projections
//    const ShapeFunctionsType CenterN = row(rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1),0);
//    const double Area = mDetJ * rGeom.IntegrationPoints(GeometryData::GI_GAUSS_1)[0].Weight();

//    // Carefully write results to nodal variables, to avoid parallelism problems
//    for (SizeType i = 0; i < NumNodes; ++i)
//    {
//        this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
//        this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += Area * CenterN[i];
//        this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
//    }

    // calculate minimum element length (used in stabilization Tau)
    array_1d<double,3> Edge(3,0.0);
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

//    // calculate minimum element height (for stabilization and Smagorinsky)
//    mElemSize = 0.0;
//    for ( unsigned int d = 0; d < TDim; ++d)
//    {
//        double hd = 1.0 / mDN_DX(1,d);
//        mElemSize += hd * hd;
//    }

//    for (unsigned int i = 1; i < NumNodes; ++i)
//    {
//        double Height = 0.0;
//        for ( unsigned int d = 0; d < TDim; ++d)
//        {
//            double hd = 1.0 / mDN_DX(i,d);
//            Height += hd * hd;
//        }
//        if (Height < mElemSize) mElemSize = Height;
//    }

//    mElemSize = sqrt(mElemSize);
}

template< unsigned int TDim >
void FractionalStep<TDim>::AddMomentumMassTerm(Matrix& rMassMatrix,
                                               const ShapeFunctionsType& rN,
                                               const double Weight)
{
    const SizeType NumNodes = TDim + 1;

    IndexType FirstRow = 0;
    IndexType FirstCol = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        for (SizeType j = 0; j < NumNodes; ++j)
        {
            const double Mij = Weight * rN[i] * rN[j];
            for (SizeType d =  0; d < TDim; ++d)
                    rMassMatrix(FirstRow+d,FirstCol+d) += Mij;
            FirstCol += TDim;
        }
        FirstRow += TDim;
        FirstCol = 0;
    }
}

template< unsigned int TDim >
void FractionalStep<TDim>::AddMomentumSystemTerms(Matrix& rLHSMatrix,
                                                  Vector& rRHSVector,
                                                  const double Density,
                                                  const array_1d<double,TDim+1>& rConvOperator,
                                                  const array_1d<double,3>& rBodyForce,
                                                  const double OldPressure,
                                                  const double TauOne,
                                                  const double TauTwo,
                                                  const array_1d<double,3>& rMomentumProjection,
                                                  const double MassProjection,
                                                  const ShapeFunctionsType& rN,
                                                  const ShapeFunctionDerivativesType& rDN_DX,
                                                  const double Weight)
{
    const SizeType NumNodes = TDim + 1;

    SizeType FirstRow = 0;
    SizeType FirstCol = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {

        // Build RHS
        for (SizeType d = 0; d < TDim; ++d)
        {
            // Body force
            double RHSi = Density * rN[i] * rBodyForce[d];
            // Pressure gradient (integrated by parts)
            RHSi += rDN_DX(i,d) * OldPressure;
            // Momentum Stabilization
            RHSi += Density * rConvOperator[i] * TauOne * ( /*rBodyForce[d] + rOldPressureGradient[d]*/ - rMomentumProjection[d] );
            // Mass Stabilization
            RHSi -= rDN_DX(i,d) * TauTwo * MassProjection;

            rRHSVector[FirstRow+d] += Weight * RHSi;
        }

        // Build LHS
        for (SizeType j = 0; j < NumNodes; ++j)
        {
            // Convective term
            double Kij = Density * rN[i] * rConvOperator[j];

            // Streamline stabilization
            Kij += Density * rConvOperator[i] * TauOne * Density * rConvOperator[j];

            Kij *= Weight;

            for (SizeType d = 0; d < TDim; ++d)
                rLHSMatrix(FirstRow + d,FirstCol + d) += Kij;

            // Mass-GLS (TauTwo) stabiliziation term
            for (SizeType m = 0; m < TDim; ++m)
                for (SizeType n = 0; n < TDim; ++n)
                    rLHSMatrix(FirstRow+m,FirstCol+n) += Weight * rDN_DX(i,m) * TauTwo * rDN_DX(j,n);

            FirstCol += TDim;
        }
        FirstRow += TDim;
        FirstCol = 0;
    }
}

template<>
void FractionalStep<2>::AddViscousTerm(MatrixType& rDampMatrix,
                                       const ShapeFunctionDerivativesType& rShapeDeriv,
                                       const double Weight)
{
    const SizeType NumNodes = 3;

    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    SizeType FirstRow(0),FirstCol(0);

    for (SizeType j = 0; j < NumNodes; ++j)
    {
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            // First Row
            rDampMatrix(FirstRow,FirstCol) += Weight * ( FourThirds * rShapeDeriv(i,0) * rShapeDeriv(j,0) + rShapeDeriv(i,1) * rShapeDeriv(j,1) );
            rDampMatrix(FirstRow,FirstCol+1) += Weight * ( nTwoThirds * rShapeDeriv(i,0) * rShapeDeriv(j,1) + rShapeDeriv(i,1) * rShapeDeriv(j,0) );

            // Second Row
            rDampMatrix(FirstRow+1,FirstCol) += Weight * ( nTwoThirds * rShapeDeriv(i,1) * rShapeDeriv(j,0) + rShapeDeriv(i,0) * rShapeDeriv(j,1) );
            rDampMatrix(FirstRow+1,FirstCol+1) += Weight * ( FourThirds * rShapeDeriv(i,1) * rShapeDeriv(j,1) + rShapeDeriv(i,0) * rShapeDeriv(j,0) );

            // Update Counter
            FirstRow += 2;
        }
        FirstRow = 0;
        FirstCol += 2;
    }
}

template <>
void FractionalStep<3>::AddViscousTerm(MatrixType& rDampMatrix,
                                       const ShapeFunctionDerivativesType& rShapeDeriv,
                                       const double Weight)
{
    const SizeType NumNodes = 4;

    const double OneThird = 1.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    unsigned int FirstRow(0),FirstCol(0);

    for (SizeType j = 0; j < NumNodes; ++j)
    {
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            // (dN_i/dx_k dN_j/dx_k)
            const double Diag =  rShapeDeriv(i,0) * rShapeDeriv(j,0) + rShapeDeriv(i,1) * rShapeDeriv(j,1) + rShapeDeriv(i,2) * rShapeDeriv(j,2);

            // First Row
            rDampMatrix(FirstRow,FirstCol) += Weight * ( OneThird * rShapeDeriv(i,0) * rShapeDeriv(j,0) + Diag );
            rDampMatrix(FirstRow,FirstCol+1) += Weight * ( nTwoThirds * rShapeDeriv(i,0) * rShapeDeriv(j,1) + rShapeDeriv(i,1) * rShapeDeriv(j,0) );
            rDampMatrix(FirstRow,FirstCol+2) += Weight * ( nTwoThirds * rShapeDeriv(i,0) * rShapeDeriv(j,2) + rShapeDeriv(i,2) * rShapeDeriv(j,0) );

            // Second Row
            rDampMatrix(FirstRow+1,FirstCol) += Weight * ( nTwoThirds * rShapeDeriv(i,1) * rShapeDeriv(j,0) + rShapeDeriv(i,0) * rShapeDeriv(j,1) );
            rDampMatrix(FirstRow+1,FirstCol+1) += Weight * ( OneThird * rShapeDeriv(i,1) * rShapeDeriv(j,1) + Diag );
            rDampMatrix(FirstRow+1,FirstCol+2) += Weight * ( nTwoThirds * rShapeDeriv(i,1) * rShapeDeriv(j,2) + rShapeDeriv(i,2) * rShapeDeriv(j,1) );

            // Third Row
            rDampMatrix(FirstRow+2,FirstCol) += Weight * ( nTwoThirds * rShapeDeriv(i,2) * rShapeDeriv(j,0) + rShapeDeriv(i,0) * rShapeDeriv(j,2) );
            rDampMatrix(FirstRow+2,FirstCol+1) += Weight * ( nTwoThirds * rShapeDeriv(i,2) * rShapeDeriv(j,1) + rShapeDeriv(i,1) * rShapeDeriv(j,2) );
            rDampMatrix(FirstRow+2,FirstCol+2) += Weight * ( OneThird * rShapeDeriv(i,2) * rShapeDeriv(j,2) + Diag );

            // Update Counter
            FirstRow += 3;
        }
        FirstRow = 0;
        FirstCol += 3;
    }
}

template <unsigned int TDim>
void FractionalStep<TDim>::CalculateTau(double &TauOne,
                                        double &TauTwo,
                                        const ProcessInfo &rCurrentProcessInfo)
{
    const ShapeFunctionsType& N = row( this->GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_1), 0);

    double Density = 0.0;
    double Viscosity = this->GetValue(VISCOSITY);
    array_1d<double,3> Velocity(3,0.0);
    array_1d<double,3> MeshVelocity(3,0.0);

    this->EvaluateInPoint(Density,DENSITY,N);
    this->EvaluateInPoint(Velocity,VELOCITY,N);
    this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);
    Velocity -= MeshVelocity;

    const double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
    const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);

    // Compute mean advective velocity norm
    double AdvVelNorm = 0.0;
    for (unsigned int d = 0; d < TDim; ++d)
        AdvVelNorm += Velocity[d] * Velocity[d];

    AdvVelNorm = sqrt(AdvVelNorm);

    TauOne = 1.0 / (Density * ( TimeFactor / DeltaTime + 4.0 * Viscosity / (mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );
    TauTwo = Density * (Viscosity + 0.5 * mElemSize * AdvVelNorm);
}

template <unsigned int TDim>
void FractionalStep<TDim>::CalculateTau(double& TauOne,
                                        double& TauTwo,
                                        const array_1d< double, 3 > & rAdvVel,
                                        const double Density,
                                        const double Viscosity,
                                        const ProcessInfo& rCurrentProcessInfo)
{
    const double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
    const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);

    // Compute mean advective velocity norm
    double AdvVelNorm = 0.0;
    for (unsigned int d = 0; d < TDim; ++d)
        AdvVelNorm += rAdvVel[d] * rAdvVel[d];

    AdvVelNorm = sqrt(AdvVelNorm);

    TauOne = 1.0 / (Density * ( TimeFactor / DeltaTime + 4.0 * Viscosity / (mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );
    TauTwo = Density * (Viscosity + 0.5 * mElemSize * AdvVelNorm);

}


template <unsigned int TDim>
void FractionalStep<TDim>::CalculateProjectionRHS(array_1d<double,3>& rMomentumRHS,
                                                  double& rMassRHS,
                                                  const ShapeFunctionsType& rN,
                                                  const ShapeFunctionDerivativesType& rDN_DX,
                                                  const double Weight)
{
    const SizeType NumNodes = TDim + 1;
    double Density;
    array_1d<double,3> BodyForce(3,0.0);
    array_1d<double,3> Velocity(3,0.0);
    array_1d<double,3> MeshVelocity(3,0.0);

    this->EvaluateInPoint(Density,DENSITY,rN);
    this->EvaluateInPoint(BodyForce,BODY_FORCE,rN);
    this->EvaluateInPoint(Velocity,VELOCITY,rN);
    this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,rN);

    array_1d<double,3> ConvVel = Velocity - MeshVelocity;
    Vector ConvOp(NumNodes);
    this->EvaluateConvection(ConvOp,ConvVel,rDN_DX);

    array_1d<double,3> Convection(3,0.0);
    for (SizeType i = 0; i < NumNodes; ++i)
        for (SizeType d = 0; d < TDim; ++d)
            Convection[d] += ConvOp[i] * this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[d];

    array_1d<double,TDim> PressureGradient(TDim,0.0);
    this->EvaluateGradientInPoint(PressureGradient,PRESSURE,rDN_DX);

    double Divergence;
    this->EvaluateDivergenceInPoint(Divergence,VELOCITY,rDN_DX);

    rMomentumRHS += Weight * ( Density * ( BodyForce - Convection) - PressureGradient );
    rMassRHS -= Weight * Divergence;
}


template <unsigned int TDim>
void FractionalStep<TDim>::CalculateProjectionRHS(VectorType& rConvTerm,
                                                  VectorType& rPresTerm,
                                                  VectorType& rDivTerm,
                                                  const ShapeFunctionsType& rN,
                                                  const ShapeFunctionDerivativesType& rDN_DX,
                                                  const double Weight)
{
    const unsigned int NumNodes = TDim + 1;
    double Density;
    array_1d<double,3> BodyForce(3,0.0);
    array_1d<double,3> Velocity(3,0.0);
    array_1d<double,3> MeshVelocity(3,0.0);

    this->EvaluateInPoint(Density,DENSITY,rN);
    this->EvaluateInPoint(BodyForce,BODY_FORCE,rN);
    this->EvaluateInPoint(Velocity,VELOCITY,rN);
    this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,rN);

    array_1d<double,3> ConvVel = Velocity - MeshVelocity;
    Vector ConvOp(NumNodes);
    this->EvaluateConvection(ConvOp,ConvVel,rDN_DX);

    array_1d<double,3> Convection(3,0.0);
    for (unsigned int i = 0; i < NumNodes; ++i)
        for (unsigned int d = 0; d < TDim; ++d)
            Convection[d] += ConvOp[i] * this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[d];

    array_1d<double,TDim> PressureGradient(TDim,0.0);
    this->EvaluateGradientInPoint(PressureGradient,PRESSURE,rDN_DX);

    double Divergence;
    this->EvaluateDivergenceInPoint(Divergence,VELOCITY,rDN_DX);

    unsigned int RowIndex = 0;

    for (unsigned int j = 0; j < NumNodes; ++j)
    {
        for (unsigned int d = 0; d < TDim; ++d)
        {
            rPresTerm[RowIndex] += Weight * rN[j] * ( Density * BodyForce[d] - PressureGradient[d] );
            rConvTerm[RowIndex] -= Weight * rN[j] * Density * Convection[d];
            ++RowIndex;
        }
        rDivTerm[j] -= Weight * rN[j] * Divergence;
    }
}

/*
 * private FractionalStep<TDim> functions
 */

template< unsigned int TDim >
void FractionalStep<TDim>::EvaluateConvection(Vector& rResult,
                                              const array_1d<double,3>& rConvVel,
                                              const ShapeFunctionDerivativesType& DN_DX)
{
    const SizeType NumNodes = TDim + 1;

    if(rResult.size() != NumNodes) rResult.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; i++)
    {
        rResult[i] = rConvVel[0]*DN_DX(i,0);
        for(SizeType k = 1; k < TDim; k++)
            rResult[i] += rConvVel[k]*DN_DX(i,k);
    }
}

//template<>
//FractionalStep<2>::ShapeFunctionsContainerType FractionalStep<2>::InitializeShapeFunctions()
//{
//    const double TwoThirds = 2.0 / 3.0;
//    const double OneSixth = 1.0 / 6.0;
//    ShapeFunctionsContainerType NContainer(3);
//    ShapeFunctionsType Temp(3,0.0);

//    // First Gauss point
//    Temp[0] = TwoThirds;
//    Temp[1] = OneSixth;
//    Temp[2] = OneSixth;
//    NContainer[0] = Temp;

//    // Second Gauss point
//    Temp[0] = OneSixth;
//    Temp[1] = TwoThirds;
//    Temp[2] = OneSixth;
//    NContainer[1] = Temp;

//    // Third Gauss point
//    Temp[0] = OneSixth;
//    Temp[1] = OneSixth;
//    Temp[2] = TwoThirds;
//    NContainer[2] = Temp;

//    return NContainer;
//}

///*
// * FractionalStep<TDim> static members
// */

//template< unsigned int TDim >
//static const FractionalStep<TDim>::ShapeFunctionsType FractionalStep<TDim>::msNg = FractionalStep<TDim>::InitializeShapeFunctions();

/*
 * Template class definition (this should allow us to compile the desired template instantiations)
 */

template class FractionalStep<2>;
template class FractionalStep<3>;

}
