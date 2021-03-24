#include "spalart_allmaras.h"
#include "utilities/math_utils.h"
#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{


Element::Pointer SpalartAllmaras::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SpalartAllmaras>(NewId, this->GetGeometry().Create(ThisNodes), pProperties, this->mIntegrationMethod);
}

Element::Pointer SpalartAllmaras::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SpalartAllmaras>(NewId, pGeom, pProperties);
}

int SpalartAllmaras::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY;

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if(ErrorCode != 0) return ErrorCode;

    // Checks on nodes
    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(VISCOSITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(MOLECULAR_VISCOSITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing MOLECULAR_VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(TURBULENT_VISCOSITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing TURBULENT_VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(TEMP_CONV_PROJ) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing TEMP_CONV_PROJ variable on solution step data for node ",this->GetGeometry()[i].Id());


        if(this->GetGeometry()[i].HasDofFor(TURBULENT_VISCOSITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing TURBULENT_VISCOSITY degree of freedom on node ",this->GetGeometry()[i].Id());
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

    return 0;

    KRATOS_CATCH("");
}

void SpalartAllmaras::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    this->InitializeElementData();

    KRATOS_CATCH( "" )
}

void SpalartAllmaras::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

    //if (FractionalStepNumber == 1)
    //{
    this->InitializeElementData();
    /*}
    else */
    if (FractionalStepNumber == 2)
    {
        const unsigned int NumNodes = this->GetGeometry().PointsNumber();
        const unsigned int NumGauss = this->GetGeometry().IntegrationPoints(this->mIntegrationMethod).size();

        const Matrix NContainer = this->GetGeometry().ShapeFunctionsValues(this->mIntegrationMethod);
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints( this->mIntegrationMethod );

        VectorType NodalArea = ZeroVector(NumNodes);
        VectorType ConvTerm = ZeroVector(NumNodes);

        for (unsigned int g = 0; g < NumGauss; g++) // Loop on Gauss points
        {
            const ShapeFunctionsType& N = row(NContainer,g);
            const ShapeDerivativesType& DN_DX = mDN_DX[g];
            const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();

            // Evaluate convective velocity on integration point
            array_1d<double,3> Velocity = ZeroVector(3);
            array_1d<double,3> MeshVelocity = ZeroVector(3);
            this->EvaluateInPoint(Velocity,VELOCITY,N);
            this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);

            // For ALE: convective velocity
            array_1d<double,3> ConvVel = Velocity - MeshVelocity;

            // Evaluate convection operator Velocity * Grad(N)
            Vector UGradN(NumNodes);
            this->EvaluateConvection(UGradN,ConvVel,DN_DX);

            // Finish evaluation of convective term as Velocity * Grad(N_i) * Viscosity_node_i
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                double Temp = GaussWeight * UGradN[i] * this->GetGeometry()[i].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
                for (unsigned int j = 0; j < NumNodes; j++)
                {
                    ConvTerm[j] += N[j] * Temp;
                    NodalArea[j] += N[j] * GaussWeight;
                }
            }
        }

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
            this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
            this->GetGeometry()[i].FastGetSolutionStepValue(TEMP_CONV_PROJ) += ConvTerm[i];
            this->GetGeometry()[i].UnSetLock();
        }
    }
    KRATOS_CATCH("")
}

void SpalartAllmaras::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
{
    // Obtain required constants
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const SizeType NumGauss = this->GetGeometry().IntegrationPoints(this->mIntegrationMethod).size();

    const Matrix NContainer = this->GetGeometry().ShapeFunctionsValues(this->mIntegrationMethod);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints( this->mIntegrationMethod );

    const double ElementSize = this->ElementSize();
    const double Tau = this->CalculateTau(ElementSize,rCurrentProcessInfo);

    // Choose between DES or URANS operation
    bool UseRealDistance = true;
    double Distance = 0.0;
    const ProcessInfo& rConstProcInfo = static_cast<const ProcessInfo&>(rCurrentProcessInfo);
    double Cdes = rConstProcInfo[C_DES];
    if ( Cdes != 0.0)
    {
        Distance = this->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE);
        for (SizeType i = 1; i < NumNodes; i++)
            Distance += this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        Distance /= double(NumNodes);

        double DesSize = Cdes * ElementSize;
        if (Distance > DesSize)
        {
            Distance = DesSize;
            UseRealDistance = false;
        }
    }


    // Initialize local contribution
    if (rLeftHandSideMatrix.size1() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

    rLeftHandSideMatrix = ZeroMatrix(NumNodes,NumNodes);

    if (rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes, false);

    rRightHandSideVector = ZeroVector(NumNodes);

    MatrixType MassMatrix = ZeroMatrix(NumNodes,NumNodes);

    // Loop on integration points
    for (SizeType g = 0; g < NumGauss; g++)
    {
        const ShapeFunctionsType& N = row(NContainer,g);
        const ShapeDerivativesType& DN_DX = mDN_DX[g];
        const double GaussWeight = mDetJ * IntegrationPoints[g].Weight();

        // Add contribution to mass matrix
        this->AddMassTerm(MassMatrix,N,GaussWeight);

        double MolecularViscosity;
        double LastViscosity;

        array_1d<double,3> Velocity = ZeroVector(3);
        array_1d<double,3> MeshVelocity = ZeroVector(3);

        this->EvaluateInPoint(MolecularViscosity,MOLECULAR_VISCOSITY,N);
        this->EvaluateInPoint(LastViscosity,TURBULENT_VISCOSITY,N);
        this->EvaluateInPoint(Velocity,VELOCITY,N);
        this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);

        // If not using DES (or if using DES but close to the wall
        if (UseRealDistance)
            this->EvaluateInPoint(Distance,DISTANCE,N);

        // For ALE: convective velocity
        array_1d<double,3> ConvVel = Velocity - MeshVelocity;

        // Evaluate convection operator Velocity * Grad(N)
        Vector UGradN(NumNodes);
        this->EvaluateConvection(UGradN,ConvVel,DN_DX);

        // Evaluate the gradient of the transported variable on the integration point
        // (used in the linearisation of the non-conservative diffusion term)
        array_1d<double,3> LastViscosityGradient = ZeroVector(3);
        for (SizeType i = 0; i < NumNodes; i++)
        {
            double NodalViscosity = this->GetGeometry()[i].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            for (SizeType d = 0; d < Dim; d++)
                LastViscosityGradient[d] += DN_DX(i,d) * NodalViscosity;
        }

        // Add convective term
        this->AddConvection(rLeftHandSideMatrix,N,UGradN,GaussWeight);

        // Add production, diffusion and destruction of eddy viscosity
        this->AddModelTerms(rLeftHandSideMatrix,MolecularViscosity,LastViscosity,LastViscosityGradient,Distance,Velocity,N,DN_DX,GaussWeight);

        // Add stabilization
        this->AddStabilization(rLeftHandSideMatrix,rRightHandSideVector,Tau,N,DN_DX,UGradN,GaussWeight);
    }

    // Add residual of previous iteration to RHS
    VectorType LastValues = ZeroVector(NumNodes);
    this->GetValuesVector(LastValues);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,LastValues);

    // Add dynamic term
    const Vector& rBDFCoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
    noalias( rLeftHandSideMatrix ) += rBDFCoeffs[0] * MassMatrix;

    VectorType TimeTerm = rBDFCoeffs[0] * LastValues;
    for(SizeType i = 1; i < rBDFCoeffs.size(); i++)
    {
        this->GetValuesVector(LastValues,i);
        noalias( TimeTerm ) += rBDFCoeffs[i] * LastValues;
    }

    noalias( rRightHandSideVector ) -= prod(MassMatrix,TimeTerm);
}

void SpalartAllmaras::CalculateRightHandSide(VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
{
    MatrixType TempMatrix;
    this->CalculateLocalSystem(TempMatrix,rRightHandSideVector,rCurrentProcessInfo);
}

void SpalartAllmaras::GetDofList(DofsVectorType &rElementalDofList, const ProcessInfo &rCurrentProcessInfo) const
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    if (rElementalDofList.size() != NumNodes)
        rElementalDofList.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(TURBULENT_VISCOSITY);
}


void SpalartAllmaras::EquationIdVector(Element::EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    if (rResult.size() != NumNodes)
        rResult.resize(NumNodes, false);

    for (SizeType i = 0; i < NumNodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(TURBULENT_VISCOSITY).EquationId();
}

void SpalartAllmaras::GetValuesVector(Vector &rValues, int Step) const
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    if(rValues.size() != NumNodes)
        rValues.resize(NumNodes);

    for(SizeType i = 0; i < NumNodes; i++)
    {
        rValues[i] = this->GetGeometry()[i].FastGetSolutionStepValue(TURBULENT_VISCOSITY,Step);
    }
}

// Protected functions: Definition of local contributions ---------------------

void SpalartAllmaras::InitializeElementData()
{
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints( this->mIntegrationMethod );

    // Initialize member variables
    mDN_DX.resize( IntegrationPoints.size() ); // Shape function derivatives container
    mDetJ = 0.00;

    GeometryType::JacobiansType J;
    J = GetGeometry().Jacobian( J, mIntegrationMethod );

    const GeometryType::ShapeFunctionsGradientsType& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients( this->mIntegrationMethod );

    // container for inverse of J
    Matrix InvJ;

    //calculating the inverse J
    for ( SizeType g = 0; g < IntegrationPoints.size(); g++ )
    {
        //calculating and storing inverse of the jacobian and the parameters needed
        MathUtils<double>::InvertMatrix( J[g], InvJ, mDetJ );

        //calculating the shape function derivatives in global coordinates
        mDN_DX[g].resize(NumNodes,Dim);
        noalias( mDN_DX[g] ) = prod( DN_De[g], InvJ );
    }
}

// Lumped version
void SpalartAllmaras::AddMassTerm(MatrixType &rMassMatrix,
                                  const ShapeFunctionsType &N,
                                  const double Weight)
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double LumpFactor = Weight / double(NumNodes);
    for(SizeType i = 0; i < NumNodes; i++)
        rMassMatrix(i,i) += LumpFactor;
}

//void SpalartAllmaras::AddMassTerm(MatrixType &rMassMatrix,
//                                  const ShapeFunctionsType N,
//                                  const double Weight)
//{
//    const SizeType NumNodes = this->GetGeometry().PointsNumber();
//    for(SizeType i = 0; i < NumNodes; i++)
//        for(SizeType j = 0; j < NumNodes; j++)
//            rMassMatrix(i,j) += Weight * N[i] * N[j];
//}

void SpalartAllmaras::AddConvection(MatrixType &rLHS,
                                    const ShapeFunctionsType &N,
                                    const Vector &UGradN,
                                    const double Weight)
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    for (SizeType j = 0; j < NumNodes; j++)
    {
        for (SizeType i = 0; i < NumNodes; i++)
        {
            rLHS(i,j) += Weight * N[i] * UGradN[j];
        }
    }
}

void SpalartAllmaras::AddModelTerms(MatrixType &rLHS,
                                    const double MolecularViscosity,
                                    const double LastEddyViscosity,
                                    const array_1d<double,3> rLastEddyViscosityGradient,
                                    const double Distance,
                                    const array_1d<double, 3> &rVelocity,
                                    const ShapeFunctionsType &N,
                                    const ShapeDerivativesType &DN_DX,
                                    const double Weight)
{
    // Constants of the Spalart-Allmaras model
    const double sigma = 2.0 / 3.0; // Prandtl number for problem variable
    //const double kappa = 0.41; // Von Karman's constant
    const double kappa_2 = 0.41*0.41;
    const double cb1 = 0.1355; // Production coefficient
    const double cb2 = 0.6220; // Coefficient for non-consistent diffusion
    const double cw1 = cb1 / (kappa_2) + (1.0 + cb2) / sigma; // Destruction coefficient
    const double cw2 = 0.3; // Used to compute fw in destruction term
    //const double cw3 = 2.0; // Used to compute fw in destruction term
    const double cw3_6 = 64.0; // Sixth power of cw3

    // Constants for variable transformation (eddy viscosity -> transported variable)
    // const double cv1 = 7.1;
    const double cv1_3 = 7.1 * 7.1 * 7.1;
    const double Xi = LastEddyViscosity / MolecularViscosity;
    const double Xi_3 = Xi*Xi*Xi;
    const double fv1 = Xi_3 / (Xi_3 + cv1_3);
    const double fv2 = 1.0 - Xi / (1.0 + Xi * fv1);

    // Modified strain rate (using rotation correction)
    const double Cprod = 2.0;
    double NormS = 0.0;
    double NormOmega = 0.0;
    this->VelocityGradientNorms(NormS,NormOmega,DN_DX);
    // S = NormOmega + Cprod * min(0,NormS-NormOmega);
    double S = NormOmega;
    if(NormS < NormOmega)
        S += Cprod * (NormS-NormOmega);
    double S_hat = S + ( fv2 * LastEddyViscosity / (kappa_2 * Distance*Distance) );
    // Numerical control on S_hat
    if (S_hat < 0.3 * NormOmega)
        S_hat = 0.3 * NormOmega;

    // Destruction function
    double r = LastEddyViscosity / (S_hat * kappa_2 * Distance*Distance);
    if (r > 10.0) r = 10.0; // Numerical control in r
    const double g = r + cw2 * (pow(r,6) - r);
    const double fw = g * pow( (1.0+cw3_6) / ( pow(g,6) + cw3_6 ) , 1.0/6.0);

    // Geometric constants
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    // Auxiliary values
    const double ProductionTerm = Weight * cb1 * S_hat;
    const double DestructionTerm = Weight * cw1 * fw * LastEddyViscosity / (Distance * Distance);
    const double Diffusivity = (MolecularViscosity + LastEddyViscosity) / sigma;

    // Add terms to matrix
    for (SizeType j = 0; j < NumNodes; j++)
    {
        // Compute intermediate value for the (linearisation of)
        // the second diffusive term cb2/sigma * Gradient(Viscosity)^2
        double tmp = rLastEddyViscosityGradient[0] * DN_DX(j,0); // Store Grad(LastViscosity):Grad(Nj)
        for (SizeType d = 1; d < Dim; d++)
            tmp += rLastEddyViscosityGradient[d] * DN_DX(j,d);
        tmp *= Weight * cb2 / sigma;

        for (SizeType i = 0; i < NumNodes; i++)
        {
            // Production and destruction
            rLHS(i,j) += (DestructionTerm - ProductionTerm) * N[i] * N[j];

            // Diffusion
            for (SizeType d = 0; d < Dim; d++)
                rLHS(i,j) += Weight * Diffusivity * DN_DX(i,d) * DN_DX(j,d);

            // Add the (linearisation of) the second diffusive term cb2/sigma * Gradient(Viscosity)^2
            rLHS(i,j) -= N[i] * tmp;
        }
    }
}

void SpalartAllmaras::AddStabilization(MatrixType &rLHS,
                                       VectorType &rRHS,
                                       const double Tau,
                                       const ShapeFunctionsType &N,
                                       const ShapeDerivativesType &DN_DX,
                                       const Vector &UGradN,
                                       const double Weight)
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    for (SizeType j = 0; j < NumNodes; j++)
    {
        const double NodalProjection = this->GetGeometry()[j].FastGetSolutionStepValue(TEMP_CONV_PROJ);
        for (SizeType i = 0; i < NumNodes; i++)
        {
            rLHS(i,j) += Weight * Tau * UGradN[i] * UGradN[j];
            rRHS[i] -= Weight * Tau * UGradN[i] * N[j] * NodalProjection;
        }
    }
}

void SpalartAllmaras::EvaluateConvection(Vector &rResult,
        const array_1d<double, 3> &rConvVel,
        const ShapeDerivativesType &DN_DX)
{
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    if(rResult.size() != NumNodes) rResult.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; i++)
    {
        rResult[i] = rConvVel[0]*DN_DX(i,0);
        for(SizeType k = 1; k < Dim; k++)
            rResult[i] += rConvVel[k]*DN_DX(i,k);
    }
}

double SpalartAllmaras::MeasureOfVorticity(const ShapeDerivativesType &DN_DX)
{
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    Matrix Omega = ZeroMatrix(Dim,Dim); // Rate of rotation matrix
    for (SizeType n = 0; n < NumNodes; n++)
    {
        const array_1d<double,3>& rNodalVelocity = this->GetGeometry()[n].FastGetSolutionStepValue(VELOCITY);
        for (SizeType i = 0; i < Dim; i++)
            for (SizeType j = 0; j < Dim; j++)
                Omega(i,j) += 0.5 * ( DN_DX(n,j)*rNodalVelocity[i] - DN_DX(n,i)*rNodalVelocity[j] );
    }

    // Compute norm of Omega
    double S = 0.0;
    for (SizeType i = 0; i < Dim; i++)
        for (SizeType j = 0; j < Dim; j++)
            S += 2.0 * Omega(i,j) * Omega(i,j);
    S = std::sqrt(S);
    return S;
}

void SpalartAllmaras::VelocityGradientNorms(double &rNormS,
        double &rNormOmega,
        const ShapeDerivativesType &DN_DX)
{
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    Matrix S = ZeroMatrix(Dim,Dim); // Rate of strain
    Matrix Omega = ZeroMatrix(Dim,Dim); // Rate of rotation

    for (SizeType n = 0; n < NumNodes; n++)
    {
        const array_1d<double,3>& rNodalVelocity = this->GetGeometry()[n].FastGetSolutionStepValue(VELOCITY);
        for (SizeType i = 0; i < Dim; i++)
            for (SizeType j = 0; j < Dim; j++)
            {
                S(i,j) += 0.5 * ( DN_DX(n,j)*rNodalVelocity[i] + DN_DX(n,i)*rNodalVelocity[j] );
                Omega(i,j) += 0.5 * ( DN_DX(n,j)*rNodalVelocity[i] - DN_DX(n,i)*rNodalVelocity[j] );
            }
    }

    // Compute norms
    rNormS = 0.0;
    rNormOmega = 0.0;

    for (SizeType i = 0; i < Dim; i++)
        for (SizeType j = 0; j < Dim; j++)
        {
            rNormS += 2.0 * S(i,j) * S(i,j);
            rNormOmega += 2.0 * Omega(i,j) * Omega(i,j);
        }
    rNormS = std::sqrt(rNormS);
    rNormOmega = std::sqrt(rNormOmega);
}

double SpalartAllmaras::CalculateTau(double ElementSize, const ProcessInfo &rCurrentProcessInfo)
{
    // Geometry
    const Element::GeometryType& rGeom = this->GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();

    // Shape functions at the element center
    const Matrix NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
    const Vector& N = row(NContainer,0);

    // Time Term
    const double TimeCoeff = rCurrentProcessInfo[BDF_COEFFICIENTS][0];

    // Diffusivity
    double MolecularViscosity = 0.0;
    double EddyViscosity = 0.0;
    this->EvaluateInPoint(MolecularViscosity,MOLECULAR_VISCOSITY,N);
    this->EvaluateInPoint(EddyViscosity,TURBULENT_VISCOSITY,N);
    const double sigma = 2.0/3.0;
    double Diffusivity = (MolecularViscosity + EddyViscosity) / sigma;

    // Convective velocity
    array_1d<double,3> Velocity = ZeroVector(3);
    array_1d<double,3> MeshVelocity = ZeroVector(3);
    this->EvaluateInPoint(Velocity,VELOCITY,N);
    this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);
    Velocity -= MeshVelocity;
    double VelNorm = Velocity[0]*Velocity[0];
    for (SizeType d = 1; d < Dim; d++)
        VelNorm += Velocity[d]*Velocity[d];
    VelNorm = sqrt(VelNorm);

    double Tmp = TimeCoeff + 4.0 * Diffusivity / (ElementSize*ElementSize) + 2.0 * VelNorm / ElementSize;
    return 1.0 / Tmp;
}

double SpalartAllmaras::ElementSize()
{
    const Element::GeometryType& rGeom = this->GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumNodes = rGeom.PointsNumber();

    // Maximum edge length
    array_1d<double,3> Edge = rGeom[1].Coordinates() - rGeom[0].Coordinates();
    double MaxLength = Edge[0]*Edge[0];
    for (SizeType d = 1; d < Dim; d++)
        MaxLength += Edge[d]*Edge[d];

    for (SizeType i = 2; i < NumNodes; i++)
        for(SizeType j = 0; j < i; j++)
        {
            Edge = rGeom[i].Coordinates() - rGeom[j].Coordinates();
            double Length = Edge[0]*Edge[0];
            for (SizeType d = 1; d < Dim; d++)
                Length += Edge[d]*Edge[d];
            if (Length > MaxLength) MaxLength = Length;
        }
    return sqrt(MaxLength);
}

}
