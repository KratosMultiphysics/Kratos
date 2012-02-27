#include "spalart_allmaras.h"
#include "utilities/math_utils.h"

namespace Kratos {


Element::Pointer SpalartAllmaras::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SpalartAllmaras(NewId, this->GetGeometry().Create(ThisNodes), pProperties, this->mIntegrationMethod));
}

int SpalartAllmaras::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if(ErrorCode != 0) return ErrorCode;

    // Check that all required variables have been registered
    if(VELOCITY.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","")
    if(MESH_VELOCITY.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check if the application was correctly registered.","")
    if(VISCOSITY.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"VISCOSITY Key is 0. Check if the application was correctly registered.","")
    if(MOLECULAR_VISCOSITY.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"MOLECULAR_VISCOSITY Key is 0. Check if the application was correctly registered.","")
    if(TURBULENT_VISCOSITY.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"TURBULENT_VISCOSITY Key is 0. Check if the application was correctly registered.","")
    if(TEMP_CONV_PROJ.Key() == 0)
        KRATOS_ERROR(std::invalid_argument,"TEMP_CONV_PROJ Key is 0. Check if the application was correctly registered.","")

    // Checks on nodes

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
            KRATOS_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
            KRATOS_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(VISCOSITY) == false)
            KRATOS_ERROR(std::invalid_argument,"missing VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(MOLECULAR_VISCOSITY) == false)
            KRATOS_ERROR(std::invalid_argument,"missing MOLECULAR_VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(TURBULENT_VISCOSITY) == false)
            KRATOS_ERROR(std::invalid_argument,"missing TURBULENT_VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(TEMP_CONV_PROJ) == false)
            KRATOS_ERROR(std::invalid_argument,"missing TEMP_CONV_PROJ variable on solution step data for node ",this->GetGeometry()[i].Id());


        if(this->GetGeometry()[i].HasDofFor(TURBULENT_VISCOSITY) == false)
            KRATOS_ERROR(std::invalid_argument,"missing TURBULENT_VISCOSITY degree of freedom on node ",this->GetGeometry()[i].Id());
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

    return 0;

    KRATOS_CATCH("");
}

void SpalartAllmaras::Initialize()
{
    KRATOS_TRY;

    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints( this->mIntegrationMethod );

    // Initialize member variables
    mDN_DX.resize( IntegrationPoints.size() ); // Shape function derivatives container
    mElementSize = 0.00;

    GeometryType::JacobiansType J;
    J = GetGeometry().Jacobian( J, mIntegrationMethod );

    const GeometryType::ShapeFunctionsGradientsType& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients( this->mIntegrationMethod );

    // containers for inverse and determinant of J
    double DetJ; // Jacobian
    Matrix InvJ;

    //calculating the inverse J
    for ( SizeType g = 0; g < IntegrationPoints.size(); g++ )
    {
        //getting informations for integration
        double IntegrationWeight = IntegrationPoints[g].Weight();

        //calculating and storing inverse of the jacobian and the parameters needed
        MathUtils<double>::InvertMatrix( J[g], InvJ, DetJ );

        //calculating the total area
        mElementSize += DetJ * IntegrationWeight;

        //calculating the shape function derivatives in global coordinates
        mDN_DX[g].resize(NumNodes,Dim);
        noalias( mDN_DX[g] ) = prod( DN_De[g], InvJ );
    }

    KRATOS_CATCH( "" )
}

void SpalartAllmaras::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY
    int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

    if (FractionalStepNumber == 2)
    {
        const SizeType NumNodes = this->GetGeometry().PointsNumber();
        const SizeType NumGauss = this->GetGeometry().IntegrationPoints(this->mIntegrationMethod).size();

        const Matrix NContainer = this->GetGeometry().ShapeFunctionsValues(this->mIntegrationMethod);
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints( this->mIntegrationMethod );

        double NodalArea = 0.0;
        double ConvTerm = 0.0;

        for (SizeType g = 0; g < NumGauss; g++) // Loop on Gauss points
        {
            const ShapeFunctionsType& N = row(NContainer,g);
            const ShapeDerivativesType& DN_DX = mDN_DX[g];
            const double GaussWeight = mElementSize * IntegrationPoints[g].Weight();

            // Evaluate convective velocity on integration point
            array_1d<double,3> Velocity(3,0.0);
            array_1d<double,3> MeshVelocity(3,0.0);
            this->EvaluateInPoint(Velocity,VELOCITY,N);
            this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);

            // For ALE: convective velocity
            array_1d<double,3> ConvVel = Velocity - MeshVelocity;

            // Evaluate convection operator Velocity * Grad(N)
            Vector UGradN(NumNodes);
            this->EvaluateConvection(UGradN,ConvVel,DN_DX);

            // Finish evaluation of convective term as Velocity * Grad(N_i) * Viscosity_node_i
            NodalArea += GaussWeight;
            for (SizeType i = 0; i < NumNodes; i++)
                ConvTerm += GaussWeight * UGradN[i] * this->GetGeometry()[i].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        }

        for (SizeType i = 0; i < NumNodes; i++)
        {
            this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
            this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea;
            this->GetGeometry()[i].FastGetSolutionStepValue(TEMP_CONV_PROJ) += ConvTerm;
            this->GetGeometry()[i].UnSetLock();
        }
    }
    KRATOS_CATCH("")
}

void SpalartAllmaras::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{
    // Obtain required constants
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const SizeType NumGauss = this->GetGeometry().IntegrationPoints(this->mIntegrationMethod).size();

    const Matrix NContainer = this->GetGeometry().ShapeFunctionsValues(this->mIntegrationMethod);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints( this->mIntegrationMethod );

    const double Tau = this->CalculateTau(rCurrentProcessInfo);

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
        const double GaussWeight = mElementSize * IntegrationPoints[g].Weight();

        // Add contribution to mass matrix
        this->AddMassTerm(MassMatrix,N,GaussWeight);

        double MolecularViscosity;
        double LastViscosity;
        double Distance; // Distance to nearest wall
        array_1d<double,3> Velocity(3,0.0);
        array_1d<double,3> MeshVelocity(3,0.0);

        this->EvaluateInPoint(MolecularViscosity,MOLECULAR_VISCOSITY,N);
        this->EvaluateInPoint(LastViscosity,TURBULENT_VISCOSITY,N);
        this->EvaluateInPoint(Distance,DISTANCE,N);
        this->EvaluateInPoint(Velocity,VELOCITY,N);
        this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);

        // For ALE: convective velocity
        array_1d<double,3> ConvVel = Velocity - MeshVelocity;

        // Evaluate convection operator Velocity * Grad(N)
        Vector UGradN(NumNodes);
        this->EvaluateConvection(UGradN,ConvVel,DN_DX);

        // Evaluate the gradient of the transported variable on the integration point
        // (used in the linearisation of the non-conservative diffusion term)
        array_1d<double,3> LastViscosityGradient(3,0.0);
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


void SpalartAllmaras::GetDofList(DofsVectorType &rElementalDofList, ProcessInfo &rCurrentProcessInfo)
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    if (rElementalDofList.size() != NumNodes)
        rElementalDofList.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(TURBULENT_VISCOSITY);
}


void SpalartAllmaras::EquationIdVector(Element::EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    if (rResult.size() != NumNodes)
        rResult.resize(NumNodes, false);

    for (SizeType i = 0; i < NumNodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(TURBULENT_VISCOSITY).EquationId();
}

void SpalartAllmaras::GetValuesVector(Vector &rValues, int Step)
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

// Lumped version
void SpalartAllmaras::AddMassTerm(MatrixType &rMassMatrix,
                                  const ShapeFunctionsType N,
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
    const double sigma = 2.0 / 3.0; // Prandtl number
    const double kappa = 0.41; // Von Karman's constant
    const double cb1 = 0.1355; // Production coefficient
    const double cb2 = 0.6220; // Coefficient for non-consistent diffusion
    const double cw1 = cb1 / (kappa*kappa) + (1.0 + cb2) / sigma; // Destruction coefficient
    const double cw2 = 0.3; // Used to compute fw in destruction term
    const double cw3 = 2.0; // Used to compute fw in destruction term

    // Constants for variable transformation (eddy viscosity -> transported variable)
    const double cv1_cube = 7.1 * 7.1 * 7.1;
    const double Xi = LastEddyViscosity / MolecularViscosity;
    const double fv1 = (Xi * Xi * Xi) / (Xi * Xi * Xi + cv1_cube);
    const double fv2 = 1.0 - Xi / (1.0 + Xi * fv1);

//    // Modified strain rate
//    const double S = this->MeasureOfVorticity(DN_DX);
//    const double S_hat = S + ( fv2 * LastEddyViscosity / (kappa*kappa * Distance*Distance) );
    // Modified strain rate (using rotation correction)
    const double Cprod = 2.0;
    double NormS = 0.0;
    double NormOmega = 0.0;
    this->VelocityGradientNorms(NormS,NormOmega,DN_DX);
    // S = NormOmega + Cprod * min(0,NormS-NormOmega);
    double S = NormOmega;
    if(NormS < NormOmega)
        S += Cprod * (NormS-NormOmega);
    double S_hat = S + ( fv2 * LastEddyViscosity / (kappa*kappa * Distance*Distance) );
    // Numerical control on S_hat
    if (S_hat < 0.3 * NormOmega)
        S_hat = 0.3 * NormOmega;

    // Destruction function
    double r = LastEddyViscosity / (S_hat * kappa*kappa * Distance*Distance);
    if (r > 10.0) r = 10.0; // Numerical control in r
    const double g = r + cw2 * (pow(r,6) - r);
    const double fw = g * pow( (1.0+pow(cw3,6)) / ( pow(g,6) + pow(cw3,6) ) , 1.0/6.0);

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
            rRHS[i] += Weight * Tau * UGradN[i] * N[j] * NodalProjection;
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

double SpalartAllmaras::CalculateTau(const ProcessInfo &rCurrentProcessInfo)
{
    // Geometry
    const Element::GeometryType& rGeom = this->GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumNodes = rGeom.PointsNumber();

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
    array_1d<double,3> Velocity(3,0.0);
    array_1d<double,3> MeshVelocity(3,0.0);
    this->EvaluateInPoint(Velocity,VELOCITY,N);
    this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);
    Velocity -= MeshVelocity;
    double VelNorm = Velocity[0]*Velocity[0];
    for (SizeType d = 1; d < Dim; d++)
        VelNorm += Velocity[d]*Velocity[d];
    VelNorm = sqrt(VelNorm);

    // Minimum element length
    array_1d<double,3> Edge(3,0.0);

    Edge = rGeom[1].Coordinates() - rGeom[0].Coordinates();
    double MinLength = Edge[0]*Edge[0];
    for (SizeType d = 0; d < Dim; d++)
        MinLength += Edge[d]*Edge[d];
    MinLength = sqrt(MinLength);

    for (SizeType i = 2; i < NumNodes; i++)
        for(SizeType j = 0; j < i; j++)
        {
            Edge = rGeom[i].Coordinates() - rGeom[j].Coordinates();
            double Length = Edge[0]*Edge[0];
            for (SizeType d = 0; d < Dim; d++)
                Length += Edge[d]*Edge[d];
            if (Length < MinLength) MinLength = Length;
        }
    MinLength = sqrt(MinLength);

    double Tmp = TimeCoeff + 4.0 * Diffusivity / (MinLength*MinLength) + 2.0 * VelNorm / MinLength;
    return 1.0 / Tmp;
}

}
