#include "stationary_stokes.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"

namespace Kratos {

template< unsigned int TDim >
Element::Pointer StationaryStokes<TDim>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared< StationaryStokes<TDim> >(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template< unsigned int TDim >
Element::Pointer StationaryStokes<TDim>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_shared<StationaryStokes>(NewId, pGeom, pProperties);
}

template< unsigned int TDim >
int StationaryStokes<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    int Result = Element::Check(rCurrentProcessInfo);

    // Check that all required variables have been registered
    if(VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","");
    if(PRESSURE.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check if the application was correctly registered.","");
    if(DENSITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check if the application was correctly registered.","");
    if(VISCOSITY.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"VISCOSITY Key is 0. Check if the application was correctly registered.","");
    if(BODY_FORCE.Key() == 0)
        KRATOS_THROW_ERROR(std::invalid_argument,"BODY_FORCE Key is 0. Check if the application was correctly registered.","");

    // Checks on nodes

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                ( this->GetGeometry().WorkingSpaceDimension() == 3 && this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false ) )
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

    return Result;
}

template< unsigned int TDim >
void StationaryStokes<TDim>::Initialize()
{
    KRATOS_TRY;

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    this->mIntegrationMethod = GeometryData::GI_GAUSS_2;

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( this->mIntegrationMethod );

    // Initialize member variables
    mDN_DX.resize( IntegrationPoints.size() ); // Shape function derivatives container
    mGaussWeight.resize( IntegrationPoints.size() ); // Integration weigths at each integration point

    GeometryType::JacobiansType J;
    J = GetGeometry().Jacobian( J, mIntegrationMethod );

    const GeometryType::ShapeFunctionsGradientsType& DNv_De = rGeom.ShapeFunctionsLocalGradients( this->mIntegrationMethod );

    // Temporary container for inverse of J
    Matrix InvJ;
    double DetJ = 0.0;

    //calculating the inverse J
    for ( unsigned int g = 0; g < IntegrationPoints.size(); g++ )
    {
        //calculating and storing inverse of the jacobian and the parameters needed
        MathUtils<double>::InvertMatrix( J[g], InvJ, DetJ );

        //calculating the shape function derivatives in global coordinates
        mDN_DX[g].resize(NumNodes,TDim, false);
        noalias( mDN_DX[g] ) = prod( DNv_De[g], InvJ );

        // Gauss point weight is stored as a fraction of the elemental area
        mGaussWeight[g] = DetJ * rGeom.IntegrationPoints(this->mIntegrationMethod)[g].Weight();
    }

    KRATOS_CATCH( "" )
}

template< unsigned int TDim >
void StationaryStokes<TDim>::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{
    // Obtain required constants
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int BlockSize = TDim + 1;
    const unsigned int LocalSize = BlockSize * NumNodes;

    const unsigned int NumGauss = rGeom.IntegrationPointsNumber();

    const Matrix NContainer = rGeom.ShapeFunctionsValues(this->mIntegrationMethod);

    // Initialize local contribution
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    rRightHandSideVector = ZeroVector(LocalSize);

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const ShapeFunctionsType& N = row(NContainer,g);
        const ShapeDerivativesType& DN_DX = mDN_DX[g];
        const double GaussWeight = mGaussWeight[g];

        double Density;
        double Viscosity;
        array_1d<double,3> BodyForce = ZeroVector(3);
        array_1d<double,3> Velocity = ZeroVector(3);

        // Interpolation
        this->EvaluateInPoint(Density,DENSITY,N,rGeom);
        this->EvaluateInPoint(Viscosity,VISCOSITY,N,rGeom);
        this->EvaluateInPoint(BodyForce,BODY_FORCE,N,rGeom);

        this->EvaluateInPoint(Velocity,VELOCITY,N,rGeom);

        double TauOne = 0.0;
        double TauTwo = 0.0;
        CalculateTau(TauOne,TauTwo,Density*Viscosity);

        // Add velocity terms in momentum equation
        this->AddMomentumTerms(rLeftHandSideMatrix,rRightHandSideVector,Density,Viscosity,BodyForce,TauTwo,N,DN_DX,GaussWeight);

        // Add velocity-pressure terms
        this->AddContinuityTerms(rLeftHandSideMatrix,rRightHandSideVector,Density,BodyForce,TauOne,N,DN_DX,GaussWeight);
    }

    // Add residual of previous iteration to RHS
    VectorType LastValues = ZeroVector(LocalSize);
    this->GetFirstDerivativesVector(LastValues);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,LastValues);
}

template< unsigned int TDim >
void StationaryStokes<TDim>::GetDofList(DofsVectorType &rElementalDofList,
                                  ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = (TDim + 1) * NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        rElementalDofList[Index++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[Index++] = rGeom[i].pGetDof(VELOCITY_Y);
        if(TDim > 2) rElementalDofList[Index++] = rGeom[i].pGetDof(VELOCITY_Z);
        rElementalDofList[Index++] = rGeom[i].pGetDof(PRESSURE);
    }
}


template< unsigned int TDim >
void StationaryStokes<TDim>::EquationIdVector(Element::EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = (TDim + 1) * NumNodes;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        if(TDim > 2) rResult[Index++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
        rResult[Index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }
}


template< unsigned int TDim >
void StationaryStokes<TDim>::GetFirstDerivativesVector(Vector &rValues, int Step)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = (TDim + 1) * NumNodes;

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
        if(TDim > 2) rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Z,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(PRESSURE,Step);
    }
}

template< unsigned int TDim >
void StationaryStokes<TDim>::AddMomentumTerms(MatrixType &rLHS,
                                              VectorType &rRHS,
                                              const double Density,
                                              const double Viscosity,
                                              const array_1d<double,3> &BodyForce,
                                              const double TauTwo,
                                              const ShapeFunctionsType &N,
                                              const ShapeDerivativesType &DN_DX,
                                              const double Weight)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = TDim+1;

    unsigned int FirstRow = 0;
    unsigned int FirstCol = 0;
    double Term_ij = 0.0;

    for(unsigned int i = 0; i < NumNodes; ++i)
    {
        // Body force
        for(unsigned int d = 0; d < TDim; ++d)
            rRHS[FirstRow + d] += Weight * N[i] * Density * BodyForce[d];

        for(unsigned int j = 0; j < NumNodes; ++j)
        {
            Term_ij = 0.0;

            // Viscous term
            for(unsigned int d = 0; d < TDim; ++d)
                Term_ij += DN_DX(i,d) * DN_DX(j,d);
            Term_ij *= Viscosity;

            Term_ij *= Density * Weight;

            for(unsigned int d = 0; d < TDim; ++d)
                rLHS(FirstRow+d,FirstCol+d) += Term_ij;

            // Stabilization
            for (unsigned int m = 0; m < TDim; ++m)
                for (unsigned int n = 0; n < TDim; ++n)
                    rLHS(FirstRow+m,FirstCol+n) += TauTwo * Weight * DN_DX(i,m) * DN_DX(j,n);

            // Update column index
            FirstCol += BlockSize;
        }
        // Update matrix indices
        FirstRow += BlockSize;
        FirstCol = 0;
    }
}

template< unsigned int TDim >
void StationaryStokes<TDim>::AddContinuityTerms(MatrixType &rLHS,
                                          VectorType &rRHS,
                                          const double Density,
                                          const array_1d<double,3> &BodyForce,
                                          const double TauOne,
                                          const ShapeFunctionsType &N,
                                          const ShapeDerivativesType &DN_DX,
                                          const double Weight)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = TDim+1;

    unsigned int FirstRow = 0;
    unsigned int FirstCol = 0;

    double DivTerm = 0.0;

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        double Qi = 0.0;
        for (unsigned int d = 0; d < TDim; ++d)
            Qi += DN_DX(i,d) * BodyForce[d];
        rRHS[FirstRow + TDim] += Weight * TauOne * Density * Qi;

        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            double Lij = 0.0;
            for (unsigned int d = 0; d < TDim; ++d)
            {
                Lij += DN_DX(i,d) * DN_DX(j,d);
                DivTerm = Weight * N[i] * DN_DX(j,d);
                rLHS(FirstRow+TDim,FirstCol + d) += DivTerm; // Divergence term
                rLHS(FirstCol + d,FirstRow+TDim) -= DivTerm; // Gradient term
            }

            rLHS(FirstRow+TDim,FirstCol+TDim) += Weight * TauOne * Lij;

            // Update column index
            FirstCol += BlockSize;
        }
        // Update matrix indices
        FirstCol = 0;
        FirstRow += BlockSize;
    }
}

template< unsigned int TDim >
void StationaryStokes<TDim>::CalculateTau(double &TauOne, double &TauTwo, const double DynViscosity)
{
    // Estimate element size
    double h = this->ElementSize();

    TauOne = (h * h) / ( 4.0 * DynViscosity );
    TauTwo = DynViscosity;
}

template< >
double StationaryStokes<2>::ElementSize()
{
    return 1.128379167 * sqrt(this->GetGeometry().DomainSize());
}

template< >
double StationaryStokes<3>::ElementSize()
{
    return 0.60046878 * pow(this->GetGeometry().DomainSize(),0.333333333333333333333);
}

// specialized class declarations

template class StationaryStokes<2>;
template class StationaryStokes<3>;


}
