//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Albert Puigferrat Perez
//                   Ignasi de Pouplana
//

// Application includes
#include "custom_elements/steady_convection_diffusion_FIC_element.hpp"
#include <math.h>

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SteadyConvectionDiffusionFICElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new SteadyConvectionDiffusionFICElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
int SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& Geom = this->GetGeometry();

    // Base class checks for positive area and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    if(Geom.DomainSize() < 1.0e-15)
        KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;

    // Verify ProcessInfo variables
    if ( CONVECTION_DIFFUSION_SETTINGS.Key() == 0 )
        KRATOS_ERROR << "CONVECTION_DIFFUSION_SETTINGS has Key zero at element " << this->Id() << std::endl;


    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // verify nodal variables and dofs

    const bool IsDefinedUnknownVariable = my_settings->IsDefinedUnknownVariable();
    if (IsDefinedUnknownVariable == false)
    {
        KRATOS_ERROR << "UnknownVariable is not defined" << std::endl;
    }

    const bool IsDefinedVelocityVariable = my_settings->IsDefinedVelocityVariable();
    if (IsDefinedVelocityVariable == false)
    {
        KRATOS_ERROR << "VelocityVariable is not defined" << std::endl;
    }

    const bool IsDefinedMeshVelocityVariable = my_settings->IsDefinedMeshVelocityVariable();
    if (IsDefinedMeshVelocityVariable == false)
    {
        KRATOS_ERROR << "MeshVelocityVariable is not defined" << std::endl;
    }

    const bool IsDefinedVolumeSourceVariable = my_settings->IsDefinedVolumeSourceVariable();
    if (IsDefinedVolumeSourceVariable == false)
    {
        KRATOS_ERROR << "VolumeSourceVariable is not defined" << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Verify properties

    const bool IsDefinedDensityVariable = my_settings->IsDefinedDensityVariable();
    if (IsDefinedDensityVariable == false)
    {
        KRATOS_ERROR << "DensityVariable is not defined" << std::endl;
    }

    const bool IsDefinedSpecificHeatVariable = my_settings->IsDefinedSpecificHeatVariable();
    if (IsDefinedSpecificHeatVariable == false)
    {
        KRATOS_ERROR << "SpecificHeatVariable is not defined" << std::endl;
    }

    const bool IsDefinedDiffusionVariable = my_settings->IsDefinedDiffusionVariable();
    if (IsDefinedDiffusionVariable == false)
    {
        KRATOS_ERROR << "DiffusionVariable is not defined" << std::endl;
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if ( TDim == 2 )
    {
        for (unsigned int i=0; i<TNumNodes; i++)
        {
            if(Geom[i].Z() != 0.0)
                KRATOS_ERROR << "Node with non-zero Z coordinate found. Id: " << Geom[i].Id() << std::endl;
        }
    }

    return ierr;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
GeometryData::IntegrationMethod SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& rGeom = this->GetGeometry();

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (rElementalDofList.size() != TNumNodes)
            rElementalDofList.resize(TNumNodes);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rElementalDofList[i] = rGeom[i].pGetDof(rUnknownVar);
        }


        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes;

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != element_size )
        rLeftHandSideMatrix.resize( element_size, element_size, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( element_size, element_size );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != element_size )
        rRightHandSideVector.resize( element_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( element_size );

    this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_THROW_ERROR(std::logic_error,"SteadyConvectionDiffusionFICElement::CalculateLeftHandSide not implemented","");

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes;

    //Resetting the RHS
    if ( rRightHandSideVector.size() != element_size )
        rRightHandSideVector.resize( element_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( element_size );

    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& rGeom = this->GetGeometry();
    const unsigned int element_size = TNumNodes;

    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    if (rResult.size() != element_size)
        rResult.resize(element_size, false);

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rResult[i] = rGeom[i].GetDof(rUnknownVar).EquationId();
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    GeometryData::IntegrationMethod ThisIntegrationMethod = GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( ThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( ThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Vector detJContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,ThisIntegrationMethod);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,Geom,Prop,CurrentProcessInfo);

    noalias(Variables.VelInter) = ZeroVector(TDim);

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNT
        noalias(Variables.GradNT) = DN_DXContainer[GPoint];

        //Compute N and Interpolated velocity
        noalias(Variables.N) = row(NContainer,GPoint);

        // TODO: This will be moved to an utility
        InterpolateVariableWithComponents(Variables.VelInter,NContainer,Variables.NodalVel,GPoint);
        this->CalculateHVector(Variables,Prop,CurrentProcessInfo);

        Variables.QSource = 0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            Variables.QSource += Variables.N[i]*Variables.NodalQSource[i];
        }

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        //Contributions to the left hand side
        this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }


    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    GeometryData::IntegrationMethod ThisIntegrationMethod = GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( ThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( ThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
    Vector detJContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,ThisIntegrationMethod);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,Geom,Prop,CurrentProcessInfo);

    noalias(Variables.VelInter) = ZeroVector(TDim);

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNT
        noalias(Variables.GradNT) = DN_DXContainer[GPoint];

        //Compute N and Interpolated velocity
        noalias(Variables.N) = row(NContainer,GPoint);

        // TODO: This will be moved to an utility
        InterpolateVariableWithComponents(Variables.VelInter,NContainer,Variables.NodalVel,GPoint);
        this->CalculateHVector(Variables,Prop,CurrentProcessInfo);

        Variables.QSource = 0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            Variables.QSource += Variables.N[i]*Variables.NodalQSource[i];
        }

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }


    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::InitializeElementVariables(ElementVariables& rVariables,
                                                                                  const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

    const GeometryType& Geom = this->GetGeometry();

    //Properties variables
    const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
    double FluidDensity = Prop[rDensityVar];

    const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
    double conductivity = Prop[rDiffusionVar];

    const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();
	double specific_heat = Prop[rSpecificHeatVar];

    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    // Compute rho*c
    rVariables.rho_dot_c = FluidDensity*specific_heat;

    // Compute DifMatrixK
    noalias(rVariables.DifMatrixK) = ZeroMatrix( TDim, TDim );
    for (unsigned int i = 0; i < TDim; i++)
    {
        rVariables.DifMatrixK(i,i) = conductivity;
    }

    //Nodal Variables
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rVariables.NodalPhi[i] = Geom[i].FastGetSolutionStepValue(rUnknownVar);

        rVariables.NodalVel[i]=ZeroVector(3);
        const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
		const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();

        rVariables.NodalVel[i] = Geom[i].FastGetSolutionStepValue(rVelocityVar) - Geom[i].FastGetSolutionStepValue(rMeshVelocityVar);
;

        const Variable<double>& rVolumeSourceVar = my_settings->GetVolumeSourceVariable();
        rVariables.NodalQSource[i] = Geom[i].FastGetSolutionStepValue(rVolumeSourceVar);
    }

    KRATOS_CATCH( "" )
}

// TODO: This will be moved to an utility
template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateHVector(ElementVariables& rVariables, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();

    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

    //Properties variables
    const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
    double conductivity = Prop[rDiffusionVar];

    double NormVel = rVariables.VelInter[0]*rVariables.VelInter[0];
    for (unsigned int d = 1; d < TDim; d++)
        NormVel += rVariables.VelInter[d]*rVariables.VelInter[d];
    NormVel = std::sqrt(NormVel);

    double Domain = rGeom.DomainSize();

    if (TDim == 2)
    {
        rVariables.lv = std::sqrt(2.0*Domain);
    }
    else
    {
        rVariables.lv = pow( (6.0*Domain/Globals::Pi) , (1.0/3.0) );
    }

    double Hu;
    if (NormVel > 1.0e-6)
    {
        // TODO: This will be moved to an utility
        Hu = this->ProjectedElementSize(rGeom,rVariables.VelInter);
        rVariables.Peclet = NormVel * Hu / (2.0 * conductivity );
        rVariables.AlphaV = (1.0/tanh(rVariables.Peclet)-1.0/rVariables.Peclet);

        // Compute HVector
        for (unsigned int i = 0; i < TDim; i++)
        {
            rVariables.HVector[i] = rVariables.AlphaV * rVariables.lv * rVariables.VelInter[i]/NormVel;
        }
    }
    else
    {
        // TODO: This will be moved to an utility
        Hu = this->AverageElementSize(rGeom);

        NormVel = 1.0e-6;

        rVariables.Peclet = NormVel * Hu / (2.0 * conductivity );
        rVariables.AlphaV = (1.0/tanh(rVariables.Peclet)-1.0/rVariables.Peclet);

        // Compute HVector
        for (unsigned int i = 0; i < TDim; i++)
        {
            rVariables.HVector[i] = rVariables.AlphaV * rVariables.lv * rVariables.VelInter[i]/NormVel;
        }
    }

}

//----------------------------------------------------------------------------------------

// Triangle2D3 version. TODO: The rest of geometries will be in a utility
template< unsigned int TDim, unsigned int TNumNodes >
double SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::ProjectedElementSize(const Geometry<Node<3> >& rGeometry,
                                                        const array_1d<double,3>& rVelocity)
{
    double Hu = 0.0;

    // Loop over edges looking for maximum 'projected' length
    array_1d<double,3> Edge(3,0.0);

    for(unsigned int i = 0; i < TNumNodes; ++i)
    {
        double lu = 0.0;
        unsigned int j = (i+1) % TNumNodes;
        Edge = rGeometry[j] - rGeometry[i];
        lu = rVelocity[0] * Edge[0];
        for (unsigned int d = 1; d < TDim; ++d)
            lu += rVelocity[d] * Edge[d];
        lu = fabs(lu);
        if(Hu < lu) Hu = lu;
    }

    if (Hu > 0.0)
    {
        double NormVel = std::sqrt(rVelocity[0]*rVelocity[0] + rVelocity[1]*rVelocity[1] + rVelocity[2]*rVelocity[2]);
        Hu /= NormVel;
    }

    return Hu;
}

//----------------------------------------------------------------------------------------

// TODO: AverageElementSize will be in a utility
// Triangle2D3 version.
template< >
double SteadyConvectionDiffusionFICElement<2,3>::AverageElementSize(const Geometry<Node<3> >& rGeometry)
{

    double x10 = rGeometry[1].X() - rGeometry[0].X();
    double y10 = rGeometry[1].Y() - rGeometry[0].Y();

    double x20 = rGeometry[2].X() - rGeometry[0].X();
    double y20 = rGeometry[2].Y() - rGeometry[0].Y();

    return std::sqrt(0.5 * (x10*y20-x20*y10) );
}

//----------------------------------------------------------------------------------------

// Quadrilateral2D4 version.
template< >
double SteadyConvectionDiffusionFICElement<2,4>::AverageElementSize(const Geometry<Node<3> >& rGeometry)
{

    double x10 = rGeometry[1].X() - rGeometry[0].X();
    double y10 = rGeometry[1].Y() - rGeometry[0].Y();

    double x30 = rGeometry[3].X() - rGeometry[0].X();
    double y30 = rGeometry[3].Y() - rGeometry[0].Y();

    return std::sqrt(x10*y30-x30*y10);
}

//----------------------------------------------------------------------------------------


// Tetrahedra3D4 version.
template<>
double SteadyConvectionDiffusionFICElement<3,4>::AverageElementSize(const Geometry<Node<3> >& rGeometry) {

    double x10 = rGeometry[1].X() - rGeometry[0].X();
    double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    double z10 = rGeometry[1].Z() - rGeometry[0].Z();

    double x20 = rGeometry[2].X() - rGeometry[0].X();
    double y20 = rGeometry[2].Y() - rGeometry[0].Y();
    double z20 = rGeometry[2].Z() - rGeometry[0].Z();

    double x30 = rGeometry[3].X() - rGeometry[0].X();
    double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    double z30 = rGeometry[3].Z() - rGeometry[0].Z();

    double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

    return pow(detJ/6.0,1./3.);
}

//----------------------------------------------------------------------------------------

// Hexahedra3D8 version.
template<>
double SteadyConvectionDiffusionFICElement<3,8>::AverageElementSize(const Geometry<Node<3> >& rGeometry) {

    double x10 = rGeometry[1].X() - rGeometry[0].X();
    double y10 = rGeometry[1].Y() - rGeometry[0].Y();
    double z10 = rGeometry[1].Z() - rGeometry[0].Z();

    double x30 = rGeometry[3].X() - rGeometry[0].X();
    double y30 = rGeometry[3].Y() - rGeometry[0].Y();
    double z30 = rGeometry[3].Z() - rGeometry[0].Z();

    double x40 = rGeometry[4].X() - rGeometry[0].X();
    double y40 = rGeometry[4].Y() - rGeometry[0].Y();
    double z40 = rGeometry[4].Z() - rGeometry[0].Z();

    double detJ = x10 * y30 * z40 - x10 * y40 * z30 + y10 * z30 * x40 - y10 * x30 * z40 + z10 * x30 * y40 - z10 * y30 * x40;
    return pow(detJ,1./3.);
}

//----------------------------------------------------------------------------------------

//TODO: This will be moved to an utility
template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::InterpolateVariableWithComponents(array_1d<double,TDim>& rVector,const Matrix& Ncontainer,
                                        const array_1d<array_1d<double,TDim>, TNumNodes>& VariableWithComponents,const unsigned int& GPoint)
{
    noalias(rVector) = ZeroVector(TDim);

    for(unsigned int i=0; i<TNumNodes; i++)
    {

        rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[i][0];
        rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[i][1];

        if(TDim>2)
        {
            rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[i][2];
        }

    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{

    this->CalculateAndAddAdvectionMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddDiffusiveMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddFICMatrix(rLeftHandSideMatrix,rVariables);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddAdvectionMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{

    noalias(rVariables.AdvMatrixAux) = rVariables.rho_dot_c*outer_prod(rVariables.N,rVariables.VelInter);

    noalias(rLeftHandSideMatrix) += prod(rVariables.AdvMatrixAux,trans(rVariables.GradNT))*rVariables.IntegrationCoefficient;

}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddDiffusiveMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{

    noalias(rVariables.DifMatrixAux) = prod(rVariables.GradNT,rVariables.DifMatrixK);

    noalias(rLeftHandSideMatrix) += prod(rVariables.DifMatrixAux,trans(rVariables.GradNT))*rVariables.IntegrationCoefficient;

}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddFICMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{

    noalias(rVariables.FICMatrixAuxOne) = rVariables.rho_dot_c*outer_prod(rVariables.HVector,rVariables.VelInter);
    noalias(rVariables.FICMatrixAuxTwo) = prod(rVariables.FICMatrixAuxOne,trans(rVariables.GradNT));

    noalias(rLeftHandSideMatrix) += 1.0/2.0*prod(rVariables.GradNT,rVariables.FICMatrixAuxTwo)*rVariables.IntegrationCoefficient;

}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{

    //Calculates -r = -K*phi+f

    //-K*Phi
    this->CalculateAndAddRHSAdvection(rRightHandSideVector, rVariables);

    this->CalculateAndAddRHSDiffusive(rRightHandSideVector, rVariables);

    this->CalculateAndAddRHSFIC(rRightHandSideVector, rVariables);

    //+f
    this->CalculateAndAddSourceForce(rRightHandSideVector, rVariables);

}
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddRHSAdvection(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{

    noalias(rVariables.AdvMatrixAux) = outer_prod(rVariables.N,rVariables.VelInter);
    noalias(rVariables.AdvMatrixAuxTwo) = prod(rVariables.AdvMatrixAux,trans(rVariables.GradNT));

    noalias(rRightHandSideVector) -= prod(rVariables.AdvMatrixAuxTwo*rVariables.IntegrationCoefficient, rVariables.NodalPhi);

}
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddRHSDiffusive(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{

    noalias(rVariables.DifMatrixAux) = prod(rVariables.GradNT,rVariables.DifMatrixK);
    noalias(rVariables.DifMatrixAuxTwo) = prod(rVariables.DifMatrixAux,trans(rVariables.GradNT));

    noalias(rRightHandSideVector) -= prod(rVariables.DifMatrixAuxTwo*rVariables.IntegrationCoefficient, rVariables.NodalPhi);

}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddRHSFIC(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{

    noalias(rVariables.FICMatrixAuxOne) = outer_prod(rVariables.HVector,rVariables.VelInter);
    noalias(rVariables.FICMatrixAuxTwo) = prod(rVariables.FICMatrixAuxOne,trans(rVariables.GradNT));
    noalias(rVariables.FICMatrixAuxThree) = prod(rVariables.GradNT,rVariables.FICMatrixAuxTwo);

    noalias(rRightHandSideVector) -= prod(1.0/2.0*rVariables.FICMatrixAuxThree*rVariables.IntegrationCoefficient, rVariables.NodalPhi);

}
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddSourceForce(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{

    noalias(rRightHandSideVector) += rVariables.N*rVariables.QSource*rVariables.IntegrationCoefficient;

}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    // Considering thickness equals 1 in 2D
    rIntegrationCoefficient = weight * detJ;
}

//----------------------------------------------------------------------------------------

template class SteadyConvectionDiffusionFICElement<2,3>;
template class SteadyConvectionDiffusionFICElement<2,4>;
template class SteadyConvectionDiffusionFICElement<3,4>;
template class SteadyConvectionDiffusionFICElement<3,8>;

} // Namespace Kratos
