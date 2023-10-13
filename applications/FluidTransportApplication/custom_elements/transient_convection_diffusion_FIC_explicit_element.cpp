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
#include "custom_elements/transient_convection_diffusion_FIC_explicit_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer TransientConvectionDiffusionFICExplicitElement<TDim,TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new TransientConvectionDiffusionFICExplicitElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer TransientConvectionDiffusionFICExplicitElement<TDim,TNumNodes>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new TransientConvectionDiffusionFICExplicitElement( NewId, pGeom, pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientConvectionDiffusionFICExplicitElement<TDim,TNumNodes>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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

    // Zero for the explicit scheme
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( element_size, element_size );

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientConvectionDiffusionFICExplicitElement<TDim, TNumNodes>::CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                        VectorType& rRightHandSideVector,
                        const ProcessInfo& rCurrentProcessInfo)
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

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    GeometryData::IntegrationMethod ThisIntegrationMethod = SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( ThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( ThisIntegrationMethod );
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);

    Vector detJContainer(NumGPoints);
    Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,ThisIntegrationMethod);

    //Element variables
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,Geom,Prop,rCurrentProcessInfo);

    Variables.IterationNumber = rCurrentProcessInfo[NL_ITERATION_NUMBER];
    BoundedMatrix<double,TNumNodes,TNumNodes> MMatrixAux = ZeroMatrix( TNumNodes, TNumNodes );

    //Loop over integration points
    for( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute GradNT
        noalias(Variables.GradNT) = DN_DXContainer[GPoint];

        //Compute N and Interpolated velocity
        noalias(Variables.N) = row(NContainer,GPoint);

        Variables.QSource = 0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            Variables.QSource += Variables.N[i]*Variables.NodalQSource[i];
        }

        noalias(Variables.VelInter) = ZeroVector(TDim);
        ElementUtilities::InterpolateVariableWithComponents(Variables.VelInter,NContainer,Variables.NodalVel,GPoint);

        noalias(Variables.DifMatrix) = ZeroMatrix( TDim, TDim );

        this->CalculateDiffusivityVariables(Variables,Prop,rCurrentProcessInfo);
        this->CalculateHVector(Variables,Prop,rCurrentProcessInfo);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, detJContainer[GPoint], integration_points[GPoint].Weight() );

        array_1d <double, TNumNodes> AuxMVector1;
        array_1d <double, TNumNodes> AuxMVector2;

        noalias(AuxMVector1) = Variables.rho_dot_c * Variables.N;
        noalias(AuxMVector2) = Variables.rho_dot_c * 0.5 * prod(Variables.GradNT,Variables.HvVector);

        //// M matrix
        BoundedMatrix<double,TNumNodes,TNumNodes> MMatrixAux1 = ZeroMatrix( TNumNodes, TNumNodes );
        BoundedMatrix<double,TNumNodes,TNumNodes> MMatrixAux2 = ZeroMatrix( TNumNodes, TNumNodes );

        MMatrixAux1 = outer_prod(AuxMVector1,Variables.N) * Variables.IntegrationCoefficient;
        MMatrixAux2 = outer_prod(AuxMVector2,Variables.N) * Variables.IntegrationCoefficient;

        // We are not considering MMatrixAux2, which is the h term
        MMatrixAux += MMatrixAux1 ;
    }

    for (unsigned int i = 0 ; i < TNumNodes ; i++ )
    {
        for (unsigned int j = 0 ; j < TNumNodes ; j ++ )
        {
            // LHS = Md lumped
            rLeftHandSideMatrix (i,i) += MMatrixAux(i,j);
        }
        // rLeftHandSideMatrix (i,i) = Geom.Area() / 3.0;

    }

    //rLeftHandSideMatrix = MMatrixAux;

    //RightHandSideVector

    const double& Theta = rCurrentProcessInfo[THETA];
    const double& DeltaTime = rCurrentProcessInfo[DELTA_TIME];

    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    //array_1d<double,TNumNodes> NodalPhiCurrent;
    //array_1d<double,TNumNodes> PrevNodalPhi;

   // const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

  //  for (unsigned int i = 0; i < TNumNodes; i++)
  //  {
  //      NodalPhiCurrent[i] = Geom[i].FastGetSolutionStepValue(rUnknownVar,0);
   //     PrevNodalPhi[i] = Geom[i].FastGetSolutionStepValue(rUnknownVar,2);
   // }

    array_1d<double,TNumNodes> aux_vector;
    //array_1d<double,TNumNodes> aux_vector_2;

    // noalias(aux_vector) = NodalPhiCurrent - 2.0 * Variables.NodalPhi + PrevNodalPhi;
    //noalias(aux_vector_2) = Variables.NodalPhi - PrevNodalPhi;

    //noalias(rRightHandSideVector) -= 1.0 / (Theta * DeltaTime) * (prod(rLeftHandSideMatrix, aux_vector) + prod(MMatrixAux, aux_vector_2));
    noalias(rRightHandSideVector) += 1.0 / (Theta * DeltaTime) * (prod(rLeftHandSideMatrix, Variables.NodalPhi));
    aux_vector = 1.0 / (Theta * DeltaTime) * (prod(rLeftHandSideMatrix, Variables.NodalPhi));

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientConvectionDiffusionFICExplicitElement<TDim,TNumNodes>::InitializeElementVariables(ElementVariables& rVariables,
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
    rVariables.rho_dot_c = FluidDensity * specific_heat;

    // Tolerance
    rVariables.LowTolerance = 1e-8;
    rVariables.HighTolerance = 1e-4;

    // Compute DifMatrixK
    noalias(rVariables.DifMatrixK) = ZeroMatrix( TDim, TDim );
    for (unsigned int i = 0; i < TDim; i++)
    {
        rVariables.DifMatrixK(i,i) = conductivity;
    }

    //Nodal Variables
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rVariables.NodalPhi[i] = Geom[i].FastGetSolutionStepValue(rUnknownVar, 1);

        for (unsigned int j = 0; j < TDim; j++)
        {
            rVariables.NodalPhiGradient [j][i] =Geom[i].FastGetSolutionStepValue(NODAL_PHI_GRADIENT)[j];
        }

        rVariables.NodalVel[i] = ZeroVector(3);
        const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
		const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();

        rVariables.NodalVel[i] = Geom[i].FastGetSolutionStepValue(rVelocityVar) - Geom[i].FastGetSolutionStepValue(rMeshVelocityVar);

        const Variable<double>& rVolumeSourceVar = my_settings->GetVolumeSourceVariable();
        rVariables.NodalQSource[i] = Geom[i].FastGetSolutionStepValue(rVolumeSourceVar);

    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void TransientConvectionDiffusionFICExplicitElement<TDim,TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{
}




//----------------------------------------------------------------------------------------

template class TransientConvectionDiffusionFICExplicitElement<2,3>;
template class TransientConvectionDiffusionFICExplicitElement<2,4>;
template class TransientConvectionDiffusionFICExplicitElement<3,4>;
template class TransientConvectionDiffusionFICExplicitElement<3,8>;

} // Namespace Kratos