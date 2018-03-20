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
    
    // const PropertiesType& Prop = this->GetProperties();
    // const GeometryType& Geom = this->GetGeometry();

    // Base class checks for positive area and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;
    
    // if(Geom.DomainSize() < 1.0e-15)
    //     KRATOS_THROW_ERROR( std::logic_error, "DomainSize < 1.0e-15 for the element ", this->Id() )    
    
    // // Verify generic variables
    // ierr = UPwElement<TDim,TNumNodes>::Check(rCurrentProcessInfo);
    // if(ierr != 0) return ierr;
    
    // // Verify specific properties
    // if ( PERMEABILITY_XX.Key() == 0 || Prop.Has( PERMEABILITY_XX ) == false || Prop[PERMEABILITY_XX] < 0.0 )
    //     KRATOS_THROW_ERROR( std::invalid_argument,"PERMEABILITY_XX has Key zero, is not defined or has an invalid value at element", this->Id() )
    // if ( PERMEABILITY_YY.Key() == 0 || Prop.Has( PERMEABILITY_YY ) == false || Prop[PERMEABILITY_YY] < 0.0 )
    //     KRATOS_THROW_ERROR( std::invalid_argument,"PERMEABILITY_YY has Key zero, is not defined or has an invalid value at element", this->Id() )
    // if ( PERMEABILITY_XY.Key() == 0 || Prop.Has( PERMEABILITY_XY ) == false || Prop[PERMEABILITY_XY] < 0.0 )
    //     KRATOS_THROW_ERROR( std::invalid_argument,"PERMEABILITY_XY has Key zero, is not defined or has an invalid value at element", this->Id() )
    // if(TDim > 2)
    // {
    //     if ( PERMEABILITY_ZZ.Key() == 0 || Prop.Has( PERMEABILITY_ZZ ) == false || Prop[PERMEABILITY_ZZ] < 0.0 )
    //         KRATOS_THROW_ERROR( std::invalid_argument,"PERMEABILITY_ZZ has Key zero, is not defined or has an invalid value at element", this->Id() )
    //     if ( PERMEABILITY_YZ.Key() == 0 || Prop.Has( PERMEABILITY_YZ ) == false || Prop[PERMEABILITY_YZ] < 0.0 )
    //         KRATOS_THROW_ERROR( std::invalid_argument,"PERMEABILITY_YZ has Key zero, is not defined or has an invalid value at element", this->Id() )
    //     if ( PERMEABILITY_ZX.Key() == 0 || Prop.Has( PERMEABILITY_ZX ) == false || Prop[PERMEABILITY_ZX] < 0.0 )
    //         KRATOS_THROW_ERROR( std::invalid_argument,"PERMEABILITY_ZX has Key zero, is not defined or has an invalid value at element", this->Id() )
    // }
    
    // // Verify the constitutive law
    // if ( CONSTITUTIVE_LAW.Key() == 0 || Prop.Has( CONSTITUTIVE_LAW ) == false )
    //     KRATOS_THROW_ERROR( std::invalid_argument, "CONSTITUTIVE_LAW has Key zero or is not defined at element ", this->Id() )
    // if ( Prop[CONSTITUTIVE_LAW] != NULL )
    // {
    //     // Verify compatibility of the element with the constitutive law
    //     ConstitutiveLaw::Features LawFeatures;
    //     Prop[CONSTITUTIVE_LAW]->GetLawFeatures(LawFeatures);
    //     bool correct_strain_measure = false;
    //     for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    //     {
    //         if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Infinitesimal)
    //             correct_strain_measure = true;
    //     }
    //     if( correct_strain_measure == false )
    //         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the element type", " StrainMeasure_Infinitesimal " );
        
    //     // Check constitutive law
    //     ierr = Prop[CONSTITUTIVE_LAW]->Check( Prop, Geom, rCurrentProcessInfo );
    // }
    // else
    //     KRATOS_THROW_ERROR( std::logic_error, "A constitutive law needs to be specified for the element ", this->Id() )
        
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

// template< unsigned int TDim, unsigned int TNumNodes >
// void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::Initialize()
// {
//     KRATOS_TRY
    
//     const PropertiesType& Prop = this->GetProperties();
//     const GeometryType& Geom = this->GetGeometry();
//     const unsigned int NumGPoints = Geom.IntegrationPointsNumber( ThisIntegrationMethod );

//     KRATOS_CATCH( "" )
// }

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (ElementalDofList.size() != TNumNodes)
            ElementalDofList.resize(TNumNodes);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int element_size = TNumNodes * (TDim + 1);
    
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

// template< unsigned int TDim, unsigned int TNumNodes >
// void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
// {
//     KRATOS_TRY
    
//     const unsigned int element_size = TNumNodes * (TDim + 1);
        
//     //Resetting the RHS
//     if ( rRightHandSideVector.size() != element_size )
//         rRightHandSideVector.resize( element_size, false );
//     noalias( rRightHandSideVector ) = ZeroVector( element_size );
    
//     this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
    
//     KRATOS_CATCH( "" )
// }

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
    }

    KRATOS_CATCH("")
}


//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::GetValuesVector( Vector& rValues, int Step )
{
    //const GeometryType& Geom = this->GetGeometry();
    const unsigned int element_size = TNumNodes * (TDim + 1);
    //unsigned int index = 0;

    if ( rValues.size() != element_size )
        rValues.resize( element_size, false );

    // for ( unsigned int i = 0; i < TNumNodes; i++ )
    // {
    //     rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
    //     rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );
    //     if ( TDim > 2 )
    //         rValues[index++] = Geom[i].FastGetSolutionStepValue( DISPLACEMENT_Z, Step );
    //     rValues[index++] = 0.0;
    // }
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
        InterpolateVariableWithComponents(Variables.VelInter,NContainer,Variables.NodalVel,GPoint);



        //Compute QSource <-- COMPROVAR
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

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::InitializeElementVariables(ElementVariables& rVariables,
                                                                                  const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo)
{   
    KRATOS_TRY
    
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

    //Properties variables    
    const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
    double FluidDensity = Prop[rDensityVar];

    const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
    double conductivity = Prop[rDiffusionVar];

    const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();
	double specific_heat = Prop[rSpecificHeatVar];

    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();


    // Compute alpha
    rVariables.alpha = conductivity / (FluidDensity*specific_heat);

    // Compute AlphaMatrix
    noalias(rVariables.AlphaMatrix) = ZeroMatrix( TDim, TDim );
    for (unsigned int i = 0; i < TDim; i++)
    {
        rVariables.AlphaMatrix(i,i) = rVariables.alpha;
    }

    //Nodal Variables

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rVariables.NodalPhi[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
        
        const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
        rVariables.NodalVel[i] = GetGeometry()[i].FastGetSolutionStepValue(rVelocityVar);

        const Variable<double>& rVolumeSourceVar = my_settings->GetVolumeSourceVariable();
        rVariables.NodalQSource[i] = GetGeometry()[i].FastGetSolutionStepValue(rVolumeSourceVar);
    }
    



    // CANVIAR AIXO
    // ElementUtilities::GetDisplacementsVector(rVariables.DisplacementVector,Geom);
    // ElementUtilities::GetVelocitiesVector(rVariables.VelocityVector,Geom);
    
    //Constitutive Law parameters
    //rVariables.N.resize(TNumNodes,false);
    //rVariables.GradNT.resize(TNumNodes,TDim,false);

    
    KRATOS_CATCH( "" )
}

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateHVector(ElementVariables& rVariables)
{
    GeometryType& rGeom = this->GetGeometry();

    double NormVel = rVariables.VelInter[0]*rVariables.VelInter[0];
    for (unsigned int d = 1; d < TDim; d++)
        NormVel += rVariables.VelInter[d]*rVariables.VelInter[d];
    NormVel = std::sqrt(NormVel);

    double Hu;
    if (NormVel > 1.0e-6)
    {
    //        Hu = ElementSizeCalculator<Dim,NumNodes>::ProjectedElementSize(rGeom,rVariables.VelInter);
        Hu = this->ProjectedElementSize(rGeom,rVariables.VelInter);
    }
    else
    {
    //        Hu = ElementSizeCalculator<Dim,NumNodes>::AverageElementSize(rGeom);
        Hu = this->AverageElementSize(rGeom);

    }

    // Compute HVector
    for (unsigned int i = 0; i < TDim; i++)
    {
        // Si la velocitat es 0 s'ha de calcular d'una altra manera
        rVariables.HVector[i] = Hu * rVariables.VelInter[i]/NormVel;
    }


}

//----------------------------------------------------------------------------------------

// Triangle2D3 version.
template< unsigned int TDim, unsigned int TNumNodes >
double SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::ProjectedElementSize(const Geometry<Node<3> >& rGeometry,
                                                        const array_1d<double,3>& rVelocity)                                                    
{

    double Hu = 0.0;

    //const unsigned int TNumNodes = 3;

    // Loop over edges looking for maximum 'projected' length
    array_1d<double,3> Edge(3,0.0);
    double lu = 0.0;
    for(unsigned int i = 0; i < TNumNodes; ++i)
    {
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

// Triangle2D3 version.
template< unsigned int TDim, unsigned int TNumNodes >
double SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::AverageElementSize(const Geometry<Node<3> >& rGeometry) 
{

    double x10 = rGeometry[1].X() - rGeometry[0].X();
    double y10 = rGeometry[1].Y() - rGeometry[0].Y();

    double x20 = rGeometry[2].X() - rGeometry[0].X();
    double y20 = rGeometry[2].Y() - rGeometry[0].Y();

    return std::sqrt(0.5 * (x10*y20-x20*y10) );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::InterpolateVariableWithComponents(array_1d<double,TDim>& rVector,const Matrix& Ncontainer, 
                                        const array_1d<array_1d<double,TDim>, TNumNodes>& VariableWithComponents,const unsigned int& GPoint)
{        
        
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // DE MOMENT ESTÀ ADAPTAT NOMÉS A 2D, CAL POSAR-HI UN IF PER FER SERVIR 3D
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    noalias(rVector) = ZeroVector(TDim);
       
    for(unsigned int i=0; i<TNumNodes; i++)
    {

        rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[i][0];
        rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[i][1];
    //    rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[i][2];

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

    noalias(rVariables.AdvMatrixAux) = outer_prod(rVariables.N,-1.0*rVariables.VelInter);

    noalias(rLeftHandSideMatrix) += prod(rVariables.AdvMatrixAux,trans(rVariables.GradNT))*rVariables.IntegrationCoefficient;
    
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddDiffusiveMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{

    noalias(rVariables.DifMatrixAux) = prod(rVariables.GradNT,rVariables.AlphaMatrix);

    noalias(rLeftHandSideMatrix) += prod(rVariables.DifMatrixAux,trans(rVariables.GradNT))*rVariables.IntegrationCoefficient;
    
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddFICMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables)
{

    noalias(rVariables.FICMatrixAuxOne) = outer_prod(rVariables.HVector,-1.0*rVariables.VelInter);
    noalias(rVariables.FICMatrixAuxTwo) = prod(rVariables.FICMatrixAuxOne,trans(rVariables.GradNT));

    noalias(rLeftHandSideMatrix) += 1/2.0*prod(rVariables.GradNT,rVariables.FICMatrixAuxTwo)*rVariables.IntegrationCoefficient;

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

    this->CalculateAndAddFICForce(rRightHandSideVector, rVariables);

}
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddRHSAdvection(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    noalias(rVariables.AdvMatrixAux) = outer_prod(rVariables.N,-1.0*rVariables.VelInter);
    noalias(rVariables.AdvMatrixAuxTwo) = prod(rVariables.AdvMatrixAux,trans(rVariables.GradNT));

    noalias(rRightHandSideVector) -= prod(rVariables.AdvMatrixAuxTwo*rVariables.IntegrationCoefficient, rVariables.NodalPhi);

}
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddRHSDiffusive(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{

    noalias(rVariables.DifMatrixAux) = prod(rVariables.GradNT,rVariables.AlphaMatrix);
    noalias(rVariables.DifMatrixAuxTwo) = prod(rVariables.DifMatrixAux,trans(rVariables.GradNT));

    noalias(rRightHandSideVector) -= prod(rVariables.DifMatrixAuxTwo*rVariables.IntegrationCoefficient, rVariables.NodalPhi);

}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddRHSFIC(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{

    noalias(rVariables.FICMatrixAuxOne) = outer_prod(rVariables.HVector,-1.0*rVariables.VelInter);
    noalias(rVariables.FICMatrixAuxTwo) = prod(rVariables.FICMatrixAuxOne,trans(rVariables.GradNT));
    noalias(rVariables.FICMatrixAuxThree) = prod(rVariables.GradNT,rVariables.FICMatrixAuxTwo);

    noalias(rRightHandSideVector) -= prod(1/2.0*rVariables.FICMatrixAuxThree*rVariables.IntegrationCoefficient, rVariables.NodalPhi);

}
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddSourceForce(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    noalias(rRightHandSideVector) += rVariables.N*rVariables.QSource*rVariables.IntegrationCoefficient;

}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateAndAddFICForce(VectorType& rRightHandSideVector, ElementVariables& rVariables)
{
    noalias(rRightHandSideVector) += 1/2.0*prod(rVariables.GradNT,rVariables.HVector*rVariables.QSource)*rVariables.IntegrationCoefficient;
    
}

//----------------------------------------------------------------------------------------

// template< unsigned int TDim, unsigned int TNumNodes >
// void SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
// {    
//     KRATOS_TRY

//     KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateRHS method for a particular element ... illegal operation!!", "" )

//     KRATOS_CATCH( "" )
// }

//----------------------------------------------------------------------------------------

template< >
void SteadyConvectionDiffusionFICElement<2,3>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    rIntegrationCoefficient = weight * detJ *1.0;
}

//----------------------------------------------------------------------------------------

template class SteadyConvectionDiffusionFICElement<2,3>;

} // Namespace Kratos
