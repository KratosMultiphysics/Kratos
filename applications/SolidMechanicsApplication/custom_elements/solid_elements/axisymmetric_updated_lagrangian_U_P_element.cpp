//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/axisymmetric_updated_lagrangian_U_P_element.hpp"
#include "solid_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymmetricUpdatedLagrangianUPElement::AxisymmetricUpdatedLagrangianUPElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacementUPElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymmetricUpdatedLagrangianUPElement::AxisymmetricUpdatedLagrangianUPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementUPElement( NewId, pGeometry, pProperties )
{
    //mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
    //mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

AxisymmetricUpdatedLagrangianUPElement::AxisymmetricUpdatedLagrangianUPElement( AxisymmetricUpdatedLagrangianUPElement const& rOther)
    :LargeDisplacementUPElement(rOther)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

AxisymmetricUpdatedLagrangianUPElement&  AxisymmetricUpdatedLagrangianUPElement::operator=(AxisymmetricUpdatedLagrangianUPElement const& rOther)
{
    LargeDisplacementUPElement::operator=(rOther);

    mDeformationGradientF0.clear();
    mDeformationGradientF0.resize(rOther.mDeformationGradientF0.size());

    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
    {
        mDeformationGradientF0[i] = rOther.mDeformationGradientF0[i];
    }

    mDeterminantF0 = rOther.mDeterminantF0;

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer AxisymmetricUpdatedLagrangianUPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new AxisymmetricUpdatedLagrangianUPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer AxisymmetricUpdatedLagrangianUPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    AxisymmetricUpdatedLagrangianUPElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );


    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());
	
	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
      }
    

    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
	NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }

    //-----------//

    if ( NewElement.mDeformationGradientF0.size() != mDeformationGradientF0.size() )
      NewElement.mDeformationGradientF0.resize(mDeformationGradientF0.size());

    for(unsigned int i=0; i<mDeformationGradientF0.size(); i++)
    {
        NewElement.mDeformationGradientF0[i] = mDeformationGradientF0[i];
    }

    NewElement.mDeterminantF0 = mDeterminantF0;

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());
        
    return Element::Pointer( new AxisymmetricUpdatedLagrangianUPElement(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

AxisymmetricUpdatedLagrangianUPElement::~AxisymmetricUpdatedLagrangianUPElement()
{
}

//************************************************************************************
//************************************************************************************


//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

  if (rVariable == DETERMINANT_F){

    const unsigned int& integration_points_number = mConstitutiveLawVector.size();

    
    for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
      {
	mDeterminantF0[PointNumber] = rValues[PointNumber];
	mConstitutiveLawVector[PointNumber]->SetValue(rVariable, rValues[PointNumber], rCurrentProcessInfo);
      }

  }
  else{

    LargeDisplacementElement::SetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

  }


}


//************************************************************************************
//************************************************************************************

//**********************************GET DOUBLE VALUE**********************************
//************************************************************************************


void AxisymmetricUpdatedLagrangianUPElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

  if (rVariable == DETERMINANT_F){

    const unsigned int& integration_points_number = mConstitutiveLawVector.size();

    if ( rValues.size() != integration_points_number )
      rValues.resize( integration_points_number );
    
    for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
      {
	rValues[PointNumber] = mDeterminantF0[PointNumber];
      }

  }
  else{

    LargeDisplacementElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

  }

}

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::Initialize()
{
    KRATOS_TRY

    LargeDisplacementElement::Initialize();

    SizeType integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    //Resize historic deformation gradient
    if ( mDeformationGradientF0.size() != integration_points_number )
        mDeformationGradientF0.resize( integration_points_number );

    if ( mDeterminantF0.size() != integration_points_number )
        mDeterminantF0.resize( integration_points_number, false );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
    {
        mDeterminantF0[PointNumber] = 1;
        mDeformationGradientF0[PointNumber] = identity_matrix<double> (3);
    }


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::InitializeElementVariables (ElementVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    const unsigned int voigt_size      = 4;
    
    rVariables.Initialize(voigt_size,dimension,number_of_nodes);

    rVariables.F.resize(3,3,false);
    rVariables.F = IdentityMatrix(3);
    rVariables.F0.resize(3,3,false);
    rVariables.F0 = IdentityMatrix(3);

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );


    //Calculate Delta Position
    rVariables.DeltaPosition = this->CalculateDeltaPosition(rVariables.DeltaPosition);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

    //stabilization factor
    double StabilizationFactor = 1.0;
    if( GetProperties().Has(STABILIZATION_FACTOR) ){
      StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
    }
    else if( rCurrentProcessInfo.Has(STABILIZATION_FACTOR) ){
      StabilizationFactor = rCurrentProcessInfo[STABILIZATION_FACTOR];
    }
    GetProperties().SetValue(STABILIZATION_FACTOR, StabilizationFactor);
    
}

////************************************************************************************
////************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::FinalizeStepVariables( ElementVariables & rVariables, const double& rPointNumber )
{ 
    //update internal (historical) variables
    mDeterminantF0[rPointNumber]         = rVariables.detF * rVariables.detF0;
    noalias(mDeformationGradientF0[rPointNumber]) = prod(rVariables.F, rVariables.F0);
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementVariables& rVariables, double& rIntegrationWeight)
{

    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;
    if ( this->GetProperties().Has( THICKNESS ) )
      IntegrationWeight /= GetProperties()[THICKNESS];

    //contributions to stiffness matrix calculated on the reference config
    rVariables.detF0   *= rVariables.detF;
    double DeterminantF = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    LargeDisplacementUPElement::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

    rVariables.detF     = DeterminantF;
    rVariables.detF0   /= rVariables.detF;
    //KRATOS_WATCH( rLeftHandSideMatrix )
}


//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{
    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;
    if ( this->GetProperties().Has( THICKNESS ) )
      IntegrationWeight /= GetProperties()[THICKNESS];

    //contribution to external forces
    rVariables.detF0   *= rVariables.detF;
    double DeterminantF = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    LargeDisplacementUPElement::CalculateAndAddRHS( rLocalSystem, rVariables, rVolumeForce, IntegrationWeight );

    rVariables.detF     = DeterminantF;
    rVariables.detF0   /= rVariables.detF;
    //KRATOS_WATCH( rRightHandSideVector )
}


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void AxisymmetricUpdatedLagrangianUPElement::CalculateKinematics(ElementVariables& rVariables,
        const double& rPointNumber)

{
    KRATOS_TRY

    //unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();
    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ(2,2);
    noalias(InvJ) = ZeroMatrix(2,2);
    MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    //std::cout<<" detJ "<<rVariables.detJ<<" Area "<<2*GetGeometry().DomainSize()<<std::endl;  
    
    //Compute cartesian derivatives [dN/dx_n]
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);

    //Calculate IntegrationPoint radius
    this->CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius, rVariables.N);

    //Current Deformation Gradient [dx_n+1/dx_n]
    CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.DeltaPosition, rVariables.CurrentRadius, rVariables.ReferenceRadius);

    //Determinant of the deformation gradient F
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
    Matrix Invj(2,2);
    noalias(Invj) = ZeroMatrix(2,2);
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ); //overwrites detJ

    //Compute cartesian derivatives [dN/dx_n+1]
    noalias(rVariables.DN_DX) = prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

    //Determinant of the Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    rVariables.F0    = mDeformationGradientF0[rPointNumber];

    //Compute the deformation matrix B
    CalculateDeformationMatrix(rVariables.B, rVariables.DN_DX, rVariables.N, rVariables.CurrentRadius);


    KRATOS_CATCH( "" )
}



//*************************COMPUTE AXYSIMMETRIC RADIUS********************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateRadius(double & rCurrentRadius,
						       double & rReferenceRadius,
						       const Vector& rN)


{

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rCurrentRadius=0;
    rReferenceRadius=0;

    if ( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            //Displacement from the reference to the current configuration
            array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;
	    array_1d<double, 3 > & CurrentPosition      = GetGeometry()[i].Coordinates();
	    array_1d<double, 3 > ReferencePosition      = CurrentPosition - DeltaDisplacement;

            rCurrentRadius   += CurrentPosition[0]*rN[i];
            rReferenceRadius += ReferencePosition[0]*rN[i];
	}
	//std::cout<<" CurrentRadius "<<rCurrentRadius<<std::endl;
    }


    if ( dimension == 3 )
    {
        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
    }

    KRATOS_CATCH( "" )
}

//*************************COMPUTE DEFORMATION GRADIENT*******************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateDeformationGradient(const Matrix& rDN_DX,
        Matrix&  rF,
        Matrix&  rDeltaPosition,
        double & rCurrentRadius,
        double & rReferenceRadius)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rF = identity_matrix<double> ( 3 );

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            rF ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
            rF ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
        }

        rF ( 2 , 2 ) = rCurrentRadius/rReferenceRadius;
    }
    else if( dimension == 3)
    {

        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}




//************************************************************************************
//************************************************************************************


void AxisymmetricUpdatedLagrangianUPElement::CalculateDeformationMatrix(Matrix& rB,
        Matrix& rDN_DX,
        Vector& rN,
        double & rCurrentRadius)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rB.clear(); //set all components to zero

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 2 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rN[i]/rCurrentRadius;
            rB( 3, index + 0 ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );

        }

    }
    else if( dimension == 3 )
    {

        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;

    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateGreenLagrangeStrain(const Matrix& rF,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

    //Right Cauchy-Green Calculation
    Matrix C ( 3, 3 );
    noalias( C ) = prod( trans( rF ), rF );

    if( dimension == 2 )
    {

        //Green Lagrange Strain Calculation
        if ( rStrainVector.size() != 4 ) rStrainVector.resize( 4, false );

        rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

        rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

        rStrainVector[2] = 0.5 * ( C( 2, 2 ) - 1.00 );

        rStrainVector[3] = C( 0, 1 ); // xy

    }
    else if( dimension == 3 )
    {

        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;

    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateAlmansiStrain(const Matrix& rF,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Left Cauchy-Green Calculation
    Matrix LeftCauchyGreen(rF.size1(), rF.size1());
    noalias(LeftCauchyGreen) = prod( rF, trans( rF ) );

    //Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen ( 3, 3 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    if( dimension == 2 )
    {

        //Almansi Strain Calculation
        if ( rStrainVector.size() != 4 ) rStrainVector.resize( 4, false );

        rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

        rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

        rStrainVector[2] = 0.5 * ( 1.00 - InverseLeftCauchyGreen( 2, 2 ) );

        rStrainVector[3] = - InverseLeftCauchyGreen( 0, 1 ); // xy

    }
    else if( dimension == 3 )
    {

        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;

    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void AxisymmetricUpdatedLagrangianUPElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
								     ElementVariables & rVariables,
								     double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int indexp = dimension;

    //VectorType Fh=rRightHandSideVector;

    double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

    //double consistent=1;

    double Coefficient = 0;
    Coefficient = this->CalculatePUCoefficient( Coefficient, rVariables ); 

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            // consistent=1;
            // if(i==j)
            //     consistent=2;

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
            rRightHandSideVector[indexp] += (1.0/BulkModulus) * rVariables.N[i] * rVariables.N[j] * Pressure * rIntegrationWeight/ (rVariables.detF0/rVariables.detF);

            //rRightHandSideVector[indexp] += consistent * (1.0/BulkModulus) * (1.0/12.0) * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D

            //std::cout<<" Pressure ["<<j<<"] : "<<Pressure<<" rhs "<<std::endl;

        }


        rRightHandSideVector[indexp] -=  Coefficient * rVariables.N[i] * rIntegrationWeight / (rVariables.detF0/rVariables.detF);

        //std::cout<< " Mpres "<<rRightHandSideVector[indexp]<<" Ppres "<<auxiliar * rVariables.N[i] * rIntegrationWeight / rVariables.detF <<std::endl;

        indexp += (dimension + 1);
    }


    // std::cout<<std::endl;
    // std::cout<<" auxiliar " <<auxiliar<<" F0 "<<rVariables.detF0<<std::endl;
    // std::cout<<" Fpres "<<rRightHandSideVector-Fh<<std::endl;

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
        ElementVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int indexp = dimension;

    // VectorType Fh=rRightHandSideVector;
    // if( this->Id() <= 4)
    //   std::cout<<" Element "<<this->Id()<<" "<<std::endl;

    //use of this variable for the complete parameter: (deffault: 4)
    double AlphaStabilization  = 1.0; 
    double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
    AlphaStabilization *= StabilizationFactor;

    const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
    const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

    double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));

    //Experimental
    // if(LameMu < rVariables.ConstitutiveMatrix(2,2))
    //   LameMu = rVariables.ConstitutiveMatrix(2,2);

    double consistent = 1;

    double FactorValue = 8.0; //JMR deffault value

    unsigned int integration_points = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    if(integration_points == 1)
        AlphaStabilization *= FactorValue/36.0;

    //use of this variable for the complete parameter:
    AlphaStabilization=(AlphaStabilization/LameMu);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
            //2D
            if(integration_points == 1)
            {

                consistent=(-1)*AlphaStabilization;
                if(i==j)
                    consistent=2*AlphaStabilization;

                rRightHandSideVector[indexp] += consistent * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF);
            }
            else
            {
                //AXISYM
                consistent = AlphaStabilization * rIntegrationWeight / (rVariables.detF0/rVariables.detF);
                // if(i==j){
                //   consistent *= ( rVariables.N[i] * rVariables.N[j] + (1.0/9.0) );
                // }
                // else{
                consistent *= ( rVariables.N[i] * rVariables.N[j] - (1.0/3.0) * rVariables.N[i] - (1.0/3.0) * rVariables.N[j] + (1.0/9.0) );
                // }

                rRightHandSideVector[indexp] += consistent * Pressure;
            }

        }


        indexp += (dimension + 1);
    }


    // if( this->Id() <= 4 ){
    //   std::cout<<std::endl;
    //   std::cout<<" FpStab "<<rRightHandSideVector-Fh<<" Stab "<<AlphaStabilization<<" LameMu "<<LameMu<<std::endl;
    // }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateAndAddKuug(MatrixType& rK,
        ElementVariables& rVariables,
        double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    int size = number_of_nodes * dimension;

    Matrix Kuu = zero_matrix<double>(size,size);

    // axisymmetric geometric matrix

    double alpha1 = 0;
    double alpha2 = 0;
    double alpha3 = 0;

    unsigned int indexi = 0;
    unsigned int indexj = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        indexj =0;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            alpha1 = rVariables.DN_DX(j,0) * ( rVariables.DN_DX(i,0) * rVariables.StressVector[0] + rVariables.DN_DX(i,1) * rVariables.StressVector[3] );
            alpha2 = rVariables.DN_DX(j,1) * ( rVariables.DN_DX(i,0) * rVariables.StressVector[3] + rVariables.DN_DX(i,1) * rVariables.StressVector[1] );
            alpha3 = rVariables.N[i] * rVariables.N[j] * rVariables.StressVector[2] * (1.0/rVariables.CurrentRadius*rVariables.CurrentRadius);

            Kuu(indexi,indexj)     = alpha1 + alpha2 + alpha3 ;
            Kuu(indexi+1,indexj+1) = alpha1 + alpha2 ;

            indexj+=2;
        }

        indexi+=2;

    }

    Kuu *= rIntegrationWeight;

    //std::cout<<std::endl;
    //std::cout<<" Kuu "<<Kuu<<std::endl;


    MatrixType Kh=rK;

    //assemble into rk the geometric uu contribution:
    indexi = 0;
    indexj = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
        {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
                {
                    rK(indexi+i,indexj+j)+=Kuu(indexi,indexj);
                    indexj++;
                }
            }
            indexi++;
        }
    }

    //std::cout<<std::endl;
    //std::cout<<" Kgeo "<<rK-Kh<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateAndAddKup (MatrixType& rK,
        ElementVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //MatrixType Kh=rK;
    //contributions to stiffness matrix calculated on the reference configuration
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexp  = dimension;
        unsigned int indexup = dimension * i + i;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rK(indexup+k,indexp) +=  rVariables.DN_DX ( i , k ) *  rVariables.N[j] * rIntegrationWeight * rVariables.detF;

                if(k==0) //axysimmetric term
                    rK(indexup+k,indexp) +=  rVariables.N[i] * rVariables.N[j] * (1.0/rVariables.CurrentRadius) * rIntegrationWeight * rVariables.detF;

            }
            indexp += (dimension + 1);
        }
    }

    // std::cout<<std::endl;
    // std::cout<<" Kup "<<rK-Kh<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateAndAddKpu (MatrixType& rK,
        ElementVariables& rVariables,
        double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //MatrixType Kh=rK;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexp = dimension;

    double DeltaCoefficient = 0;
    DeltaCoefficient = this->CalculatePUDeltaCoefficient( DeltaCoefficient, rVariables ); 

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            int indexup= dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rK(indexp,indexup+k) +=  DeltaCoefficient * rVariables.N[i] * rVariables.DN_DX ( j , k ) * rIntegrationWeight * rVariables.detF;

                //std::cout<<" value ("<<indexp<<","<<indexup+k<<") "<<(2*detF) * rN[i] * rDN_DX ( j , k ) * rIntegrationWeight<<std::endl;
                if(k==0) //axysimmetric term
                    rK(indexp,indexup+k) +=  DeltaCoefficient * rVariables.N[i] * rVariables.N[j] * (1.0/rVariables.CurrentRadius) * rIntegrationWeight * rVariables.detF;

            }
        }
        indexp += (dimension + 1);
    }


    // std::cout<<std::endl;
    // std::cout<<" Kpu "<<rK-Kh<<std::endl;


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateAndAddKpp (MatrixType& rK,
        ElementVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    //repasar

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

    //MatrixType Kh=rK;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexpi = dimension;
    //double consistent = 1.0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexpj = dimension;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            // consistent=1;
            // if(indexpi==indexpj)
            //     consistent=2;

            rK(indexpi,indexpj)  -= ((1.0)/(BulkModulus)) * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight / (rVariables.detF0/rVariables.detF);
            //rK(indexpi,indexpj)  -= consistent * ((1.0)/(BulkModulus)) * (1.0/12.0) * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D

            indexpj += (dimension + 1);
        }

        indexpi += (dimension + 1);
    }

    // std::cout<<std::endl;
    // std::cout<<" Kpp "<<rK-Kh<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateAndAddKppStab (MatrixType& rK,
        ElementVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    //repasar

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // MatrixType Kh=rK;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexpi = dimension;

    //use of this variable for the complete parameter: (deffault: 4)    
    double AlphaStabilization  = 1.0; 
    double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
    AlphaStabilization *= StabilizationFactor;

    const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
    const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

    double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));

    //Experimental
    // if(LameMu < rVariables.ConstitutiveMatrix(2,2))
    //   LameMu = rVariables.ConstitutiveMatrix(2,2);

    double consistent = 1.0;

    double FactorValue = 8.0; //JMR deffault value

    unsigned int integration_points = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    if(integration_points == 1)
        AlphaStabilization *= FactorValue/36.0;

    //use of this variable for the complete parameter:
    AlphaStabilization=(AlphaStabilization/LameMu);
    

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexpj = dimension;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
	    //AXISYM
            if(integration_points == 1)
            {

                consistent=(-1)*AlphaStabilization;
                if(indexpi==indexpj)
                    consistent=2*AlphaStabilization;

                rK(indexpi,indexpj)  -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF);

            }
            else
            {

	        consistent  = AlphaStabilization * rIntegrationWeight / (rVariables.detF0/rVariables.detF);
                // if(indexpi==indexpj){
                //   consistent *= ( rVariables.N[i] * rVariables.N[j] + (1.0/9.0) );
                // }
                // else{
                consistent *= ( rVariables.N[i] * rVariables.N[j] - (1.0/3.0) * rVariables.N[i] - (1.0/3.0) * rVariables.N[j] + (1.0/9.0) );
                // }

                rK(indexpi,indexpj)  -= consistent;;
            }

            indexpj += (dimension + 1);
        }

        indexpi += (dimension + 1);
    }

    // std::cout<<std::endl;
    // std::cout<<" KppStab "<<rK-Kh<<" alha stab "<<GetProperties()[STABILIZATION_FACTOR]<<std::endl;

    KRATOS_CATCH( "" )
}



//************************************CALCULATE TOTAL MASS****************************
//************************************************************************************

double& AxisymmetricUpdatedLagrangianUPElement::CalculateTotalMass( double& rTotalMass, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //Compute the Volume Change acumulated:
    ElementVariables Variables;
    this->InitializeElementVariables(Variables,rCurrentProcessInfo);

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    rTotalMass = 0;
    //reading integration points
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
	//compute element kinematics
	this->CalculateKinematics(Variables,PointNumber);
	
	//getting informations for integration
        double IntegrationWeight = Variables.detJ * integration_points[PointNumber].Weight();

	//compute point volume change
	double PointVolumeChange = 0;
	PointVolumeChange = this->CalculateVolumeChange( PointVolumeChange, Variables );
	
	rTotalMass += PointVolumeChange * GetProperties()[DENSITY] * 2.0 * 3.141592654 * Variables.CurrentRadius * IntegrationWeight;

      }

    return rTotalMass;


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //lumped
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

    if ( rMassMatrix.size1() != MatSize )
        rMassMatrix.resize( MatSize, MatSize, false );

    noalias(rMassMatrix) = ZeroMatrix( MatSize, MatSize );

    // Not Lumped Mass Matrix (numerical integration):

    //reading integration points
    IntegrationMethod CurrentIntegrationMethod = mThisIntegrationMethod; //GeometryData::GI_GAUSS_2; //GeometryData::GI_GAUSS_1;

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( CurrentIntegrationMethod  );

    ElementVariables Variables;
    this->InitializeElementVariables(Variables,rCurrentProcessInfo);

    
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
   
      //compute element kinematics
      this->CalculateKinematics( Variables, PointNumber );

      //getting informations for integration
      double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ * 2.0 * 3.141592654 * Variables.CurrentRadius;


      //compute point volume change
      double PointVolumeChange = 0;
      PointVolumeChange = this->CalculateVolumeChange( PointVolumeChange, Variables );
	
      double CurrentDensity = PointVolumeChange * GetProperties()[DENSITY];

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      	{
      	  unsigned int indexupi = dimension * i + i;

      	  for ( unsigned int j = 0; j < number_of_nodes; j++ )
      	    {
      	      unsigned int indexupj = dimension * j + j;

      	      for ( unsigned int k = 0; k < dimension; k++ )
      		{
      		  rMassMatrix( indexupi+k , indexupj+k ) += Variables.N[i] * Variables.N[j] * CurrentDensity * IntegrationWeight;
      		}
      	    }
      	}

    }

    // Lumped Mass Matrix:

    // double TotalMass = 0;

    // this->CalculateTotalMass( TotalMass, rCurrentProcessInfo );

    // Vector LumpFact(number_of_nodes);    
    // noalias(LumpFact) = ZeroVector(number_of_nodes);
    
    // LumpFact = GetGeometry().LumpingFactors( LumpFact );
    
    // for ( unsigned int i = 0; i < number_of_nodes; i++ )
    // 	{
    // 	  double temp = LumpFact[i] * TotalMass;
    
    // 	  unsigned int indexup = dimension * i + i;
    
    // 	  for ( unsigned int j = 0; j < dimension; j++ )
    // 	    {
    // 	      rMassMatrix( indexup+j , indexup+j ) += temp;
    // 	    }
    // 	}

    // std::cout<<std::endl;
    // std::cout<<" Mass Matrix "<<rMassMatrix<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymmetricUpdatedLagrangianUPElement::GetHistoricalVariables( ElementVariables& rVariables, const double& rPointNumber )
{
    LargeDisplacementElement::GetHistoricalVariables(rVariables,rPointNumber);

    //Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    rVariables.F0    = mDeformationGradientF0[rPointNumber];

    rVariables.CurrentRadius = rVariables.ReferenceRadius;
}

//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& AxisymmetricUpdatedLagrangianUPElement::CalculateVolumeChange( double& rVolumeChange, ElementVariables& rVariables )
{
    KRATOS_TRY
      
    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
  
int AxisymmetricUpdatedLagrangianUPElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = LargeDisplacementUPElement::Check(rCurrentProcessInfo);     

    return ErrorCode;

    KRATOS_CATCH( "" );
}
  
//************************************************************************************
//************************************************************************************
  
void AxisymmetricUpdatedLagrangianUPElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementUPElement )
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
}

void AxisymmetricUpdatedLagrangianUPElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementUPElement )
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
}


} // Namespace Kratos


