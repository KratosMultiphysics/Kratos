//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/axisym_spatial_lagrangian_U_P_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymSpatialLagrangianUPElement::AxisymSpatialLagrangianUPElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacementUPElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymSpatialLagrangianUPElement::AxisymSpatialLagrangianUPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementUPElement( NewId, pGeometry, pProperties )
{
    //mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    //mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
    mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

AxisymSpatialLagrangianUPElement::AxisymSpatialLagrangianUPElement( AxisymSpatialLagrangianUPElement const& rOther)
    :LargeDisplacementUPElement(rOther)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

AxisymSpatialLagrangianUPElement&  AxisymSpatialLagrangianUPElement::operator=(AxisymSpatialLagrangianUPElement const& rOther)
{
    LargeDisplacementUPElement::operator=(rOther);

    mDeformationGradientF0.clear();
    mDeformationGradientF0.resize(rOther.mDeformationGradientF0.size());

    for(unsigned int i=0; i<<mConstitutiveLawVector.size(); i++)
    {
        mDeformationGradientF0[i] = rOther.mDeformationGradientF0[i];
    }

    mDeterminantF0 = rOther.mDeterminantF0;

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer AxisymSpatialLagrangianUPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new AxisymSpatialLagrangianUPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer AxisymSpatialLagrangianUPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    AxisymSpatialLagrangianUPElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );


    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());
	
	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
      }
    

    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
	NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }

    //-----------//

    if ( NewElement.mDeformationGradientF0.size() != mDeformationGradientF0.size() )
      NewElement.mDeformationGradientF0.resize(mDeformationGradientF0.size());

    for(unsigned int i=0; i<<mDeformationGradientF0.size(); i++)
    {
        NewElement.mDeformationGradientF0[i] = mDeformationGradientF0[i];
    }

    NewElement.mDeterminantF0 = mDeterminantF0;

        
    return Element::Pointer( new AxisymSpatialLagrangianUPElement(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

AxisymSpatialLagrangianUPElement::~AxisymSpatialLagrangianUPElement()
{
}

//************************************************************************************
//************************************************************************************


//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

void AxisymSpatialLagrangianUPElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

  if (rVariable == DETERMINANT_F){

    const unsigned int& integration_points_number = mConstitutiveLawVector.size();

    
    for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
      {
	mDeterminantF0[PointNumber] = rValues[PointNumber];
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


void AxisymSpatialLagrangianUPElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
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

void AxisymSpatialLagrangianUPElement::Initialize()
{
    KRATOS_TRY

    LargeDisplacementUPElement::Initialize();

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

void AxisymSpatialLagrangianUPElement::InitializeGeneralVariables (GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();

    rVariables.detF  = 1;

    rVariables.detF0 = 1;

    rVariables.DomainSize = 1;

    rVariables.B.resize( 4 , number_of_nodes * 2 );

    rVariables.F.resize( 3, 3 );

    rVariables.F0.resize( 3, 3 );

    rVariables.ConstitutiveMatrix.resize( 4, 4 );

    rVariables.StrainVector.resize( 4 );

    rVariables.StressVector.resize( 4 );

    rVariables.DN_DX.resize( number_of_nodes, 2 );

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );


    //Calculate Delta Position
    rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );


}

////************************************************************************************
////************************************************************************************

void AxisymSpatialLagrangianUPElement::FinalizeStepVariables( GeneralVariables & rVariables, const double& rPointNumber )
{ 
    //update internal (historical) variables
    mDeterminantF0[rPointNumber]         = rVariables.detF0 ;
    mDeformationGradientF0[rPointNumber] = rVariables.F0;
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void AxisymSpatialLagrangianUPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
{

    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

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

void AxisymSpatialLagrangianUPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{
    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

    //contribution to external forces
    rVariables.detF0   *= rVariables.detF;
    double DeterminantF = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    LargeDisplacementUPElement::CalculateAndAddRHS( rLocalSystem, rVariables, rVolumeForce, IntegrationWeight );

    rVariables.detF     = DeterminantF;
    rVariables.detF0   /= rVariables.detF;
    //KRATOS_WATCH( rRightHandSideVector )
}



//************************************CALCULATE TOTAL MASS****************************
//************************************************************************************

double& AxisymSpatialLagrangianUPElement::CalculateTotalMass( double& rTotalMass )
{
    KRATOS_TRY

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( GeometryData::GI_GAUSS_1 );
    const unsigned int PointNumber = 0;
    Vector N = row(Ncontainer , PointNumber);

    double CurrentRadius = 0;
    double ReferenceRadius = 0;
    CalculateRadius (CurrentRadius, ReferenceRadius, N);

    rTotalMass = GetGeometry().DomainSize() * GetProperties()[DENSITY] * 2.0 * 3.141592654 * CurrentRadius;

    return rTotalMass;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymSpatialLagrangianUPElement::SetGeneralVariables(GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues,
        const int & rPointNumber)
{
    LargeDisplacementElement::SetGeneralVariables(rVariables,rValues,rPointNumber);

    //Set extra options for the contitutive law
    Flags &ConstitutiveLawOptions=rValues.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::FINAL_CONFIGURATION);

}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void AxisymSpatialLagrangianUPElement::CalculateKinematics(GeneralVariables& rVariables,
        const double& rPointNumber)

{
    KRATOS_TRY

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();
    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    //Step domain size
    rVariables.DomainSize = rVariables.detJ;

    //Compute cartesian derivatives [dN/dx_n]
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);

    //Calculate IntegrationPoint radius
    CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius, rVariables.N);

    //Current Deformation Gradient [dx_n+1/dx_n]
    CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.DeltaPosition, rVariables.CurrentRadius, rVariables.ReferenceRadius);

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ); //overwrites detJ

    //Compute cartesian derivatives [dN/dx_n+1]
    rVariables.DN_DX = prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

    //Determinant of the Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    rVariables.F0    = mDeformationGradientF0[rPointNumber];

    //Compute the deformation matrix B
    CalculateDeformationMatrix(rVariables.B, rVariables.DN_DX, rVariables.N, rVariables.CurrentRadius);


    KRATOS_CATCH( "" )
}



//*************************COMPUTE AXYSIMMETRIC RADIUS********************************
//************************************************************************************

void AxisymSpatialLagrangianUPElement::CalculateRadius(double & rCurrentRadius,
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
            array_1d<double, 3 > & ReferencePosition    = GetGeometry()[i].Coordinates();
            array_1d<double, 3 > CurrentPosition        = ReferencePosition + DeltaDisplacement;

            rCurrentRadius   += CurrentPosition[0]*rN[i];
            rReferenceRadius += ReferencePosition[0]*rN[i];
            //std::cout<<" node "<<i<<" -> DeltaDisplacement : "<<DeltaDisplacement<<std::endl;
        }
    }


    if ( dimension == 3 )
    {
        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
    }

    KRATOS_CATCH( "" )
}

//*************************COMPUTE DEFORMATION GRADIENT*******************************
//************************************************************************************

void AxisymSpatialLagrangianUPElement::CalculateDeformationGradient(const Matrix& rDN_DX,
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

        KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}




//************************************************************************************
//************************************************************************************


void AxisymSpatialLagrangianUPElement::CalculateDeformationMatrix(Matrix& rB,
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

        KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymSpatialLagrangianUPElement::CalculateGreenLagrangeStrain(const Matrix& rF,
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

        KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymSpatialLagrangianUPElement::CalculateAlmansiStrain(const Matrix& rF,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Left Cauchy-Green Calculation
    Matrix LeftCauchyGreen = prod( rF, trans( rF ) );

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

        KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void AxisymSpatialLagrangianUPElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int indexp = dimension;

    VectorType Fh=rRightHandSideVector;

    double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

    //double consistent=1;

    //double auxiliar = 0.5*(rVariables.detF0*rVariables.detF0-1)/rVariables.detF0); //(J²-1)
    double auxiliar = (std::log(rVariables.detF0)/rVariables.detF0);  //(ln(J))

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            // consistent=1;
            // if(i==j)
            //     consistent=2;

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
            rRightHandSideVector[indexp] += (1.0/BulkModulus) * rVariables.N[i] * rVariables.N[j] * Pressure * rIntegrationWeight/ rVariables.detF0;

            //rRightHandSideVector[indexp] += consistent * (1.0/BulkModulus) * (1.0/12.0) * Pressure * rIntegrationWeight / rVariables.detF0 ; //2D

            //std::cout<<" Pressure ["<<j<<"] : "<<Pressure<<" rhs "<<std::endl;

        }

        //rRightHandSideVector[indexp] -=  auxiliar * rVariables.N[i] * rIntegrationWeight / rVariables.detF;

        rRightHandSideVector[indexp] -=  auxiliar * rVariables.N[i] * rIntegrationWeight / rVariables.detF0;

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

void AxisymSpatialLagrangianUPElement::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
        GeneralVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int indexp = dimension;

    VectorType Fh=rRightHandSideVector;

    const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
    const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

    double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));

    //Experimental
    // if(LameMu < rVariables.ConstitutiveMatrix(2,2))
    //   LameMu = rVariables.ConstitutiveMatrix(2,2);

    //use of this variable for the complete parameter:
    double AlphaStabilization = 4.0; //GetProperties()[STABILIZATION];

    unsigned int integration_points = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    if(integration_points == 1)
        AlphaStabilization = 1.0/4.5;

    AlphaStabilization=(AlphaStabilization/LameMu);


    double consistent = 1;

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

                rRightHandSideVector[indexp] += consistent * Pressure * rIntegrationWeight / rVariables.detF0;
            }
            else
            {
                //AXISYM
                consistent = AlphaStabilization * rIntegrationWeight / rVariables.detF0;
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


    // std::cout<<std::endl;
    // std::cout<<" FpStab "<<rRightHandSideVector-Fh<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void AxisymSpatialLagrangianUPElement::CalculateAndAddKuug(MatrixType& rK,
        GeneralVariables& rVariables,
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

void AxisymSpatialLagrangianUPElement::CalculateAndAddKup (MatrixType& rK,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    MatrixType Kh=rK;
    //contributions to stiffness matrix calculated on the reference configuration
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexp  = dimension;
        unsigned int indexup = dimension * i + i;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rK(indexup+k,indexp) +=  rVariables.DN_DX ( i , k ) *  rVariables.N[j] * rIntegrationWeight;

                if(k==0) //axysimmetric term
                    rK(indexup+k,indexp) +=  rVariables.N[i] * rVariables.N[j] * (1.0/rVariables.CurrentRadius) * rIntegrationWeight;

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

void AxisymSpatialLagrangianUPElement::CalculateAndAddKpu (MatrixType& rK,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    MatrixType Kh=rK;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexp = dimension;

    //double auxiliar = (rVariables.detF0*rVariables.detF0 + 1)/(rVariables.detF0*rVariables.detF0); //(J²-1)
    double auxiliar = (1.0-std::log(rVariables.detF0))/(rVariables.detF0*rVariables.detF0);   //(ln(J))


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            int indexup= dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rK(indexp,indexup+k) +=  auxiliar  * rVariables.N[i] * rVariables.DN_DX ( j , k ) * rIntegrationWeight;

                //std::cout<<" value ("<<indexp<<","<<indexup+k<<") "<<(2*detF) * rN[i] * rDN_DX ( j , k ) * rIntegrationWeight<<std::endl;
                if(k==0) //axysimmetric term
                    rK(indexp,indexup+k) +=  auxiliar  * rVariables.N[i] * rVariables.N[j] * (1.0/rVariables.CurrentRadius) * rIntegrationWeight;

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

void AxisymSpatialLagrangianUPElement::CalculateAndAddKpp (MatrixType& rK,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    //repasar

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

    MatrixType Kh=rK;

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

            rK(indexpi,indexpj)  -= ((1.0)/(BulkModulus)) * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight / rVariables.detF0;
            //rK(indexpi,indexpj)  -= consistent * ((1.0)/(BulkModulus)) * (1.0/12.0) * rIntegrationWeight / rVariables.detF0; //2D

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

void AxisymSpatialLagrangianUPElement::CalculateAndAddKppStab (MatrixType& rK,
        GeneralVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    //repasar

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    MatrixType Kh=rK;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexpi = dimension;
    double consistent = 1.0;

    const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
    const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

    double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));

    //Experimental
    // if(LameMu < rVariables.ConstitutiveMatrix(2,2))
    //   LameMu = rVariables.ConstitutiveMatrix(2,2);


    //use of this variable for the complete parameter:
    double AlphaStabilization =4.0; //GetProperties()[STABILIZATION];

    unsigned int integration_points = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    if(integration_points == 1)
        AlphaStabilization = 1.0/4.5;

    AlphaStabilization=(AlphaStabilization/LameMu);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexpj = dimension;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            if(integration_points == 1)
            {

                consistent=(-1)*AlphaStabilization;
                if(indexpi==indexpj)
                    consistent=2*AlphaStabilization;

                rK(indexpi,indexpj)  -= consistent * rIntegrationWeight / rVariables.detF0;

            }
            else
            {
                //AXISYM
                consistent  = AlphaStabilization * rIntegrationWeight / rVariables.detF0;
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
    // std::cout<<" KppStab "<<rK-Kh<<std::endl;

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************
void AxisymSpatialLagrangianUPElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementUPElement )
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
}

void AxisymSpatialLagrangianUPElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementUPElement )
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
}


} // Namespace Kratos


