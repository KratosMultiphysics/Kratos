//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
// 

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/small_displacement_beam_element.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

#include "solid_mechanics_application_variables.h"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  SmallDisplacementBeamElement::SmallDisplacementBeamElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : BeamElement(NewId, pGeometry)
  {
  }

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  SmallDisplacementBeamElement::SmallDisplacementBeamElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : BeamElement(NewId, pGeometry, pProperties)
  {
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  SmallDisplacementBeamElement::SmallDisplacementBeamElement( SmallDisplacementBeamElement const& rOther)
    : BeamElement(rOther)
  {
  }

  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer SmallDisplacementBeamElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
  {
    return Element::Pointer(new SmallDisplacementBeamElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  SmallDisplacementBeamElement::~SmallDisplacementBeamElement()
  {
  }
 

  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    const unsigned int voigt_size      = dimension * (dimension +1) * 0.5;
    
    rVariables.Initialize(voigt_size,dimension,number_of_nodes);
    
    //Compute Section Properties:
    this->CalculateSectionProperties(rVariables.Section);

    rVariables.Length = GetGeometry().Length();

    if(rVariables.Length == 0.00)
      KRATOS_ERROR << "Zero length found in element #" << this->Id() << std::endl;

    
    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    //rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

    //Calculate Delta Position
    rVariables.DeltaPosition = this->CalculateTotalDeltaPosition(rVariables.DeltaPosition);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

    KRATOS_CATCH( "" )
  }
  
  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateKinematics(ElementVariables& rVariables, const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().size();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    //point number
    rVariables.PointNumber = rPointNumber;
      
    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);
    
    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    Vector Jacobian = ZeroVector(3);
    
    //Calculating the jacobian [dx_0/d£]
    for( unsigned int i=0; i<dimension; i++)
      {
	Jacobian[i] = rVariables.J[rPointNumber](i,0);
      }

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_0]
    rVariables.detJ  =  norm_2(Jacobian);
    rVariables.DN_DX = DN_De[rPointNumber] * 1.0/rVariables.detJ; 

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.N, rVariables.DN_DX);
    
    
    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************
  void SmallDisplacementBeamElement::CalculateDeformationMatrix(Matrix& rB, const Vector& rN, const Matrix& rDN_DX)
  {
    KRATOS_TRY
      
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int voigt_size = dimension * (dimension +1) * 0.5;

    if ( rB.size1() != voigt_size || rB.size2() != dimension*number_of_nodes )
      rB.resize(voigt_size, ( (dimension-1) * 3 ) * number_of_nodes, false );
    
    if( dimension == 2 )
      {

	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	    unsigned int index = 2 * i;

	    rB( 0, index + 0 ) = rDN_DX( i, 0 );
	    rB( 1, index + 1 ) = rDN_DX( i, 0 );
	    rB( 1, index + 2 ) = rN[i];
	    rB( 2, index + 2 ) = rDN_DX( i, 0 );

	  }

      }
    else if( dimension == 3 )
      {
	unsigned int index = 0;
	
	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	    index = (dimension-1) * 3 ) * i;

	    rB( 0, index + 0 ) = rDN_DX( i, 0 );
	    rB( 1, index + 1 ) = rDN_DX( i, 0 );
	    rB( 2, index + 2 ) = rDN_DX( i, 0 );

	    rB( 1, index + 5 ) = -rN[i];
	    rB( 2, index + 4 ) = rN[i];

	    rB( 3, index + 3 ) = rDN_DX( i, 0 );
	    rB( 4, index + 4 ) = rDN_DX( i, 0 );
	    rB( 5, index + 5 ) = rDN_DX( i, 0 );

	  }

      }
    else
      {
	KRATOS_ERROR << " something is wrong with the dimension strain matrix " << std::endl;

      }

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateConstitutiveMatrix(ElementVariables& rVariables)
  {
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 ){

      if( rConstitutiveMatrix.size1() != 3 || rConstitutiveMatrix.size2() != 3)
	rConstitutiveMatrix.resize(3,3, false);

      rConstitutiveMatrix = ZeroMatrix(3,3);

      if(GetProperties().Has(LOCAL_CONSTITUTIVE_MATRIX)){

	//Axis local E2
	rConstitutiveMatrix = GetProperties()[LOCAL_CONSTITUTIVE_MATRIX];

      }
      else{
      
	const double PoissonCoefficient = GetProperties()[POISSON_RATIO];
	const double YoungModulus       = GetProperties()[YOUNG_MODULUS];
	const double ShearModulus       = YoungModulus*0.5/(1.0 + PoissonCoefficient);

	//Axis local E2
	rConstitutiveMatrix( 0, 0 ) = ShearModulus * rVariables.Section.Area;         //local vertical axis
	rConstitutiveMatrix( 1, 1 ) = YoungModulus * rVariables.Section.Area;         //local beam axis
	rConstitutiveMatrix( 2, 2 ) = YoungModulus * rVariables.Section.Inertia_z;    //local horizontal axis
	
      }
      
    }
    else{

      if( rConstitutiveMatrix.size1() != 6 || rConstitutiveMatrix.size2() != 6)
	rConstitutiveMatrix.resize(6,6, false);

      rConstitutiveMatrix = ZeroMatrix(6,6);

      if(GetProperties().Has(LOCAL_CONSTITUTIVE_MATRIX)){

	//Axis local E1
	rConstitutiveMatrix = GetProperties()[LOCAL_CONSTITUTIVE_MATRIX];

      }
      else{
      
	const double PoissonCoefficient = GetProperties()[POISSON_RATIO];
	const double YoungModulus       = GetProperties()[YOUNG_MODULUS];
	const double ShearModulus       = YoungModulus*0.5/(1.0 + PoissonCoefficient);


	//Axis local E1
	rConstitutiveMatrix( 0, 0 ) = YoungModulus * rVariables.Section.Area;          //local beam axis		
	rConstitutiveMatrix( 1, 1 ) = ShearModulus * rVariables.Section.Area;          //local vertial axis
	rConstitutiveMatrix( 2, 2 ) = ShearModulus * rVariables.Section.Area;          //local horizontal axis

	rConstitutiveMatrix( 3, 3 ) = ShearModulus * rVariables.Section.Polar_Inertia; //local torsion
	rConstitutiveMatrix( 4, 4 ) = YoungModulus * rVariables.Section.Inertia_y;     //local vertial axis
	rConstitutiveMatrix( 5, 5 ) = YoungModulus * rVariables.Section.Inertia_z;     //local horizontal axis
      }

    }
    
    //std::cout<<" ConstitutiveMatrix "<<rConstitutiveMatrix<<std::endl;

    KRATOS_CATCH( "" )
  }
  
  //************************************************************************************
  //************************************************************************************


  void SmallDisplacementBeamElement::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
							       ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //size needed
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * ( dimension * 2 );

    //initialize local transformation/rotation matrix
    rLocalSystem.RotationMatrix = ZeroMatrix(MatSize,MatSize);

    if ( rLocalSystem.CalculationFlags.Is(SmallDisplacementBeamElement::COMPUTE_LHS_MATRIX) || rLocalSystem.CalculationFlags.Is(SmallDisplacementBeamElement::COMPUTE_RHS_VECTOR) ) {
    
      //Local to Global Transformation Matrix
      this->CalculateTransformationMatrix(rLocalSystem.RotationMatrix);
    }

    //reading integration points (in fact is the two nodes beam element, only one integration point)
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );


    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    
    //auxiliary terms
    Vector VolumeForce;

    //(in fact is the two nodes beam element, only one integration point)
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {

	Vector N = row( Ncontainer, PointNumber);

	if ( rLocalSystem.CalculationFlags.Is(SmallDisplacementBeamElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
	  {
	    this->CalculateAndAddLHS( rLocalSystem );
	  }

	if ( rLocalSystem.CalculationFlags.Is(SmallDisplacementBeamElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
	  {
	    //contribution to external forces
	    VolumeForce  = this->CalculateVolumeForce( VolumeForce, N );

	    this->CalculateAndAddRHS( rLocalSystem , VolumeForce );
	  }

      }


    KRATOS_CATCH( "" )
      }

  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
							 ElementVariables& rVariables,
							 double& rIntegrationWeight)
  {
    KRATOS_TRY

    MatrixType LocalStiffnessMatrix = rLeftHandSideMatrix;
    this->CalculateLocalStiffnessMatrix(Matrix& LocalStiffnessMatrix)

    std::cout<<" Direct Stiffness Matrix  "<<LocalStiffnessMatrix<<std::endl;
    
    //contributions to stiffness matrix calculated on the reference config
    noalias( rLeftHandSideMatrix ) += prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) ); //to be optimized to remove the temporary

    std::cout<<" Stiffness Matrix Timoshenko "<<rLeftHandSideMatrix<<std::endl;  
    
    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, Vector& VolumeForce)
  {
    KRATOS_TRY

    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 
  
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    const unsigned int MatSize         = rRightHandSideVector.size();

    //Calculate Body Force
    VectorType BodyForceVector = ZeroVector(MatSize);
    this->CalculateGlobalBodyForce(BodyForceVector, VolumeForce);

    //Displacements and Rotations Vector
    Vector LocalVector = ZeroVector(MatSize);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	int index = i * ( dimension * 2 );
	LocalVector[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X );
	LocalVector[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y );
	LocalVector[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z );

	LocalVector[index + 3] = GetGeometry()[i].GetSolutionStepValue( ROTATION_X );
	LocalVector[index + 4] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Y );
	LocalVector[index + 5] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z );
      }

    //std::cout<<" LocalVector "<<LocalVector<<std::endl;

    //Stiffness Matrix
    Matrix GlobalMatrix = ZeroMatrix(MatSize);
    if ( rLocalSystem.CalculationFlags.Is(SmallDisplacementBeamElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
      {
	GlobalMatrix = rLocalSystem.GetLeftHandSideMatrix();
      }
    else
      {
	LocalSystemComponents LocalSystem;
	LocalSystem.SetLeftHandSideMatrix(GlobalMatrix);
	LocalSystem.RotationMatrix = rLocalSystem.RotationMatrix;
	this->CalculateAndAddLHS(LocalSystem);
      }

    //Force Vector
    noalias(rRightHandSideVector) += BodyForceVector;
    noalias(rRightHandSideVector) -= prod(GlobalMatrix, LocalVector);

    //std::cout<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateLocalStiffnessMatrix(Matrix& LocalStiffnessMatrix)
  {
    KRATOS_TRY

    const double PoissonCoefficient = GetProperties()[POISSON_RATIO];
    const double YoungModulus       = GetProperties()[YOUNG_MODULUS];
    const double ShearModulus       = YoungModulus*0.5/(1.0 + PoissonCoefficient);

    const double L   = mLength;
    const double LL  = mLength * mLength;
    const double LLL = mLength * mLength * mLength;

    double const EA  = mSection.Area       * YoungModulus;
    double const EIz = mSection.Inertia_z  * YoungModulus;
    double const EIy = mSection.Inertia_y  * YoungModulus;
    double const JG  = mSection.Polar_Inertia * ShearModulus;

    LocalStiffnessMatrix(0,0)	=   (EA)/(L);
    LocalStiffnessMatrix(6,0)	=  -(EA)/(L);

    LocalStiffnessMatrix(1,1)	=   (12*EIz)/(LLL);
    LocalStiffnessMatrix(5,1)	=   (6*EIz)/(LL);
    LocalStiffnessMatrix(7,1)	=  -(12*EIz)/(LLL);
    LocalStiffnessMatrix(11,1)	=   (6*EIz)/(LL);

    LocalStiffnessMatrix(2,2)	=   (12*EIy)/(LLL);
    LocalStiffnessMatrix(4,2)	=  -(6*EIy)/(LL);
    LocalStiffnessMatrix(8,2)	=  -(12*EIy)/(LLL);
    LocalStiffnessMatrix(10,2)	=  -(6*EIy)/(LL);

    LocalStiffnessMatrix(3,3)	=   (JG)/L;
    LocalStiffnessMatrix(9,3)	=  -(JG)/L;

    LocalStiffnessMatrix(2,4)	=  -(6*EIy)/(LL);
    LocalStiffnessMatrix(4,4)	=   (4*EIy)/L;
    LocalStiffnessMatrix(8,4)	=   (6*EIy)/(LL);
    LocalStiffnessMatrix(10,4)	=   (2*EIy)/L;

    LocalStiffnessMatrix(1,5)	=   (6*EIz)/(LL);
    LocalStiffnessMatrix(5,5)	=   (4*EIz)/L;
    LocalStiffnessMatrix(7,5)	=  -(6*EIz)/(LL);
    LocalStiffnessMatrix(11,5)	=   (2*EIz)/L;

    LocalStiffnessMatrix(0,6)	=  -(EA)/( L);
    LocalStiffnessMatrix(6,6)	=   (EA)/( L);

    LocalStiffnessMatrix(1,7)	=  -(12*EIz)/(LLL);
    LocalStiffnessMatrix(5,7)	=  -(6*EIz)/(LL);
    LocalStiffnessMatrix(7,7)	=   (12*EIz)/(LLL);
    LocalStiffnessMatrix(11,7)	=  -(6*EIz)/(LL);

    LocalStiffnessMatrix(2,8)	=  -(12*EIy)/(LLL);
    LocalStiffnessMatrix(4,8)	=   (6*EIy)/(LL);
    LocalStiffnessMatrix(8,8)	=   (12*EIy)/(LLL);
    LocalStiffnessMatrix(10,8)	=   (6*EIy)/(LL);

    LocalStiffnessMatrix(3,9)	=  -(JG)/L;
    LocalStiffnessMatrix(9,9)	=   (JG)/L;

    LocalStiffnessMatrix(2,10)	=  -(6*EIy)/(LL);
    LocalStiffnessMatrix(4,10)	=   (2*EIy)/L;
    LocalStiffnessMatrix(8,10)	=   (6*EIy)/(LL);
    LocalStiffnessMatrix(10,10) =   (4*EIy)/L;

    LocalStiffnessMatrix(1,11)	=   (6*EIz)/(LL);
    LocalStiffnessMatrix(5,11)	=   (2*EIz)/L;
    LocalStiffnessMatrix(7,11)	=  -(6*EIz)/(LL);
    LocalStiffnessMatrix(11,11) =   (4*EIz)/L;


    KRATOS_CATCH( "" )

  }

  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::MapLocalToGlobal(ElementVariables& rVariables, MatrixType& rMatrix)
  {
    KRATOS_TRY
      
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    // calculate rotation matrix from spatial frame to material frame
    this->CalculateTransformationMatrix(rVariables.CurrentRotationMatrix);

    if( dimension == 2 )
      BeamMathUtilsType::MapLocalToGlobal2D(rVariables.CurrentRotationMatrix, rMatrix);
    else
      BeamMathUtilsType::MapLocalToGlobal3D(rVariables.CurrentRotationMatrix, rMatrix);
    
    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::MapLocalToGlobal(ElementVariables& rVariables, VectorType& rVector)
  {
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // calculate rotation matrix from spatial frame to material frame
    this->CalculateTransformationMatrix(rVariables.CurrentRotationMatrix);
   
    if( dimension == 2 )
      BeamMathUtilsType::MapLocalToGlobal2D(rVariables.CurrentRotationMatrix, rVector);
    else
      BeamMathUtilsType::MapLocalToGlobal3D(rVariables.CurrentRotationMatrix, rVector);
    
    KRATOS_CATCH( "" )
  }

  //*****************************************************************************
  //*****************************************************************************

  void SmallDisplacementBeamElement::CalculateTransformationMatrix(Matrix& rRotationMatrix)
  {
    KRATOS_TRY
      
    this->CalculateLocalAxesMatrix(rRotationMatrix);

    QuaternionType CurrentLocalQuaternion = QuaternionType::FromRotationMatrix( LocalTransformationMatrix );

    CurrentRotationQuaternion = mInitialLocalQuaternion.conjugate() * RotationQuaternion;

    CurrentRotationQuaternion.ToRotationMatrix(rRotationMatrix);
    
    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateLocalBodyForce( Vector& rLocalForceVector, Vector& rVolumeForce )

  {
    KRATOS_TRY

    // External force vector calculation
    // external forces are uniformely distributed
    // must be located somewhere else
      
    Vector Weight = rVolumeForce * mSection.Area;

    array_1d<double, 6 > ReferenceCoordinates;
    ReferenceCoordinates(0)= GetGeometry()[0].X0();
    ReferenceCoordinates(1)= GetGeometry()[0].Y0();
    ReferenceCoordinates(2)= GetGeometry()[0].Z0();
    ReferenceCoordinates(3)= GetGeometry()[1].X0();
    ReferenceCoordinates(4)= GetGeometry()[1].Y0();
    ReferenceCoordinates(5)= GetGeometry()[1].Z0();

    Vector LocalBeamAxis;
    LocalBeamAxis.resize(3,false);

    for (unsigned int i=0; i<3; i++)
      {
        LocalBeamAxis[i] = ReferenceCoordinates[i+3] - ReferenceCoordinates[i];
      }

    array_1d<double, 12 > Cargas_X = ZeroVector(12);
    array_1d<double, 12 > Cargas_Y = ZeroVector(12);
    array_1d<double, 12 > Cargas_Z = ZeroVector(12);

    array_1d<double, 2 > Load;

    Vector NormalPlaneProjection;
    NormalPlaneProjection.resize(3,false);

    double alpha =  0.00;
    double sign  =  1.00;

    double sinus;
    double cosinus;

    //Load X : longitudinal direction
    //***********************************
    if(Weight[0]!=0.00)
      {
        NormalPlaneProjection[0] = 0.00;
        NormalPlaneProjection[1] = LocalBeamAxis[1] ;
        NormalPlaneProjection[2] = LocalBeamAxis[2] ;

        if (LocalBeamAxis[0]<0)
	  {
            sign =-1.00;
	  }
        if( norm_2(NormalPlaneProjection)==0 || norm_2(LocalBeamAxis)==0  )
	  {
            alpha = sign * PI * 0.5;
	  }
        else
	  {
            alpha = inner_prod(NormalPlaneProjection,LocalBeamAxis)/(norm_2(LocalBeamAxis)*norm_2(NormalPlaneProjection));
            alpha = sign*acos(alpha);
	  }

        sinus   = sin(alpha);
        cosinus = cos(alpha);;
        if(fabs(sinus) < 1E-7) sinus = 0.00;
        if(fabs(cosinus) < 1E-7) cosinus = 0.00;

	// self weight only
        Load[0]= Weight[0]*sinus;        // Axial load
        Load[1]= Weight[0]*cosinus;      // Global Vertical load

        Cargas_X[0]=   Load[0]*mLength*0.5;                  // Load in X;
        Cargas_X[1]=   -(Load[1]*mLength)*0.5;               // Load in Y; gravity
        Cargas_X[2]=   0.00;				     // Load in Z
        Cargas_X[3]=   0.00;		                     // Torsional Moment X;
        Cargas_X[4]=   0.00;		                     // Moment Y;
        Cargas_X[5]=  -(Load[1])*mLength*mLength/12.00;      // Moment Z;
        Cargas_X[6]=   Load[0]*mLength*0.5;
        Cargas_X[7]=  -(Load[1])*mLength*0.5;
        Cargas_X[8]=   0.00;
        Cargas_X[9]=   0.00;
        Cargas_X[10]=  0.00;
        Cargas_X[11]=  (Load[1])*mLength*mLength/12.00;


        noalias(rLocalForceVector)  = Cargas_X;

      }


    //Load Z : horizontal section direction 
    //***********************************
    if(Weight[2]!=0.00)
      {
        NormalPlaneProjection[0] = LocalBeamAxis[0] ;
        NormalPlaneProjection[1] = LocalBeamAxis[1] ;
        NormalPlaneProjection[2] = 0.00;

        if (LocalBeamAxis[2]<0)
	  {
            sign = 1.00;
	  }
        if( norm_2(NormalPlaneProjection)==0 || norm_2(LocalBeamAxis)==0  )
	  {
            alpha = sign * PI * 0.5;
	  }
        else
	  {
            alpha = inner_prod(NormalPlaneProjection,LocalBeamAxis)/(norm_2(LocalBeamAxis)*norm_2(NormalPlaneProjection));
            alpha = sign*acos(alpha);
	  }

        sinus = sin(alpha);
        cosinus = cos(alpha);

        if(fabs(sinus) < 1E-7) sinus = 0.00;
        if(fabs(cosinus) < 1E-7) cosinus = 0.00;

	// self weight only
        Load[0]= Weight[2]*sinus;         // Axial load
        Load[1]= Weight[2]*cosinus;       // Global Vertical load

        Cargas_Z[0]=  -(Load[0]*mLength)*0.5;              // Fuerza en X;
        Cargas_Z[1]=   0.00;                               // Fuerza en Y;
        Cargas_Z[2]=  -(Load[1]*mLength)*0.5;              // Fuerza en Z; gravity	  		 
	Cargas_Z[3]=   0.00;			           // Torsional Moment X;
        Cargas_Z[4]=  -(Load[1])*mLength*mLength/12.00;    // Moment Y;
        Cargas_Z[5]=   0.00;                               // Moment Z;
        Cargas_Z[6]=  -(Load[0]*mLength)*0.5;
        Cargas_Z[7]=   0.00;
        Cargas_Z[8]=  -(Load[1])*mLength*0.5;
        Cargas_Z[9]=   0.00;
        Cargas_Z[10]=  (Load[1])*mLength*mLength/12.00;
        Cargas_Z[11]=  0.00;

        noalias(rLocalForceVector)  = Cargas_Z;

      }

    //Load Y : vertical section direction 
    //***********************************
    if(Weight[1]!=0.00)
      {
        NormalPlaneProjection = ZeroVector(3);

        NormalPlaneProjection[0] = LocalBeamAxis[0] ;
        NormalPlaneProjection[1] = 0.00 ;
        NormalPlaneProjection[2] = LocalBeamAxis[2];

        if (LocalBeamAxis[1]<0)
	  {
            sign =-1.00;
	  }
        if( norm_2(NormalPlaneProjection)==0 || norm_2(LocalBeamAxis)==0  )
	  {
            alpha = sign * PI * 0.5;
	  }
        else
	  {
            alpha = inner_prod(NormalPlaneProjection,LocalBeamAxis)/(norm_2(LocalBeamAxis)*norm_2(NormalPlaneProjection));
            alpha = sign*acos(alpha);
	  }

        sinus = sin(alpha);
        cosinus = cos(alpha);

        if(fabs(sinus) < 1E-7) sinus = 0.00;
        if(fabs(cosinus) < 1E-7) cosinus = 0.00;

        // self weight only
        Load[0]= Weight[1]*sinus;         // Axial load
        Load[1]= Weight[1]*cosinus;       // Global Vertical load

        Cargas_Y[0]=   -(Load[0]*mLength)*0.5;              // Fuerza en X;
        Cargas_Y[1]=   -(Load[1]*mLength)*0.5;              // Fuerza en Y; gravity
        Cargas_Y[2]=    0.00;		                    // Fuerza en Z
	Cargas_Y[3]=    0.00;			            // Torsional Moment X;
        Cargas_Y[4]=    0.00;				    // Moment Y;
        Cargas_Y[5]=   -(Load[1])*mLength*mLength/12.00;    // Moment Z;
        Cargas_Y[6]=   -Load[0]*mLength*0.5;
        Cargas_Y[7]=   -(Load[1])*mLength*0.5;
        Cargas_Y[8]=    0.00;
        Cargas_Y[9]=    0.00;
        Cargas_Y[10]=   0.00;
        Cargas_Y[11]=   (Load[1])*mLength*mLength/12.00;

        noalias(rLocalForceVector)  = Cargas_Y;
      }

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateGlobalBodyForce( Vector& rGlobalForceVector, Vector& rVolumeForce )

  {
    KRATOS_TRY

    // External force vector calculation
    // external forces are uniformely distributed
    // must be located somewhere else
       
    Vector Weight = rVolumeForce * mSection.Area;

    array_1d<double, 6 > ReferenceCoordinates;

    Vector LocalBeamAxis;
    LocalBeamAxis.resize(3,false);


    ReferenceCoordinates(0)= GetGeometry()[0].X0();
    ReferenceCoordinates(1)= GetGeometry()[0].Y0();
    ReferenceCoordinates(2)= GetGeometry()[0].Z0();
    ReferenceCoordinates(3)= GetGeometry()[1].X0();
    ReferenceCoordinates(4)= GetGeometry()[1].Y0();
    ReferenceCoordinates(5)= GetGeometry()[1].Z0();


    for (unsigned int i=0; i<3; i++)
      {
        LocalBeamAxis[i] = ReferenceCoordinates[i+3] - ReferenceCoordinates[i];
      }


    rGlobalForceVector[0]=   Weight[0]*mLength*0.5;  // Load in X;
    rGlobalForceVector[1]=   Weight[1]*mLength*0.5;  // Load in Y; gravity
    rGlobalForceVector[2]=   Weight[2]*mLength*0.5;  // Load in Z
    rGlobalForceVector[3]=   0.00;		     // Torsional Moment X;
    rGlobalForceVector[4]=   0.00;		     // Moment Y;
    rGlobalForceVector[5]=   0.00;                   // Moment Z; 
    rGlobalForceVector[6]=   Weight[0]*mLength*0.5;  // Load in X; 
    rGlobalForceVector[7]=   Weight[1]*mLength*0.5;  // Load in Y; gravity 
    rGlobalForceVector[8]=   Weight[2]*mLength*0.5;  // Load in Z; 
    rGlobalForceVector[9]=   0.00;                   // Torsional Moment X;
    rGlobalForceVector[10]=  0.00;                   // Moment Y;
    rGlobalForceVector[11]=  0.00;                   // Moment Z;
    

    KRATOS_CATCH( "" )

  }

  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
  {

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension * 2 );

    if(rMassMatrix.size1() != MatSize)
      rMassMatrix.resize (MatSize, MatSize, false);

    rMassMatrix = ZeroMatrix( MatSize, MatSize );

    double TotalMass = 0;
    TotalMass = ( mSection.Area * mLength ) * GetProperties()[DENSITY];

    Vector LumpFact = ZeroVector(number_of_nodes); 

    LumpFact = GetGeometry().LumpingFactors(LumpFact);

    for( unsigned int i=0; i < number_of_nodes; i++ )
      {
        double temp = LumpFact[i] * TotalMass;

        for( unsigned int j=0; j < dimension; j++ )
	  {
 	    unsigned int index = i * (dimension * 2) + j;

            rMassMatrix(index,index) = temp;
	  }

      }

    KRATOS_CATCH( "" )

  }

  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateOnIntegrationPoints(  const Variable<array_1d<double, 3 > >& rVariable,
								    std::vector< array_1d<double, 3 > >& rOutput, 
								    const ProcessInfo& rCurrentProcessInfo )
  {

    KRATOS_TRY

      const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    const unsigned int dimension                  = GetGeometry().WorkingSpaceDimension();
   
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    
    Vector Stress;
    std::vector<Vector> Load(dimension);
    
    //auxiliary terms
    int factor = 1;
    
    //(in fact is the two nodes beam element, only one integration point)
    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
      {
	
	Vector N = row( Ncontainer, PointNumber);
	
	//contribution to external forces
	Vector VolumeForce;
	VolumeForce = this->CalculateVolumeForce( VolumeForce, N );

	this->CalculateLocalNodalStress(Stress, VolumeForce);

        //Transform to Global Stress
	Matrix Rotation = ZeroMatrix(12,12);
	this->CalculateTransformationMatrix(Rotation);
        Stress = prod(Rotation, Stress);

	//std::cout<<" Stress "<<Stress<<std::endl;
    

	//dangerous:
	for(unsigned int i = 0; i<Stress.size(); i++)
	  {
	    if( std::fabs(Stress[i])< 1E-6) Stress[i] = 0.00;
	  }
	
	
	double x_tolerance     = GetGeometry()[1].X0() - GetGeometry()[0].X0();
	double y_tolerance     = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
	const double tolerance = 1E-6;

	if(fabs(x_tolerance)>tolerance)
	  {
	    if(GetGeometry()[1].X0() > GetGeometry()[0].X0())
	      factor = 1;
	  }
	else if(fabs(y_tolerance)>tolerance)
	  {
	    if(GetGeometry()[1].Y0() > GetGeometry()[0].Y0())
	      factor = 1;
	  }
	else
	  {
	    factor = 1; //-1;
	  }
	
	// for( unsigned int i= 0; i< dimension; i++ )
	//   CalculateDistributedBodyForce(i, Load[i], VolumeForce);

	if( Load[PointNumber].size() != 3 )
	  Load[PointNumber].resize(3, false);

	Load[PointNumber] = ZeroVector(3);
	//std::cout<<" VolumeForce "<<VolumeForce<<std::endl;
      }


    //it is written in 3 integration points instead of 1
    const unsigned int&  write_points_number = GetGeometry().IntegrationPointsNumber( GetIntegrationMethod() );
      
    if ( rOutput.size() != write_points_number )
      rOutput.resize( write_points_number );
  
    //only moment in z axis (global ?)
    if(rVariable==MOMENT)
      {
	//internal moment in not using the given load
   
	/// Punto Inical
	rOutput[0][0] = factor * CalculateInternalMoment(Stress[3], Stress[9], Load[0][1], 0.25);  //Stress[3];
	rOutput[0][1] = factor * CalculateInternalMoment(Stress[4], Stress[10], Load[0][1], 0.25);  //Stress[4];
	rOutput[0][2] = factor * CalculateInternalMoment(Stress[5], Stress[11], Load[0][1], 0.25);  //Stress[5];
        //rOutput[0][2] = factor * CalculateInternalMoment(Stress[5], Stress[1], Load[0][1], mLength * 0.25); 


        rOutput[1][0] = factor * CalculateInternalMoment(Stress[3], Stress[9], Load[0][1],  0.5);
        rOutput[1][1] = factor * CalculateInternalMoment(Stress[4], Stress[10], Load[0][1], 0.5);
        rOutput[1][2] = factor * CalculateInternalMoment(Stress[5], Stress[11], Load[0][1], 0.5);
        //rOutput[1][2] = factor * CalculateInternalMoment(Stress[5], Stress[1], Load[0][1], mLength * 0.5);


        rOutput[2][0] = factor * CalculateInternalMoment(Stress[3], Stress[9], Load[0][1], 0.75);
        rOutput[2][1] = factor * CalculateInternalMoment(Stress[4], Stress[10], Load[0][1], 0.75);
        rOutput[2][2] = factor * CalculateInternalMoment(Stress[5], Stress[11], Load[0][1], 0.75);
        //rOutput[2][2] = factor * CalculateInternalMoment(Stress[5], Stress[1], Load[0][1], 3.00 * mLength * 0.25);

      }

    //only force in x and y axis (global?)
    if(rVariable==FORCE)
      {
   
        rOutput[0][0] = factor * CalculateInternalAxil(Stress[0], Load[0][0], mLength * 0.25);
        rOutput[0][1] = factor * CalculateInternalShear(Stress[1], Load[0][1], mLength * 0.25);
        rOutput[0][2] = factor * CalculateInternalShear(Stress[2], Load[0][2], mLength * 0.25);

        rOutput[1][0] = factor * CalculateInternalAxil(Stress[0], Load[0][0], mLength * 0.5);
        rOutput[1][1] = factor * CalculateInternalShear(Stress[1], Load[0][1], mLength * 0.5);
        rOutput[1][2] = factor * CalculateInternalShear(Stress[2], Load[0][2], mLength * 0.5);

        rOutput[2][0] = factor * CalculateInternalAxil(Stress[0], Load[0][0], mLength * 0.75);
        rOutput[2][1] = factor * CalculateInternalShear(Stress[1], Load[0][1], mLength * 0.75);
        rOutput[2][2] = factor * CalculateInternalShear(Stress[2], Load[0][2], mLength * 0.75);
      }

    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************


  double SmallDisplacementBeamElement::CalculateInternalMoment(const double& M1, const double& M2, const double& ShearLoad, const double& X)
  {
    //return Mo - Vo*X + 0.5 * ShearLoad * X * X;
    return M1*(1-X) - M2*X;
  }

  //************************************************************************************
  //************************************************************************************


  double SmallDisplacementBeamElement::CalculateInternalShear(const double& Q, const double& ShearLoad, const double& X)
  {
    return  -Q + ShearLoad * X;
  }

  //************************************************************************************
  //************************************************************************************


  double SmallDisplacementBeamElement::CalculateInternalAxil(const double& N, const double& AxialLoad, const double& X)
  {
    return  -N + AxialLoad * X;

  }


  //************************************************************************************
  //************************************************************************************


  void SmallDisplacementBeamElement::CalculateLocalNodalStress(Vector& Stress, Vector& rVolumeForce)
  {
    KRATOS_TRY

      array_1d<double, 12 > CurrentDisplacement;
    CurrentDisplacement(0)		=   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X);
    CurrentDisplacement(1)		=   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y);
    CurrentDisplacement(2)		=   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z);
    CurrentDisplacement(3)		=   GetGeometry()[0].GetSolutionStepValue(ROTATION_X);
    CurrentDisplacement(4)		=   GetGeometry()[0].GetSolutionStepValue(ROTATION_Y);
    CurrentDisplacement(5)		=   GetGeometry()[0].GetSolutionStepValue(ROTATION_Z);
    CurrentDisplacement(6)		=   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_X);
    CurrentDisplacement(7)		=   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Y);
    CurrentDisplacement(8)		=   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Z);
    CurrentDisplacement(9)		=   GetGeometry()[1].GetSolutionStepValue(ROTATION_X);
    CurrentDisplacement(10)	        =   GetGeometry()[1].GetSolutionStepValue(ROTATION_Y);
    CurrentDisplacement(11)	        =   GetGeometry()[1].GetSolutionStepValue(ROTATION_Z);

    Matrix Rotation = ZeroMatrix(12,12);
    this->CalculateTransformationMatrix(Rotation);

    Matrix LocalMatrix = ZeroMatrix(12,12);
    this->CalculateLocalStiffnessMatrix(LocalMatrix);

    array_1d<double, 12 > LocalDisplacement;
    noalias(LocalDisplacement) = prod(Matrix(trans(Rotation)), CurrentDisplacement);

    Vector LocalForceVector  = ZeroVector(12);
    this->CalculateLocalBodyForce( LocalForceVector, rVolumeForce );

    if( Stress.size() != 12 )
      Stress.resize(12, false);

    // K·u - fext = fint; where fint is named: Local Stress
    noalias(Stress) = prod(LocalMatrix, LocalDisplacement); 
    Stress -= LocalForceVector; 


    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  /**
   * This function provides the place to perform checks on the completeness of the input.
   * It is designed to be called only once (or anyway, not often) typically at the beginning
   * of the calculations, so to verify that nothing is missing from the input
   * or that no common error is found.
   * @param rCurrentProcessInfo
   */
  int  SmallDisplacementBeamElement::Check(const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    BeamElement::Check(rCurrentProcessInfo);

    return 0;

    KRATOS_CATCH( "" )
  }




} // Namespace Kratos


