//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:         Nelson Maireni $
//   Last modified by:    $Co-Author:         JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/small_displacement_beam_element_3D2N.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

#include "solid_mechanics_application_variables.h"


namespace Kratos
{

/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( SmallDisplacementBeamElement3D2N, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( SmallDisplacementBeamElement3D2N, COMPUTE_LHS_MATRIX,                 1 );


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementBeamElement3D2N::SmallDisplacementBeamElement3D2N(IndexType NewId,GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************


SmallDisplacementBeamElement3D2N::SmallDisplacementBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    KRATOS_TRY

    mThisIntegrationMethod = GeometryData::GI_GAUSS_1;

    //DO NOT ADD DOFS HERE!!!

    KRATOS_CATCH( "" )

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

SmallDisplacementBeamElement3D2N::SmallDisplacementBeamElement3D2N( SmallDisplacementBeamElement3D2N const& rOther)
    :Element(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
{
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer SmallDisplacementBeamElement3D2N::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
  return Kratos::make_shared<SmallDisplacementBeamElement3D2N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementBeamElement3D2N::~SmallDisplacementBeamElement3D2N()
{
}


//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

SmallDisplacementBeamElement3D2N::IntegrationMethod  SmallDisplacementBeamElement3D2N::GetIntegrationMethod() const
{
  //return mThisIntegrationMethod;
  return GeometryData::GI_GAUSS_3; //writing results with 3 integration points....
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo&
                             CurrentProcessInfo)
{
    ElementalDofList.resize(0);

    for ( SizeType i = 0; i < GetGeometry().size(); i++ )
      {
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));

	ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
      }
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::EquationIdVector(EquationIdVectorType& rResult,
                                   ProcessInfo& CurrentProcessInfo)
{

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size    = number_of_nodes * ( dimension * 2 );

    if ( rResult.size() != dofs_size )
        rResult.resize( dofs_size, false );

    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      int index = i * ( dimension * 2 );
      rResult[index]   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
      rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
      rResult[index+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();

      rResult[index+3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
      rResult[index+4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
      rResult[index+5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
    }

}


//*********************************DISPLACEMENT***************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::GetValuesVector(Vector& rValues, int Step)
{

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = number_of_nodes * ( dimension * 2 );

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      int index = i * ( dimension * 2 );
      rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
      rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
      rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );

      rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ROTATION_X, Step );
      rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Y, Step );
      rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z, Step );
    }

}

//************************************VELOCITY****************************************
//************************************************************************************

//************************************************************************************
//************************************************************************************
void SmallDisplacementBeamElement3D2N::GetFirstDerivativesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = number_of_nodes * ( dimension * 2 );

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      int index = i * ( dimension * 2 );
      rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
      rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
      rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );

      rValues[index + 3] = 0.0;
      rValues[index + 4] = 0.0;
      rValues[index + 5] = 0.0;
    }


    KRATOS_CATCH( "" )
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::GetSecondDerivativesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = number_of_nodes * ( dimension * 2 );

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      int index = i * ( dimension * 2 );
      rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
      rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
      rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );

      rValues[index + 3] = 0.0;
      rValues[index + 4] = 0.0;
      rValues[index + 5] = 0.0;
    }


    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

//**********************************GET VECTOR VALUE**********************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::GetValueOnIntegrationPoints( const Variable<array_1d<double, 3 > >& rVariable,
								    std::vector< array_1d<double, 3 > >& rValues,
								    const ProcessInfo& rCurrentProcessInfo )
{

    this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);

}


//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::Initialize()
{
    KRATOS_TRY

    //Compute Setction Properties:
    this->CalculateSectionProperties();

    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
					   VectorType& rRightHandSideVector,
					   Flags& rCalculationFlags)

{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * ( dimension * 2 );

    if ( rCalculationFlags.Is(SmallDisplacementBeamElement3D2N::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }



    //resizing as needed the RHS
    if ( rCalculationFlags.Is(SmallDisplacementBeamElement3D2N::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
	    rRightHandSideVector.resize( MatSize, false );

	rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS

    }

    KRATOS_CATCH( "" )
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::CalculateSectionProperties()

{
    KRATOS_TRY

    if( GetProperties().Has(CROSS_SECTION_AREA) ){
        mSection.Area = GetProperties()[CROSS_SECTION_AREA];
    }
    else{
        mSection.Area = GetValue(CROSS_SECTION_AREA);
    }


    if( GetProperties().Has(LOCAL_INERTIA_TENSOR) )
    {
        Matrix& inertia = GetProperties()[LOCAL_INERTIA_TENSOR];
        mSection.Inertia_z = inertia(0,0);
        mSection.Inertia_y = inertia(1,1);
        mSection.Polar_Inertia = inertia(0,1);
    }
    else
    {
        Matrix& inertia = GetValue(LOCAL_INERTIA_TENSOR);
        mSection.Inertia_z = inertia(0,0);
        mSection.Inertia_y = inertia(1,1);
        mSection.Polar_Inertia = inertia(0,1);
    }

    mLength = GetGeometry().Length();

    if(mLength == 0.00)
        KRATOS_THROW_ERROR(std::invalid_argument, "Zero length found in elemnet #", this->Id());


    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************


void SmallDisplacementBeamElement3D2N::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
								 ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //size needed
    const SizeType number_of_nodes  = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * ( dimension * 2 );

    //initialize local transformation/rotation matrix
    rLocalSystem.RotationMatrix     = ZeroMatrix(MatSize,MatSize);

    if ( rLocalSystem.CalculationFlags.Is(SmallDisplacementBeamElement3D2N::COMPUTE_LHS_MATRIX) || rLocalSystem.CalculationFlags.Is(SmallDisplacementBeamElement3D2N::COMPUTE_RHS_VECTOR) ) {

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

	Vector N = matrix_row<const Matrix>( Ncontainer, PointNumber);

	if ( rLocalSystem.CalculationFlags.Is(SmallDisplacementBeamElement3D2N::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
	  {
	    this->CalculateAndAddLHS( rLocalSystem );
	  }

	if ( rLocalSystem.CalculationFlags.Is(SmallDisplacementBeamElement3D2N::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
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

void SmallDisplacementBeamElement3D2N::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem)
{
    KRATOS_TRY

    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

    unsigned int MatSize = rLeftHandSideMatrix.size1();

    //Initialize Local Matrices
    Matrix LocalStiffnessMatrix = ZeroMatrix(MatSize,MatSize);

    //Local Stiffness Matrix
    CalculateLocalStiffnessMatrix(LocalStiffnessMatrix);

    Matrix aux_matrix   = ZeroMatrix(MatSize,MatSize);
    noalias(aux_matrix) = prod(rLocalSystem.RotationMatrix, LocalStiffnessMatrix);


    //Stiffness Matrix
    noalias(rLeftHandSideMatrix) = prod(aux_matrix,Matrix(trans(rLocalSystem.RotationMatrix)));

    // std::cout<<" rLocalSystem.RotationMatrix "<<rLocalSystem.RotationMatrix<<std::endl;
    // std::cout<<" rLeftHandSideMatrix "<<rLeftHandSideMatrix<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, Vector& VolumeForce)
{
    KRATOS_TRY

    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    const unsigned int MatSize         = rRightHandSideVector.size();

    //Calculate Body Force
    VectorType BodyForceVector = ZeroVector(MatSize);
    this->CalculateGlobalBodyForce(BodyForceVector, VolumeForce);

    //Displacements and Rotations Vector
    Vector LocalVector = ZeroVector(MatSize);
    for ( SizeType i = 0; i < number_of_nodes; i++ )
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
    Matrix GlobalMatrix(MatSize,MatSize);
    noalias(GlobalMatrix) = ZeroMatrix(MatSize,MatSize);
    if ( rLocalSystem.CalculationFlags.Is(SmallDisplacementBeamElement3D2N::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
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

void  SmallDisplacementBeamElement3D2N::CalculateRightHandSide(VectorType& rRightHandSideVector,
					  ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(SmallDisplacementBeamElement3D2N::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(SmallDisplacementBeamElement3D2N::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector,  LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(RightHandSideVector);

    //Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(SmallDisplacementBeamElement3D2N::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(SmallDisplacementBeamElement3D2N::COMPUTE_LHS_MATRIX);

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::CalculateLocalStiffnessMatrix(Matrix& LocalStiffnessMatrix)
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


//*****************************************************************************
//*****************************************************************************


void SmallDisplacementBeamElement3D2N::CalculateTransformationMatrix(Matrix& rRotationMatrix)

{

    KRATOS_TRY
    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int       size            = number_of_nodes * dimension;
    unsigned int       MatSize         = 2 * size;


    Vector GlobalY = ZeroVector(3);
    GlobalY[1]=1.0;

    Vector GlobalZ = ZeroVector(3);
    GlobalZ[2]=1.0;

    Vector DirectionVectorX     = ZeroVector(3);
    Vector ReferenceCoordinates = ZeroVector(size);


    ReferenceCoordinates[0] = GetGeometry()[0].X();
    ReferenceCoordinates[1] = GetGeometry()[0].Y();
    ReferenceCoordinates[2] = GetGeometry()[0].Z();

    int k = number_of_nodes - 1 ;

    ReferenceCoordinates[3] = GetGeometry()[k].X();
    ReferenceCoordinates[4] = GetGeometry()[k].Y();
    ReferenceCoordinates[5] = GetGeometry()[k].Z();

    for( SizeType i = 0; i < dimension; i++ )
    {
      DirectionVectorX[i]  = (ReferenceCoordinates[i+3] - ReferenceCoordinates[i]);
    }

    // local x-axis (e1_local) is the beam axis  (in GID is e3_local)
    double VectorNorm = MathUtils<double>::Norm(DirectionVectorX);
    if( VectorNorm != 0)
      DirectionVectorX /= VectorNorm;

    // local y-axis (e2_local) (in GID is e1_local)
    Vector DirectionVectorY = ZeroVector(3);

    double tolerance = 1.0/64.0;
    if(fabs(DirectionVectorX[0])< tolerance && fabs(DirectionVectorX[1])< tolerance){
      MathUtils<double>::CrossProduct(DirectionVectorY, GlobalY, DirectionVectorX);
    }
    else{
      MathUtils<double>::CrossProduct(DirectionVectorY, GlobalZ, DirectionVectorX);
    }

    VectorNorm = MathUtils<double>::Norm(DirectionVectorY);
    if( VectorNorm != 0)
      DirectionVectorY /= VectorNorm;

    // local z-axis (e3_local) (in GID is e2_local)
    Vector DirectionVectorZ;
    MathUtils<double>::CrossProduct(DirectionVectorZ, DirectionVectorX,DirectionVectorY);

    VectorNorm = MathUtils<double>::Norm(DirectionVectorZ);
    if( VectorNorm != 0 )
      DirectionVectorZ /= VectorNorm;


    //Transformation matrix T = [e1_local, e2_local, e3_local]
    Matrix AuxRotationMatrix;

    if( AuxRotationMatrix.size1() != dimension )
      AuxRotationMatrix.resize(dimension, dimension, false);

    // std::cout<<" Xlocal "<<DirectionVectorX<<std::endl;
    // std::cout<<" Ylocal "<<DirectionVectorY<<std::endl;
    // std::cout<<" Zlocal "<<DirectionVectorZ<<std::endl;

    for (SizeType i=0; i<dimension; i++)
      {
	AuxRotationMatrix(i,0) = DirectionVectorX[i];  // column distribution
	AuxRotationMatrix(i,1) = DirectionVectorY[i];
	AuxRotationMatrix(i,2) = DirectionVectorZ[i];
      }


    if( rRotationMatrix.size1() != MatSize )
      rRotationMatrix.resize(MatSize, MatSize, false);

    rRotationMatrix = ZeroMatrix(MatSize,MatSize);

    //Building the rotation matrix for the local element matrix
    for (unsigned int kk=0; kk < MatSize; kk += dimension)
    {
        for (SizeType i=0; i<dimension; i++)
        {
            for(SizeType j=0; j<dimension; j++)
            {
	      rRotationMatrix(i+kk,j+kk) = AuxRotationMatrix(i,j);
            }
        }
    }


    KRATOS_CATCH( "" )

}

//************************************CALCULATE VOLUME ACCELERATION*******************
//************************************************************************************

Vector&  SmallDisplacementBeamElement3D2N::CalculateVolumeForce( Vector& rVolumeForce, const Vector &rN)
{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    rVolumeForce = ZeroVector(dimension);
    for ( SizeType j = 0; j < number_of_nodes; j++ )
    {
        if( GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION) ) //temporary, will be checked once at the beginning only
            rVolumeForce += rN[j] * GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    rVolumeForce *= GetProperties()[DENSITY];

    return rVolumeForce;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::CalculateLocalBodyForce( Vector& rLocalForceVector, Vector& rVolumeForce )

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
            alpha = sign * 0.5 * Globals::Pi;
        }
        else
        {
            alpha = inner_prod(NormalPlaneProjection,LocalBeamAxis)/(norm_2(LocalBeamAxis)*norm_2(NormalPlaneProjection));
            alpha = sign*acos(alpha);
        }

        sinus   = sin(alpha);
        cosinus = cos(alpha);
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
            alpha = sign * 0.5 * Globals::Pi;
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
            alpha = sign * 0.5 * Globals::Pi;
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

void SmallDisplacementBeamElement3D2N::CalculateGlobalBodyForce( Vector& rGlobalForceVector, Vector& rVolumeForce )

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

void SmallDisplacementBeamElement3D2N::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension * 2 );

    if(rMassMatrix.size1() != MatSize)
        rMassMatrix.resize (MatSize, MatSize, false);

    rMassMatrix = ZeroMatrix( MatSize, MatSize );

    double TotalMass = 0;
    TotalMass = ( mSection.Area * mLength ) * GetProperties()[DENSITY];

    Vector LumpFact = ZeroVector(number_of_nodes);

    LumpFact = GetGeometry().LumpingFactors(LumpFact);

    for( SizeType i=0; i < number_of_nodes; i++ )
    {
        double temp = LumpFact[i] * TotalMass;

        for( SizeType j=0; j < dimension; j++ )
        {
 	    unsigned int index = i * (dimension * 2) + j;

            rMassMatrix(index,index) = temp;
        }

    }

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************


void SmallDisplacementBeamElement3D2N::CalculateDistributedBodyForce(const int Direction, Vector& Load, Vector& rVolumeForce)
{
    KRATOS_TRY

    Vector Weight = mSection.Area * rVolumeForce ;

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

    if( Load.size() != 3 )
      Load.resize(3, false);

    double alpha  =  0.00;
    double sign  =  1.00;

    double sinus;
    double cosinus;

    Vector NormalPlaneProjection;
    NormalPlaneProjection.resize(3,false);

    if(Direction==0)  // 0->x, 1->y, 2->z
    {
	NormalPlaneProjection[0] = 0.00;
	NormalPlaneProjection[1] = LocalBeamAxis[1];
	NormalPlaneProjection[2] = LocalBeamAxis[2];

	if (LocalBeamAxis[0]<0)
	  {
	    sign =-1.00;
	  }
	if( norm_2(NormalPlaneProjection)==0 || norm_2(LocalBeamAxis)==0  )
	  {
	    alpha = sign * 0.5 * Globals::Pi;
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
	Load[0]= Weight[0]*sinus;         // Axial load
	Load[1]= Weight[0]*cosinus;       // Global Vertical load
    }

    if(Direction==1) // 0->x, 1->y, 2->z
    {
        NormalPlaneProjection = ZeroVector(3);
        NormalPlaneProjection[0] = LocalBeamAxis[0];
        NormalPlaneProjection[1] = 0.00 ;
        NormalPlaneProjection[2] = LocalBeamAxis[2];

        if (LocalBeamAxis[1]<0)
        {
            sign =-1.00;
        }
        if( norm_2(NormalPlaneProjection)==0 || norm_2(LocalBeamAxis)==0  )
        {
            alpha = sign * 0.5 * Globals::Pi;
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
    }

    if(Direction==2) // 0->x, 1->y, 2->z
    {
        NormalPlaneProjection[0] = LocalBeamAxis[0] ;
        NormalPlaneProjection[1] = LocalBeamAxis[1] ;
        NormalPlaneProjection[2] = 0.00;

        if (LocalBeamAxis[2]<0)
        {
            sign =-1.00;
        }
        if( norm_2(NormalPlaneProjection)==0 || norm_2(LocalBeamAxis)==0  )
        {
            alpha = sign * 0.5 * Globals::Pi;
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
    }

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void SmallDisplacementBeamElement3D2N::CalculateOnIntegrationPoints(  const Variable<array_1d<double, 3 > >& rVariable,
								      std::vector< array_1d<double, 3 > >& rOutput,
								      const ProcessInfo& rCurrentProcessInfo )
{

   KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    const SizeType dimension                      = GetGeometry().WorkingSpaceDimension();

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

    Vector Stress;
    std::vector<Vector> Load(dimension);

    //auxiliary terms
    int factor = 1;

    //(in fact is the two nodes beam element, only one integration point)
    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
      {

	Vector N = matrix_row<const Matrix>( Ncontainer, PointNumber);

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

	// for( SizeType i= 0; i< dimension; i++ )
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


double SmallDisplacementBeamElement3D2N::CalculateInternalMoment(const double& M1, const double& M2, const double& ShearLoad, const double& X)
{
    //return Mo - Vo*X + 0.5 * ShearLoad * X * X;
    return M1*(1-X) - M2*X;
}

//************************************************************************************
//************************************************************************************


double SmallDisplacementBeamElement3D2N::CalculateInternalShear(const double& Q, const double& ShearLoad, const double& X)
{
    return  -Q + ShearLoad * X;
}

//************************************************************************************
//************************************************************************************


double SmallDisplacementBeamElement3D2N::CalculateInternalAxil(const double& N, const double& AxialLoad, const double& X)
{
    return  -N + AxialLoad * X;

}


//************************************************************************************
//************************************************************************************


void SmallDisplacementBeamElement3D2N::CalculateLocalNodalStress(Vector& Stress, Vector& rVolumeForce)
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
int  SmallDisplacementBeamElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = Element::Check(rCurrentProcessInfo);

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(ROTATION);
    KRATOS_CHECK_VARIABLE_KEY(ANGULAR_VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ANGULAR_ACCELERATION);

    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(CROSS_SECTION_AREA);
    KRATOS_CHECK_VARIABLE_KEY(LOCAL_INERTIA_TENSOR);
    //KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION);

    // Check that the element nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
      {
	// Nodal data
	Node<3> &rNode = this->GetGeometry()[i];
	KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rNode);
	KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION,rNode);

	// Nodal dofs
	KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,rNode);
	KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,rNode);
	if( rCurrentProcessInfo[SPACE_DIMENSION] == 3)
	  KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z,rNode);

	KRATOS_CHECK_DOF_IN_NODE(ROTATION_Z,rNode);
	if( rCurrentProcessInfo[SPACE_DIMENSION] == 3){
	  KRATOS_CHECK_DOF_IN_NODE(ROTATION_X,rNode);
	  KRATOS_CHECK_DOF_IN_NODE(ROTATION_Y,rNode);
	}
      }

    //verify that the area is given by properties
    if ( this->GetProperties().Has(CROSS_SECTION_AREA) == false )
      {
        if( this->GetProperties()[CROSS_SECTION_AREA] == 0.0 )
	  KRATOS_ERROR << "CROSS_SECTION_AREA not provided for this element " << this->Id() << std::endl;
      }

    if ( this->GetProperties().Has(LOCAL_INERTIA_TENSOR) == false )
      {
        if( this->GetProperties()[LOCAL_INERTIA_TENSOR](1,1) == 0.0 )
	  KRATOS_ERROR << "LOCAL_INERTIA_TENSOR not provided for this element " << this->Id() << std::endl;

	// if ( LOCAL_CONSTITUTIVE_MATRIX.Key() == 0)
	//   KRATOS_ERROR <<  "LOCAL_CONSTITUTIVE_MATRIX has Key zero! (check if the application is correctly registered)" << std::endl;

      }

    if (GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size()!=2 )
    {
      KRATOS_THROW_ERROR( std::invalid_argument, "This element works only in 3D and with 2 noded linear elements", "")
    }


    return ErrorCode;


    KRATOS_CATCH( "" )
}




} // Namespace Kratos


