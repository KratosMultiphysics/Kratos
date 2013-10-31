//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/kratos_flags.h"
#include "custom_conditions/contact_domain_condition.hpp"

#include "pfem_solid_mechanics_application.h"


//#include <omp.h>

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ContactDomainCondition::ContactDomainCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
  //DO NOT ADD DOFS HERE!!!
   this->Set(CONTACT);
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ContactDomainCondition::ContactDomainCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
  this->Set(CONTACT);
  mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ContactDomainCondition::ContactDomainCondition( ContactDomainCondition const& rOther)
    :Condition(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    ,mContactVariables(rOther.mContactVariables)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ContactDomainCondition&  ContactDomainCondition::operator=(ContactDomainCondition const& rOther)
{
    Condition::operator=(rOther);

    mThisIntegrationMethod = rOther.mThisIntegrationMethod;

    mConstitutiveLawVector.clear();
    mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
    {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
    }

    mContactVariables = rOther.mContactVariables;

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Condition::Pointer ContactDomainCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer(new ContactDomainCondition( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************


ContactDomainCondition::~ContactDomainCondition()
{
}



//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

ContactDomainCondition::IntegrationMethod ContactDomainCondition::GetIntegrationMethod()
{
    return mThisIntegrationMethod;
}


//************************************************************************************
//************************************************************************************

void ContactDomainCondition::GetDofList( DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo )
{
    rConditionalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

        if ( GetGeometry().WorkingSpaceDimension() == 3 )
        {
            rConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

    //ADD MASTER NODE
    unsigned int vsize=GetValue(MASTER_NODES).size();
    Element::NodeType&    MasterNode   = GetValue(MASTER_NODES)[vsize-1];
    
    rConditionalDofList.push_back( MasterNode.pGetDof( DISPLACEMENT_X ) );
    rConditionalDofList.push_back( MasterNode.pGetDof( DISPLACEMENT_Y ) );
    if ( GetGeometry().WorkingSpaceDimension() == 3 )
        rConditionalDofList.push_back( MasterNode.pGetDof( DISPLACEMENT_Z ) );
}

//************************************************************************************
//************************************************************************************

void ContactDomainCondition::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    int number_of_nodes = GetGeometry().size();
    int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int mat_sizes = (number_of_nodes + 1)  * dimension;

    if ( rResult.size() != mat_sizes )
        rResult.resize( mat_sizes, false );

    for ( int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension;
        rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if ( dimension == 3 )
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

    //ADD MASTER NODE
    int index = number_of_nodes * dimension;
    unsigned int vsize=GetValue(MASTER_NODES).size();
    Element::NodeType&    MasterNode   = GetValue(MASTER_NODES)[vsize-1];
 
    rResult[index]   = MasterNode.GetDof( DISPLACEMENT_X ).EquationId();
    rResult[index+1] = MasterNode.GetDof( DISPLACEMENT_Y ).EquationId();
    if ( dimension == 3 )
        rResult[index+2] = MasterNode.GetDof( DISPLACEMENT_Z ).EquationId();

}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void ContactDomainCondition::GetValuesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int mat_size = (number_of_nodes + 1) * dimension;

    if ( rValues.size() != mat_size ) rValues.resize( mat_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        rValues[index] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
    }


    //ADD MASTER NODE
    unsigned int index = number_of_nodes * dimension;
    unsigned int vsize=GetValue(MASTER_NODES).size();
    Element::NodeType&    MasterNode   = GetValue(MASTER_NODES)[vsize-1];

    rValues[index] = MasterNode.GetSolutionStepValue( DISPLACEMENT_X, Step );
    rValues[index+1] = MasterNode.GetSolutionStepValue( DISPLACEMENT_Y, Step );

    if ( dimension == 3 )
        rValues[index+1] = MasterNode.GetSolutionStepValue( DISPLACEMENT_Z, Step );

}


//************************************VELOCITY****************************************
//************************************************************************************

void ContactDomainCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int mat_size = (number_of_nodes + 1) * dimension;

    if ( rValues.size() != mat_size ) rValues.resize( mat_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        rValues[index] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
    }

    //ADD MASTER NODE
    unsigned int index = number_of_nodes * dimension;
    unsigned int vsize=GetValue(MASTER_NODES).size();
    Element::NodeType&    MasterNode   = GetValue(MASTER_NODES)[vsize-1];

    rValues[index] = MasterNode.GetSolutionStepValue( VELOCITY_X, Step );
    rValues[index+1] = MasterNode.GetSolutionStepValue( VELOCITY_Y, Step );

    if ( dimension == 3 )
        rValues[index+1] = MasterNode.GetSolutionStepValue( VELOCITY_Z, Step );

}



//*********************************ACCELERATION***************************************
//************************************************************************************

void ContactDomainCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int mat_size = (number_of_nodes + 1) * dimension;

    if ( rValues.size() != mat_size ) rValues.resize( mat_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        rValues[index] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
    }

    //ADD MASTER NODE
    unsigned int index = number_of_nodes * dimension;
    unsigned int vsize=GetValue(MASTER_NODES).size();
    Element::NodeType&    MasterNode   = GetValue(MASTER_NODES)[vsize-1];

    rValues[index] = MasterNode.GetSolutionStepValue( ACCELERATION_X, Step );
    rValues[index+1] = MasterNode.GetSolutionStepValue( ACCELERATION_Y, Step );

    if ( dimension == 3 )
        rValues[index+1] = MasterNode.GetSolutionStepValue( ACCELERATION_Z, Step );


}


//************************************************************************************
//************************************************************************************

//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

void ContactDomainCondition::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }
}

//*********************************SET VECTOR VALUE***********************************
//************************************************************************************

void ContactDomainCondition::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

    for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable,rValues[PointNumber], rCurrentProcessInfo );
    }

}


//*********************************SET MATRIX VALUE***********************************
//************************************************************************************

void ContactDomainCondition::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable,rValues[PointNumber], rCurrentProcessInfo );
    }

}

//*********************************GET DOUBLE VALUE***********************************
//************************************************************************************

void ContactDomainCondition::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
        rValues.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(), false );

    for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
}


//**********************************GET VECTOR VALUE**********************************
//************************************************************************************

void ContactDomainCondition::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int& size = GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();

    if ( rValues.size() != size )
        rValues.resize( size );

    if ( rVariable == PK2_STRESS_TENSOR || rVariable == CAUCHY_STRESS_TENSOR)
    {
        for ( unsigned int PointNumber = 0;
                PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
                PointNumber++ )
        {
            rValues[PointNumber] =
                mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
        }
    }


    if ( rVariable == INTERNAL_VARIABLES )
    {
        for ( unsigned int PointNumber = 0;
                PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
                PointNumber++ )
        {
            rValues[PointNumber] =
                mConstitutiveLawVector[PointNumber]->GetValue( INTERNAL_VARIABLES, rValues[PointNumber] );

        }
    }

}

//***********************************GET MATRIX VALUE*********************************
//************************************************************************************

void ContactDomainCondition::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }

    if ( rVariable == PK2_STRESS_TENSOR || rVariable == CAUCHY_STRESS_TENSOR)
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }


}



//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void ContactDomainCondition::Initialize()
{
    KRATOS_TRY

	    std::cout<<" The position update on the iteration requires a modification in the condition "<<std::endl;

    KRATOS_CATCH( "" )
}



////************************************************************************************
////************************************************************************************

void ContactDomainCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{

    //0.- Initialize Iteration Counter
    mContactVariables.IterationCounter = 0;

    //1.- Clear nodal contact forces
    ClearNodalForces();
   
    //2.-Set Master Element Geometry: Master Elements and Nodes
    this->SetMasterGeometry();

    //3.-Get ConstitutiveLaw from the selected Master Element
    ElementType& MasterElement = mContactVariables.GetMasterElement();
    
    MasterElement.GetValueOnIntegrationPoints( CONSTITUTIVE_LAW_POINTER, mConstitutiveLawVector, rCurrentProcessInfo );
    
    //4.- Clear possible residual forces from the mesh refining and interpolation:
    ClearMasterElementNodalForces( MasterElement );


    //5.- Calculate Contact Factor (stabilization or penalty)
    this->CalculateContactFactor( rCurrentProcessInfo );

    //Previous Gap Calculation
    this->CalculatePreviousGap();

}

////************************************************************************************
////************************************************************************************

void ContactDomainCondition::InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
{
  //0.- Clear nodal contact forces
  ClearNodalForces();

  CurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] = 0;
  CurrentProcessInfo[NUMBER_OF_STICK_CONTACTS]  = 0;
  CurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS]   = 0;

}
//************************************************************************************
//************************************************************************************

void ContactDomainCondition::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
  KRATOS_TRY

  CurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] = 0;
  CurrentProcessInfo[NUMBER_OF_STICK_CONTACTS]  = 0;
  CurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS]   = 0;

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void ContactDomainCondition::ClearNodalForces()
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
	GetGeometry()[i].SetLock();
	VectorType & ContactForceNormal  = GetGeometry()[i].FastGetSolutionStepValue(FORCE_CONTACT_NORMAL);
	ContactForceNormal.clear();

	VectorType & ContactForceTangent  = GetGeometry()[i].FastGetSolutionStepValue(FORCE_CONTACT_TANGENT);
	ContactForceTangent.clear();
	GetGeometry()[i].UnSetLock();
    }


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void ContactDomainCondition::ClearMasterElementNodalForces(ElementType& rMasterElement)
{
    KRATOS_TRY
      
    
    //------------------------------------//
    const unsigned int number_of_nodes = rMasterElement.GetGeometry().PointsNumber();
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
	rMasterElement.GetGeometry()[i].SetLock();
	VectorType & ContactForceNormal   = rMasterElement.GetGeometry()[i].FastGetSolutionStepValue(FORCE_CONTACT_NORMAL);
	VectorType & ContactForceTangent  = rMasterElement.GetGeometry()[i].FastGetSolutionStepValue(FORCE_CONTACT_TANGENT);


	ContactForceNormal.clear();
	ContactForceTangent.clear();
	rMasterElement.GetGeometry()[i].UnSetLock();
    }
    //------------------------------------//
    KRATOS_CATCH( "" )
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//*****************************COMPUTE AVERAGE VALUE ON GAUSS POINT*******************
//************************************************************************************

void ContactDomainCondition::SetContactIntegrationVariable ( Vector & rContactVariable, std::vector<Vector> & rMasterVariables, const unsigned int& rPointNumber )
{
    
	//option if master element has one integration point:
	unsigned int master_integration_points_number  = rMasterVariables.size();
	unsigned int contact_integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
	if( master_integration_points_number == contact_integration_points_number ){
		rContactVariable = rMasterVariables[rPointNumber];
	}
	else{     //option if master element has several integration points:
		rContactVariable=rMasterVariables[0];
		for(unsigned int i=1; i<master_integration_points_number; i++)
		{
			rContactVariable += rMasterVariables[i];
		}
		rContactVariable *= (1.0/double(master_integration_points_number));
	}

}



//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void ContactDomainCondition::CalculateKinematics( GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo, const unsigned int& rPointNumber )
{
    KRATOS_TRY

    ElementType&  MasterElement  = mContactVariables.GetMasterElement();
    GeometryType& MasterGeometry = mContactVariables.GetMasterGeometry();
      
    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = MasterGeometry.ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    //Get integration points number
    unsigned int integration_points_number = MasterGeometry.IntegrationPointsNumber( MasterElement.GetIntegrationMethod() );

    unsigned int voigtsize = 3;
    int dimension = GetGeometry().WorkingSpaceDimension();
    if( dimension == 3 )
        voigtsize  = 6;

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);


    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ );

    //Compute cartesian derivatives
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber] , Invj );   


    //Get Current DeformationGradient
    std::vector<Matrix> DeformationGradientVector ( integration_points_number );
    DeformationGradientVector[rPointNumber]=identity_matrix<double>( dimension );
    MasterElement.GetValueOnIntegrationPoints(DEFORMATION_GRADIENT,DeformationGradientVector,rCurrentProcessInfo);
    rVariables.F = DeformationGradientVector[rPointNumber];

    rVariables.detF = MathUtils<double>::Det(rVariables.F);

    //Get Current Stress
    std::vector<Vector> StressVector ( integration_points_number );
    StressVector[rPointNumber]=ZeroVector(voigtsize);
    MasterElement.GetValueOnIntegrationPoints(PK2_STRESS_VECTOR,StressVector,rCurrentProcessInfo);
    
    SetContactIntegrationVariable( rVariables.StressVector, StressVector, rPointNumber );


    //std::cout<<" StressVector "<<rVariables.StressVector<<std::endl;
    
    //Get Current Strain
    std::vector<Matrix> StrainTensor ( integration_points_number );
    StrainTensor[rPointNumber]=ZeroMatrix(dimension,dimension);
    MasterElement.GetValueOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_TENSOR,StrainTensor,rCurrentProcessInfo);
    std::vector<Vector> StrainVector ( integration_points_number );
    for(unsigned int i=1; i<integration_points_number; i++)
    {
	    StrainVector[i] = MathUtils<double>::StrainTensorToVector( StrainTensor[i], voigtsize );
    }

    SetContactIntegrationVariable( rVariables.StrainVector, StrainVector, rPointNumber );


    //Get Current Constitutive Matrix
    std::vector<Matrix> ConstitutiveMatrix(mConstitutiveLawVector.size());   
    MasterElement.CalculateOnIntegrationPoints(CONSTITUTIVE_MATRIX,ConstitutiveMatrix,rCurrentProcessInfo);

    rVariables.ConstitutiveMatrix = ConstitutiveMatrix[rPointNumber];

    //Calculate Explicit Lagrange Multipliers or Penalty Factors
    this->CalculateExplicitFactors( rVariables, rCurrentProcessInfo );


    KRATOS_CATCH( "" )
}


//***********************COMPUTE LOCAL SYSTEM CONTRIBUTIONS***************************
//************************************************************************************


void ContactDomainCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{

    //calculation flags
    Flags CalculationFlags;
    CalculationFlags.Set(ContactDomainUtilities::COMPUTE_LHS_MATRIX);
    CalculationFlags.Set(ContactDomainUtilities::COMPUTE_RHS_VECTOR);


    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, CalculationFlags );

    //Calculate elemental system
    this->CalculateConditionalSystem( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculationFlags );

}

//************************************************************************************
//************************************************************************************


void ContactDomainCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    Flags CalculationFlags;
    CalculationFlags.Set(ContactDomainUtilities::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, CalculationFlags );

    //Calculate elemental system
    this->CalculateConditionalSystem( LeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculationFlags );

}


//************************************************************************************
//************************************************************************************


void ContactDomainCondition::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    Flags CalculationFlags;
    CalculationFlags.Set(ContactDomainUtilities::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, CalculationFlags );

    //Calculate elemental system
    this->CalculateConditionalSystem( rLeftHandSideMatrix, RightHandSideVector, rCurrentProcessInfo, CalculationFlags );

}

//************************************************************************************
//************************************************************************************

void ContactDomainCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
						      VectorType& rRightHandSideVector,
						      Flags& rCalculationFlags)
{

  const unsigned int number_of_nodes = GetGeometry().size();
  const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

  //resizing as needed the LHS
  unsigned int mat_size = (number_of_nodes + 1) * dimension;

  if ( rCalculationFlags.Is(ContactDomainUtilities::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
      if ( rLeftHandSideMatrix.size1() != mat_size )
	rLeftHandSideMatrix.resize( mat_size, mat_size, false );

      noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }


  //resizing as needed the RHS
  if ( rCalculationFlags.Is(ContactDomainUtilities::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
      if ( rRightHandSideVector.size() != mat_size )
	rRightHandSideVector.resize( mat_size, false );

      rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS
    }
}



//************************************************************************************
//************************************************************************************

void ContactDomainCondition::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType & MasterGeometry = mContactVariables.GetMasterGeometry();

    const unsigned int number_of_nodes = MasterGeometry.size();
    const unsigned int dimension       = MasterGeometry.WorkingSpaceDimension();

    unsigned int voigtsize = 3;
    if( dimension == 3 )
    {
        voigtsize  = 6;
    }

    rVariables.F.resize( dimension, dimension );

    rVariables.ConstitutiveMatrix.resize( voigtsize, voigtsize );

    rVariables.StressVector.resize( voigtsize );

    rVariables.DN_DX.resize( number_of_nodes, dimension );

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

}


//************************************************************************************
//************************************************************************************

void ContactDomainCondition::CalculateRelativeVelocity (GeneralVariables& rVariables, VectorType & TangentVelocity)
{
    //if current tangent is not previously computed, do it here.
    rVariables.Contact.CurrentSurface.Tangent = this->CalculateCurrentTangent( rVariables.Contact.CurrentSurface.Tangent );


   if(double(inner_prod(rVariables.Contact.CurrentSurface.Tangent,mContactVariables.ReferenceSurface.Tangent))<0) //to give the correct direction
        rVariables.Contact.CurrentSurface.Tangent*=-1;

   //std::cout<<" Normal ["<<this->Id()<<"] :"<<rVariables.Contact.CurrentSurface.Normal<<std::endl;
   //std::cout<<" Tangent ["<<this->Id()<<"] :"<<rVariables.Contact.CurrentSurface.Tangent<<std::endl;

    // (Tangent vector previously computed)
    const int number_of_nodes = GetGeometry().size();

    //compute relative velocities
    int slave=mContactVariables.slaves[0];
    VectorType  CurrentVelocity;
    for (int i = 0; i < number_of_nodes; i++ )
    {
        //Displacement from the reference to the current configuration
        CurrentVelocity  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        if(i!=slave)
	    CurrentVelocity *=(-1)*(1.0/double(number_of_nodes - 1));

        TangentVelocity+=CurrentVelocity;
    }

    //Relative tangent movement of the slave if the master is fixed (the direction is implicit in the method)
    TangentVelocity =  rVariables.Contact.CurrentSurface.Tangent*(inner_prod(TangentVelocity,rVariables.Contact.CurrentSurface.Tangent));

    //Filter for low velocities (relatives to dynamic waves)
    CurrentVelocity.clear();
    CalculateRelativeDisplacement(rVariables, CurrentVelocity);

    if(norm_2(TangentVelocity)<1e-2*(norm_2(CurrentVelocity)/norm_2(TangentVelocity)))
    {
      TangentVelocity.clear();
    }

    //std::cout<<" TangentVelocity ["<<this->Id()<<"] :"<<TangentVelocity<<std::endl;

}

//************************************************************************************
//************************************************************************************


void ContactDomainCondition::CalculateRelativeDisplacement (GeneralVariables& rVariables, VectorType & TangentDisplacement)
{

    // (Tangent vector previously computed)
    const int number_of_nodes = GetGeometry().size();

    //compute relative displacements
    int slave=mContactVariables.slaves[0];
    VectorType CurrentDisplacement;
    for (int i = 0; i < number_of_nodes; i++ )
    {
        //Displacement from the reference to the current configuration
        CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        if(i!=slave)
	    CurrentDisplacement *=(-1)*(1.0/double(number_of_nodes - 1));

        TangentDisplacement+=CurrentDisplacement;
    }

    //Relative tangent movement of the slave if the master is fixed (the direction is implicit in the method)
    TangentDisplacement = rVariables.Contact.CurrentSurface.Tangent*(inner_prod(TangentDisplacement,rVariables.Contact.CurrentSurface.Tangent));

}

//************************************************************************************
//************************************************************************************


void ContactDomainCondition::CalculateFrictionCoefficient (GeneralVariables& rVariables, const VectorType & TangentVelocity)
{
    //---FRICTION LAW in function of the relative sliding velocity ---//

    //Addicional constitutive parameter  paramC
    //which describes how fast the static coefficient approaches the dynamic:
    double paramC=0.1;

    //Addicional constitutive parameter  paramE
    //regularization parameter (->0, classical Coulomb law)
    double paramE=0.01;

    double Velocity=norm_2(TangentVelocity);

    double muStatic  = 0.3;
    double muDynamic = 0.2;
    muStatic  =GetProperties()[MU_STATIC];
    muDynamic =GetProperties()[MU_DYNAMIC];
   

    if(Velocity!=0)
    {


        rVariables.Contact.FrictionCoefficient=muDynamic+(muStatic-muStatic)*exp((-1)*paramC*fabs(Velocity));


        //if (TangentVelocity.modulus()>paramE){

        //1.- square root regularization
        rVariables.Contact.FrictionCoefficient*=fabs(Velocity)/sqrt((Velocity*Velocity)+(paramE*paramE));//square root

        //2.-hyperbolic regularization
        //FrictionCoefficient*=tanh(fabs(Velocity)/paramE);

        //}


    }
    else
    {
	rVariables.Contact.FrictionCoefficient=muStatic;
    }

    //Activate or deactivate friction (simulation type)

    if( rVariables.Contact.Options.IsNot(ContactDomainUtilities::COMPUTE_FRICTION_FORCES) ){
	    rVariables.Contact.FrictionCoefficient = 0;
    }

    //std::cout<<" friction coefficient "<<rVariables.Contact.FrictionCoefficient<<std::endl;

}

//************************************************************************************
//************************************************************************************

void ContactDomainCondition::CalculateConditionalSystem( MatrixType& rLeftHandSideMatrix,
							 VectorType& rRightHandSideVector,
							 ProcessInfo& rCurrentProcessInfo,
							 Flags& rCalculationFlags )
{
    KRATOS_TRY

    //std::cout<<"//******** CONTACT ELEMENT "<<this->Id()<<" ********// "<<std::endl;
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

    
    //SET TANGENT DIRECTION-RELATIVE VELOCITY AND FRICTION PARAMETERS

    //Initialize friction parameter
    Variables.Contact.FrictionCoefficient   =0;

    Variables.Contact.Options.Set(ContactDomainUtilities::COMPUTE_FRICTION_STIFFNESS);
    if( GetProperties()[FRICTION_ACTIVE] == 1 )
	    Variables.Contact.Options.Set(ContactDomainUtilities::COMPUTE_FRICTION_FORCES);

    //Tangent velocity and stablish friction parameter
    VectorType TangentVelocity (3,0.0);

    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        //Calculate Element Kinematics
        this->CalculateKinematics( Variables, rCurrentProcessInfo, PointNumber );
	
	//Calculate Functions for the Contact Tangent Matrix construction
	this->CalculateDomainShapeN(Variables);

	//Calculate Relative Velocity
	this->CalculateRelativeVelocity    ( Variables, TangentVelocity);

	//Calculate Friction Coefficient
	this->CalculateFrictionCoefficient ( Variables, TangentVelocity);

        double IntegrationWeight =0.5 * Variables.Contact.ReferenceBase[0].L;  //all components are multiplied by this
        IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

        if(Variables.Contact.Options.Is(ACTIVE))
        {
	  rCalculationFlags.Set(ContactDomainUtilities::COMPUTE_LHS_MATRIX,true); //take a look on strategy and impose it

	  if ( rCalculationFlags.Is(ContactDomainUtilities::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
	    {
	      //contributions to stiffness matrix calculated on the reference config
	      this->CalculateAndAddLHS ( rLeftHandSideMatrix, Variables, IntegrationWeight );
	    }

	  if ( rCalculationFlags.Is(ContactDomainUtilities::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
	    {
	      //contribution to contact forces
	      this->CalculateAndAddRHS ( rRightHandSideVector, Variables, IntegrationWeight );

	    }


        }
    }


    KRATOS_CATCH( "" )
}


//***********************************************************************************
//************************************************************************************

double& ContactDomainCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    return rIntegrationWeight;
}


//************************************************************************************
//************************************************************************************

void ContactDomainCondition::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
{
  //contributions to stiffness matrix calculated on the reference config
  this->CalculateAndAddKm( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

  //KRATOS_WATCH(rLeftHandSideMatrix)
}


//************************************************************************************
//************************************************************************************

void ContactDomainCondition::CalculateAndAddRHS(VectorType& rRightHandSideVector, GeneralVariables& rVariables, double& rIntegrationWeight)
{
  //contribution to contact forces
  this->CalculateAndAddContactForces(rRightHandSideVector, rVariables, rIntegrationWeight);

  //KRATOS_WATCH(rRightHandSideVector)
}

//************************************************************************************
//************************************************************************************

inline void ContactDomainCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
								 GeneralVariables &rVariables,
								 double& rIntegrationWeight)
{
    KRATOS_TRY

    //contributions to stiffness matrix calculated on the reference config
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int size      = rRightHandSideVector.size()/dimension;

    Vector Nforce;
    Vector Tforce;

    VectorType NormalForce  (3,0.0);
    VectorType TangentForce (3,0.0);

    unsigned int index=0;
    for (unsigned int ndi=0; ndi<size; ndi++)
    {
	Nforce=ZeroVector( dimension );
        Tforce=ZeroVector( dimension );

	NormalForce.clear();
	TangentForce.clear();

        for (unsigned int i=0; i<dimension; i++)
        {
            //NORMAL FORCE
   	    this->CalculateNormalForce(Nforce[i],rVariables,ndi,i);

            //TANGENT FORCE
            if(rVariables.Contact.Options.Is(NOT_SLIP))
            {
		this->CalculateTangentStickForce(Tforce[i],rVariables,ndi,i);
            }
            else
            {
                this->CalculateTangentSlipForce(Tforce[i],rVariables,ndi,i);
            }

	    rRightHandSideVector[index] -=(Nforce[i] + Tforce[i]);
 
	    NormalForce[i]  -= (Nforce[i])*rIntegrationWeight;
	    TangentForce[i] -= (Tforce[i])*rIntegrationWeight;

	    index++;
        }
	

	if(ndi<GetGeometry().PointsNumber())
	{
		GetGeometry()[ndi].SetLock();

		VectorType & ForceContactNormal  = GetGeometry()[ndi].FastGetSolutionStepValue(FORCE_CONTACT_NORMAL);

		VectorType & ForceContactTangent = GetGeometry()[ndi].FastGetSolutionStepValue(FORCE_CONTACT_TANGENT);


		//Only tangent force
		ForceContactTangent += (TangentForce);

		//Only normal force
		ForceContactNormal  += (NormalForce);

		GetGeometry()[ndi].UnSetLock();
    		
	}

    }


    rRightHandSideVector  *=  rIntegrationWeight;


    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************


void ContactDomainCondition::CalculateAndAddKm(MatrixType& rLeftHandSideMatrix,
					       GeneralVariables& rVariables,
					       double& rIntegrationWeight)
  
{ 
    KRATOS_TRY

    //contributions to stiffness matrix calculated on the reference config
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int size      = rLeftHandSideMatrix.size1()/dimension;
    double kcont=0;

    for (unsigned int ndi=0; ndi<size; ndi++)
    {      
        for (unsigned int ndj=0; ndj<size; ndj++)
        {
            for (unsigned int i=0; i<dimension; i++)
            {
                for (unsigned int j=0; j<dimension; j++)
                {
		    kcont=0;
		    this->CalcContactStiffness(kcont,rVariables,ndi,ndj,i,j);
		    rLeftHandSideMatrix(ndi*2+i,ndj*2+j)+=kcont;
                }
            }
        }
    }

    rLeftHandSideMatrix *= rIntegrationWeight;

    // std::cout<<std::endl;
    // std::cout<<" Kcontact "<<rLeftHandSideMatrix<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void ContactDomainCondition::CalculateOnIntegrationPoints( const Variable<double>& rVariable, Vector& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rOutput.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
        rOutput.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(), false );

    for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
}


//************************************************************************************
//************************************************************************************

void ContactDomainCondition::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    unsigned int StrainSize;

    if ( GetGeometry().WorkingSpaceDimension() == 2 )
    {
        StrainSize = 3;
    }
    else
    {
        StrainSize = 6;
    }

    Vector StrainVector( StrainSize );

    if ( rVariable == PK2_STRESS_TENSOR )
    {
        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        {
            if ( rOutput[ii].size() != StrainVector.size() )
                rOutput[ii].resize( StrainVector.size(), false );

            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable , rOutput[ii] );
        }
    }
    else
    {
        if ( rOutput.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
            rOutput.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );

        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
    }

}

//************************************************************************************
//************************************************************************************

void ContactDomainCondition::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int StrainSize;

    if ( dimension == 2 )
        StrainSize = 3;
    else
        StrainSize = 6;

 
    if ( rVariable == PK2_STRESS_TENSOR )
    {

        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {

            if ( rOutput[PointNumber].size2() != StrainSize )
                rOutput[PointNumber].resize( 1 ,  StrainSize , false );

	    //Get Current Stress
	    ElementType& MasterElement = mContactVariables.GetMasterElement();

	    std::vector<Vector> StressVector ( mConstitutiveLawVector.size() );
	    StressVector[PointNumber]=ZeroVector(StrainSize);

	    MasterElement.GetValueOnIntegrationPoints(PK2_STRESS_VECTOR,StressVector,rCurrentProcessInfo);

            for ( unsigned int ii = 0; ii < StressVector[PointNumber].size(); ii++ )
            {
                rOutput[PointNumber]( 0, ii ) = StressVector[PointNumber][ii];
            }


        }
    }
    if ( rVariable == CAUCHY_STRESS_TENSOR )
    {

        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            if ( rOutput[PointNumber].size2() != StrainSize)
                rOutput[PointNumber].resize( 1 , StrainSize , false );

            Matrix StressMatrix ( dimension, dimension );
            StressMatrix = mConstitutiveLawVector[PointNumber]->GetValue( rVariable , StressMatrix );

            Vector StressVector ( StrainSize );
            StressVector = MathUtils<double>::StressTensorToVector( StressMatrix );

            for ( unsigned int ii = 0; ii < StressVector.size(); ii++ )
            {
                rOutput[PointNumber]( 0, ii ) = StressVector[ii];
            }


        }

    }

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
int  ContactDomainCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();



    //verify that the variables are correctly initialized

    if ( VELOCITY.Key() == 0 )
        KRATOS_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" );

    if ( DISPLACEMENT.Key() == 0 )
        KRATOS_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" );

    if ( ACCELERATION.Key() == 0 )
        KRATOS_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" );

    if ( THICKNESS.Key() == 0 )
        KRATOS_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "" );

    //verify that the dofs exist
    for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            KRATOS_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() );

        if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
            KRATOS_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() );
    }

    //verify that the constitutive law exists
    if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
    {
        KRATOS_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() );
    }


    //verify that the constitutive law has the correct dimension
    if ( dimension == 2 )
    {
        if ( this->GetProperties().Has( THICKNESS ) == false )
            KRATOS_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() );

        if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 3 )
            KRATOS_ERROR( std::logic_error, "wrong constitutive law used. This is a 2D element! expected strain size is 3 (el id = ) ", this->Id() );
    }
    else
    {
        if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6 )
            KRATOS_ERROR( std::logic_error, "wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ", this->Id() );
    }

    //check constitutive law
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
    {
        return mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
    }

    //check if it is in the XY plane for 2D case


    return 0;

    KRATOS_CATCH( "" );
}

//Note: in the restart the contact mesh is generated from the begining

void ContactDomainCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    // int IntMethod = (int)mThisIntegrationMethod;
    // rSerializer.save("IntegrationMethod",IntMethod);
    // rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    // rSerializer.save("ContactVariables",mContactVariables);
}

void ContactDomainCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    // int IntMethod;
    // rSerializer.load("IntegrationMethod",IntMethod);
    // mThisIntegrationMethod = IntegrationMethod(IntMethod);
    // rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    // rSerializer.load("ContactVariables",mContactVariables);
}



} // Namespace Kratos


