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
#include "custom_conditions/force_load_condition.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

/**
 * Flags related to the condition computation
 */
KRATOS_CREATE_LOCAL_FLAG( ForceLoadCondition, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( ForceLoadCondition, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( ForceLoadCondition, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( ForceLoadCondition, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );


//***********************************************************************************
//***********************************************************************************
ForceLoadCondition::ForceLoadCondition()
    : Condition()
{
  //DO NOT CALL IT: only needed for Register and Serialization!!!
}


//***********************************************************************************
//***********************************************************************************
ForceLoadCondition::ForceLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
ForceLoadCondition::ForceLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{

    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
ForceLoadCondition::ForceLoadCondition( ForceLoadCondition const& rOther )
    : Condition(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
      
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer ForceLoadCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new ForceLoadCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer ForceLoadCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  std::cout<<" Call base class FORCE LOAD CONDITION Clone "<<std::endl;
  
  ForceLoadCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

  return Condition::Pointer( new ForceLoadCondition(NewCondition) );
}


//***********************************************************************************
//***********************************************************************************
ForceLoadCondition::~ForceLoadCondition()
{
}

//************* GETTING METHODS

//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::GetDofList(DofsVectorType& rConditionDofList,
				    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rConditionDofList.resize(0);
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
	if( dimension == 3 )
	  rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }


    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::EquationIdVector(EquationIdVectorType& rResult,
					  ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int condition_size        = number_of_nodes * dimension;

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        int index = i * dimension;
        rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
	if( dimension == 3)
	  rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::GetValuesVector(Vector& rValues, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size  = number_of_nodes * dimension;

    if ( rValues.size() != condition_size ) 
      rValues.resize( condition_size, false );

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int index = i * dimension;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
    }
}

//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size    = number_of_nodes * dimension;

    if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
    }
}


//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size    = number_of_nodes * dimension;

    if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
    }

}


//************************************************************************************
//************************************************************************************
void ForceLoadCondition::ClearNodalForces()
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	if( GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_FORCE) && GetGeometry()[i].SolutionStepsDataHas(INTERNAL_FORCE) ){
	  
	  array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
	  array_1d<double, 3 > & InternalForce = GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);
  
	  GetGeometry()[i].SetLock();
	  ExternalForce.clear();
	  InternalForce.clear();
	  GetGeometry()[i].UnSetLock();

	}

      }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::AddExplicitContribution(const VectorType& rRHS, 
						 const Variable<VectorType>& rRHSVariable,
						 Variable<array_1d<double,3> >& rDestinationVariable, 
						 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == EXTERNAL_FORCES_VECTOR && rDestinationVariable == EXTERNAL_FORCE )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		ExternalForce[j] += rRHS[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		ForceResidual[j] += rRHS[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    KRATOS_CATCH( "" )
}

//************* STARTING - ENDING  METHODS
//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::Initialize()
{
    KRATOS_TRY

    mEnergy = 0;

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    ClearNodalForces();

    mEnergy = 0;
 
    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    mEnergy = 0;

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
						  VectorType& rRightHandSideVector,
						  Flags& rCalculationFlags)

{

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension;

    if ( rCalculationFlags.Is(ForceLoadCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(ForceLoadCondition::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
	    rRightHandSideVector.resize( MatSize, false );
      
	rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
	  
    }
}


//************************************************************************************
//************************************************************************************

void ForceLoadCondition::InitializeGeneralVariables(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{

    //const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    //const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rVariables.DomainSize = 1;

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void ForceLoadCondition::CalculateKinematics(GeneralVariables& rVariables,
					     const double& rPointNumber)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateKinematics method for a force load condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

Vector& ForceLoadCondition::CalculateVectorForce(Vector& rVectorForce, GeneralVariables& rVariables)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateVectorForce method for a force load condition ... illegal operation!!", "" )


    return rVectorForce;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void ForceLoadCondition::CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
						  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create and initialize condition variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //force terms
    Vector VectorForce;

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

        //calculating weights for integration on the "reference configuration" Jacobian respect to the reference configuration
	//take in account in a linear element (Jacobian=2*Area) this is the relation

        double IntegrationWeight = Variables.Jacobian * integration_points[PointNumber].Weight();

        IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

	//std::cout<<" Variables.Jacobian "<<Variables.Jacobian<<" Weight "<<integration_points[PointNumber].Weight()<<" / "<<std::endl;

       
	//calculation of the force and the pressure loads
	VectorForce = this->CalculateVectorForce( VectorForce, Variables );

        if ( rLocalSystem.CalculationFlags.Is(ForceLoadCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
        {
            //contributions to stiffness matrix calculated on the reference config
	    this->CalculateAndAddLHS ( rLocalSystem, Variables, IntegrationWeight );
        }

        if ( rLocalSystem.CalculationFlags.Is(ForceLoadCondition::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
        {
            //contribution to external forces 
	    this->CalculateAndAddRHS ( rLocalSystem, Variables, VectorForce, IntegrationWeight );
        }

    }

    KRATOS_CATCH( "" )
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void ForceLoadCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
{

  //contributions of the stiffness matrix calculated on the reference configuration
  if( rLocalSystem.CalculationFlags.Is( ForceLoadCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
    {
      std::vector<MatrixType>& rLeftHandSideMatrices = rLocalSystem.GetLeftHandSideMatrices();
      const std::vector< Variable< MatrixType > >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables();

      for( unsigned int i=0; i<rLeftHandSideVariables.size(); i++ )
	{
	  bool calculated = false;
	  
	  if( rLeftHandSideVariables[i] == GEOMETRIC_STIFFNESS_MATRIX ){
	    // operation performed: add Kg to the rLefsHandSideMatrix
	    this->CalculateAndAddKuug( rLeftHandSideMatrices[i], rVariables, rIntegrationWeight );
	    calculated = true;
	  }

	  if(calculated == false)
	    {
	      KRATOS_THROW_ERROR( std::logic_error, " CONDITION can not supply the required local system variable: ",rLeftHandSideVariables[i] )
	    }

	}
    } 
  else{
    
    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix(); 

    // operation performed: add Kg to the rLefsHandSideMatrix
    this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    //KRATOS_WATCH( rLeftHandSideMatrix )
  }

}


//************************************************************************************
//************************************************************************************

void ForceLoadCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVectorForce, double& rIntegrationWeight)
{
    //contribution of the internal and external forces
    if( rLocalSystem.CalculationFlags.Is( ForceLoadCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
    {

      std::vector<VectorType>& rRightHandSideVectors = rLocalSystem.GetRightHandSideVectors();
      const std::vector< Variable< VectorType > >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables();
      for( unsigned int i=0; i<rRightHandSideVariables.size(); i++ )
	{
	  bool calculated = false;
	  if( rRightHandSideVariables[i] == EXTERNAL_FORCES_VECTOR ){
	    // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
	    this->CalculateAndAddExternalForces( rRightHandSideVectors[i], rVariables, rVectorForce, rIntegrationWeight );
	    calculated = true;
	  }

	  if( rRightHandSideVariables[i] == CONTACT_FORCES_VECTOR ){
	    // operation performed: rRightHandSideVector += ContactForce*IntToReferenceWeight
	    rRightHandSideVectors[i] += ZeroVector( rRightHandSideVectors[i].size() );
	    calculated = true;
	  }
	  
	  if(calculated == false)
	    {
	      KRATOS_THROW_ERROR( std::logic_error, " CONDITION can not supply the required local system variable: ",rRightHandSideVariables[i] )
	    }

	}
    }
    else{
      
      VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 

      // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
      this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVectorForce, rIntegrationWeight );

      //KRATOS_WATCH( rRightHandSideVector )

    }
    

}

//***********************************************************************************
//************************************************************************************

double& ForceLoadCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    return rIntegrationWeight;
}


//************************************************************************************
//************************************************************************************

void ForceLoadCondition::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(ForceLoadCondition::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(RightHandSideVector);

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );

    //KRATOS_WATCH( rRightHandSideVector )

}

//************************************************************************************
//************************************************************************************

void ForceLoadCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(ForceLoadCondition::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );

}

//************************************************************************************
//************************************************************************************

void ForceLoadCondition::CalculateRightHandSide( std::vector< VectorType >& rRightHandSideVectors, const std::vector< Variable< VectorType > >& rRHSVariables, ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(ForceLoadCondition::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(ForceLoadCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    if( rRHSVariables.size() != rRightHandSideVectors.size() )
      rRightHandSideVectors.resize(rRHSVariables.size());
    
    for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
      {
	this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVectors[i], LocalSystem.CalculationFlags );
      }

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );


}



//************************************************************************************
//************************************************************************************

void ForceLoadCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(ForceLoadCondition::COMPUTE_LHS_MATRIX);
    LocalSystem.CalculationFlags.Set(ForceLoadCondition::COMPUTE_RHS_VECTOR);

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );

    //KRATOS_WATCH( rLeftHandSideMatrix )
    //KRATOS_WATCH( rRightHandSideVector )

}


//************************************************************************************
//************************************************************************************

void ForceLoadCondition::CalculateLocalSystem( std::vector< MatrixType >& rLeftHandSideMatrices,
					       const std::vector< Variable< MatrixType > >& rLHSVariables,
					       std::vector< VectorType >& rRightHandSideVectors,
					       const std::vector< Variable< VectorType > >& rRHSVariables,
					       ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(ForceLoadCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
    LocalSystem.CalculationFlags.Set(ForceLoadCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);


    //Initialize sizes for the system components:
    if( rLHSVariables.size() != rLeftHandSideMatrices.size() )
      rLeftHandSideMatrices.resize(rLHSVariables.size());

    if( rRHSVariables.size() != rRightHandSideVectors.size() )
      rRightHandSideVectors.resize(rRHSVariables.size());
    
    LocalSystem.CalculationFlags.Set(ForceLoadCondition::COMPUTE_LHS_MATRIX);
    for( unsigned int i=0; i<rLeftHandSideMatrices.size(); i++ )
      {
	//Note: rRightHandSideVectors.size() > 0
	this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0], LocalSystem.CalculationFlags );
      }

    LocalSystem.CalculationFlags.Set(ForceLoadCondition::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(ForceLoadCondition::COMPUTE_LHS_MATRIX,false);

    for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
      {
	//Note: rLeftHandSideMatrices.size() > 0
    	this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags );
      }
    LocalSystem.CalculationFlags.Set(ForceLoadCondition::COMPUTE_LHS_MATRIX,true);


    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrices(rLeftHandSideMatrices);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetLeftHandSideVariables(rLHSVariables);
    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );

}


//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rMassMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rDampingMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
					     GeneralVariables& rVariables,
					     double& rIntegrationWeight)

{
    KRATOS_TRY



    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
						       GeneralVariables& rVariables,
						       Vector& rVectorForce,
						       double& rIntegrationWeight)

{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Energy Calculation:
    Vector CurrentValueVector = ZeroVector(dimension); 
    Vector Displacements = ZeroVector(dimension);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	//current displacements to compute energy
	CurrentValueVector = GetCurrentValue( DISPLACEMENT, CurrentValueVector, i );

	Displacements += rVariables.N[i] * Displacements;
      }
    //------

    Vector ForceVector = ZeroVector(3);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = dimension * i;
	
        for ( unsigned int j = 0; j < dimension; j++ )
        {
	  rRightHandSideVector[index + j] += rVariables.N[i] * rVectorForce[j] * rIntegrationWeight;
        }

	ForceVector += rVariables.N[i] * rVectorForce * rIntegrationWeight;
    }


    mEnergy += inner_prod( ForceVector, Displacements );


    //KRATOS_WATCH( rRightHandSideVector )

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void ForceLoadCondition::GetNodalDeltaMovements(Vector& rValues, const int& rNode)
{
  unsigned int dimension = GetGeometry().WorkingSpaceDimension();

  if( rValues.size() != dimension )
    rValues.resize(dimension);

  rValues = ZeroVector(dimension);
  
  Vector CurrentValueVector = ZeroVector(3);
  CurrentValueVector = GetCurrentValue( DISPLACEMENT, CurrentValueVector, rNode );

  Vector PreviousValueVector = ZeroVector(3);
  CurrentValueVector = GetPreviousValue( DISPLACEMENT, CurrentValueVector, rNode );


  rValues[0] = CurrentValueVector[0] - PreviousValueVector[0];
  rValues[1] = CurrentValueVector[1] - PreviousValueVector[1];

  if( dimension == 3 )
    rValues[2] = CurrentValueVector[2] - PreviousValueVector[2];
	
  //take imposed values away
  // if( (GetGeometry()[rNode].pGetDof(DISPLACEMENT_X))->IsFixed() )
  //   rValues[0] = 0;
  // if( (GetGeometry()[rNode].pGetDof(DISPLACEMENT_Y))->IsFixed() )
  //   rValues[1] = 0;

  // if( dimension == 3 )
  //   if( (GetGeometry()[rNode].pGetDof(DISPLACEMENT_Z))->IsFixed() )
  //     rValues[2] = 0;

}


//************************************************************************************
//************************************************************************************

Vector& ForceLoadCondition::GetCurrentValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode)
{
    KRATOS_TRY

    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    array_1d<double,3> ArrayValue;
    ArrayValue = GetGeometry()[rNode].FastGetSolutionStepValue( rVariable );
    
    if( rValue.size() != dimension )
      rValue.resize(dimension, false);

    for( unsigned int i=0; i<dimension; i++ )
      {
	rValue[i] = ArrayValue[i];
      }
   

    return rValue;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

Vector& ForceLoadCondition::GetPreviousValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode)
{
    KRATOS_TRY

    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    array_1d<double,3> ArrayValue;
    ArrayValue = GetGeometry()[rNode].FastGetSolutionStepValue( rVariable, 1 );
    
    if( rValue.size() != dimension )
      rValue.resize(dimension, false);

    for( unsigned int i=0; i<dimension; i++ )
      {
	rValue[i] = ArrayValue[i];
      }
   

    return rValue;

    KRATOS_CATCH( "" )
}

//*********************************GET DOUBLE VALUE***********************************
//************************************************************************************

void ForceLoadCondition::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
						      std::vector<double>& rValues,
						      const ProcessInfo& rCurrentProcessInfo )
{ 
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

//************************************************************************************
//************************************************************************************

void ForceLoadCondition::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    unsigned int integration_points = 0;
    if( integration_points_number == 0 )
      integration_points = 1;

    if ( rOutput.size() != integration_points )
        rOutput.resize( integration_points, false );


    if ( rVariable == EXTERNAL_ENERGY )
    {
      
      //reading integration points
      for ( unsigned int PointNumber = 0; PointNumber < integration_points; PointNumber++ )
        {
	  rOutput[PointNumber] = mEnergy; //fabs(mEnergy);
	}
    }

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************


int ForceLoadCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************

void ForceLoadCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
}

void ForceLoadCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
}


} // Namespace Kratos.
