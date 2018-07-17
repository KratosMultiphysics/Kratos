//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/kratos_flags.h"
#include "custom_conditions/thermal_contact_domain_condition.hpp"

#include "contact_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ThermalContactDomainCondition::ThermalContactDomainCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
    this->Set(THERMAL);
    this->Set(CONTACT);
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ThermalContactDomainCondition::ThermalContactDomainCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    this->Set(THERMAL);
    this->Set(CONTACT);
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ThermalContactDomainCondition::ThermalContactDomainCondition( ThermalContactDomainCondition const& rOther)
    :Condition(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    ,mContactVariables(rOther.mContactVariables)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ThermalContactDomainCondition&  ThermalContactDomainCondition::operator=(ThermalContactDomainCondition const& rOther)
{
    Condition::operator=(rOther);

    mThisIntegrationMethod = rOther.mThisIntegrationMethod;

    mConstitutiveLawVector.clear();
    mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

    for(unsigned int i=0; i<<mConstitutiveLawVector.size(); i++)
    {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
    }

    mContactVariables = rOther.mContactVariables;


    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Condition::Pointer ThermalContactDomainCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared<ThermalContactDomainCondition>(NewId, GetGeometry().Create( ThisNodes ), pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************


Condition::Pointer ThermalContactDomainCondition::Clone( IndexType NewId, NodesArrayType const& ThisNodes ) const
{
  return this->Create(NewId, ThisNodes, pGetProperties());
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************


ThermalContactDomainCondition::~ThermalContactDomainCondition()
{
}



//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

ThermalContactDomainCondition::IntegrationMethod ThermalContactDomainCondition::GetIntegrationMethod()
{
    return mThisIntegrationMethod;
}


//************************************************************************************
//************************************************************************************

void ThermalContactDomainCondition::GetDofList( DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo )
{
    rConditionalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rConditionalDofList.push_back( GetGeometry()[i].pGetDof( TEMPERATURE ) );
    }
}

//************************************************************************************
//************************************************************************************

void ThermalContactDomainCondition::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    int number_of_nodes = GetGeometry().size();
    unsigned int elementdimension = number_of_nodes;

   if ( rResult.size() != elementdimension )
        rResult.resize( elementdimension, false );

    for ( int i = 0; i < number_of_nodes; i++ )
    {
        rResult[i] = GetGeometry()[i].GetDof( TEMPERATURE ).EquationId();
    }


}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void ThermalContactDomainCondition::GetValuesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes;

    if ( rValues.size() != MatSize ) rValues.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        rValues[i] = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE, Step );
    }

}


//************************************VELOCITY****************************************
//************************************************************************************

void ThermalContactDomainCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes;

    if ( rValues.size() != MatSize ) rValues.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        rValues[i] = 0;
    }
}



//*********************************ACCELERATION***************************************
//************************************************************************************

void ThermalContactDomainCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes;

    if ( rValues.size() != MatSize ) rValues.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        rValues[i] = 0;
    }
}


//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

void ThermalContactDomainCondition::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

    for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable,rValues[PointNumber], rCurrentProcessInfo );
    }

}

//*********************************SET VECTOR VALUE***********************************
//************************************************************************************

void ThermalContactDomainCondition::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
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

void ThermalContactDomainCondition::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
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

void ThermalContactDomainCondition::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
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

void ThermalContactDomainCondition::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int& size = GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();

    if ( rValues.size() != size )
        rValues.resize( size );

    for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
        {
            rValues[PointNumber] =
                mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
        }


}

//***********************************GET MATRIX VALUE*********************************
//************************************************************************************

void ThermalContactDomainCondition::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

}



//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void ThermalContactDomainCondition::Initialize()
{
    KRATOS_TRY

	    std::cout<<" The position update on the iteration requires a modification in the condition "<<std::endl;

    KRATOS_CATCH( "" )
}



////************************************************************************************
////************************************************************************************

void ThermalContactDomainCondition::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
{

    //0.- Initialize Iteration Counter

    //1.-Set Master Element Geometry
    this->SetMasterGeometry();

    //3.-Get ConstitutiveLaw from the selected Master Element
    ElementType& MasterElement = mContactVariables.GetMasterElement();

    MasterElement.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,mConstitutiveLawVector,CurrentProcessInfo);


    //4.- Calculate Contact Factor (stabilization or penalty)
    CalculateHeatConductivity();

}

////************************************************************************************
////************************************************************************************

void ThermalContactDomainCondition::InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
{
}


//************************************************************************************
//************************************************************************************

void ThermalContactDomainCondition::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//****************************COMPUTE THERMAL CONDUCTIVITY****************************
//************************************************************************************


void ThermalContactDomainCondition::CalculateHeatConductivity()
{
    //Initilialize Tau for the stabilization
    double alpha_stab=1;

    // unsigned int vsize=GetValue(MASTER_ELEMENTS).size();
    // Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS)[vsize-1];
    Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS).back();

    //Look at the nodes, get the slave and get the Emin

    //Contact face segment node1-node2
    unsigned int slave=mContactVariables.slaves.back();

    double Kslave=GetGeometry()[slave].GetValue(NEIGHBOUR_ELEMENTS)[0].GetProperties()[HEAT_CONDUCTIVITY];
    double Kmin  =MasterElement.GetProperties()[HEAT_CONDUCTIVITY];

    if(Kmin>Kslave)
	Kmin=Kslave;

    if(Kmin==0){
      Kmin = 1e-4;
    }
    else{
      Kmin*= 1e-4;
    }

    //std::cout<<" Kslave "<<Kslave<<" Kmin "<<Kmin<<std::endl;

    mContactVariables.StabilizationFactor= alpha_stab/Kmin;

    //this must be supplied by properties:
    double HeatTransferCoeffitient = 5.0e6;
    mContactVariables.StabilizationFactor= alpha_stab * HeatTransferCoeffitient;


}



//***********************COMPUTE LOCAL SYSTEM CONTRIBUTIONS***************************
//************************************************************************************


void ThermalContactDomainCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
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


void ThermalContactDomainCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
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


void ThermalContactDomainCondition::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
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

void ThermalContactDomainCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
							     VectorType& rRightHandSideVector,
							     Flags& rCalculationFlags)
{

  const unsigned int number_of_nodes = GetGeometry().size();

  //resizing as needed the LHS
  unsigned int MatSize = (number_of_nodes);

  if ( rCalculationFlags.Is(ContactDomainUtilities::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
      if ( rLeftHandSideMatrix.size1() != MatSize )
	rLeftHandSideMatrix.resize( MatSize, MatSize, false );

      noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


  //resizing as needed the RHS
  if ( rCalculationFlags.Is(ContactDomainUtilities::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
      if ( rRightHandSideVector.size() != MatSize )
	rRightHandSideVector.resize( MatSize, false );

      rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
    }


}


//************************************************************************************
//************************************************************************************

void ThermalContactDomainCondition::CalculateConditionalSystem( MatrixType& rLeftHandSideMatrix,
								VectorType& rRightHandSideVector,
								ProcessInfo& rCurrentProcessInfo,
								Flags& rCalculationFlags )
{
    KRATOS_TRY

    GeneralVariables Variables;

    for ( unsigned int PointNumber = 0; PointNumber < 1 ; PointNumber++ )
    {
	this->CalculateKinematics(Variables, rCurrentProcessInfo, PointNumber);

        double IntegrationWeight = 1 ;  //all components are multiplied by this
        IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

        if(Variables.Options.Is(ACTIVE))
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

//************************************************************************************
//************************************************************************************

void ThermalContactDomainCondition::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
{
  //contributions to stiffness matrix calculated on the reference config
  this->CalculateAndAddThermalKm( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

  //KRATOS_WATCH(rLeftHandSideMatrix)
}


//************************************************************************************
//************************************************************************************

void ThermalContactDomainCondition::CalculateAndAddRHS(VectorType& rRightHandSideVector, GeneralVariables& rVariables, double& rIntegrationWeight)
{
  //contribution to contact forces
  this->CalculateAndAddThermalContactForces(rRightHandSideVector, rVariables, rIntegrationWeight);

  //KRATOS_WATCH(rRightHandSideVector)
}

//***********************************************************************************
//************************************************************************************

double& ThermalContactDomainCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    return rIntegrationWeight;
}



//************************************************************************************
//************************************************************************************

inline void ThermalContactDomainCondition::CalculateAndAddThermalContactForces(VectorType& rRightHandSideVector,
									       GeneralVariables& rVariables,
									       double& rIntegrationWeight)
{
    KRATOS_TRY

    //contributions to stiffness matrix calculated on the reference config
    unsigned int size=rRightHandSideVector.size();

    Vector ThermalConductionForce = ZeroVector(3);
    Vector ThermalFrictionForce   = ZeroVector(3);

    for (unsigned int ndi=0; ndi<size; ndi++)
    {
        this->CalculateThermalFrictionForce(ThermalFrictionForce[ndi], rVariables, ndi);

	this->CalculateThermalConductionForce(ThermalConductionForce[ndi], rVariables, ndi);


	rRightHandSideVector[ndi] -=(ThermalConductionForce[ndi] + ThermalFrictionForce[ndi]);
    }

    rRightHandSideVector  *=  rIntegrationWeight;

    // std::cout<<std::endl;
    // std::cout<<" ThermalConductionForce "<<ThermalConductionForce*rIntegrationWeight<<" tauk "<<mContactVariables.StabilizationFactor<<" gap Th "<<rVariables.ThermalGap<<" Projections Vector "<< rVariables.ProjectionsVector<<std::endl;
    // std::cout<<" ThermalFrictionForce "<<ThermalFrictionForce*rIntegrationWeight<<std::endl;
    // std::cout<<" Ftherm_contact "<<rRightHandSideVector<<std::endl;


    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************


void ThermalContactDomainCondition::CalculateAndAddThermalKm(MatrixType& rLeftHandSideMatrix,
							     GeneralVariables& rVariables,
							     double& rIntegrationWeight)

{    KRATOS_TRY

    //contributions to stiffness matrix calculated on the reference config
    unsigned int size=rLeftHandSideMatrix.size1();
    Vector ThermalProjection = ZeroVector(3);
    for (unsigned int ndi=0; ndi<size; ndi++)
      {
	ThermalProjection[ndi] =  rVariables.ProjectionsVector[ndi];
      }

    rLeftHandSideMatrix  =  outer_prod(ThermalProjection,ThermalProjection);
    rLeftHandSideMatrix *=  mContactVariables.StabilizationFactor;
    rLeftHandSideMatrix *=  rIntegrationWeight;

    // std::cout<<std::endl;
    // std::cout<<" Kthermcontact "<<rLeftHandSideMatrix<<std::endl;

    KRATOS_CATCH( "" )
}




//************************************************************************************
//************************************************************************************

void ThermalContactDomainCondition::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY


    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void ThermalContactDomainCondition::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void ThermalContactDomainCondition::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY


    KRATOS_CATCH( "" )
}




//************************************************************************************
//************************************************************************************

void ThermalContactDomainCondition::CalculateRelativeVelocity(GeneralVariables& rVariables, PointType & TangentVelocity, ProcessInfo& rCurrentProcessInfo)
{
    //if current tangent is not previously computed, do it here.
    rVariables.CurrentSurface.Tangent = this->CalculateCurrentTangent( rVariables.CurrentSurface.Tangent );

    if(double(inner_prod(rVariables.CurrentSurface.Tangent,rVariables.ReferenceSurface.Tangent))<0) //to give the correct direction
        rVariables.CurrentSurface.Tangent*=-1;


    // (Tangent vector previously computed)
    const int number_of_nodes = GetGeometry().size();

    //compute relative velocities
    int slave=mContactVariables.slaves[0];
    PointType  CurrentVelocity;
    for (int i = 0; i < number_of_nodes; i++ )
    {
        //Current velocity
        CurrentVelocity  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        if(i!=slave)
	    CurrentVelocity *=(-1)*(1.0/double(number_of_nodes - 1));

        TangentVelocity+=CurrentVelocity;
    }

    //Relative tangent movement of the slave if the master is fixed (the direction is implicit in the method)
    TangentVelocity =  rVariables.CurrentSurface.Tangent*(inner_prod(TangentVelocity,rVariables.CurrentSurface.Tangent));

    //Filter for low velocities (relatives to dynamic waves)
    CurrentVelocity.clear();
    CalculateRelativeDisplacement(rVariables, CurrentVelocity,rCurrentProcessInfo);

    if( norm_2(TangentVelocity)>0 ){

      if(norm_2(TangentVelocity)<1e-2*(norm_2(CurrentVelocity)/norm_2(TangentVelocity)))
	{
	  TangentVelocity.clear();
	}
    }
    else{

      TangentVelocity = CurrentVelocity;
    }
}

//************************************************************************************
//************************************************************************************


void ThermalContactDomainCondition::CalculateRelativeDisplacement(GeneralVariables& rVariables, PointType & TangentDisplacement, ProcessInfo& rCurrentProcessInfo)
{

    // (Tangent vector previously computed)
    const int number_of_nodes = GetGeometry().size();

    //compute relative displacements
    int slave=mContactVariables.slaves[0];
    PointType CurrentDisplacement;
    for (int i = 0; i < number_of_nodes; i++ )
    {
        //Displacement from the reference to the current configuration
        CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        if(i!=slave)
	    CurrentDisplacement *=(-1)*(1.0/double(number_of_nodes - 1));

        TangentDisplacement+=CurrentDisplacement;
    }

    //Relative tangent movement of the slave if the master is fixed (the direction is implicit in the method)
    TangentDisplacement = rVariables.CurrentSurface.Tangent*(inner_prod(TangentDisplacement,rVariables.CurrentSurface.Tangent));

    TangentDisplacement /= rCurrentProcessInfo[DELTA_TIME];
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
int  ThermalContactDomainCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH( "" );
}


//Note: in the restart the contact mesh is generated from the begining

void ThermalContactDomainCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    // int IntMethod = int(mThisIntegrationMethod);
    // rSerializer.save("IntegrationMethod",IntMethod);
    // rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    // rSerializer.save("ContactVariables",mContactVariables);
}

void ThermalContactDomainCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    // int IntMethod;
    // rSerializer.load("IntegrationMethod",IntMethod);
    // mThisIntegrationMethod = IntegrationMethod(IntMethod);
    // rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    // rSerializer.load("ContactVariables",mContactVariables);
}



} // Namespace Kratos
