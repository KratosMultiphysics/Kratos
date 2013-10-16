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
#include "includes/define.h"
#include "includes/element.h"
#include "custom_conditions/contact_domain_condition.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "pfem_solid_mechanics_application.h"


//#include <omp.h>

namespace Kratos
{

  KRATOS_CREATE_LOCAL_FLAG( ContactDomainCondition, PENALTY,                          1 );
  KRATOS_CREATE_LOCAL_FLAG( ContactDomainCondition, COMPUTE_FRICTION_FORCES,          2 );
  KRATOS_CREATE_LOCAL_FLAG( ContactDomainCondition, COMPUTE_FRICTION_STIFFNESS,       3 );


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ContactDomainCondition::ContactDomainCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ContactDomainCondition::ContactDomainCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ContactDomainCondition::ContactDomainCondition( ContactDomainCondition const& rOther)
    :Condition(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
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

    for(unsigned int i=0; i<<mConstitutiveLawVector.size(); i++)
    {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
    }

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
    //Element::NodeType&    MasterNode   = GetValue(MASTER_NODES).back();
    
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
    unsigned int dim2 = (number_of_nodes + 1)  * dimension;

    if ( rResult.size() != dim2 )
        rResult.resize( dim2, false );

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
    //Element::NodeType&    MasterNode   = GetValue(MASTER_NODES).back();

    rResult[index]   = MasterNode.GetDof( DISPLACEMENT_X ).EquationId();
    rResult[index+1] = MasterNode.GetDof( DISPLACEMENT_Y ).EquationId();
    if ( dimension == 3 )
        rResult[index+2] = MasterNode.GetDof( DISPLACEMENT_Z ).EquationId();

    // std::cout<<" ID "<<this->Id()<<std::endl;
    // for(unsigned int i=0; i<rResult.size(); i++)
    // 	std::cout<<" c2d Equation Id "<<rResult[i]<<std::endl;

}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void ContactDomainCondition::GetValuesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = (number_of_nodes + 1) * dimension;

    if ( rValues.size() != MatSize ) rValues.resize( MatSize, false );

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
    //NodeType& MasterNode=GetValue(MASTER_NODES).back();
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
    unsigned int MatSize = (number_of_nodes + 1) * dimension;

    if ( rValues.size() != MatSize ) rValues.resize( MatSize, false );

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
    //NodeType& MasterNode=GetValue(MASTER_NODES).back();
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
    unsigned int MatSize = (number_of_nodes + 1) * dimension;

    if ( rValues.size() != MatSize ) rValues.resize( MatSize, false );

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
    //NodeType& MasterNode=GetValue(MASTER_NODES).back();

    rValues[index] = MasterNode.GetSolutionStepValue( ACCELERATION_X, Step );
    rValues[index+1] = MasterNode.GetSolutionStepValue( ACCELERATION_Y, Step );

    if ( dimension == 3 )
        rValues[index+1] = MasterNode.GetSolutionStepValue( ACCELERATION_Z, Step );


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

 
    KRATOS_CATCH( "" )
}



////************************************************************************************
////************************************************************************************

void ContactDomainCondition::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
    
    ClearNodalForces();
   
    //Set Master Element Geometry
    SetMasterGeometry();

    unsigned int vsize=GetValue(MASTER_ELEMENTS).size();
    Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS)[vsize-1];
    
    MasterElement.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW_POINTER,mConstitutiveLawVector,CurrentProcessInfo);

    //Clear possible residual forces from the mesh refining and interpolation:
    //------------------------------------//
    const unsigned int number_of_nodes = MasterElement.GetGeometry().PointsNumber();
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
	
	VectorType & ContactForceNormal  = MasterElement.GetGeometry()[i].FastGetSolutionStepValue(FORCE_CONTACT_NORMAL);

	VectorType & ContactForceTangent  = MasterElement.GetGeometry()[i].FastGetSolutionStepValue(FORCE_CONTACT_TANGENT);

	//KRATOS_WATCH(ContactForce)
	ContactForceNormal.clear();
	ContactForceTangent.clear();
    }
    //------------------------------------//


    mVariables.Contact.IterationCounter = 0;

    //Calculate Tau Stab
    CalculateTauStab( CurrentProcessInfo );

    //Previous Gap Calculation
    CalcPreviousGap();

}

////************************************************************************************
////************************************************************************************

void ContactDomainCondition::InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
{
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
	
	VectorType & ContactForceNormal  = GetGeometry()[i].FastGetSolutionStepValue(FORCE_CONTACT_NORMAL);
	ContactForceNormal.clear();

	VectorType & ContactForceTangent  = GetGeometry()[i].FastGetSolutionStepValue(FORCE_CONTACT_TANGENT);
	ContactForceTangent.clear();
    }


    KRATOS_CATCH( "" )
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//**********************************COMPUTE TAU STAB**********************************
//************************************************************************************


void ContactDomainCondition::CalculateTauStab( ProcessInfo& rCurrentProcessInfo )
{
    //Initilialize Tau for the stabilization
    double alpha_stab=0.1;
    alpha_stab=GetProperties()[TAU_STAB];

    unsigned int vsize=GetValue(MASTER_ELEMENTS).size();
    Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS)[vsize-1];
  
    //Look at the nodes, get the slave and get the Emin

    //Contact face segment node1-node2
    unsigned int slave=mVariables.slaves.back();

    // std::cout<<" slave "<<slave<<std::endl;
    
    // std::cout<<" vertices for the contact element "<<std::endl;
    // std::cout<<" (1): ["<<GetGeometry()[0].Id()<<"] "<<GetGeometry()[0]<<std::endl;
    // std::cout<<" (2): ["<<GetGeometry()[1].Id()<<"] "<<GetGeometry()[1]<<std::endl;
    // std::cout<<" (3): ["<<GetGeometry()[2].Id()<<"] "<<GetGeometry()[2]<<std::endl;

    // std::cout<<" vertices for the contact element "<<std::endl;
    // std::cout<<" (1): ["<<GetGeometry()[0].Id()<<"] "<<GetGeometry()[0]<<std::endl;
    // std::cout<<" (2): ["<<GetGeometry()[1].Id()<<"] "<<GetGeometry()[1]<<std::endl;
    // std::cout<<" (3): ["<<GetGeometry()[2].Id()<<"] "<<GetGeometry()[2]<<std::endl;

    // std::cout<<" contact slave "<<std::endl;
    // GetGeometry()[slave].GetValue(NEIGHBOUR_ELEMENTS)[0].PrintInfo(std::cout);
    // std::cout<<std::endl;

    double Eslave=GetGeometry()[slave].GetValue(NEIGHBOUR_ELEMENTS)[0].GetProperties()[YOUNG_MODULUS];
    double Emin  =MasterElement.GetProperties()[YOUNG_MODULUS];

    //STANDARD OPTION
    if(Emin>Eslave)
	Emin=Eslave;

    mVariables.Contact.Tau=alpha_stab/Emin;


    //EXPERIMENTAL OPTION
    // const GeometryType::IntegrationPointsArrayType& integration_points = MasterElement.GetGeometry().IntegrationPoints( mThisIntegrationMethod );
   
    // //Get Current ConstitutiveMatrix
    // int size = integration_points.size();
    // std::vector<Matrix> ConstitutiveMatrix(size);   
    // MasterElement.CalculateOnIntegrationPoints(CONSTITUTIVE_MATRIX,ConstitutiveMatrix,rCurrentProcessInfo);

    // //Calc Norm of the Constitutive tensor:
    // double Cnorm=0;
    // for(int i=0; i<size; i++){
    //   for(int j=0; j<size; j++){
    // 	Cnorm += ConstitutiveMatrix[0](i,j)*ConstitutiveMatrix[0](i,j);
    //   }
    // }

    // Cnorm = sqrt(Cnorm)*0.5;

       
    // if(Emin>Cnorm){
    //   //std::cout<<" --Tau Stab A "<<mVariables.Contact.Tau<<" -- Tau Stab B "<<alpha_stab/Cnorm<<std::endl;
    //   mVariables.Contact.Tau=alpha_stab/Cnorm;
    // }
    

    
}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void ContactDomainCondition::CalculateKinematics(const double& rPointNumber,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType::ShapeFunctionsGradientsType& DN_De = mpMasterGeometry->ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    unsigned int dimension = mpMasterGeometry->WorkingSpaceDimension();
    Matrix J ( dimension , dimension);
    J = mpMasterGeometry->Jacobian( J, rPointNumber , mThisIntegrationMethod );

    Matrix InvJ;

    //Calculating the inverse of the jacobian and the parameters needed
    MathUtils<double>::InvertMatrix( J, InvJ, mVariables.detJ);

    //Compute cartesian derivatives
    noalias( mVariables.DN_DX ) = prod( DN_De[rPointNumber] , InvJ );   

    unsigned int vsize=GetValue(MASTER_ELEMENTS).size();
    Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS)[vsize-1];

    //Current Deformation Gradient and Stresses
    //CalculateDeformationGradient (mVariables.DN_DX, mVariables.F );

    //Get current DeformationGradient
    std::vector<Matrix> DeformationGradientVector (mConstitutiveLawVector.size());
    DeformationGradientVector[rPointNumber]=identity_matrix<double>( 2 );
    
    MasterElement.GetValueOnIntegrationPoints(DEFORMATION_GRADIENT,DeformationGradientVector,rCurrentProcessInfo);
    mVariables.F=DeformationGradientVector[rPointNumber];
    //Get Current Stress
    std::vector<Vector> StressVector (mConstitutiveLawVector.size());   
    StressVector[rPointNumber]=ZeroVector(3);
    MasterElement.GetValueOnIntegrationPoints(PK2_STRESS_VECTOR,StressVector,rCurrentProcessInfo);
    mVariables.StressVector=StressVector[rPointNumber];

    //std::cout<<" StressVector "<<mVariables.StressVector<<std::endl;
    
    //Get Current Strain
    std::vector<Matrix> StrainTensor (mConstitutiveLawVector.size());   
    StrainTensor[rPointNumber]=ZeroMatrix(3,3);
    MasterElement.GetValueOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_TENSOR,StrainTensor,rCurrentProcessInfo);
    mVariables.StrainVector=MathUtils<double>::StrainTensorToVector( StrainTensor[rPointNumber] );


    if(mVariables.Contact.Options.Is(PENALTY)){
      //Calculate Current PenaltyParameters
      CalcPenaltyParameters(rCurrentProcessInfo);
    }
    else{
      //Calculate Current LagrangeMultipliers
      CalcMultipliers(rCurrentProcessInfo);
    }


    KRATOS_CATCH( "" )
}


//*************************COMPUTE DEFORMATION GRADIENT*******************************
//************************************************************************************

void ContactDomainCondition::CalculateDeformationGradient(const Matrix& rDN_DX,
        Matrix& rF )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();


    if ( dimension == 2 )
    {
        rF=identity_matrix<double> ( 2 );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            //Displacement from the reference to the current configuration
            VectorType & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            VectorType & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            VectorType DeltaDisplacement=CurrentDisplacement-PreviousDisplacement;

            rF ( 0 , 0 ) += DeltaDisplacement[0]*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += DeltaDisplacement[0]*rDN_DX ( i , 1 );
            rF ( 1 , 0 ) += DeltaDisplacement[1]*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += DeltaDisplacement[1]*rDN_DX ( i , 1 );

        }

    }


    if ( dimension == 3 )
    {

        rF=identity_matrix<double> ( 3 );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            //Displacement from the reference to the current configuration
            VectorType & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            VectorType & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            VectorType DeltaDisplacement=CurrentDisplacement-PreviousDisplacement;

            rF ( 0 , 0 ) += DeltaDisplacement[0]*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += DeltaDisplacement[0]*rDN_DX ( i , 1 );
            rF ( 0 , 2 ) += DeltaDisplacement[0]*rDN_DX ( i , 2 );
            rF ( 1 , 0 ) += DeltaDisplacement[1]*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += DeltaDisplacement[1]*rDN_DX ( i , 1 );
            rF ( 1 , 2 ) += DeltaDisplacement[1]*rDN_DX ( i , 2 );
            rF ( 2 , 0 ) += DeltaDisplacement[0]*rDN_DX ( i , 0 );
            rF ( 2 , 1 ) += DeltaDisplacement[0]*rDN_DX ( i , 1 );
            rF ( 2 , 2 ) += DeltaDisplacement[0]*rDN_DX ( i , 2 );
        }

    }

    KRATOS_CATCH( "" )
}




//***********************COMPUTE LOCAL SYSTEM CONTRIBUTIONS***************************
//************************************************************************************


void ContactDomainCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag  = true;

    CalculateConditionalSystem( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );


}

//************************************************************************************
//************************************************************************************


void ContactDomainCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateConditionalSystem( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}


//************************************************************************************
//************************************************************************************


void ContactDomainCondition::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0);
}



//************************************************************************************
//************************************************************************************

void ContactDomainCondition::SetMasterGeometry()

{
    unsigned int vsize=GetValue(MASTER_ELEMENTS).size();
    Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS)[vsize-1];
    //Element::ElementType& MasterElement= GetValue(MASTER_ELEMENTS).back();

    vsize=GetValue(MASTER_NODES).size();
    Element::NodeType&    MasterNode   = GetValue(MASTER_NODES)[vsize-1];
    //Element::NodeType&    MasterNode   = GetValue(MASTER_NODES).back();

     int  slave=-1;
    for(unsigned int i=0; i<MasterElement.GetGeometry().PointsNumber(); i++)
    {
        if(MasterNode.Id()==MasterElement.GetGeometry()[i].Id())
        {
	    slave=i;
	}
    }


    if(slave>=0)
    {

        NodesArrayType vertex;
	mVariables.order.resize(GetGeometry().PointsNumber(),false);
	
        for(unsigned int i=0; i<GetGeometry().PointsNumber(); i++)
        {
            bool iset=false;
            for(unsigned int j=0; j<GetGeometry().PointsNumber(); j++)
            {

                if(GetGeometry()[i].Id()==MasterElement.GetGeometry()[j].Id())
                {
                    //vertex[i]=GetGeometry()[i];
		    //vertex.push_back( GetGeometry()[i] );
   		    //mVariables.nodes.push_back(i);
                    mVariables.order[i] = j;
 		    //std::cout<<" order ["<<i<<"] = "<<j<<std::endl;
                    iset=true;
		    break;
                }

            }

            if(iset==false)
            {
		//vertex[i]=GetGeometry()[i];
                //vertex[i]=MasterNode;
		//vertex.push_back( MasterNode );
         	mVariables.order[i] = slave;
		//std::cout<<" order ["<<i<<"] = "<<slave<<std::endl;
                mVariables.slaves.push_back(i);
            }
        }



	//Permute
	std::vector<unsigned int> permute (5);
	
	permute[0]=0; 
	permute[1]=1;
	permute[2]=2;
	permute[3]=0;
	permute[4]=1;
	
	//reorder counter-clock-wise
	mVariables.nodes.push_back(permute[mVariables.slaves.back()+1]);
	mVariables.nodes.push_back(permute[mVariables.slaves.back()+2]);

	mpMasterGeometry= MasterElement.pGetGeometry();
       
    }
    else
    {
        KRATOS_ERROR( std::invalid_argument, "MASTERNODE do not belongs to MASTER ELEMENT", "" );

    }

    
}

//************************************************************************************
//************************************************************************************


void ContactDomainCondition::CalcPreviousGap() //prediction of the lagrange multiplier
{
 
    //Contact face segment node1-node2
    unsigned int node1=mVariables.nodes[0];
    unsigned int node2=mVariables.nodes[1];
    unsigned int slave=mVariables.slaves.back();

    //Get Reference Normal
    mVariables.Contact.ReferenceSurface.Normal=GetValue(NORMAL);

    //std::cout<<" Got Normal ["<<this->Id()<<"] "<<mVariables.Contact.ReferenceSurface.Normal<<std::endl;

    //Set Reference Tangent
    mVariables.Contact.ReferenceSurface.Tangent=ComputeFaceTangent(mVariables.Contact.ReferenceSurface.Tangent,mVariables.Contact.ReferenceSurface.Normal);

    //4.- Compute Effective Gaps: (g^eff=g_n3+2*Tau*tn=2*Tau*NormalMultiplier)

    //a.- Recover tensils from previous step: (it must be done once for each time step)

    //compare to auxiliar variables stored in the contact nodes to restore the LocalTensils
    //from the previous configuration

    // unsigned int vsize=GetValue(MASTER_NODES).size();
    // Element::NodeType&    MasterNode   = GetValue(MASTER_NODES)[vsize-1];
    // Element::NodeType&    MasterNode   = GetValue(MASTER_NODES).back();

    Condition::Pointer MasterCondition = GetValue(MASTER_CONDITION);

    // std::cout<<" Master Condition : "<<MasterCondition->Id()<<" Contact "<<std::endl;

    //Get previous mechanics stored in the master node
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Matrix StressMatrix =zero_matrix<double> (dimension);
    Matrix F =zero_matrix<double> (dimension);

    Vector StressVector;
    StressVector = MasterCondition->GetValue(CAUCHY_STRESS_VECTOR);  //it means that has been stored
    F            = MasterCondition->GetValue(DEFORMATION_GRADIENT);  //it means that has been stored

    // std::cout<<" StressVector "<<StressVector<<std::endl;
    // std::cout<<" F "<<F<<std::endl;

    StressMatrix = MathUtils<double>::StressVectorToTensor( StressVector );

    //we are going to need F here from Cn-1 to Cn
    // F0 from C0 to Cn is need for the stress recovery on domain elements

    double detF =MathUtils<double>::Det(F);

    //b.- Compute the 1srt Piola Kirchhoff stress tensor  (P=J*CauchyStress*F^-T)
    mConstitutiveLawVector[0]->TransformStresses(StressMatrix,F,detF,ConstitutiveLaw::StressMeasure_Cauchy,ConstitutiveLaw::StressMeasure_PK1);

    //Compute the tension (or traction) vector T=P*N (in the Reference configuration)

    //c.- Transform to 3 components
    Matrix StressMatrix3D= zero_matrix<double> ( 3 );
    for(unsigned int i=0; i<2; i++)
    {
	for(unsigned int j=0; j<2; j++)
	{
	    StressMatrix3D(i,j)=StressMatrix(i,j);
	}
    }



    //c.- Compute (n-1) normals, tangents and relative displacements from historic mX on boundaries

    //Previous normal and tangent:  n_n-1,t_n-1
    //Previous Position
    VectorType PS  =  GetGeometry()[slave].Coordinates() - (GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,2) );
    VectorType P1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) -GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,2) );
    VectorType P2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,2) );

    //Compute Previous Normal
    mVariables.Contact.PreStepSurface.Normal=ComputeFaceNormal(mVariables.Contact.PreStepSurface.Normal,P1,P2);

    //std::cout<<" Pre Normal ["<<this->Id()<<"] "<<mVariables.Contact.PreStepSurface.Normal<<std::endl;

    if((inner_prod(mVariables.Contact.PreStepSurface.Normal,mVariables.Contact.ReferenceSurface.Normal))<0) //to give the correct direction
        mVariables.Contact.PreStepSurface.Normal*=-1;

    mVariables.Contact.PreStepSurface.Normal /= norm_2(mVariables.Contact.PreStepSurface.Normal);       //to be unitary

    if(!(norm_2(mVariables.Contact.PreStepSurface.Normal)))
    {
        mVariables.Contact.PreStepSurface.Normal=mVariables.Contact.ReferenceSurface.Normal;
    }

    //Set Previous Tangent
    mVariables.Contact.PreStepSurface.Tangent=ComputeFaceTangent(mVariables.Contact.PreStepSurface.Tangent,mVariables.Contact.PreStepSurface.Normal);
  
    //Traction vector
    mVariables.Contact.TractionVector=prod(StressMatrix3D,mVariables.Contact.PreStepSurface.Normal);



    //Reference normal: n_n,t_n  -> mVariables.Contact.ReferenceSurface.Normal / mVariables.Contact.ReferenceSurface.Tangent

    //d.- Compute A_n-1,B_n-1,L_n-1

    //A_n-1, B_n-1, L_n-1:
    BaseLengths PreviousBase;
    CalcBaseDistances(PreviousBase,P1,P2,PS,mVariables.Contact.PreStepSurface.Normal);

    //std::cout<<" L :"<<PreviousBase.L<<" A :"<<PreviousBase.A<<" B :"<<PreviousBase.B<<std::endl;

    //complete the computation of the stabilization gap
    mVariables.Contact.Tau*=PreviousBase.L;

    //std::cout<<" Tau "<<mVariables.Contact.Tau<<std::endl;

    //e.-obtain the (g_N)3 and (g_T)3 for the n-1 configuration

    VectorType DS  =  GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,2);
    VectorType D1  =  GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,2);
    VectorType D2  =  GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,2);


    mVariables.Contact.PreStepNormalGap = 0;
    mVariables.Contact.PreStepTangentialGap = 0;


    mVariables.Contact.PreStepNormalGap = inner_prod((PS-P1),mVariables.Contact.PreStepSurface.Normal);
    mVariables.Contact.PreStepTangentialGap = mVariables.Contact.PreStepNormalGap;

    mVariables.Contact.PreStepNormalGap*= inner_prod(mVariables.Contact.ReferenceSurface.Normal,mVariables.Contact.PreStepSurface.Normal);

    mVariables.Contact.PreStepNormalGap+=inner_prod(mVariables.Contact.ReferenceSurface.Normal,(D1*(-PreviousBase.A/PreviousBase.L)));
    mVariables.Contact.PreStepNormalGap+=inner_prod(mVariables.Contact.ReferenceSurface.Normal,(D2*(-PreviousBase.B/PreviousBase.L)));
    mVariables.Contact.PreStepNormalGap+=inner_prod(mVariables.Contact.ReferenceSurface.Normal,DS);


    mVariables.Contact.PreStepTangentialGap*= inner_prod(mVariables.Contact.ReferenceSurface.Tangent,mVariables.Contact.PreStepSurface.Normal);

    mVariables.Contact.PreStepTangentialGap+=inner_prod(mVariables.Contact.ReferenceSurface.Tangent,(D1*(-PreviousBase.A/PreviousBase.L)));
    mVariables.Contact.PreStepTangentialGap+=inner_prod(mVariables.Contact.ReferenceSurface.Tangent,(D2*(-PreviousBase.B/PreviousBase.L)));
    mVariables.Contact.PreStepTangentialGap+=inner_prod(mVariables.Contact.ReferenceSurface.Tangent,DS);


    //d_n-1=X_n - X_n-1

    //f.- get total effective gap as: gap_n^eff=gap_n+(PreviousTimeStep/CurrentTimeStep)*(gap_n-gap_n-1)

    //gap_n-1 (in function of the n-1 position of hte other node) gap_n-1=(g_N)3_n-1+2*Tau*tn_n-1

    double NormalTensil=0,TangentTensil=0;
    
    
    //g.- Compute normal component of the tension vector:   (tn=n·P·N)
    NormalTensil=inner_prod(mVariables.Contact.PreStepSurface.Normal,mVariables.Contact.TractionVector);

    //h.- Compute tangent component of the tension vector:  (tt=t·P·N)
    TangentTensil=inner_prod(mVariables.Contact.PreStepSurface.Tangent,mVariables.Contact.TractionVector);


    mVariables.Contact.PreStepNormalGap +=2*mVariables.Contact.Tau*NormalTensil;
    mVariables.Contact.PreStepTangentialGap +=2*mVariables.Contact.Tau*TangentTensil;
   
    //restore the computation of the stabilization gap
    mVariables.Contact.Tau/=PreviousBase.L;

    // std::cout<<"ConditionID:  "<<this->Id()<<" -> Previous Tractions [tN:"<<NormalTensil<<", tT:"<<TangentTensil<<"] "<<std::endl; 
    // std::cout<<" Previous Gaps [gN:"<<mVariables.Contact.PreStepNormalGap<<", gT:"<<mVariables.Contact.PreStepTangentialGap<<"] "<<std::endl; 
}

//************************************************************************************
//************************************************************************************



inline void ContactDomainCondition:: CalcBaseDistances (BaseLengths& Base,VectorType& P1,VectorType& P2,VectorType& PS,VectorType& Normal)
{

    Base.L=norm_2(P2-P1);


    //if the normal points from the side to the slave node:
    VectorType Projection= PS-P1;
    Projection-=Normal*(inner_prod(Projection,Normal));
    Projection+=P1;

    double sign=1;
    VectorType Pro1 = Projection-P1; 
    VectorType Pro2 = P2-P1;

    if(double(inner_prod(Pro2,Pro1))<0)
        sign*=(-1);

    //signed distance to node 1
    Base.B= sign*norm_2(P1-Projection);

    sign=1;
    Pro1 = Projection-P2;
    Pro2 = P1-P2;
    if(inner_prod(Pro2,Pro1)<0)
        sign*=(-1);

    //signed distance to node 2
    Base.A= sign*norm_2(P2-Projection);

}


//************************************************************************************
//************************************************************************************
inline VectorType & ContactDomainCondition::ComputeFaceNormal(VectorType &Normal, VectorType& P1, VectorType &P2)
{

    Normal.clear();
    Normal[0] =    P2[1] - P1[1];
    Normal[1] = - (P2[0] - P1[0]);
    Normal[2] =    0.00;

    if(norm_2(Normal)!=0)
	Normal/=norm_2(Normal);

    return Normal;
}

//************************************************************************************
//************************************************************************************


inline VectorType & ContactDomainCondition::ComputeFaceTangent(VectorType &Tangent ,VectorType& P1, VectorType &P2)
{

    Tangent.clear();
    Tangent[0] =    (P2[0] - P1[0]);
    Tangent[1] =   -(P2[1] - P1[1]);
    Tangent[2] =    0.00;

    if(norm_2(Tangent)!=0)
	Tangent/=norm_2(Tangent);

    return Tangent;

}



//************************************************************************************
//************************************************************************************


inline VectorType & ContactDomainCondition::CalcCurrentTangent ( VectorType &Tangent )
{

  unsigned int node1=mVariables.nodes[0];
  unsigned int node2=mVariables.nodes[1];

  VectorType P1  =  GetGeometry()[node1].Coordinates() + ( GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
  VectorType P2  =  GetGeometry()[node2].Coordinates() + ( GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );
  
  //Set Reference Tangent
  Tangent=ComputeFaceTangent(Tangent,P1,P2);
  
  return Tangent;

}

//************************************************************************************
//************************************************************************************


inline VectorType & ContactDomainCondition::ComputeFaceTangent(VectorType &Tangent ,VectorType& Normal)
{

    Tangent.clear();
    Tangent[0] =  - Normal[1];
    Tangent[1] =    Normal[0];
    Tangent[2] =    0.00;

    if(norm_2(Tangent)!=0)
	Tangent/=norm_2(Tangent);

    return Tangent;

}

//************************************************************************************
//************************************************************************************

inline double ContactDomainCondition::CalculateVol(const double x0, const double y0,
						     const double x1, const double y1,
						     const double x2, const double y2)
{
  return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
}

//************************************************************************************
//************************************************************************************

inline bool ContactDomainCondition::CalculatePosition(const double x0, const double y0,
							const double x1, const double y1,
							const double x2, const double y2,
							const double xc, const double yc)
{
  double area = CalculateVol(x0,y0,x1,y1,x2,y2);

  //std::cout<<" Area "<<area<<std::endl;
	    
  if(area < 1e-15)
    {
      //KRATOS_ERROR(std::logic_error,"element with zero area found","");
      std::cout<<"element with zero area found: "<<area<<" position ("<<x0<<", "<<y0<<") ("<<x1<<", "<<y1<<") ("<<x2<<", "<<y2<<") "<<std::endl;
    }

  VectorType N;

  N[0] = CalculateVol(x1,y1,x2,y2,xc,yc)  / area;
  N[1] = CalculateVol(x2,y2,x0,y0,xc,yc)  / area;
  N[2] = CalculateVol(x0,y0,x1,y1,xc,yc)  / area;

  double tol = 1e-3;
  double upper_limit = 1.0+tol;
  double lower_limit = -tol;

  if(N[0] >= lower_limit && N[1] >= lower_limit && N[2] >= lower_limit && N[0] <= upper_limit && N[1] <= upper_limit && N[2] <= upper_limit) //if the xc yc is inside the triangle
    return true;

  return false;
}

//************************************************************************************
//************************************************************************************

inline bool ContactDomainCondition::CalculateObtuseAngle(const double x0, const double y0,
							   const double x1, const double y1,
							   const double xc, const double yc)
{

  double side0 = sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) ); //master side
  double side1 = sqrt( (x0-xc)*(x0-xc) + (y0-yc)*(y0-yc) );
  double side2 = sqrt( (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1) );

  double cos_angle = 0; 
  double aux       = (2*side1*side0);
  if(aux!=0)
    cos_angle = ((side1*side1) + (side0*side0) - (side2*side2)) / aux;

  if(cos_angle<(-0.1))
    return true;

  aux       = (2*side2*side0);
  if(aux!=0)
    cos_angle = ((side2*side2) + (side0*side0) - (side1*side1)) / aux;

  if(cos_angle<(-0.1))
    return true;

  return false;
}

//************************************************************************************
//************************************************************************************


inline bool ContactDomainCondition::CheckFictiousContacts()
{

  bool real_contact = false;

  //Contact face segment node1-node2
  unsigned int node1=mVariables.nodes[0];
  unsigned int node2=mVariables.nodes[1];    
  unsigned int slave=mVariables.slaves.back();

  double offset_factor = mVariables.Contact.NormalGap; 
  
  VectorType PS  =  GetGeometry()[slave].Coordinates() + ( GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) );

  VectorType Normal =GetGeometry()[slave].FastGetSolutionStepValue(NORMAL); 
  double  Shrink             =1;//GetGeometry()[slave].FastGetSolutionStepValue(SHRINK_FACTOR);   
  VectorType Offset =GetGeometry()[slave].FastGetSolutionStepValue(OFFSET);   
  offset_factor = norm_2(Offset);

  //modify slave position projection following slave normal
  double Sx1 = PS[0]+Normal[0]*Shrink*offset_factor;
  double Sy1 = PS[1]+Normal[1]*Shrink*offset_factor;

  //modify slave position projection following master normal
  double Mx1 = PS[0]-mVariables.Contact.CurrentSurface.Normal[0]*Shrink*offset_factor;
  double My1 = PS[1]-mVariables.Contact.CurrentSurface.Normal[1]*Shrink*offset_factor;

  //Domain neighbours:
      

  //Check slave node inside the contacting domain:

  //node1:
  WeakPointerVector<Element >& rNeighbours_n1 = GetGeometry()[node1].GetValue(NEIGHBOUR_ELEMENTS);
  //node2:
  WeakPointerVector<Element >& rNeighbours_n2 = GetGeometry()[node2].GetValue(NEIGHBOUR_ELEMENTS);

  unsigned int NumberOfNeighbours_n1 = rNeighbours_n1.size();
  unsigned int NumberOfNeighbours_n2 = rNeighbours_n2.size();

  bool is_inside_a = false;
  //following slave normal projection of the slave Sx1 and Sy1
  for(unsigned int i = 0; i < NumberOfNeighbours_n1; i++)
    {
      GeometryType::PointsArrayType& vertices=rNeighbours_n1[i].GetGeometry().Points();
  
      is_inside_a = CalculatePosition( vertices[0].X(), vertices[0].Y(),
				       vertices[1].X(), vertices[1].Y(),
				       vertices[2].X(), vertices[2].Y(),
				       Sx1, Sy1);

      if(is_inside_a)
	break;
    }

  if(!is_inside_a){
  
    for(unsigned int i = 0; i < NumberOfNeighbours_n2; i++)
      {
	GeometryType::PointsArrayType& vertices=rNeighbours_n2[i].GetGeometry().Points();
      
	is_inside_a = CalculatePosition( vertices[0].X(), vertices[0].Y(),
					 vertices[1].X(), vertices[1].Y(),
					 vertices[2].X(), vertices[2].Y(),
					 Sx1, Sy1);
      
	if(is_inside_a)
	  break;
      }

  }

  bool is_inside_b = false;
  //Check projection of the slave node inside the contacting domain:
  //following master normal projection of the slave Mx1 and My1
  for(unsigned int i = 0; i < NumberOfNeighbours_n1; i++)
    {
      GeometryType::PointsArrayType& vertices=rNeighbours_n1[i].GetGeometry().Points();
  
      is_inside_b = CalculatePosition( vertices[0].X(), vertices[0].Y(),
				       vertices[1].X(), vertices[1].Y(),
				       vertices[2].X(), vertices[2].Y(),
				       Mx1, My1);

      if(is_inside_b)
	break;
    }

  if(!is_inside_b){
  

    for(unsigned int i = 0; i < NumberOfNeighbours_n2; i++)
      {
	GeometryType::PointsArrayType& vertices=rNeighbours_n2[i].GetGeometry().Points();
      
	is_inside_b = CalculatePosition( vertices[0].X(), vertices[0].Y(),
					 vertices[1].X(), vertices[1].Y(),
					 vertices[2].X(), vertices[2].Y(),
					 Mx1, My1);
      
	if(is_inside_b)
	  break;
      }

  }

  if(is_inside_a && is_inside_b)
    real_contact = true; //if the slave node is inside of the domain --> real contact
  else
    real_contact = false;

  if( real_contact == false && (is_inside_b || is_inside_a) ){ //following the master normal is in.
    std::cout<<" THERE IS a SERIOUS DOUBT IN A FICTIOUS CONTACT "<<this->Id()<<std::endl;

    VectorType P1  =  GetGeometry()[node1].Coordinates() + ( GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    VectorType P2  =  GetGeometry()[node2].Coordinates() + ( GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );

    bool is_obtuse = CalculateObtuseAngle( P1[0], P1[1],
					   P2[0], P2[1],
					   PS[0], PS[1] );

    if(!is_obtuse){
      real_contact=true;
    }
    else{
      std::cout<<" BUT IT IS OBTUSE --> FICTIOUS "<<std::endl;
    }
      
  }

  double projection = inner_prod(Normal,mVariables.Contact.CurrentSurface.Normal);

  if(real_contact==false && fabs(projection)>0.707){
    real_contact = true;
    std::cout<<" NORMALS say that this is a REAL CONTACT "<<std::endl;
  }

  if(real_contact==false){
    std::cout<<" S normal ("<<Sx1<<","<<Sy1<<")"<<std::endl;
    std::cout<<" P normal ("<<Mx1<<","<<My1<<")"<<std::endl;
    std::cout<<" Current Normal "<<mVariables.Contact.CurrentSurface.Normal<<" Slave normal "<<Normal<<std::endl;
    std::cout<<" Shrink "<<Shrink<<" offset_factor "<<offset_factor<<" Shrink*offset_factor "<<Shrink*offset_factor<<std::endl;
  }

 
  return real_contact;
}

//************************************************************************************
//************************************************************************************

void ContactDomainCondition:: CalcMultipliers(ProcessInfo& rCurrentProcessInfo)
{

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Contact face segment node1-node2
    unsigned int node1=mVariables.nodes[0];
    unsigned int node2=mVariables.nodes[1];    
    unsigned int slave=mVariables.slaves.back();

    // std::cout<<" ************ CONTACT ELEMENT "<<this->Id()<<" ************* "<<std::endl;
    // std::cout<<std::endl;
    
    // unsigned int vsize=GetValue(MASTER_ELEMENTS).size();
    // Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS)[vsize-1];
    
    // std::cout<<" master element "<<MasterElement.Id()<<std::endl;
    // std::cout<<" Elastic Modulus "<<MasterElement.GetProperties()[YOUNG_MODULUS]<<std::endl;

    // std::cout<<" Nodes 1,2,S "<<node1<<" "<<node2<<" "<<slave<<std::endl;

    // std::cout<<" Nodes ID  "<<GetGeometry()[0].Id()<<" "<<GetGeometry()[1].Id()<<" "<<GetGeometry()[2].Id()<<std::endl;

    // std::cout<<" Master Nodes ID  "<<(*mpMasterGeometry)[0].Id()<<" "<<(*mpMasterGeometry)[1].Id()<<" "<<(*mpMasterGeometry)[2].Id()<<std::endl;

    //1.- Compute tension vector:  (must be updated each iteration)
    Matrix StressMatrix ( dimension, dimension );

    // std::cout<<" 2nd PK stress "<<mVariables.StressVector<<std::endl;

    //a.- Assign initial 2nd Piola Kirchhoff stress:
    StressMatrix=MathUtils<double>::StressVectorToTensor( mVariables.StressVector );

    //b.- Compute the 1srt Piola Kirchhoff stress tensor  (P=F·S)
    StressMatrix=prod(mVariables.F,StressMatrix);

    //c.- Transform to 3 components
    Matrix StressMatrix3D= zero_matrix<double> ( 3 );
    for(unsigned int i=0; i<2; i++)
      {
	for(unsigned int j=0; j<2; j++)
	  {
	    StressMatrix3D(i,j)=StressMatrix(i,j);
	  }
      }
   
    // std::cout<<" 1st PK stress "<<StressMatrix3D<<std::endl;

    //c.- Compute the tension (or traction) vector T=P*N (in the Reference configuration)
    mVariables.Contact.TractionVector=prod(StressMatrix3D,mVariables.Contact.ReferenceSurface.Normal);

    //std::cout<<" Reference Normal "<<mVariables.Contact.ReferenceSurface.Normal<<std::endl;
    //std::cout<<" Traction  Vector "<<mVariables.Contact.TractionVector<<std::endl;


    //d.- Compute the Current Normal and Tangent

    VectorType PS  =  GetGeometry()[slave].Coordinates() + ( GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) );
    VectorType P1  =  GetGeometry()[node1].Coordinates() + ( GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    VectorType P2  =  GetGeometry()[node2].Coordinates() + ( GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );

    //compute the current normal vector
    mVariables.Contact.CurrentSurface.Normal=ComputeFaceNormal(mVariables.Contact.CurrentSurface.Normal,P1,P2);

    //std::cout<<" Current Normal "<<mVariables.Contact.CurrentSurface.Normal<<std::endl;

    if(double(inner_prod(mVariables.Contact.CurrentSurface.Normal,mVariables.Contact.ReferenceSurface.Normal))<0) //to give the correct direction
        mVariables.Contact.CurrentSurface.Normal*=-1;

    mVariables.Contact.CurrentSurface.Normal /= norm_2(mVariables.Contact.CurrentSurface.Normal);  //to be unitary

    if(!(norm_2(mVariables.Contact.CurrentSurface.Normal)))
        mVariables.Contact.CurrentSurface.Normal=mVariables.Contact.ReferenceSurface.Normal;


    //compute the current tangent vector
    //mVariables.Contact.CurrentSurface.Tangent=ComputeFaceTangent(mVariables.Contact.CurrentSurface.Tangent,P1,P2);
    mVariables.Contact.CurrentSurface.Tangent=ComputeFaceTangent(mVariables.Contact.CurrentSurface.Tangent,mVariables.Contact.CurrentSurface.Normal);



    //std::cout<<" reference face  normal  "<<mVariables.Contact.ReferenceSurface.Normal<<std::endl;
    //std::cout<<" reference face  tangent  "<<mVariables.Contact.ReferenceSurface.Tangent<<std::endl;


    //std::cout<<" current face  normal  "<<mVariables.Contact.CurrentSurface.Normal<<std::endl;
    //std::cout<<" current face  tangent  "<<mVariables.Contact.CurrentSurface.Tangent<<std::endl;

    //Current normal:   mVariables.Contact.ReferenceSurface.Normal

    //2.- Compute normal component of the tension vector:   (tn=n·P·N)
    mVariables.Contact.CurrentTensil.Normal=inner_prod(mVariables.Contact.CurrentSurface.Normal,mVariables.Contact.TractionVector);

    //3.- Compute tangent component of the tension vector:  (tt=t·P·N)
    mVariables.Contact.CurrentTensil.Tangent=inner_prod(mVariables.Contact.CurrentSurface.Tangent,mVariables.Contact.TractionVector);

    //4.- Compute Effective Gaps: (g^eff=g_n3+2*Tau*tn=2*Tau*NormalMultiplier)

    //Reference normal: n_n,t_n  -> mVariables.Contact.ReferenceSurface.Normal / mVariables.Contact.Tangent
    //Current normal:   n,t      -> mVariables.Contact.CurrentSurface.Normal /  mVariables.Contact.CurrentSurface.Tangent

    //d.- Compute A_n,B_n,L_n
    mVariables.Contact.ReferenceBase.resize(1);
    mVariables.Contact.CurrentBase.resize(1);

    //a, b, l:
    CalcBaseDistances (mVariables.Contact.CurrentBase[0],P1,P2,PS,mVariables.Contact.CurrentSurface.Normal);

    //Write Current Positions:
    // std::cout<<" Current position node 1 "<<P1<<std::endl;
    // std::cout<<" Current position node 2 "<<P2<<std::endl;
    // std::cout<<" Current position node s "<<PS<<std::endl;


    //A, B, L:

    PS =  GetGeometry()[slave].Coordinates();
    P1 =  GetGeometry()[node1].Coordinates();
    P2 =  GetGeometry()[node2].Coordinates();

    CalcBaseDistances (mVariables.Contact.ReferenceBase[0],P1,P2,PS,mVariables.Contact.ReferenceSurface.Normal);


    //complete the computation of the stabilization gap
    mVariables.Contact.Tau*=mVariables.Contact.ReferenceBase[0].L;

    //e.-obtain the (g_N)3 and (g_T)3 for the n configuration
    //Write Current Positions:
    // std::cout<<" Reference position node 1 "<<P1<<std::endl;
    // std::cout<<" Reference position node 2 "<<P2<<std::endl;
    // std::cout<<" Reference position node s "<<PS<<std::endl;

    double ReferenceGapN = inner_prod((PS - P1),mVariables.Contact.ReferenceSurface.Normal);
    //std::cout<<" Reference GAP "<<ReferenceGapN<<std::endl;
    
    double ReferenceGapT = ReferenceGapN;

    double H = ReferenceGapN;

    VectorType DS  =  GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1);
    VectorType D1  =  GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1);
    VectorType D2  =  GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1);

    //(g_N)3
    ReferenceGapN*=inner_prod(mVariables.Contact.CurrentSurface.Normal,mVariables.Contact.ReferenceSurface.Normal);

    ReferenceGapN+=inner_prod(mVariables.Contact.CurrentSurface.Normal,(D1*(-mVariables.Contact.ReferenceBase[0].A/mVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapN+=inner_prod(mVariables.Contact.CurrentSurface.Normal,(D2*(-mVariables.Contact.ReferenceBase[0].B/mVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapN+=inner_prod(mVariables.Contact.CurrentSurface.Normal,DS);

    //(g_T)3
    ReferenceGapT*=inner_prod(mVariables.Contact.CurrentSurface.Tangent,mVariables.Contact.ReferenceSurface.Normal);

    ReferenceGapT+=inner_prod(mVariables.Contact.CurrentSurface.Tangent,(D1*(-mVariables.Contact.ReferenceBase[0].A/mVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapT+=inner_prod(mVariables.Contact.CurrentSurface.Tangent,(D2*(-mVariables.Contact.ReferenceBase[0].B/mVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapT+=inner_prod(mVariables.Contact.CurrentSurface.Tangent,DS);

 
    //Write Displacements:
    // std::cout<<" displacement node 1 "<<D1<<std::endl;
    // std::cout<<" displacement node 2 "<<D2<<std::endl;
    // std::cout<<" displacement node s "<<DS<<std::endl;

    // std::cout<<" L :"<<mVariables.Contact.ReferenceBase[0].L<<" A :"<<mVariables.Contact.ReferenceBase[0].A<<" B :"<<mVariables.Contact.ReferenceBase[0].B<<std::endl;
    // std::cout<<" gN3 ref "<<ReferenceGapN<<std::endl;


     // std::cout<<" current   normal  "<<mVariables.Contact.CurrentSurface.Normal<<std::endl;
     // std::cout<<" reference normal  "<<mVariables.Contact.ReferenceSurface.Normal<<std::endl;


    mVariables.Contact.NormalGap=ReferenceGapN; //(g_N)3 -- needed in the Kcont1 computation
    mVariables.Contact.TangentialGap=ReferenceGapT; //(g_T)3 -- needed in the Kcont1 computation

    //f.- get total effective gap as: gap_n^eff=gap_n+(PreviousTimeStep/CurrentTimeStep)*(gap_n-gap_n-1)

    //gap_n   (in function of the n position of the other node) gap_n=(g_N)3+2*Tau*tn_n

    ReferenceGapN+=2*mVariables.Contact.Tau*mVariables.Contact.CurrentTensil.Normal;
    ReferenceGapT+=2*mVariables.Contact.Tau*mVariables.Contact.CurrentTensil.Tangent;

    mVariables.Contact.TangentialGapsign=1;

    if((mVariables.Contact.TangentialGap)<0)
    {
        mVariables.Contact.TangentialGapsign*=(-1);
    }


    if(H==0) mVariables.Contact.TangentialGapsign=0;
    
    //CORRECTION: to skip change on diagonals problems //convergence problems !!! take care;
    if(mVariables.Contact.NormalGap>0 && ReferenceGapN<0)
      {
	//look at the magnitud
	if(fabs(mVariables.Contact.NormalGap) > 2*fabs(ReferenceGapN)){
	  ReferenceGapN = 0; //mVariables.Contact.NormalGap +  mVariables.Contact.Tau*mVariables.Contact.CurrentTensil.Normal;
	}

      }

    //std::cout<<" Tensil "<<mVariables.Contact.CurrentTensil.ReferenceSurface.Normal<<" Tau "<<mVariables.Contact.Tau<<" product "<<2*mVariables.Contact.Tau*mVariables.Contact.CurrentTensil.ReferenceSurface.Normal<<std::endl;
    //std::cout<<" gN3 ref total "<<ReferenceGap<<std::endl;

    //5.- Compute (Lagrange) Multipliers

    //From effective gaps set active contact domain:

    double EffectiveGapN = ReferenceGapN;
    double EffectiveGapT = ReferenceGapT;

    double CurrentTimeStep  = rCurrentProcessInfo[DELTA_TIME];
    double PreviousTimeStep = rCurrentProcessInfo[PREVIOUS_DELTA_TIME];
    
    
    if(mVariables.Contact.PreStepNormalGap!=0 && mVariables.Contact.IterationCounter<1){
	EffectiveGapN+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapN-mVariables.Contact.PreStepNormalGap);
	// std::cout<<" Effective prediction first iteration +:"<<(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapN-mVariables.Contact.PreStepNormalGap)<<" PreStepNormalGap "<<mVariables.Contact.PreStepNormalGap<<" ReferenceGapN "<<ReferenceGapN<<std::endl;
	// std::cout<<" EffectiveGapN "<<EffectiveGapN<<" PreviousEffectiveGapN "<<ReferenceGapN<<std::endl;
    }

    //only in the first iteration:
    //mVariables.Contact.PreStepNormalGap=ReferenceGapN;

    if(mVariables.Contact.PreStepTangentialGap!=0 && mVariables.Contact.IterationCounter<1){
      EffectiveGapT+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapT-mVariables.Contact.PreStepTangentialGap);
      // std::cout<<" Effective prediction first iteration +:"<<(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapT-mVariables.Contact.PreStepTangentialGap)<<" PreStepTangentialGap "<<mVariables.Contact.PreStepTangentialGap<<" ReferenceGapT "<<ReferenceGapT<<std::endl;
      // 	std::cout<<" EffectiveGapT "<<EffectiveGapT<<" PreviousEffectiveGapT "<<ReferenceGapT<<std::endl;

    }

    // std::cout<<"ConditionID:  "<<this->Id()<<" -> Previous Gap [gN:"<<ReferenceGapN<<", gT:"<<ReferenceGapT<<"] "<<std::endl; 
    // std::cout<<" -> Effective Gap [gN:"<<EffectiveGapN<<", gT:"<<EffectiveGapT<<"] "<<std::endl; 


    //only in the first iteration:
    //mVariables.Contact.PreStepTangentialGap=ReferenceGapT;

    //std::cout<<" PreTime "<<Time.PreStep<<" Time "<<Time.Step<<std::endl;

    //CHECK IF THE ELEMENT IS ACTIVE:

    mVariables.Contact.Options.Set(SLIP,false);

    // if( mVariables.Contact.NormalGap>0 && EffectiveGapN<0 && mVariables.Contact.NormalGap>mVariables.Contact.ReferenceBase[0].L*0.001){
    //   EffectiveGapN =  mVariables.Contact.NormalGap;
    // }

    //CORRECTION: to skip tip contact elements problems:
    
    bool check_fictious_geometry = false;
    
    //Check ORTHOGONAL FACES in contact
    if(check_fictious_geometry ==true){
      VectorType& SlaveNormal  =  GetGeometry()[slave].FastGetSolutionStepValue(NORMAL);
      double orthogonal = inner_prod(SlaveNormal,mVariables.Contact.CurrentSurface.Normal);

      if(EffectiveGapN<=0 && fabs(orthogonal)<=1){

	bool real_contact = CheckFictiousContacts();

	if(!real_contact || fabs(orthogonal)<=0.25){
	  EffectiveGapN = 1; //not active element: geometrically wrong
	  std::cout<<" DISABLE ContactElement "<<this->Id()<<" real contact "<<real_contact<<" geometrically orthogonal "<<orthogonal<<std::endl;
	}
      }
    }
    //Check ORTHOGONAL FACES in contact
    
   
    if(EffectiveGapN<=0)   //if(EffectiveGap<0){
    {

        mVariables.Contact.Options.Set(ACTIVE);  //normal contact active

        if(fabs(EffectiveGapT)<=mVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN))
        {
	    mVariables.Contact.Options.Set(SLIP,false); //contact stick case active
	    rCurrentProcessInfo[NUMBER_OF_STICK_CONTACTS] += 1;
        }
        else
        {

	    mVariables.Contact.Options.Set(SLIP,true);  //contact slip  case active
	    rCurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS] += 1;
        }

	rCurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] += 1;
	
	this->Set(ContactDomainCondition::ACTIVE);
    }
    else
    {
	mVariables.Contact.Options.Set(ACTIVE,false); //normal contact not active
	this->Reset(ContactDomainCondition::ACTIVE);
    }


    //mVariables.Contact.Options.Set(SLIP,false); //impose stick

    //From total current gap compute multipliers:

    //mVariables.Contact.NormalMultiplier = EffectiveGap*(1./(2.0*mVariables.Contact.Tau)); //posible computation of the Lagrange Multiplier
    mVariables.Contact.NormalMultiplier =mVariables.Contact.CurrentTensil.Normal;
    mVariables.Contact.NormalMultiplier+=mVariables.Contact.NormalGap*(1./(2.0*mVariables.Contact.Tau));

    mVariables.Contact.TangentialMultiplier =mVariables.Contact.CurrentTensil.Tangent;
    mVariables.Contact.TangentialMultiplier+=mVariables.Contact.TangentialGap*(1./(2.0*mVariables.Contact.Tau));


    if(mVariables.Contact.TangentialMultiplier<0)  //add the sign of the Lagrange Multiplier
    {
        mVariables.Contact.TangentialGapsign*=(-1);
    }

    //std::cout<<" Effective GapN "<<EffectiveGapN<<" NormalMultiplier "<<mVariables.Contact.NormalMultiplier<<" CurrentTensil.N "<<mVariables.Contact.CurrentTensil.Normal<<" GapN "<<mVariables.Contact.NormalGap<<" ReferenceGapN "<<ReferenceGapN<<" Tau "<<mVariables.Contact.Tau<<" iteration "<<mVariables.Contact.IterationCounter<<std::endl;

    // std::cout<<" current face  normal  "<<mVariables.Contact.CurrentSurface.Normal<<std::endl;
    // std::cout<<" current face  tangent  "<<mVariables.Contact.CurrentSurface.Tangent<<std::endl;

    // if(mVariables.Contact.Options.Is(ACTIVE) && mVariables.Contact.NormalGap>0)
    //   std::cout<<" Condition ["<<this->Id()<<"]:  Effective GapN "<<EffectiveGapN<<" NormalMultiplier "<<mVariables.Contact.NormalMultiplier<<" CurrentTensil.N "<<mVariables.Contact.CurrentTensil.Normal<<" GapN "<<mVariables.Contact.NormalGap<<" ReferenceGapN "<<ReferenceGapN<<" Tau "<<mVariables.Contact.Tau<<" iteration "<<mVariables.Contact.IterationCounter<<std::endl;
      

    if(mVariables.Contact.IterationCounter < 1)
      mVariables.Contact.IterationCounter += 1;

    //restore computation of the stabilization gap
    mVariables.Contact.Tau/=mVariables.Contact.ReferenceBase[0].L;

}


//************************************************************************************
//************************************************************************************

void ContactDomainCondition:: CalcPenaltyParameters(ProcessInfo& rCurrentProcessInfo)
{


    //Contact face segment node1-node2
    unsigned int node1=mVariables.nodes[0];
    unsigned int node2=mVariables.nodes[1];    
    unsigned int slave=mVariables.slaves.back();

    // std::cout<<" ************ CONTACT ELEMENT "<<this->Id()<<" ************* "<<std::endl;
    // std::cout<<std::endl;
    
    // unsigned int vsize=GetValue(MASTER_ELEMENTS).size();
    // Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS)[vsize-1];
    
    // std::cout<<" master element "<<MasterElement.Id()<<std::endl;
    // std::cout<<" Elastic Modulus "<<MasterElement.GetProperties()[YOUNG_MODULUS]<<std::endl;

    // std::cout<<" Nodes 1,2,S "<<node1<<" "<<node2<<" "<<slave<<std::endl;

    // std::cout<<" Nodes ID  "<<GetGeometry()[0].Id()<<" "<<GetGeometry()[1].Id()<<" "<<GetGeometry()[2].Id()<<std::endl;

    // std::cout<<" Master Nodes ID  "<<(*mpMasterGeometry)[0].Id()<<" "<<(*mpMasterGeometry)[1].Id()<<" "<<(*mpMasterGeometry)[2].Id()<<std::endl;


    //1.- Compute the Current Normal and Tangent

    VectorType PS  =  GetGeometry()[slave].Coordinates() + ( GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) );
    VectorType P1  =  GetGeometry()[node1].Coordinates() + ( GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    VectorType P2  =  GetGeometry()[node2].Coordinates() + ( GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );

    //compute the current normal vector
    mVariables.Contact.CurrentSurface.Normal=ComputeFaceNormal(mVariables.Contact.CurrentSurface.Normal,P1,P2);

    //std::cout<<" Current Normal "<<mVariables.Contact.CurrentSurface.Normal<<std::endl;

    if(double(inner_prod(mVariables.Contact.CurrentSurface.Normal,mVariables.Contact.ReferenceSurface.Normal))<0) //to give the correct direction
        mVariables.Contact.CurrentSurface.Normal*=-1;

    mVariables.Contact.CurrentSurface.Normal /= norm_2(mVariables.Contact.CurrentSurface.Normal);  //to be unitary

    if(!(norm_2(mVariables.Contact.CurrentSurface.Normal)))
        mVariables.Contact.CurrentSurface.Normal=mVariables.Contact.ReferenceSurface.Normal;


    //compute the current tangent vector
    //mVariables.Contact.CurrentSurface.Tangent=ComputeFaceTangent(mVariables.Contact.CurrentSurface.Tangent,P1,P2);
    mVariables.Contact.CurrentSurface.Tangent=ComputeFaceTangent(mVariables.Contact.CurrentSurface.Tangent,mVariables.Contact.CurrentSurface.Normal);



    //std::cout<<" reference face  normal  "<<mVariables.Contact.ReferenceSurface.Normal<<std::endl;
    //std::cout<<" reference face  tangent  "<<mVariables.Contact.ReferenceSurface.Tangent<<std::endl;


    //std::cout<<" current face  normal  "<<mVariables.Contact.CurrentSurface.Normal<<std::endl;
    //std::cout<<" current face  tangent  "<<mVariables.Contact.CurrentSurface.Tangent<<std::endl;

    //Current normal:   mVariables.Contact.ReferenceSurface.Normal

    //Reference normal: n_n,t_n  -> mVariables.Contact.ReferenceSurface.Normal / mVariables.Contact.Tangent
    //Current normal:   n,t      -> mVariables.Contact.CurrentSurface.Normal /  mVariables.Contact.CurrentSurface.Tangent

    //d.- Compute A_n,B_n,L_n
    mVariables.Contact.ReferenceBase.resize(1);
    mVariables.Contact.CurrentBase.resize(1);

    //a, b, l:
    CalcBaseDistances (mVariables.Contact.CurrentBase[0],P1,P2,PS,mVariables.Contact.CurrentSurface.Normal);

    //Write Current Positions:
    // std::cout<<" Current position node 1 "<<P1<<std::endl;
    // std::cout<<" Current position node 2 "<<P2<<std::endl;
    // std::cout<<" Current position node s "<<PS<<std::endl;


    //A, B, L:

    PS =  GetGeometry()[slave].Coordinates();
    P1 =  GetGeometry()[node1].Coordinates();
    P2 =  GetGeometry()[node2].Coordinates();

    CalcBaseDistances (mVariables.Contact.ReferenceBase[0],P1,P2,PS,mVariables.Contact.ReferenceSurface.Normal);


    //complete the computation of the stabilization gap
    mVariables.Contact.Tau*=mVariables.Contact.ReferenceBase[0].L;

    //e.-obtain the (g_N)3 and (g_T)3 for the n configuration
    //Write Current Positions:
    // std::cout<<" Reference position node 1 "<<P1<<std::endl;
    // std::cout<<" Reference position node 2 "<<P2<<std::endl;
    // std::cout<<" Reference position node s "<<PS<<std::endl;

    double ReferenceGapN = inner_prod((PS - P1),mVariables.Contact.ReferenceSurface.Normal);
    //std::cout<<" Reference GAP "<<ReferenceGapN<<std::endl;
    
    double ReferenceGapT = ReferenceGapN;

    double H = ReferenceGapN;

    VectorType DS  =  GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1);
    VectorType D1  =  GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1);
    VectorType D2  =  GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1);

    //(g_N)3
    ReferenceGapN*=inner_prod(mVariables.Contact.CurrentSurface.Normal,mVariables.Contact.ReferenceSurface.Normal);

    ReferenceGapN+=inner_prod(mVariables.Contact.CurrentSurface.Normal,(D1*(-mVariables.Contact.ReferenceBase[0].A/mVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapN+=inner_prod(mVariables.Contact.CurrentSurface.Normal,(D2*(-mVariables.Contact.ReferenceBase[0].B/mVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapN+=inner_prod(mVariables.Contact.CurrentSurface.Normal,DS);

    //(g_T)3
    ReferenceGapT*=inner_prod(mVariables.Contact.CurrentSurface.Tangent,mVariables.Contact.ReferenceSurface.Normal);

    ReferenceGapT+=inner_prod(mVariables.Contact.CurrentSurface.Tangent,(D1*(-mVariables.Contact.ReferenceBase[0].A/mVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapT+=inner_prod(mVariables.Contact.CurrentSurface.Tangent,(D2*(-mVariables.Contact.ReferenceBase[0].B/mVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapT+=inner_prod(mVariables.Contact.CurrentSurface.Tangent,DS);

 
    //Write Displacements:
    // std::cout<<" displacement node 1 "<<D1<<std::endl;
    // std::cout<<" displacement node 2 "<<D2<<std::endl;
    // std::cout<<" displacement node s "<<DS<<std::endl;

    // std::cout<<" L :"<<mVariables.Contact.ReferenceBase[0].L<<" A :"<<mVariables.Contact.ReferenceBase[0].A<<" B :"<<mVariables.Contact.ReferenceBase[0].B<<std::endl;
    // std::cout<<" gN3 ref "<<ReferenceGapN<<std::endl;


     // std::cout<<" current   normal  "<<mVariables.Contact.CurrentSurface.Normal<<std::endl;
     // std::cout<<" reference normal  "<<mVariables.Contact.ReferenceSurface.Normal<<std::endl;


    mVariables.Contact.NormalGap=ReferenceGapN; //(g_N)3 -- needed in the Kcont1 computation
    mVariables.Contact.TangentialGap=ReferenceGapT; //(g_T)3 -- needed in the Kcont1 computation

    //f.- get total effective gap as: gap_n^eff=gap_n+(PreviousTimeStep/CurrentTimeStep)*(gap_n-gap_n-1)

    //gap_n   (in function of the n position of the other node)

    mVariables.Contact.TangentialGapsign=1;

    if((mVariables.Contact.TangentialGap)<0)
    {
        mVariables.Contact.TangentialGapsign*=(-1);
    }

    if(H==0) mVariables.Contact.TangentialGapsign=0;

    //5.- Compute (Lagrange) Multipliers


    // double CurrentTimeStep  = rCurrentProcessInfo[DELTA_TIME];
    // double PreviousTimeStep = rCurrentProcessInfo[PREVIOUS_DELTA_TIME];
        

    mVariables.Contact.NormalMultiplier =  H * (0.5/mVariables.Contact.Tau); //(1/(2*Tau)) is now KN, the penalty parameter
    mVariables.Contact.TangentialMultiplier = 0;
    //std::cout<<" PreTime "<<Time.PreStep<<" Time "<<Time.Step<<std::endl;

    //CHECK IF THE ELEMENT IS ACTIVE:

    mVariables.Contact.Options.Set(SLIP,false); 

    if(mVariables.Contact.NormalMultiplier<0)   //if(EffectiveGap<0){
    {

	 mVariables.Contact.Options.Set(ACTIVE,true); //normal contact active

        // if(fabs(EffectiveGapT)<=mVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN))
        // {
        //     mVariables.Contact.Options.Set(SLIP,false);  //contact stick case active
	//     rCurrentProcessInfo[NUMBER_OF_STICK_CONTACTS] += 1;
        // }
        // else
        // {
        //     mVariables.Contact.Options.Set(SLIP,true);  //contact slip  case active
	//     rCurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS] += 1;
        // }
	mVariables.Contact.Options.Set(COMPUTE_FRICTION_STIFFNESS); 
	mVariables.Contact.Options.Set(COMPUTE_FRICTION_FORCES);
	//rCurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] += 1;

    }
    else
    {
	mVariables.Contact.Options.Set(ACTIVE,false); //normal contact not active
    }


    mVariables.Contact.Options.Set(SLIP,false); //impose stick


    if(mVariables.Contact.TangentialMultiplier<0)  //add the sign of the Lagrange Multiplier
    {
        mVariables.Contact.TangentialGapsign*=(-1);
    }

    if(mVariables.Contact.IterationCounter < 1)
      mVariables.Contact.IterationCounter += 1;

    
    //now NormalMultiplier is the penalty parameter:
    mVariables.Contact.NormalMultiplier = (0.5/mVariables.Contact.Tau);

    //restore computation of the stabilization gap
    mVariables.Contact.Tau/=mVariables.Contact.ReferenceBase[0].L;

}





//************************************************************************************
//************************************************************************************


void ContactDomainCondition::CalcDomainShapeN()
{

    unsigned int ndi=mVariables.nodes[0];
    unsigned int ndj=mVariables.nodes[1];
    unsigned int ndk=mVariables.slaves[0];
    unsigned int ndr=3;

    //Set discrete variations of the shape function on the normal and tangent directions:


    Matrix DN_DX = mVariables.DN_DX;
    for(unsigned int i=0;i<DN_DX.size1();i++)
      {
	for(unsigned int j=0;j<DN_DX.size2();j++)
	  {
	    mVariables.DN_DX(i,j) = DN_DX(mVariables.order[i],j);
	  }
      }

    // std::cout<<" DN_DX "<<DN_DX<<std::endl;
    // std::cout<<" Variables_DN_DX "<<mVariables.DN_DX<<std::endl;

    //dN_dn:

    mVariables.Contact.dN_dn=ZeroVector(4);

    mVariables.Contact.dN_dn[ndi]=(-1)*mVariables.Contact.CurrentBase[0].A/mVariables.Contact.CurrentBase[0].L;
    mVariables.Contact.dN_dn[ndj]=(-1)*mVariables.Contact.CurrentBase[0].B/mVariables.Contact.CurrentBase[0].L;
    mVariables.Contact.dN_dn[ndk]= 1.0;


    //dN_drn:

    mVariables.Contact.dN_drn=ZeroVector(4);

    mVariables.Contact.dN_drn[ndi]=(-1)*mVariables.Contact.ReferenceBase[0].A/mVariables.Contact.ReferenceBase[0].L;
    mVariables.Contact.dN_drn[ndj]=(-1)*mVariables.Contact.ReferenceBase[0].B/mVariables.Contact.ReferenceBase[0].L;
    mVariables.Contact.dN_drn[ndk]= 1.0;


    //dN_dt:

    mVariables.Contact.dN_dt=ZeroVector(4);

    mVariables.Contact.dN_dt[ndi]=   1.0/mVariables.Contact.CurrentBase[0].L;
    mVariables.Contact.dN_dt[ndj]=(-1.0)/mVariables.Contact.CurrentBase[0].L;

    // std::cout<<" dN_dn "<<mVariables.Contact.dN_dn<<std::endl;
    // std::cout<<" dN_drn"<<mVariables.Contact.dN_drn<<std::endl;
    // std::cout<<" dN_dt "<<mVariables.Contact.dN_dt<<std::endl;
  
    //TsigmaP :
    mVariables.Contact.Tsigma.resize(4);
    FSigmaP(mVariables.Contact.Tsigma,mVariables.Contact.CurrentSurface.Tangent,ndi,ndj,ndk,ndr);

    // for(unsigned int i=0; i<mVariables.Contact.Tsigma.size(); i++)
    //   std::cout<<" i: "<<i<<" Tsigma "<<mVariables.Contact.Tsigma[i]<<std::endl;


    //NsigmaP :
    mVariables.Contact.Nsigma.resize(4);
    FSigmaP(mVariables.Contact.Nsigma,mVariables.Contact.CurrentSurface.Normal,ndi,ndj,ndk,ndr);

    // for(unsigned int i=0; i<mVariables.Contact.Nsigma.size(); i++)
    //   std::cout<<" i: "<<i<<" Nsigma "<<mVariables.Contact.Nsigma[i]<<std::endl;


    mVariables.DN_DX = DN_DX;

}


//************************************************************************************
//************************************************************************************

void ContactDomainCondition::FSigmaP(std::vector<Vector > &SigmaP, VectorType& DirVector,unsigned int &ndi,unsigned int &ndj,unsigned int &ndk,unsigned int &ndr)
{
    //Computation with the ndi and storage to ndj

    // std::cout<<" DirVector "<<DirVector<<std::endl;
    // std::cout<<" StressVector "<<mVariables.StressVector<<std::endl;

    FSigmaPnd(SigmaP,DirVector,ndi,ndi);

    FSigmaPnd(SigmaP,DirVector,ndj,ndj);

    SigmaP[ndk]=ZeroVector(2);

    FSigmaPnd(SigmaP,DirVector,ndk,ndr);

}

//************************************************************************************
//************************************************************************************


void ContactDomainCondition::FSigmaPnd(std::vector<Vector > &SigmaP, VectorType& DirVector,unsigned int &ndi,unsigned int &ndj)
{
    //Computation with the ndi and storage to ndj
    SigmaP[ndj]=ZeroVector(2);

    //part1:
    SigmaP[ndj][0]= DirVector[0]*mVariables.Contact.ReferenceSurface.Normal[0]*(mVariables.StressVector[0]*mVariables.DN_DX(ndi,0)+mVariables.StressVector[2]*mVariables.DN_DX(ndi,1))+ DirVector[0]*mVariables.Contact.ReferenceSurface.Normal[1]*(mVariables.StressVector[1]*mVariables.DN_DX(ndi,1)+mVariables.StressVector[2]*mVariables.DN_DX(ndi,0));

       
    SigmaP[ndj][1]= DirVector[1]*mVariables.Contact.ReferenceSurface.Normal[1]*(mVariables.StressVector[1]*mVariables.DN_DX(ndi,1)+mVariables.StressVector[2]*mVariables.DN_DX(ndi,0))+ DirVector[1]*mVariables.Contact.ReferenceSurface.Normal[0]*(mVariables.StressVector[0]*mVariables.DN_DX(ndi,0)+mVariables.StressVector[2]*mVariables.DN_DX(ndi,1));

    //part2:
    std::vector<Vector> FD(4);

    FD[0].resize(4);
    FD[0][0]=(mVariables.F(0,0)*mVariables.ConstitutiveMatrix(0,0)+mVariables.F(0,1)*mVariables.ConstitutiveMatrix(2,0));
    FD[0][1]=(mVariables.F(0,0)*mVariables.ConstitutiveMatrix(0,1)+mVariables.F(0,1)*mVariables.ConstitutiveMatrix(2,1));
    FD[0][2]=(mVariables.F(0,0)*mVariables.ConstitutiveMatrix(0,2)+mVariables.F(0,1)*mVariables.ConstitutiveMatrix(2,2));
    FD[0][3]=(mVariables.F(0,0)*mVariables.ConstitutiveMatrix(0,2)+mVariables.F(0,1)*mVariables.ConstitutiveMatrix(2,2));

    FD[1].resize(4);
    FD[1][0]=(mVariables.F(1,1)*mVariables.ConstitutiveMatrix(1,0)+mVariables.F(1,0)*mVariables.ConstitutiveMatrix(2,0));
    FD[1][1]=(mVariables.F(1,1)*mVariables.ConstitutiveMatrix(1,1)+mVariables.F(1,0)*mVariables.ConstitutiveMatrix(2,1));
    FD[1][2]=(mVariables.F(1,1)*mVariables.ConstitutiveMatrix(1,2)+mVariables.F(1,0)*mVariables.ConstitutiveMatrix(2,2));
    FD[1][3]=(mVariables.F(1,1)*mVariables.ConstitutiveMatrix(1,2)+mVariables.F(1,0)*mVariables.ConstitutiveMatrix(2,2));

    FD[2].resize(4);
    FD[2][0]=(mVariables.F(0,1)*mVariables.ConstitutiveMatrix(1,0)+mVariables.F(0,0)*mVariables.ConstitutiveMatrix(2,0));
    FD[2][1]=(mVariables.F(0,1)*mVariables.ConstitutiveMatrix(1,1)+mVariables.F(0,0)*mVariables.ConstitutiveMatrix(2,1));
    FD[2][2]=(mVariables.F(0,1)*mVariables.ConstitutiveMatrix(1,2)+mVariables.F(0,0)*mVariables.ConstitutiveMatrix(2,2));
    FD[2][3]=(mVariables.F(0,1)*mVariables.ConstitutiveMatrix(1,2)+mVariables.F(0,0)*mVariables.ConstitutiveMatrix(2,2));

    FD[3].resize(4);
    FD[3][0]=(mVariables.F(1,0)*mVariables.ConstitutiveMatrix(0,0)+mVariables.F(1,1)*mVariables.ConstitutiveMatrix(2,0));
    FD[3][1]=(mVariables.F(1,0)*mVariables.ConstitutiveMatrix(0,1)+mVariables.F(1,1)*mVariables.ConstitutiveMatrix(2,1));
    FD[3][2]=(mVariables.F(1,0)*mVariables.ConstitutiveMatrix(0,2)+mVariables.F(1,1)*mVariables.ConstitutiveMatrix(2,2));
    FD[3][3]=(mVariables.F(1,0)*mVariables.ConstitutiveMatrix(0,2)+mVariables.F(1,1)*mVariables.ConstitutiveMatrix(2,2));

    std::vector<Vector> FDB(4);

    FDB[0]=ZeroVector(2);
    FDB[1]=ZeroVector(2);
    FDB[2]=ZeroVector(2);
    FDB[3]=ZeroVector(2);

    for(int i=0; i<4; i++)
    {
        for(int j=0; j<2; j++)
        {
            FDB[i][j]=FD[i][0]*mVariables.F(j,0)*mVariables.DN_DX(ndi,0)+FD[i][1]*mVariables.F(j,1)*mVariables.DN_DX(ndi,1)+
                      (FD[i][2]+FD[i][3])*(0.5)*(mVariables.F(j,0)*mVariables.DN_DX(ndi,1)+mVariables.F(j,1)*mVariables.DN_DX(ndi,0));

        }

    }

    // std::cout<<" FBD [0] "<<FDB[0]<<std::endl;
    // std::cout<<" FBD [1] "<<FDB[1]<<std::endl;
    // std::cout<<" FBD [2] "<<FDB[2]<<std::endl;
    // std::cout<<" FBD [3] "<<FDB[3]<<std::endl;

    for(int i=0; i<2; i++)
    {
        SigmaP[ndj][i]+=DirVector[0]*mVariables.Contact.ReferenceSurface.Normal[0]*(FDB[0][i])+
                        DirVector[1]*mVariables.Contact.ReferenceSurface.Normal[1]*(FDB[1][i])+
                        DirVector[0]*mVariables.Contact.ReferenceSurface.Normal[1]*(FDB[2][i])+
                        DirVector[1]*mVariables.Contact.ReferenceSurface.Normal[0]*(FDB[3][i]);
    }



}

//************************************************************************************
//************************************************************************************

void ContactDomainCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
							VectorType& rRightHandSideVector,
							bool CalculateStiffnessMatrixFlag,
							bool CalculateResidualVectorFlag )
{

  const unsigned int number_of_nodes = GetGeometry().size();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

  //resizing as needed the LHS
  unsigned int MatSize = (number_of_nodes + 1) * dimension;

  if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
      if ( rLeftHandSideMatrix.size1() != MatSize )
	rLeftHandSideMatrix.resize( MatSize, MatSize, false );

      noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


  //resizing as needed the RHS
  if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
      if ( rRightHandSideVector.size() != MatSize )
	rRightHandSideVector.resize( MatSize, false );

      rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
    }
}


//************************************************************************************
//************************************************************************************

void ContactDomainCondition::CalculateConditionalSystem( MatrixType& rLeftHandSideMatrix,
							   VectorType& rRightHandSideVector,
							   ProcessInfo& rCurrentProcessInfo,
							   bool CalculateStiffnessMatrixFlag,
							   bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    unsigned int vsize=GetValue(MASTER_ELEMENTS).size();
    Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS)[vsize-1];

    //std::cout<<"//******** CONTACT ELEMENT "<<this->Id()<<" ********// "<<std::endl;

    const unsigned int number_of_nodes = mpMasterGeometry->size();
    const unsigned int dimension       = mpMasterGeometry->WorkingSpaceDimension();

    unsigned int StrainSize;

    if ( dimension == 2 )
        StrainSize = 3;
    else
        StrainSize = 6;

    mVariables.F.resize( dimension, dimension );

    mVariables.ConstitutiveMatrix.resize( StrainSize, StrainSize );

    mVariables.StressVector.resize( StrainSize );

    mVariables.DN_DX.resize( number_of_nodes, dimension );

    
    InitializeSystemMatrices(rLeftHandSideMatrix,rRightHandSideVector,CalculateStiffnessMatrixFlag,CalculateResidualVectorFlag);

    //SET TANGENT DIRECTION-RELATIVE VELOCITY AND FRICTION PARAMETERS

    if( GetProperties()[PENALTY_CONTACT] == 1 )
	    mVariables.Contact.Options.Set(PENALTY);

    //Initialize friction parameter
    mVariables.Contact.FrictionCoefficient   =0;

    mVariables.Contact.Options.Set(COMPUTE_FRICTION_STIFFNESS);
    if( GetProperties()[FRICTION_ACTIVE] == 1 )
	    mVariables.Contact.Options.Set(COMPUTE_FRICTION_FORCES);



    //Tangent velocity and stablish friction parameter
    VectorType TangentVelocity (3,0.0);
    CalcRelativeVelocity    (TangentVelocity);
    CalcFrictionCoefficient (TangentVelocity);


    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = mpMasterGeometry->IntegrationPoints( mThisIntegrationMethod );
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

    mVariables.detF =MathUtils<double>::Det(mVariables.F);

    //Get Current ConstitutiveMatrix
    std::vector<Matrix> ConstitutiveMatrix(integration_points.size());   
    MasterElement.CalculateOnIntegrationPoints(CONSTITUTIVE_MATRIX,ConstitutiveMatrix,rCurrentProcessInfo);

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        CalculateKinematics(PointNumber,rCurrentProcessInfo);
	
	//set standart parameters
	mVariables.N=row(Ncontainer , PointNumber);

	//set constitutive matrix
	mVariables.ConstitutiveMatrix = ConstitutiveMatrix[PointNumber];

	//std::cout<<" ConstitutiveMatrix "<<mVariables.ConstitutiveMatrix<<std::endl;

	//Calculate Functions for the Tangent construction
	CalcDomainShapeN();

        double IntegrationWeight =0.5*mVariables.Contact.ReferenceBase[0].L;  //all components are multiplied by this
        if ( dimension == 2 ) IntegrationWeight *=  MasterElement.GetProperties()[THICKNESS];

        if(mVariables.Contact.Options.Is(ACTIVE))
        {

            //if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            //{
                //contributions to stiffness matrix calculated on the reference config

                // operation performed: add Km to the rLefsHandSideMatrix
	      if(mVariables.Contact.Options.Is(PENALTY)){
		std::cout<<" PenaltyStiffness Computed "<<std::endl;
                CalculateAndAddPenaltyKm( rLeftHandSideMatrix, IntegrationWeight );
	      }
	      else{
		CalculateAndAddKm( rLeftHandSideMatrix, IntegrationWeight );
	      }
  	    //}

	    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
            {
                // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
	      if(mVariables.Contact.Options.Is(PENALTY)){
		std::cout<<" PenaltyForces Computed "<<std::endl;
                CalculateAndAddContactPenaltyForces(rRightHandSideVector,IntegrationWeight);
	      }
	      else{
		CalculateAndAddContactForces(rRightHandSideVector,IntegrationWeight);
	      }

            }
        }
    }


    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

inline void ContactDomainCondition::CalcNormalForce (double &F,unsigned int& ndi,unsigned int& idir)
{    
    F = mVariables.Contact.NormalMultiplier*mVariables.Contact.dN_dn[ndi]*mVariables.Contact.CurrentSurface.Normal[idir];

}

//************************************************************************************
//************************************************************************************


inline void ContactDomainCondition::CalcTangentStickForce (double &F,unsigned int& ndi,unsigned int& idir)
{

	if( mVariables.Contact.Options.Is(COMPUTE_FRICTION_FORCES) )
	{
		F=mVariables.Contact.TangentialMultiplier*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir]);
	}
	else{
		F=0.0;
	}

}

//************************************************************************************
//************************************************************************************


inline void ContactDomainCondition::CalcTangentSlipForce (double &F,unsigned int& ndi,unsigned int& idir)
{

	if( mVariables.Contact.Options.Is(COMPUTE_FRICTION_FORCES) )
	{
		F=mVariables.Contact.NormalMultiplier*(mVariables.Contact.FrictionCoefficient*mVariables.Contact.TangentialGapsign)*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir]);

	}
	else{
		F=0.0;
	}



}

//************************************************************************************
//************************************************************************************

inline void ContactDomainCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
								   double& rIntegrationWeight)
{
    KRATOS_TRY

    //contributions to stiffness matrix calculated on the reference config
    unsigned int dimension=2;
    unsigned int size=rRightHandSideVector.size()/dimension;

    Vector Nforce;
    Vector Tforce;

    // std::cout<<" size "<<size<<std::endl;
    // for (unsigned int ndi=0; ndi<3; ndi++)
    //   std::cout<<" Node "<<GetGeometry()[ndi].Id()<<std::endl;

    VectorType NormalForce (3,0.0);
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
	    CalcNormalForce(Nforce[i],ndi,i);
            //TANGENT FORCE
            if(mVariables.Contact.Options.Is(NOT_SLIP))
            {
                CalcTangentStickForce(Tforce[i],ndi,i);
            }
            else
            {
                CalcTangentSlipForce(Tforce[i],ndi,i);
            }

	    rRightHandSideVector[index] -=(Nforce[i] + Tforce[i]);
 
	    NormalForce[i]  -= (Nforce[i])*rIntegrationWeight;
	    TangentForce[i] -= (Tforce[i])*rIntegrationWeight;

	    index++;
        }
	

	if(ndi<GetGeometry().PointsNumber())
	  {
	    VectorType & ForceContactNormal  = GetGeometry()[ndi].FastGetSolutionStepValue(FORCE_CONTACT_NORMAL);

	    VectorType & ForceContactTangent = GetGeometry()[ndi].FastGetSolutionStepValue(FORCE_CONTACT_TANGENT);

	    //Total force
	    //ForceContact+=ContactForce;	    	

	    //Only tangent force
	    ForceContactTangent += (TangentForce);

	    //Only normal force
	    ForceContactNormal  += (NormalForce);

	    //Normal Direction
	    //ForceContactTangent += ((-1)*mVariables.Contact.NormalMultiplier*mVariables.Contact.dN_dn[ndi]*mVariables.Contact.CurrentSurface.Normal*rIntegrationWeight);	 

	    // if(mVariables.Contact.NormalMultiplier>0)
	    //   std::cout<<" Element "<<this->Id()<<" positive NormalMultiplier "<<mVariables.Contact.NormalMultiplier<<std::endl;
	    		
	  }

	// if(ndi<this->GetGeometry().size())
	//   std::cout<<" Node:"<<this->GetGeometry()[ndi].Id()<<"[ Nforce "<<Nforce*rIntegrationWeight<<" Tforce "<<Tforce*rIntegrationWeight<<"]"<<std::endl;

	// std::cout<<" Nforce "<<Nforce*rIntegrationWeight<<std::endl;
	// std::cout<<" Tforce "<<Tforce*rIntegrationWeight<<std::endl;

    }


    rRightHandSideVector  *=  rIntegrationWeight;

    // std::cout<<std::endl;
    // std::cout<<" Fcontact "<<rRightHandSideVector<<std::endl;
    

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

inline void ContactDomainCondition::CalcNormalPenaltyForce (double &F,unsigned int& ndi,unsigned int& idir)
{    
    F = mVariables.Contact.gapN*mVariables.Contact.NormalMultiplier*mVariables.Contact.dN_dn[ndi]*mVariables.Contact.CurrentSurface.Normal[idir];

}

//************************************************************************************
//************************************************************************************


inline void ContactDomainCondition::CalcTangentStickPenaltyForce (double &F,unsigned int& ndi,unsigned int& idir)
{

	if( mVariables.Contact.Options.Is(COMPUTE_FRICTION_FORCES) )
	{
		F=mVariables.Contact.TangentialMultiplier*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir]);
	}
	else{
		F=0.0;
	}


    

}

//************************************************************************************
//************************************************************************************


inline void ContactDomainCondition::CalcTangentSlipPenaltyForce (double &F,unsigned int& ndi,unsigned int& idir)
{

	if( mVariables.Contact.Options.Is(COMPUTE_FRICTION_FORCES) )
	{
		F=mVariables.Contact.NormalMultiplier*(mVariables.Contact.FrictionCoefficient*mVariables.Contact.TangentialGapsign)*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir]);

	}
	else{
		F=0.0;
	}



}

//************************************************************************************
//************************************************************************************

inline void ContactDomainCondition::CalculateAndAddContactPenaltyForces(VectorType& rRightHandSideVector,
									  double& rIntegrationWeight)
{
    KRATOS_TRY

    //contributions to stiffness matrix calculated on the reference config
    unsigned int dimension=2;
    unsigned int size=rRightHandSideVector.size()/dimension;

    Vector Nforce;
    Vector Tforce;

    // std::cout<<" size "<<size<<std::endl;
    // for (unsigned int ndi=0; ndi<3; ndi++)
    //   std::cout<<" Node "<<GetGeometry()[ndi].Id()<<std::endl;

    VectorType NormalForce (3,0.0);
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
	    CalcNormalPenaltyForce(Nforce[i],ndi,i);
            //TANGENT FORCE
            if(mVariables.Contact.Options.Is(NOT_SLIP))
            {
                CalcTangentStickForce(Tforce[i],ndi,i);
            }
            else
            {
                CalcTangentSlipForce(Tforce[i],ndi,i);
            }

	    rRightHandSideVector[index] -=(Nforce[i] + Tforce[i]);
 
	    NormalForce[i]  -= (Nforce[i])*rIntegrationWeight;
	    TangentForce[i] -= (Tforce[i])*rIntegrationWeight;

	    index++;
        }
	

	if(ndi<GetGeometry().PointsNumber())
	  {
	    VectorType & ForceContactNormal  = GetGeometry()[ndi].FastGetSolutionStepValue(FORCE_CONTACT_NORMAL);
	    VectorType & ForceContactTangent = GetGeometry()[ndi].FastGetSolutionStepValue(FORCE_CONTACT_TANGENT);

	    //Total force
	    //ForceContact+=ContactForce;	    	

	    //Only tangent force
	    ForceContactTangent += (TangentForce);

	    //Only normal force
	    ForceContactNormal  += (NormalForce);

	    //Normal Direction
	    //ForceContactTangent += ((-1)*mVariables.Contact.NormalMultiplier*mVariables.Contact.dN_dn[ndi]*mVariables.Contact.CurrentSurface.Normal*rIntegrationWeight);	 

	    // if(mVariables.Contact.NormalMultiplier>0)
	    //   std::cout<<" Element "<<this->Id()<<" positive NormalMultiplier "<<mVariables.Contact.NormalMultiplier<<std::endl;
	    		
	  }

	// if(ndi<this->GetGeometry().size())
	//   std::cout<<" Node:"<<this->GetGeometry()[ndi].Id()<<"[ Nforce "<<Nforce*rIntegrationWeight<<" Tforce "<<Tforce*rIntegrationWeight<<"]"<<std::endl;

    }


    rRightHandSideVector  *=  rIntegrationWeight;

    // std::cout<<std::endl;
    // std::cout<<" Fcontact "<<rRightHandSideVector<<std::endl;
    

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void ContactDomainCondition::CalcContactStiffness (double &Kcont,unsigned int& ndi,unsigned int& ndj,unsigned int& idir,unsigned int& jdir)
{
    Kcont=0;

    //Normal contact contribution:
    //KI:
    Kcont= (0.5/(mVariables.Contact.Tau*mVariables.Contact.ReferenceBase[0].L))*mVariables.Contact.dN_dn[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dn[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]-
           mVariables.Contact.NormalMultiplier*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]+
                                   mVariables.Contact.dN_dn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]+
                                   mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dn[ndj]*mVariables.Contact.CurrentSurface.Tangent[jdir])-
           mVariables.Contact.CurrentTensil.Tangent*(mVariables.Contact.dN_dn[ndi]*mVariables.Contact.CurrentSurface.Normal[idir])*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir];

    // std::cout<<" ndi "<<ndi<<" ndj "<<ndj<<" i "<<idir<<" j "<<jdir<<std::endl;
    // std::cout<<" Kg "<<Kcont;
    // //std::cout<<" constant "<<constant<<" Kg "<<Kcont*constant;
    // double K1=Kcont;

    //KII:
    Kcont+= mVariables.Contact.dN_dn[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.Nsigma[ndj][jdir];


    //Stick contact contribution:
    if(mVariables.Contact.Options.Is(NOT_SLIP))
    {
	//std::cout<<" + stick ";
        if(mVariables.Contact.Options.Is(COMPUTE_FRICTION_STIFFNESS))
        {
	    //std::cout<<"(mu_on)";
            //KI:
            Kcont+= (0.5/(mVariables.Contact.Tau*mVariables.Contact.ReferenceBase[0].L))*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir])*
                    (mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]+mVariables.Contact.dN_drn[ndj]*mVariables.Contact.CurrentSurface.Tangent[jdir])+ mVariables.Contact.TangentialMultiplier*(mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dn[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]+ mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]- mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]- mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Tangent[jdir])+ mVariables.Contact.CurrentTensil.Normal*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir])*(mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]);

            //KII:
            Kcont+= (mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir])*mVariables.Contact.Tsigma[ndj][jdir];

            //(he_a*ddN_dt_a(in)*ne_a(idime) + hddN_dn_n(in)*te_a(idime))*raux
        }
    }
    else
    {
        //Slip contact contribution:
	//std::cout<<" + slip ";
        if(mVariables.Contact.Options.Is(COMPUTE_FRICTION_STIFFNESS))
        {
	    //std::cout<<"(mu_on)";
            //KI:
            Kcont+= (mVariables.Contact.FrictionCoefficient*mVariables.Contact.TangentialGapsign)*((0.5/(mVariables.Contact.Tau*mVariables.Contact.ReferenceBase[0].L))*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir])* (mVariables.Contact.dN_dn[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir])+
                    mVariables.Contact.NormalMultiplier*(mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dn[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]+ mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]- mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]- mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Tangent[jdir])- mVariables.Contact.CurrentTensil.Tangent*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir])*(mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]));


            //KII:
            Kcont+= (mVariables.Contact.FrictionCoefficient*mVariables.Contact.TangentialGapsign)*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir])*mVariables.Contact.Nsigma[ndj][jdir];

        }
    }

    //std::cout<<" Ks "<<Kcont-K1<<std::endl;

}


//************************************************************************************
//************************************************************************************

void ContactDomainCondition::CalcContactPenaltyStiffness (double &Kcont,unsigned int& ndi,unsigned int& ndj,unsigned int& idir,unsigned int& jdir)
{
    Kcont=0;

    
    //Normal contact penalty contribution:
    //KI:
    double KN = mVariables.Contact.NormalMultiplier;

    Kcont= KN*( (mVariables.Contact.dN_dn[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dn[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir])-
		(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dn[ndj]*mVariables.Contact.CurrentSurface.Tangent[jdir])-
		(mVariables.Contact.gapN*mVariables.Contact.dN_dn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir])-
		((mVariables.Contact.gapN*mVariables.Contact.gapN)*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]));

    // std::cout<<" ndi "<<ndi<<" ndj "<<ndj<<" i "<<idir<<" j "<<jdir<<std::endl;
    // std::cout<<" Kg "<<Kcont;
    // //std::cout<<" constant "<<constant<<" Kg "<<Kcont*constant;
    // double K1=Kcont;


    //Stick contact contribution:
    // if(mVariables.Contact.Options.Is(NOT_SLIP))
    // {
    // 	//std::cout<<" + stick ";
    //     if(mVariables.Contact.Options.Is(COMPUTE_FRICTION_STIFFNESS))
    //     {
    // 	    //std::cout<<"(mu_on)";
    //         //KI:
    //         Kcont+= (0.5/(mVariables.Contact.Tau*mVariables.Contact.ReferenceBase[0].L))*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir])*
    //                 (mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]+mVariables.Contact.dN_drn[ndj]*mVariables.Contact.CurrentSurface.Tangent[jdir])+ mVariables.Contact.TangentialMultiplier*(mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dn[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]+ mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]- mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]- mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Tangent[jdir])+ mVariables.Contact.CurrentTensil.Normal*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir])*(mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]);

    //         //KII:
    //         Kcont+= (mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir])*mVariables.Contact.Tsigma[ndj][jdir];

    //         //(he_a*ddN_dt_a(in)*ne_a(idime) + hddN_dn_n(in)*te_a(idime))*raux
    //     }
    // }
    // else
    // {
    //     //Slip contact contribution:
    // 	//std::cout<<" + slip ";
    //     if(mVariables.Contact.Options.Is(COMPUTE_FRICTION_STIFFNESS))
    //     {
    // 	    //std::cout<<"(mu_on)";
    //         //KI:
    //         Kcont+= (mVariables.Contact.FrictionCoefficient*mVariables.Contact.TangentialGapsign)*((0.5/(mVariables.Contact.Tau*mVariables.Contact.ReferenceBase[0].L))*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir])* (mVariables.Contact.dN_dn[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir])+
    //                 mVariables.Contact.NormalMultiplier*(mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dn[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]+ mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]- mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]- mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]*mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Tangent[jdir])- mVariables.Contact.CurrentTensil.Tangent*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir])*(mVariables.Contact.dN_dt[ndj]*mVariables.Contact.CurrentSurface.Normal[jdir]));


    //         //KII:
    //         Kcont+= (mVariables.Contact.FrictionCoefficient*mVariables.Contact.TangentialGapsign)*(mVariables.Contact.gapN*mVariables.Contact.dN_dt[ndi]*mVariables.Contact.CurrentSurface.Normal[idir]+mVariables.Contact.dN_drn[ndi]*mVariables.Contact.CurrentSurface.Tangent[idir])*mVariables.Contact.Nsigma[ndj][jdir];

    //     }
    // }

    //std::cout<<" Ks "<<Kcont-K1<<std::endl;

}

//************************************************************************************
//************************************************************************************


void ContactDomainCondition::CalculateAndAddKm(MatrixType& rK,
						 double& rIntegrationWeight)
  
{ 
    KRATOS_TRY

    //contributions to stiffness matrix calculated on the reference config
    unsigned int dimension=2;
    unsigned int size=rK.size1()/dimension;
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
		    CalcContactStiffness(kcont,ndi,ndj,i,j);
		    rK(ndi*2+i,ndj*2+j)+=kcont;
                }
            }
        }
    }

    rK *= rIntegrationWeight;

    // std::cout<<std::endl;
    // std::cout<<" Kcontact "<<rK<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************


void ContactDomainCondition::CalculateAndAddPenaltyKm(MatrixType& rK,
							double& rIntegrationWeight)
							
{    
    KRATOS_TRY

    //contributions to stiffness matrix calculated on the reference config
    unsigned int dimension=2;
    unsigned int size=rK.size1()/dimension;
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
		    CalcContactPenaltyStiffness(kcont,ndi,ndj,i,j);
		    rK(ndi*2+i,ndj*2+j)+=kcont;
                }
            }
        }
    }

    rK *= rIntegrationWeight;

    // std::cout<<std::endl;
    // std::cout<<" Kcontact "<<rK<<std::endl;

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

        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {

            if ( rOutput[PointNumber].size2() != StrainSize )
                rOutput[PointNumber].resize( 1 ,  StrainSize , false );


            mVariables.detF =MathUtils<double>::Det(mVariables.F);

            Matrix StressMatrix ( dimension, dimension );
            StressMatrix = mConstitutiveLawVector[PointNumber]->GetValue( rVariable , StressMatrix );

            StressMatrix = mConstitutiveLawVector[PointNumber]->TransformStresses(StressMatrix,mVariables.F,mVariables.detF,ConstitutiveLaw::StressMeasure_Cauchy,ConstitutiveLaw::StressMeasure_PK2);

            Vector StressVector ( StrainSize );
            StressVector = MathUtils<double>::StressTensorToVector( StressMatrix );

            for ( unsigned int ii = 0; ii < StressVector.size(); ii++ )
            {
                rOutput[PointNumber]( 0, ii ) = StressVector[ii];
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

inline void ContactDomainCondition::CalcRelativeVelocity (VectorType & TangentVelocity)
{
    //if current tangent is not previously computed, do it here.
    mVariables.Contact.CurrentSurface.Tangent = CalcCurrentTangent( mVariables.Contact.CurrentSurface.Tangent );


   if(double(inner_prod(mVariables.Contact.CurrentSurface.Tangent,mVariables.Contact.ReferenceSurface.Tangent))<0) //to give the correct direction
        mVariables.Contact.CurrentSurface.Tangent*=-1;

   //std::cout<<" Normal ["<<this->Id()<<"] :"<<mVariables.Contact.CurrentSurface.Normal<<std::endl;
   //std::cout<<" Tangent ["<<this->Id()<<"] :"<<mVariables.Contact.CurrentSurface.Tangent<<std::endl;

    // (Tangent vector previously computed)
    const int number_of_nodes = GetGeometry().size();

    //compute relative velocities
    int slave=mVariables.slaves[0];
    VectorType  CurrentVelocity;
    for (int i = 0; i < number_of_nodes; i++ )
    {
        //Displacement from the reference to the current configuration
        CurrentVelocity  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        if(i!=slave)
            CurrentVelocity *=(-0.5);

        TangentVelocity+=CurrentVelocity;
    }

    //Relative tangent movement of the slave if the master is fixed (the direction is implicit in the method)
    TangentVelocity =  mVariables.Contact.CurrentSurface.Tangent*(inner_prod(TangentVelocity,mVariables.Contact.CurrentSurface.Tangent));

    //Filter for low velocities (relatives to dynamic waves)
    CurrentVelocity.clear();
    CalcRelativeDisplacement(CurrentVelocity);

    if(norm_2(TangentVelocity)<1e-2*(norm_2(CurrentVelocity)/norm_2(TangentVelocity)))
    {
      TangentVelocity.clear();
    }

    //std::cout<<" TangentVelocity ["<<this->Id()<<"] :"<<TangentVelocity<<std::endl;

}

//************************************************************************************
//************************************************************************************


inline void ContactDomainCondition::CalcRelativeDisplacement (VectorType & TangentDisplacement)
{

    // (Tangent vector previously computed)
    const int number_of_nodes = GetGeometry().size();

    //compute relative displacements
    int slave=mVariables.slaves[0];
    VectorType CurrentDisplacement;
    for (int i = 0; i < number_of_nodes; i++ )
    {
        //Displacement from the reference to the current configuration
        CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        if(i!=slave)
            CurrentDisplacement *=(-0.5);

        TangentDisplacement+=CurrentDisplacement;
    }

    //Relative tangent movement of the slave if the master is fixed (the direction is implicit in the method)
    TangentDisplacement = mVariables.Contact.CurrentSurface.Tangent*(inner_prod(TangentDisplacement,mVariables.Contact.CurrentSurface.Tangent));

}

//************************************************************************************
//************************************************************************************


void ContactDomainCondition::CalcFrictionCoefficient (const VectorType & TangentVelocity)
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


        mVariables.Contact.FrictionCoefficient=muDynamic+(muStatic-muStatic)*exp((-1)*paramC*fabs(Velocity));


        //if (TangentVelocity.modulus()>paramE){

        //1.- square root regularization
        mVariables.Contact.FrictionCoefficient*=fabs(Velocity)/sqrt((Velocity*Velocity)+(paramE*paramE));//square root

        //2.-hyperbolic regularization
        //FrictionCoefficient*=tanh(fabs(Velocity)/paramE);

        //}


    }
    else
    {
	mVariables.Contact.FrictionCoefficient=muStatic;
    }

    //Activate or deactivate friction (simulation type)

    if( mVariables.Contact.Options.IsNot(COMPUTE_FRICTION_FORCES) ){
	    mVariables.Contact.FrictionCoefficient = 0;
    }

    //std::cout<<" friction coefficient "<<mVariables.Contact.FrictionCoefficient<<std::endl;

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


void ContactDomainCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    int IntMethod = (int)mThisIntegrationMethod;
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("MasterGeometry",mpMasterGeometry);
    //rSerializer.save("StandardVariables",mVariables);
}

void ContactDomainCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("MasterGeometry",mpMasterGeometry);
    //rSerializer.load("StandardVariables",mVariables);

}



} // Namespace Kratos


