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
#include "custom_conditions/axisym_contact_domain_penalty_2D_condition.hpp"

#include "contact_mechanics_application_variables.h"


namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymContactDomainPenalty2DCondition::AxisymContactDomainPenalty2DCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : ContactDomainPenalty2DCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymContactDomainPenalty2DCondition::AxisymContactDomainPenalty2DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : ContactDomainPenalty2DCondition( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

AxisymContactDomainPenalty2DCondition::AxisymContactDomainPenalty2DCondition( AxisymContactDomainPenalty2DCondition const& rOther)
    :ContactDomainPenalty2DCondition(rOther)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

AxisymContactDomainPenalty2DCondition&  AxisymContactDomainPenalty2DCondition::operator=(AxisymContactDomainPenalty2DCondition const& rOther)
{
    ContactDomainPenalty2DCondition::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Condition::Pointer AxisymContactDomainPenalty2DCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared<AxisymContactDomainPenalty2DCondition>(NewId, GetGeometry().Create( ThisNodes ), pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer AxisymContactDomainPenalty2DCondition::Clone( IndexType NewId, NodesArrayType const& ThisNodes ) const
{
  return this->Create(NewId, ThisNodes, pGetProperties());
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************


AxisymContactDomainPenalty2DCondition::~AxisymContactDomainPenalty2DCondition()
{
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//************************************************************************************
//************************************************************************************

void AxisymContactDomainPenalty2DCondition::InitializeConditionVariables (ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType & MasterGeometry = mContactVariables.GetMasterGeometry();

    const unsigned int number_of_nodes = MasterGeometry.size();
    const unsigned int dimension       = MasterGeometry.WorkingSpaceDimension();

    unsigned int voigtsize = 4;

    rVariables.F.resize( dimension, dimension );

    rVariables.ConstitutiveMatrix.resize( voigtsize, voigtsize );

    rVariables.StressVector.resize( voigtsize );

    rVariables.DN_DX.resize( number_of_nodes, dimension );

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    // UL
    //Calculate Delta Position
    //rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    //rVariables.J = MasterGeometry.Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

     // SL
    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = MasterGeometry.Jacobian( rVariables.j, mThisIntegrationMethod );

}

//*********************************COMPUTE RADIUS*************************************
//************************************************************************************

void AxisymContactDomainPenalty2DCondition::CalculateRadius(double & rCurrentRadius,
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
            //std::cout<<" node "<<i<<" -> DeltaDisplacement : "<<DeltaDisplacement<<std::endl;
        }
    }


    if ( dimension == 3 )
    {
      std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
    }

    KRATOS_CATCH( "" )
}


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void AxisymContactDomainPenalty2DCondition::CalculateKinematics( ConditionVariables& rVariables, ProcessInfo& rCurrentProcessInfo, const unsigned int& rPointNumber )
{
    KRATOS_TRY


    ElementType&  MasterElement  = mContactVariables.GetMasterElement();
    GeometryType& MasterGeometry = mContactVariables.GetMasterGeometry();

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = MasterGeometry.ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    //Get integration points number
    unsigned int integration_points_number = MasterGeometry.IntegrationPointsNumber( MasterElement.GetIntegrationMethod() );

    unsigned int voigtsize = 4;
    int dimension = GetGeometry().WorkingSpaceDimension();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);


    //Calculate radius
    rVariables.CurrentRadius = 0;
    rVariables.ReferenceRadius = 0;
    CalculateRadius( rVariables.CurrentRadius, rVariables.ReferenceRadius, rVariables.N );

   // UL
    // //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    // Matrix Invj;
    // MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    // //Compute cartesian derivatives [dN/dx_n]
    // noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber] , InvJ );

    // SL
    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ ); //overwrites detJ

    //Compute cartesian derivatives [dN/dx_n+1]
    rVariables.DN_DX = prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

    //Get Current DeformationGradient
    std::vector<Matrix> DeformationGradientVector ( integration_points_number );
    DeformationGradientVector[rPointNumber]=identity_matrix<double>( dimension );
    MasterElement.GetValueOnIntegrationPoints(DEFORMATION_GRADIENT,DeformationGradientVector,rCurrentProcessInfo);
    rVariables.F = DeformationGradientVector[rPointNumber];

    rVariables.detF = MathUtils<double>::Det(rVariables.F);

    //Get Current Stress
    std::vector<Vector> StressVector ( integration_points_number );
    StressVector[rPointNumber]=ZeroVector(voigtsize);
    //MasterElement.GetValueOnIntegrationPoints(PK2_STRESS_VECTOR,StressVector,rCurrentProcessInfo);
    MasterElement.GetValueOnIntegrationPoints(CAUCHY_STRESS_VECTOR,StressVector,rCurrentProcessInfo);

    // for( unsigned int i=0; i<StressVector.size(); i++)
    //   {
    // 	StressVector[i] = mConstitutiveLawVector[rPointNumber]->TransformStresses(StressVector[i], rVariables.F, rVariables.detF, ConstitutiveLaw::StressMeasure_Cauchy, ConstitutiveLaw::StressMeasure_PK2);
    //   }

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


//************************************************************************************
//************************************************************************************

void AxisymContactDomainPenalty2DCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
{
  ElementType&  MasterElement  = mContactVariables.GetMasterElement();
  // UL
  //double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.ReferenceRadius;
  //if ( MasterElement.GetProperties().Has(THICKNESS) )
  //   rIntegrationWeight /= MasterElement.GetProperties()[THICKNESS];

  // SL
  double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;
  if ( MasterElement.GetProperties().Has(THICKNESS) )
      rIntegrationWeight /= MasterElement.GetProperties()[THICKNESS];

  ContactDomainCondition::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

  //KRATOS_WATCH( rLeftHandSideMatrix )
}


//************************************************************************************
//************************************************************************************

void AxisymContactDomainPenalty2DCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
{
  ElementType&  MasterElement  = mContactVariables.GetMasterElement();
  // UL
  //double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.ReferenceRadius;
  //if ( MasterElement.GetProperties().Has(THICKNESS) )
  //   rIntegrationWeight /= MasterElement.GetProperties()[THICKNESS];

  // SL
  double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius;
  if ( MasterElement.GetProperties().Has(THICKNESS) )
      rIntegrationWeight /= MasterElement.GetProperties()[THICKNESS];

  ContactDomainCondition::CalculateAndAddRHS( rLocalSystem, rVariables, IntegrationWeight );

  //KRATOS_WATCH( rRightHandSideVector )
}

//************************************************************************************
//************************************************************************************


void AxisymContactDomainPenalty2DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ContactDomainPenalty2DCondition )
}

void AxisymContactDomainPenalty2DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ContactDomainPenalty2DCondition )
}



} // Namespace Kratos
