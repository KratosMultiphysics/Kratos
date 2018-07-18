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
#include "custom_elements/thermal_elements/axisymmetric_thermal_element.hpp"

#include "solid_mechanics_application_variables.h"

//#include <omp.h>

namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymmetricThermalElement::AxisymmetricThermalElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : ThermalElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

AxisymmetricThermalElement::AxisymmetricThermalElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : ThermalElement( NewId, pGeometry, pProperties )
{
  //const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  //mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
  //mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

AxisymmetricThermalElement::AxisymmetricThermalElement( AxisymmetricThermalElement const& rOther)
    :ThermalElement(rOther)
{
}



//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

AxisymmetricThermalElement&  AxisymmetricThermalElement::operator=(AxisymmetricThermalElement const& rOther)
{
    ThermalElement::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer AxisymmetricThermalElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared<AxisymmetricThermalElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

AxisymmetricThermalElement::~AxisymmetricThermalElement()
{
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void AxisymmetricThermalElement::CalculateKinematics(GeneralVariables& rVariables,
        const double& rPointNumber)

{
    KRATOS_TRY

    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    Matrix InvJ;
    //Calculating the inverse of the jacobian and the parameters needed
    MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    //Compute cartesian derivatives
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber] , InvJ );

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber );

    //Calculate IntegrationPoint radius
    CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius, rVariables.N);

    Matrix Invj;
    //Calculating the inverse of the jacobian and the parameters needed
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ); //overwrites detJ

    //Compute cartesian derivatives
    rVariables.DN_DX = prod( DN_De[rPointNumber] , Invj ); //overwrites DX now is the current position



    KRATOS_CATCH( "" )
}


//*************************COMPUTE AXYSIMMETRIC RADIUS********************************
//************************************************************************************
void AxisymmetricThermalElement::CalculateRadius(double & rCurrentRadius,
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






//***********************COMPUTE LOCAL SYSTEM CONTRIBUTIONS***************************
//************************************************************************************

void AxisymmetricThermalElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
{
  double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

  //contributions to stiffness matrix calculated on the reference config

  ThermalElement::CalculateAndAddLHS( rLeftHandSideMatrix, rVariables, IntegrationWeight );

}


//************************************************************************************
//************************************************************************************

void AxisymmetricThermalElement::CalculateAndAddRHS(VectorType& rRightHandSideVector, GeneralVariables& rVariables, double& rHeatSource, double& rIntegrationWeight)
{

    double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

    //contribution to external forces

    ThermalElement::CalculateAndAddRHS( rRightHandSideVector, rVariables, rHeatSource, IntegrationWeight );

}



void AxisymmetricThermalElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ThermalElement );
}

void AxisymmetricThermalElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ThermalElement );
}





} // Namespace Kratos


