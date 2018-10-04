//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:               MCaicedo $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/solid_elements/small_displacement_bbar_element.hpp"
#include "includes/constitutive_law.h"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementBbarElement::SmallDisplacementBbarElement()
  : SmallDisplacementElement()
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementBbarElement::SmallDisplacementBbarElement(IndexType NewId,
                                                           GeometryType::Pointer pGeometry)
    : SmallDisplacementElement(NewId, pGeometry)
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementBbarElement::SmallDisplacementBbarElement(IndexType NewId,
                                                           GeometryType::Pointer pGeometry,
                                                           PropertiesType::Pointer pProperties)
    : SmallDisplacementElement(NewId, pGeometry, pProperties)
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

SmallDisplacementBbarElement::SmallDisplacementBbarElement(SmallDisplacementBbarElement const& rOther)
    : SmallDisplacementElement(rOther)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

SmallDisplacementBbarElement& SmallDisplacementBbarElement::operator=(SmallDisplacementBbarElement const& rOther)
{
    SmallDisplacementElement::operator=(rOther);

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer SmallDisplacementBbarElement::Create(IndexType NewId,
                                                      NodesArrayType const& rThisNodes,
                                                      PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared< SmallDisplacementBbarElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer SmallDisplacementBbarElement::Clone(IndexType NewId,
                                                     NodesArrayType const& rThisNodes) const
{
    SmallDisplacementBbarElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if (NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size())
    {
        NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

        if (NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber())
            KRATOS_THROW_ERROR(std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size());
    }

    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
    }

    //-----------//

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Kratos::make_shared< SmallDisplacementBbarElement >(NewElement);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementBbarElement::~SmallDisplacementBbarElement()
{
}


//************************************************************************************
//************************************************************************************


void SmallDisplacementBbarElement::InitializeElementData(ElementDataType& rVariables,
                                                              const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    SmallDisplacementElement::InitializeElementData(rVariables, rCurrentProcessInfo);

    // calculate volumetric deformation matrix (stored in rVariables.H)
    CalculateVolumetricDeformationMatrix(rVariables);

    KRATOS_CATCH( "" )
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


void SmallDisplacementBbarElement::CalculateKinematics(ElementDataType& rVariables, const double& rPointNumber)
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

    //Compute cartesian derivatives  [dN/dx_n]
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber] , InvJ );

    //Set Shape Functions Values for this integration point
    noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    // Compute the deformation matrix B_bar
    this->CalculateDeformationMatrixBbar(rVariables.B, rVariables.H, rVariables.DN_DX);

    // Compute infinitessimal B_bar strain
    this->CalculateInfinitesimalStrainBbar(rVariables.StrainVector, rVariables.H, rVariables.DN_DX);

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementBbarElement::CalculateInfinitesimalStrainBbar(Vector& rStrainVector,
								    const Matrix& rH,
								    const Matrix& rDN_DX)
{
    KRATOS_TRY

    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    Matrix J (dimension,dimension);

    //Displacement gradient J [dU/dx_n]
    this->CalculateDisplacementGradient( J, rDN_DX );

    //Compute infinitessimal strain
    this->CalculateInfinitesimalStrain( J, rStrainVector );

    //Add Bbar terms:

    double VolumetricStrain = 0;
    for( SizeType i=0; i<dimension; i++ )
      {
	VolumetricStrain -= rStrainVector[i];
      }
    VolumetricStrain /= double(dimension);


    Vector Displacements;
    this->GetValuesVector(Displacements,0);

    double DeformationVolume = 0;
    for(unsigned int i =0; i<Displacements.size(); i++)
       DeformationVolume += rH(i,0) * Displacements[i];


    VolumetricStrain += DeformationVolume;

    for( SizeType i=0; i<dimension; i++ )
      {
	rStrainVector[i] += VolumetricStrain;
      }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBbarElement::CalculateDeformationMatrixBbar(Matrix& rB,
                                                                  const Matrix& rH,
                                                                  const Matrix& rDN_DX)
{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    // Compute deformation matrix
    const GeometryType& rGeometry = GetGeometry();
    ElementUtilities::CalculateLinearDeformationMatrix(rB,rGeometry,rDN_DX);

    if (dimension == 2)
    {
        unsigned int index = 0;
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
            index = 2 * i;

            rB( 0, index + 0 ) += rH( index  , 0 ) - 0.5 * rDN_DX( i, 0 );
	    rB( 0, index + 1 ) += rH( index+1, 0 ) - 0.5 * rDN_DX( i, 1 );
            rB( 1, index + 1 ) += rH( index+1, 0 ) - 0.5 * rDN_DX( i, 1 );
	    rB( 1, index + 0 ) += rH( index  , 0 ) - 0.5 * rDN_DX( i, 0 );
        }
    }
    else if (dimension == 3)
    {
      double athird = 1.0/3.0;
      unsigned int index = 0;

      for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
	  index = 3 * i;

	  rB( 0, index + 0 ) += rH( index  , 0 ) - athird * ( rDN_DX( i, 0 ) );
	  rB( 0, index + 1 ) += rH( index+1, 0 ) - athird * ( rDN_DX( i, 1 ) );
	  rB( 0, index + 2 ) += rH( index+2, 0 ) - athird * ( rDN_DX( i, 2 ) );

	  rB( 1, index + 1 ) += rH( index+1, 0 ) - athird * ( rDN_DX( i, 1 ) );
	  rB( 1, index + 0 ) += rH( index  , 0 ) - athird * ( rDN_DX( i, 0 ) );
	  rB( 1, index + 2 ) += rH( index+2, 0 ) - athird * ( rDN_DX( i, 2 ) );

	  rB( 2, index + 2 ) += rH( index+2, 0 ) - athird * ( rDN_DX( i, 2 ) );
	  rB( 2, index + 0 ) += rH( index  , 0 ) - athird * ( rDN_DX( i, 0 ) );
	  rB( 2, index + 1 ) += rH( index+1, 0 ) - athird * ( rDN_DX( i, 1 ) );
	}

    }
    else
    {
        KRATOS_THROW_ERROR(std::invalid_argument,
                           "something is wrong with the dimension", "")
    }


    KRATOS_CATCH("")
}

void SmallDisplacementBbarElement::CalculateVolumetricDeformationMatrix(ElementDataType& rVariables)
{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    rVariables.H.resize(dimension * number_of_nodes, 1, false);
    rVariables.H.clear();

    // reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    double GeometrySize = 0.0;

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
      {
	const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

	Matrix InvJ;
	MathUtils<double>::InvertMatrix(rVariables.J[PointNumber], InvJ, rVariables.detJ);
	noalias(rVariables.DN_DX) = prod(DN_De[PointNumber], InvJ);

	double IntegrationWeight = integration_points[PointNumber].Weight() * rVariables.detJ;
	GeometrySize += IntegrationWeight;

	unsigned int index = 0;
	for (SizeType i = 0; i < number_of_nodes; i++)
	  {
	    for (SizeType j= 0; j < dimension; j++)
	      {
		index = i * dimension + j;
		rVariables.H(index,0) += rVariables.DN_DX(i,j) * IntegrationWeight;
	      }
	  }
      }

    rVariables.H /= (GeometrySize*dimension);

    KRATOS_CATCH("")
}

void SmallDisplacementBbarElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementElement)
}

void SmallDisplacementBbarElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacementElement)
}

} // Namespace Kratos
