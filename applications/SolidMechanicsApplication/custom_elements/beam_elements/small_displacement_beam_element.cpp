//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/small_displacement_beam_element.hpp"

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
    return Kratos::make_shared<SmallDisplacementBeamElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  SmallDisplacementBeamElement::~SmallDisplacementBeamElement()
  {
  }


  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::InitializeElementData(ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
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

  void SmallDisplacementBeamElement::CalculateKinematics(ElementDataType& rVariables, const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //point number
    rVariables.PointNumber = rPointNumber;

    //Set Shape Functions Values for this integration point
    noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    Vector Jacobian(3);
    noalias(Jacobian) = ZeroVector(3);

    //Calculating the jacobian [dx_0/d£]
    for( SizeType i=0; i<dimension; i++)
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

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int voigt_size = dimension * (dimension +1) * 0.5;

    if ( rB.size1() != voigt_size || rB.size2() != dimension*number_of_nodes )
      rB.resize(voigt_size, ( (dimension-1) * 3 ) * number_of_nodes, false );

    if( dimension == 2 )
      {

	for ( SizeType i = 0; i < number_of_nodes; i++ )
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

	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    index = ((dimension-1) * 3 ) * i;

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

  void SmallDisplacementBeamElement::CalculateConstitutiveMatrix(ElementDataType& rVariables)
  {
    KRATOS_TRY

    this->CalculateMaterialConstitutiveMatrix(rVariables.ConstitutiveMatrix, rVariables);

    //std::cout<<" ConstitutiveMatrix "<<rConstitutiveMatrix<<std::endl;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************
  void SmallDisplacementBeamElement::CalculateStressResultants(ElementDataType& rVariables,
							       const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    Vector NodalDofs;
    this->GetValuesVector(NodalDofs,0);

    rVariables.StrainVector = prod( rVariables.B, NodalDofs );

    Matrix ConstitutiveMatrix;
    this->CalculateMaterialConstitutiveMatrix(ConstitutiveMatrix, rVariables);

    //Reference Stress Vector
    rVariables.StressVector = prod( ConstitutiveMatrix, rVariables.StrainVector );


    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
							 ElementDataType& rVariables,
							 double& rIntegrationWeight)
  {
    KRATOS_TRY

    MatrixType LocalStiffnessMatrix = rLeftHandSideMatrix;
    this->CalculateLocalStiffnessMatrix(LocalStiffnessMatrix,rVariables);

    std::cout<<" Direct Stiffness Matrix  "<<LocalStiffnessMatrix<<std::endl;

    //contributions to stiffness matrix calculated on the reference config
    noalias( rLeftHandSideMatrix ) += prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) ); //to be optimized to remove the temporary

    std::cout<<" Stiffness Matrix Timoshenko "<<rLeftHandSideMatrix<<std::endl;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
								   ElementDataType & rVariables,
								   double& rIntegrationWeight)
  {
    KRATOS_TRY

    //contributions to the internal force vector calculated on the reference config
    noalias( rRightHandSideVector ) += rIntegrationWeight * prod( trans( rVariables.B ),  rVariables.StressVector); //to be optimized to remove the temporary

    std::cout<<" Internal Forces Timoshenko "<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateLocalStiffnessMatrix(Matrix& LocalStiffnessMatrix, ElementDataType& rVariables)
  {
    KRATOS_TRY

    const double PoissonCoefficient = GetProperties()[POISSON_RATIO];
    const double YoungModulus       = GetProperties()[YOUNG_MODULUS];
    const double ShearModulus       = YoungModulus*0.5/(1.0 + PoissonCoefficient);

    const double L   = rVariables.Length;
    const double LL  = rVariables.Length * rVariables.Length;
    const double LLL = rVariables.Length * rVariables.Length * rVariables.Length;

    double const EA  = rVariables.Section.Area * YoungModulus;
    double const EIz = rVariables.Section.Inertia_z  * YoungModulus;
    double const EIy = rVariables.Section.Inertia_y  * YoungModulus;
    double const JG  = rVariables.Section.Polar_Inertia * ShearModulus;

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

  void SmallDisplacementBeamElement::MapLocalToGlobal(ElementDataType& rVariables, MatrixType& rMatrix)
  {
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

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

  void SmallDisplacementBeamElement::MapLocalToGlobal(ElementDataType& rVariables, VectorType& rVector)
  {
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

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

    QuaternionType CurrentLocalQuaternion = QuaternionType::FromRotationMatrix( rRotationMatrix );

    CurrentLocalQuaternion = mInitialLocalQuaternion.conjugate() * CurrentLocalQuaternion;

    CurrentLocalQuaternion.ToRotationMatrix(rRotationMatrix);

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

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = BeamElement::Check(rCurrentProcessInfo);

    return ErrorCode;

    KRATOS_CATCH( "" )
  }




} // Namespace Kratos


