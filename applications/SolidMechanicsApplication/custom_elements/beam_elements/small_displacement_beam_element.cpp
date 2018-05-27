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
    return Element::Pointer(new SmallDisplacementBeamElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  SmallDisplacementBeamElement::~SmallDisplacementBeamElement()
  {
  }
 

  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::InitializeElementData(ElementDataPointerType& pVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType& dimension       = this->Dimension();
    const unsigned int voigt_size      = dimension * (dimension +1) * 0.5;
    
    pVariables->Initialize(voigt_size,dimension,number_of_nodes);
    
    //Compute Section Properties:
    this->CalculateSectionProperties(pVariables->Section);

    pVariables->Length = GetGeometry().Length();

    if(pVariables->Length == 0.00)
      KRATOS_ERROR << "Zero length found in element #" << this->Id() << std::endl;

    
    //set variables including all integration points values

    //reading shape functions
    pVariables->SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    pVariables->SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    //pVariables->j = GetGeometry().Jacobian( pVariables->j, mThisIntegrationMethod );

    //Calculate Delta Position
    pVariables->DeltaPosition = this->CalculateTotalDeltaPosition(pVariables->DeltaPosition);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    pVariables->J = GetGeometry().Jacobian( pVariables->J, mThisIntegrationMethod, pVariables->DeltaPosition );

    KRATOS_CATCH( "" )
  }
  
  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateKinematics(ElementDataPointerType& pVariables, const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    const SizeType& dimension       = this->Dimension();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = pVariables->GetShapeFunctions();

    //point number
    pVariables->PointNumber = rPointNumber;
      
    //Set Shape Functions Values for this integration point
    noalias(pVariables->N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);
    
    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = pVariables->GetShapeFunctionsGradients();

    Vector Jacobian(3);
    noalias(Jacobian) = ZeroVector(3);
    
    //Calculating the jacobian [dx_0/d£]
    for( SizeType i=0; i<dimension; i++)
      {
	Jacobian[i] = pVariables->J[rPointNumber](i,0);
      }

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_0]
    pVariables->detJ  =  norm_2(Jacobian);
    pVariables->DN_DX = DN_De[rPointNumber] * 1.0/pVariables->detJ; 

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(pVariables->B, pVariables->N, pVariables->DN_DX);
    
    
    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************
  void SmallDisplacementBeamElement::CalculateDeformationMatrix(Matrix& rB, const Vector& rN, const Matrix& rDN_DX)
  {
    KRATOS_TRY
      
    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType& dimension       = this->Dimension();
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

  void SmallDisplacementBeamElement::CalculateConstitutiveMatrix(ElementDataPointerType& pVariables)
  {
    KRATOS_TRY

    this->CalculateMaterialConstitutiveMatrix(pVariables->ConstitutiveMatrix, pVariables);
    
    //std::cout<<" ConstitutiveMatrix "<<rConstitutiveMatrix<<std::endl;

    KRATOS_CATCH( "" )
  }
  
  //************************************************************************************
  //************************************************************************************
  void SmallDisplacementBeamElement::CalculateStressResultants(ElementDataPointerType& pVariables,
							       const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    Vector NodalDofs;
    this->GetValuesVector(NodalDofs,0);
    
    pVariables->StrainVector = prod( pVariables->B, NodalDofs );

    Matrix ConstitutiveMatrix;
    this->CalculateMaterialConstitutiveMatrix(ConstitutiveMatrix, pVariables);

    //Reference Stress Vector
    pVariables->StressVector = prod( ConstitutiveMatrix, pVariables->StrainVector );


    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
							 ElementDataPointerType& pVariables,
							 double& rIntegrationWeight)
  {
    KRATOS_TRY

    MatrixType LocalStiffnessMatrix = rLeftHandSideMatrix;
    this->CalculateLocalStiffnessMatrix(LocalStiffnessMatrix,pVariables);

    std::cout<<" Direct Stiffness Matrix  "<<LocalStiffnessMatrix<<std::endl;
    
    //contributions to stiffness matrix calculated on the reference config
    noalias( rLeftHandSideMatrix ) += prod( trans( pVariables->B ),  rIntegrationWeight * Matrix( prod( pVariables->ConstitutiveMatrix, pVariables->B ) ) ); //to be optimized to remove the temporary

    std::cout<<" Stiffness Matrix Timoshenko "<<rLeftHandSideMatrix<<std::endl;  
    
    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
								   ElementDataPointerType & pVariables,
								   double& rIntegrationWeight)
  {
    KRATOS_TRY

    //contributions to the internal force vector calculated on the reference config
    noalias( rRightHandSideVector ) += rIntegrationWeight * prod( trans( pVariables->B ),  pVariables->StressVector); //to be optimized to remove the temporary
    
    std::cout<<" Internal Forces Timoshenko "<<rRightHandSideVector<<std::endl; 
      
    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::CalculateLocalStiffnessMatrix(Matrix& LocalStiffnessMatrix, ElementDataPointerType& pVariables)
  {
    KRATOS_TRY

    const double PoissonCoefficient = GetProperties()[POISSON_RATIO];
    const double YoungModulus       = GetProperties()[YOUNG_MODULUS];
    const double ShearModulus       = YoungModulus*0.5/(1.0 + PoissonCoefficient);

    const double L   = pVariables->Length;
    const double LL  = pVariables->Length * pVariables->Length;
    const double LLL = pVariables->Length * pVariables->Length * pVariables->Length;

    double const EA  = pVariables->Section.Area * YoungModulus;
    double const EIz = pVariables->Section.Inertia_z  * YoungModulus;
    double const EIy = pVariables->Section.Inertia_y  * YoungModulus;
    double const JG  = pVariables->Section.Polar_Inertia * ShearModulus;

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

  void SmallDisplacementBeamElement::MapLocalToGlobal(ElementDataPointerType& pVariables, MatrixType& rMatrix)
  {
    KRATOS_TRY
      
    const SizeType& dimension = this->Dimension();
    
    // calculate rotation matrix from spatial frame to material frame
    this->CalculateTransformationMatrix(pVariables->CurrentRotationMatrix);

    if( dimension == 2 )
      BeamMathUtilsType::MapLocalToGlobal2D(pVariables->CurrentRotationMatrix, rMatrix);
    else
      BeamMathUtilsType::MapLocalToGlobal3D(pVariables->CurrentRotationMatrix, rMatrix);
    
    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void SmallDisplacementBeamElement::MapLocalToGlobal(ElementDataPointerType& pVariables, VectorType& rVector)
  {
    KRATOS_TRY

    const SizeType& dimension = this->Dimension();

    // calculate rotation matrix from spatial frame to material frame
    this->CalculateTransformationMatrix(pVariables->CurrentRotationMatrix);
   
    if( dimension == 2 )
      BeamMathUtilsType::MapLocalToGlobal2D(pVariables->CurrentRotationMatrix, rVector);
    else
      BeamMathUtilsType::MapLocalToGlobal3D(pVariables->CurrentRotationMatrix, rVector);
    
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


