//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CONSTITUTIVE_MODEL_UTILITIES )
#define  KRATOS_CONSTITUTIVE_MODEL_UTILITIES

// System includes
#include <cmath>
#include <limits>
#include <algorithm>

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "geometries/geometry.h"
#include "includes/node.h"

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
   */

  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) ConstitutiveModelUtilities
  {
  public:

    ///@name Type Definitions
    ///@{
    typedef BoundedMatrix<double,3,3>    MatrixType;
    typedef array_1d<double,6>            VectorType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConstitutiveModelUtilities() {};

    /// Destructor.
    ~ConstitutiveModelUtilities() {};

    ///@}
    ///@name Operators
    ///@{


    /**
     * Takes a matrix 2x2 and transforms it to a 3x3 adding a 3rd row and a 3rd column with a 1 in the diagonal
     * if the matrix passed is 3D is does nothing
     * if the matrix passed is bigger or smaller throws an error
     * @param rF : the DeformationGradientF in 2D / 3D
     * @param rF3D : the DeformationGradientF in 3D
     */
    static inline MatrixType& DeformationGradientTo3D(const MatrixType& rF, MatrixType& rF3D)
    {
      KRATOS_TRY

      for(unsigned int i=0; i<rF.size1(); i++)
	for(unsigned int j=0; j<rF.size2(); j++)
	  rF3D(i,j) = rF(i,j);

      if (rF.size1() == 2 && rF.size2() == 2)
	{
	  rF3D( 0 , 2 ) = 0.0;
	  rF3D( 1 , 2 ) = 0.0;

	  rF3D( 2 , 0 ) = 0.0;
	  rF3D( 2 , 1 ) = 0.0;

	  rF3D( 2 , 2 ) = 1.0;
	}
      else if(rF.size1() != 3 && rF.size2() != 3)
	{
	  KRATOS_ERROR << "Matrix Dimensions are not correct" << std::endl;
	}

      return rF3D;

      KRATOS_CATCH(" ")
    }

    /**
     * Takes a matrix 2x2 and transforms it to a 3x3 adding a 3rd row and a 3rd column with a 0 in the diagonal
     * if the matrix passed is 3D is does nothing
     * if the matrix passed is bigger or smaller throws an error
     * @param rL : the VelocityGradient in 2D / 3D
     * @param rL3D : the VelocityGradient in 3D
     */
    static inline MatrixType& VelocityGradientTo3D(const MatrixType& rL, MatrixType& rL3D)
    {
      KRATOS_TRY

      for(unsigned int i=0; i<rL.size1(); i++)
	for(unsigned int j=0; j<rL.size2(); j++)
	  rL3D(i,j) = rL(i,j);

      if (rL.size1() == 2 && rL.size2() == 2)
	{
	  rL3D( 0 , 2 ) = 0.0;
	  rL3D( 1 , 2 ) = 0.0;

	  rL3D( 2 , 0 ) = 0.0;
	  rL3D( 2 , 1 ) = 0.0;

	  rL3D( 2 , 2 ) = 0.0;
	}
      else if(rL.size1() != 3 && rL.size2() != 3)
	{
	  KRATOS_ERROR << "Matrix Dimensions are not correct" << std::endl;
	}

      return rL3D;

      KRATOS_CATCH(" ")
    }


    /**
     * Computes the RightCauchyGreen (C=FT*F) given the DeformationGradientF
     * @param rDeformationGradientF input matrix
     * @param rRightCauchyGreen output matrix
     * correct dimensions for the input/output is needed
     */
    static inline void CalculateRightCauchyGreen( const MatrixType& rDeformationGradientF,
						  MatrixType& rRightCauchyGreen )
    {
      noalias( rRightCauchyGreen ) = prod( trans(rDeformationGradientF), rDeformationGradientF );
    }


    /**
     * Computes the LeftCauchyGreen (b=F*FT) given the DeformationGradientF
     * @param rDeformationGradientF input matrix
     * @param rRightCauchyGreen output matrix
     * correct dimensions for the input/output is needed
     */
    static inline void CalculateLeftCauchyGreen( const Matrix & rDeformationGradientF,
						 Matrix& rLeftCauchyGreen )
    {
      noalias( rLeftCauchyGreen ) = prod( rDeformationGradientF, trans(rDeformationGradientF) );
    }


    /**
     * Computes the InverseLeftCauchyGreen (invb=invFT*invF) given the DeformationGradientF
     * @param rDeformationGradientF input matrix
     * @param rRightCauchyGreen output matrix
     * correct dimensions for the input/output is needed
     */
    static inline void CalculateInverseLeftCauchyGreen( const Matrix & rDeformationGradientF,
							Matrix& rInverseLeftCauchyGreen )
    {
      Matrix LeftCauchyGreen;
      CalculateLeftCauchyGreen( rDeformationGradientF, LeftCauchyGreen );

      //calculating the inverse
      double det_b=0;
      MathUtils<double>::InvertMatrix( LeftCauchyGreen, rInverseLeftCauchyGreen, det_b);
    }

        /**
     * Computes the GreenLagrangeStrain E= 0.5*(C-1) given the RightCauchyGreenTensor
     * @param rRightCauchyGreen input matrix 3x3
     * @param rStrainVector output vector
     * correct dimensions for the input/output is needed
     */
    static inline void RightCauchyToGreenLagrangeStrain( const MatrixType& rRightCauchyGreen,
							 MatrixType& rStrainMatrix )
    {

      rStrainMatrix = rRightCauchyGreen;

      rStrainMatrix(0,0) -= 1.0;
      rStrainMatrix(1,1) -= 1.0;
      rStrainMatrix(2,2) -= 1.0;

      rStrainMatrix *= 0.5;
    }

    /**
     * Computes the AlmansiStrain e = 0.5*(1-inv(b)) given the InverseCauchyGreenTensor
     * @param rInverseLeftCauchyGreen input matrix 3x3
     * @param rStrainVector output vector
     * correct dimensions for the input/output is needed
     */
    static inline void LeftCauchyToAlmansiStrain( const MatrixType& rLeftCauchyGreen,
              MatrixType& rStrainMatrix )
    {
      double I3 = 0;

      InvertMatrix3( rLeftCauchyGreen, rStrainMatrix, I3 );

      rStrainMatrix *= (-0.5);

      rStrainMatrix(0,0) += 0.5;
      rStrainMatrix(1,1) += 0.5;
      rStrainMatrix(2,2) += 0.5;
   }

    /**
     * Computes the GreenLagrangeStrain E= 0.5*(C-1) given the RightCauchyGreenTensor
     * @param rRightCauchyGreen input matrix 3x3
     * @param rStrainVector output vector
     * correct dimensions for the input/output is needed
     */
    static inline void RightCauchyToGreenLagrangeStrain( const Matrix& rRightCauchyGreen,
							 Vector& rStrainVector )
    {

      if( rStrainVector.size() == 6 ){

	rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
	rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
	rStrainVector[2] = 0.5 * ( rRightCauchyGreen( 2, 2 ) - 1.00 );
	rStrainVector[3] = rRightCauchyGreen( 0, 1 ); // xy
	rStrainVector[4] = rRightCauchyGreen( 1, 2 ); // yz
	rStrainVector[5] = rRightCauchyGreen( 0, 2 ); // xz

      }
      else if( rStrainVector.size() == 4 ){

	rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
        rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
        rStrainVector[2] = 0.5 * ( rRightCauchyGreen( 2, 2 ) - 1.00 );
        rStrainVector[3] = rRightCauchyGreen( 0, 1 ); // xy
      }
      else if( rStrainVector.size() == 3){

	rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
	rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
	rStrainVector[2] = rRightCauchyGreen( 0, 1 ); // xy

      }
      else{
        KRATOS_ERROR << "Strain Vector dimensions are not correct" << std::endl;
      }

    }

    /**
     * Computes the GreenLagrangeStrain E= 0.5*(FT*F-1) given the DeformationGradientF
     * @param rDeformationGradientF input matrix 3x3
     * @param rStrainVector output vector
     * correct dimensions for the input/output is needed
     */
    static inline void CalculateGreenLagrangeStrain( const Matrix& rDeformationGradientF,
						     Vector& rStrainVector )
    {

      //E= 0.5*(FT*F-1) or E = 0.5*(C-1)
      MatrixType RightCauchyGreen;
      CalculateRightCauchyGreen( rDeformationGradientF, RightCauchyGreen );


      RightCauchyToGreenLagrangeStrain(RightCauchyGreen,rStrainVector);

      // Matrix StrainMatrix(3,3);
      // CalculateGreenLagrangeStrain(rDeformationGradientF, StrainMatrix);
      // noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector( StrainMatrix, rStrainVector.size() );

    }



    /**
     * Computes the GreenLagrangeStrain E= 0.5*(FT*F-1) given the DeformationGradientF
     * @param rDeformationGradientF input matrix 3x3
     * @param rStrainMatrix output matrix
     */
    static inline void CalculateGreenLagrangeStrain( const MatrixType & rDeformationGradientF,
						     MatrixType& rStrainMatrix )
    {

      CalculateRightCauchyGreen( rDeformationGradientF, rStrainMatrix );

      rStrainMatrix(0,0) -= 1;
      rStrainMatrix(1,1) -= 1;
      rStrainMatrix(2,2) -= 1;

      rStrainMatrix *= 0.5;

    }

    /**
     * Computes the AlmansiStrain e = 0.5*(1-inv(b)) given the InverseCauchyGreenTensor
     * @param rInverseLeftCauchyGreen input matrix 3x3
     * @param rStrainVector output vector
     * correct dimensions for the input/output is needed
     */
    static inline void LeftCauchyToAlmansiStrain( const Matrix& rLeftCauchyGreen,
						  Vector& rStrainVector )
    {
      double I3 = 0;
      Matrix InverseLeftCauchyGreen;
      MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, I3 );
      InverseLeftCauchyToAlmansiStrain( InverseLeftCauchyGreen, rStrainVector );
    }

    /**
     * Computes the AlmansiStrain e = 0.5*(1-inv(b)) given the InverseCauchyGreenTensor
     * @param rInverseLeftCauchyGreen input matrix 3x3
     * @param rStrainVector output vector
     * correct dimensions for the input/output is needed
     */
    static inline void InverseLeftCauchyToAlmansiStrain( const Matrix & rInverseLeftCauchyGreen,
							 Vector& rStrainVector )
    {

      if( rStrainVector.size() == 6 ){

	rStrainVector[0] = 0.5 * (  1.00 - rInverseLeftCauchyGreen( 0, 0 ) );
	rStrainVector[1] = 0.5 * (  1.00 - rInverseLeftCauchyGreen( 1, 1 ) );
	rStrainVector[2] = 0.5 * (  1.00 - rInverseLeftCauchyGreen( 2, 2 ) );
	rStrainVector[3] = - rInverseLeftCauchyGreen( 0, 1 ); // xy
	rStrainVector[4] = - rInverseLeftCauchyGreen( 1, 2 ); // yz
	rStrainVector[5] = - rInverseLeftCauchyGreen( 0, 2 ); // xz

      }
      else if( rStrainVector.size() == 4 ){

	rStrainVector[0] = 0.5 * (  1.00 - rInverseLeftCauchyGreen( 0, 0 ) );
	rStrainVector[1] = 0.5 * (  1.00 - rInverseLeftCauchyGreen( 1, 1 ) );
	rStrainVector[2] = 0.5 * (  1.00 - rInverseLeftCauchyGreen( 2, 2 ) );
	rStrainVector[3] = - rInverseLeftCauchyGreen( 0, 1 ); // xy

      }
      else if( rStrainVector.size() == 3){

	rStrainVector[0] = 0.5 * (  1.00 - rInverseLeftCauchyGreen( 0, 0 ) );
	rStrainVector[1] = 0.5 * (  1.00 - rInverseLeftCauchyGreen( 1, 1 ) );
	rStrainVector[2] = - rInverseLeftCauchyGreen( 0, 1 ); // xy

      }
      else{
        KRATOS_ERROR << "Strain Vector dimensions are not correct" << std::endl;
      }

      // Matrix StrainMatrix(3,3);
      // CalculateAlmansiStrain(rDeformationGradientF, StrainMatrix);
      // noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector( StrainMatrix, rStrainVector.size() );

    }


    /**
     * Computes the AlmansiStrain e = 0.5*(1-invFT*invF) given the DeformationGradientF
     * @param rDeformationGradientF input matrix 3x3
     * @param rStrainVector output vector
     * correct dimensions for the input/output is needed
     */
    static inline void CalculateAlmansiStrain( const Matrix & rDeformationGradientF,
					       Vector& rStrainVector )
    {

      // e = 0.5*(1-invFT*invF) or e = 0.5*(1-inv(b))
      Matrix InverseLeftCauchyGreen(3,3);
      CalculateInverseLeftCauchyGreen( rDeformationGradientF, InverseLeftCauchyGreen );

      InverseLeftCauchyToAlmansiStrain( InverseLeftCauchyGreen, rStrainVector );

      // Matrix StrainMatrix(3,3);
      // CalculateAlmansiStrain(rDeformationGradientF, StrainMatrix);
      // noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector( StrainMatrix, rStrainVector.size() );

    }


    /**
     * Computes the AlmansiStrain e = 0.5*(1-invFT*invF) given the DeformationGradientF
     * @param rDeformationGradientF input matrix 3x3
     * @param rStrainMatrix output matrix
     */
    static inline void CalculateAlmansiStrain( const Matrix & rDeformationGradientF,
					       Matrix& rStrainMatrix )
    {

      CalculateInverseLeftCauchyGreen( rDeformationGradientF, rStrainMatrix );

      rStrainMatrix(0,0) -= 1;
      rStrainMatrix(1,1) -= 1;
      rStrainMatrix(2,2) -= 1;

      rStrainMatrix *= -0.5;

    }


    /**
     * Transforms a given 3D Constitutive Tensor to VoigtSize Constitutive Matrix:
     * in the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (VoigtSize*VoightSize) Vector
     * @param rTensor the given second order tensor in matrix form
     * @param rMatrix the corresponding second order tensor in voigt size matrix form
     */

    static inline Matrix& ConstitutiveTensorToMatrix(const BoundedMatrix<double,6,6>& rTensor, Matrix& rMatrix)
    {
        KRATOS_TRY;

	if( rMatrix.size1() == 6 ){

	  rMatrix = rTensor;

	}
	else if( rMatrix.size1() == 4 ){

	  for(unsigned int i=0; i<3; i++)
	    {
	      for(unsigned int j=0; j<3; j++)
		{
		  rMatrix(i,j) = rTensor(i,j);
		}
	    }

	  rMatrix(3,3) = rTensor(3,3);

	}
	else if( rMatrix.size1() == 3){

	  for(unsigned int i=0; i<2; i++)
	    {
	      for(unsigned int j=0; j<2; j++)
		{
		  rMatrix(i,j) = rTensor(i,j);
		}
	    }

	  rMatrix(2,2) = rTensor(3,3);

	}
	else{
	  KRATOS_ERROR << "Constitutive Matrix dimensions are not correct" << std::endl;
	}

	return rMatrix;

        KRATOS_CATCH("");
     }


    /**
     * Transforms a given Vector to a non symmetric 3D Tensor:
     * in the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (9*1) Vector
     * in the 2D case: from a second order tensor (3*3) Matrix  to a corresponing (4*1) Vector
     * @param rTensor the given symmetric second order stress tensor
     * @return the corresponding stress tensor in vector form
     */
    static inline MatrixType& VectorToTensor(const Vector& rVector, MatrixType& rTensor)
    {
        KRATOS_TRY;

        // vector2D = [ a00, a11, a01, a10, a21, a02 ]
        // vector3D = [ a00, a11, a22, a01, a12, a20, a10, a21, a02 ]

	if (rVector.size() == 4)
        {
   	    rTensor(0,0) = rVector[0];
	    rTensor(0,1) = rVector[2];
	    rTensor(0,2) = 0.0;

	    rTensor(1,0) = rVector[3];
	    rTensor(1,1) = rVector[1];
	    rTensor(1,2) = 0.0;

	    rTensor(2,0) = 0.0;
	    rTensor(2,1) = 0.0;
	    rTensor(2,2) = 0.0;
        }
        else if (rVector.size() == 9)
        {
	    rTensor(0,0) = rVector[0];
	    rTensor(0,1) = rVector[3];
	    rTensor(0,2) = rVector[8];

	    rTensor(1,0) = rVector[6];
	    rTensor(1,1) = rVector[1];
	    rTensor(1,2) = rVector[4];

	    rTensor(2,0) = rVector[5];
	    rTensor(2,1) = rVector[7];
	    rTensor(2,2) = rVector[2];
        }
        else{
          KRATOS_ERROR << " VectorToTensor transform Vector Size not correct : " << rVector.size() <<std::endl;
        }

        return rTensor;

        KRATOS_CATCH("");
    }

    /**
     * Transforms a given non symmetric Tensor to a Vector:
     * in the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (9*1) Vector
     * in the 2D case: from a second order tensor (2*2) Matrix  to a corresponing (4*1) Vector
     * @param rTensor the given symmetric second order stress tensor
     * @return the corresponding stress tensor in vector form
     */

    static inline Vector& TensorToVector(const MatrixType& rTensor, Vector& rVector)
    {
        KRATOS_TRY;

        // vector2D = [ a00, a11, a01, a10, a21, a02 ]
        // vector3D = [ a00, a11, a22, a01, a12, a20, a10, a21, a02 ]

        if (rVector.size() == 4)
        {
	    rVector[0] = rTensor(0,0);
            rVector[1] = rTensor(1,1);
            rVector[2] = rTensor(0,1);
            rVector[3] = rTensor(1,0);
        }
        else if (rVector.size() == 9)
        {
            rVector[0] = rTensor(0,0);
            rVector[1] = rTensor(1,1);
            rVector[2] = rTensor(2,2);
            rVector[3] = rTensor(0,1);
            rVector[4] = rTensor(1,2);
            rVector[5] = rTensor(2,0);
            rVector[6] = rTensor(1,0);
            rVector[7] = rTensor(2,1);
            rVector[8] = rTensor(0,2);
        }
        else{
          KRATOS_ERROR << " TensorToVector transform Vector Size not correct : " << rVector.size() <<std::endl;
        }

        return rVector;

        KRATOS_CATCH("");
    }

    /**
     * Transforms a given 3D symmetric Tensor from Voigt notation to Matrix notation
     * in the 3D case: from a second order tensor (6*1) Vector to a corresponing (3*3) Matrix
     * @param rVector the given symmetric second order tensor in vector form
     * @param rMatrix the corresponding second order tensor in matrix form
     */
    static inline MatrixType& VectorToSymmetricTensor(const array_1d<double,6>& rVector, MatrixType& rMatrix)
    {
        KRATOS_TRY;

	rMatrix(0,0) = rVector[0];
	rMatrix(0,1) = rVector[3];
	rMatrix(0,2) = rVector[5];

	rMatrix(1,0) = rVector[3];
	rMatrix(1,1) = rVector[1];
	rMatrix(1,2) = rVector[4];

	rMatrix(2,0) = rVector[5];
	rMatrix(2,1) = rVector[4];
	rMatrix(2,2) = rVector[2];

        return rMatrix;

        KRATOS_CATCH("");
    }


    /**
     * Transforms a given 3D symmetric Tensor to Voigt Notation:
     * in the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     * @param rMatrix the given symmetric second order tensor in matrix form
     * @param rVector the corresponding second order tensor in vector form
     */

    static inline void SymmetricTensorToVector(const MatrixType& rMatrix, array_1d<double,6>& rVector)
    {
        KRATOS_TRY;

	rVector[0]= rMatrix(0,0);
	rVector[1]= rMatrix(1,1);
	rVector[2]= rMatrix(2,2);

	rVector[3]= rMatrix(0,1);
	rVector[4]= rMatrix(1,2);
	rVector[5]= rMatrix(2,0);

        KRATOS_CATCH("");
     }


    /**
     * Transforms a given 3D symmetric Tensor from Voigt notation to Matrix notation
     * in the 3D case: from a second order tensor (6*1) Vector to a corresponing (3*3) Matrix
     * @param rVector the given symmetric second order tensor in vector form
     * @param rMatrix the corresponding second order tensor in matrix form
     */
    static inline MatrixType& StrainVectorToTensor(const array_1d<double,6>& rVector, MatrixType& rMatrix)
    {
        KRATOS_TRY;

	rMatrix(0,0) = rVector[0];
	rMatrix(0,1) = 0.5*rVector[3];
	rMatrix(0,2) = 0.5*rVector[5];
	rMatrix(1,0) = 0.5*rVector[3];
	rMatrix(1,1) = rVector[1];
	rMatrix(1,2) = 0.5*rVector[4];
	rMatrix(2,0) = 0.5*rVector[5];
	rMatrix(2,1) = 0.5*rVector[4];
	rMatrix(2,2) = rVector[2];

        return rMatrix;

        KRATOS_CATCH("");
    }


    /**
     * Transforms a given 3D symmetric Tensor to Voigt Notation:
     * in the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     * @param rMatrix the given symmetric second order tensor in matrix form
     * @param rVector the corresponding second order tensor in vector form
     */

    static inline void StrainTensorToVector(const MatrixType& rMatrix, array_1d<double,6>& rVector)
    {
        KRATOS_TRY;

	rVector[0]= rMatrix(0,0);
	rVector[1]= rMatrix(1,1);
	rVector[2]= rMatrix(2,2);
	rVector[3]= 2.0*rMatrix(0,1);
	rVector[4]= 2.0*rMatrix(1,2);
	rVector[5]= 2.0*rMatrix(0,2);

        KRATOS_CATCH("");
     }

    /**
     * Transforms a given symmetric Strain Tensor to Voigt Notation:
     * in the 3D  case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     * in the 2Da case: from a second order tensor (3*3) Matrix  to a corresponing (4*1) Vector
     * in the 2D  case: from a second order tensor (3*3) Matrix  to a corresponing (3*1) Vector
     * @param rStrainTensor the given symmetric second order stress tensor
     * @return the corresponding stress tensor in vector form
     */

    static inline MatrixType& StrainVectorToTensor(const Vector& rStrainVector, MatrixType& rStrainTensor)
    {
        KRATOS_TRY;

	if (rStrainVector.size() == 3)
        {
   	    rStrainTensor(0,0) = rStrainVector[0];
	    rStrainTensor(0,1) = 0.5*rStrainVector[2];
	    rStrainTensor(0,2) = 0.0;

	    rStrainTensor(1,0) = 0.5*rStrainVector[2];
	    rStrainTensor(1,1) = rStrainVector[1];
	    rStrainTensor(1,2) = 0.0;

	    rStrainTensor(2,0) = 0.0;
	    rStrainTensor(2,1) = 0.0;
	    rStrainTensor(2,2) = 0.0;
        }
        else if (rStrainVector.size() == 4)
        {
   	    rStrainTensor(0,0) = rStrainVector[0];
	    rStrainTensor(0,1) = 0.5*rStrainVector[3];
	    rStrainTensor(0,2) = 0.0;

	    rStrainTensor(1,0) = 0.5*rStrainVector[3];
	    rStrainTensor(1,1) = rStrainVector[1];
	    rStrainTensor(1,2) = 0.0;

	    rStrainTensor(2,0) = 0.0;
	    rStrainTensor(2,1) = 0.0;
	    rStrainTensor(2,2) = rStrainVector[2];
        }
        else if (rStrainVector.size() == 6)
        {
	    rStrainTensor(0,0) = rStrainVector[0];
	    rStrainTensor(0,1) = 0.5*rStrainVector[3];
	    rStrainTensor(0,2) = 0.5*rStrainVector[5];

	    rStrainTensor(1,0) = 0.5*rStrainVector[3];
	    rStrainTensor(1,1) = rStrainVector[1];
	    rStrainTensor(1,2) = 0.5*rStrainVector[4];

	    rStrainTensor(2,0) = 0.5*rStrainVector[5];
	    rStrainTensor(2,1) = 0.5*rStrainVector[4];
	    rStrainTensor(2,2) = rStrainVector[2];
        }
        else{
          KRATOS_ERROR << "Unexpected voigt size: " << rStrainVector.size() << std::endl;
        }

        return rStrainTensor;

        KRATOS_CATCH("");
    }

    /**
     * Transforms a given symmetric Strain Tensor to Voigt Notation:
     * in the 3D  case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     * in the 2Da case: from a second order tensor (3*3) Matrix  to a corresponing (4*1) Vector
     * in the 2D  case: from a second order tensor (3*3) Matrix  to a corresponing (3*1) Vector
     * @param rStrainTensor the given symmetric second order stress tensor
     * @return the corresponding stress tensor in vector form
     */

    static inline Vector& StrainTensorToVector(const MatrixType& rStrainTensor, Vector& rStrainVector)
    {
        KRATOS_TRY;

	if (rStrainVector.size() == 3)
        {
          rStrainVector[0] = rStrainTensor(0,0);
          rStrainVector[1] = rStrainTensor(1,1);
          rStrainVector[2] = 2.0*rStrainTensor(0,1);
        }
        else if (rStrainVector.size() == 4)
        {
          rStrainVector[0] = rStrainTensor(0,0);
          rStrainVector[1] = rStrainTensor(1,1);
          rStrainVector[2] = rStrainTensor(2,2);
          rStrainVector[3] = 2.0*rStrainTensor(0,1);
        }
        else if (rStrainVector.size() == 6)
        {
          rStrainVector[0] = rStrainTensor(0,0);
          rStrainVector[1] = rStrainTensor(1,1);
          rStrainVector[2] = rStrainTensor(2,2);
          rStrainVector[3] = 2.0*rStrainTensor(0,1);
          rStrainVector[4] = 2.0*rStrainTensor(1,2);
          rStrainVector[5] = 2.0*rStrainTensor(0,2);
        }
        else{
          KRATOS_ERROR << "Unexpected voigt size: " << rStrainVector.size() << std::endl;
        }


        return rStrainVector;


        KRATOS_CATCH("");
    }

    /**
     * Transforms a given symmetric Stress Tensor to Voigt Notation:
     * in the 3D  case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     * in the 2Da case: from a second order tensor (3*3) Matrix  to a corresponing (4*1) Vector
     * in the 2D  case: from a second order tensor (3*3) Matrix  to a corresponing (3*1) Vector
     * @param rStressTensor the given symmetric second order stress tensor
     * @return the corresponding stress tensor in vector form
     */

    static inline MatrixType& StressVectorToTensor(const Vector& rStressVector, MatrixType& rStressTensor)
    {
        KRATOS_TRY;

	if (rStressVector.size() == 3)
        {
          rStressTensor(0,0) = rStressVector[0];
          rStressTensor(0,1) = rStressVector[2];
          rStressTensor(0,2) = 0.0;

          rStressTensor(1,0) = rStressVector[2];
          rStressTensor(1,1) = rStressVector[1];
          rStressTensor(1,2) = 0.0;

          rStressTensor(2,0) = 0.0;
          rStressTensor(2,1) = 0.0;
          rStressTensor(2,2) = 0.0;
        }
        else if (rStressVector.size() == 4)
        {
          rStressTensor(0,0) = rStressVector[0];
          rStressTensor(0,1) = rStressVector[3];
          rStressTensor(0,2) = 0.0;

          rStressTensor(1,0) = rStressVector[3];
          rStressTensor(1,1) = rStressVector[1];
          rStressTensor(1,2) = 0.0;

          rStressTensor(2,0) = 0.0;
          rStressTensor(2,1) = 0.0;
          rStressTensor(2,2) = rStressVector[2];
        }
        else if (rStressVector.size() == 6)
        {
          rStressTensor(0,0) = rStressVector[0];
          rStressTensor(0,1) = rStressVector[3];
          rStressTensor(0,2) = rStressVector[5];

          rStressTensor(1,0) = rStressVector[3];
          rStressTensor(1,1) = rStressVector[1];
          rStressTensor(1,2) = rStressVector[4];

          rStressTensor(2,0) = rStressVector[5];
          rStressTensor(2,1) = rStressVector[4];
          rStressTensor(2,2) = rStressVector[2];
        }
        else{
          KRATOS_ERROR << "Unexpected voigt size: " << rStressVector.size() << std::endl;
        }

        return rStressTensor;

        KRATOS_CATCH("");
    }

    /**
     * Transforms a given symmetric Stress Tensor to Voigt Notation:
     * in the 3D  case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     * in the 2Da case: from a second order tensor (3*3) Matrix  to a corresponing (4*1) Vector
     * in the 2D  case: from a second order tensor (3*3) Matrix  to a corresponing (3*1) Vector
     * @param rStressTensor the given symmetric second order stress tensor
     * @return the corresponding stress tensor in vector form
     */

    static inline Vector& StressTensorToVector(const MatrixType& rStressTensor, Vector& rStressVector)
    {
        KRATOS_TRY;

        if (rStressVector.size() == 3)
        {
	    rStressVector[0] = rStressTensor(0,0);
            rStressVector[1] = rStressTensor(1,1);
            rStressVector[2] = rStressTensor(0,1);
        }
        else if (rStressVector.size() == 4)
        {
	    rStressVector[0] = rStressTensor(0,0);
            rStressVector[1] = rStressTensor(1,1);
            rStressVector[2] = rStressTensor(2,2);
            rStressVector[3] = rStressTensor(0,1);
        }
        else if (rStressVector.size() == 6)
        {
            rStressVector[0] = rStressTensor(0,0);
            rStressVector[1] = rStressTensor(1,1);
            rStressVector[2] = rStressTensor(2,2);
            rStressVector[3] = rStressTensor(0,1);
            rStressVector[4] = rStressTensor(1,2);
            rStressVector[5] = rStressTensor(0,2);
        }

        return rStressVector;

        KRATOS_CATCH("");
    }


    /**
     * Computes the Stress Norm
     * @param rStressMatrix input matrix 3x3
     * @param rStressNorm output double
     */
    static inline double& CalculateStressNorm( const MatrixType& rStressMatrix,
					       double& rStressNorm )
    {

      rStressNorm =  sqrt((rStressMatrix(0,0)*rStressMatrix(0,0))+(rStressMatrix(1,1)*rStressMatrix(1,1))+(rStressMatrix(2,2)*rStressMatrix(2,2))+
			  (rStressMatrix(0,1)*rStressMatrix(0,1))+(rStressMatrix(0,2)*rStressMatrix(0,2))+(rStressMatrix(1,2)*rStressMatrix(1,2))+
			  (rStressMatrix(1,0)*rStressMatrix(1,0))+(rStressMatrix(2,0)*rStressMatrix(2,0))+(rStressMatrix(2,1)*rStressMatrix(2,1)));


      return rStressNorm;

    }


    /**
     * Computes the Characteristic Size of a Geometry
     * @param rDomainGeomtry input ElementGeometry
     * @param rCharacteristicSize output Size of the Geometry
     */
    static inline double& CalculateCharacteristicSize( const Geometry<Node<3> >& rDomainGeometry,
						       double& rCharacteristicSize )
    {

      if( rDomainGeometry.WorkingSpaceDimension() == 2)
	//rCharacteristicSize is the diameter of a circle with the same area as the element
	rCharacteristicSize = sqrt(4.0*rDomainGeometry.Area()/Globals::Pi);

      if( rDomainGeometry.WorkingSpaceDimension() == 3)
	//rCharacteristicSize is the diameter of a sphere with the same volume as the element
	rCharacteristicSize = pow((6.0*rDomainGeometry.Volume()/Globals::Pi),0.33333333333333);

      return rCharacteristicSize;

    }


    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculates perturbations in each direction of a given vector.
     * @param InputVector the given vector used to obtain the perturbations.
     */

    static inline void ComputePerturbationVector( Vector& rPerturbationVector, const Vector& InputVector )
    {
        const unsigned int VSize = InputVector.size();
        if(rPerturbationVector.size() != VSize)
            rPerturbationVector.resize(VSize,false);

        const double MinTol = 1.0e-10;
        const double MaxTol = 1.0e-5;

        //Maximum and minimum vector components
        double max_component = fabs(InputVector[0]) , min_component = fabs(InputVector[0]);

        for( unsigned int i=1; i<VSize; i++ )
        {
            if( fabs(InputVector[i]) < min_component )
            {
                min_component = fabs(InputVector[i]);
            }
            else if( fabs(InputVector[i]) > max_component )
            {
                max_component = fabs(InputVector[i]);
            }
        }

        double aux = min_component*MaxTol;

        if( aux < (max_component*MinTol) )
        {
            aux = max_component*MinTol;
        }

        //PerturbationVector
        for( unsigned int i=0; i<VSize; i++ )
        {
            if( fabs(InputVector[i]) > 1.0e-20 ) // different from zero
            {
                rPerturbationVector[i] = InputVector[i]*MaxTol;
            }
            else if( InputVector[i] >= 0.0 )
            {
                rPerturbationVector[i] = aux;
            }
            else
            {
                rPerturbationVector[i] = -aux;
            }
        }
    }


    /**
     * calculates the eigenvectiors using a direct method.
     * @param A the given square matrix the eigenvalues are to be calculated.
     * WARNING only valid symmetric 3*3 Matrices
     */


    static inline Vector EigenValuesDirectMethod(const Matrix& A)
    {
        // Given a real symmetric 3x3 matrix A, compute the eigenvalues
        int dim= A.size1();
        Vector Result(dim);
	noalias(Result) = ZeroVector(dim);

        const double p1 = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2);
        if (p1 == 0)
        {//A is diagonal.
            Result[0] = A(0,0);
            Result[1] = A(1,1);
            Result[2] = A(2,2);
            return Result;
        }

        const double q = (A(0,0) + A(1,1) + A(2,2)) / 3.0;
        const double p2 = (A(0,0) - q) * (A(0,0) - q) + (A(1,1) - q) * (A(1,1) - q) + (A(2,2) - q) * (A(2,2) - q) + 2.0 * p1;
        const double p = sqrt(p2 / 6.0);

        Matrix B(3,3);
        const double inv_p = 1.0/p;

        // B = (1 / p) * (A - q * I)  where  I is the identity matrix
        B(0,0) = inv_p * (A(0,0) - q);
        B(1,1) = inv_p * (A(1,1) - q);
        B(2,2) = inv_p * (A(2,2) - q);
        B(0,1) = inv_p * A(0,1);
        B(1,0) = inv_p * A(1,0);
        B(0,2) = inv_p * A(0,2);
        B(2,0) = inv_p * A(2,0);
        B(1,2) = inv_p * A(1,2);
        B(2,1) = inv_p * A(2,1);

        //r = det(B) / 2
        double r = 0.5 * ( B(0,0)*B(1,1)*B(2,2) + B(0,1)*B(1,2)*B(2,0) + B(1,0)*B(2,1)*B(0,2) - B(2,0)*B(1,1)*B(0,2) - B(1,0)*B(0,1)*B(2,2) - B(0,0)*B(2,1)*B(1,2) );

        // In exact arithmetic for a symmetric matrix  -1 <= r <= 1
        // but computation error can leave it slightly outside this range.
        double phi = 0.0;
        if (r <= -1) { phi = Globals::Pi / 3.0; }
        else if (r >= 1) { phi = 0.0; }
        else { phi = acos(r) / 3.0;}

        // the eigenvalues satisfy eig3 <= eig2 <= eig1
        Result[0] = q + 2.0 * p * cos(phi);
        Result[2] = q + 2.0 * p * cos(phi + (2.0*Globals::Pi/3.0));
        Result[1] = 3.0 * q - Result[0] - Result[2];     //% since trace(A) = eig1 + eig2 + eig3

        return Result;
    }


    /**
     * It inverts matrices of order 3 //VERIFIED!!!
     * @param InputMatrix: Is the input matrix (unchanged at output)
     * @return InvertedMatrix: Is the inverse of the input matrix
     * @return InputMatrixDet: Is the determinant of the input matrix
     */

    static void InvertMatrix3(const MatrixType& InputMatrix,
			      MatrixType& InvertedMatrix,
			      double& InputMatrixDet
			      )
    {
        KRATOS_TRY;


        // Filling the inverted matrix with the algebraic complements
        // First column
        InvertedMatrix(0,0) = InputMatrix(1,1)*InputMatrix(2,2) - InputMatrix(1,2)*InputMatrix(2,1);
        InvertedMatrix(1,0) = -InputMatrix(1,0)*InputMatrix(2,2) + InputMatrix(1,2)*InputMatrix(2,0);
        InvertedMatrix(2,0) = InputMatrix(1,0)*InputMatrix(2,1) - InputMatrix(1,1)*InputMatrix(2,0);

        // Second column
        InvertedMatrix(0,1) = -InputMatrix(0,1)*InputMatrix(2,2) + InputMatrix(0,2)*InputMatrix(2,1);
        InvertedMatrix(1,1) = InputMatrix(0,0)*InputMatrix(2,2) - InputMatrix(0,2)*InputMatrix(2,0);
        InvertedMatrix(2,1) = -InputMatrix(0,0)*InputMatrix(2,1) + InputMatrix(0,1)*InputMatrix(2,0);

        // Third column
        InvertedMatrix(0,2) = InputMatrix(0,1)*InputMatrix(1,2) - InputMatrix(0,2)*InputMatrix(1,1);
        InvertedMatrix(1,2) = -InputMatrix(0,0)*InputMatrix(1,2) + InputMatrix(0,2)*InputMatrix(1,0);
        InvertedMatrix(2,2) = InputMatrix(0,0)*InputMatrix(1,1) - InputMatrix(0,1)*InputMatrix(1,0);

        // Calculation of determinant (of the input matrix)
        InputMatrixDet = InputMatrix(0,0)*InvertedMatrix(0,0) + InputMatrix(0,1)*InvertedMatrix(1,0) + InputMatrix(0,2)*InvertedMatrix(2,0);

        // Finalizing the calculation of the inverted matrix
        InvertedMatrix *= ( 1.0 / InputMatrixDet );

        KRATOS_CATCH("")
    }


    /**
     * Calculates perturbations in each direction of a given vector.
     * @param InputVector the given vector used to obtain the perturbations.
     */

    static inline void ComputePerturbationVector( VectorType& rPerturbationVector, const VectorType& InputVector )
    {
        const unsigned int VSize = InputVector.size();
        if(rPerturbationVector.size() != VSize)
            rPerturbationVector.resize(VSize,false);

        const double MinTol = 1.0e-10;
        const double MaxTol = 1.0e-5;

        //Maximum and minimum vector components
        double max_component = fabs(InputVector[0]) , min_component = fabs(InputVector[0]);

        for( unsigned int i=1; i<VSize; i++ )
        {
            if( fabs(InputVector[i]) < min_component )
            {
                min_component = fabs(InputVector[i]);
            }
            else if( fabs(InputVector[i]) > max_component )
            {
                max_component = fabs(InputVector[i]);
            }
        }

        double aux = min_component*MaxTol;

        if( aux < (max_component*MinTol) )
        {
            aux = max_component*MinTol;
        }

        //PerturbationVector
        for( unsigned int i=0; i<VSize; i++ )
        {
            if( fabs(InputVector[i]) > 1.0e-20 ) // different from zero
            {
                rPerturbationVector[i] = InputVector[i]*MaxTol;
            }
            else if( InputVector[i] >= 0.0 )
            {
                rPerturbationVector[i] = aux;
            }
            else
            {
                rPerturbationVector[i] = -aux;
            }
        }
    }


    /**
     * Computes fourth order unit tensor
     * @param rValue output double
     * @param a index for the fourth order tensor
     * @param b index for the fourth order tensor
     * @param c index for the fourth order tensor
     * @param d index for the fourth order tensor
     */
    static inline double& CalculateFourthOrderUnitTensor( double& rValue,
							  const unsigned int& a, const unsigned int& b,
							  const unsigned int& c, const unsigned int& d )
    {
      MatrixType IdentityMatrix;
      noalias(IdentityMatrix) = identity_matrix<double>(3);

      rValue = CalculateFourthOrderUnitTensor(IdentityMatrix,rValue,a,b,c,d);

      return rValue;
    }

    /**
     * Computes fourth order unit tensor
     * @param rIdentityMatrix input tensor identity matrix 3x3
     * @param rValue output double
     * @param a index for the fourth order tensor
     * @param b index for the fourth order tensor
     * @param c index for the fourth order tensor
     * @param d index for the fourth order tensor
     */
    static inline double& CalculateFourthOrderUnitTensor( const MatrixType& rIdentityMatrix, double& rValue,
							  const unsigned int& a, const unsigned int& b,
							  const unsigned int& c, const unsigned int& d )
    {
      rValue = 0.5*(rIdentityMatrix(a,c)*rIdentityMatrix(b,d)+rIdentityMatrix(a,d)*rIdentityMatrix(b,c));

      return rValue;
    }


    /**
     * Computes fourth order unit tensor
     * @param rMatrix input tensor identity matrix 3x3
     * @param rValue output double
     * @param a index for the fourth order tensor
     * @param b index for the fourth order tensor
     * @param c index for the fourth order tensor
     * @param d index for the fourth order tensor
     */
    static inline double& CalculateFourthOrderTensor(const MatrixType& rMatrix,
						     double& rValue,
						     const double& a,
						     const double& b,
						     const double& c,
						     const double& d)
    {
	rValue = 0.5*(rMatrix(a,c)*rMatrix(b,d)+rMatrix(a,d)*rMatrix(b,c));

	return rValue;
    }

    /**
     * Computes the Square Tensor Derivative
     * @param rMatrix input tensor matrix 3x3
     * @param rValue output double
     * @param a index for the fourth order tensor
     * @param b index for the fourth order tensor
     * @param c index for the fourth order tensor
     * @param d index for the fourth order tensor
     */
    static inline double& CalculateSquareTensorDerivative( const MatrixType& rMatrix, double& rValue,
							   const unsigned int& a, const unsigned int& b,
							   const unsigned int& c, const unsigned int& d )
    {
      MatrixType IdentityMatrix;
      noalias(IdentityMatrix) = identity_matrix<double>(3);

      rValue = CalculateSquareTensorDerivative(rMatrix,IdentityMatrix,rValue,a,b,c,d);

      return rValue;
    }

    /**
     * Computes the Square Tensor Derivative
     * @param rMatrix input tensor matrix 3x3
     * @param rIdentityMatrix input tensor identity matrix 3x3
     * @param rValue output double
     * @param a index for the fourth order tensor
     * @param b index for the fourth order tensor
     * @param c index for the fourth order tensor
     * @param d index for the fourth order tensor
     */
    static inline double& CalculateSquareTensorDerivative( const MatrixType& rMatrix, const MatrixType& rIdentityMatrix, double& rValue,
							   const unsigned int& a, const unsigned int& b,
							   const unsigned int& c, const unsigned int& d )
    {
      rValue = 0.5*(rIdentityMatrix(a,c)*rMatrix(d,b)+rIdentityMatrix(a,d)*rMatrix(c,b)+rIdentityMatrix(b,d)*rMatrix(a,c)+rIdentityMatrix(c,b)*rMatrix(a,d));

      return rValue;
    }


    /**
     * Computes the FourthOrder Tensor Product
     * @param rMatrixA input tensor matrix 3x3
     * @param rMatrixB input tensor matrix 3x3
     * @param rValue output double
     * @param a index for the fourth order tensor
     * @param b index for the fourth order tensor
     * @param c index for the fourth order tensor
     * @param d index for the fourth order tensor
     */
    static inline double& CalculateFourthOrderTensorProduct( const MatrixType& rMatrixA, const MatrixType& rMatrixB, double& rValue,
							     const unsigned int& a, const unsigned int& b,
							     const unsigned int& c, const unsigned int& d )
    {
      rValue = rMatrixA(a,b)*rMatrixB(c,d);

      return rValue;
    }

     /**
     * Computes the FourthOrder Tensor Product
     * @param rVectorA input vector 3
     * @param rVectorB input vector 3
     * @param rValue output double
     * @param a index for the fourth order tensor
     * @param b index for the fourth order tensor
     * @param c index for the fourth order tensor
     * @param d index for the fourth order tensor
     */
    static inline double& CalculateFourthOrderTensorProduct( const array_1d<double,3>& rVectorA, const array_1d<double,3>& rVectorB, double& rValue,
							     const unsigned int& a, const unsigned int& b,
							     const unsigned int& c, const unsigned int& d )
    {
      rValue = (rVectorA[a]*rVectorA[b])*(rVectorB[c]*rVectorB[d]);

      return rValue;
    }


     /**
     * Computes the FourthOrder Tensor Product
     * @param rMatrixA input tensor matrix 3x3
     * @param rVectorB input vector 3
     * @param rValue output double
     * @param a index for the fourth order tensor
     * @param b index for the fourth order tensor
     * @param c index for the fourth order tensor
     * @param d index for the fourth order tensor
     */
    static inline double& CalculateFourthOrderTensorProduct( const MatrixType& rMatrixA, const array_1d<double,3>& rVectorB, double& rValue,
							     const unsigned int& a, const unsigned int& b,
							     const unsigned int& c, const unsigned int& d )
    {
      rValue = rMatrixA(a,b)*(rVectorB[c]*rVectorB[d]);

      return rValue;
    }

     /**
     * Computes the FourthOrder Tensor Product
     * @param rVectorA input vector 3
     * @param rMatrixB input tensor matrix 3x3
     * @param rValue output double
     * @param a index for the fourth order tensor
     * @param b index for the fourth order tensor
     * @param c index for the fourth order tensor
     * @param d index for the fourth order tensor
     */
    static inline double& CalculateFourthOrderTensorProduct( const array_1d<double,3>& rVectorA, const MatrixType& rMatrixB, double& rValue,
							     const unsigned int& a, const unsigned int& b,
							     const unsigned int& c, const unsigned int& d )
    {
      rValue = (rVectorA[a]*rVectorA[b])*rMatrixB(c,d);

      return rValue;
    }


     /**
     * Checks if two doubles are equal
     * @param rA input double
     * @param rB input double
     */
    static inline bool AreEqual( const double& rA, const double& rB )
    {
	//different sign
	if( rA*rB < 0 )
	    return false;

	//different order of magnitude
	double absDiff = std::fabs(rA - rB);
	if(absDiff <= std::numeric_limits<double>::epsilon())
	{
	    return true;
	}

	//similar order of magnitude
	double maxAbs  = std::max(std::fabs(rA), std::fabs(rB));
	return (absDiff/maxAbs) < 1E-8;

    }



    /**
     * Methods to transform Constitutive Matrices:
     * @param rConstitutiveMatrix the constitutive matrix
     * @param rF the DeformationGradientF matrix between the configurations
     */

    /**
     * This method performs a pull-back of the constitutive matrix
     */
    static inline void PullBackConstitutiveMatrix( Matrix& rConstitutiveMatrix,
                                                   const Matrix & rF )
    {
      Matrix OriginalConstitutiveMatrix = rConstitutiveMatrix;

      rConstitutiveMatrix.clear();

      Matrix InverseF ( 3, 3 );
      double detF = 0;
      MathUtils<double>::InvertMatrix( rF, InverseF, detF);

      ConstitutiveMatrixTransformation( rConstitutiveMatrix, OriginalConstitutiveMatrix, InverseF );
    }


    /**
     * This method performs a push-forward of the constitutive matrix
     */
    static inline void PushForwardConstitutiveMatrix( Matrix& rConstitutiveMatrix,
                                                      const Matrix & rF )
    {
      Matrix OriginalConstitutiveMatrix = rConstitutiveMatrix;

      rConstitutiveMatrix.clear();

      ConstitutiveMatrixTransformation( rConstitutiveMatrix, OriginalConstitutiveMatrix, rF );
    }


    /**
     * This method performs a pull-back or a push-forward between two constitutive matrices
     */
    static inline void ConstitutiveMatrixTransformation ( Matrix& rConstitutiveMatrix,
                                                          const Matrix& rOriginalConstitutiveMatrix,
                                                          const Matrix & rF )
    {
      unsigned int size = rOriginalConstitutiveMatrix.size1();
      if(  size == 6 )
      {
        const unsigned int IndexVoigt3D6C [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };

        for(unsigned int i=0; i<6; i++)
        {
          for(unsigned int j=0; j<6; j++)
          {
            rConstitutiveMatrix( i, j ) = TransformConstitutiveComponent(rConstitutiveMatrix( i, j ), rOriginalConstitutiveMatrix, rF,
                                                                         IndexVoigt3D6C[i][0], IndexVoigt3D6C[i][1], IndexVoigt3D6C[j][0], IndexVoigt3D6C[j][1]);
          }

        }
      }
      else if( size == 4 )
      {

        const unsigned int IndexVoigt2D4C [4][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1} };

        for(unsigned int i=0; i<4; i++)
        {
          for(unsigned int j=0; j<4; j++)
          {
            rConstitutiveMatrix( i, j ) = TransformConstitutiveComponent(rConstitutiveMatrix( i, j ), rOriginalConstitutiveMatrix, rF,
                                                                         IndexVoigt2D4C[i][0], IndexVoigt2D4C[i][1], IndexVoigt2D4C[j][0], IndexVoigt2D4C[j][1]);
          }

        }
      }
      else if( size == 3 )
      {

        const unsigned int IndexVoigt2D3C [3][2] = { {0, 0}, {1, 1}, {0, 1} };

        for(unsigned int i=0; i<3; i++)
        {
          for(unsigned int j=0; j<3; j++)
          {
            rConstitutiveMatrix( i, j ) = TransformConstitutiveComponent(rConstitutiveMatrix( i, j ), rOriginalConstitutiveMatrix, rF,
                                                                         IndexVoigt2D3C[i][0], IndexVoigt2D3C[i][1], IndexVoigt2D3C[j][0], IndexVoigt2D3C[j][1]);
          }

        }
      }


    }


    /**
     * This method performs a pull-back or a push-forward between two constitutive tensor components
     */
    static inline double& TransformConstitutiveComponent(double & rCabcd,
                                                         const Matrix & rConstitutiveMatrix,
                                                         const Matrix & rF,
                                                         const unsigned int& a, const unsigned int& b,
                                                         const unsigned int& c, const unsigned int& d)

    {

      rCabcd = 0;
      double Cijkl=0;

      unsigned int dimension = rF.size1();

      //Cabcd
      for(unsigned int j=0; j<dimension; j++)
      {
        for(unsigned int l=0; l<dimension; l++)
        {
          for(unsigned int k=0; k<dimension; k++)
          {
            for(unsigned int i=0; i<dimension; i++)
            {
              //Cijkl
              rCabcd +=rF(a,i)*rF(b,j)*rF(c,k)*rF(d,l)*GetConstitutiveComponent(Cijkl,rConstitutiveMatrix,i,j,k,l);
            }
          }
        }
      }

      return rCabcd;

    }

    /**
     * This method gets the constitutive tensor components
     * from a consitutive matrix supplied in voigt notation
     */
    static inline double& GetConstitutiveComponent(double & rCabcd,
                                                   const Matrix& rConstitutiveMatrix,
                                                   const unsigned int& a, const unsigned int& b,
                                                   const unsigned int& c, const unsigned int& d)
    {
      // matrix indices
      unsigned int k=0, l= 0;

      unsigned int size = rConstitutiveMatrix.size1();

      if( size == 3 )
      {

        const unsigned int IndexVoigt2D3C [3][2] = { {0, 0}, {1, 1}, {0, 1} };

        //index k
        for(unsigned int i=0; i<3; i++)
        {
          if( a == b )
          {
            if( IndexVoigt2D3C[i][0] == a && IndexVoigt2D3C[i][1] == b )
            {
              k = i;
              break;
            }
          }
          else
          {
            if( (IndexVoigt2D3C[i][0] == a && IndexVoigt2D3C[i][1] == b) ||
                (IndexVoigt2D3C[i][1] == a && IndexVoigt2D3C[i][0] == b) )
            {
              k = i;
              break;
            }
          }
        }

        //index l
        for(unsigned int i=0; i<3; i++)
        {
          if( c == d )
          {
            if( IndexVoigt2D3C[i][0] == c && IndexVoigt2D3C[i][1] == d )
            {
              l = i;
              break;
            }
          }
          else
          {
            if( (IndexVoigt2D3C[i][0] == c && IndexVoigt2D3C[i][1] == d) ||
                (IndexVoigt2D3C[i][1] == c && IndexVoigt2D3C[i][0] == d) )
            {
              l = i;
              break;
            }
          }
        }


      }
      else if( size == 4 )
      {

        const unsigned int IndexVoigt2D4C [4][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1} };
        //index k
        for(unsigned int i=0; i<4; i++)
        {
          if( a == b )
          {
            if( IndexVoigt2D4C[i][0] == a && IndexVoigt2D4C[i][1] == b )
            {
              k = i;
              break;
            }
          }
          else
          {
            if( (IndexVoigt2D4C[i][0] == a && IndexVoigt2D4C[i][1] == b) ||
                (IndexVoigt2D4C[i][1] == a && IndexVoigt2D4C[i][0] == b) )
            {
              k = i;
              break;
            }
          }
        }

        //index l
        for(unsigned int i=0; i<4; i++)
        {
          if( c == d )
          {
            if( IndexVoigt2D4C[i][0] == c && IndexVoigt2D4C[i][1] == d )
            {
              l = i;
              break;
            }
          }
          else
          {
            if( (IndexVoigt2D4C[i][0] == c && IndexVoigt2D4C[i][1] == d) ||
                (IndexVoigt2D4C[i][1] == c && IndexVoigt2D4C[i][0] == d) )
            {
              l = i;
              break;
            }
          }
        }

      }
      else if( size == 6 )
      {

        const unsigned int IndexVoigt3D6C [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };

        //index k
        for(unsigned int i=0; i<6; i++)
        {
          if( a == b )
          {
            if( IndexVoigt3D6C[i][0] == a && IndexVoigt3D6C[i][1] == b )
            {
              k = i;
              break;
            }
          }
          else
          {
            if( (IndexVoigt3D6C[i][0] == a && IndexVoigt3D6C[i][1] == b) ||
                (IndexVoigt3D6C[i][1] == a && IndexVoigt3D6C[i][0] == b) )
            {
              k = i;
              break;
            }
          }
        }

        //index l
        for(unsigned int i=0; i<6; i++)
        {
          if( c == d )
          {
            if( IndexVoigt3D6C[i][0] == c && IndexVoigt3D6C[i][1] == d )
            {
              l = i;
              break;
            }
          }
          else
          {
            if( (IndexVoigt3D6C[i][0] == c && IndexVoigt3D6C[i][1] == d) ||
                (IndexVoigt3D6C[i][1] == c && IndexVoigt3D6C[i][0] == d) )
            {
              l = i;
              break;
            }
          }
        }
      }

      rCabcd = rConstitutiveMatrix(k,l);

      return rCabcd;
    }



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{


    ///@}


  protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

  private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}


  }; // Class ConstitutiveModelUtilities

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{
  ///@}

  ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_CONSTITUTIVE_MODEL_UTILITIES defined

