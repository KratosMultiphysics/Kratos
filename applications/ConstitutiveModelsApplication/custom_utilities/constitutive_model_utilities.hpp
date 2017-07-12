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
    typedef bounded_matrix<double,3,3>    MatrixType;
	
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
    
    static inline Matrix& ConstitutiveTensorToMatrix(const bounded_matrix<double,6,6>& rTensor, Matrix& rMatrix)
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
	  
	  rMatrix(2,2) = rTensor(2,2);

	}
	else{
	  KRATOS_ERROR << "Constitutive Matrix dimensions are not correct" << std::endl;
	}
        
	return rMatrix;
	
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
    
    static inline array_1d<double,6>& SymmetricTensorToVector(const MatrixType& rMatrix, array_1d<double,6>& rVector)
    {
        KRATOS_TRY;
        
	rVector[0]= rMatrix(0,0);
	rVector[1]= rMatrix(1,1);
	rVector[2]= rMatrix(2,2);
	
	rVector[3]= rMatrix(0,1);
	rVector[4]= rMatrix(1,2);
	rVector[5]= rMatrix(2,0);

        return rVector;
        
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
    
    static inline array_1d<double,6>& StrainTensorToVector(const MatrixType& rMatrix, array_1d<double,6>& rVector)
    {
        KRATOS_TRY;
        
	rVector[0]= rMatrix(0,0);
	rVector[1]= rMatrix(1,1);
	rVector[2]= rMatrix(2,2);
	rVector[3]= 2.0*rMatrix(0,1);
	rVector[4]= 2.0*rMatrix(1,2);
	rVector[5]= 2.0*rMatrix(0,2);

        return rVector;
        
        KRATOS_CATCH("");
     }
    
    /**
     * Transforms a given symmetric Stress Tensor to Voigt Notation:
     * in the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     * in the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (4*1) Vector
     * in the 2D case: from a second order tensor (3*3) Matrix  to a corresponing (3*1) Vector
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
	rCharacteristicSize = sqrt(4.0*rDomainGeometry.Area()/KRATOS_M_PI);

      if( rDomainGeometry.WorkingSpaceDimension() == 3)
	//rCharacteristicSize is the diameter of a sphere with the same volume as the element
	rCharacteristicSize = pow((6.0*rDomainGeometry.Volume()/KRATOS_M_PI),0.33333333333333);

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
        if (r <= -1) { phi = KRATOS_M_PI / 3.0; }
        else if (r >= 1) { phi = 0.0; }
        else { phi = acos(r) / 3.0;}
            
        // the eigenvalues satisfy eig3 <= eig2 <= eig1
        Result[0] = q + 2.0 * p * cos(phi);
        Result[2] = q + 2.0 * p * cos(phi + (2.0*KRATOS_M_PI/3.0));
        Result[1] = 3.0 * q - Result[0] - Result[2];     //% since trace(A) = eig1 + eig2 + eig3   

        return Result;
    }


    /**
     * It inverts matrices of order 3 //VERIFIED!!!
     * @param InputMatrix: Is the input matrix (unchanged at output)
     * @return InvertedMatrix: Is the inverse of the input matrix
     * @return InputMatrixDet: Is the determinant of the input matrix
     */
    
    static void InvertMatrix3(
        const MatrixType& InputMatrix,
        MatrixType& InvertedMatrix,
        double& InputMatrixDet
        ) 
    {
        KRATOS_TRY;
        
        if(InvertedMatrix.size1() != 3 || InvertedMatrix.size2() != 3)
        {
            InvertedMatrix.resize(3,3,false);
        }

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
        InvertedMatrix /= InputMatrixDet;
        
        KRATOS_CATCH("")
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

