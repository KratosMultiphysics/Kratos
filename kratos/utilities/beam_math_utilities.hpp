//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Main author:      Josep Maria Carbonell
//  coming from       StringDynamicsApplication
//
//

#if !defined(KRATOS_BEAM_MATH_UTILITIES)
#define KRATOS_BEAM_MATH_UTILITIES

// System includes
#include <cmath>

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "utilities/quaternion.h"
#include "geometries/point.h"


namespace Kratos
{

template<class TDataType> class BeamMathUtils
{
public:
  ///name Type Definitions
  ///@{
  
  typedef Matrix                          MatrixType;
  typedef Vector                          VectorType;

  typedef Quaternion<double>          QuaternionType;   
  typedef MathUtils<TDataType>         MathUtilsType;
  typedef BeamMathUtils<TDataType> BeamMathUtilsType; 

  //typedef bounded_vector<double, 3>      PointType;
  typedef array_1d<double, 3>              PointType;

  ///@}
  ///name Math Utilities for beams
  ///@{

  //************************************************************************************
  //************************************************************************************

  /**
   * Transform a vector from the reference to the current local frame (MATERIAL frame for a beam)
   * @param rQuaternion: Quaternion representing the rotation from the reference to the current local frames
   * @param rVector: Vector to be rotated
   * the rotated vector rVector is returned. 
   */
  static inline VectorType& MapToCurrentLocalFrame(QuaternionType& rQuaternion, VectorType& rVector)
  {
    KRATOS_TRY

    //(rQuaternion.conjugate()).RotateVector3(rVariable); 
    // precision problems due to a rest included in the rotation

    //vector value :  v' = QT * v

    Matrix RotationMatrix;
    (rQuaternion.conjugate()).ToRotationMatrix(RotationMatrix);
    
    rVector = prod(RotationMatrix,rVector);

    return rVector;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  /**
   * Transform a vector from the reference to the current local frame (MATERIAL frame for a beam)
   * @param rQuaternion: Quaternion representing the rotation from the reference to the current local frames
   * @param rVector: Vector to be rotated
   * the rotated vector rVector is returned. 
   */
  static inline PointType& MapToCurrentLocalFrame(QuaternionType& rQuaternion, PointType& rVector)
  {
    KRATOS_TRY

    //(rQuaternion.conjugate()).RotateVector3(rVariable); 
    // precision problems due to a rest included in the rotation

    //vector value :  v' = QT * v

    Matrix RotationMatrix;
    (rQuaternion.conjugate()).ToRotationMatrix(RotationMatrix);
    
    rVector = prod(RotationMatrix,rVector);

    return rVector;
      
    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************


  /**
   * Transform a vector from the current to the reference local frame (SPATIAL frame for a beam)
   * @param rQuaternion: Quaternion representing the rotation from the reference to the current local frames
   * @param rVector: Vector to be rotated
   * the rotated vector rVariable is returned.
   */
  static inline VectorType& MapToReferenceLocalFrame(QuaternionType& rQuaternion, VectorType& rVector)
  {
    KRATOS_TRY

    //rQuaternion.RotateVector3(rVariable); 
    // precision problems due to a rest included in the rotation
      
    //vector value :  v = Q * v'

    Matrix RotationMatrix;
    rQuaternion.ToRotationMatrix(RotationMatrix);
    
    rVector = prod(RotationMatrix,rVector);

    return rVector;

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************


  /**
   * Transform a vector from the current to the reference local frame (SPATIAL frame for a beam)
   * @param rQuaternion: Quaternion representing the rotation from the reference to the current local frames
   * @param rVector: Vector to be rotated
   * the rotated vector rVariable is returned.
   */
  static inline PointType& MapToReferenceLocalFrame(QuaternionType& rQuaternion, PointType& rVector)
  {
    KRATOS_TRY
      
    //rQuaternion.RotateVector3(rVariable); 
    // precision problems due to a rest included in the rotation
      
    //vector value :  v = Q * v'

    Matrix RotationMatrix;
    rQuaternion.ToRotationMatrix(RotationMatrix);
    
    rVector = prod(RotationMatrix,rVector);

    return rVector;
    
    KRATOS_CATCH( "" )

  }

  ///@}


  //************************************************************************************
  //************************************************************************************

  /**
   * Transform a rotation vector to a rotation matrix with the Exponential transform
   * @param rVector: rotation vector (input parameter)
   * @param rExponentialTensor: rotation matrix (output parameter)
   */
  static inline void ExponentialTransform(const VectorType& rVector, MatrixType& rExponentialTensor)
  {  
    KRATOS_TRY

    //Initialize Local Matrices
    if( rExponentialTensor.size1() != 3 )
      rExponentialTensor.resize(3, 3, false);

    QuaternionType QuaternionValue = QuaternionType::FromRotationVector(rVector);
    QuaternionValue.ToRotationMatrix(rExponentialTensor);


    KRATOS_CATCH( "" )
  } 

  //************************************************************************************
  //************************************************************************************

  /**
   * Transform a rotation matrix to a rotation vector with the inverse of the Exponential transform
   * @param rCayleyTensor: rotation matrix (input parameter)
   * @param rVector: rotation vector (output parameter)
   */
  static inline void InverseExponentialTransform(const MatrixType& rExponentialTensor, VectorType& rVector)
  {  
    KRATOS_TRY
    
    //Initialize Local Matrices
    if( rVector.size() != 3 )
      rVector.resize(3,false);
    
    QuaternionType QuaternionValue = QuaternionType::FromRotationMatrix(rExponentialTensor);
    QuaternionValue.ToRotationVector(rVector);

    KRATOS_CATCH( "" )
  } 

  //************************************************************************************
  //************************************************************************************

  /**
   * Transform a rotation vector to a rotation matrix with the Cayley transform
   * @param rVector: rotation vector (input parameter)
   * @param rCayleyTensor: rotation matrix (output parameter)
   */
  static inline void CayleyTransform(const VectorType& rVector, MatrixType& rCayleyTensor)
  {  
    KRATOS_TRY

    //Initialize Local Matrices
    if( rCayleyTensor.size1() != 3 )
      rCayleyTensor.resize(3, 3, false);

    rCayleyTensor(0,0) = 0.5 * (4.0+rVector[0]*rVector[0]-rVector[1]*rVector[1]-rVector[2]*rVector[2]);
    rCayleyTensor(1,1) = 0.5 * (4.0-rVector[0]*rVector[0]+rVector[1]*rVector[1]-rVector[2]*rVector[2]);
    rCayleyTensor(2,2) = 0.5 * (4.0-rVector[0]*rVector[0]-rVector[1]*rVector[1]+rVector[2]*rVector[2]);
  
    rCayleyTensor(0,1) = (rVector[0]*rVector[1]-2.0*rVector[2]);
    rCayleyTensor(0,2) = (rVector[0]*rVector[2]+2.0*rVector[1]);
    rCayleyTensor(1,2) = (rVector[1]*rVector[2]-2.0*rVector[0]);
  
    rCayleyTensor(1,0) = (rVector[0]*rVector[1]+2.0*rVector[2]);
    rCayleyTensor(2,0) = (rVector[0]*rVector[2]-2.0*rVector[1]);
    rCayleyTensor(2,1) = (rVector[1]*rVector[2]+2.0*rVector[0]);
  
    rCayleyTensor *= (2.0 / (4.0+rVector[0]*rVector[0]+rVector[1]*rVector[1]+rVector[2]*rVector[2]));
  
    KRATOS_CATCH( "" )
  } 

  //***************************************************************************
  //***************************************************************************

  /**
   * Transform a rotation tensor to a rotation vector with the inverse of the Cayley transform
   * @param rCayleyTensor: rotation matrix (input parameter)
   * @param rVector: rotation vector (output parameter)
   */
  static inline void InverseCayleyTransform(const MatrixType& rCayleyTensor, VectorType& rVector)
  {  
    KRATOS_TRY
    
    //Initialize Local Matrices
    if( rVector.size() != 3 )
      rVector.resize(3,false);
    
    Matrix Identity = IdentityMatrix(3);     
    Matrix CayI     = rCayleyTensor + Identity;
    Matrix CayII    = rCayleyTensor - Identity; 

    Matrix InvertCayI;
    double det;
    MathUtilsType::InvertMatrix3(CayI, InvertCayI, det);
 
    Matrix SkewSymTensor = 2 * prod( InvertCayI, CayII );
      
    BeamMathUtilsType::SkewSymmetricTensorToVector( SkewSymTensor, rVector);
     
    KRATOS_CATCH( "" )
  } 

  //************************************************************************************
  //************************************************************************************

  /**
   * Transform a rotation vector to a skew symmetric matrix
   * @param rVector: rotation vector (input parameter)
   * @param rSkewSymmetricTensor: skew symmetric matrix (output parameter)
   */
  static inline void VectorToSkewSymmetricTensor(const VectorType& rVector, MatrixType& rSkewSymmetricTensor)
  {
    KRATOS_TRY

    //Initialize Local Matrices
    if( rSkewSymmetricTensor.size1() != 3 )
      rSkewSymmetricTensor.resize(3, 3, false);
    
    rSkewSymmetricTensor = ZeroMatrix(3,3);

    rSkewSymmetricTensor( 0, 1 ) = -rVector[2];
    rSkewSymmetricTensor( 0, 2 ) =  rVector[1];
    rSkewSymmetricTensor( 1, 2 ) = -rVector[0];

    rSkewSymmetricTensor( 1, 0 ) =  rVector[2];
    rSkewSymmetricTensor( 2, 0 ) = -rVector[1];
    rSkewSymmetricTensor( 2, 1 ) =  rVector[0];

    
    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************
  /**
   * Transform a skew symmetric matrix to a rotation vector
   * @param rSkewSymmetricTensor: skew symmetric matrix (input parameter)
   * @param rVector: rotation vector (output parameter)
   */
  static inline void SkewSymmetricTensorToVector(const MatrixType& rSkewSymmetricTensor, VectorType& rVector)
  { 
    KRATOS_TRY

    //Initialize Local Matrices
    if( rVector.size() != 3 )
      rVector.resize(3,false);
         
    rVector[0] = 0.5 * ( rSkewSymmetricTensor( 2, 1 ) - rSkewSymmetricTensor( 1, 2 ) );
    rVector[1] = 0.5 * ( rSkewSymmetricTensor( 0, 2 ) - rSkewSymmetricTensor( 2, 0 ) );
    rVector[2] = 0.5 * ( rSkewSymmetricTensor( 1, 0 ) - rSkewSymmetricTensor( 0, 1 ) );
    
         
    KRATOS_CATCH( "" )
   }

  //************************************************************************************
  //************************************************************************************
  /**
   * Get kroneckrDelta :: two order tensor that gives the value of the identity
   * @param i : coefficient i
   * @param j : coefficient j
   * Returns a double with the value
   */
  static inline double KroneckerDelta(int i, int j)
  { 
    KRATOS_TRY

    if( i==j )
      return 1;
    else
      return 0;
         
    KRATOS_CATCH( "" )
   }


  //************************************************************************************
  //************************************************************************************
  /**
   * Get Levi Civita Epsilon :: third order tensor that gives the value of a vectorial product
   * @param i : coefficient i
   * @param j : coefficient j
   * @param k : coefficient k
   * Returns a double with the sign
   */
  static inline double LeviCivitaEpsilon(int i, int j, int k)
  { 
    KRATOS_TRY

    if( i==j || j==k || k==i )
	return 0;
  
    if( (i==1 && j==2 && k==3) || (i==2 && j==3 && k==1) || (i==3 && j==1 && k==2) )
      return 1;

    if( (i==3 && j==2 && k==1) || (i==1 && j==3 && k==2) || (i==2 && j==1 && k==3) )
      return -1;

    return 0;
         
    KRATOS_CATCH( "" )
   }

  //************************************************************************************
  //************************************************************************************

  /**
   * Add a Nodal vector to a Elemental vector
   * @param rInputVector: LocalNodalVector (input parameter)
   * @param rOutputVector: LocalElementalVector (output parameter)
   * @param InitialRow: InitialRowNumber, initial index of the OutputVector
   * note the initialization of the outputvector must be done previously to the call of the method
   */
  static inline void AddVector(const VectorType& rInputVector, VectorType& rOutputVector, unsigned int InitialRow)
  {
    KRATOS_TRY
    
    for(unsigned int i=0; i<rInputVector.size(); i++)
      rOutputVector[InitialRow+i] += rInputVector[i];

    KRATOS_CATCH("")
   }

  //************************************************************************************
  //************************************************************************************

  /**
   * Substract a Nodal vector from an Elemental vector
   * @param rInputVector: LocalNodalVector (input parameter)
   * @param rOutputVector: LocalElementalVector (output parameter)
   * @param InitialRow: InitialRowNumber, initial index of the OutputVector
   * note the initialization of the outputvector must be done previously to the call of the method
   */
  static inline void SubstractVector(const VectorType& rInputVector, VectorType& rOutputVector, unsigned int InitialRow)
  {
    KRATOS_TRY

    for(unsigned int i=0; i<rInputVector.size(); i++)
      rOutputVector[InitialRow+i] -= rInputVector[i];

    KRATOS_CATCH("")
  }

  //*****************************************************************************
  //*****************************************************************************

  /**
   * Map a Matrix expressed on the Local frame of an element to a Global frame expression 
   * @param rLocalToGlobalMatrix: transformation matrix from local to global frame
   * @param rMatrix: matrix to be transformed (output parameter)
   * note the initialization of the Matrices must be done previously to the call of the method
   * return value :  A = Q * A' * QT
   */
  static inline void MapLocalToGlobal3D(const MatrixType& rLocalToGlobalMatrix, MatrixType& rMatrix)
  {

    KRATOS_TRY

    unsigned int MatSize = rMatrix.size1();

    Matrix AuxiliarRotationMatrix = ZeroMatrix(MatSize,MatSize);
 
    //Building the rotation matrix for the local element matrix
    for (unsigned int kk=0; kk < MatSize; kk += 3)
    {
        for (unsigned int i=0; i<3; i++)
        {
            for(unsigned int j=0; j<3; j++)
            {
	      AuxiliarRotationMatrix(i+kk,j+kk) = rLocalToGlobalMatrix(i,j);
            }
        }
    }

    //Rotate Local Stiffness Matrix
    Matrix aux_matrix   = ZeroMatrix(MatSize,MatSize);
    noalias(aux_matrix) = prod(AuxiliarRotationMatrix, rMatrix);

    //Stiffness Matrix
    rMatrix = ZeroMatrix(MatSize,MatSize);
    noalias(rMatrix) = prod(aux_matrix,trans(AuxiliarRotationMatrix));
         

    KRATOS_CATCH( "" )

  }

  //*****************************************************************************
  //*****************************************************************************

  /**
   * Map a Vector expressed on the Local frame of an element to a Global frame expression 
   * @param rLocalToGlobalMatrix: transformation matrix from local to global frame
   * @param rVector: vector to be transformed (output parameter)
   * note the initialization of the Matrices must be done previously to the call of the method
   */
  static inline void MapLocalToGlobal3D(const MatrixType& rLocalToGlobalMatrix, VectorType& rVector)
  {

    KRATOS_TRY

    unsigned int MatSize = rVector.size();

    Matrix AuxiliarRotationMatrix = ZeroMatrix(MatSize,MatSize);
 
    //Building the rotation matrix for the local element matrix
    for (unsigned int kk=0; kk < MatSize; kk += 3)
    {
        for (unsigned int i=0; i<3; i++)
        {
            for(unsigned int j=0; j<3; j++)
            {
	      AuxiliarRotationMatrix(i+kk,j+kk) = rLocalToGlobalMatrix(i,j);
            }
        }
    }

    rVector = prod(AuxiliarRotationMatrix, rVector);
          
    KRATOS_CATCH( "" )

  }


  //****************GID DEFINITION OF THE AUTOMATIC LOCAL AXES******************
  //*****************************************************************************

  /**
   * Deffault expression for GID local beam axes:: E3 is considered the local beam axial direction
   * @param rLocalZ: Local Beam axis vector (input parameter)
   * @param rRotationMatrix: transformation matrix from local to global frame (output parameter)
   */
  static inline  void CalculateLocalAxesMatrix(const VectorType& rLocalZ, MatrixType& rRotationMatrix)
  {

    KRATOS_TRY

    VectorType LocalX = ZeroVector(3);
    VectorType LocalY = ZeroVector(3);
    VectorType LocalZ = rLocalZ;

    BeamMathUtilsType::CalculateLocalAxesVectors(LocalZ,LocalX,LocalY);
        
    //Transformation matrix T = [e1_local, e2_local, e3_local] 
    if( rRotationMatrix.size1() != 3 )
      rRotationMatrix.resize(3, 3, false);
    
    //Building the rotation matrix
    for (unsigned int i=0; i<3; i++)
      {
    	rRotationMatrix(i,0) = LocalX[i];  // column distribution
    	rRotationMatrix(i,1) = LocalY[i];
    	rRotationMatrix(i,2) = LocalZ[i];	
      }
    
    KRATOS_CATCH( "" )

  }

  //*****************************************************************************
  //*****************************************************************************

  /**
   * Deffault expression for GID local beam axes:: E3 is considered the local beam axial direction
   * @param rLocalZ: Local Beam axis director vector E3 (input parameter) (output parameter)
   * @param rLocalX: Local Beam axis director vector E1 (output parameter)
   * @param rLocalY: Local Beam axis director vector E2 (output parameter)
   */
  static inline  void CalculateLocalAxesVectors(VectorType& rLocalZ, VectorType& rLocalX, VectorType& rLocalY)
  {

    KRATOS_TRY

    VectorType GlobalY = ZeroVector(3);
    GlobalY[1]=1.0;

    VectorType GlobalZ = ZeroVector(3);
    GlobalZ[2]=1.0;

    // local z-axis (e3_local) is the beam axis
    double VectorNorm = MathUtilsType::Norm(rLocalZ);
    if( VectorNorm != 0)
      rLocalZ /= VectorNorm;
    
    // local x-axis (e1_local)  
    double tolerance = 1.0/64.0;
    if(fabs(rLocalZ[0])< tolerance && fabs(rLocalZ[1])< tolerance){
      rLocalX = MathUtilsType::CrossProduct(GlobalY, rLocalZ);
    }
    else{
      rLocalX = MathUtilsType::CrossProduct(GlobalZ, rLocalZ);
    }
    
    VectorNorm = MathUtilsType::Norm(rLocalX);
    if( VectorNorm != 0)
      rLocalX /= VectorNorm;
    
    // local y-axis (e2_local)
    rLocalY = MathUtilsType::CrossProduct(rLocalZ,rLocalX);
    
    VectorNorm = MathUtilsType::Norm(rLocalY);
    if( VectorNorm != 0 )
      rLocalY /= VectorNorm;
        
    
    KRATOS_CATCH( "" )
  }

  //*****************************************************************************
  //*****************************************************************************

  /**
   * Check if a Tensor is orthogonal
   * @param rTensor: transformation matrix to be checked
   * @param Name: Name of the checked tensor to be verbosed
   * @param Verbose: Set verbosity for the check
   */
  static inline bool CheckOrthogonality( const Matrix& rTensor, std::string Name, bool Verbose )
  {  
    KRATOS_TRY
    	
    double Determinant    = MathUtilsType::Det(rTensor);
    Matrix IdentitySearch = prod( rTensor, trans(rTensor) );

    if( Determinant > 1.0001 || Determinant < 0.9999 ){

      if( Verbose ){
	std::cout<<Name<<std::endl;
	std::cout<<" Warning Matrix is not orthogonal "<<Determinant<<std::endl;
	std::cout<<" Matrix: "<<rTensor<<std::endl;
	std::cout<<" Identity not matches "<<IdentitySearch<<std::endl;
      }

      return false;
    }
  
    return true;

    KRATOS_CATCH( "" )
  }



private:
};// class BeamMathUtils
}
#endif /* KRATOS_BEAM_MATH_UTILITIESS defined */
