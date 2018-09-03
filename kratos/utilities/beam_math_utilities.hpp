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

  ///@}
  ///name Math Utilities for beams
  ///@{

  //************************************************************************************
  //************************************************************************************

  /**
   * Transform a vector from the reference to the current local frame (MATERIAL frame for a beam)
   * @param rQuaternion Quaternion representing the rotation from the reference to the current local frames
   * @param rVector Vector to be rotated
   * the rotated vector rVector is returned. 
   */
  template<class TVector3>
  static inline TVector3& MapToCurrentLocalFrame(QuaternionType& rQuaternion, TVector3& rVector)  
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
   * @param rQuaternion Quaternion representing the rotation from the reference to the current local frames
   * @param rVector Vector to be rotated
   * the rotated vector rVariable is returned.
   */
  template<class TVector3>
  static inline TVector3& MapToReferenceLocalFrame(QuaternionType& rQuaternion, TVector3& rVector)
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
   * Transform a rotation vector to a rotation matrix with the Exponential transform
   * @param rVector rotation vector (input parameter)
   * @param rExponentialTensor rotation matrix (output parameter)
   */
  template<class TVector3>
  static inline void ExponentialTransform(const TVector3& rVector, MatrixType& rExponentialTensor)
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
   * @param rCayleyTensor rotation matrix (input parameter)
   * @param rVector rotation vector (output parameter)
   */
  template<class TVector3>
  static inline void InverseExponentialTransform(const MatrixType& rExponentialTensor, TVector3& rVector)
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
   * @param rVector rotation vector (input parameter)
   * @param rCayleyTensor rotation matrix (output parameter)
   */
  template<class TVector3>
  static inline void CayleyTransform(const TVector3& rVector, MatrixType& rCayleyTensor)
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
   * @param rCayleyTensor rotation matrix (input parameter)
   * @param rVector rotation vector (output parameter)
   */
  template<class TVector3>
  static inline void InverseCayleyTransform(const MatrixType& rCayleyTensor, TVector3& rVector)
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
   * @param rVector rotation vector (input parameter)
   * @param rSkewSymmetricTensor skew symmetric matrix (output parameter)
   */
  template<class TVector3, class TMatrix3>
  static inline void VectorToSkewSymmetricTensor(const TVector3& rVector, TMatrix3& rSkewSymmetricTensor)
  {
    KRATOS_TRY

    //Initialize Local Matrices
    if( rSkewSymmetricTensor.size1() != 3 )
      rSkewSymmetricTensor.resize(3, 3, false);

    rSkewSymmetricTensor( 0, 0 ) = 0.0;
    rSkewSymmetricTensor( 1, 1 ) = 0.0;
    rSkewSymmetricTensor( 2, 2 ) = 0.0;
    
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
   * @param rSkewSymmetricTensor skew symmetric matrix (input parameter)
   * @param rVector rotation vector (output parameter)
   */
  template<class TMatrix3, class TVector3>
  static inline void SkewSymmetricTensorToVector(const TMatrix3& rSkewSymmetricTensor, TVector3& rVector)
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
   * @param i  coefficient i
   * @param j  coefficient j
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
   * @param i  coefficient i
   * @param j  coefficient j
   * @param k  coefficient k
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
   * @param rInputVector LocalNodalVector (input parameter)
   * @param rOutputVector LocalElementalVector (output parameter)
   * @param InitialRow InitialRowNumber, initial index of the OutputVector
   * note the initialization of the outputvector must be done previously to the call of the method
   */
  static inline void AddVector(const VectorType& rInputVector,VectorType& rOutputVector,const unsigned int InitialRow)
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
   * @param rInputVector LocalNodalVector (input parameter)
   * @param rOutputVector LocalElementalVector (output parameter)
   * @param InitialRow InitialRowNumber, initial index of the OutputVector
   * note the initialization of the outputvector must be done previously to the call of the method
   */
  static inline void SubstractVector(const VectorType& rInputVector,VectorType& rOutputVector,const unsigned int InitialRow)
  {
    KRATOS_TRY

    for(unsigned int i=0; i<rInputVector.size(); i++)
      rOutputVector[InitialRow+i] -= rInputVector[i];

    KRATOS_CATCH("")
  }


  /**
   * "InputMatrix" is ADDED to "Destination" matrix starting from
   * InitialRow and InitialCol of the destination matrix
   * "Destination" is assumed to be able to contain the "input matrix"
   * (no check is performed on the bounds)
   * @return rDestination: The matric destination
   * @param rInputMatrix The input matrix to be computed
   * @param rInitialRow The initial row to compute
   * @param rInitialCol The initial column to compute
   */
  
  template<class TMatrix, class TInputMatrix>
  static inline void  AddMatrix(TMatrix& rDestination,const TInputMatrix& rInputMatrix,const unsigned int rInitialRow,const unsigned int rInitialCol)
  {
    KRATOS_TRY
    for(unsigned int i = 0; i < rInputMatrix.size1(); i++)
      {
	for(unsigned int j = 0; j < rInputMatrix.size2(); j++)
	  {
	    rDestination(rInitialRow+i, rInitialCol+j) += rInputMatrix(i,j);
	  }
      }
    KRATOS_CATCH("")
  }
    
  /**
   *  "InputMatrix" is SUBTRACTED to "Destination" matrix starting from
   * InitialRow and InitialCol of the destination matrix
   * "Destination" is assumed to be able to contain the "input matrix"
   * (no check is performed on the bounds)
   * @return rDestination: The matric destination
   * @param rInputMatrix The input matrix to be computed
   * @param rInitialRow The initial row to compute
   * @param rInitialCol The initial column to compute
   */

  template<class TMatrix, class TInputMatrix>
  static inline void  SubtractMatrix(TMatrix& rDestination,const TInputMatrix& rInputMatrix,const unsigned int rInitialRow,const unsigned int rInitialCol)
  {
    KRATOS_TRY;
        
    for(unsigned int i = 0; i<rInputMatrix.size1(); i++)
      {
	for(unsigned int j = 0; j<rInputMatrix.size2(); j++)
	  {
	    rDestination(rInitialRow+i, rInitialCol+j) -= rInputMatrix(i,j);
	  }
      }
        
    KRATOS_CATCH("");
  }

  //*****************************************************************************
  //*****************************************************************************

  /**
   * Map a Matrix expressed on the Local frame of an element to a Global frame expression 
   * @param rLocalToGlobalQuaternion transformation quaternion from local to global frame
   * @param rMatrix matrix to be transformed (output parameter)
   * note the initialization of the Matrices must be done previously to the call of the method
   * return value :  A = Q * A' * QT
   */
  static inline void MapLocalToGlobal2D(const QuaternionType& rLocalToGlobalQuaternion, MatrixType& rMatrix)
  {

    KRATOS_TRY

    MatrixType LocalToGlobalMatrix;
    rLocalToGlobalQuaternion.ToRotationMatrix(LocalToGlobalMatrix);

    MapLocalToGlobal2D(LocalToGlobalMatrix,rMatrix);

    KRATOS_CATCH( "" )

  }
  
  //*****************************************************************************
  //*****************************************************************************

  /**
   * Map a Vector expressed on the Local frame of an element to a Global frame expression 
   * @param rLocalToGlobalQuaternion transformation quaternion from local to global frame
   * @param rVector vector to be transformed (output parameter)
   * note the initialization of the Matrices must be done previously to the call of the method
   */  
  static inline void MapLocalToGlobal2D(const QuaternionType& rLocalToGlobalQuaternion, VectorType& rVector)
  {

    KRATOS_TRY

    MatrixType LocalToGlobalMatrix;
    rLocalToGlobalQuaternion.ToRotationMatrix(LocalToGlobalMatrix);
  
    MapLocalToGlobal2D(LocalToGlobalMatrix,rVector);

    KRATOS_CATCH( "" )

  }

  
  //*****************************************************************************
  //*****************************************************************************

  /**
   * Map a Matrix expressed on the Local frame of an element to a Global frame expression 
   * @param rLocalToGlobalMatrix transformation matrix from local to global frame
   * @param rMatrix matrix to be transformed (output parameter)
   * note the initialization of the Matrices must be done previously to the call of the method
   * return value :  A = Q * A' * QT
   */
  static inline void MapLocalToGlobal2D(const MatrixType& rLocalToGlobalMatrix, MatrixType& rMatrix)
  {

    KRATOS_TRY

    unsigned int MatSize = rMatrix.size1();

    Matrix AuxiliarRotationMatrix(MatSize,MatSize);
    noalias(AuxiliarRotationMatrix) = ZeroMatrix(MatSize,MatSize);
 
    //Building the rotation matrix for the local element matrix
    for (unsigned int kk=0; kk < MatSize; kk += 2)
    {
        for (unsigned int i=0; i<2; i++)
        {
            for(unsigned int j=0; j<2; j++)
            {
	      AuxiliarRotationMatrix(i+kk,j+kk) = rLocalToGlobalMatrix(i,j);
            }
        }
    }

    //Rotate Local Stiffness Matrix
    Matrix aux_matrix(MatSize,MatSize);
    noalias(aux_matrix) = prod(AuxiliarRotationMatrix, rMatrix);

    //Stiffness Matrix
    noalias(rMatrix) = prod(aux_matrix,trans(AuxiliarRotationMatrix));
         

    KRATOS_CATCH( "" )

  }

  //*****************************************************************************
  //*****************************************************************************

  /**
   * Map a Vector expressed on the Local frame of an element to a Global frame expression 
   * @param rLocalToGlobalMatrix transformation matrix from local to global frame
   * @param rVector vector to be transformed (output parameter)
   * note the initialization of the Matrices must be done previously to the call of the method
   */
  static inline void MapLocalToGlobal2D(const MatrixType& rLocalToGlobalMatrix, VectorType& rVector)
  {

    KRATOS_TRY

    unsigned int MatSize = rVector.size();

    Matrix AuxiliarRotationMatrix(MatSize,MatSize);
    noalias(AuxiliarRotationMatrix) = ZeroMatrix(MatSize,MatSize);
 
    //Building the rotation matrix for the local element matrix
    for (unsigned int kk=0; kk < MatSize; kk += 2)
    {
        for (unsigned int i=0; i<2; i++)
        {
            for(unsigned int j=0; j<2; j++)
            {
	      AuxiliarRotationMatrix(i+kk,j+kk) = rLocalToGlobalMatrix(i,j);
            }
        }
    }

    rVector = prod(AuxiliarRotationMatrix, rVector);
          
    KRATOS_CATCH( "" )

  }

  //*****************************************************************************
  //*****************************************************************************

  /**
   * Map a Matrix expressed on the Local frame of an element to a Global frame expression 
   * @param rLocalToGlobalQuaternion transformation quaternion from local to global frame
   * @param rMatrix matrix to be transformed (output parameter)
   * note the initialization of the Matrices must be done previously to the call of the method
   * return value :  A = Q * A' * QT
   */
  static inline void MapLocalToGlobal3D(const QuaternionType& rLocalToGlobalQuaternion, MatrixType& rMatrix)
  {

    KRATOS_TRY

    MatrixType LocalToGlobalMatrix(3,3);
    rLocalToGlobalQuaternion.ToRotationMatrix(LocalToGlobalMatrix);
  
    MapLocalToGlobal3D(LocalToGlobalMatrix,rMatrix);

    KRATOS_CATCH( "" )

  }
  
  //*****************************************************************************
  //*****************************************************************************

  /**
   * Map a Vector expressed on the Local frame of an element to a Global frame expression 
   * @param rLocalToGlobalQuaternion transformation quaternion from local to global frame
   * @param rVector vector to be transformed (output parameter)
   * note the initialization of the Matrices must be done previously to the call of the method
   */
  static inline void MapLocalToGlobal3D(const QuaternionType& rLocalToGlobalQuaternion, VectorType& rVector)
  {

    KRATOS_TRY

    MatrixType LocalToGlobalMatrix(3,3);
    rLocalToGlobalQuaternion.ToRotationMatrix(LocalToGlobalMatrix);
  
    MapLocalToGlobal3D(LocalToGlobalMatrix,rVector);

    KRATOS_CATCH( "" )

  }
  
  //*****************************************************************************
  //*****************************************************************************

  /**
   * Map a Matrix expressed on the Local frame of an element to a Global frame expression 
   * @param rLocalToGlobalMatrix transformation matrix from local to global frame
   * @param rMatrix matrix to be transformed (output parameter)
   * note the initialization of the Matrices must be done previously to the call of the method
   * return value :  A = Q * A' * QT
   */
  static inline void MapLocalToGlobal3D(const MatrixType& rLocalToGlobalMatrix, MatrixType& rMatrix)
  {

    KRATOS_TRY

    unsigned int MatSize = rMatrix.size1();

    Matrix AuxiliarRotationMatrix(MatSize,MatSize);
    noalias(AuxiliarRotationMatrix) = ZeroMatrix(MatSize,MatSize);

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
    Matrix aux_matrix(MatSize,MatSize);
    noalias(aux_matrix) = prod(AuxiliarRotationMatrix, rMatrix);

    //Stiffness Matrix
    noalias(rMatrix) = prod(aux_matrix,trans(AuxiliarRotationMatrix));
         

    KRATOS_CATCH( "" )

  }

  //*****************************************************************************
  //*****************************************************************************

  /**
   * Map a Vector expressed on the Local frame of an element to a Global frame expression 
   * @param rLocalToGlobalMatrix transformation matrix from local to global frame
   * @param rVector vector to be transformed (output parameter)
   * note the initialization of the Matrices must be done previously to the call of the method
   */
  static inline void MapLocalToGlobal3D(const MatrixType& rLocalToGlobalMatrix, VectorType& rVector)
  {

    KRATOS_TRY

    unsigned int MatSize = rVector.size();

    Matrix AuxiliarRotationMatrix(MatSize,MatSize);
    noalias(AuxiliarRotationMatrix) = ZeroMatrix(MatSize,MatSize);

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
   * Deffault expression for GID local beam axes:: E1 is considered the local beam axial direction
   * @param rLocalX Local Beam axis vector (input parameter)
   * @param rRotationMatrix transformation matrix from local to global frame (output parameter)
   */
  template<class TVector3>
  static inline  void CalculateLocalAxesMatrix(const TVector3& rLocalX, MatrixType& rRotationMatrix)
  {

    KRATOS_TRY

    TVector3 LocalX = rLocalX;
    TVector3 LocalY = ZeroVector(3);
    TVector3 LocalZ = ZeroVector(3);

    BeamMathUtilsType::CalculateLocalAxesVectors(LocalX,LocalY,LocalZ);
        
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
   * Deffault expression for GID local beam axes:: E1 is considered the local beam axial direction
   * @param rLocalX Local Beam axis director vector E1 (input parameter) (output parameter)
   * @param rLocalY Local Beam axis director vector E2 (output parameter)
   * @param rLocalZ Local Beam axis director vector E3 (output parameter)
   */
  template<class TVector3>
  static inline  void CalculateLocalAxesVectors(TVector3& rLocalX, TVector3& rLocalY, TVector3& rLocalZ)
  {

    KRATOS_TRY

    TVector3 GlobalY = ZeroVector(3);
    GlobalY[1]=1.0;

    TVector3 GlobalZ = ZeroVector(3);
    GlobalZ[2]=1.0;

    // local x-axis (e1_local) is the beam axis
    double VectorNorm = MathUtilsType::Norm(rLocalX);
    if( VectorNorm != 0)
      rLocalX /= VectorNorm;
    
    // local y-axis (e2_local)  
    double tolerance = 1.0/64.0;
    if(fabs(rLocalX[0])< tolerance && fabs(rLocalX[1])< tolerance){
      MathUtilsType::CrossProduct(rLocalY,GlobalY,rLocalX);
    }
    else{
      MathUtilsType::CrossProduct(rLocalY,GlobalZ,rLocalX);
    }
    
    VectorNorm = MathUtilsType::Norm(rLocalY);
    if( VectorNorm != 0)
      rLocalY /= VectorNorm;
    
    // local z-axis (e3_local)
    MathUtilsType::CrossProduct(rLocalZ,rLocalX,rLocalY);
    
    VectorNorm = MathUtilsType::Norm(rLocalZ);
    if( VectorNorm != 0 )
      rLocalZ /= VectorNorm;
        
    
    KRATOS_CATCH( "" )
  }


  //****************GID DEFINITION OF THE CUSTOM LOCAL AXES***********************
  //*****************************************************************************

  /**
   * Deffault expression for GID local beam axes:: E1 is considered the local beam axial direction
   * @param rLocalX Local Beam axis vector (input parameter)
   * @param rLocalY Local axis 2 vector (input parameter)
   * @param rRotationMatrix transformation matrix from local to global frame (output parameter)
   */
  template<class TVector3>
  static inline  void CalculateLocalAxesMatrix(const TVector3& rLocalX, const TVector3& rLocalY, MatrixType& rRotationMatrix)
  {

    KRATOS_TRY

    TVector3 LocalX = rLocalX;
    TVector3 LocalY = rLocalY;
    TVector3 LocalZ = ZeroVector(3);

    BeamMathUtilsType::CalculateLocalAxisVector(LocalX,LocalY,LocalZ);
        
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
   * Deffault expression for GID local beam axes:: E1 is considered the local beam axial direction
   * @param rLocalX Local Beam axis director vector E1 (input parameter) (output parameter)
   * @param rLocalY Local Beam axis director vector E2 (input parameter) (output parameter)
   * @param rLocalZ Local Beam axis director vector E3 (output parameter)
   */
  template<class TVector3>
  static inline  void CalculateLocalAxisVector(TVector3& rLocalX, TVector3& rLocalY, TVector3& rLocalZ)
  {
    KRATOS_TRY

    // local x-axis (e1_local) is the beam axis
    double VectorNorm = MathUtilsType::Norm(rLocalX);
    if( VectorNorm != 0)
      rLocalX /= VectorNorm;
    
    // local y-axis (e2_local) supplied
    VectorNorm = MathUtilsType::Norm(rLocalY);
    if( VectorNorm != 0)
      rLocalY /= VectorNorm;
    
    // local z-axis (e3_local)
    MathUtilsType::CrossProduct(rLocalZ,rLocalX,rLocalY);
    
    VectorNorm = MathUtilsType::Norm(rLocalZ);
    if( VectorNorm != 0 )
      rLocalZ /= VectorNorm;
            
    KRATOS_CATCH( "" )
  }

  //*****************************************************************************
  //*****************************************************************************

  /**
   * Check if a Tensor is orthogonal
   * @param rTensor transformation matrix to be checked
   * @param Name Name of the checked tensor to be verbosed
   * @param Verbose Set verbosity for the check
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
