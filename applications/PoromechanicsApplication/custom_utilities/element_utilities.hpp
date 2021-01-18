//
//   Project Name:        KratosPoromechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ELEMENT_UTILITIES_H_INCLUDED )
#define  KRATOS_ELEMENT_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class ElementUtilities
{
 public:

  ///@name Type Definitions
  ///@{

  //definition of the size type
  typedef std::size_t SizeType;

  ///definition of node type (default is: Node<3>)
  typedef Node<3> NodeType;

  ///definition of the geometry type with given NodeType
  typedef Geometry<NodeType> GeometryType;
  ///@}


  /**
   * @brief Calculate Delta Position
   * @param rDeltaPosition, matrix storing the displacement or position increment, returned parameter
   * @param rGeometry, geometry where the gradient is calculated
   */
  static inline void CalculateDeltaPosition(Matrix& rDeltaPosition, const GeometryType& rGeometry)
  {

    const SizeType number_of_nodes = rGeometry.PointsNumber();
    const SizeType dimension  = rGeometry.WorkingSpaceDimension();

    if(rDeltaPosition.size1() != number_of_nodes || rDeltaPosition.size2() !=  dimension)
      rDeltaPosition.resize(number_of_nodes,dimension,false);

    //noalias(rDeltaPosition) = ZeroMatrix(number_of_nodes,dimension);

    if( rGeometry[0].SolutionStepsDataHas(STEP_DISPLACEMENT) )
    {
      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        const array_1d<double, 3 > & CurrentStepDisplacement = rGeometry[i].FastGetSolutionStepValue(STEP_DISPLACEMENT,0);

        for ( SizeType j = 0; j < dimension; j++ )
        {
          rDeltaPosition(i,j) = CurrentStepDisplacement[j];
        }

      }
    }
    else{

      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        const array_1d<double, 3 > & CurrentDisplacement  = rGeometry[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 > & PreviousDisplacement = rGeometry[i].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( SizeType j = 0; j < dimension; j++ )
        {
          rDeltaPosition(i,j) = CurrentDisplacement[j]-PreviousDisplacement[j];
        }
      }

    }

  }

  /**
   * @brief Calculate Total Delta Position
   * @param rDeltaPosition, matrix storing the displacement or position increment from origin, returned parameter
   * @param rGeometry, geometry where the gradient is calculated
   */
  static inline void CalculateTotalDeltaPosition(Matrix & rDeltaPosition, const GeometryType& rGeometry)
  {
    const SizeType number_of_nodes = rGeometry.PointsNumber();
    const SizeType dimension  = rGeometry.WorkingSpaceDimension();

    if(rDeltaPosition.size1() != number_of_nodes || rDeltaPosition.size2() !=  dimension)
      rDeltaPosition.resize(number_of_nodes,dimension,false);

    //noalias(rDeltaPosition) = ZeroMatrix(number_of_nodes,dimension);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      const array_1d<double, 3 > & CurrentDisplacement  = rGeometry[i].FastGetSolutionStepValue(DISPLACEMENT);

      for ( SizeType j = 0; j < dimension; j++ )
      {
        rDeltaPosition(i,j) = CurrentDisplacement[j];
      }
    }

  }

  /**
   * @brief Calculate Norm of stresses.VelocityGradient
   * @param rVelocityGradient, matrix form of the velocity gradient, returned parameter
   * @param rGeometry, geometry where the gradient is calculated
   * @param rDN_DX, shape functions derivatives
   * @param Alpha, parameter to change the step calculation [0,1]
   */
  static inline void CalculateVelocityGradient(Matrix& rVelocityGradient, const GeometryType& rGeometry,
                                               const Matrix& rDN_DX, const double Alpha = 1.0)
  {

    const SizeType number_of_nodes  = rGeometry.PointsNumber();
    const SizeType dimension        = rGeometry.WorkingSpaceDimension();

    if( rVelocityGradient.size1() != dimension || rVelocityGradient.size2() != dimension )
      rVelocityGradient.resize(dimension,dimension);

    noalias(rVelocityGradient) = ZeroMatrix(dimension,dimension);

    if( Alpha != 1.0 ){

      if( dimension == 2 )
      {
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          const array_1d<double,3>& rPreviousVelocity = rGeometry[i].FastGetSolutionStepValue(VELOCITY,1);
          const array_1d<double,3>& rCurrentVelocity  = rGeometry[i].FastGetSolutionStepValue(VELOCITY);
          rVelocityGradient ( 0 , 0 ) += (rCurrentVelocity[0] * Alpha + rPreviousVelocity[0]* (1.0-Alpha))*rDN_DX ( i , 0 );
          rVelocityGradient ( 0 , 1 ) += (rCurrentVelocity[0] * Alpha + rPreviousVelocity[0]* (1.0-Alpha))*rDN_DX ( i , 1 );
          rVelocityGradient ( 1 , 0 ) += (rCurrentVelocity[1] * Alpha + rPreviousVelocity[1]* (1.0-Alpha))*rDN_DX ( i , 0 );
          rVelocityGradient ( 1 , 1 ) += (rCurrentVelocity[1] * Alpha + rPreviousVelocity[1]* (1.0-Alpha))*rDN_DX ( i , 1 );
        }
      }
      else if( dimension == 3)
      {
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          const array_1d<double,3>& rPreviousVelocity = rGeometry[i].FastGetSolutionStepValue(VELOCITY,1);
          const array_1d<double,3>& rCurrentVelocity  = rGeometry[i].FastGetSolutionStepValue(VELOCITY);
          rVelocityGradient ( 0 , 0 ) += (rCurrentVelocity[0] * Alpha + rPreviousVelocity[0]* (1.0-Alpha))*rDN_DX ( i , 0 );
          rVelocityGradient ( 0 , 1 ) += (rCurrentVelocity[0] * Alpha + rPreviousVelocity[0]* (1.0-Alpha))*rDN_DX ( i , 1 );
          rVelocityGradient ( 0 , 2 ) += (rCurrentVelocity[0] * Alpha + rPreviousVelocity[0]* (1.0-Alpha))*rDN_DX ( i , 2 );
          rVelocityGradient ( 1 , 0 ) += (rCurrentVelocity[1] * Alpha + rPreviousVelocity[1]* (1.0-Alpha))*rDN_DX ( i , 0 );
          rVelocityGradient ( 1 , 1 ) += (rCurrentVelocity[1] * Alpha + rPreviousVelocity[1]* (1.0-Alpha))*rDN_DX ( i , 1 );
          rVelocityGradient ( 1 , 2 ) += (rCurrentVelocity[1] * Alpha + rPreviousVelocity[1]* (1.0-Alpha))*rDN_DX ( i , 2 );
          rVelocityGradient ( 2 , 0 ) += (rCurrentVelocity[2] * Alpha + rPreviousVelocity[2]* (1.0-Alpha))*rDN_DX ( i , 0 );
          rVelocityGradient ( 2 , 1 ) += (rCurrentVelocity[2] * Alpha + rPreviousVelocity[2]* (1.0-Alpha))*rDN_DX ( i , 1 );
          rVelocityGradient ( 2 , 2 ) += (rCurrentVelocity[2] * Alpha + rPreviousVelocity[2]* (1.0-Alpha))*rDN_DX ( i , 2 );
        }
      }

    }
    else{

      if( dimension == 2 )
      {
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          const array_1d<double,3>& rCurrentVelocity = rGeometry[i].FastGetSolutionStepValue(VELOCITY);
          rVelocityGradient ( 0 , 0 ) += rCurrentVelocity[0]*rDN_DX ( i , 0 );
          rVelocityGradient ( 0 , 1 ) += rCurrentVelocity[0]*rDN_DX ( i , 1 );
          rVelocityGradient ( 1 , 0 ) += rCurrentVelocity[1]*rDN_DX ( i , 0 );
          rVelocityGradient ( 1 , 1 ) += rCurrentVelocity[1]*rDN_DX ( i , 1 );
        }

      }
      else if( dimension == 3)
      {
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          const array_1d<double,3>& rCurrentVelocity = rGeometry[i].FastGetSolutionStepValue(VELOCITY);
          rVelocityGradient ( 0 , 0 ) += rCurrentVelocity[0]*rDN_DX ( i , 0 );
          rVelocityGradient ( 0 , 1 ) += rCurrentVelocity[0]*rDN_DX ( i , 1 );
          rVelocityGradient ( 0 , 2 ) += rCurrentVelocity[0]*rDN_DX ( i , 2 );
          rVelocityGradient ( 1 , 0 ) += rCurrentVelocity[1]*rDN_DX ( i , 0 );
          rVelocityGradient ( 1 , 1 ) += rCurrentVelocity[1]*rDN_DX ( i , 1 );
          rVelocityGradient ( 1 , 2 ) += rCurrentVelocity[1]*rDN_DX ( i , 2 );
          rVelocityGradient ( 2 , 0 ) += rCurrentVelocity[2]*rDN_DX ( i , 0 );
          rVelocityGradient ( 2 , 1 ) += rCurrentVelocity[2]*rDN_DX ( i , 1 );
          rVelocityGradient ( 2 , 2 ) += rCurrentVelocity[2]*rDN_DX ( i , 2 );
        }
      }
    }

  }

  /**
   * @brief Calculate the Deformation Gradient Tensor
   * @param rVelocityGradient, matrix form for the deformation gradient, returned parameter
   * @param rGeometry, geometry where the gradient is calculated
   * @param rDN_DX, shape functions derivatives
   * @param rDeltaPosition, matrix containing increment of position
   */
  static inline void CalculateDeformationGradient(Matrix& rDeformationGradient, const GeometryType& rGeometry,
                                                  const Matrix& rDN_DX, const Matrix& rDeltaPosition)
  {
    const SizeType number_of_nodes = rGeometry.PointsNumber();
    const SizeType dimension       = rGeometry.WorkingSpaceDimension();

    if( rDeformationGradient.size1() != dimension || rDeformationGradient.size2() != dimension )
      rDeformationGradient.resize(dimension, dimension, false);

    noalias(rDeformationGradient) = IdentityMatrix(dimension);

    if( dimension == 2 )
    {
      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        rDeformationGradient ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
        rDeformationGradient ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
        rDeformationGradient ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
        rDeformationGradient ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
      }

    }
    else if( dimension == 3)
    {
      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        rDeformationGradient ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
        rDeformationGradient ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
        rDeformationGradient ( 0 , 2 ) += rDeltaPosition(i,0)*rDN_DX ( i , 2 );
        rDeformationGradient ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
        rDeformationGradient ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
        rDeformationGradient ( 1 , 2 ) += rDeltaPosition(i,1)*rDN_DX ( i , 2 );
        rDeformationGradient ( 2 , 0 ) += rDeltaPosition(i,2)*rDN_DX ( i , 0 );
        rDeformationGradient ( 2 , 1 ) += rDeltaPosition(i,2)*rDN_DX ( i , 1 );
        rDeformationGradient ( 2 , 2 ) += rDeltaPosition(i,2)*rDN_DX ( i , 2 );
      }
    }
    else
    {
      KRATOS_ERROR << " something is wrong with the dimension when computing the Deformation Gradient " << std::endl;
    }

  }

  /**
   * @brief Calculate the VelocityGradient vector (no voigt form)
   * @param rVelocityGradient, vector form of the non symmetric velocity gradient, returned parameter
   * @param rGeometry, geometry where the gradient is calculated
   * @param rDN_DX, shape functions derivatives
   * @param Alpha, parameter to change the step calculation [0,1]
   */
  static inline void CalculateVelocityGradientVector(Vector& rVelocityGradient, const GeometryType& rGeometry,
                                                     const Matrix& rDN_DX, const double Alpha = 1.0)
  {

    const SizeType number_of_nodes  = rGeometry.PointsNumber();
    const SizeType dimension        = rGeometry.WorkingSpaceDimension();

    if( rVelocityGradient.size() != dimension*dimension )
      rVelocityGradient.resize(dimension*dimension);

    noalias(rVelocityGradient) = ZeroVector(dimension * dimension);

    array_1d<double,3> Velocity;
    if( dimension == 2 )
    {
      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        const array_1d<double,3>& rPreviousVelocity = rGeometry[i].FastGetSolutionStepValue(VELOCITY,1);
        const array_1d<double,3>& rCurrentVelocity  = rGeometry[i].FastGetSolutionStepValue(VELOCITY);

        Velocity = rCurrentVelocity * Alpha + rPreviousVelocity * (1.0-Alpha);

        rVelocityGradient[0] += Velocity[0]*rDN_DX ( i , 0 );
        rVelocityGradient[1] += Velocity[1]*rDN_DX ( i , 1 );
        rVelocityGradient[2] += Velocity[0]*rDN_DX ( i , 1 );
        rVelocityGradient[3] += Velocity[1]*rDN_DX ( i , 0 );
      }
    }
    else if( dimension == 3)
    {

      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        const array_1d<double,3>& rPreviousVelocity = rGeometry[i].FastGetSolutionStepValue(VELOCITY,1);
        const array_1d<double,3>& rCurrentVelocity  = rGeometry[i].FastGetSolutionStepValue(VELOCITY);

        Velocity = rCurrentVelocity * Alpha + rPreviousVelocity * (1.0-Alpha);

        rVelocityGradient[0] += Velocity[0]*rDN_DX ( i , 0 );
        rVelocityGradient[1] += Velocity[1]*rDN_DX ( i , 1 );
        rVelocityGradient[2] += Velocity[2]*rDN_DX ( i , 2 );

        rVelocityGradient[3] += Velocity[0]*rDN_DX ( i , 1 );
        rVelocityGradient[4] += Velocity[1]*rDN_DX ( i , 2 );
        rVelocityGradient[5] += Velocity[2]*rDN_DX ( i , 0 );

        rVelocityGradient[6] += Velocity[1]*rDN_DX ( i , 0 );
        rVelocityGradient[7] += Velocity[2]*rDN_DX ( i , 1 );
        rVelocityGradient[8] += Velocity[0]*rDN_DX ( i , 2 );
      }

    }
    else
    {
      KRATOS_ERROR << " something is wrong with the dimension when computing symmetric velocity gradient " << std::endl;
    }


  }

  /**
   * @brief Calculate the symmetric VelocityGradient vector
   * @param rVelocityGradientMatrix, matrix form of the velocity gradient
   * @param rSymmetricVelocityGradientVector, vector form of the symmetric velocity gradient, returned parameter
   */
  static inline void CalculateSymmetricVelocityGradientVector(const Matrix& rVelocityGradientMatrix,
                                                              Vector& rSymmetricVelocityGradientVector,
                                                              const SizeType& rDimension)
  {

    if( rDimension == 2 )
    {
      if ( rSymmetricVelocityGradientVector.size() != 3 ) rSymmetricVelocityGradientVector.resize( 3, false );

      rSymmetricVelocityGradientVector[0] = rVelocityGradientMatrix( 0, 0 );
      rSymmetricVelocityGradientVector[1] = rVelocityGradientMatrix( 1, 1 );
      rSymmetricVelocityGradientVector[2] = 0.5 * (rVelocityGradientMatrix( 0, 1 ) + rVelocityGradientMatrix( 1, 0 )); // xy

    }
    else if( rDimension == 3 )
    {
      if ( rSymmetricVelocityGradientVector.size() != 6 ) rSymmetricVelocityGradientVector.resize( 6, false );

      rSymmetricVelocityGradientVector[0] = rVelocityGradientMatrix( 0, 0 );
      rSymmetricVelocityGradientVector[1] = rVelocityGradientMatrix( 1, 1 );
      rSymmetricVelocityGradientVector[2] = rVelocityGradientMatrix( 2, 2 );
      rSymmetricVelocityGradientVector[3] = 0.5 * ( rVelocityGradientMatrix( 0, 1 ) + rVelocityGradientMatrix( 1, 0 ) ); // xy
      rSymmetricVelocityGradientVector[4] = 0.5 * ( rVelocityGradientMatrix( 1, 2 ) + rVelocityGradientMatrix( 2, 1 ) ); // yz
      rSymmetricVelocityGradientVector[5] = 0.5 * ( rVelocityGradientMatrix( 0, 2 ) + rVelocityGradientMatrix( 2, 0 ) ); // xz

    }
    else
    {
      KRATOS_ERROR << " something is wrong with the dimension symmetric velocity gradient " << std::endl;
    }

  }

  /**
   * @brief Calculate the skew-symmetric VelocityGradient vector
   * @param rVelocityGradientMatrix, matrix form of the velocity gradient
   * @param rSkewSymmetricVelocityGradientVector, vector form of the symmetric velocity gradient, returned parameter
   */
  static inline void CalculateSkewSymmetricVelocityGradientVector(const Matrix& rVelocityGradientMatrix,
                                                                  Vector& rSkewSymmetricVelocityGradientVector,
                                                                  const SizeType& rDimension)
  {

    if( rDimension == 2 )
    {
      if ( rSkewSymmetricVelocityGradientVector.size() != 3 ) rSkewSymmetricVelocityGradientVector.resize( 3, false );

      rSkewSymmetricVelocityGradientVector[0] = 0.0;
      rSkewSymmetricVelocityGradientVector[1] = 0.0;
      rSkewSymmetricVelocityGradientVector[2] = 0.5 * (rVelocityGradientMatrix( 0, 1 ) - rVelocityGradientMatrix( 1, 0 )); // xy

    }
    else if( rDimension == 3 )
    {
      if ( rSkewSymmetricVelocityGradientVector.size() != 6 ) rSkewSymmetricVelocityGradientVector.resize( 6, false );

      rSkewSymmetricVelocityGradientVector[0] = 0.0;
      rSkewSymmetricVelocityGradientVector[1] = 0.0;
      rSkewSymmetricVelocityGradientVector[2] = 0.0;
      rSkewSymmetricVelocityGradientVector[3] = 0.5 * ( rVelocityGradientMatrix( 0, 1 ) - rVelocityGradientMatrix( 1, 0 ) ); // xy
      rSkewSymmetricVelocityGradientVector[4] = 0.5 * ( rVelocityGradientMatrix( 1, 2 ) - rVelocityGradientMatrix( 2, 1 ) ); // yz
      rSkewSymmetricVelocityGradientVector[5] = 0.5 * ( rVelocityGradientMatrix( 0, 2 ) - rVelocityGradientMatrix( 2, 0 ) ); // xz

    }
    else
    {
      KRATOS_ERROR << " something is wrong with the dimension skew symmetric velocity gradient " << std::endl;
    }

  }


  /**
   * @brief Calculate Linear deformation matrix BL
   * @param rDeformationMatrix, matrix form, returned parameter
   * @param rGeometry, geometry where the gradient is calculated
   * @param rDN_DX, shape functions derivatives
   */
  static inline void CalculateLinearDeformationMatrix(Matrix& rDeformationMatrix, const GeometryType& rGeometry, const Matrix& rDN_DX)
  {

    const SizeType number_of_nodes  = rGeometry.PointsNumber();
    const SizeType dimension        = rGeometry.WorkingSpaceDimension();
    unsigned int voigt_size         = dimension * (dimension +1) * 0.5;

    if ( rDeformationMatrix.size1() != voigt_size || rDeformationMatrix.size2() != dimension*number_of_nodes )
      rDeformationMatrix.resize(voigt_size, dimension*number_of_nodes, false );


    if( dimension == 2 )
    {

      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        unsigned int index = 2 * i;

        rDeformationMatrix( 0, index + 0 ) = rDN_DX( i, 0 );
        rDeformationMatrix( 0, index + 1 ) = 0.0;
        rDeformationMatrix( 1, index + 0 ) = 0.0;
        rDeformationMatrix( 1, index + 1 ) = rDN_DX( i, 1 );
        rDeformationMatrix( 2, index + 0 ) = rDN_DX( i, 1 );
        rDeformationMatrix( 2, index + 1 ) = rDN_DX( i, 0 );
      }

    }
    else if( dimension == 3 )
    {
      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        unsigned int index = 3 * i;

        rDeformationMatrix( 0, index + 0 ) = rDN_DX( i, 0 );
        rDeformationMatrix( 0, index + 1 ) = 0.0;
        rDeformationMatrix( 0, index + 2 ) = 0.0;

        rDeformationMatrix( 1, index + 0 ) = 0.0;
        rDeformationMatrix( 1, index + 1 ) = rDN_DX( i, 1 );
        rDeformationMatrix( 1, index + 2 ) = 0.0;

        rDeformationMatrix( 2, index + 0 ) = 0.0;
        rDeformationMatrix( 2, index + 1 ) = 0.0;
        rDeformationMatrix( 2, index + 2 ) = rDN_DX( i, 2 );

        rDeformationMatrix( 3, index + 0 ) = rDN_DX( i, 1 );
        rDeformationMatrix( 3, index + 1 ) = rDN_DX( i, 0 );
        rDeformationMatrix( 3, index + 2 ) = 0.0;

        rDeformationMatrix( 4, index + 0 ) = 0.0;
        rDeformationMatrix( 4, index + 1 ) = rDN_DX( i, 2 );
        rDeformationMatrix( 4, index + 2 ) = rDN_DX( i, 1 );

        rDeformationMatrix( 5, index + 0 ) = rDN_DX( i, 2 );
        rDeformationMatrix( 5, index + 1 ) = 0.0;
        rDeformationMatrix( 5, index + 2 ) = rDN_DX( i, 0 );

      }
    }
    else
    {
      KRATOS_ERROR << " something is wrong with the dimension when computing linear DeformationMatrix " << std::endl;
    }

  }


  /**
   * @brief Calculate Norm of stresses.
   * @param rStressVector, the stress tensor in voigt form
   * @return StressNorm, the norm of stresses
   */
  static inline double CalculateStressNorm(const Vector& rStressVector)
  {
    Matrix LocalStressTensor  = MathUtils<double>::StressVectorToTensor(rStressVector); //reduced dimension stress tensor

    Matrix StressTensor(3,3); //3D stress tensor
    noalias(StressTensor) = ZeroMatrix(3,3);
    for(unsigned int i=0; i<LocalStressTensor.size1(); i++)
    {
      for(unsigned int j=0; j<LocalStressTensor.size2(); j++)
      {
        StressTensor(i,j) = LocalStressTensor(i,j);
      }
    }

    double StressNorm =  ((StressTensor(0,0)*StressTensor(0,0))+(StressTensor(1,1)*StressTensor(1,1))+(StressTensor(2,2)*StressTensor(2,2))+
                          (StressTensor(0,1)*StressTensor(0,1))+(StressTensor(0,2)*StressTensor(0,2))+(StressTensor(1,2)*StressTensor(1,2))+
                          (StressTensor(1,0)*StressTensor(1,0))+(StressTensor(2,0)*StressTensor(2,0))+(StressTensor(2,1)*StressTensor(2,1)));

    StressNorm = sqrt(StressNorm);

    return StressNorm;

  };


  /**
   * @brief Calculate VonMises stress.
   * @param rStressVector, the stress tensor in voigt form
   * @return VonMisesStress, the von mises stress
   */
  static inline double CalculateVonMises(const Vector& rStressVector)
  {
    Matrix LocalStressTensor  = MathUtils<double>::StressVectorToTensor(rStressVector); //reduced dimension stress tensor

    Matrix StressTensor(3,3); //3D stress tensor
    noalias(StressTensor) = ZeroMatrix(3,3);
    for(unsigned int i=0; i<LocalStressTensor.size1(); i++)
    {
      for(unsigned int j=0; j<LocalStressTensor.size2(); j++)
      {
        StressTensor(i,j) = LocalStressTensor(i,j);
      }
    }

    //in general coordinates:
    double SigmaEquivalent =  (0.5)*((StressTensor(0,0)-StressTensor(1,1))*((StressTensor(0,0)-StressTensor(1,1)))+
                                     (StressTensor(1,1)-StressTensor(2,2))*((StressTensor(1,1)-StressTensor(2,2)))+
                                     (StressTensor(2,2)-StressTensor(0,0))*((StressTensor(2,2)-StressTensor(0,0)))+
                                     6*(StressTensor(0,1)*StressTensor(1,0)+StressTensor(1,2)*StressTensor(2,1)+StressTensor(2,0)*StressTensor(0,2)));

    if( SigmaEquivalent < 0 )
      SigmaEquivalent = 0;

    SigmaEquivalent = sqrt(SigmaEquivalent);

    return SigmaEquivalent;
  }

  /**
   * @brief Calculate VonMises stress.
   * @param rStressVector, the stress tensor in voigt form
   * @return VonMisesStress, the von mises stress
   */
  static inline double CalculateVonMisesUsingPrincipalStresses(const Vector& rStressVector)
  {

    //in principal stresses:

    Matrix LocalStressTensor  = MathUtils<double>::StressVectorToTensor(rStressVector); //reduced dimension stress tensor

    Matrix StressTensor(3,3); //3D stress tensor
    noalias(StressTensor) = ZeroMatrix(3,3);
    for(unsigned int i=0; i<LocalStressTensor.size1(); i++)
    {
      for(unsigned int j=0; j<LocalStressTensor.size2(); j++)
      {
        StressTensor(i,j) = LocalStressTensor(i,j);
      }
    }


    double tolerance  = 1e-10;
    double zero       = 1e-10;
    double NormStress = 0.00;

    CheckZeroDiagonalComponents (StressTensor);

    Vector PrincipalStress(3);
    noalias(PrincipalStress) = ZeroVector(3);

    NormStress =SolidMechanicsMathUtilities<double>::NormTensor(StressTensor);

    Vector MainStresses(3);
    noalias(MainStresses) = ZeroVector(3);

    bool main_tensor = CheckPrincipalStresses( StressTensor );

    if(!main_tensor)
    {

      if(NormStress>1e-6)
      {
        MainStresses = SolidMechanicsMathUtilities<double>::EigenValues(StressTensor,tolerance,zero);
      }
      else
      {
        noalias(MainStresses) = ZeroVector(3);
      }

    }
    else
    {
      noalias(MainStresses) = ZeroVector(3);
      for(unsigned int i=0; i<StressTensor.size1(); i++)
        MainStresses[i]=StressTensor(i,i);
    }


    for(unsigned int i=0; i<MainStresses.size(); i++)
      PrincipalStress[i]=MainStresses[i];


    double SigmaEquivalent =  (0.5)*((PrincipalStress[0]-PrincipalStress[1])*(PrincipalStress[0]-PrincipalStress[1]) +
                                     (PrincipalStress[1]-PrincipalStress[2])*(PrincipalStress[1]-PrincipalStress[2]) +
                                     (PrincipalStress[2]-PrincipalStress[0])*(PrincipalStress[2]-PrincipalStress[0]));


    SigmaEquivalent = sqrt(SigmaEquivalent);

    return SigmaEquivalent;
  }





 protected:

  /**
   * @brief Check and correct diagonal terms in the stress tensor
   * @param rStressTensor, the stress tensor in matrix form
   */
  static inline void  CheckZeroDiagonalComponents (Matrix& StressTensor)
  {
    // No null diagonal terms are accepted in the eigenvalue calculation
    for(unsigned int i=0; i<StressTensor.size1(); i++)
    {
      if (fabs(StressTensor(i,i))<1e-10)
      {
        StressTensor(i,i) = 1e-10;
      }
    }
  }

  /**
   * @brief Check no zero diagonal terms in the diagonalized stress tensor
   * @param rStressTensor, the stress tensor in matrix form
   * @return bool, if zero principal stresses are detected
   */
  static inline bool CheckPrincipalStresses (Matrix& StressTensor)
  {
    // No null diagonal terms are accepted in the eigenvalue calculation
    bool main = true;
    for(unsigned int i=0; i<StressTensor.size1(); i++)
    {
      for(unsigned int j=0; j<StressTensor.size2(); j++)
      {
        if(i!=j)
        {
          if (fabs(StressTensor(i,j))>1e-10)
          {
            main = false;
          }
        }
      }
    }

    return main;
  }

};
} // namespace Kratos.

#endif // KRATOS_ELELMENT_UTILITIES_H_INCLUDED
