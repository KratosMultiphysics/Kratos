//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#ifndef KRATOS_MPM_STRESS_PRINCIPAL_INVARIANTS_UTILITY
#define KRATOS_MPM_STRESS_PRINCIPAL_INVARIANTS_UTILITY

// System includes
#include <cmath>
#include <set>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{

   class MPMStressPrincipalInvariantsUtility
   {

      public:

            typedef Matrix MatrixType;

            typedef Vector VectorType;

            typedef unsigned int IndexType;

            typedef unsigned int SizeType;

            // Sort Principal Stresses, Strains, and Directions to magnitude order
            static inline void SortPrincipalStress(Vector& rPrincipalStress, Vector& rMainStrain, Matrix& rMainDirections)
            {
                  // Create Copy
                  Vector principal_direction_1 = ZeroVector(3);
                  Vector principal_direction_2 = ZeroVector(3);
                  Vector principal_direction_3 = ZeroVector(3);

                  for(unsigned int i=0; i<3; i++)
                  {
                        principal_direction_1(i) = rMainDirections(0,i);
                        principal_direction_2(i) = rMainDirections(1,i);
                        principal_direction_3(i) = rMainDirections(2,i);
                  }

                  // Reorder and swap
                  if(rPrincipalStress[0]<rPrincipalStress[1])
                  {
                        std::swap(rPrincipalStress[0],rPrincipalStress[1]);
                        std::swap(rMainStrain[0],rMainStrain[1]);
                        Vector temp_vector = principal_direction_1;
                        principal_direction_1 = principal_direction_2;
                        principal_direction_2 = temp_vector;
                  }

                  if(rPrincipalStress[1]<rPrincipalStress[2])
                  {
                        std::swap(rPrincipalStress[1],rPrincipalStress[2]);
                        std::swap(rMainStrain[1],rMainStrain[2]);
                        Vector temp_vector = principal_direction_2;
                        principal_direction_2 = principal_direction_3;
                        principal_direction_3 = temp_vector;
                  }

                  if(rPrincipalStress[0]<rPrincipalStress[1])
                  {
                        std::swap(rPrincipalStress[0],rPrincipalStress[1]);
                        std::swap(rMainStrain[0],rMainStrain[1]);
                        Vector temp_vector = principal_direction_1;
                        principal_direction_1 = principal_direction_2;
                        principal_direction_2 = temp_vector;
                  }

                  // Copy back to original matrix
                  for(unsigned int i=0; i<3; i++)
                  {
                        rMainDirections(i,0) = principal_direction_1(i);
                        rMainDirections(i,1) = principal_direction_2(i);
                        rMainDirections(i,2) = principal_direction_3(i);
                  }
            }

            // Calculate invariants I1 and J2 of a tensor
            static inline void CalculateTensorInvariants( const Vector& rVector, double& rI1, double& rJ2 )
            {
                  // Purely calculate invariant I1 = trace(rVector)
                  rI1 = 0;
                  for (unsigned int i = 0; i < 3; ++i)
                        rI1 += rVector[i];

                  // Purely calculate invariant J2
                  rJ2 = 0;
                  for (unsigned int i = 0; i < 3; i++)
                        rJ2 += std::pow( rVector[i] - rI1/3.0, 2);

                  if ( rVector.size() == 6 )
                  {
                        for (unsigned int i = 3; i < 6; i++)
                              rJ2 += 2.0 * std::pow( rVector[i], 2);
                  }
                  rJ2 /= 2.0;

            }

            // Calculate invariants I1, J2, and J3 of a tensor
            static inline void CalculateTensorInvariants( const Vector& rVector, double& rI1, double& rJ2, double& rJ3 )
            {
                  CalculateTensorInvariants( rVector, rI1, rJ2 );

                  // Purely calculate invariant J3
                  rJ3 = 0;
                  Vector s_vector = rVector;
                  for (unsigned int i = 0; i < 3; ++i)
                        s_vector[i] -= rI1/3.0;
                  if (s_vector.size() == 3)
                  {
                        rJ3 = s_vector[0] * s_vector[1] * s_vector[2]; // J_3 = det(tensor)
                  }
                  else if (s_vector.size() == 6)
                  {
                        Matrix tensor = MathUtils<double>::StressVectorToTensor( s_vector );
                        rJ3 = MathUtils<double>::Det(tensor); // J_3 = det(tensor)
                  }
            }

            // Calculate derivatives of invariants dI1/dtensor, dJ2/dtensor, and dJ3/dtensor
            static inline void CalculateTensorInvariantsDerivatives( const Vector& rVector, Vector& rDI1, Vector& rDJ2, Vector& rDJ3 )
            {
                  double i_1, j_2, j_3;
                  CalculateTensorInvariants(rVector, i_1, j_2, j_3);

                  // dI1/dtensor
                  rDI1 = ZeroVector(rVector.size());
                  for (unsigned int i = 0; i < 3; i++)
                        rDI1[i] = 1.0;

                  // dJ2/dtensor
                  rDJ2 = rVector;
                  for (unsigned int i = 0; i < 3; ++i)
                        rDJ2[i] -= i_1/3.0;

                  // dJ3/dtensor
                  Matrix s_tensor = MathUtils<double>::StressVectorToTensor( rDJ2 );
                  Matrix t_tensor = prod(s_tensor,s_tensor);
                  for (unsigned int i = 0; i < 3; ++i)
                        t_tensor(i,i) -= 2.0/3.0 * j_2;
                  rDJ3 = MathUtils<double>::StrainTensorToVector( t_tensor, rVector.size());
            }

            // Calculate second derivatives of invariants d2I1/d2tensor, d2J2/d2tensor, and d2J3/d2tensor
            static inline void CalculateTensorInvariantsSecondDerivatives( const Vector& rVector, Matrix& rD2I1, Matrix& rD2J2, Matrix& rD2J3 )
            {
                  double i_1, j_2;
                  CalculateTensorInvariants(rVector, i_1, j_2);

                  // change given vector to tensor
                  const Matrix aux_tensor = MathUtils<double>::StressVectorToTensor( rVector );

                  // d2I1/d2tensor
                  rD2I1 = ZeroMatrix(3);

                  // d2J2/d2tensor
                  rD2J2 = IdentityMatrix(3);

                  // d2J3/d2tensor
                  rD2J3  = aux_tensor;
                  for (unsigned int i = 0; i < 3; ++i)
                        rD2J3(i,i) -= i_1/3.0;
                  rD2J3 *= 2.0;

            }

            // Calculate stress invariants p (volumetric equivalent) and q (deviatoric equivalent)
            static inline void CalculateStressInvariants( const Vector& rStress, double& rMeanStressP, double& rDeviatoricQ)
            {
                  CalculateTensorInvariants(rStress, rMeanStressP, rDeviatoricQ);

                  // Volumetric Equivalent (WARNING, you might want to check the sign. P is defined positive here)
                  rMeanStressP = rMeanStressP / 3.0;

                  // Deviatoric Equivalent
                  rDeviatoricQ = std::sqrt(3 * rDeviatoricQ);
            }

            // Calculate the third stress invariants lode angle (we are using positive sine definition)
            static inline void CalculateThirdStressInvariant(const Vector& rStress, double& rLodeAngle)
            {
                  double i_1, j_2, j_3;
                  CalculateTensorInvariants( rStress, i_1, j_2, j_3);

                  // Lode Angle
                  rLodeAngle   = -j_3 / 2.0 * std::pow(3.0/j_2, 1.5);
                  double epsilon = 1.0e-9;
                  if ( std::abs(j_2) < epsilon ) {                                               // if j_2 is 0
                        rLodeAngle = GetPI() / 6.0;
                  }
                  else if ( std::abs( rLodeAngle ) > (1.0 - epsilon) ) {                        // if current rLodeAngle magnitude is larger than 1.0
                        rLodeAngle = ( GetPI() / 6.0 ) * rLodeAngle / std::abs(rLodeAngle);
                  }
                  else {                                                                        // otherwise
                        rLodeAngle = std::asin( rLodeAngle ) / 3.0;
                  }
            }

            // Calculate stress invariants p, q, and lode angle
            static inline void CalculateStressInvariants( const Vector& rStress, double& rMeanStressP, double& rDeviatoricQ, double& rLodeAngle)
            {
                  // Calculate first two stress invariant
                  CalculateStressInvariants( rStress, rMeanStressP, rDeviatoricQ);

                  // Lode Angle
                  CalculateThirdStressInvariant(rStress, rLodeAngle);
            }

            // Calculate stress invariant derivatives dp/dsigma (volumetric) and dq/dsigma (deviatoric)
            static inline void CalculateDerivativeVectors( const Vector rStress, Vector& rC1, Vector& rC2)
            {
                  // Calculate stress invariants
                  double i_1, j_2;
                  CalculateTensorInvariants( rStress, i_1, j_2);

                  // dP/dstress (WARNING, you might want to check the sign. dP is defined positive here)
                  rC1 = ZeroVector(rStress.size());
                  for (unsigned int i = 0; i < 3; i++)
                        rC1[i] = 1.0/3.0;

                  // dQ/dstress
                  rC2 = ZeroVector(rStress.size());
                  if ( std::abs(j_2) > 1E-9) {
                        for (unsigned int i = 0; i < 3; i++)
                              rC2[i] = rStress[i] - i_1/3.0;

                        rC2 *= 3.0 / (2.0 * std::sqrt(3 * j_2));
                  }

            }

            // Calculate stress invariant derivatives dp/dsigma (volumetric), dq/dsigma (deviatoric), and dlodeangle/dsigma
            static inline void CalculateDerivativeVectors( const Vector rStress, Vector& rC1, Vector& rC2, Vector& rC3)
            {
                  double i_1, j_2, j_3;
                  CalculateTensorInvariants(rStress, i_1, j_2, j_3);

                  Vector di_1, dj_2, dj_3;
                  CalculateTensorInvariantsDerivatives(rStress, di_1, dj_2, dj_3);

                  // Compute dP/dstress and dQ/dstress
                  CalculateDerivativeVectors(rStress, rC1, rC2);

                  // Compute dLodeAngle/dstress
                  double lode_angle;
                  CalculateThirdStressInvariant(rStress, lode_angle);

                  // Compute dLodeAngle/dstress
                  rC3  = dj_3 - (3.0/2.0 * j_3/j_2) * dj_2;
                  rC3 *= - std::sqrt(3.0) / (2.0 * std::cos(3*lode_angle) * std::pow(j_2, 1.5));
            }

            // Calculate second stress invariant derivatives d2p/d2sigma (volumetric) and d2q/d2sigma (deviatoric)
            static inline void CalculateSecondDerivativeMatrices( const Vector rStress, Matrix& r2C1, Matrix& r2C2)
            {
                  // Calculate stress invariants
                  double i_1, j_2;
                  CalculateTensorInvariants( rStress, i_1, j_2);

                  // d2P/d2stress
                  r2C1 = ZeroMatrix(3);

                  // d2Q/d2stress
                  r2C2 = IdentityMatrix(3);
                  if ( std::abs(j_2) > 1E-9) r2C2 *= 3.0 / (2.0 * std::sqrt(3 * j_2));

            }

            // Calculate second stress invariant derivatives d2p/d2sigma (volumetric), d2q/d2sigma (deviatoric), and d2lodeangle/d2sigma
            static inline void CalculateSecondDerivativeMatrices( const Vector rStress, Matrix& r2C1, Matrix& r2C2, Matrix& r2C3)
            {
                  double i_1, j_2, j_3;
                  CalculateTensorInvariants(rStress, i_1, j_2, j_3);

                  Matrix d2i_1, d2j_2, d2j_3;
                  CalculateTensorInvariantsSecondDerivatives(rStress, d2i_1, d2j_2, d2j_3);

                  // Compute d2P/d2stress and d2Q/d2stress
                  CalculateSecondDerivativeMatrices(rStress, r2C1, r2C2);

                  // Compute dLodeAngle/dstress
                  double lode_angle;
                  CalculateThirdStressInvariant(rStress, lode_angle);

                  // Compute dLodeAngle/dstress
                  r2C3  = d2j_3 - (3.0/2.0 * j_3/j_2) * d2j_2;
                  r2C3 *= - std::sqrt(3.0) / (2.0 * std::cos(3*lode_angle) * std::pow(j_2, 1.5));
            }

            static double GetPI()
            {
                  return std::atan(1.0)*4.0;
            }

   }; // end Class MPMStressPrincipalInvariantsUtility

} // end namespace Kratos

#endif // KRATOS_MPM_STRESS_PRINCIPAL_INVARIANTS_UTILITY