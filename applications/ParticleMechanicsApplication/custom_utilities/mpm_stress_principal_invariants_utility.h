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

            // Tolerance
            static constexpr double tolerance = 1.0e-9;


            /**
             * @brief Sort Principal Stresses, Strains, and Directions to magnitude order.
             * @param[in/out] rPrincipalStress The principal stresses to be sorted in descending magnitude order.
             * @param[in/out] rMainStrain The principal strain corresponding to the stresses.
             * @param[in/out] rMainDirections The eigen directions corresponding to the stresses.
             */
            static inline void SortPrincipalStress(Vector& rPrincipalStress, Vector& rMainStrain, Matrix& rMainDirections)
            {
                  // Create Copy
                  Vector principal_direction_1 = ZeroVector(3);
                  Vector principal_direction_2 = ZeroVector(3);
                  Vector principal_direction_3 = ZeroVector(3);

                  for(unsigned int i=0; i<3; i++)
                  {
                        principal_direction_1[i] = rMainDirections(0,i);
                        principal_direction_2[i] = rMainDirections(1,i);
                        principal_direction_3[i] = rMainDirections(2,i);
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
                        rMainDirections(i,0) = principal_direction_1[i];
                        rMainDirections(i,1) = principal_direction_2[i];
                        rMainDirections(i,2) = principal_direction_3[i];
                  }
            }


            /**
             * @brief Calculate invariants I1 and J2 of a tensor.
             * @param[in] rVector Input principal tensor.
             * @param[out] rI1 First stress tensor invariant.
             * @param[out] rJ2 Second stress deviator tensor invariant.
             */
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


            /**
             * @brief Calculate invariants I1, J2, and J3 of a tensor.
             * @param[in] rVector Input principal tensor.
             * @param[out] rI1 First stress tensor invariant.
             * @param[out] rJ2 Second stress deviator tensor invariant.
             * @param[out] rJ3 Third stress deviator tensor invariant.
             */
            static inline void CalculateTensorInvariants( const Vector& rVector, double& rI1, double& rJ2, double& rJ3 )
            {
                  CalculateTensorInvariants( rVector, rI1, rJ2 );

                  // Purely calculate invariant J3
                  rJ3 = 0;
                  Vector s_vector = rVector;
                  for (unsigned int i = 0; i < 3; ++i)
                        s_vector[i] -= rI1/3.0;
                  Matrix tensor = PrincipalVectorToMatrix( s_vector , rVector.size());
                  rJ3 = MathUtils<double>::Det(tensor); // J_3 = det(tensor)
            }


            /**
             * @brief Calculate derivatives of invariants dI1/dtensor, dJ2/dtensor, and dJ3/dtensor.
             * @param[in] rVector Input principal tensor.
             * @param[out] rDI1 First derivative of first stress tensor invariant.
             * @param[out] rDJ2 First derivative of second stress deviator tensor invariant.
             * @param[out] rDJ3 First derivative of third stress deviator tensor invariant.
             */
            static inline void CalculateTensorInvariantsDerivatives( const Vector& rVector, Vector& rDI1, Vector& rDJ2, Vector& rDJ3 )
            {
                  double i_1, j_2, j_3;
                  CalculateTensorInvariants(rVector, i_1, j_2, j_3);

                  // dI1/dtensor
                  rDI1 = ZeroVector(rVector.size());
                  for (unsigned int i = 0; i < 3; ++i)
                        rDI1[i] = 1.0;

                  // dJ2/dtensor
                  rDJ2 = ZeroVector(rVector.size());
                  rDJ2 = rVector;
                  for (unsigned int i = 0; i < 3; ++i)
                        rDJ2[i] -= i_1/3.0;

                  // dJ3/dtensor
                  rDJ3 = ZeroVector(rVector.size());
                  Matrix s_tensor = PrincipalVectorToMatrix( rDJ2, rVector.size() );
                  Matrix t_tensor = prod(s_tensor,s_tensor);
                  for (unsigned int i = 0; i < 3; ++i)
                        t_tensor(i,i) -= 2.0/3.0 * j_2;
                  rDJ3 = PrincipalMatrixtoVector( t_tensor, rVector.size());

            }


            /**
             * @brief Calculate second derivatives of invariants d2I1/d2tensor, d2J2/d2tensor, and d2J3/d2tensor.
             * @param[in] rVector Input principal tensor.
             * @param[out] rD2I1 Second derivative of first stress tensor invariant.
             * @param[out] rD2J2 Second derivative of second stress deviator tensor invariant.
             * @param[out] rD2J3 Second derivative of third stress deviator tensor invariant.
             */
            static inline void CalculateTensorInvariantsSecondDerivatives( const Vector& rVector, Matrix& rD2I1, Matrix& rD2J2, Matrix& rD2J3 )
            {
                  // This function should return a forth-order tensor and thus the current usage is limited only to rVector.size() = 3
                  if (rVector.size() == 3)
                  {
                        double i_1, j_2;
                        CalculateTensorInvariants(rVector, i_1, j_2);

                        // d2I1/d2tensor
                        rD2I1 = ZeroMatrix(3,3);

                        // d2J2/d2tensor = delta_mp * delta_nq - 1/3 delta_pq * delta_mn
                        rD2J2 = ZeroMatrix(3,3);
                        for (unsigned int i = 0; i < 3; ++i)
                              for (unsigned int j = 0; j < 3; ++j)
                              {
                                    if (i == j) rD2J2(i,j) = 2.0/3.0;
                                    else rD2J2(i,j) = -1.0/3.0;
                              }

                        // d2J3/d2tensor
                        rD2J3 = ZeroMatrix(3,3);
                        Vector s_vector = rVector;
                        for (unsigned int i = 0; i < 3; ++i)
                              s_vector[i] -= i_1/3.0;

                        for (unsigned int i = 0; i < 3; ++i)
                              for (unsigned int j = 0; j < 3; ++j)
                              {
                                    if (i == j) rD2J3(i,j) = 2.0/3.0 * s_vector[i];
                                    else rD2J3(i,j) = - 2.0/3.0* (s_vector[i] + s_vector[j]);
                              }
                  }
                  else
                  {
                        KRATOS_ERROR <<  "The given vector dimension is wrong! Given: " << rVector.size() << "; Expected: 3" << std::endl;
                  }

            }


            /**
             * @brief Calculate stress invariants p (volumetric equivalent) and q (deviatoric equivalent).
             * @param[in] rStress Input principal stress tensor.
             * @param[out] rMeanStressP Hydrostatic mean stress.
             * @param[out] rDeviatoricQ Deviatoric stress component.
             */
            static inline void CalculateStressInvariants( const Vector& rStress, double& rMeanStressP, double& rDeviatoricQ)
            {
                  CalculateTensorInvariants(rStress, rMeanStressP, rDeviatoricQ);

                  // Volumetric Equivalent (WARNING, you might want to check the sign. P is defined positive here)
                  rMeanStressP = rMeanStressP / 3.0;

                  // Deviatoric Equivalent
                  rDeviatoricQ = std::sqrt(3 * rDeviatoricQ);
            }


            /**
             * @brief Calculate the third stress invariants lode angle (we are using positive sine definition).
             * @param[in] rStress Input principal tensor.
             * @param[out] rLodeAngle Third stress invariant direction used for non-circular octahedral profile yield surface.
             */
            static inline void CalculateThirdStressInvariant(const Vector& rStress, double& rLodeAngle)
            {
                  double i_1, j_2, j_3;
                  CalculateTensorInvariants( rStress, i_1, j_2, j_3);

                  // Lode Angle
                  if ( std::abs(j_2) < tolerance ) // if j_2 is 0
                        j_2 = tolerance;

                  rLodeAngle   = -j_3 / 2.0 * std::pow(3.0/j_2, 1.5);
                  if ( std::abs( rLodeAngle ) > 1.0 )                                           // if current rLodeAngle magnitude is larger than 1.0
                        rLodeAngle = ( Globals::Pi / 6.0 ) * rLodeAngle / std::abs(rLodeAngle);
                  else                                                                          // otherwise
                        rLodeAngle = std::asin( rLodeAngle ) / 3.0;

            }


            /**
             * @brief Calculate stress invariants p, q, and lode angle.
             * @param[in] rStress Input principal stress tensor.
             * @param[out] rMeanStressP Hydrostatic mean stress.
             * @param[out] rDeviatoricQ Deviatoric stress component.
             * @param[out] rLodeAngle Third stress invariant direction used for non-circular octahedral profile yield surface.
             */
            static inline void CalculateStressInvariants( const Vector& rStress, double& rMeanStressP, double& rDeviatoricQ, double& rLodeAngle)
            {
                  // Calculate first two stress invariant
                  CalculateStressInvariants( rStress, rMeanStressP, rDeviatoricQ);

                  // Lode Angle
                  CalculateThirdStressInvariant(rStress, rLodeAngle);
            }


            /**
             * @brief Calculate stress invariant derivatives dp/dsigma (volumetric) and dq/dsigma (deviatoric).
             * @param[in] rStress Input principal stress tensor.
             * @param[out] rC1 Hydrostatic mean stress derivative with respect to principal stresses.
             * @param[out] rC2 Deviatoric stress derivative with respect to principal stresses.
             */
            static inline void CalculateDerivativeVectors( const Vector rStress, Vector& rC1, Vector& rC2)
            {
                  // Calculate stress invariants
                  double p, q;
                  CalculateStressInvariants( rStress, p, q);

                  // dP/dstress (WARNING, you might want to check the sign. dP is defined positive here)
                  rC1 = ZeroVector(rStress.size());
                  for (unsigned int i = 0; i < 3; i++)
                        rC1[i] = 1.0/3.0;

                  // dQ/dstress
                  rC2 = ZeroVector(rStress.size());
                  if ( std::abs(q) > tolerance) {
                        rC2 = rStress;
                        for (unsigned int i = 0; i < 3; i++)
                              rC2[i] -= p;

                        rC2 *= 3.0 / (2.0 * q);
                  }

            }


            /**
             * @brief Calculate stress invariant derivatives dp/dsigma (volumetric), dq/dsigma (deviatoric), and dlodeangle/dsigma.
             * @param[in] rStress Input principal stress tensor.
             * @param[out] rC1 Hydrostatic mean stress derivative with respect to principal stresses.
             * @param[out] rC2 Deviatoric stress derivative with respect to principal stresses.
             * @param[out] rC3 Derivative of lode angle with respect to principal stresses.
             */
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
                  rC3  = ZeroVector(rStress.size());
                  if (std::abs(j_2) > tolerance)
                  {
                        rC3  = dj_3 - (3.0/2.0 * j_3/j_2) * dj_2;
                        rC3 *= - std::sqrt(3.0) / (2.0 * std::cos(3*lode_angle) * std::pow(j_2, 1.5));
                  }
            }


            /**
             * @brief Calculate second stress invariant derivatives d2p/d2sigma (volumetric) and d2q/d2sigma (deviatoric).
             * @param[in] rStress Input principal stress tensor.
             * @param[out] r2C1 Hydrostatic mean stress second derivative with respect to principal stresses.
             * @param[out] r2C2 Deviatoric stress second derivative with respect to principal stresses.
             */
            static inline void CalculateSecondDerivativeMatrices( const Vector rStress, Matrix& r2C1, Matrix& r2C2)
            {
                  // This function should return a forth-order tensor and thus the current usage is limited only to rStress.size() = 3
                  if (rStress.size() == 3)
                  {
                        // Calculate stress invariants
                        double p, q;
                        CalculateStressInvariants( rStress, p, q);

                        // d2P/d2stress
                        r2C1 = ZeroMatrix(3,3);

                        // d2Q/d2stress
                        r2C2 = ZeroMatrix(3,3);
                        if (std::abs(q) > tolerance)
                        {
                              Vector s_vector = rStress;
                              for (unsigned int i = 0; i < 3; ++i)
                                    s_vector[i] -= p;

                              for (unsigned int i = 0; i < 3; ++i)
                                    for (unsigned int j = 0; j < 3; ++j)
                                    {
                                          if (i == j) r2C2(i,j) = 1.0 / q;
                                          else r2C2(i,j) = - 1.0/2.0/q;

                                          r2C2(i,j) -= (9.0/4.0/std::pow(q, 3) * s_vector[i] * s_vector[j]);
                                    }
                        }

                  }
                  else
                  {
                        KRATOS_ERROR <<  "The given vector dimension is wrong! Given: " << rStress.size() << "; Expected: 3" << std::endl;
                  }

            }


            /**
             * @brief Calculate second stress invariant derivatives d2p/d2sigma (volumetric), d2q/d2sigma (deviatoric), and d2lodeangle/d2sigma.
             * @param[in] rStress Input principal stress tensor.
             * @param[out] r2C1 Hydrostatic mean stress second derivative with respect to principal stresses.
             * @param[out] r2C2 Deviatoric stress second derivative with respect to principal stresses.
             * @param[out] r2C3 Second Derivative of lode angle with respect to principal stresses.
             */
            static inline void CalculateSecondDerivativeMatrices( const Vector rStress, Matrix& r2C1, Matrix& r2C2, Matrix& r2C3)
            {
                  // This function should return a forth-order tensor and thus the current usage is limited only to rStress.size() = 3
                  if (rStress.size() == 3)
                  {
                        // Compute d2P/d2stress and d2Q/d2stress
                        CalculateSecondDerivativeMatrices(rStress, r2C1, r2C2);

                        // Compute d2LodeAngle/d2stress
                        double i_1, j_2, j_3;
                        CalculateTensorInvariants(rStress, i_1, j_2, j_3);

                        Vector di_1, dj_2, dj_3;
                        CalculateTensorInvariantsDerivatives(rStress, di_1, dj_2, dj_3);

                        Matrix d2i_1, d2j_2, d2j_3;
                        CalculateTensorInvariantsSecondDerivatives(rStress, d2i_1, d2j_2, d2j_3);

                        Vector dp, dq, dlode_angle;
                        CalculateDerivativeVectors(rStress, dp, dq, dlode_angle);

                        double lode_angle;
                        CalculateThirdStressInvariant(rStress, lode_angle);

                        r2C3 = ZeroMatrix(3,3);
                        if (std::abs(j_2) > tolerance)
                        {
                              if (std::abs(j_3) < tolerance) j_3 = tolerance;

                              // To reduce complication, we simplify into two components: d(LodeAngle)/dstress = AT * J3' + AS * J2'
                              const double AT = - std::sqrt(3) / (2.0 * j_2 * std::sqrt(j_2) * std::cos(3 * lode_angle));
                              const double AS = AT * (- 1.50 * j_3 / j_2);

                              // Then we need to compute the derivative of each component with respect to stress
                              // 1. d(AT * J3')/dstress = dAT * J3' + AT * J3''
                              Matrix first_component = ZeroMatrix(3,3);
                              Vector aux_dAT = ZeroVector(3);
                              aux_dAT = AT * (-3.0/2.0/j_2 * dj_2 + 3.0 * std::tan(3*lode_angle) * dlode_angle);
                              first_component = MathUtils<double>::TensorProduct3(aux_dAT, dj_3);
                              first_component += AT * d2j_3;

                              // 2. d(AS * J2')/dstress = dAS * J2' + AS * J2''
                              Matrix second_component = ZeroMatrix(3,3);
                              Vector aux_dAS = ZeroVector(3);
                              aux_dAS = AS * (-5.0/2.0/j_2 * dj_2 + 3.0 * std::tan(3*lode_angle) * dlode_angle + 1.0/j_3 * dj_3);
                              second_component = MathUtils<double>::TensorProduct3(aux_dAS, dj_2);
                              second_component += AS * d2j_2;

                              // Sum the two components up
                              r2C3 = first_component + second_component;
                        }
                  }
                  else
                  {
                        KRATOS_ERROR <<  "The given vector dimension is wrong! Given: " << rStress.size() << "; Expected: 3" << std::endl;
                  }
            }


            /**
             * @brief Transform Stress or Principal stress tensor to Matrix form (fully assuming 3D space).
             * @param[in] rVector Vector tensor to be transformed.
             * @param[in] rSize Matrix size, 3 for only principal stresses, 6 for normal stresses.
             * @return transformed matrix.
             */
            static Matrix PrincipalVectorToMatrix(const Vector& rVector, const unsigned int& rSize)
            {
                  Matrix matrix = ZeroMatrix(3,3);
                  if (rSize == 3)
                  {
                        for(unsigned int i=0; i<3; ++i)
                              matrix(i,i) = rVector[i];
                  }
                  else if (rSize == 6)
                  {
                        matrix = MathUtils<double>::StressVectorToTensor( rVector );
                  }
                  return matrix;
            }


            /**
             * @brief Transform Stress or Principal stress tensor to Vector form (fully assuming 3D space).
             * @param[in] rMatrix Matrix tensor to be transformed.
             * @param[in] rSize Matrix size, 3 for only principal stresses, 6 for normal stresses.
             * @return transformed vector.
             */
            static Vector PrincipalMatrixtoVector(const Matrix& rMatrix, const unsigned int& rSize)
            {
                  Vector vector = ZeroVector(3);
                  if (rSize == 3)
                  {
                        for(unsigned int i=0; i<3; ++i)
                              vector[i] = rMatrix(i,i);
                  }
                  else if (rSize == 6)
                  {
                        vector = MathUtils<double>::StressTensorToVector( rMatrix, rSize );
                  }
                  return vector;
            }


            static double CalculateMatrixDoubleContraction(const Matrix& rInput)
            {
                double result = 0.0;
                KRATOS_ERROR_IF(rInput.size1() != rInput.size2()) 
                    << "CalculateMatrixDoubleContraction only takes square matrices\n";
                for (size_t i = 0; i < rInput.size1(); ++i)
                {
                    for (size_t j = 0; j < rInput.size2(); ++j)
                    {
                        result += rInput(i, j) * rInput(i, j);
                    }
                }
                return result;
            }


            static double CalculateMatrixTrace(const Matrix& rInput)
            {
                KRATOS_ERROR_IF(rInput.size1() != rInput.size2()) << "Can only calculate trace of square matrices";
                double trace = 0.0;
                for (size_t i = 0; i < rInput.size1(); ++i) trace += rInput(i, i);
                return trace;
            }
   }; // end Class MPMStressPrincipalInvariantsUtility

} // end namespace Kratos

#endif // KRATOS_MPM_STRESS_PRINCIPAL_INVARIANTS_UTILITY