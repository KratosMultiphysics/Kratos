//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                 LlMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_STRESS_INVARIANTS_UTILITIES)
#define KRATOS_STRESS_INVARIANTS_UTILITIES


// System includes
#include <cmath>
#include <set>


// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"

#include "constitutive_model_utilities.hpp"

namespace Kratos
{


   class StressInvariantsUtilities
   {

      public:

         typedef BoundedMatrix<double,3,3> MatrixType;

         typedef array_1d<double, 6> VectorType;

         typedef unsigned int IndexType;

         typedef unsigned int SizeType;

         static inline void CalculateStressInvariants( const MatrixType& rStressMatrix, double & rI1, double & rJ2)
         {
            Vector StressVector = ZeroVector(6);
            ConstitutiveModelUtilities::StressTensorToVector( rStressMatrix, StressVector);
            CalculateStressInvariants( StressVector, rI1, rJ2);
         }

         static inline void CalculateStressInvariants( const Vector& rStress, double& I1, double& J2)
         {

            I1 = 0;
            for (unsigned int i = 0; i < 3; ++i)
               I1 += rStress(i);
            I1 /= 3.0;

            J2 = 0;
            for (unsigned int i = 0; i < 3; i++)
               J2 += pow( rStress(i) - I1, 2);

            for (unsigned int i = 3; i < 6; i++)
               J2 += 2.0*pow( rStress(i), 2);

            J2 = sqrt( 0.5 * J2 );

         }

         static inline void CalculateStressInvariants( const MatrixType& rStressMatrix, double & rI1, double & rJ2, double & rLode)
         {
            Vector StressVector = ZeroVector(6);
            ConstitutiveModelUtilities::StressTensorToVector( rStressMatrix, StressVector);
            CalculateStressInvariants( StressVector, rI1, rJ2, rLode);
         }
         static inline void CalculateStressInvariants( const Vector& rStress, double& I1, double& J2, double& Lode)
         {
            CalculateStressInvariants( rStress, I1, J2);
            Matrix StressTensor  = MathUtils<double>::StressVectorToTensor( rStress);

            for (unsigned int i = 0; i < 3; ++i)
               StressTensor(i,i) -= I1;

            Lode = MathUtils<double>::Det(StressTensor);
            Lode = 3.0 * sqrt(3.0) / 2.0 * Lode / pow( J2, 3);

            double epsi = 1.0e-9;
            if ( fabs( Lode ) > 1.0-epsi) {
               Lode = -30.0*Globals::Pi / 180.0 * Lode / fabs(Lode);
            }
            else if ( J2 < 10.0*epsi) {
               Lode = 30.0*Globals::Pi / 180.0;
            }
            else {
               Lode = std::asin( -Lode) / 3.0;
            }

            //std::cout << " THISLODE " << Lode << " Stress;atrox " << StressTensor << " Stress " << rStress <<
         }

         static inline void CalculateDerivativeVectors( const MatrixType& rStressMatrix, VectorType& C1, VectorType & C2)
         {
            Vector StressVector = ZeroVector(6);
            ConstitutiveModelUtilities::StressTensorToVector( rStressMatrix, StressVector);
            CalculateDerivativeVectors( StressVector, C1, C2);
         }

         static inline void CalculateDerivativeVectors( const Vector & rStress, VectorType& C1, VectorType& C2)
         {

            C1.clear();
            for (unsigned int i = 0; i < 3; i++)
               C1(i) = 1.0/3.0;

            double I1, J2;
            CalculateStressInvariants( rStress, I1, J2);

            C2.clear();

            for (unsigned int i = 0; i < 3; i++)
               C2(i) = rStress(i) - I1;

            for (unsigned int i = 3; i < 6; i++)
               C2(i) = 2.0*rStress(i);

            if ( J2 > 1E-5) {
               C2 /= 2.0 * J2;
            }
            else {
               C2 = ZeroVector(6);
            }

         }

         static inline void CalculateDerivativeVectors( const Vector rStress, Vector& C1, Vector& C2, Vector& C3)
         {
            // COPY
            C1 = ZeroVector(6);
            for (unsigned int i = 0; i < 3; i++)
               C1(i) = 1.0/3.0;

            double I1, J2;
            CalculateStressInvariants( rStress, I1, J2);

            C2 = ZeroVector(6);

            for (unsigned int i = 0; i < 3; i++)
               C2(i) = rStress(i) - I1;

            for (unsigned int i = 3; i < 6; i++)
               C2(i) = 2.0*rStress(i);

            if ( J2 > 1E-5) {
               C2 /= 2.0 * J2;
            }
            else {
               C2 = ZeroVector(6);
               C3 = ZeroVector(6);
               return;
            }

            C3 = ZeroVector(6);

            Vector ShearStress = rStress;
            for (int i = 0; i < 3; i++)
               ShearStress(i) -= I1;

            C3(0) = ShearStress(1)*ShearStress(2) - pow( ShearStress(4), 2);
            C3(1) = ShearStress(2)*ShearStress(0) - pow( ShearStress(5), 2);
            C3(2) = ShearStress(0)*ShearStress(1) - pow( ShearStress(3), 2);

            C3(3) = 2.0 * ( ShearStress(4)*ShearStress(5) - ShearStress(2)*ShearStress(3));
            C3(4) = 2.0 * ( ShearStress(5)*ShearStress(3) - ShearStress(0)*ShearStress(4));
            C3(5) = 2.0 * ( ShearStress(3)*ShearStress(4) - ShearStress(1)*ShearStress(5));

            for (unsigned int i = 0; i < 3; ++i)
               C3(i) += pow(J2, 2) / 3.0;

            Matrix Aux, ShearStressM;
            ShearStressM = MathUtils<double>::StressVectorToTensor( ShearStress);
            Aux = prod(  ShearStressM, ShearStressM);

            for (unsigned int i = 0; i < 3; i++)
               Aux(i,i) -= 1.0/3.0 * 2.0 * pow( J2 , 2);

            C3 = MathUtils<double>::StrainTensorToVector( Aux, 6);


         }

   }; // end Class StressInvariantsUtilities

} // end namespace Kratos

#endif // KRATOS_STRESS_INVARIANTS_UTILITIES
