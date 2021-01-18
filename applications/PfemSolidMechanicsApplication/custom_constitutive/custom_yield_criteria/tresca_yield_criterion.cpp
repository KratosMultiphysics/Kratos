//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_constitutive/custom_yield_criteria/tresca_yield_criterion.hpp"

#include "pfem_solid_mechanics_application_variables.h"

// ROUNDED TRESCA YIELD CRITERION (Sloan & Booker, 1986 )

// This is the C2 approximation ... (Abbo, Lyamin, Sloan, Hambleton,) 
// Maybe I can convert it to Mohr Coulomb (the current one may not work)
// the functions of invariants must disappear

namespace Kratos
{


   //*******************************CONSTRUCTOR******************************************
   //************************************************************************************
   TrescaYieldCriterion::TrescaYieldCriterion()
      :YieldCriterion()
   {

   }

   //*****************************INITIALIZATION CONSTRUCTOR*****************************
   //************************************************************************************

   TrescaYieldCriterion::TrescaYieldCriterion(HardeningLawPointer pHardeningLaw)
      :YieldCriterion(pHardeningLaw)
   {

   }

   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   TrescaYieldCriterion& TrescaYieldCriterion::operator=(TrescaYieldCriterion const& rOther)
   {
      YieldCriterion::operator=(rOther);
      return *this;
   }

   //*******************************COPY CONSTRUCTOR*************************************
   //************************************************************************************

   TrescaYieldCriterion::TrescaYieldCriterion(TrescaYieldCriterion const& rOther)
      :YieldCriterion(rOther)
   {

   }


   //********************************DESTRUCTOR******************************************
   //************************************************************************************

   TrescaYieldCriterion::~TrescaYieldCriterion()
   {
   }



   //************************* CALCULATE YIELD FUNCTION  ******************
   //**********************************************************************

   double& TrescaYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const double& rAlpha)
   {
      // That is, undrained shear strength
      double YieldStress = this->GetHardeningLaw().GetProperties()[YIELD_STRESS];
      TrescaStressInvariants StressInvariants;
      this->CalculateStressInvariants(rStressVector, StressInvariants );


      double LodeCut = this->GetSmoothingLodeAngle();

      if ( fabs(StressInvariants.LodeAngle) < LodeCut ) {

         rStateFunction = StressInvariants.J2InvSQ * std::cos( StressInvariants.LodeAngle );

      }
      else {

         this->CalculateSmoothingInvariants( StressInvariants);

         rStateFunction = StressInvariants.J2InvSQ * ( StressInvariants.A +  StressInvariants.B * std::sin( 3.0* StressInvariants.LodeAngle )  + StressInvariants.C * pow( std::sin( 3.0*StressInvariants.LodeAngle), 2 )  );
      }

      //rStateFunction = rStateFunction/YieldStress - 1.0;
      rStateFunction = rStateFunction - 1.0*YieldStress;
      return rStateFunction; 

   }

   //************************* YIELD FUNCTION DERIVATIVE ******************
   //**********************************************************************
   void TrescaYieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rYieldFunctionD, const double& rAlpha)
   {
      //double YieldStress = this->GetHardeningLaw().GetProperties()[YIELD_STRESS];

      TrescaStressInvariants StressInvariants;
      this->CalculateStressInvariants(rStressVector, StressInvariants );


      double LodeCut = this->GetSmoothingLodeAngle();

      double C2;
      double C3;

      if ( fabs(StressInvariants.LodeAngle) < LodeCut ) {

         C2 =  std::cos( StressInvariants.LodeAngle) * ( 1.0 - std::tan( StressInvariants.LodeAngle) * std::tan( 3.0 * StressInvariants.LodeAngle) );
         // WRONG??
         C2 =  std::cos( StressInvariants.LodeAngle) * ( 1.0 + std::tan( StressInvariants.LodeAngle) * std::tan( 3.0 * StressInvariants.LodeAngle) );

         C3 = -sqrt( 3.0)/2.0 * std::sin( StressInvariants.LodeAngle) / std::cos( 3.0* StressInvariants.LodeAngle) ;  // also wrong?
         C3 /= pow(StressInvariants.J2InvSQ, 2);

         C3 = sqrt(3.0) / 2.0 * std::sin( StressInvariants.LodeAngle) / std::cos( 3.0* StressInvariants.LodeAngle);
         C3 /= pow( StressInvariants.J2InvSQ, 2);

         if ( fabs( StressInvariants.J2InvSQ) < 1E-6)
            C3 = 0.0;
      }
      else {

         this->CalculateSmoothingInvariants( StressInvariants);

         C2 = StressInvariants.A + 2.0 * StressInvariants.B * std::sin ( 3.0 * StressInvariants.LodeAngle) ;

         C3 = 3.0 * sqrt(3.0)  / 2.0  * StressInvariants.B / pow(StressInvariants.J2InvSQ, 2);

         // ARA EL C2Continuous

         C2 = StressInvariants.A - 2.0 * StressInvariants.B * std::sin( 3.0 * StressInvariants.LodeAngle ) - 5.0 * StressInvariants.C * pow( std::sin( 3.0 * StressInvariants.LodeAngle ) , 2 );

         C3 = StressInvariants.B + 2.0 * StressInvariants.C * std::sin( 3.0 * StressInvariants.LodeAngle ) ;
         C3 *= -3.0 * sqrt( 3.0) / ( 2.0 * pow( StressInvariants.J2InvSQ, 2 ) );

         if ( fabs( StressInvariants.J2InvSQ) < 1E-6) {
            C3 = 0.0;
         }

      } 

      Vector C2Vector, C3Vector; 
      ComputeC2andC3Vector( rStressVector, StressInvariants, C2Vector, C3Vector);

      rYieldFunctionD = C2*C2Vector + C3*C3Vector;

     return;  
      std::cout << " " << std::endl;
      
      Vector  StressVector = ZeroVector(6);
      StressVector(0) = 1.0;
      StressVector(1) = 2.0;
      StressVector(2) = 3.0;

      this->CalculateStressInvariants( StressVector, StressInvariants); 
      ComputeC2andC3Vector( StressVector, StressInvariants, C2Vector, C3Vector);
      std::cout << " StressInv: p " << StressInvariants.MeanStress <<  " J2 " << StressInvariants.J2InvSQ << " LODE " << StressInvariants.LodeAngle << std::endl;
      std::cout << " C2 Vector " << C2Vector << std::endl;
      std::cout << " C3 Vector " << C3Vector << std::endl;

   }


   // YIELD FUNCION SECOND DERIVATIVE
/*   void TrescaYieldCriterion::CalculateYieldFunctionSecondDerivative(const Vector& rStressVector, Matrix& rYieldFunctionDD, const double& rAlpha)
   {
      rYieldFunctionDD = ZeroMatrix(1,1);

      TrescaStressInvariants StressInvariants;
      this->CalculateStressInvariants(rStressVector, StressInvariants );
      //std::cout << " LODE ANGLE " << StressInvariants.LodeAngle * 180.0 / GetPI() << std::endl;
      // this part is a copy

      double LodeCut = this->GetSmoothingLodeAngle();

      double C2;
      double C3;

      if ( fabs(StressInvariants.LodeAngle) < LodeCut ) {

         C2 =  std::cos( StressInvariants.LodeAngle) * ( 1.0 - std::tan( StressInvariants.LodeAngle) * std::tan( 3.0 * StressInvariants.LodeAngle) );
         // WRONG??
         C2 =  std::cos( StressInvariants.LodeAngle) * ( 1.0 + std::tan( StressInvariants.LodeAngle) * std::tan( 3.0 * StressInvariants.LodeAngle) );

         C3 = -sqrt(3.0) /2.0 * std::sin( StressInvariants.LodeAngle) / std::cos( 3.0* StressInvariants.LodeAngle) ;  // also wrong?
         C3 /= pow(StressInvariants.J2InvSQ, 2);

         C3 = sqrt(3.0) / 2.0 * std::sin( StressInvariants.LodeAngle) / std::cos( 3.0* StressInvariants.LodeAngle);
         C3 /= pow( StressInvariants.J2InvSQ, 2);

         if ( fabs( StressInvariants.J2InvSQ) < 1E-6)
            C3 = 0.0;
      }
      else {

         this->CalculateSmoothingInvariants( StressInvariants);

         C2 = StressInvariants.A + 2.0 * StressInvariants.B * std::sin ( 3.0 * StressInvariants.LodeAngle) ;

         C3 = 3.0 * sqrt(3.0)  / 2.0  * StressInvariants.B / pow(StressInvariants.J2InvSQ, 2);

         // ARA EL C2Continuous

         C2 = StressInvariants.A - 2.0 * StressInvariants.B * std::sin( 3.0 * StressInvariants.LodeAngle ) - 5.0 * StressInvariants.C * pow( std::sin( 3.0 * StressInvariants.LodeAngle ) , 2 );

         C3 = StressInvariants.B + 2.0 * StressInvariants.C * std::sin( 3.0 * StressInvariants.LodeAngle ) ;
         C3 *= -3.0 * sqrt( 3.0) / ( 2.0 * pow( StressInvariants.J2InvSQ, 2 ) );

         if ( fabs( StressInvariants.J2InvSQ) < 1E-6) {
            C3 = 0.0;
         }
      } 

      // here begins new things
      Vector C2Vector, C3Vector, LodeDerivative;
      Matrix C2Matrix, C3Matrix; 
      ComputeC2andC3VectorDD( rStressVector, StressInvariants, C2Vector, C3Vector, C2Matrix, C3Matrix, LodeDerivative);

      Vector C2V, C3V;
      if ( fabs( StressInvariants.LodeAngle) < LodeCut)
      {
         C2V = LodeDerivative;
         double LD = StressInvariants.LodeAngle;

         //C2V *= -std::sin(LD)  + std::cos( LD) * std::tan(3.0*LD) + 3.0 * std::sin( LD) / pow( cos( 3.0* LD) , 2.0) ; // wrong
         C2V *= -std::sin(LD) - std::sin(LD)*std::tan(LD)*std::tan(3.0*LD)  + std::tan(3.0*LD) / std::cos(LD) + 3.0 * std::sin( LD) / pow( cos( 3.0* LD) , 2) ; 


         // MALAMENT CASI SEGUR
         Vector Aux1, Aux2; 
         Aux1 = LodeDerivative; 
         Aux2 = C2Vector;

         Aux1 = C2Vector; 
         Aux1 *= -2.0 / StressInvariants.J2InvSQ * std::sin( LD);

         Aux2 = LodeDerivative; 
         Aux2 *= + std::cos( LD )  + 3.0 * std::sin ( LD) * std::tan( 3.0*LD) ;

         C3V = Aux1 + Aux2;
         C3V *=  + sqrt(3)  / 2.0 / pow ( StressInvariants.J2InvSQ, 2) / std::cos( 3.0* LD );
         

      }
      else {
         this->CalculateSmoothingInvariants( StressInvariants);
         C2V = LodeDerivative;
         C2V *= -6.0 * StressInvariants.B * std::cos( 3.0 * StressInvariants.LodeAngle) - 15.0 * StressInvariants.C * std::sin( 6.0 * StressInvariants.LodeAngle) ;


         Vector Aux1, Aux2; 

         Aux1 = LodeDerivative; 
         Aux1 *= -3.0 * StressInvariants.C * std::cos( 3.0 * StressInvariants.LodeAngle) ; 

         Aux2 = C2Vector;
         Aux2 *= StressInvariants.B  +  2.0 * StressInvariants.C * std::sin( 3.0 * StressInvariants.LodeAngle) ; 
         Aux2 /= StressInvariants.J2InvSQ ;

         C3V = Aux1 + Aux2;
         C3V *=  + sqrt(3) * 3.0 / pow ( StressInvariants.J2InvSQ, 2);

         if ( fabs( StressInvariants.J2InvSQ) < 1e-6)
            C3V = ZeroVector(6);

      }

      rYieldFunctionDD = C2 * C2Matrix + C3*C3Matrix;

      Matrix Tonti = rYieldFunctionDD; 
      for (unsigned int i = 0; i < 6; ++i) {
         for (unsigned int j = 0; j < 6; j++) {
            rYieldFunctionDD(i,j) += (C2V(i)* C2Vector(j) ); // + C3V(i)*C3Vector(j) ) / 2.0;
            rYieldFunctionDD(i,j) += (C3V(i)* C3Vector(j) ); //+ C3V(i)*C3Vector(j) ) / 2.0;
         }
      }
     
      return;


      // ESCRIURE TOTA LA PUTA DERIVADA
      std::cout << " C2Matrix " << std::endl;
      for (unsigned int i = 0; i < 6; i++) {
         for (unsigned int j = 0; j < 6; j++) {
            std::cout << C2Matrix(i,j) * C2 << " ";
         }
         std::cout << std::endl; 
      }
      std::cout << " C3Matrix " << std::endl;
      for (unsigned int i = 0; i < 6; i++) {
         for (unsigned int j = 0; j < 6; j++) {
            std::cout << C3Matrix(i,j) * C3 << " ";
         }
         std::cout <<   std::endl; 
      }

      Matrix B1 = ZeroMatrix(6,6);
      Matrix B2 = B1; 

      for (unsigned int i = 0; i < 6; i++) {
         for (unsigned int j = 0 ; j < 6; ++j) {
            B1(i,j) = C2V(i) * C2Vector(j) ; //+ C2V(j) * C2Vector(i);
            B2(i,j) = C3V(i) * C3Vector(j) ; //+ C3V(j) * C3Vector(i);
         }
      }

      std::cout << " B1Matrix " << std::endl;
      for (unsigned int i = 0; i < 6; i++) {
         for (unsigned int j = 0; j < 6; j++) {
            std::cout << B1(i,j) << " ";
         }
         std::cout << std::endl; 
      }

      std::cout << " B2Matrix " << std::endl;
      for (unsigned int i = 0; i < 6; i++) {
         for (unsigned int j = 0; j < 6; j++) {
            std::cout << B2(i,j) << " ";
         }
         std::cout <<   std::endl; 
      }


   }*/


   void TrescaYieldCriterion::ComputeC2andC3Vector( const Vector&  rStressVector, TrescaStressInvariants& rStressInvariants, Vector& C2Vector, Vector & C3Vector)
   {

      Vector ShearStress = rStressVector;
      for (unsigned int i = 0; i < 3; ++i )
         ShearStress(i) -= rStressInvariants.MeanStress;

      C2Vector = ShearStress;
      for (unsigned int i = 3; i < 6; ++i )
         C2Vector(i) *= 2.0;

      C2Vector /= 2.0 * rStressInvariants.J2InvSQ;

      if ( fabs( rStressInvariants.J2InvSQ) < 1E-6)
         C2Vector = ZeroVector(6);

      C3Vector = ZeroVector(6);


      C3Vector(0) = ShearStress(1)*ShearStress(2) - pow( ShearStress(4), 2); 
      C3Vector(1) = ShearStress(2)*ShearStress(0) - pow( ShearStress(5), 2); 
      C3Vector(2) = ShearStress(0)*ShearStress(1) - pow( ShearStress(3), 2); 

      C3Vector(3) = 2.0 * ( ShearStress(4)*ShearStress(5) - ShearStress(2)*ShearStress(3));
      C3Vector(4) = 2.0 * ( ShearStress(5)*ShearStress(3) - ShearStress(0)*ShearStress(4));
      C3Vector(5) = 2.0 * ( ShearStress(3)*ShearStress(4) - ShearStress(1)*ShearStress(5));

      for (unsigned int i = 0; i < 3; ++i)
         C3Vector(i) += pow(rStressInvariants.J2InvSQ, 2) / 3.0;

      Matrix Aux, ShearStressM;
      ShearStressM = MathUtils<double>::StressVectorToTensor( ShearStress);
      Aux = prod(  ShearStressM, ShearStressM);

      for (unsigned int i = 0; i < 3; i++)
         Aux(i,i) -= 1.0/3.0 * 2.0*pow( rStressInvariants.J2InvSQ, 2);

      C3Vector = MathUtils<double>::StrainTensorToVector( Aux, 6);




   }

   void TrescaYieldCriterion::ComputeC2andC3VectorDD( const Vector&  rStressVector, TrescaStressInvariants& rStressInvariants, Vector& C2Vector, Vector & C3Vector, Matrix& C2Matrix, Matrix& C3Matrix, Vector & rLodeDerivative)
   {
      ComputeC2andC3Vector( rStressVector, rStressInvariants, C2Vector, C3Vector);

      C3Matrix = ZeroMatrix(6,6);
      rLodeDerivative = ZeroVector(6);

      Vector ShearStress = rStressVector; 
      for (unsigned int i = 0; i < 3; i++)
         ShearStress(i) -= rStressInvariants.MeanStress;


      C2Matrix = ZeroMatrix(6,6);
      double times;
      for (unsigned int i = 0; i < 6; i++) {
         for (unsigned int j = 0; j < 6; j++) {
            times = 0.5;
            if ( i > 2)
               times *= 2;
            if (j > 2)
               times *= 2;
            C2Matrix(i,j) = -  times * ShearStress(i) * ShearStress(j) / pow( rStressInvariants.J2InvSQ, 2) ;
         }
      }

      for (unsigned int i = 0; i < 3; i++)
      {
         C2Matrix(i,i) += 1.0;
         for (unsigned int j = 0; j < 3; j++)
            C2Matrix(i,j) -= 1.0/3.0;
      }

      for (unsigned int i = 3; i < 6; i++)
         C2Matrix(i,i) += 2.0;

      C2Matrix /= 2.0 * rStressInvariants.J2InvSQ;

      C3Matrix = ZeroMatrix(6,6);

      
      // FIRST QUAT
      C3Matrix(0,0) = ShearStress(0) - ShearStress(1) - ShearStress(2);
      C3Matrix(1,1) = ShearStress(1) - ShearStress(0) - ShearStress(2);
      C3Matrix(2,2) = ShearStress(2) - ShearStress(1) - ShearStress(0);

      C3Matrix(1,0) = 2 * ShearStress(2);
      C3Matrix(2,0) = 2 * ShearStress(1);
      C3Matrix(2,1) = 2 * ShearStress(0);

      //
      C3Matrix(3,0) = 2 * ShearStress(3);
      C3Matrix(3,1) = 2 * ShearStress(3);
      C3Matrix(3,2) = -4 * ShearStress(3);

      C3Matrix(4,0) = -4 * ShearStress(4);
      C3Matrix(4,1) = 2 * ShearStress(4);
      C3Matrix(4,2) = 2 * ShearStress(4);

      C3Matrix(5,0) = 2 * ShearStress(5);
      C3Matrix(5,1) = -4 * ShearStress(5);
      C3Matrix(5,2) = 2 * ShearStress(5);

      C3Matrix(3,3) = -6 * ShearStress(2);

      C3Matrix(4,3) = 6 * ShearStress(4);
      C3Matrix(4,4) = -6 * ShearStress(0);

      C3Matrix(5,3) = 6 * ShearStress(5);
      C3Matrix(5,4) = 6 * ShearStress(3);
      C3Matrix(5,5) = -6 * ShearStress(1);  // not shure but who cares

      for (unsigned int i = 0; i < 6; i++) {
         for (unsigned int j = i+1; j < 6; j++)
         {
            C3Matrix(i,j) = C3Matrix(j,i);
         }
      }
      C3Matrix /= 3.0;

      rLodeDerivative = ZeroVector(6);
      Matrix ShearStressM = MathUtils<double>::StressVectorToTensor( ShearStress);
      double J3 = MathUtils<double>::Det( ShearStressM);

      rLodeDerivative = C3Vector - 3.0*J3 / rStressInvariants.J2InvSQ * C2Vector; 
      

      rLodeDerivative *= - sqrt(3) / 2.0 / pow( rStressInvariants.J2InvSQ, 3) / std::cos( 3.0 * rStressInvariants.LodeAngle);  // un dos un tres no se ni idea



   }

   // *********** CALCULATE SMOOTHING (OF THE CORNER ) CONSTANTS  **********
   // **********************************************************************
   void TrescaYieldCriterion::CalculateSmoothingInvariants(  TrescaStressInvariants& rStressInvariants)
   {
      // NO SE PERQUE TINC LA B DE SIGNE CONTRARI, QUE COSAS...

      double SignedSmoothing = this->GetSmoothingLodeAngle();
      double SmoothingAngle = this->GetSmoothingLodeAngle();

      if ( rStressInvariants.LodeAngle < 0.0)
         SignedSmoothing = -SmoothingAngle;

      /*
      //rStressInvariants.A = std::cos( SignedSmoothing ) + std::sin( 3.0*SignedSmoothing ) * std::tan( 3.0*SignedSmoothing ) / 3.0;
      rStressInvariants.A = std::cos( SignedSmoothing ) * ( 3.0 + std::tan( SignedSmoothing ) * std::tan(3.0 * SignedSmoothing) ) / 3.0;
      rStressInvariants.A = std::cos( SmoothingAngle ) *  ( 3.0 + std::tan(SmoothingAngle) * std::tan(3.0*SmoothingAngle) ) / 3.0;

      rStressInvariants.B = std::sin( SignedSmoothing ) / ( 3.0 * std::cos( SignedSmoothing * 3.0 ) );
      rStressInvariants.B = std::sin( SignedSmoothing ) / ( 3.0 * std::cos( 3.0*SmoothingAngle) );
      */
      // VALE ARA EL C2

      double Denom = 18.0 * pow ( std::cos( 3.0 * SmoothingAngle) , 3 );
      
      rStressInvariants.C = - std::cos( 3.0*SmoothingAngle) * std::cos( SmoothingAngle) - 3.0*std::sin ( 3.0 * SmoothingAngle) * std::sin( SmoothingAngle);

      rStressInvariants.B = std::sin( 6.0 * SignedSmoothing) * std::cos( SmoothingAngle) - 6.0 * std::cos( 6.0 * SmoothingAngle) * std::sin( SignedSmoothing) ;

      rStressInvariants.B /= Denom;
      rStressInvariants.C /= Denom;

      rStressInvariants.A = - rStressInvariants.B * std::sin( 3.0 * SignedSmoothing) - rStressInvariants.C * pow( std::sin( 3.0 * SmoothingAngle), 2 ) + std::cos( SmoothingAngle) ;

   }

   // ***************** COMPUTE STRESS INVARIANTS ******************************
   // **************************************************************************
   void TrescaYieldCriterion::CalculateStressInvariants( const Vector& rStressVector, TrescaStressInvariants& rStressInvariants)
   {

      Matrix StressMatrix;
      StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector);

      rStressInvariants.MeanStress = 0.0;
      for (unsigned int i = 0; i < 3; ++i)
         rStressInvariants.MeanStress += StressMatrix(i,i)/3.0;


      rStressInvariants.J2InvSQ = 0.0;
      for (unsigned int i = 0; i <3; ++i) 
         rStressInvariants.J2InvSQ += pow( StressMatrix(i,i) - rStressInvariants.MeanStress, 2);

      rStressInvariants.J2InvSQ += 2.0*pow( StressMatrix(0,1) , 2);
      rStressInvariants.J2InvSQ += 2.0*pow( StressMatrix(0,2) , 2);
      rStressInvariants.J2InvSQ += 2.0*pow( StressMatrix(1,2) , 2);

      rStressInvariants.J2InvSQ = sqrt(  rStressInvariants.J2InvSQ / 2.0);


      // THIRD INVARIANT COMPUTATION
      rStressInvariants.LodeAngle = 0.0;

      for (unsigned int i = 0; i < 3 ; ++i ) 
         StressMatrix(i,i) -= rStressInvariants.MeanStress;

      rStressInvariants.LodeAngle  = MathUtils<double>::Det( StressMatrix);
      rStressInvariants.LodeAngle *= 3.0*sqrt(3.0) / 2.0;    
      rStressInvariants.LodeAngle /= pow( rStressInvariants.J2InvSQ, 3);

      // SHOULD BE LESS THAN ONE TO COMPUTE INVERSE SINUS
      double Epsi = 10e-5;
      if ( fabs(rStressInvariants.LodeAngle) > 1.0-Epsi ) { 
         rStressInvariants.LodeAngle = -30.0*GetPI() / 180.0 * rStressInvariants.LodeAngle / fabs(rStressInvariants.LodeAngle) ;
      }
      else {
         rStressInvariants.LodeAngle  = std::asin( -rStressInvariants.LodeAngle) / 3.0;
      }
      if (fabs(rStressInvariants.J2InvSQ) < 0.001)  
         rStressInvariants.LodeAngle = 0.0;
      // RADIANTS 

   }





   // DEFINE SOME CONSTANTS.
   double TrescaYieldCriterion::GetSmoothingLodeAngle()
   {
      //return 29.8*GetPI()/180.0;
      //return 22.8*GetPI()/180.0;
      return 29.0*GetPI()/180.0;
   }

   double TrescaYieldCriterion::GetPI()
   {
      return 3.14159265359;
   }


   void TrescaYieldCriterion::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YieldCriterion )
   }

   void TrescaYieldCriterion::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion )
   }


}  // namespace Kratos.
