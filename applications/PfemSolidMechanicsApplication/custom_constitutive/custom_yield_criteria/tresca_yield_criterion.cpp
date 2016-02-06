// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_constitutive/custom_yield_criteria/tresca_yield_criterion.hpp"

#include "pfem_solid_mechanics_application_variables.h"

// ROUNDED TRESCA YIELD CRITERION (Sloan & Booker, 1986 )

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

        rStateFunction = StressInvariants.J2InvSQ * ( StressInvariants.A - StressInvariants.B * std::sin( 3.0* StressInvariants.LodeAngle ) );
   }

   //rStateFunction = rStateFunction/YieldStress - 1.0;
   rStateFunction = rStateFunction - 1.0*YieldStress;
   return rStateFunction; 

}


void TrescaYieldCriterion::CalculateSmoothingInvariants(  TrescaStressInvariants& rStressInvariants)
{
   // NO SE PERQUE TINC LA B DE SIGNE CONTRARI, QUE COSAS...

    double SignedSmoothing = this->GetSmoothingLodeAngle();
    double SmoothingAngle = this->GetSmoothingLodeAngle();

    if ( rStressInvariants.LodeAngle < 0.0)
           SignedSmoothing = -SmoothingAngle;

    //rStressInvariants.A = std::cos( SignedSmoothing ) + std::sin( 3.0*SignedSmoothing ) * std::tan( 3.0*SignedSmoothing ) / 3.0;
    rStressInvariants.A = std::cos( SignedSmoothing ) * ( 3.0 + std::tan( SignedSmoothing ) * std::tan(3.0 * SignedSmoothing) ) / 3.0;
    rStressInvariants.A = std::cos( SmoothingAngle ) *  ( 3.0 + std::tan(SmoothingAngle) * std::tan(3.0*SmoothingAngle) ) / 3.0;

    rStressInvariants.B = std::sin( SignedSmoothing ) / ( 3.0 * std::cos( SignedSmoothing * 3.0 ) );
    rStressInvariants.B = std::sin( SignedSmoothing ) / ( 3.0 * std::cos( 3.0*SmoothingAngle) );

}

void TrescaYieldCriterion::CalculateStressInvariants( const Vector& rStressVector, TrescaStressInvariants& rStressInvariants)
{

         Matrix StressMatrix;
         StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector);

         rStressInvariants.MeanStress = 0.0;
         for (unsigned int i = 0; i < 3; ++i)
              rStressInvariants.MeanStress += StressMatrix(i,i)/3.0;


         rStressInvariants.J2InvSQ = 0.0;
         for (unsigned int i = 0; i <3; ++i) 
              rStressInvariants.J2InvSQ += pow( StressMatrix(i,i) - rStressInvariants.MeanStress, 2.0);

         rStressInvariants.J2InvSQ += 2.0*pow( StressMatrix(0,1) , 2.0);
         rStressInvariants.J2InvSQ += 2.0*pow( StressMatrix(0,2) , 2.0);
         rStressInvariants.J2InvSQ += 2.0*pow( StressMatrix(1,2) , 2.0);
         
         rStressInvariants.J2InvSQ = pow( rStressInvariants.J2InvSQ / 2.0, 1.0/2.0);


         // THIRD INVARIANT COMPUTATION
         rStressInvariants.LodeAngle = 0.0;

         for (unsigned int i = 0; i < 3 ; ++i ) 
              StressMatrix(i,i) -= rStressInvariants.MeanStress;

          rStressInvariants.LodeAngle  = MathUtils<double>::Det( StressMatrix);
          rStressInvariants.LodeAngle *= 3.0*pow( 3.0, 1.0/2.0) / 2.0;    
          rStressInvariants.LodeAngle /= pow( rStressInvariants.J2InvSQ, 3.0);

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



// SOLO ME QUEDA ESTO, y....



double TrescaYieldCriterion::GetSmoothingLodeAngle()
{
    return 29.8*GetPI()/180.0;
    //return 22.8*GetPI()/180.0;
    //return 10.0*GetPI()/180.0;
}


double TrescaYieldCriterion::GetPI()
{
   return 3.14159265359;
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

        C3 = -pow( 3.0, 1.0/2.0)/2.0 * std::sin( StressInvariants.LodeAngle) / std::cos( 3.0* StressInvariants.LodeAngle) ;
        C3 /= pow(StressInvariants.J2InvSQ, 2.0);
     }
     else {

        this->CalculateSmoothingInvariants( StressInvariants);
    
        C2 = StressInvariants.A + 2.0 * StressInvariants.B * std::sin ( 3.0 * StressInvariants.LodeAngle) ;

        C3 = 3.0 / 2.0 * pow ( 3.0, 1.0/2.0) * StressInvariants.B / pow(StressInvariants.J2InvSQ, 2.0);
     } 

     Vector ShearStress = rStressVector;
     for (unsigned int i = 0; i < 3; ++i )
         ShearStress(i) -= StressInvariants.MeanStress;

     Vector C2Vector = ShearStress;
     for (unsigned int i = 3; i < 6; ++i )
          C2Vector(i) *= 2.0;

     C2Vector /= 2.0 * StressInvariants.J2InvSQ;

     Vector C3Vector = ZeroVector(6);


     // FALTER TERMES
     C3Vector(0) = ShearStress(1)*ShearStress(2) - pow( ShearStress(4), 2.0); 
     C3Vector(1) = ShearStress(2)*ShearStress(0) - pow( ShearStress(5), 2.0); 
     C3Vector(2) = ShearStress(0)*ShearStress(1) - pow( ShearStress(3), 2.0); 

     C3Vector(3) = 2.0 * ( ShearStress(4)*ShearStress(5) - ShearStress(2)*ShearStress(3));
     C3Vector(4) = 2.0 * ( ShearStress(5)*ShearStress(3) - ShearStress(0)*ShearStress(4));
     C3Vector(5) = 2.0 * ( ShearStress(3)*ShearStress(4) - ShearStress(1)*ShearStress(5));

     for (unsigned int i = 0; i < 3; ++i)
         C3Vector(i) += pow(StressInvariants.J2InvSQ, 2.0) / 3.0;

     rYieldFunctionD = C2*C2Vector + C3*C3Vector;
     //rYieldFunctionD /= YieldStress;
 
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
