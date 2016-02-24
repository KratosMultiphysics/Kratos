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
#include "custom_constitutive/custom_flow_rules/mohr_coulomb_explicit_plastic_flow_rule.hpp"

#include "pfem_solid_mechanics_application_variables.h"


namespace Kratos
{



//************ CONSTRUCTOR ***********
MohrCoulombExplicitFlowRule::MohrCoulombExplicitFlowRule()
   :TrescaExplicitFlowRule()
{
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

MohrCoulombExplicitFlowRule::MohrCoulombExplicitFlowRule(YieldCriterionPointer pYieldCriterion)
	:TrescaExplicitFlowRule(pYieldCriterion)
{
   
}

//********* ASSIGMENT OPERATOR
MohrCoulombExplicitFlowRule& MohrCoulombExplicitFlowRule::operator=(MohrCoulombExplicitFlowRule const& rOther)
{
	TrescaExplicitFlowRule::operator=(rOther);
	return *this;

}



//********** COPY CONSTRUCTOR *********
MohrCoulombExplicitFlowRule::MohrCoulombExplicitFlowRule(MohrCoulombExplicitFlowRule const& rOther)
      :TrescaExplicitFlowRule(rOther)
{
}

//*******   CLONE ********
FlowRule::Pointer MohrCoulombExplicitFlowRule::Clone() const
{
  FlowRule::Pointer p_clone(new MohrCoulombExplicitFlowRule(*this));
  return p_clone;
}



// ********** DESTRUCTOR **************
MohrCoulombExplicitFlowRule::~MohrCoulombExplicitFlowRule()
{
}

// THE IMPORTANT PART, THE DERIVATIVES
void MohrCoulombExplicitFlowRule::CalculatePlasticPotentialDerivatives(const Vector& rStressVector, Vector& rFirstDerivative, Matrix& rSecondDerivative)
{
     //double YieldStress = mpYieldCriterion->GetHardeningLaw().GetProperties()[YIELD_STRESS];
 
   rFirstDerivative = ZeroVector(1);
   rSecondDerivative = ZeroMatrix(1);

   //double Cohesion = this->GetHardeningLaw().GetProperties()[COHESION];
   double FrictionAngle = mpYieldCriterion->GetHardeningLaw().GetProperties()[INTERNAL_DILATANCY_ANGLE];
   FrictionAngle *= GetPI() / 180.0;
 
   MCStressInvariants StressInvariants;
   MCSmoothingConstants SmoothingConstants;

   this->CalculateStressInvariants(rStressVector, StressInvariants );  // ME SIRVEN LOS DE ANTES ( de aquella manera)


   double LodeCut = this->GetSmoothingLodeAngle();
   double HiperbolicA = this->GetSmoothingHiperbolic();
 
   double K, Kderivative;
   double Alpha;
   double C1, C2, C3;


   C1 = std::sin( FrictionAngle);
   if ( fabs(StressInvariants.LodeAngle) < LodeCut ) {

        K = std::cos( StressInvariants.LodeAngle) - 1.0/sqrt(3.0) * std::sin(FrictionAngle) * std::sin( StressInvariants.LodeAngle );
        Kderivative = -std::sin( StressInvariants.LodeAngle) - 1.0/sqrt(3.0) * std::sin(FrictionAngle)*std::cos( StressInvariants.LodeAngle); 
  
        C2 = K - std::tan( 3.0*StressInvariants.LodeAngle ) * Kderivative;
        C3 = - sqrt(3.0) / ( 2.0 * pow( StressInvariants.J2InvSQ, 2.0) * std::cos( 3.0*StressInvariants.LodeAngle) ) * Kderivative;

   }
   else {

        this->CalculateSmoothingConstants( SmoothingConstants, StressInvariants);

        K = SmoothingConstants.A + SmoothingConstants.B * std::sin( 3.0*StressInvariants.LodeAngle) ;
        Kderivative = 3.0*SmoothingConstants.B * std::cos( 3.0 * StressInvariants.LodeAngle);


        C2 = K - 3.0 * SmoothingConstants.B * std::sin(3.0*StressInvariants.LodeAngle);
        C3 =   - 3.0 * sqrt(3.0) * SmoothingConstants.B / ( 2.0* pow(StressInvariants.J2InvSQ, 2.0) ); // to avoid dividing by cos(3teta)
   }

   Alpha = pow ( StressInvariants.J2InvSQ * K, 2.0) + pow( HiperbolicA * std::sin( FrictionAngle), 2.0);
   Alpha = K*StressInvariants.J2InvSQ / sqrt( Alpha);

   C2 *= Alpha;
   C3 *= Alpha;

     Vector C1Vector = ZeroVector(6);
     for (unsigned int i = 0; i < 3; ++i)
          C1Vector(i) = 1.0/3.0;

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

     rFirstDerivative = C1 * C1Vector  + C2*C2Vector + C3*C3Vector;

/*     std::cout << " AND FINALLY " << rYieldFunctionD << std::endl;
     
     Matrix StressMatrix = MathUtils<double>::StressVectorToTensor(rStressVector); 
     Matrix EigenVectors;
     Vector EigenValues;
    
     SolidMechanicsMathUtilities<double>::EigenVectors( StressMatrix, EigenVectors, EigenValues);

     std::cout << "     StressVECTOR " << rStressVector << std::endl;
     std::cout << " EIGENV " <<EigenValues << " VECTORS " << EigenVectors << std::endl;

     std::cout << " FIRST C1 " << C1 << " VECTOR " << C1Vector << std::endl;
     std::cout << C1*C1Vector << std::endl;
     std::cout << " FIRST C2 " << C2 << " VECTOR " << C2Vector << std::endl;
     std::cout << C2*C2Vector << std::endl;
     std::cout << " FIRST C3 " << C3 << " VECTOR " << C3Vector << std::endl;
     std::cout << C3*C3Vector << std::endl;
     std::cout << " alpHA " << Alpha << std::endl;
     std::cout << std::endl;
*/

}

void MohrCoulombExplicitFlowRule::CalculateSmoothingConstants( MCSmoothingConstants& rSmoothingConstants, const MCStressInvariants& rStressInvariants)
{
    double SmoothingAngle = this->GetSmoothingLodeAngle();
    double FrictionAngle = mpYieldCriterion->GetHardeningLaw().GetProperties()[INTERNAL_DILATANCY_ANGLE];
    FrictionAngle *= GetPI() / 180.0;

    double Sign = 1.0;
    if ( rStressInvariants.LodeAngle < 0.0)
           Sign = -1.0;

    rSmoothingConstants.A = 3.0 +  std::tan(SmoothingAngle) * std::tan(3.0*SmoothingAngle) + Sign * (std::tan( 3.0*SmoothingAngle) - 3.0*std::tan(SmoothingAngle)) * std::sin( FrictionAngle) / sqrt(3.0);
    rSmoothingConstants.A *= (1.0/3.0) * std::cos( SmoothingAngle );

    rSmoothingConstants.B = -1.0/ ( 3.0*std::cos( SmoothingAngle) ) * ( Sign * std::sin(SmoothingAngle) + std::sin(FrictionAngle)*std::cos(SmoothingAngle) / sqrt(3.0));
    rSmoothingConstants.B = -1.0 * ( Sign* std::sin(SmoothingAngle) + std::sin(FrictionAngle)*std::cos(SmoothingAngle) / sqrt(3.0) ) / ( 3.0*std::cos(3.0*SmoothingAngle) );

}


// CALCULATE THE INVARIANTS. (NOTHING NEEDS TO BE CHANGED, the are invariants)
void MohrCoulombExplicitFlowRule::CalculateStressInvariants( const Vector& rStressVector, MCStressInvariants& rStressInvariants)
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


double MohrCoulombExplicitFlowRule::GetSmoothingLodeAngle()
{
    return 29.0*GetPI()/180.0;
    //return 22.8*GetPI()/180.0;
    //return 10.0*GetPI()/180.0;
}


double MohrCoulombExplicitFlowRule::GetPI()
{
   return 3.14159265359;
}

double MohrCoulombExplicitFlowRule::GetSmoothingHiperbolic()
{
   return 10.0;
}


void MohrCoulombExplicitFlowRule::ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH)
{
   
     rH = 0.0;

}
 
void MohrCoulombExplicitFlowRule::save( Serializer& rSerializer) const 
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
}

void MohrCoulombExplicitFlowRule::load( Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )

}

} //end namespace kratos
