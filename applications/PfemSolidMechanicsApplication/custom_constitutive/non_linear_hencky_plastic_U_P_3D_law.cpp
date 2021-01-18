//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/non_linear_hencky_plastic_U_P_3D_law.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{



// CONSTRUCTORS
NonLinearHenckyElasticPlasticUP3DLaw::NonLinearHenckyElasticPlasticUP3DLaw()
   : NonLinearHenckyElasticPlastic3DLaw()
{

}

NonLinearHenckyElasticPlasticUP3DLaw::NonLinearHenckyElasticPlasticUP3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
   : NonLinearHenckyElasticPlastic3DLaw()
{

}

// ************* COPY CONSTRUCTOR ******************
NonLinearHenckyElasticPlasticUP3DLaw::NonLinearHenckyElasticPlasticUP3DLaw(const NonLinearHenckyElasticPlasticUP3DLaw&  rOther)
  : NonLinearHenckyElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NonLinearHenckyElasticPlasticUP3DLaw::Clone() const
{
    NonLinearHenckyElasticPlasticUP3DLaw::Pointer p_clone(new NonLinearHenckyElasticPlasticUP3DLaw(*this));
    return p_clone;
}

NonLinearHenckyElasticPlasticUP3DLaw::~NonLinearHenckyElasticPlasticUP3DLaw()
{
}




void NonLinearHenckyElasticPlasticUP3DLaw::CorrectDomainPressure(Matrix& rStressMatrix, const MaterialResponseVariables& rElasticVariables)
{

    double MeanPressure = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
        MeanPressure += rStressMatrix(i,i);

    MeanPressure /=3.0;
    //if ( fabs(MeanPressure) > 1.0E-4)
    //   std::cout << " UNCORRECTED PRESSURE " << MeanPressure << std::endl;

    for (unsigned int i = 0; i < 3; ++i)
        rStressMatrix(i,i) -= MeanPressure;


    double Pressure = 0;
    GetDomainPressure( Pressure, rElasticVariables);

    for (unsigned int i = 0; i < 3; ++i)
        rStressMatrix(i,i) += Pressure * rElasticVariables.DeterminantF;

    //std::cout << " THIS DET " << rElasticVariables.DeterminantF << std::endl;
}

void NonLinearHenckyElasticPlasticUP3DLaw::GetDomainPressure( double& rPressure, const MaterialResponseVariables& rElasticVariables)
{

    rPressure = 0.0;
    const GeometryType&  DomainGeometry =  rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues  =  rElasticVariables.GetShapeFunctionsValues();

    const unsigned int number_of_nodes  =  DomainGeometry.size();

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        rPressure += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(PRESSURE); //NOOOOO
    }

}



void NonLinearHenckyElasticPlasticUP3DLaw::CalculateElastoPlasticTangentMatrix( const FlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen,const double& rAlpha, Matrix& rElastoPlasticTangentMatrix, const MaterialResponseVariables& rElasticVariables)
{

     mpFlowRule->ComputeElastoPlasticTangentMatrix( rReturnMappingVariables,  rNewElasticLeftCauchyGreen, rAlpha, rElastoPlasticTangentMatrix);


     Matrix StressTensor = mpFlowRule->ComputeKirchhoffStressMatrix( rNewElasticLeftCauchyGreen);
     CorrectDomainPressure( StressTensor, rElasticVariables);
     Matrix ExtraMatrix = this->CalculateExtraMatrix( StressTensor);

     rElastoPlasticTangentMatrix += ExtraMatrix;


     // ADDING THE K TERMS


     double Young = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
     double Nu = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];

     double K = Young / 3.0 / (1.0 - 2.0*Nu);

     for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3 ; ++j) {
           rElastoPlasticTangentMatrix(i,j)  -= K;
         }
     }

     double Pressure;
     //GetDomainPressure( Pressure, rElasticVariables);
     GetDomainPressure( Pressure, rElasticVariables);

     Pressure *= rElasticVariables.DeterminantF;


     Matrix DeviatoricTensor = ZeroMatrix(6,6);
     for (unsigned int i = 0; i < 6; i++)
        DeviatoricTensor(i,i) = 1.0;

     for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
           DeviatoricTensor(i,j) = -1.0/3.0;
        }
     }

     //rElastoPlasticTangentMatrix = prod( DeviatoricTensor, rElastoPlasticTangentMatrix);


     Matrix FourthOrderIdentity = ZeroMatrix(6,6);
     for (unsigned int i = 0; i<3; ++i)
        FourthOrderIdentity(i,i) = 1.0;

     for (unsigned int i = 3; i<6; ++i)
        FourthOrderIdentity(i,i) = 0.50;
        // VOIGT NOTATION AND NOT KELVIN

     Matrix IdentityCross = ZeroMatrix(6,6);
     for (unsigned int i = 0; i<3; ++i) {
          for (unsigned int j = 0; j<3; ++j) {
             IdentityCross(i,j) = 1.0;
          }
     }

     rElastoPlasticTangentMatrix += Pressure* ( IdentityCross - 2.0 * FourthOrderIdentity);

     double det =  MathUtils<double>::Det( mElasticLeftCauchyGreen);
     mElasticLeftCauchyGreen /= pow( det, 1/3);
}

// I HAVE TO DELETE THE DEFINITION FROM THE HPP IN ORDER TO REMOVE THE TO FOLLOWING FUNCTIONS; THAT ARE COPY PASTE; THEY ARE ALREADY DEFINED BY INHERITANCE !!!!
void NonLinearHenckyElasticPlasticUP3DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = 0.5 * ( rRightCauchyGreen( 2, 2 ) - 1.00 );
    rStrainVector[3] = rRightCauchyGreen( 0, 1 ); // xy
    rStrainVector[4] = rRightCauchyGreen( 1, 2 ); // yz
    rStrainVector[5] = rRightCauchyGreen( 0, 2 ); // xz

}



//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************
void NonLinearHenckyElasticPlasticUP3DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{

    // e = 0.5*(1-invFT*invF) or e = 0.5*(1-inv(b))

    //Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen ( 3, 3 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 2, 2 ) );
    rStrainVector[3] = - InverseLeftCauchyGreen( 0, 1 ); // xy
    rStrainVector[4] = - InverseLeftCauchyGreen( 1, 2 ); // yz
    rStrainVector[5] = - InverseLeftCauchyGreen( 0, 2 ); // xz
}




void NonLinearHenckyElasticPlasticUP3DLaw::CalculateOnlyDeviatoricPart( Matrix& rIncrementalDeformationGradient)
{
     Matrix Aux = rIncrementalDeformationGradient;

     double det = 0.0;
     det  = Aux(0,0) * Aux(1,1) * Aux(2,2);
     det += Aux(1,0) * Aux(2,1) * Aux(0,2);
     det += Aux(2,0) * Aux(0,1) * Aux(1,2);
     det -= Aux(0,2) * Aux(1,1) * Aux(2,0);
     det -= Aux(1,2) * Aux(2,1) * Aux(0,0);
     det -= Aux(2,2) * Aux(0,1) * Aux(1,0);

     det = pow( det, 1.0/3.0);

     rIncrementalDeformationGradient /= det;

}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void NonLinearHenckyElasticPlasticUP3DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( FINITE_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );
	rFeatures.mOptions.Set( U_P_LAW );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}

Matrix& NonLinearHenckyElasticPlasticUP3DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
{
   if ( rThisVariable == KIRCHHOFF_STRESS_TENSOR )
   {
      Matrix StressMatrix;
      StressMatrix = mpFlowRule->ComputeKirchhoffStressMatrix( mElasticLeftCauchyGreen);

      // Check if p == 0 and correct
      double MeanStress = 0;
      for (unsigned int i = 0; i < 3; ++i)
         MeanStress += StressMatrix(i,i);
      MeanStress /= 3.0;
      for (unsigned int i = 0; i < 3; ++i)
         StressMatrix(i,i) -= MeanStress;


      rValue = StressMatrix;

   }
   else {
      rValue = NonLinearHenckyElasticPlastic3DLaw::GetValue( rThisVariable, rValue);
   }
   return rValue;
}

void NonLinearHenckyElasticPlasticUP3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
{

   if ( rThisVariable == ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS)
   {
      Vector ThisVector = rValue;
      double MeanStress = 0;
      for (unsigned int i = 0; i < 3; i++)
         MeanStress += ThisVector(i);
      MeanStress /= 3.0;
      for (unsigned int i = 0; i < 3; i++)
         ThisVector(i) -= MeanStress;

      NonLinearHenckyElasticPlastic3DLaw::SetValue( rThisVariable, ThisVector, rCurrentProcessInfo);

   }
   else {
      NonLinearHenckyElasticPlastic3DLaw::SetValue( rThisVariable, rValue, rCurrentProcessInfo);
   }

}



} //end namespace kratos
