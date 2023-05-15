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
#include "custom_constitutive/non_linear_hencky_plastic_U_P_axisym_2D_law.hpp"

#include "pfem_solid_mechanics_application_variables.h"

//Molt important, el tema de constructors... etc
namespace Kratos
{


NonLinearHenckyElasticPlasticUPAxisym2DLaw::NonLinearHenckyElasticPlasticUPAxisym2DLaw()
   : NonLinearHenckyElasticPlasticUP3DLaw()
{

}

NonLinearHenckyElasticPlasticUPAxisym2DLaw::NonLinearHenckyElasticPlasticUPAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
   : NonLinearHenckyElasticPlasticUP3DLaw()
{

}



// ************* COPY CONSTRUCTOR ******************
NonLinearHenckyElasticPlasticUPAxisym2DLaw::NonLinearHenckyElasticPlasticUPAxisym2DLaw(const NonLinearHenckyElasticPlasticUPAxisym2DLaw&  rOther)
  : NonLinearHenckyElasticPlasticUP3DLaw(rOther)
{

}


//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NonLinearHenckyElasticPlasticUPAxisym2DLaw::Clone() const
{
    NonLinearHenckyElasticPlasticUPAxisym2DLaw::Pointer p_clone(new NonLinearHenckyElasticPlasticUPAxisym2DLaw(*this));
    return p_clone;
}

NonLinearHenckyElasticPlasticUPAxisym2DLaw::~NonLinearHenckyElasticPlasticUPAxisym2DLaw()
{
}

//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void NonLinearHenckyElasticPlasticUPAxisym2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = 0.5 * ( rRightCauchyGreen( 2, 2 ) - 1.00 );
    rStrainVector[3] = rRightCauchyGreen( 0, 1 );

}


//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void NonLinearHenckyElasticPlasticUPAxisym2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{

    // e= 0.5*(1-invbT*invb)
    Matrix InverseLeftCauchyGreen ( rLeftCauchyGreen.size1() , rLeftCauchyGreen.size2() );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector.clear();
    rStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 2, 2 ) );
    rStrainVector[3] = -InverseLeftCauchyGreen( 0, 1 );


}


Matrix NonLinearHenckyElasticPlasticUPAxisym2DLaw::SetConstitutiveMatrixToAppropiateDimension(const Matrix& rElastoPlasticTangentMatrix)
{
     Matrix Result = ZeroMatrix(4,4);


     for (unsigned int i = 0; i < 4; ++i)  {
         for (unsigned int j = 0; j < 4; ++j) {
             Result( i, j) = 1.0*rElastoPlasticTangentMatrix( i, j);
         }
     }

     double det =  MathUtils<double>::Det( mElasticLeftCauchyGreen);
     mElasticLeftCauchyGreen /= pow( det, 1/3);
     return Result;

}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void NonLinearHenckyElasticPlasticUPAxisym2DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( AXISYMMETRIC_LAW );
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



} //end namespace kratos
