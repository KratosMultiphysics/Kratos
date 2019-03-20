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
#include "custom_constitutive/non_linear_hencky_plastic_axisym_2D_law.hpp"

#include "pfem_solid_mechanics_application_variables.h"

//Molt important, el tema de constructors... etc
namespace Kratos
{


NonLinearHenckyElasticPlasticAxisym2DLaw::NonLinearHenckyElasticPlasticAxisym2DLaw()
   : NonLinearHenckyElasticPlastic3DLaw()
{

}

NonLinearHenckyElasticPlasticAxisym2DLaw::NonLinearHenckyElasticPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
   : NonLinearHenckyElasticPlastic3DLaw()
{

}



// ************* COPY CONSTRUCTOR ******************
NonLinearHenckyElasticPlasticAxisym2DLaw::NonLinearHenckyElasticPlasticAxisym2DLaw(const NonLinearHenckyElasticPlasticAxisym2DLaw&  rOther)
  : NonLinearHenckyElasticPlastic3DLaw(rOther)
{

}

ConstitutiveLaw::Pointer NonLinearHenckyElasticPlasticAxisym2DLaw::Clone() const
{
    NonLinearHenckyElasticPlasticAxisym2DLaw::Pointer p_clone(new NonLinearHenckyElasticPlasticAxisym2DLaw(*this));
    return p_clone;
}

NonLinearHenckyElasticPlasticAxisym2DLaw::~NonLinearHenckyElasticPlasticAxisym2DLaw()
{
}

//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void NonLinearHenckyElasticPlasticAxisym2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
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

void NonLinearHenckyElasticPlasticAxisym2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
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


Matrix NonLinearHenckyElasticPlasticAxisym2DLaw::SetConstitutiveMatrixToAppropiateDimension(const Matrix& rElastoPlasticTangentMatrix)
{
     Matrix Result = ZeroMatrix(4,4);


     for (unsigned int i = 0; i < 4; ++i)  {
         for (unsigned int j = 0; j < 4; ++j) {
             Result( i, j) = rElastoPlasticTangentMatrix( i, j);
         }
     }

     return Result;

}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void NonLinearHenckyElasticPlasticAxisym2DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( AXISYMMETRIC_LAW );
	rFeatures.mOptions.Set( FINITE_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}

} //end namespace kratos
