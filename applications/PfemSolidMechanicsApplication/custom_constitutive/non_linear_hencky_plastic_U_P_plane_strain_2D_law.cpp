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
#include "custom_constitutive/non_linear_hencky_plastic_U_P_plane_strain_2D_law.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{



// CONSTRUCTORS
NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw::NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw()
   : NonLinearHenckyElasticPlasticUP3DLaw()
{

}

NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw::NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
   : NonLinearHenckyElasticPlasticUP3DLaw()
{

}

// ************* COPY CONSTRUCTOR ******************
NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw::NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw(const NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw&  rOther)
  : NonLinearHenckyElasticPlasticUP3DLaw(rOther)
{

}

ConstitutiveLaw::Pointer NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw::Clone() const
{
    NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw::Pointer p_clone(new NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw(*this));
    return p_clone;
}


NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw::~NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw()
{
}






// I HAVE TO DELETE THE DEFINITION FROM THE HPP IN ORDER TO REMOVE THE TO FOLLOWING FUNCTIONS; THAT ARE COPY PASTE; THEY ARE ALREADY DEFINED BY INHERITANCE !!!!
void NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = rRightCauchyGreen( 0, 1 );

}



//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************
void NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{

    // e = 0.5*(1-invFT*invF) or e = 0.5*(1-inv(b))

    //Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen ( rLeftCauchyGreen.size1() , rLeftCauchyGreen.size2() );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector.clear();
    rStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = -InverseLeftCauchyGreen( 0, 1 );
}


Matrix NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw::SetConstitutiveMatrixToAppropiateDimension(const Matrix& rElastoPlasticTangentMatrix)
{
     Matrix Result = ZeroMatrix(3,3);

     Result(0, 0) = rElastoPlasticTangentMatrix(0, 0);
     Result(0, 1) = rElastoPlasticTangentMatrix(0, 1);
     Result(1, 0) = rElastoPlasticTangentMatrix(1, 0);
     Result(1, 1) = rElastoPlasticTangentMatrix(1, 1);


     Result(2, 0) = rElastoPlasticTangentMatrix(3, 0);
     Result(2, 1) = rElastoPlasticTangentMatrix(3, 1);
     Result(2, 2) = rElastoPlasticTangentMatrix(3, 3);


     Result(0, 2) = rElastoPlasticTangentMatrix(0, 3);
     Result(1, 2) = rElastoPlasticTangentMatrix(1, 3);


     return Result;

}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************


void NonLinearHenckyElasticPlasticUPPlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
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
