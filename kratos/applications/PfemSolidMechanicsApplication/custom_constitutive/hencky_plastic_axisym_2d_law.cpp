//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_plastic_axisym_2d_law.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticAxisym2DLaw::HenckyElasticPlasticAxisym2DLaw()
    : HenckyElasticPlastic3DLaw()
{

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticAxisym2DLaw::HenckyElasticPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : HenckyElasticPlastic3DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyElasticPlasticAxisym2DLaw::HenckyElasticPlasticAxisym2DLaw(const HenckyElasticPlasticAxisym2DLaw& rOther)
    : HenckyElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyElasticPlasticAxisym2DLaw::Clone() const
{
    HenckyElasticPlasticAxisym2DLaw::Pointer p_clone(new HenckyElasticPlasticAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticAxisym2DLaw::~HenckyElasticPlasticAxisym2DLaw()
{
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HenckyElasticPlasticAxisym2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
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

void HenckyElasticPlasticAxisym2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
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




void HenckyElasticPlasticAxisym2DLaw::SetConstitutiveMatrixToAppropiateDimension(Matrix&  rPrincipalTangentMatrix)
{

   Matrix Auxiliar = ZeroMatrix(4);
 
   for (unsigned int i = 0; i < 4; ++i) {
        for (unsigned int j = 0; j < 4; ++j) {
            Auxiliar(i,j) = rPrincipalTangentMatrix(i,j);
        }
    }

    rPrincipalTangentMatrix = Auxiliar;


}







} // Namespace Kratos

