// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Malik Ali Dawi
//                   Ruben Zorrilla
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_plane_stress_2d.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/
HyperElasticIsotropicKirchhoffPlaneStress2D::HyperElasticIsotropicKirchhoffPlaneStress2D()
    : BaseType() 
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

HyperElasticIsotropicKirchhoffPlaneStress2D::HyperElasticIsotropicKirchhoffPlaneStress2D(const HyperElasticIsotropicKirchhoffPlaneStress2D& rOther)
    : BaseType(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer HyperElasticIsotropicKirchhoffPlaneStress2D::Clone() const 
{
    return Kratos::make_shared<HyperElasticIsotropicKirchhoffPlaneStress2D>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

HyperElasticIsotropicKirchhoffPlaneStress2D::~HyperElasticIsotropicKirchhoffPlaneStress2D() 
{
};

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void HyperElasticIsotropicKirchhoffPlaneStress2D::GetLawFeatures(Features& rFeatures) 
{

    //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRESS_LAW );
    rFeatures.mOptions.Set( FINITE_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_GreenLagrange);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicKirchhoffPlaneStress2D::CalculateConstitutiveMatrixPK2(
    Matrix& rConstitutiveMatrix,
    const double YoungModulus,
    const double PoissonCoefficient
    ) 
{
    rConstitutiveMatrix.clear();
    rConstitutiveMatrix = ZeroMatrix(3,3);

    double c1 = YoungModulus / (1.0 - PoissonCoefficient*PoissonCoefficient);
    double c2 = c1 * PoissonCoefficient;
    double c3 = 0.5*YoungModulus / (1 + PoissonCoefficient);

    rConstitutiveMatrix(0,0) = c1;
    rConstitutiveMatrix(0,1) = c2;
    rConstitutiveMatrix(1,0) = c2;
    rConstitutiveMatrix(1,1) = c1;
    rConstitutiveMatrix(2,2) = c3;
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicKirchhoffPlaneStress2D::CalculateGreenLagrangianStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector
    ) 
{
    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // E = 0.5*(inv(C) - I)
    Matrix C_tensor = prod(trans(F),F);

    rStrainVector[0] = 0.5 * ( C_tensor( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( C_tensor( 1, 1 ) - 1.00 );
    rStrainVector[2] = C_tensor( 0, 1 ); // xy
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticIsotropicKirchhoffPlaneStress2D::CalculateAlmansiStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector
    ) 
{

    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // e = 0.5*(1-inv(B))
    Matrix B_tensor = prod(F,trans(F));

    //Calculating the inverse of the jacobian
    Matrix inverse_B_tensor ( 2, 2 );
    double aux_det_b = 0;
    MathUtils<double>::InvertMatrix( B_tensor, inverse_B_tensor, aux_det_b);

    rStrainVector[0] = 0.5 * ( 1.00 - inverse_B_tensor( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.00 - inverse_B_tensor( 1, 1 ) );
    rStrainVector[2] = - inverse_B_tensor( 0, 1 ); // xy ??? yz
}

} // Namespace Kratos
