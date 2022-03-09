// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alireza Taherzadeh Fard
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/cohesive_law.h"

#include "constitutive_laws_application_variables.h"

#include "utilities/math_utils.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

CohesiveLaw::CohesiveLaw()
    : ElasticIsotropic3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

CohesiveLaw::CohesiveLaw(const CohesiveLaw& rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer CohesiveLaw::Clone() const
{
    return Kratos::make_shared<CohesiveLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

CohesiveLaw::~CohesiveLaw()
{
}


void  CohesiveLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    const Vector& r_strain_vector = rValues.GetStrainVector();

    Vector& r_stress_vector = rValues.GetStressVector();

    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E  = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    const double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
    const double c2 = c1 * ( 1 - NU );
    const double c3 = c1 * NU;
    const double c4 = c1 * 0.5 * ( 1 - 2 * NU );

    r_constitutive_matrix( 0, 0 ) = c2;
    r_constitutive_matrix( 0, 1 ) = c3;
    r_constitutive_matrix( 0, 2 ) = c3;
    r_constitutive_matrix( 1, 0 ) = c3;
    r_constitutive_matrix( 1, 1 ) = c2;
    r_constitutive_matrix( 1, 2 ) = c3;
    r_constitutive_matrix( 2, 0 ) = c3;
    r_constitutive_matrix( 2, 1 ) = c3;
    r_constitutive_matrix( 2, 2 ) = c2;
    r_constitutive_matrix( 3, 3 ) = c4;
    r_constitutive_matrix( 4, 4 ) = c4;
    r_constitutive_matrix( 5, 5 ) = c4;

    double H_variable = mH;
    double rD_variable  = mrD;
    double D_variable = mD;
    double Hmax_variable = mHmax;
   
    
    // We compute H
    H_variable = std::min(((r_strain_vector[2]-0.02) / (0.042-0.02)),1.0);
    
    // We store the maximum amount of H over time
    if(H_variable > Hmax_variable){
        Hmax_variable = H_variable;
    }
    
    // We compute rD
    rD_variable = std::max(0.0,Hmax_variable);

    // We compute damage
    D_variable = (rD_variable*0.042) / ((rD_variable*0.042)+((1-rD_variable)*0.042));
    

    
    r_stress_vector[0] = (c2 * r_strain_vector[0] + c3 * r_strain_vector[1] + c3 * r_strain_vector[2]);
    r_stress_vector[1] = (c3 * r_strain_vector[0] + c2 * r_strain_vector[1] + c3 * r_strain_vector[2]);
    r_stress_vector[2] = (-1)*D_variable*(c3 * r_strain_vector[0] + c3 * r_strain_vector[1] + c2 * r_strain_vector[2]);
    r_stress_vector[3] = (c4 * r_strain_vector[3]);
    r_stress_vector[4] = (c4 * r_strain_vector[4]);
    r_stress_vector[5] = (c4 * r_strain_vector[5]);

    // rStressVector[0] = c2 * rStrainVector[0] + c3 * rStrainVector[1] + c3 * rStrainVector[2];
    // rStressVector[1] = c3 * rStrainVector[0] + c2 * rStrainVector[1] + c3 * rStrainVector[2];
    // rStressVector[2] = c3 * rStrainVector[0] + c3 * rStrainVector[1] + c2 * rStrainVector[2];
    // rStressVector[3] = c4 * rStrainVector[3];
    // rStressVector[4] = c4 * rStrainVector[4];
    // rStressVector[5] = c4 * rStrainVector[5];
    

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void CohesiveLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{

const Vector& r_strain_vector = rValues.GetStrainVector();

    Vector& r_stress_vector = rValues.GetStressVector();

    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E  = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    const double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
    const double c2 = c1 * ( 1 - NU );
    const double c3 = c1 * NU;
    const double c4 = c1 * 0.5 * ( 1 - 2 * NU );

    r_constitutive_matrix( 0, 0 ) = c2;
    r_constitutive_matrix( 0, 1 ) = c3;
    r_constitutive_matrix( 0, 2 ) = c3;
    r_constitutive_matrix( 1, 0 ) = c3;
    r_constitutive_matrix( 1, 1 ) = c2;
    r_constitutive_matrix( 1, 2 ) = c3;
    r_constitutive_matrix( 2, 0 ) = c3;
    r_constitutive_matrix( 2, 1 ) = c3;
    r_constitutive_matrix( 2, 2 ) = c2;
    r_constitutive_matrix( 3, 3 ) = c4;
    r_constitutive_matrix( 4, 4 ) = c4;
    r_constitutive_matrix( 5, 5 ) = c4;

    double H_variable = mH;
    double rD_variable  = mrD;
    double D_variable = mD;
    double Hmax_variable = mHmax;

    // We compute H
    H_variable = std::min(((r_strain_vector[2]-0.02) / (0.042-0.02)),1.0);
    
    // We store the maximum amount of H over time
    if(H_variable > Hmax_variable){
        Hmax_variable = H_variable;
    }
    
    // We compute rD
    rD_variable = std::max(0.0,Hmax_variable);

    // We compute damage
    D_variable = (rD_variable*0.042) / ((rD_variable*0.042)+((1-rD_variable)*0.042));
    

    
    r_stress_vector[0] = (c2 * r_strain_vector[0] + c3 * r_strain_vector[1] + c3 * r_strain_vector[2]);
    r_stress_vector[1] = (c3 * r_strain_vector[0] + c2 * r_strain_vector[1] + c3 * r_strain_vector[2]);
    r_stress_vector[2] = (-1)*D_variable*(c3 * r_strain_vector[0] + c3 * r_strain_vector[1] + c2 * r_strain_vector[2]);
    r_stress_vector[3] = (c4 * r_strain_vector[3]);
    r_stress_vector[4] = (c4 * r_strain_vector[4]);
    r_stress_vector[5] = (c4 * r_strain_vector[5]);

    mH = H_variable;
    mrD = rD_variable;
    mD = D_variable;
    mHmax = Hmax_variable;
    
    KRATOS_WATCH(mD);
    
}

} // Namespace Kratos
