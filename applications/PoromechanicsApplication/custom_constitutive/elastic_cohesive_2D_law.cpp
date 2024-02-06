//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana and Danilo Cavalcanti
//

// Application includes
#include "custom_constitutive/elastic_cohesive_2D_law.hpp"

namespace Kratos
{

void ElasticCohesive2DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
	rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	//rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the spacedimension
	rFeatures.mSpaceDimension = 2;

	//Set the strain size
	rFeatures.mStrainSize = 2;
}

//----------------------------------------------------------------------------------------

void ElasticCohesive2DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                                                ConstitutiveLawVariables& rVariables,
                                                                Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    double cp = 1.0;
    // Penalization coefficient, in case it is a compression
    if(StrainVector[1] < 0.0){
        cp = rVariables.PenaltyStiffness;
    }

    // Fill the constitutive matrix
    noalias(rConstitutiveMatrix) = ZeroMatrix(2,2);
    rConstitutiveMatrix(0,0) = rVariables.ShearStiffness;
    rConstitutiveMatrix(1,1) = rVariables.NormalStiffness * cp;

}

//----------------------------------------------------------------------------------------

void ElasticCohesive2DLaw::ComputeStressVector(Vector& rStressVector,
                                                ConstitutiveLawVariables& rVariables,
                                                Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

// ------------------------------------------------------ NEW ------------------------------------------------------
    const Vector& N = rValues.GetShapeFunctionsValues();
    const Element::GeometryType& geometry = rValues.GetElementGeometry();
    const unsigned int number_of_nodes = geometry.size();

    // const unsigned int dimension = rValues.GetSpaceDimension();
    // const unsigned int strainSize = rValues.GetStrainSize();

    const int dimension = 2; // Temporary solution. Need to generalise
    const int strainSize = 3; // Temporary solution. Need to generalise

    // Create the necessary components for the initialisation of the stresses
    Vector nodal_initial_stress_vector(strainSize);
    Matrix nodal_initial_stress_tensor(dimension,dimension);

    Vector gp_initial_stress_vector(strainSize);
    noalias(gp_initial_stress_vector) = ZeroVector(strainSize);
    
    for (unsigned int i = 0; i < number_of_nodes; i++) {
        const Matrix& r_initial_stress_tensor = geometry[i].GetSolutionStepValue(INITIAL_STRESS_TENSOR);
        for(unsigned int j=0; j < dimension; j++) {
            for(unsigned int k=0; k < dimension; k++) {
                nodal_initial_stress_tensor(j,k) = r_initial_stress_tensor(j,k);
            }
        }
        noalias(nodal_initial_stress_vector) = MathUtils<double>::StressTensorToVector(nodal_initial_stress_tensor);

        for(unsigned int j=0; j < strainSize; j++){
            gp_initial_stress_vector[j] += N[i] * nodal_initial_stress_vector[j];
        }
    }

    //Define mid-plane points for quadrilateral_interface_2d_4
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    noalias(pmid0) = 0.5 * (geometry.GetPoint( 0 ) + geometry.GetPoint( 3 ));
    noalias(pmid1) = 0.5 * (geometry.GetPoint( 1 ) + geometry.GetPoint( 2 ));

    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = pmid1 - pmid0;
    double inv_norm_x = 1.0/norm_2(Vx);
    Vx[0] *= inv_norm_x; // cos
    Vx[1] *= inv_norm_x; // sin

    // Define and assign the rotation matrix
    Matrix RotationInterface(dimension,strainSize);
    RotationInterface(0,0) = -Vx[1]*Vx[0];
    RotationInterface(0,1) = Vx[1]*Vx[0];
    RotationInterface(0,2) = Vx[0]*Vx[0] - Vx[1]*Vx[1];
    RotationInterface(1,0) = Vx[1]*Vx[1];
    RotationInterface(1,1) = Vx[0]*Vx[0];
    RotationInterface(1,2) = -2*Vx[1]*Vx[0];

    // Local stress vector in local coordinates
    Vector LocalInitialStresses(dimension); 
    noalias(LocalInitialStresses) = prod(RotationInterface, trans(gp_initial_stress_vector));

// -----------------------------------------------------------------------------------------------------------------
    double cp = 1.0;

    // Penalization coefficient, in case it is a compression
    if(StrainVector[1] < 0.0){
        cp = rVariables.PenaltyStiffness;
    }

    rStressVector[0] = LocalInitialStresses[0] + StrainVector[0] * rVariables.ShearStiffness;
    rStressVector[1] = LocalInitialStresses[1] + StrainVector[1] * rVariables.NormalStiffness * cp;

}

} // Namespace Kratos
