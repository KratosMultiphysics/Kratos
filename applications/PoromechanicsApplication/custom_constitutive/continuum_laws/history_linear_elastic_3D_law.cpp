//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// Application includes
#include "custom_constitutive/continuum_laws/history_linear_elastic_3D_law.hpp"

namespace Kratos
{
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void HistoryLinearElastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{
    KRATOS_TRY

    Flags& Options = rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();

    Vector& StrainVector = rValues.GetStrainVector();
    Vector& StressVector = rValues.GetStressVector();

    //1.- Lame constants
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
	    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
            Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
            this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );
            this->AddInitialStresses(rValues,StressVector);
	    } else {
            Matrix ConstitutiveMatrix( StrainVector.size() ,StrainVector.size());
            noalias(ConstitutiveMatrix) = ZeroMatrix( StrainVector.size() ,StrainVector.size());
            this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
            this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );
            this->AddInitialStresses(rValues,StressVector);
	    }
    } else if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
        Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
	    this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
	}

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void HistoryLinearElastic3DLaw::AddInitialStresses( Parameters& rValues, Vector& rStressVector )
{
    const Vector& N = rValues.GetShapeFunctionsValues();
    const Element::GeometryType& geometry = rValues.GetElementGeometry();
    const unsigned int number_of_nodes = geometry.size();

    unsigned int voigt_size = GetStrainSize();
    unsigned int dimension = WorkingSpaceDimension();

    Vector nodal_initial_stress_vector(voigt_size);
    Matrix nodal_initial_stress_tensor(dimension,dimension);

    Vector gp_initial_stress_vector(voigt_size);
    noalias(gp_initial_stress_vector) = ZeroVector(voigt_size);

    for (unsigned int i = 0; i < number_of_nodes; i++) {
        const Matrix& r_initial_stress_tensor = geometry[i].GetSolutionStepValue(INITIAL_STRESS_TENSOR);
        for(unsigned int j=0; j < dimension; j++) {
            for(unsigned int k=0; k < dimension; k++) {
                nodal_initial_stress_tensor(j,k) = r_initial_stress_tensor(j,k);
            }
        }
        noalias(nodal_initial_stress_vector) = MathUtils<double>::StressTensorToVector(nodal_initial_stress_tensor);

        for(unsigned int j=0; j < voigt_size; j++) {
            gp_initial_stress_vector[j] += N[i] * nodal_initial_stress_vector[j];
        }
    }

    noalias(rStressVector) += gp_initial_stress_vector;
}

} // Namespace Kratos
