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
#include "custom_constitutive/interface_element_laws/elastic_cohesive_3D_law.hpp"

namespace Kratos
{

void ElasticCohesive3DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	//rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the spacedimension
	rFeatures.mSpaceDimension = 3;

	//Set the strain size
	rFeatures.mStrainSize = 3;
}

//----------------------------------------------------------------------------------------

int ElasticCohesive3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo) const
{

    // Verify Properties variables

    if(rMaterialProperties.Has(NORMAL_STIFFNESS)) {
        KRATOS_ERROR_IF(rMaterialProperties[NORMAL_STIFFNESS] <= 0.0) << "NORMAL_STIFFNESS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "NORMAL_STIFFNESS not defined" << std::endl;
    }

    if(rMaterialProperties.Has(SHEAR_STIFFNESS)) {
        KRATOS_ERROR_IF(rMaterialProperties[SHEAR_STIFFNESS] <= 0.0) << "SHEAR_STIFFNESS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "SHEAR_STIFFNESS not defined" << std::endl;
    }

    if(rMaterialProperties.Has(PENALTY_STIFFNESS)) {
        KRATOS_ERROR_IF(rMaterialProperties[PENALTY_STIFFNESS] <= 0.0) << "PENALTY has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "PENALTY not defined" << std::endl;
    }

    return 0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ElasticCohesive3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    //Check
    rValues.CheckAllParameters();

    //Initialize main variables
    Flags& Options = rValues.GetOptions();
    ConstitutiveLawVariables Variables;
    this->InitializeConstitutiveLawVariables(Variables,rValues);

    // Compute the elastic constitutive matrix
    Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
    this->ComputeConstitutiveMatrix(rConstitutiveMatrix,Variables,rValues);

    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {
        if(Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS))
        {
            // COMPUTE_CONSTITUTIVE_TENSOR
            Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

            this->ComputeConstitutiveMatrix(rConstitutiveMatrix,Variables,rValues);
        }
        else
        {
            // COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
            Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            Vector& rStressVector = rValues.GetStressVector();

            this->ComputeConstitutiveMatrix(rConstitutiveMatrix,Variables,rValues);
            this->ComputeStressVector(rStressVector,Variables,rValues);
        }
    }
    else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
    {
        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();

        this->ComputeStressVector(rStressVector,Variables,rValues);
    }
}

//----------------------------------------------------------------------------------------
// This is not being used, can we exclude it?
void ElasticCohesive3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    if(rValues.GetProcessInfo()[IS_CONVERGED]==true) //Convergence is achieved. Save equilibrium state variable
    {
        rValues.CheckAllParameters();

        ConstitutiveLawVariables Variables;
        this->InitializeConstitutiveLawVariables(Variables,rValues);
        
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ElasticCohesive3DLaw::InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables,
                                                                Parameters& rValues)

{
    const Properties& MaterialProperties = rValues.GetMaterialProperties();

    rVariables.ShearStiffness   = MaterialProperties[SHEAR_STIFFNESS];
    rVariables.NormalStiffness  = MaterialProperties[NORMAL_STIFFNESS];
    rVariables.PenaltyStiffness = MaterialProperties[PENALTY_STIFFNESS];
}

//----------------------------------------------------------------------------------------

void ElasticCohesive3DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                                        ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    double cp = 1.0;
    // Penalization coefficient, in case it is a compression
    if(StrainVector[2] < 1.0e-20){
        cp = rVariables.PenaltyStiffness;
    }

    // Fill the constitutive matrix
    noalias(rConstitutiveMatrix) = ZeroMatrix(3,3);
    rConstitutiveMatrix(0,0) = rVariables.ShearStiffness;
    rConstitutiveMatrix(1,1) = rVariables.ShearStiffness;
    rConstitutiveMatrix(2,2) = rVariables.NormalStiffness * cp;
    
}

//----------------------------------------------------------------------------------------

void ElasticCohesive3DLaw::ComputeStressVector(Vector& rStressVector,
                                                ConstitutiveLawVariables& rVariables,
                                                Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    double cp = 1.0;
    // Penalization coefficient, in case it is a compression
    if(StrainVector[2] < 1.0e-20){
        cp = rVariables.PenaltyStiffness;
    }

    rStressVector[0] = StrainVector[0] * rVariables.ShearStiffness;
    rStressVector[1] = StrainVector[1] * rVariables.ShearStiffness;
    rStressVector[2] = StrainVector[2] * rVariables.NormalStiffness * cp;

}

} // Namespace Kratos
