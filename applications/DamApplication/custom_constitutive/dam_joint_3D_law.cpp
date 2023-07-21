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
#include "custom_constitutive/dam_joint_3D_law.hpp"

namespace Kratos
{

void DamJoint3DLaw::GetLawFeatures(Features& rFeatures)
{
    KRATOS_TRY

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

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

int DamJoint3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Verify Properties variables
    if(rMaterialProperties.Has(YOUNG_MODULUS)) {
        KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "YOUNG_MODULUS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "YOUNG_MODULUS not defined" << std::endl;
    }

    if(rMaterialProperties.Has(YIELD_STRESS)) {
        KRATOS_ERROR_IF(rMaterialProperties[YIELD_STRESS] < 0.0) << "YIELD_STRESS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "YIELD_STRESS not defined" << std::endl;
    }

    return 0;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

void DamJoint3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    KRATOS_TRY

    mStateVariable = 0.0;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void DamJoint3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    KRATOS_TRY

    //Check
    rValues.CheckAllParameters();

    //Initialize main variables
    Flags& Options = rValues.GetOptions();
    ConstitutiveLawVariables Variables;
    this->InitializeConstitutiveLawVariables(Variables,rValues);

    this->ComputeEquivalentStrain(Variables,rValues);
    this->CheckLoadingFunction(Variables,rValues);

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

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

void DamJoint3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    KRATOS_TRY

    if(rValues.GetProcessInfo()[IS_CONVERGED]==true) //Convergence is achieved. Save equilibrium state variable
    {
        rValues.CheckAllParameters();

        ConstitutiveLawVariables Variables;
        this->InitializeConstitutiveLawVariables(Variables,rValues);

        this->ComputeEquivalentStrain(Variables,rValues);
        this->CheckLoadingFunction(Variables,rValues);

        if(Variables.LoadingFlag)
        {
            mStateVariable = Variables.EquivalentStrain;
            if(mStateVariable > 1.0)
                mStateVariable = 1.0;
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double& DamJoint3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    KRATOS_TRY

    if( rThisVariable == DAMAGE_VARIABLE || rThisVariable == STATE_VARIABLE )
    {
        rValue = mStateVariable;
    }

    return rValue;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

void DamJoint3DLaw::SetValue(const Variable<double>& rThisVariable, const double& rValue,
                             const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if (rThisVariable == STATE_VARIABLE)
    {
        mStateVariable = rValue;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void DamJoint3DLaw::InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables,
                                                       Parameters& rValues)

{
    KRATOS_TRY

    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    rVariables.YieldStress = MaterialProperties[YIELD_STRESS];
    rVariables.YoungModulus = MaterialProperties[YOUNG_MODULUS];

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void DamJoint3DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                            Parameters& rValues)
{}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void DamJoint3DLaw::CheckLoadingFunction(ConstitutiveLawVariables& rVariables,
                                         Parameters& rValues)
{}

//----------------------------------------------------------------------------------------

void DamJoint3DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                              ConstitutiveLawVariables& rVariables,
                                              Parameters& rValues)
{
    KRATOS_TRY

    if(rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY)) // No contact between interfaces
    {
        rConstitutiveMatrix(0,0) = rVariables.YieldStress;
        rConstitutiveMatrix(1,1) = rVariables.YieldStress;
        rConstitutiveMatrix(2,2) = rVariables.YieldStress;

        rConstitutiveMatrix(0,1) = 0.0;
        rConstitutiveMatrix(0,2) = 0.0;
        rConstitutiveMatrix(1,2) = 0.0;
        rConstitutiveMatrix(1,0) = 0.0;
        rConstitutiveMatrix(2,0) = 0.0;
        rConstitutiveMatrix(2,1) = 0.0;
    }

    else // Contact between interfaces
    {
        rConstitutiveMatrix(0,0) = rVariables.YieldStress;
        rConstitutiveMatrix(1,1) = rVariables.YieldStress;
        rConstitutiveMatrix(2,2) = rVariables.YoungModulus;

        rConstitutiveMatrix(0,1) = 0.0;
        rConstitutiveMatrix(0,2) = 0.0;
        rConstitutiveMatrix(1,2) = 0.0;
        rConstitutiveMatrix(1,0) = 0.0;
        rConstitutiveMatrix(2,0) = 0.0;
        rConstitutiveMatrix(2,1) = 0.0;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

void DamJoint3DLaw::ComputeStressVector(Vector& rStressVector,
                                        ConstitutiveLawVariables& rVariables,
                                        Parameters& rValues)
{
    KRATOS_TRY

    const Vector& StrainVector = rValues.GetStrainVector();

    if(rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY)) // No contact between interfaces
    {
        rStressVector[0] = rVariables.YieldStress * StrainVector[0];
        rStressVector[1] = rVariables.YieldStress * StrainVector[1];
        rStressVector[2] = rVariables.YieldStress * StrainVector[2];
    }
    else // Contact between interfaces
    {
        rStressVector[0] = rVariables.YieldStress * StrainVector[0];
        rStressVector[1] = rVariables.YieldStress * StrainVector[1];
        rStressVector[2] = rVariables.YoungModulus * StrainVector[2];
    }

    KRATOS_CATCH("")
}

} // Namespace Kratos
