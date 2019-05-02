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
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"

namespace Kratos
{

void BilinearCohesive3DLaw::GetLawFeatures(Features& rFeatures)
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

int BilinearCohesive3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo)
{
    // Verify ProcessInfo variables
    KRATOS_CHECK_VARIABLE_KEY(IS_CONVERGED);

    // Verify Properties variables
    KRATOS_CHECK_VARIABLE_KEY(CRITICAL_DISPLACEMENT);
    if(rMaterialProperties.Has(CRITICAL_DISPLACEMENT)) {
        KRATOS_ERROR_IF(rMaterialProperties[CRITICAL_DISPLACEMENT] <= 0.0) << "CRITICAL_DISPLACEMENT has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "CRITICAL_DISPLACEMENT not defined" << std::endl;
    }

    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
    if(rMaterialProperties.Has(YOUNG_MODULUS)) {
        KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "YOUNG_MODULUS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "YOUNG_MODULUS not defined" << std::endl;
    }

    KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS);
    if(rMaterialProperties.Has(YIELD_STRESS)) {
        KRATOS_ERROR_IF(rMaterialProperties[YIELD_STRESS] < 0.0) << "YIELD_STRESS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "YIELD_STRESS not defined" << std::endl;
    }

    KRATOS_CHECK_VARIABLE_KEY(FRICTION_COEFFICIENT);
    if(rMaterialProperties.Has(FRICTION_COEFFICIENT)) {
        KRATOS_ERROR_IF(rMaterialProperties[FRICTION_COEFFICIENT] < 0.0) << "FRICTION_COEFFICIENT has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "FRICTION_COEFFICIENT not defined" << std::endl;
    }

    KRATOS_CHECK_VARIABLE_KEY(DAMAGE_THRESHOLD);
    if(rMaterialProperties.Has(DAMAGE_THRESHOLD)) {
        const double& damage_threshold = rMaterialProperties[DAMAGE_THRESHOLD];
        const bool check = static_cast<bool>((damage_threshold <= 0.0) || (damage_threshold > 1.0));
        KRATOS_ERROR_IF(check) << "DAMAGE_THRESHOLD has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "DAMAGE_THRESHOLD not defined" << std::endl;
    }

    return 0;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    mStateVariable = rMaterialProperties[DAMAGE_THRESHOLD];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
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
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
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
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double& BilinearCohesive3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if( rThisVariable == DAMAGE_VARIABLE || rThisVariable == STATE_VARIABLE )
    {
        rValue = mStateVariable;
    }

    return rValue;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                        const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == STATE_VARIABLE)
    {
        mStateVariable = rValue;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables,
                                                                Parameters& rValues)

{
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    rVariables.CriticalDisplacement = MaterialProperties[CRITICAL_DISPLACEMENT];
    rVariables.DamageThreshold = MaterialProperties[DAMAGE_THRESHOLD];
    rVariables.YieldStress = MaterialProperties[YIELD_STRESS];

    rVariables.YoungModulus = MaterialProperties[YOUNG_MODULUS];
    rVariables.FrictionCoefficient = MaterialProperties[FRICTION_COEFFICIENT];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        rVariables.EquivalentStrain = std::sqrt(StrainVector[0]*StrainVector[0]+
                                                StrainVector[1]*StrainVector[1]+
                                                StrainVector[2]*StrainVector[2])/rVariables.CriticalDisplacement;
    }
    else // Contact between interfaces
    {
        rVariables.EquivalentStrain = std::sqrt(StrainVector[0]*StrainVector[0]+
                                                StrainVector[1]*StrainVector[1])/rVariables.CriticalDisplacement;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::CheckLoadingFunction(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues)
{
    rVariables.LoadingFlag = false;
    rVariables.LoadingFunction = 0.0;

    if(rVariables.EquivalentStrain >= mStateVariable)
    {
        rVariables.LoadingFlag = true;
        rVariables.LoadingFunction = 1.0;
    }
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                                        ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        if(rVariables.LoadingFlag) // Loading
        {
            rConstitutiveMatrix(0,0) = rVariables.YieldStress/((1.0-rVariables.DamageThreshold)*rVariables.CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                        StrainVector[0]*StrainVector[0]/(rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
            rConstitutiveMatrix(1,1) = rVariables.YieldStress/((1.0-rVariables.DamageThreshold)*rVariables.CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                        StrainVector[1]*StrainVector[1]/(rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
            rConstitutiveMatrix(2,2) = rVariables.YieldStress/((1.0-rVariables.DamageThreshold)*rVariables.CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                        StrainVector[2]*StrainVector[2]/(rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );

            rConstitutiveMatrix(0,1) = -rVariables.YieldStress*StrainVector[0]*StrainVector[1]/( (1.0-rVariables.DamageThreshold)*
                                        rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
            rConstitutiveMatrix(0,2) = -rVariables.YieldStress*StrainVector[0]*StrainVector[2]/( (1.0-rVariables.DamageThreshold)*
                                        rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
            rConstitutiveMatrix(1,2) = -rVariables.YieldStress*StrainVector[1]*StrainVector[2]/( (1.0-rVariables.DamageThreshold)*
                                        rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
            rConstitutiveMatrix(1,0) = rConstitutiveMatrix(0,1);
            rConstitutiveMatrix(2,0) = rConstitutiveMatrix(0,2);
            rConstitutiveMatrix(2,1) = rConstitutiveMatrix(1,2);
        }
        else // Unloading
        {
            rConstitutiveMatrix(0,0) = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold);
            rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
            rConstitutiveMatrix(2,2) = rConstitutiveMatrix(0,0);

            rConstitutiveMatrix(0,1) = 0.0;
            rConstitutiveMatrix(0,2) = 0.0;
            rConstitutiveMatrix(1,2) = 0.0;
            rConstitutiveMatrix(1,0) = 0.0;
            rConstitutiveMatrix(2,0) = 0.0;
            rConstitutiveMatrix(2,1) = 0.0;
        }
    }
    else // Contact between interfaces
    {
        if(rVariables.LoadingFlag) // Loading
        {
            rConstitutiveMatrix(0,0) = rVariables.YieldStress/((1.0-rVariables.DamageThreshold)*rVariables.CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                        StrainVector[0]*StrainVector[0]/(rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
            rConstitutiveMatrix(1,1) = rVariables.YieldStress/((1.0-rVariables.DamageThreshold)*rVariables.CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                        StrainVector[1]*StrainVector[1]/(rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
            rConstitutiveMatrix(2,2) = rVariables.YoungModulus/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);

            rConstitutiveMatrix(0,1) = -rVariables.YieldStress*StrainVector[0]*StrainVector[1]/( (1.0-rVariables.DamageThreshold)*
                                        rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
            if(StrainVector[0] > 1.0e-20)
            {
                rConstitutiveMatrix(0,2) = -rVariables.YieldStress*StrainVector[0]*StrainVector[2]/( (1.0-rVariables.DamageThreshold)*
                                            rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) -
                                            rVariables.YoungModulus*rVariables.FrictionCoefficient/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
            }
            else if(StrainVector[0] < -1.0e-20)
            {
                rConstitutiveMatrix(0,2) = -rVariables.YieldStress*StrainVector[0]*StrainVector[2]/( (1.0-rVariables.DamageThreshold)*
                                            rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) +
                                            rVariables.YoungModulus*rVariables.FrictionCoefficient/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
            }
            else
            {
                rConstitutiveMatrix(0,2) = 0.0;
            }
            if(StrainVector[1] > 1.0e-20)
            {
                rConstitutiveMatrix(1,2) = -rVariables.YieldStress*StrainVector[1]*StrainVector[2]/( (1.0-rVariables.DamageThreshold)*
                                            rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) -
                                            rVariables.YoungModulus*rVariables.FrictionCoefficient/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
            }
            else if(StrainVector[1] < -1.0e-20)
            {
                rConstitutiveMatrix(1,2) = -rVariables.YieldStress*StrainVector[1]*StrainVector[2]/( (1.0-rVariables.DamageThreshold)*
                                            rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) +
                                            rVariables.YoungModulus*rVariables.FrictionCoefficient/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
            }
            else
            {
                rConstitutiveMatrix(1,2) = 0.0;
            }
            rConstitutiveMatrix(1,0) = rConstitutiveMatrix(0,1);
            rConstitutiveMatrix(2,0) = 0.0;
            rConstitutiveMatrix(2,1) = 0.0;
        }
        else // Unloading
        {
            rConstitutiveMatrix(0,0) = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold);
            rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
            rConstitutiveMatrix(2,2) = rVariables.YoungModulus/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);

            rConstitutiveMatrix(0,1) = 0.0;
            if(StrainVector[0] > 1.0e-20)
            {
                rConstitutiveMatrix(0,2) = -rVariables.YoungModulus*rVariables.FrictionCoefficient/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
            }
            else if(StrainVector[0] < -1.0e-20)
            {
                rConstitutiveMatrix(0,2) = rVariables.YoungModulus*rVariables.FrictionCoefficient/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
            }
            else
            {
                rConstitutiveMatrix(0,2) = 0.0;
            }
            if(StrainVector[1] > 1.0e-20)
            {
                rConstitutiveMatrix(1,2) = -rVariables.YoungModulus*rVariables.FrictionCoefficient/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
            }
            else if(StrainVector[1] < -1.0e-20)
            {
                rConstitutiveMatrix(1,2) = rVariables.YoungModulus*rVariables.FrictionCoefficient/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
            }
            else
            {
                rConstitutiveMatrix(1,2) = 0.0;
            }
            rConstitutiveMatrix(1,0) = 0.0;
            rConstitutiveMatrix(2,0) = 0.0;
            rConstitutiveMatrix(2,1) = 0.0;
        }
    }
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeStressVector(Vector& rStressVector,
                                                ConstitutiveLawVariables& rVariables,
                                                Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        rStressVector[0] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold) * StrainVector[0];
        rStressVector[1] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold) * StrainVector[1];
        rStressVector[2] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold) * StrainVector[2];
    }
    else // Contact between interfaces
    {
        // Note: StrainVector[2] < 0.0
        rStressVector[2] = rVariables.YoungModulus/(rVariables.DamageThreshold*rVariables.CriticalDisplacement)*StrainVector[2];

        if(StrainVector[0] > 1.0e-20)
        {
            rStressVector[0] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold)*StrainVector[0] - rVariables.FrictionCoefficient*rStressVector[2];
        }
        else if(StrainVector[0] < -1.0e-20)
        {
            rStressVector[0] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold)*StrainVector[0] + rVariables.FrictionCoefficient*rStressVector[2];
        }
        else
        {
            rStressVector[0] = 0.0;
        }

        if(StrainVector[1] > 1.0e-20)
        {
            rStressVector[1] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold)*StrainVector[1] - rVariables.FrictionCoefficient*rStressVector[2];
        }
        else if(StrainVector[1] < -1.0e-20)
        {
            rStressVector[1] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold)*StrainVector[1] + rVariables.FrictionCoefficient*rStressVector[2];
        }
        else
        {
            rStressVector[1] = 0.0;
        }
    }
}

} // Namespace Kratos
