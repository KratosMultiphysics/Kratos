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
#include "custom_constitutive/exponential_cohesive_3D_law.hpp"

namespace Kratos
{

int ExponentialCohesive3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo)
{
    // Verify ProcessInfo variables
    KRATOS_CHECK_VARIABLE_KEY(IS_CONVERGED);

    // Verify Properties variables
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

    KRATOS_CHECK_VARIABLE_KEY(FRACTURE_ENERGY);
    if(rMaterialProperties.Has(FRACTURE_ENERGY)) {
        KRATOS_ERROR_IF(rMaterialProperties[FRACTURE_ENERGY] <= 0.0) << "FRACTURE_ENERGY has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "FRACTURE_ENERGY not defined" << std::endl;
    }

    KRATOS_CHECK_VARIABLE_KEY(SHEAR_FRACTURE_ENERGY);
    if(rMaterialProperties.Has(SHEAR_FRACTURE_ENERGY)) {
        KRATOS_ERROR_IF(rMaterialProperties[SHEAR_FRACTURE_ENERGY] < 0.0) << "SHEAR_FRACTURE_ENERGY has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "SHEAR_FRACTURE_ENERGY not defined" << std::endl;
    }

    return 0;
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    mStateVariable = 1.0e-12;
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
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
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double& ExponentialCohesive3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if(rThisVariable == DAMAGE_VARIABLE)
    {
        double Damage = 0.0;

        // TODO: Calculate damage

        rValue = Damage;
    }
    else if(rThisVariable == STATE_VARIABLE )
    {
        rValue = mStateVariable;
    }

    return rValue;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables,
                                                                    Parameters& rValues)

{
    // const GeometryType& ElementGeometry = rValues.GetElementGeometry();
    // const unsigned int Dim = ElementGeometry.WorkingSpaceDimension();
    const Vector& StrainVector = rValues.GetStrainVector();
    const Properties& MaterialProperties = rValues.GetMaterialProperties();

    rVariables.YoungModulus = MaterialProperties[YOUNG_MODULUS];
    rVariables.YieldStress = MaterialProperties[YIELD_STRESS];
    rVariables.FractureEnergy = MaterialProperties[FRACTURE_ENERGY];
    rVariables.ShearFractureEnergy = MaterialProperties[SHEAR_FRACTURE_ENERGY];

    rVariables.CriticalDisplacement = rVariables.FractureEnergy / (std::exp(1.0) * rVariables.YieldStress);
    double ShearStrain2 = StrainVector[0]*StrainVector[0]+StrainVector[1]*StrainVector[1];
    double PositiveNormalStrain = this->MacaulayBrackets(StrainVector[2]);
    double TotalStrain2 = PositiveNormalStrain*PositiveNormalStrain+ShearStrain2;
    rVariables.ModeMixingRatio = ShearStrain2 / TotalStrain2;
    rVariables.CurveFittingParameter = 1.0; // TODO

    // TODO: seguir

    rVariables.CompressionMatrix.resize(3,3);
    noalias(rVariables.CompressionMatrix) = ZeroMatrix(3,3);


    rVariables.WeightMatrix.resize(3,3);
    noalias(rVariables.WeightMatrix) = ZeroMatrix(3,3);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double ExponentialCohesive3DLaw::MacaulayBrackets(const double& Value)
{
    if(Value > 0.0)
        return Value;
    else
        return 0.0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues)
{
    // TODO

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

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                                        ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues)
{
    // TODO

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

void ExponentialCohesive3DLaw::ComputeStressVector(Vector& rStressVector,
                                                ConstitutiveLawVariables& rVariables,
                                                Parameters& rValues)
{
    // TODO

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
