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
#include "custom_constitutive/bilinear_cohesive_2D_law.hpp"

namespace Kratos
{

void BilinearCohesive2DLaw::GetLawFeatures(Features& rFeatures)
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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    rVariables.EquivalentStrain = std::sqrt(StrainVector[0]*StrainVector[0]+StrainVector[1]*StrainVector[1])/rVariables.CriticalDisplacement;

    rVariables.LoadingFlag = false;
    if(rVariables.EquivalentStrain >= mStateVariable)
    {
        rVariables.LoadingFlag = true;
    }
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeEquivalentStrainContact(ConstitutiveLawVariables& rVariables,
                                                            Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    rVariables.EquivalentStrain = std::abs(StrainVector[0])/rVariables.CriticalDisplacement;

    rVariables.LoadingFlag = false;
    if(rVariables.EquivalentStrain >= mStateVariable)
    {
        rVariables.LoadingFlag = true;
    }
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeConstitutiveMatrixLoading(Matrix& rConstitutiveMatrix,
                                                                ConstitutiveLawVariables& rVariables,
                                                                Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    rConstitutiveMatrix(0,0) = rVariables.YieldStress/((1.0-rVariables.DamageThreshold)*rVariables.CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[0]*StrainVector[0]/(rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(1,1) = rVariables.YieldStress/((1.0-rVariables.DamageThreshold)*rVariables.CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[1]*StrainVector[1]/(rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );

    rConstitutiveMatrix(0,1) = -rVariables.YieldStress*StrainVector[0]*StrainVector[1]/( (1.0-rVariables.DamageThreshold)*
                                rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    rConstitutiveMatrix(1,0) = rConstitutiveMatrix(0,1);
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeConstitutiveMatrixContactLoading(Matrix& rConstitutiveMatrix,
                                                                        ConstitutiveLawVariables& rVariables,
                                                                        Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    rConstitutiveMatrix(0,0) = rVariables.YieldStress/((1.0-rVariables.DamageThreshold)*rVariables.CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[0]*StrainVector[0]/(rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(1,1) = rVariables.YoungModulus/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);

    if(StrainVector[0] > 1.0e-20)
    {
        rConstitutiveMatrix(0,1) = -rVariables.YieldStress*StrainVector[0]*StrainVector[1]/( (1.0-rVariables.DamageThreshold)*
                                    rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) -
                                    rVariables.YoungModulus*rVariables.FrictionCoefficient/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
    }
    else if(StrainVector[0] < -1.0e-20)
    {
        rConstitutiveMatrix(0,1) = -rVariables.YieldStress*StrainVector[0]*StrainVector[1]/( (1.0-rVariables.DamageThreshold)*
                                    rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) +
                                    rVariables.YoungModulus*rVariables.FrictionCoefficient/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
    }
    else
    {
        rConstitutiveMatrix(0,1) = 0.0;
    }

    rConstitutiveMatrix(1,0) = 0.0;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeConstitutiveMatrixUnloading(Matrix& rConstitutiveMatrix,
                                                                ConstitutiveLawVariables& rVariables,
                                                                Parameters& rValues)
{
    rConstitutiveMatrix(0,0) = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold);
    rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);

    rConstitutiveMatrix(0,1) = 0.0;
    rConstitutiveMatrix(1,0) = 0.0;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeConstitutiveMatrixContactUnloading(Matrix& rConstitutiveMatrix,
                                                                        ConstitutiveLawVariables& rVariables,
                                                                        Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    rConstitutiveMatrix(0,0) = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold);
    rConstitutiveMatrix(1,1) = rVariables.YoungModulus/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);

    if(StrainVector[0] > 0.0)
    {
        rConstitutiveMatrix(0,1) = -rVariables.YoungModulus*rVariables.FrictionCoefficient/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
    }
    else if(StrainVector[0] < 0.0)
    {
        rConstitutiveMatrix(0,1) = rVariables.YoungModulus*rVariables.FrictionCoefficient/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
    }
    else
    {
        rConstitutiveMatrix(0,1) = 0.0;
    }

    rConstitutiveMatrix(1,0) = 0.0;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeStressVector(Vector& rStressVector,
                                                ConstitutiveLawVariables& rVariables,
                                                Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    rStressVector[0] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold) * StrainVector[0];
    rStressVector[1] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold) * StrainVector[1];
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeStressVectorContact(Vector& rStressVector,
                                                        ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    // Note: StrainVector[1] < 0.0
    rStressVector[1] = rVariables.YoungModulus/(rVariables.DamageThreshold*rVariables.CriticalDisplacement)*StrainVector[1];

    if(StrainVector[0] > 0.0)
    {
        rStressVector[0] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold)*StrainVector[0] - rVariables.FrictionCoefficient*rStressVector[1];
    }
    else if(StrainVector[0] < 0.0)
    {
        rStressVector[0] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold)*StrainVector[0] + rVariables.FrictionCoefficient*rStressVector[1];
    }
    else
    {
        rStressVector[0] = 0.0;
    }

    //TODO

    // rStressVector[1] = -0.5*rVariables.YoungModulus*std::exp(-StrainVector[1]); // Kn = 0.5E

    // // rStressVector[1] = rVariables.YoungModulus/(rVariables.DamageThreshold*rVariables.CriticalDisplacement)*StrainVector[1];

    // rStressVector[0] = 0.5*rVariables.YoungModulus*StrainVector[0]; // Kt = 0.5E

    // double Tauc = std::abs(rStressVector[1])*rVariables.FrictionCoefficient;
    // double Taum = (1.0-mStateVariable)*rVariables.YieldStress/std::sqrt(3.0); // n=(1.0-mStateVariable)?

    // // ShearStrength = min(Tauc,Taum)
    // double ShearStrength = Tauc;
    // if (Taum < ShearStrength)
    // {
    //     ShearStrength = Taum;
    // }

    // double TangentialStress = std::abs(rStressVector[0]);

    // if(TangentialStress >= ShearStrength)
    // {
    //     if(StrainVector[0] > 1.0e-20)
    //     {
    //         rStressVector[0] = ShearStrength;
    //     }
    //     else if(StrainVector[0] < -1.0e-20)
    //     {
    //         rStressVector[0] = -1.0*ShearStrength;
    //     }
    //     else
    //     {
    //         rStressVector[0] = 0.0;
    //     }
    // }

}

} // Namespace Kratos
