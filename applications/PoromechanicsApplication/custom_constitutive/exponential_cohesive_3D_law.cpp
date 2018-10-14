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

void ExponentialCohesive3DLaw::GetLawFeatures(Features& rFeatures)
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

int ExponentialCohesive3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo)
{
    // Verify ProcessInfo variables
    if ( IS_CONVERGED.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"IS_CONVERGED Key is 0. Check if all applications were correctly registered.", "" )

    // Verify Properties variables
    if(CRITICAL_DISPLACEMENT.Key() == 0 || rMaterialProperties.Has( CRITICAL_DISPLACEMENT ) == false || rMaterialProperties[CRITICAL_DISPLACEMENT] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"CRITICAL_DISPLACEMENT has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties.Has( YOUNG_MODULUS ) == false || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
        KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(YIELD_STRESS.Key() == 0 || rMaterialProperties.Has( YIELD_STRESS ) == false || rMaterialProperties[YIELD_STRESS] < 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"YIELD_STRESS has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(FRICTION_COEFFICIENT.Key() == 0 || rMaterialProperties.Has( FRICTION_COEFFICIENT ) == false || rMaterialProperties[FRICTION_COEFFICIENT] < 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"FRICTION_COEFFICIENT has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    const double& DamageThreshold = rMaterialProperties[DAMAGE_THRESHOLD];
    if(DAMAGE_THRESHOLD.Key() == 0 || rMaterialProperties.Has( DAMAGE_THRESHOLD ) == false || DamageThreshold<=0.0 || DamageThreshold > 1.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DAMAGE_THRESHOLD has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )

    return 0;
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    mStateVariable = rMaterialProperties[DAMAGE_THRESHOLD];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    //Check
    rValues.CheckAllParameters();

    //Initialize main variables
    Vector& rStrainVector = rValues.GetStrainVector();
    double EquivalentStrain;

    //Material properties
    Flags& Options = rValues.GetOptions();
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    const double& CriticalDisplacement = MaterialProperties[CRITICAL_DISPLACEMENT];
    const double& DamageThreshold = MaterialProperties[DAMAGE_THRESHOLD];
    const double& YieldStress = MaterialProperties[YIELD_STRESS];

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        this->ComputeEquivalentStrain(EquivalentStrain,rStrainVector,CriticalDisplacement);

        if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
        {
            if(Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS))
            {
                // COMPUTE_CONSTITUTIVE_TENSOR
                Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

                if(EquivalentStrain >= mStateVariable) //Loading
                {
                    this->ComputeConstitutiveMatrixLoading(rConstitutiveMatrix,rStrainVector,YieldStress,DamageThreshold,CriticalDisplacement);
                }
                else //Unloading
                {
                    this->ComputeConstitutiveMatrixUnloading(rConstitutiveMatrix,YieldStress,DamageThreshold,CriticalDisplacement);
                }
            }
            else
            {
                // COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
                Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                Vector& rStressVector = rValues.GetStressVector();

                if(EquivalentStrain >= mStateVariable) //Loading
                {
                    this->ComputeConstitutiveMatrixLoading(rConstitutiveMatrix,rStrainVector,YieldStress,DamageThreshold,CriticalDisplacement);
                }
                else //Unloading
                {
                    this->ComputeConstitutiveMatrixUnloading(rConstitutiveMatrix,YieldStress,DamageThreshold,CriticalDisplacement);
                }
                this->ComputeStressVector(rStressVector,rStrainVector,YieldStress,DamageThreshold,CriticalDisplacement);
            }
        }
        else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
        {
            // COMPUTE_STRESS
            Vector& rStressVector = rValues.GetStressVector();

            this->ComputeStressVector(rStressVector,rStrainVector,YieldStress,DamageThreshold,CriticalDisplacement);
        }
    }
    else  // Contact between interfaces
    {
        const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
        const double& FrictionCoefficient = MaterialProperties[FRICTION_COEFFICIENT];

        this->ComputeEquivalentStrainContact(EquivalentStrain,rStrainVector,CriticalDisplacement);

        if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
        {
            if(Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS))
            {
                // COMPUTE_CONSTITUTIVE_TENSOR
                Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

                if(EquivalentStrain >= mStateVariable) //Loading
                {
                    //TODO
                    this->ComputeConstitutiveMatrixContactLoading(rConstitutiveMatrix,rStrainVector,YoungModulus,FrictionCoefficient,YieldStress,DamageThreshold,CriticalDisplacement);
                    // this->ComputeConstitutiveMatrixContactUnloading(rConstitutiveMatrix,rStrainVector,YoungModulus,FrictionCoefficient,YieldStress,DamageThreshold,CriticalDisplacement);
                }
                else //Unloading
                {
                    this->ComputeConstitutiveMatrixContactUnloading(rConstitutiveMatrix,rStrainVector,YoungModulus,FrictionCoefficient,YieldStress,DamageThreshold,CriticalDisplacement);
                }
            }
            else
            {
                // COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
                Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                Vector& rStressVector = rValues.GetStressVector();

                if(EquivalentStrain >= mStateVariable) //Loading
                {
                    //TODO
                    this->ComputeConstitutiveMatrixContactLoading(rConstitutiveMatrix,rStrainVector,YoungModulus,FrictionCoefficient,YieldStress,DamageThreshold,CriticalDisplacement);
                    // this->ComputeConstitutiveMatrixContactUnloading(rConstitutiveMatrix,rStrainVector,YoungModulus,FrictionCoefficient,YieldStress,DamageThreshold,CriticalDisplacement);
                }
                else //Unloading
                {
                    this->ComputeConstitutiveMatrixContactUnloading(rConstitutiveMatrix,rStrainVector,YoungModulus,FrictionCoefficient,YieldStress,DamageThreshold,CriticalDisplacement);
                }
                this->ComputeStressVectorContact(rStressVector,rStrainVector,YoungModulus,FrictionCoefficient,YieldStress,DamageThreshold,CriticalDisplacement);
            }
        }
        else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
        {
            // COMPUTE_STRESS
            Vector& rStressVector = rValues.GetStressVector();

            this->ComputeStressVectorContact(rStressVector,rStrainVector,YoungModulus,FrictionCoefficient,YieldStress,DamageThreshold,CriticalDisplacement);
        }
    }
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    if(rValues.GetProcessInfo()[IS_CONVERGED]==true) //Convergence is achieved. Save equilibrium state variable
    {
        //Check
        rValues.CheckAllParameters();

        //Initialize main variables
        Vector& rStrainVector = rValues.GetStrainVector();
        double EquivalentStrain;

        //Material properties
        const double& CriticalDisplacement = rValues.GetMaterialProperties()[CRITICAL_DISPLACEMENT];

        if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
        {
            this->ComputeEquivalentStrain(EquivalentStrain,rStrainVector,CriticalDisplacement);
        }
        else // Contact between interfaces
        {
            this->ComputeEquivalentStrainContact(EquivalentStrain,rStrainVector,CriticalDisplacement);
        }

        if(EquivalentStrain >= mStateVariable)
        {
            mStateVariable = EquivalentStrain;
            if(mStateVariable > 1.0) mStateVariable = 1.0;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double& ExponentialCohesive3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if( rThisVariable == DAMAGE_VARIABLE || rThisVariable == STATE_VARIABLE )
    {
        rValue = mStateVariable;
    }

    return rValue;
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                                        const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == STATE_VARIABLE)
    {
        mStateVariable = rValue;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeEquivalentStrain(double& rEquivalentStrain,const Vector& StrainVector,const double& CriticalDisplacement)
{
    rEquivalentStrain = sqrt(StrainVector[0]*StrainVector[0]+StrainVector[1]*StrainVector[1]+StrainVector[2]*StrainVector[2])/CriticalDisplacement;
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeEquivalentStrainContact(double& rEquivalentStrain,const Vector& StrainVector,const double& CriticalDisplacement)
{
    rEquivalentStrain = sqrt(StrainVector[0]*StrainVector[0]+StrainVector[1]*StrainVector[1])/CriticalDisplacement;
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeConstitutiveMatrixLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& YieldStress,
                                                                        const double& DamageThreshold,const double& CriticalDisplacement)
{
    rConstitutiveMatrix(0,0) = YieldStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[0]*StrainVector[0]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(1,1) = YieldStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[1]*StrainVector[1]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(2,2) = YieldStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[2]*StrainVector[2]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );

    rConstitutiveMatrix(0,1) = -YieldStress*StrainVector[0]*StrainVector[1]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    rConstitutiveMatrix(0,2) = -YieldStress*StrainVector[0]*StrainVector[2]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    rConstitutiveMatrix(1,2) = -YieldStress*StrainVector[1]*StrainVector[2]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    rConstitutiveMatrix(1,0) = rConstitutiveMatrix(0,1);
    rConstitutiveMatrix(2,0) = rConstitutiveMatrix(0,2);
    rConstitutiveMatrix(2,1) = rConstitutiveMatrix(1,2);
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeConstitutiveMatrixContactLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                                            const double& YieldStress,const double& DamageThreshold,const double& CriticalDisplacement)
{
    rConstitutiveMatrix(0,0) = YieldStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[0]*StrainVector[0]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(1,1) = YieldStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[1]*StrainVector[1]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(2,2) = YoungModulus/(DamageThreshold*CriticalDisplacement);

    rConstitutiveMatrix(0,1) = -YieldStress*StrainVector[0]*StrainVector[1]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    if(StrainVector[0] > 1.0e-20)
    {
        rConstitutiveMatrix(0,2) = -YieldStress*StrainVector[0]*StrainVector[2]/( (1.0-DamageThreshold)*
                                    CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) -
                                    YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    }
    else if(StrainVector[0] < -1.0e-20)
    {
        rConstitutiveMatrix(0,2) = -YieldStress*StrainVector[0]*StrainVector[2]/( (1.0-DamageThreshold)*
                                    CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) +
                                    YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    }
    else
    {
        rConstitutiveMatrix(0,2) = 0.0;
    }
    if(StrainVector[1] > 1.0e-20)
    {
        rConstitutiveMatrix(1,2) = -YieldStress*StrainVector[1]*StrainVector[2]/( (1.0-DamageThreshold)*
                                    CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) -
                                    YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    }
    else if(StrainVector[1] < -1.0e-20)
    {
        rConstitutiveMatrix(1,2) = -YieldStress*StrainVector[1]*StrainVector[2]/( (1.0-DamageThreshold)*
                                    CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) +
                                    YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    }
    else
    {
        rConstitutiveMatrix(1,2) = 0.0;
    }
    rConstitutiveMatrix(1,0) = rConstitutiveMatrix(0,1);
    rConstitutiveMatrix(2,0) = 0.0;
    rConstitutiveMatrix(2,1) = 0.0;
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeConstitutiveMatrixUnloading(Matrix& rConstitutiveMatrix,const double& YieldStress,
                                                                        const double& DamageThreshold,const double& CriticalDisplacement)
{
    rConstitutiveMatrix(0,0) = YieldStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold);
    rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
    rConstitutiveMatrix(2,2) = rConstitutiveMatrix(0,0);

    rConstitutiveMatrix(0,1) = 0.0;
    rConstitutiveMatrix(0,2) = 0.0;
    rConstitutiveMatrix(1,2) = 0.0;
    rConstitutiveMatrix(1,0) = 0.0;
    rConstitutiveMatrix(2,0) = 0.0;
    rConstitutiveMatrix(2,1) = 0.0;
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeConstitutiveMatrixContactUnloading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                                            const double& YieldStress,const double& DamageThreshold,const double& CriticalDisplacement)
{
    rConstitutiveMatrix(0,0) = YieldStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold);
    rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
    rConstitutiveMatrix(2,2) = YoungModulus/(DamageThreshold*CriticalDisplacement);

    rConstitutiveMatrix(0,1) = 0.0;
    if(StrainVector[0] > 1.0e-20)
    {
        rConstitutiveMatrix(0,2) = -YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    }
    else if(StrainVector[0] < -1.0e-20)
    {
        rConstitutiveMatrix(0,2) = YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    }
    else
    {
        rConstitutiveMatrix(0,2) = 0.0;
    }
    if(StrainVector[1] > 1.0e-20)
    {
        rConstitutiveMatrix(1,2) = -YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    }
    else if(StrainVector[1] < -1.0e-20)
    {
        rConstitutiveMatrix(1,2) = YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    }
    else
    {
        rConstitutiveMatrix(1,2) = 0.0;
    }
    rConstitutiveMatrix(1,0) = 0.0;
    rConstitutiveMatrix(2,0) = 0.0;
    rConstitutiveMatrix(2,1) = 0.0;
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeStressVector(Vector& rStressVector,const Vector& StrainVector,const double& YieldStress,
                                                        const double& DamageThreshold,const double& CriticalDisplacement)
{
    rStressVector[0] = YieldStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold) * StrainVector[0];
    rStressVector[1] = YieldStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold) * StrainVector[1];
    rStressVector[2] = YieldStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold) * StrainVector[2];
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeStressVectorContact(Vector& rStressVector,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                            const double& YieldStress,const double& DamageThreshold,const double& CriticalDisplacement)
{
    // Note: StrainVector[2] < 0.0
    rStressVector[2] = YoungModulus/(DamageThreshold*CriticalDisplacement)*StrainVector[2];

    if(StrainVector[0] > 1.0e-20)
    {
        rStressVector[0] = YieldStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold)*StrainVector[0] - FrictionCoefficient*rStressVector[2];
    }
    else if(StrainVector[0] < -1.0e-20)
    {
        rStressVector[0] = YieldStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold)*StrainVector[0] + FrictionCoefficient*rStressVector[2];
    }
    else
    {
        rStressVector[0] = 0.0;
    }

    if(StrainVector[1] > 1.0e-20)
    {
        rStressVector[1] = YieldStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold)*StrainVector[1] - FrictionCoefficient*rStressVector[2];
    }
    else if(StrainVector[1] < -1.0e-20)
    {
        rStressVector[1] = YieldStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold)*StrainVector[1] + FrictionCoefficient*rStressVector[2];
    }
    else
    {
        rStressVector[1] = 0.0;
    }
}

} // Namespace Kratos
