//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"

namespace Kratos
{

//Default Constructor
BilinearCohesive3DLaw::BilinearCohesive3DLaw() : ConstitutiveLaw() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
BilinearCohesive3DLaw::BilinearCohesive3DLaw(const BilinearCohesive3DLaw& rOther) : ConstitutiveLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
BilinearCohesive3DLaw::~BilinearCohesive3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
    if(CRITICAL_DISPLACEMENT.Key() == 0 || rMaterialProperties.Has( CRITICAL_DISPLACEMENT ) == false || rMaterialProperties[CRITICAL_DISPLACEMENT] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"CRITICAL_DISPLACEMENT has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(YOUNG_MODULUS.Key() == 0 || rMaterialProperties.Has( YOUNG_MODULUS ) == false || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
        KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(RESIDUAL_STRESS.Key() == 0 || rMaterialProperties.Has( RESIDUAL_STRESS ) == false || rMaterialProperties[RESIDUAL_STRESS] < 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"RESIDUAL_STRESS has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if(FRICTION_COEFFICIENT.Key() == 0 || rMaterialProperties.Has( FRICTION_COEFFICIENT ) == false || rMaterialProperties[FRICTION_COEFFICIENT] < 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"FRICTION_COEFFICIENT has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    const double& DamageThreshold = rMaterialProperties[DAMAGE_THRESHOLD];
    if(DAMAGE_THRESHOLD.Key() == 0 || rMaterialProperties.Has( DAMAGE_THRESHOLD ) == false || DamageThreshold<=0.0 || DamageThreshold > 1.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DAMAGE_THRESHOLD has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )

    return 0;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer BilinearCohesive3DLaw::Clone() const
{
    BilinearCohesive3DLaw::Pointer p_clone(new BilinearCohesive3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    mStateVariable = rMaterialProperties[DAMAGE_THRESHOLD];
    mStateVariableEquilibrium = mStateVariable;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    //Check
    rValues.CheckAllParameters();

    //Initialize main variables
    Vector& rStrainVector = rValues.GetStrainVector();
    Vector& rStressVector = rValues.GetStressVector();
    double EffectiveDisplacement;
    
    //Material properties
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    const double& CriticalDisplacement = MaterialProperties[CRITICAL_DISPLACEMENT];
    const double& DamageThreshold = MaterialProperties[DAMAGE_THRESHOLD];
    const double& ResidualStress = MaterialProperties[RESIDUAL_STRESS];
        
    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        this->ComputeEffectiveDisplacement(EffectiveDisplacement,rStrainVector,CriticalDisplacement);

        if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR) ) //Computation of ConstitutiveMatrix and StressVector. WARNING: it is not general
        {
            Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                        
            if(EffectiveDisplacement >= mStateVariable) //Loading
            {
                mStateVariable = EffectiveDisplacement;
                if(mStateVariable > 1.0) mStateVariable = 1.0;
                
                this->ComputeConstitutiveMatrixLoading(rConstitutiveMatrix,rStrainVector,ResidualStress,DamageThreshold,CriticalDisplacement);
            }
            else //Unloading
            {
                this->ComputeConstitutiveMatrixUnloading(rConstitutiveMatrix,ResidualStress,DamageThreshold,CriticalDisplacement);
            }
            this->ComputeStressVector(rStressVector,rStrainVector,ResidualStress,DamageThreshold,CriticalDisplacement);
        }
        else //Computation of StressVector only
        {
            if(EffectiveDisplacement >= mStateVariable)
            {
                mStateVariable = EffectiveDisplacement;
                if(mStateVariable > 1.0) mStateVariable = 1.0;
            }
            this->ComputeStressVector(rStressVector,rStrainVector,ResidualStress,DamageThreshold,CriticalDisplacement);
        }
    }
    else  // Contact between interfaces
    {
        const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
        const double& FrictionCoefficient = MaterialProperties[FRICTION_COEFFICIENT];
        
        this->ComputeEffectiveDisplacementContact(EffectiveDisplacement,rStrainVector,CriticalDisplacement);
        
        if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR) ) //Computation of ConstitutiveMatrix and StressVector. WARNING: it is not general
        {
            Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            
            if(EffectiveDisplacement >= mStateVariable) //Loading
            {
                mStateVariable = EffectiveDisplacement;
                if(mStateVariable > 1.0) mStateVariable = 1.0;
                
                this->ComputeConstitutiveMatrixContactLoading(rConstitutiveMatrix,rStrainVector,YoungModulus,FrictionCoefficient,ResidualStress,DamageThreshold,CriticalDisplacement);
            }
            else //Unloading
            {
                this->ComputeConstitutiveMatrixContactUnloading(rConstitutiveMatrix,YoungModulus,FrictionCoefficient,ResidualStress,DamageThreshold,CriticalDisplacement);
            }
            this->ComputeStressVectorContact(rStressVector,rStrainVector,YoungModulus,FrictionCoefficient,ResidualStress,DamageThreshold,CriticalDisplacement);
        }
        else //Computation of StressVector only
        {
            if(EffectiveDisplacement >= mStateVariable)
            {
                mStateVariable = EffectiveDisplacement;
                if(mStateVariable > 1.0) mStateVariable = 1.0;
            }
            this->ComputeStressVectorContact(rStressVector,rStrainVector,YoungModulus,FrictionCoefficient,ResidualStress,DamageThreshold,CriticalDisplacement);
        }
    }
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    if(rValues.GetProcessInfo()[NO_CONVERGENCE]==0) //Convergence is achieved. Save equilibrium state variable
    {
        //Check
        rValues.CheckAllParameters();

        //Initialize main variables
        Vector& rStrainVector = rValues.GetStrainVector();
        double EffectiveDisplacement;
        
        //Material properties
        const double& CriticalDisplacement = rValues.GetMaterialProperties()[CRITICAL_DISPLACEMENT];
            
        if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
        {
            this->ComputeEffectiveDisplacement(EffectiveDisplacement,rStrainVector,CriticalDisplacement);
        }
        else // Contact between interfaces
        {            
            this->ComputeEffectiveDisplacementContact(EffectiveDisplacement,rStrainVector,CriticalDisplacement);
        }
            
        if(EffectiveDisplacement >= mStateVariable)
        {
            mStateVariable = EffectiveDisplacement;
            if(mStateVariable > 1.0) mStateVariable = 1.0;
        }
        
        mStateVariableEquilibrium = mStateVariable;
    }
    else // No convergence is achieved. Restore state variable to equilibrium
    {
        mStateVariable = mStateVariableEquilibrium;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double& BilinearCohesive3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if( rThisVariable == DAMAGE_VARIABLE )
    {
        rValue = mStateVariable;
    }
        
    return rValue;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeEffectiveDisplacement(double& rEffectiveDisplacement,const Vector& StrainVector,const double& CriticalDisplacement)
{
    rEffectiveDisplacement = sqrt(StrainVector[0]*StrainVector[0]+StrainVector[1]*StrainVector[1]+StrainVector[2]*StrainVector[2])/CriticalDisplacement;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeEffectiveDisplacementContact(double& rEffectiveDisplacement,const Vector& StrainVector,const double& CriticalDisplacement)
{
    rEffectiveDisplacement = sqrt(StrainVector[0]*StrainVector[0]+StrainVector[1]*StrainVector[1])/CriticalDisplacement;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeConstitutiveMatrixLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& ResidualStress,
                                                                        const double& DamageThreshold,const double& CriticalDisplacement)
{
    rConstitutiveMatrix(0,0) = ResidualStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[0]*StrainVector[0]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(1,1) = ResidualStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[1]*StrainVector[1]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(2,2) = ResidualStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[2]*StrainVector[2]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );

    rConstitutiveMatrix(0,1) = -ResidualStress*StrainVector[0]*StrainVector[1]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    rConstitutiveMatrix(0,2) = -ResidualStress*StrainVector[0]*StrainVector[2]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    rConstitutiveMatrix(1,2) = -ResidualStress*StrainVector[1]*StrainVector[2]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    rConstitutiveMatrix(1,0) = rConstitutiveMatrix(0,1);
    rConstitutiveMatrix(2,0) = rConstitutiveMatrix(0,2);
    rConstitutiveMatrix(2,1) = rConstitutiveMatrix(1,2);
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeConstitutiveMatrixContactLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                                            const double& ResidualStress,const double& DamageThreshold,const double& CriticalDisplacement)
{
    rConstitutiveMatrix(0,0) = ResidualStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[0]*StrainVector[0]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(1,1) = ResidualStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[1]*StrainVector[1]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(2,2) = YoungModulus/(DamageThreshold*CriticalDisplacement);

    rConstitutiveMatrix(0,1) = -ResidualStress*StrainVector[0]*StrainVector[1]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    rConstitutiveMatrix(0,2) = -ResidualStress*StrainVector[0]*StrainVector[2]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) +
                                YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    rConstitutiveMatrix(1,2) = -ResidualStress*StrainVector[1]*StrainVector[2]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) +
                                YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    rConstitutiveMatrix(1,0) = rConstitutiveMatrix(0,1);
    rConstitutiveMatrix(2,0) = 0.0;
    rConstitutiveMatrix(2,1) = 0.0;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeConstitutiveMatrixUnloading(Matrix& rConstitutiveMatrix,const double& ResidualStress,
                                                                        const double& DamageThreshold,const double& CriticalDisplacement)
{
    rConstitutiveMatrix(0,0) = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold);
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

void BilinearCohesive3DLaw::ComputeConstitutiveMatrixContactUnloading(Matrix& rConstitutiveMatrix,const double& YoungModulus,const double& FrictionCoefficient,
                                                                            const double& ResidualStress,const double& DamageThreshold,const double& CriticalDisplacement)
{
    rConstitutiveMatrix(0,0) = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold);
    rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
    rConstitutiveMatrix(2,2) = YoungModulus/(DamageThreshold*CriticalDisplacement);

    rConstitutiveMatrix(0,1) = 0.0;
    rConstitutiveMatrix(0,2) = YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    rConstitutiveMatrix(1,2) = YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    rConstitutiveMatrix(1,0) = 0.0; 
    rConstitutiveMatrix(2,0) = 0.0;
    rConstitutiveMatrix(2,1) = 0.0;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeStressVector(Vector& rStressVector,const Vector& StrainVector,const double& ResidualStress,
                                                        const double& DamageThreshold,const double& CriticalDisplacement)
{
    rStressVector[0] = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold) * StrainVector[0];
    rStressVector[1] = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold) * StrainVector[1];
    rStressVector[2] = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold) * StrainVector[2];
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeStressVectorContact(Vector& rStressVector,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                            const double& ResidualStress,const double& DamageThreshold,const double& CriticalDisplacement)
{
    if(StrainVector[0] > 0.0)
    {
        rStressVector[0] = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold)*StrainVector[0] -
                            FrictionCoefficient*YoungModulus/(DamageThreshold*CriticalDisplacement)*StrainVector[2];
    }
    else
    {
        rStressVector[0] = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold)*StrainVector[0] +
                            FrictionCoefficient*YoungModulus/(DamageThreshold*CriticalDisplacement)*StrainVector[2];
    }

    if(StrainVector[1] > 0.0)
    {
        rStressVector[1] = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold)*StrainVector[1] -
                            FrictionCoefficient*YoungModulus/(DamageThreshold*CriticalDisplacement)*StrainVector[2];
    }
    else
    {
        rStressVector[1] = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold)*StrainVector[1] +
                            FrictionCoefficient*YoungModulus/(DamageThreshold*CriticalDisplacement)*StrainVector[2];
    }
                    
    rStressVector[2] = YoungModulus/(DamageThreshold*CriticalDisplacement)*StrainVector[2];
}

} // Namespace Kratos
