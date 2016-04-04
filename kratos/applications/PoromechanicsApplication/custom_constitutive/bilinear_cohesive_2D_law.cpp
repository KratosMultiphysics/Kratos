//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/bilinear_cohesive_2D_law.hpp"

namespace Kratos
{

//Default Constructor
BilinearCohesive2DLaw::BilinearCohesive2DLaw() : BilinearCohesive3DLaw() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
BilinearCohesive2DLaw::BilinearCohesive2DLaw(const BilinearCohesive2DLaw& rOther) : BilinearCohesive3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
BilinearCohesive2DLaw::~BilinearCohesive2DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer BilinearCohesive2DLaw::Clone() const
{
    BilinearCohesive2DLaw::Pointer p_clone(new BilinearCohesive2DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeEffectiveDisplacement(double& rEffectiveDisplacement,const Vector& StrainVector,const double& CriticalDisplacement)
{
    rEffectiveDisplacement = sqrt(StrainVector[0]*StrainVector[0]+StrainVector[1]*StrainVector[1])/CriticalDisplacement;    
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeEffectiveDisplacementContact(double& rEffectiveDisplacement,const Vector& StrainVector,const double& CriticalDisplacement)
{
    rEffectiveDisplacement = fabs(StrainVector[0])/CriticalDisplacement;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeConstitutiveMatrixLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& ResidualStress,
                                                                        const double& DamageThreshold,const double& CriticalDisplacement)
{
    rConstitutiveMatrix(0,0) = ResidualStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[0]*StrainVector[0]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(1,1) = ResidualStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[1]*StrainVector[1]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );

    rConstitutiveMatrix(0,1) = -ResidualStress*StrainVector[0]*StrainVector[1]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    rConstitutiveMatrix(1,0) = rConstitutiveMatrix(0,1);
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeConstitutiveMatrixContactLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                                            const double& ResidualStress,const double& DamageThreshold,const double& CriticalDisplacement)
{
    rConstitutiveMatrix(0,0) = ResidualStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[0]*StrainVector[0]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(1,1) = YoungModulus/(DamageThreshold*CriticalDisplacement);

    rConstitutiveMatrix(0,1) = -ResidualStress*StrainVector[0]*StrainVector[1]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable ) +
                                YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    rConstitutiveMatrix(1,0) = 0.0;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeConstitutiveMatrixUnloading(Matrix& rConstitutiveMatrix,const double& ResidualStress,
                                                                        const double& DamageThreshold,const double& CriticalDisplacement)
{
    rConstitutiveMatrix(0,0) = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold);
    rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);

    rConstitutiveMatrix(0,1) = 0.0;
    rConstitutiveMatrix(1,0) = 0.0;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeConstitutiveMatrixContactUnloading(Matrix& rConstitutiveMatrix,const double& YoungModulus,const double& FrictionCoefficient,
                                                                            const double& ResidualStress,const double& DamageThreshold,const double& CriticalDisplacement)
{
    rConstitutiveMatrix(0,0) = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold);
    rConstitutiveMatrix(1,1) = YoungModulus/(DamageThreshold*CriticalDisplacement);

    rConstitutiveMatrix(0,1) = YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
    rConstitutiveMatrix(1,0) = 0.0;
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeStressVector(Vector& rStressVector,const Vector& StrainVector,const double& ResidualStress,
                                                        const double& DamageThreshold,const double& CriticalDisplacement)
{
    rStressVector[0] = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold) * StrainVector[0];
    rStressVector[1] = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold) * StrainVector[1];
}

//----------------------------------------------------------------------------------------

void BilinearCohesive2DLaw::ComputeStressVectorContact(Vector& rStressVector,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                            const double& ResidualStress,const double& DamageThreshold,const double& CriticalDisplacement)
{
    if(StrainVector[0] > 0.0)
    {
        rStressVector[0] = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold)*StrainVector[0] -
                            FrictionCoefficient*YoungModulus/(DamageThreshold*CriticalDisplacement)*StrainVector[1];
    }
    else
    {
        rStressVector[0] = ResidualStress/(CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-DamageThreshold)*StrainVector[0] +
                            FrictionCoefficient*YoungModulus/(DamageThreshold*CriticalDisplacement)*StrainVector[1];
    }
                            
    rStressVector[1] = YoungModulus/(DamageThreshold*CriticalDisplacement)*StrainVector[1];
}

} // Namespace Kratos
