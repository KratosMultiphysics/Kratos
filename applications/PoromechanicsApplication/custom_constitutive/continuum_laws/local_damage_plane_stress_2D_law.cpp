//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/continuum_laws/local_damage_plane_stress_2D_law.hpp"

namespace Kratos
{

//Default Constructor
LocalDamagePlaneStress2DLaw::LocalDamagePlaneStress2DLaw() : LocalDamagePlaneStrain2DLaw() {}

//----------------------------------------------------------------------------------------

//Second Constructor
LocalDamagePlaneStress2DLaw::LocalDamagePlaneStress2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : LocalDamagePlaneStrain2DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
LocalDamagePlaneStress2DLaw::LocalDamagePlaneStress2DLaw(const LocalDamagePlaneStress2DLaw& rOther) : LocalDamagePlaneStrain2DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
LocalDamagePlaneStress2DLaw::~LocalDamagePlaneStress2DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LocalDamagePlaneStress2DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
	rFeatures.mOptions.Set( PLANE_STRESS_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	//rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer LocalDamagePlaneStress2DLaw::Clone() const
{
    LocalDamagePlaneStress2DLaw::Pointer p_clone(new LocalDamagePlaneStress2DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LocalDamagePlaneStress2DLaw::CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,
        const double& YoungModulus,
        const double& PoissonCoefficient )
{
    rLinearElasticMatrix.clear();
    
    // Plane stress constitutive matrix
    rLinearElasticMatrix ( 0 , 0 ) = (YoungModulus)/(1.0-PoissonCoefficient*PoissonCoefficient);
    rLinearElasticMatrix ( 1 , 1 ) = rLinearElasticMatrix ( 0 , 0 );

    rLinearElasticMatrix ( 2 , 2 ) = rLinearElasticMatrix ( 0 , 0 )*(1.0-PoissonCoefficient)*0.5;

    rLinearElasticMatrix ( 0 , 1 ) = rLinearElasticMatrix ( 0 , 0 )*PoissonCoefficient;
    rLinearElasticMatrix ( 1 , 0 ) = rLinearElasticMatrix ( 0 , 1 );
}

} // Namespace Kratos
