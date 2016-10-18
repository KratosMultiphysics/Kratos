//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/restore_simo_ju_plane_strain_2D_law.hpp"

namespace Kratos
{

//Default Constructor
RestoreSimoJuPlaneStrain2DLaw::RestoreSimoJuPlaneStrain2DLaw() : RestoreSimoJu3DLaw() {}

//----------------------------------------------------------------------------------------

//Second Constructor
RestoreSimoJuPlaneStrain2DLaw::RestoreSimoJuPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : RestoreSimoJu3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
RestoreSimoJuPlaneStrain2DLaw::RestoreSimoJuPlaneStrain2DLaw(const RestoreSimoJuPlaneStrain2DLaw& rOther) : RestoreSimoJu3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
RestoreSimoJuPlaneStrain2DLaw::~RestoreSimoJuPlaneStrain2DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void RestoreSimoJuPlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
	rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
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

ConstitutiveLaw::Pointer RestoreSimoJuPlaneStrain2DLaw::Clone() const
{
    RestoreSimoJuPlaneStrain2DLaw::Pointer p_clone(new RestoreSimoJuPlaneStrain2DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void RestoreSimoJuPlaneStrain2DLaw::CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry )
{
    //rCharacteristicSize is the diameter of a circle with the same area as the element
    rCharacteristicSize = sqrt(4.0*DomainGeometry.Area()/KRATOS_M_PI);
}

//----------------------------------------------------------------------------------------

void RestoreSimoJuPlaneStrain2DLaw::CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,
        const double& YoungModulus,
        const double& PoissonCoefficient )
{
    rLinearElasticMatrix.clear();
    
    // Plane strain constitutive matrix
    rLinearElasticMatrix ( 0 , 0 ) = (YoungModulus*(1.0-PoissonCoefficient)/((1.0+PoissonCoefficient)*(1.0-2.0*PoissonCoefficient)));
    rLinearElasticMatrix ( 1 , 1 ) = rLinearElasticMatrix ( 0 , 0 );

    rLinearElasticMatrix ( 2 , 2 ) = rLinearElasticMatrix ( 0 , 0 )*(1.0-2.0*PoissonCoefficient)/(2.0*(1.0-PoissonCoefficient));

    rLinearElasticMatrix ( 0 , 1 ) = rLinearElasticMatrix ( 0 , 0 )*PoissonCoefficient/(1.0-PoissonCoefficient);
    rLinearElasticMatrix ( 1 , 0 ) = rLinearElasticMatrix ( 0 , 1 );
}

} // Namespace Kratos
