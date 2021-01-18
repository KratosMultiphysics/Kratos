//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2017 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_constitutive/thermal_nonlocal_damage_plane_strain_2D_law.hpp"

namespace Kratos
{

//Default Constructor
ThermalNonlocalDamagePlaneStrain2DLaw::ThermalNonlocalDamagePlaneStrain2DLaw() : ThermalNonlocalDamage3DLaw() {}

//----------------------------------------------------------------------------------------

//Second Constructor
ThermalNonlocalDamagePlaneStrain2DLaw::ThermalNonlocalDamagePlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : ThermalNonlocalDamage3DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalNonlocalDamagePlaneStrain2DLaw::ThermalNonlocalDamagePlaneStrain2DLaw(const ThermalNonlocalDamagePlaneStrain2DLaw& rOther) : ThermalNonlocalDamage3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
ThermalNonlocalDamagePlaneStrain2DLaw::~ThermalNonlocalDamagePlaneStrain2DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalNonlocalDamagePlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
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

ConstitutiveLaw::Pointer ThermalNonlocalDamagePlaneStrain2DLaw::Clone() const
{
    ThermalNonlocalDamagePlaneStrain2DLaw::Pointer p_clone(new ThermalNonlocalDamagePlaneStrain2DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalNonlocalDamagePlaneStrain2DLaw::CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,
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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalNonlocalDamagePlaneStrain2DLaw::CalculateThermalStrain(Vector& rThermalStrainVector, const MaterialResponseVariables& ElasticVariables, double & rNodalReferenceTemperature)
{
    KRATOS_TRY

    //1.-Temperature from nodes
    const GeometryType& DomainGeometry = ElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues = ElasticVariables.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();

    double Temperature = 0.0;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        Temperature += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(TEMPERATURE);
    }

    //Identity vector
    if(rThermalStrainVector.size()!=3)
        rThermalStrainVector.resize(3,false);
    rThermalStrainVector[0] = 1.0;
    rThermalStrainVector[1] = 1.0;
    rThermalStrainVector[2] = 0.0;

    // Delta T
    double DeltaTemperature = Temperature - rNodalReferenceTemperature;

    //Thermal strain vector
    for(unsigned int i = 0; i < 3; i++)
        rThermalStrainVector[i] *= ElasticVariables.LameMu * ElasticVariables.ThermalExpansionCoefficient * DeltaTemperature;

    KRATOS_CATCH( "" )
}

} // Namespace Kratos
