//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

/* Project includes */
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress_nodal.hpp"

namespace Kratos
{

//Default Constructor
ThermalLinearElastic2DPlaneStressNodal::ThermalLinearElastic2DPlaneStressNodal() : ThermalLinearElastic2DPlaneStrainNodal() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalLinearElastic2DPlaneStressNodal::ThermalLinearElastic2DPlaneStressNodal(const ThermalLinearElastic2DPlaneStressNodal& rOther) : ThermalLinearElastic2DPlaneStrainNodal(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
ThermalLinearElastic2DPlaneStressNodal::~ThermalLinearElastic2DPlaneStressNodal() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer ThermalLinearElastic2DPlaneStressNodal::Clone() const
{
    ThermalLinearElastic2DPlaneStressNodal::Pointer p_clone(new ThermalLinearElastic2DPlaneStressNodal(*this));
    return p_clone;
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void  ThermalLinearElastic2DPlaneStressNodal::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,
        const double &rYoungModulus,
        const double &rPoissonCoefficient )
{
    rConstitutiveMatrix.clear();

    //plane stress constitutive matrix:
    rConstitutiveMatrix ( 0 , 0 ) = (rYoungModulus)/(1.0-rPoissonCoefficient*rPoissonCoefficient);
    rConstitutiveMatrix ( 1 , 1 ) = rConstitutiveMatrix ( 0 , 0 );

    rConstitutiveMatrix ( 2 , 2 ) = rConstitutiveMatrix ( 0 , 0 )*(1.0-rPoissonCoefficient)*0.5;

    rConstitutiveMatrix ( 0 , 1 ) = rConstitutiveMatrix ( 0 , 0 )*rPoissonCoefficient;
    rConstitutiveMatrix ( 1 , 0 ) = rConstitutiveMatrix ( 0 , 1 );
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic2DPlaneStressNodal::GetLawFeatures(Features& rFeatures)
{
        //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRESS_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic2DPlaneStressNodal::CalculateThermalStrain( Vector& rThermalStrainVector, const MaterialResponseVariables& rElasticVariables, double & rTemperature, double & rNodalReferenceTemperature)
{
    KRATOS_TRY


    //Identity vector
    rThermalStrainVector.resize(3,false);
    rThermalStrainVector[0] = 1.0;
    rThermalStrainVector[1] = 1.0;
    rThermalStrainVector[2] = 0.0;

    // Delta T
    double DeltaTemperature = rTemperature - rNodalReferenceTemperature;

    //Thermal strain vector
    for(unsigned int i = 0; i < 3; i++)
        rThermalStrainVector[i] *= rElasticVariables.ThermalExpansionCoefficient * DeltaTemperature;

    KRATOS_CATCH( "" )
}

} // Namespace Kratos
