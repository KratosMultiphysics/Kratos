//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

/* Project includes */
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain_nodal.hpp"

namespace Kratos
{

//Default Constructor
ThermalLinearElastic2DPlaneStrainNodal::ThermalLinearElastic2DPlaneStrainNodal() : ThermalLinearElastic3DLawNodal() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalLinearElastic2DPlaneStrainNodal::ThermalLinearElastic2DPlaneStrainNodal(const ThermalLinearElastic2DPlaneStrainNodal& rOther) : ThermalLinearElastic3DLawNodal(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
ThermalLinearElastic2DPlaneStrainNodal::~ThermalLinearElastic2DPlaneStrainNodal() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer ThermalLinearElastic2DPlaneStrainNodal::Clone() const
{
    ThermalLinearElastic2DPlaneStrainNodal::Pointer p_clone(new ThermalLinearElastic2DPlaneStrainNodal(*this));
    return p_clone;
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void  ThermalLinearElastic2DPlaneStrainNodal::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,
        const double &rYoungModulus,
        const double &rPoissonCoefficient )
{
    rConstitutiveMatrix.clear();

    //plane strain constitutive matrix:
    rConstitutiveMatrix ( 0 , 0 ) = (rYoungModulus*(1.0-rPoissonCoefficient)/((1.0+rPoissonCoefficient)*(1.0-2.0*rPoissonCoefficient)));
    rConstitutiveMatrix ( 1 , 1 ) = rConstitutiveMatrix ( 0 , 0 );

    rConstitutiveMatrix ( 2 , 2 ) = rConstitutiveMatrix ( 0 , 0 )*(1.0-2.0*rPoissonCoefficient)/(2.0*(1.0-rPoissonCoefficient));

    rConstitutiveMatrix ( 0 , 1 ) = rConstitutiveMatrix ( 0 , 0 )*rPoissonCoefficient/(1.0-rPoissonCoefficient);
    rConstitutiveMatrix ( 1 , 0 ) = rConstitutiveMatrix ( 0 , 1 );
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic2DPlaneStrainNodal::GetLawFeatures(Features& rFeatures)
{
        //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
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

void ThermalLinearElastic2DPlaneStrainNodal::CalculateThermalStrain( Vector& rThermalStrainVector, const MaterialResponseVariables& rElasticVariables, double & rTemperature, double & rNodalReferenceTemperature)
{
    KRATOS_TRY


    //Identity vector
    rThermalStrainVector.resize(3,false);
    rThermalStrainVector[0] = 1.0;
    rThermalStrainVector[1] = 1.0;
    rThermalStrainVector[2] = 0.0;

    // Delta T
    double DeltaTemperature = rTemperature - rNodalReferenceTemperature;

    //Thermal strain vector // LameMu = (1 + poisson)
    for(unsigned int i = 0; i < 3; i++)
        rThermalStrainVector[i] *= rElasticVariables.LameMu * rElasticVariables.ThermalExpansionCoefficient * DeltaTemperature;

    KRATOS_CATCH( "" )
}

} // Namespace Kratos
