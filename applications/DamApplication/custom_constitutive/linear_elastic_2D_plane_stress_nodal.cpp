//
//   Project Name:   
//   Last modified by:    $Author:     
//   Date:                $Date:     
//   Revision:            $Revision:     
//

/* Project includes */
#include "custom_constitutive/linear_elastic_2D_plane_stress_nodal.hpp"

namespace Kratos
{

//Default Constructor
LinearElastic2DPlaneStressNodal::LinearElastic2DPlaneStressNodal() : LinearElastic2DPlaneStrainNodal() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
LinearElastic2DPlaneStressNodal::LinearElastic2DPlaneStressNodal(const LinearElastic2DPlaneStressNodal& rOther) : LinearElastic2DPlaneStrainNodal(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
LinearElastic2DPlaneStressNodal::~LinearElastic2DPlaneStressNodal() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer LinearElastic2DPlaneStressNodal::Clone() const
{
    LinearElastic2DPlaneStressNodal::Pointer p_clone(new LinearElastic2DPlaneStressNodal(*this));
    return p_clone;
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void  LinearElastic2DPlaneStressNodal::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,
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

void LinearElastic2DPlaneStressNodal::GetLawFeatures(Features& rFeatures)
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


} // Namespace Kratos
