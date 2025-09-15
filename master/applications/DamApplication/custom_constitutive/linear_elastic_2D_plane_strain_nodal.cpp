//
//   Project Name:   
//   Last modified by:    $Author:     
//   Date:                $Date:     
//   Revision:            $Revision:     
//

/* Project includes */
#include "custom_constitutive/linear_elastic_2D_plane_strain_nodal.hpp"

namespace Kratos
{

//Default Constructor
LinearElastic2DPlaneStrainNodal::LinearElastic2DPlaneStrainNodal() : LinearElastic3DLawNodal() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
LinearElastic2DPlaneStrainNodal::LinearElastic2DPlaneStrainNodal(const LinearElastic2DPlaneStrainNodal& rOther) : LinearElastic3DLawNodal(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
LinearElastic2DPlaneStrainNodal::~LinearElastic2DPlaneStrainNodal() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer LinearElastic2DPlaneStrainNodal::Clone() const
{
    LinearElastic2DPlaneStrainNodal::Pointer p_clone(new LinearElastic2DPlaneStrainNodal(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LinearElastic2DPlaneStrainNodal::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{

    // e= 0.5*(1-invbT*invb)
    Matrix InverseLeftCauchyGreen ( 2 , 2 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector.clear();
    rStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = -InverseLeftCauchyGreen( 0, 1 );


}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void  LinearElastic2DPlaneStrainNodal::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,
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

void LinearElastic2DPlaneStrainNodal::GetLawFeatures(Features& rFeatures)
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

} // Namespace Kratos