//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

/* Project includes */
#include "custom_constitutive/linear_elastic_2D_plane_strain_nodal.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/nodal_young_modulus_utilities.h"

namespace Kratos
{

//Default Constructor
LinearElastic2DPlaneStrainNodal::LinearElastic2DPlaneStrainNodal() : BaseType() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
LinearElastic2DPlaneStrainNodal::LinearElastic2DPlaneStrainNodal(const LinearElastic2DPlaneStrainNodal& rOther) : BaseType(rOther) {}

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

void LinearElastic2DPlaneStrainNodal::CalculateElasticMatrix(
    ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = NodalYoungModulusUtilities::InterpolatedYoungModulus(
        rValues.GetElementGeometry(), rValues.GetShapeFunctionsValues());
    const double NU = r_material_properties.GetValue(POISSON_RATIO,
        rValues.GetElementGeometry(), rValues.GetShapeFunctionsValues(), rValues.GetProcessInfo());
    ConstitutiveLawUtilities<3>::CalculateElasticMatrixPlaneStrain(rConstitutiveMatrix, E, NU);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LinearElastic2DPlaneStrainNodal::CalculatePK2Stress(
    const ConstitutiveLaw::StrainVectorType& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = NodalYoungModulusUtilities::InterpolatedYoungModulus(
        rValues.GetElementGeometry(), rValues.GetShapeFunctionsValues());
    const double NU = r_material_properties.GetValue(POISSON_RATIO,
        rValues.GetElementGeometry(), rValues.GetShapeFunctionsValues(), rValues.GetProcessInfo());
    ConstitutiveLawUtilities<3>::CalculatePK2StressFromStrainPlaneStrain(rStressVector, rStrainVector, E, NU);
}

} // Namespace Kratos