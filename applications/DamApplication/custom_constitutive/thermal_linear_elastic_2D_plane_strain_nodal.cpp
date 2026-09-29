//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

/* Project includes */
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain_nodal.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/nodal_young_modulus_utilities.h"

namespace Kratos
{

//Default Constructor
ThermalLinearElastic2DPlaneStrainNodal::ThermalLinearElastic2DPlaneStrainNodal() : BaseType() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalLinearElastic2DPlaneStrainNodal::ThermalLinearElastic2DPlaneStrainNodal(const ThermalLinearElastic2DPlaneStrainNodal& rOther) : BaseType(rOther) {}

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

void ThermalLinearElastic2DPlaneStrainNodal::CalculateElasticMatrix(
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

void ThermalLinearElastic2DPlaneStrainNodal::CalculatePK2Stress(
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