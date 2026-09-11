//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

/* Project includes */
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress_nodal.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/nodal_young_modulus_utilities.h"

namespace Kratos
{

//Default Constructor
ThermalLinearElastic2DPlaneStressNodal::ThermalLinearElastic2DPlaneStressNodal() : BaseType() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalLinearElastic2DPlaneStressNodal::ThermalLinearElastic2DPlaneStressNodal(const ThermalLinearElastic2DPlaneStressNodal& rOther) : BaseType(rOther) {}

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

void ThermalLinearElastic2DPlaneStressNodal::CalculateElasticMatrix(
    ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = NodalYoungModulusUtilities::InterpolatedYoungModulus(
        rValues.GetElementGeometry(), rValues.GetShapeFunctionsValues());
    const double NU = r_material_properties.GetValue(POISSON_RATIO,
        rValues.GetElementGeometry(), rValues.GetShapeFunctionsValues(), rValues.GetProcessInfo());
    ConstitutiveLawUtilities<3>::CalculateElasticMatrixPlaneStress(rConstitutiveMatrix, E, NU);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic2DPlaneStressNodal::CalculatePK2Stress(
    const ConstitutiveLaw::StrainVectorType& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = NodalYoungModulusUtilities::InterpolatedYoungModulus(
        rValues.GetElementGeometry(), rValues.GetShapeFunctionsValues());
    const double NU = r_material_properties.GetValue(POISSON_RATIO,
        rValues.GetElementGeometry(), rValues.GetShapeFunctionsValues(), rValues.GetProcessInfo());
    ConstitutiveLawUtilities<3>::CalculatePK2StressFromStrainPlaneStress(rStressVector, rStrainVector, E, NU);
}

} // Namespace Kratos