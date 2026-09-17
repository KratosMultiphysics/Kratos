//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

/* Project includes */
#include "custom_constitutive/linear_elastic_3D_law_nodal.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/nodal_young_modulus_utilities.h"

namespace Kratos
{

//Default Constructor
LinearElastic3DLawNodal::LinearElastic3DLawNodal() : BaseType() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
LinearElastic3DLawNodal::LinearElastic3DLawNodal(const LinearElastic3DLawNodal& rOther) : BaseType(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
LinearElastic3DLawNodal::~LinearElastic3DLawNodal() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer LinearElastic3DLawNodal::Clone() const
{
    LinearElastic3DLawNodal::Pointer p_clone(new LinearElastic3DLawNodal(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LinearElastic3DLawNodal::CalculateElasticMatrix(
    ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = NodalYoungModulusUtilities::InterpolatedYoungModulus(
        rValues.GetElementGeometry(), rValues.GetShapeFunctionsValues());
    const double NU = r_material_properties.GetValue(POISSON_RATIO,
        rValues.GetElementGeometry(), rValues.GetShapeFunctionsValues(), rValues.GetProcessInfo());
    ConstitutiveLawUtilities<6>::CalculateElasticMatrix(rConstitutiveMatrix, E, NU);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LinearElastic3DLawNodal::CalculatePK2Stress(
    const ConstitutiveLaw::StrainVectorType& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = NodalYoungModulusUtilities::InterpolatedYoungModulus(
        rValues.GetElementGeometry(), rValues.GetShapeFunctionsValues());
    const double NU = r_material_properties.GetValue(POISSON_RATIO,
        rValues.GetElementGeometry(), rValues.GetShapeFunctionsValues(), rValues.GetProcessInfo());
    ConstitutiveLawUtilities<6>::CalculatePK2StressFromStrain(rStressVector, rStrainVector, E, NU);
}

} // Namespace Kratos