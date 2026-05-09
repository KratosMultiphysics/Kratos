// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

// Project includes
#include "custom_elements/geo_linear_timoshenko_beam_element_2D2N.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "includes/constitutive_law.h"

namespace Kratos
{
void GeoLinearTimoshenkoBeamElement2D2N::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (mConstitutiveLawVector.empty()) return;

    const auto& r_geometry = GetGeometry();
    const auto& r_props    = GetProperties();
    const auto  length     = CalculateLength();
    const auto  Phi        = StructuralMechanicsElementUtilities::CalculatePhi(r_props, length);
    const auto& integration_points = IntegrationPoints(GetIntegrationMethod());
    const auto  strain_size        = mConstitutiveLawVector[0]->GetStrainSize();

    VectorType nodal_values(GetDoFsPerNode() * r_geometry.size());
    GetNodalValuesVector(nodal_values);

    VectorType strain_vector(strain_size);
    VectorType stress_vector(strain_size);
    strain_vector.clear();
    stress_vector.clear();

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rCurrentProcessInfo);
    cl_values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);

    for (IndexType ip = 0; ip < integration_points.size(); ++ip) {
        if (!mConstitutiveLawVector[ip]->RequiresFinalizeMaterialResponse()) continue;
        CalculateGeneralizedStrainsVector(strain_vector, length, Phi, integration_points[ip].X(), nodal_values);
        mConstitutiveLawVector[ip]->FinalizeMaterialResponsePK2(cl_values);
    }

    KRATOS_CATCH("")
}

} // namespace Kratos
