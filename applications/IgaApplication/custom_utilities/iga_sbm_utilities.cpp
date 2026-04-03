//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// Project includes
#include "custom_utilities/iga_sbm_utilities.h"
#include "custom_conditions/gap_sbm_contact_condition.h"

#include <cmath>

namespace Kratos
{
namespace
{
class GapSbmContactConditionAccessor : public GapSbmContactCondition
{
public:
    GapSbmContactConditionAccessor(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : GapSbmContactCondition(NewId, pGeometry)
    {
    }

    void ComputeTaylorExpansionContributionWithDim(
        const GeometryType& rGeometry,
        const Vector& rDistanceVector,
        const SizeType BasisFunctionsOrder,
        const SizeType Dim,
        Vector& rHsumVector)
    {
        mDim = static_cast<unsigned int>(Dim);
        ComputeTaylorExpansionContribution(
            rGeometry,
            rDistanceVector,
            BasisFunctionsOrder,
            rHsumVector);
    }
};
} // namespace

void IgaSbmUtilities::GetDeformedPosition(
    const Condition& rCondition,
    array_1d<double, 3>& rPointDeformedCoordinates)
{
    KRATOS_ERROR_IF_NOT(rCondition.Has(NEIGHBOUR_GEOMETRIES))
        << "::[IgaSbmUtilities]:: Condition #" << rCondition.Id()
        << " missing NEIGHBOUR_GEOMETRIES." << std::endl;
    const auto& r_neighbour_geometries = rCondition.GetValue(NEIGHBOUR_GEOMETRIES);
    KRATOS_ERROR_IF(r_neighbour_geometries.empty())
        << "::[IgaSbmUtilities]:: Condition #" << rCondition.Id()
        << " has empty NEIGHBOUR_GEOMETRIES." << std::endl;

    const GeometryType& r_reference_geom = *r_neighbour_geometries[0];

    const array_1d<double, 3> skin_coordinates = rCondition.GetGeometry().Center().Coordinates();

    Vector distance_vector(3);
    noalias(distance_vector) = skin_coordinates - r_reference_geom.Center().Coordinates();

    const auto& r_DN_De = r_reference_geom.ShapeFunctionsLocalGradients(
        r_reference_geom.GetDefaultIntegrationMethod());
    const SizeType dim = r_DN_De[0].size2();
    KRATOS_ERROR_IF(dim != 2)
        << "::[IgaSbmUtilities]:: GetDeformedPosition currently supports only 2D. "
        << "Condition #" << rCondition.Id() << " has dimension " << dim << "." << std::endl;

    SizeType basis_functions_order = std::sqrt(r_DN_De[0].size1()) - 1;
    basis_functions_order *= 2;

    Vector H_sum_vec;
    GapSbmContactConditionAccessor accessor(0, rCondition.pGetGeometry());
    accessor.ComputeTaylorExpansionContributionWithDim(
        r_reference_geom,
        distance_vector,
        basis_functions_order,
        dim,
        H_sum_vec);

    Vector displacement_coefficient;
    GetSolutionCoefficientVector(DISPLACEMENT, r_reference_geom, displacement_coefficient);

    const SizeType number_of_cp = r_reference_geom.size();
    Matrix H_sum = ZeroMatrix(dim, dim * number_of_cp);
    for (IndexType i = 0; i < number_of_cp; ++i) {
        for (IndexType idim = 0; idim < dim; ++idim) {
            H_sum(idim, dim * i + idim) = H_sum_vec(i);
        }
    }

    const Vector displacement = prod(H_sum, displacement_coefficient);

    rPointDeformedCoordinates = skin_coordinates;
    rPointDeformedCoordinates[0] += displacement[0];
    rPointDeformedCoordinates[1] += displacement[1];
}

void IgaSbmUtilities::GetSolutionCoefficientVector(
    const Variable<array_1d<double, 3>>& rVariable,
    const GeometryType& rReferenceGeometry,
    Vector& rValues)
{
    const SizeType number_of_control_points = rReferenceGeometry.size();
    const SizeType mat_size = number_of_control_points * 2;

    if (rValues.size() != mat_size) {
        rValues.resize(mat_size, false);
    }

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        const array_1d<double, 3>& r_value = rReferenceGeometry[i].GetSolutionStepValue(rVariable);
        const IndexType index = i * 2;
        rValues[index] = r_value[0];
        rValues[index + 1] = r_value[1];
    }
}

} // namespace Kratos
