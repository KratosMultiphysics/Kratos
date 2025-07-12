// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Anne van de Graaf
//
#include "interface_stress_state.h"
#include "geo_mechanics_application_constants.h"
#include "includes/serializer.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <type_traits>

namespace
{

using namespace Kratos;

Matrix CalculateBMatrix(const Vector& rN, const Geometry<Node>& rGeometry, const std::vector<std::size_t>& rComponentOrder)
{
    KRATOS_ERROR_IF(rN.empty())
        << "Shape function values are empty. Therefore, the B matrix can not be computed.\n";
    KRATOS_ERROR_IF_NOT(rN.size() == rGeometry.size() / 2)
        << "The number of shape functions should be equal to the number of node pairs. Therefore, "
           "the B matrix can not be computed.\n";

    auto result =
        Matrix{ZeroMatrix{rComponentOrder.size(), rGeometry.WorkingSpaceDimension() * rGeometry.size()}};

    const auto number_of_u_dofs_per_side = result.size2() / 2;
    for (unsigned int i = 0; i < rGeometry.size() / 2; ++i) {
        // Use the order in which the degrees of freedom at any node must be processed to compute
        // the normal component(s) first and then the tangential component(s)
        for (unsigned int j = 0; j < rComponentOrder.size(); ++j) {
            result(j, i * rGeometry.WorkingSpaceDimension() + rComponentOrder[j]) = -rN[i];
            result(j, i * rGeometry.WorkingSpaceDimension() + rComponentOrder[j] + number_of_u_dofs_per_side) =
                rN[i];
        }
    }

    return result;
}

} // namespace

namespace Kratos
{

Matrix Line2DInterfaceStressState::CalculateBMatrix(const Matrix&, const Vector& rN, const Geometry<Node>& rGeometry) const
{
    // Define the order in which the degrees of freedom at any node must be processed to compute the
    // normal component first and then the tangential component
    const auto component_order = std::vector<std::size_t>{1, 0};
    return ::CalculateBMatrix(rN, rGeometry, component_order);
}

Vector Line2DInterfaceStressState::CalculateGreenLagrangeStrain(const Matrix&) const
{
    KRATOS_ERROR << "For line interfaces, it is not possible to calculate the Green-Lagrange "
                    "strain.\n";
}

std::unique_ptr<StressStatePolicy> Line2DInterfaceStressState::Clone() const
{
    return std::make_unique<Line2DInterfaceStressState>();
}

const Vector& Line2DInterfaceStressState::GetVoigtVector() const { return VoigtVectorInterface2D; }

SizeType Line2DInterfaceStressState::GetVoigtSize() const { return VOIGT_SIZE_2D_INTERFACE; }

SizeType Line2DInterfaceStressState::GetStressTensorSize() const
{
    KRATOS_ERROR << "For line interfaces, the stress tensor size is not implemented.\n";
}

Vector Line2DInterfaceStressState::DefineInterfaceVoigtVector()
{
    Vector result{VOIGT_SIZE_2D_INTERFACE};
    result <<= 1.0, 0.0;

    return result;
}

const Vector Line2DInterfaceStressState::VoigtVectorInterface2D =
    Line2DInterfaceStressState::DefineInterfaceVoigtVector();

void Line2DInterfaceStressState::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void Line2DInterfaceStressState::load(Serializer&)
{
    // No data members to be loaded (yet)
}

static_assert(!std::is_copy_constructible_v<Line2DInterfaceStressState>);
static_assert(!std::is_copy_assignable_v<Line2DInterfaceStressState>);
static_assert(std::is_move_constructible_v<Line2DInterfaceStressState>);
static_assert(std::is_move_assignable_v<Line2DInterfaceStressState>);

Matrix SurfaceInterfaceStressState::CalculateBMatrix(const Matrix&, const Vector& rN, const Geometry<Node>& rGeometry) const
{
    // Define the order in which the degrees of freedom at any node must be processed to compute the
    // normal component first and then the two tangential components
    const auto component_order = std::vector<std::size_t>{2, 0, 1};
    return ::CalculateBMatrix(rN, rGeometry, component_order);
}

Vector SurfaceInterfaceStressState::CalculateGreenLagrangeStrain(const Matrix&) const
{
    KRATOS_ERROR << "For surface interfaces, it is not possible to calculate the Green-Lagrange "
                    "strain.\n";
}

std::unique_ptr<StressStatePolicy> SurfaceInterfaceStressState::Clone() const
{
    return std::make_unique<SurfaceInterfaceStressState>();
}

const Vector& SurfaceInterfaceStressState::GetVoigtVector() const { return VoigtVectorInterface3D; }

SizeType SurfaceInterfaceStressState::GetVoigtSize() const { return VOIGT_SIZE_3D_INTERFACE; }

SizeType SurfaceInterfaceStressState::GetStressTensorSize() const
{
    KRATOS_ERROR << "For surface interfaces, the stress tensor size is not implemented.\n";
}

Vector SurfaceInterfaceStressState::DefineInterfaceVoigtVector()
{
    Vector result{VOIGT_SIZE_3D_INTERFACE};
    result <<= 1.0, 0.0, 0.0;

    return result;
}

const Vector SurfaceInterfaceStressState::VoigtVectorInterface3D =
    SurfaceInterfaceStressState::DefineInterfaceVoigtVector();

void SurfaceInterfaceStressState::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void SurfaceInterfaceStressState::load(Serializer&)
{
    // No data members to be loaded (yet)
}

static_assert(!std::is_copy_constructible_v<SurfaceInterfaceStressState>);
static_assert(!std::is_copy_assignable_v<SurfaceInterfaceStressState>);
static_assert(std::is_move_constructible_v<SurfaceInterfaceStressState>);
static_assert(std::is_move_assignable_v<SurfaceInterfaceStressState>);

} // namespace Kratos
