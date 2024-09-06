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
//
#include "line_interface_element.h"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/equation_of_motion_utilities.h"
#include "interface_stress_state.h"
#include "lobatto_integration_scheme.h"

namespace Kratos
{

LineInterfaceElement::LineInterfaceElement() = default;

LineInterfaceElement::LineInterfaceElement(IndexType                                      NewId,
                                           Geometry<GeometricalObject::NodeType>::Pointer pGeometry,
                                           Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties),
      mIntegrationScheme(std::make_unique<LobattoIntegrationScheme>(GetGeometry().PointsNumber() / 2)),
      mStressStatePolicy(std::make_unique<InterfaceStressState>())
{
}

void LineInterfaceElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

void LineInterfaceElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix = ZeroMatrix{GetDofs().size(), GetDofs().size()};

    auto       shape_function_values_at_integration_points = std::vector<Vector>{};
    const auto number_of_node_pairs                        = GetGeometry().PointsNumber() / 2;
    for (const auto& integration_point : mIntegrationScheme->GetIntegrationPoints()) {
        auto shape_function_values = Vector{number_of_node_pairs};
        for (auto i = std::size_t{0}; i < number_of_node_pairs; ++i) {
            shape_function_values[i] = GetGeometry().ShapeFunctionValue(i, integration_point);
        }
        shape_function_values_at_integration_points.push_back(shape_function_values);
    }

    auto       b_matrices      = std::vector<Matrix>{};
    const auto dummy_gradients = Matrix{};
    for (auto i = std::size_t{0}; i < mIntegrationScheme->GetNumberOfIntegrationPoints(); ++i) {
        b_matrices.push_back(mStressStatePolicy->CalculateBMatrix(
            dummy_gradients, shape_function_values_at_integration_points[i], GetGeometry()));
    }

    const auto& properties         = GetProperties();
    auto        p_constitutive_law = properties.GetValue(CONSTITUTIVE_LAW);
    auto        law_parameters     = ConstitutiveLaw::Parameters{};
    law_parameters.SetMaterialProperties(GetProperties());
    auto constitutive_matrix = Matrix{};
    p_constitutive_law->CalculateValue(law_parameters, CONSTITUTIVE_MATRIX, constitutive_matrix);
    const auto constitutive_matrices =
        std::vector<Matrix>{mIntegrationScheme->GetNumberOfIntegrationPoints(), constitutive_matrix};

    auto determinants_of_jacobian = std::vector<double>{};
    for (const auto& integration_point : mIntegrationScheme->GetIntegrationPoints()) {
        determinants_of_jacobian.push_back(GetGeometry().DeterminantOfJacobian(integration_point));
    }

    auto integration_coefficients = std::vector<double>{};
    auto index = std::size_t{0};
    for (const auto& integration_point : mIntegrationScheme->GetIntegrationPoints()) {
        integration_coefficients.push_back(mStressStatePolicy->CalculateIntegrationCoefficient(integration_point, determinants_of_jacobian[index], GetGeometry()));
    }

    rLeftHandSideMatrix = GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
        b_matrices, constitutive_matrices, integration_coefficients);
}

Element::Pointer LineInterfaceElement::Create(IndexType               NewId,
                                              const NodesArrayType&   rNodes,
                                              PropertiesType::Pointer pProperties) const
{
    return Create(NewId, this->GetGeometry().Create(rNodes), pProperties);
}

Element::Pointer LineInterfaceElement::Create(IndexType               NewId,
                                              GeometryType::Pointer   pGeometry,
                                              PropertiesType::Pointer pProperties) const
{
    return {make_intrusive<LineInterfaceElement>(NewId, pGeometry, pProperties)};
}

void LineInterfaceElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    rElementalDofList = GetDofs();
}

Element::DofsVectorType LineInterfaceElement::GetDofs() const
{
    // At this point we only look at the U dofs, so we leave the water pressure nodes empty.
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(GetGeometry(), Geometry<Node>(),
                                                      GetGeometry().WorkingSpaceDimension());
}

} // namespace Kratos
