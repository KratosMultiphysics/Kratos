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
#include "custom_utilities/geometry_utilities.h"
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
    auto evaluate_shape_function_values = [&geometry = GetGeometry()](const auto& rIntegrationPoint) {
        auto result = Vector{};
        geometry.ShapeFunctionsValues(result, rIntegrationPoint);
        return result;
    };
    const auto& r_integration_points = mIntegrationScheme->GetIntegrationPoints();
    auto        shape_function_values_at_integration_points = std::vector<Vector>{};
    std::transform(r_integration_points.begin(), r_integration_points.end(),
                   std::back_inserter(shape_function_values_at_integration_points),
                   evaluate_shape_function_values);

    auto       b_matrices      = std::vector<Matrix>{};
    const auto dummy_gradients = Matrix{};
    for (auto i = std::size_t{0}; i < mIntegrationScheme->GetNumberOfIntegrationPoints(); ++i) {
        const Matrix b_matrix = prod(
            GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(GetGeometry(), r_integration_points[i]),
            mStressStatePolicy->CalculateBMatrix(
                dummy_gradients, shape_function_values_at_integration_points[i], GetGeometry()));

        b_matrices.emplace_back(b_matrix);
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
    auto index                    = std::size_t{0};
    for (const auto& integration_point : mIntegrationScheme->GetIntegrationPoints()) {
        integration_coefficients.push_back(mStressStatePolicy->CalculateIntegrationCoefficient(
            integration_point, determinants_of_jacobian[index], GetGeometry()));
    }

    rLeftHandSideMatrix = GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
        b_matrices, constitutive_matrices, integration_coefficients);
}

void LineInterfaceElement::CalculateRightHandSide(Element::VectorType& rRightHandSideVector,
                                                  const ProcessInfo&   rCurrentProcessInfo)
{
    rRightHandSideVector = ZeroVector{GetDofs().size()};
}

void LineInterfaceElement::CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                                        std::vector<ConstitutiveLaw::Pointer>& rOutput,
                                                        const ProcessInfo&)
{
    KRATOS_ERROR_IF_NOT(rVariable == CONSTITUTIVE_LAW)
        << "Cannot calculate on integration points: got unexpected variable " << rVariable.Name() << "\n";

    rOutput = mConstitutiveLaws;
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

void LineInterfaceElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    Element::Initialize(rCurrentProcessInfo);

    for (auto i = std::size_t{0}; i < mIntegrationScheme->GetNumberOfIntegrationPoints(); ++i) {
        mConstitutiveLaws.push_back(GetProperties()[CONSTITUTIVE_LAW]->Clone());
    }
}

Element::DofsVectorType LineInterfaceElement::GetDofs() const
{
    // At this point we only look at the U dofs, so we leave the water pressure nodes empty.
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(GetGeometry(), Geometry<Node>(),
                                                      GetGeometry().WorkingSpaceDimension());
}

} // namespace Kratos
