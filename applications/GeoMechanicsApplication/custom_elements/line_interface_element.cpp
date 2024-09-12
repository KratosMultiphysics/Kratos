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
    rLeftHandSideMatrix = GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
        CalculateLocalBMatricesAtIntegrationPoints(),
        CalculateConstitutiveMatricesAtIntegrationPoints(), CalculateIntegrationCoefficients());
}

void LineInterfaceElement::CalculateRightHandSide(Element::VectorType& rRightHandSideVector,
                                                  const ProcessInfo&   rCurrentProcessInfo)
{
    const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();

    auto relative_displacements = CalculateRelativeDisplacements(local_b_matrices);

    auto tractions = std::vector<Vector>{};
    for (auto i = std::size_t{0}; i < mConstitutiveLaws.size(); ++i) {
        auto law_parameters = ConstitutiveLaw::Parameters{};
        law_parameters.SetStrainVector(relative_displacements[i]);
        auto traction = Vector{};
        law_parameters.SetStressVector(traction);
        law_parameters.SetMaterialProperties(GetProperties());
        mConstitutiveLaws[i]->CalculateMaterialResponseCauchy(law_parameters);
        tractions.push_back(traction);
    }

    const auto integration_coefficients = CalculateIntegrationCoefficients();

    rRightHandSideVector = -GeoEquationOfMotionUtilities::CalculateInternalForceVector(
        local_b_matrices, tractions, integration_coefficients);
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
    BaseType::Initialize(rCurrentProcessInfo);

    const auto unused_shape_function_values = Vector{};
    for (auto i = std::size_t{0}; i < mIntegrationScheme->GetNumberOfIntegrationPoints(); ++i) {
        mConstitutiveLaws.push_back(GetProperties()[CONSTITUTIVE_LAW]->Clone());
        mConstitutiveLaws.back()->InitializeMaterial(GetProperties(), GetGeometry(), unused_shape_function_values);
    }
}

Element::DofsVectorType LineInterfaceElement::GetDofs() const
{
    // At this point we only look at the U dofs, so we leave the water pressure nodes empty.
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(GetGeometry(), Geometry<Node>(),
                                                      GetGeometry().WorkingSpaceDimension());
}

std::vector<Matrix> LineInterfaceElement::CalculateLocalBMatricesAtIntegrationPoints() const
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

    auto       result          = std::vector<Matrix>{};
    const auto dummy_gradients = Matrix{};
    for (auto i = std::size_t{0}; i < mIntegrationScheme->GetNumberOfIntegrationPoints(); ++i) {
        const Matrix b_matrix = prod(
            GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(GetGeometry(), r_integration_points[i]),
            mStressStatePolicy->CalculateBMatrix(
                dummy_gradients, shape_function_values_at_integration_points[i], GetGeometry()));

        result.emplace_back(b_matrix);
    }

    return result;
}

std::vector<double> LineInterfaceElement::CalculateIntegrationCoefficients() const
{
    auto determinants_of_jacobian = std::vector<double>{};
    for (const auto& integration_point : mIntegrationScheme->GetIntegrationPoints()) {
        determinants_of_jacobian.push_back(GetGeometry().DeterminantOfJacobian(integration_point));
    }

    auto integration_coefficients = std::vector<double>{};
    auto index                    = size_t{0};
    for (const auto& integration_point : mIntegrationScheme->GetIntegrationPoints()) {
        integration_coefficients.push_back(mStressStatePolicy->CalculateIntegrationCoefficient(
            integration_point, determinants_of_jacobian[index], GetGeometry()));
        ++index;
    }

    return integration_coefficients;
}

std::vector<Matrix> LineInterfaceElement::CalculateConstitutiveMatricesAtIntegrationPoints()
{
    auto get_constitutive_matrix = [&properties = GetProperties()](const auto& p_law) {
        auto result         = Matrix{};
        auto law_parameters = ConstitutiveLaw::Parameters{};
        law_parameters.SetMaterialProperties(properties);
        p_law->CalculateValue(law_parameters, CONSTITUTIVE_MATRIX, result);
        return result;
    };
    auto constitutive_matrices = std::vector<Matrix>{};
    std::transform(mConstitutiveLaws.begin(), mConstitutiveLaws.end(),
                   std::back_inserter(constitutive_matrices), get_constitutive_matrix);
    return constitutive_matrices;
}

std::vector<Vector> LineInterfaceElement::CalculateRelativeDisplacements(const std::vector<Matrix>& rLocalBMatrices) const
{
    const auto dofs = Geo::DofUtilities::ExtractUPwDofsFromNodes(
        GetGeometry(), Geometry<Node>(), GetGeometry().WorkingSpaceDimension());
    auto nodal_displacements = Vector{dofs.size()};
    std::transform(dofs.begin(), dofs.end(), nodal_displacements.begin(),
                   [](auto p_dof) { return p_dof->GetSolutionStepValue(); });

    auto result = std::vector<Vector>{};
    for (const auto& r_b : rLocalBMatrices) {
        result.emplace_back(prod(r_b, nodal_displacements));
    }

    return result;
}

} // namespace Kratos
