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
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/equation_of_motion_utilities.h"
#include "custom_utilities/geometry_utilities.h"
#include "interface_stress_state.h"
#include "lobatto_integration_scheme.h"

namespace
{

using namespace Kratos;

// The functions in this anonymous namespace are not specific to interface elements and
// will be moved to a utility module.
Vector CalculateDeterminantsOfJacobiansAtIntegrationPoints(const Geo::IntegrationPointVectorType& rIntegrationPoints,
                                                           const Geometry<Node>& rGeometry)
{
    auto result = Vector(rIntegrationPoints.size());
    std::transform(rIntegrationPoints.begin(), rIntegrationPoints.end(), result.begin(),
                   [&rGeometry](const auto& rIntegrationPoint) {
        return rGeometry.DeterminantOfJacobian(rIntegrationPoint);
    });

    return result;
}

std::vector<Matrix> CalculateConstitutiveMatricesAtIntegrationPoints(
    const std::vector<ConstitutiveLaw::Pointer>& rConstitutiveLaws, const Properties& rProperties)
{
    auto get_constitutive_matrix = [&rProperties](const auto& p_constitutive_law) {
        auto result         = Matrix{};
        auto law_parameters = ConstitutiveLaw::Parameters{};
        law_parameters.SetMaterialProperties(rProperties);
        p_constitutive_law->CalculateValue(law_parameters, CONSTITUTIVE_MATRIX, result);
        return result;
    };
    auto result = std::vector<Matrix>{};
    result.reserve(rConstitutiveLaws.size());
    std::transform(rConstitutiveLaws.begin(), rConstitutiveLaws.end(), std::back_inserter(result),
                   get_constitutive_matrix);

    return result;
}

} // namespace

namespace Kratos
{

LineInterfaceElement::LineInterfaceElement(IndexType NewId,
                                           const Geometry<GeometricalObject::NodeType>::Pointer& rGeometry,
                                           const Properties::Pointer& rProperties)
    : Element(NewId, rGeometry, rProperties),
      mIntegrationScheme(std::make_unique<LobattoIntegrationScheme>(GetGeometry().PointsNumber() / 2)),
      mStressStatePolicy(std::make_unique<InterfaceStressState>())
{
}

LineInterfaceElement::LineInterfaceElement(IndexType NewId, const GeometryType::Pointer& rGeometry)
    : Element(NewId, rGeometry),
      mIntegrationScheme(std::make_unique<LobattoIntegrationScheme>(GetGeometry().PointsNumber() / 2)),
      mStressStatePolicy(std::make_unique<InterfaceStressState>())
{
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
    return make_intrusive<LineInterfaceElement>(NewId, pGeometry, pProperties);
}

void LineInterfaceElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

void LineInterfaceElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo&)
{
    // Currently, the left-hand side matrix only includes the stiffness matrix. In the future, it
    // will also include water pressure contributions and coupling terms.
    rLeftHandSideMatrix = GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
        CalculateLocalBMatricesAtIntegrationPoints(),
        CalculateConstitutiveMatricesAtIntegrationPoints(), CalculateIntegrationCoefficients());
}

void LineInterfaceElement::CalculateRightHandSide(Element::VectorType& rRightHandSideVector, const ProcessInfo&)
{
    // Currently, the right-hand side only includes the internal force vector. In the future, it
    // will also include water pressure contributions and coupling terms.
    const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();
    const auto relative_displacements = CalculateRelativeDisplacementsAtIntegrationPoints(local_b_matrices);
    const auto tractions = CalculateTractionsAtIntegrationPoints(relative_displacements);
    const auto integration_coefficients = CalculateIntegrationCoefficients();
    rRightHandSideVector = -GeoEquationOfMotionUtilities::CalculateInternalForceVector(
        local_b_matrices, tractions, integration_coefficients);
}

void LineInterfaceElement::CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                                                VectorType&        rRightHandSideVector,
                                                const ProcessInfo& rCurrentProcessInfo)
{
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void LineInterfaceElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                        std::vector<Vector>&    rOutput,
                                                        const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == STRAIN) {
        const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();
        rOutput = CalculateRelativeDisplacementsAtIntegrationPoints(local_b_matrices);
    } else if (rVariable == CAUCHY_STRESS_VECTOR) {
        const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();
        const auto relative_displacements = CalculateRelativeDisplacementsAtIntegrationPoints(local_b_matrices);
        rOutput = CalculateTractionsAtIntegrationPoints(relative_displacements);
    } else {
        Element::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

void LineInterfaceElement::CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                                        std::vector<ConstitutiveLaw::Pointer>& rOutput,
                                                        const ProcessInfo&)
{
    KRATOS_ERROR_IF_NOT(rVariable == CONSTITUTIVE_LAW)
        << "Cannot calculate on integration points: got unexpected variable " << rVariable.Name() << "\n";

    rOutput = mConstitutiveLaws;
}

void LineInterfaceElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

void LineInterfaceElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    Element::Initialize(rCurrentProcessInfo);

    const auto shape_function_values_at_integration_points =
        GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(
            mIntegrationScheme->GetIntegrationPoints(), GetGeometry());

    mConstitutiveLaws.clear();
    for (const auto& r_shape_function_values : shape_function_values_at_integration_points) {
        mConstitutiveLaws.push_back(GetProperties()[CONSTITUTIVE_LAW]->Clone());
        mConstitutiveLaws.back()->InitializeMaterial(GetProperties(), GetGeometry(), r_shape_function_values);
    }
}

int LineInterfaceElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    int error = Element::Check(rCurrentProcessInfo);
    if (error != 0) return error;

    if (this->IsActive()) {
        KRATOS_ERROR_IF(mIntegrationScheme->GetNumberOfIntegrationPoints() != mConstitutiveLaws.size())
            << "Number of integration points (" << mIntegrationScheme->GetNumberOfIntegrationPoints()
            << ") and constitutive laws (" << mConstitutiveLaws.size() << ") do not match.\n";

        for (const auto& r_constitutive_law : mConstitutiveLaws) {
            error = r_constitutive_law->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);
            if (error != 0) return error;
        }
    }

    return 0;
}

Element::DofsVectorType LineInterfaceElement::GetDofs() const
{
    const auto no_Pw_geometry_yet = Geometry<Node>{};
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(GetGeometry(), no_Pw_geometry_yet,
                                                      GetGeometry().WorkingSpaceDimension());
}

std::vector<Matrix> LineInterfaceElement::CalculateLocalBMatricesAtIntegrationPoints() const
{
    const auto& r_integration_points = mIntegrationScheme->GetIntegrationPoints();
    const auto  shape_function_values_at_integration_points =
        GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(r_integration_points, GetGeometry());

    auto result = std::vector<Matrix>{};
    result.reserve(shape_function_values_at_integration_points.size());
    auto calculate_local_b_matrix = [&r_geometry = GetGeometry(), p_policy = mStressStatePolicy.get()](
                                        const auto& rShapeFunctionValuesAtIntegrationPoint,
                                        const auto& rIntegrationPoint) {
        // For interface elements, the shape function gradients are not used, since these are
        // non-continuum elements. Therefore, we pass an empty matrix.
        const auto dummy_gradients = Matrix{};
        return Matrix{prod(
            GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(r_geometry, rIntegrationPoint),
            p_policy->CalculateBMatrix(dummy_gradients, rShapeFunctionValuesAtIntegrationPoint, r_geometry))};
    };
    std::transform(shape_function_values_at_integration_points.begin(),
                   shape_function_values_at_integration_points.end(), r_integration_points.begin(),
                   std::back_inserter(result), calculate_local_b_matrix);

    return result;
}

std::vector<double> LineInterfaceElement::CalculateIntegrationCoefficients() const
{
    const auto determinants_of_jacobian = CalculateDeterminantsOfJacobiansAtIntegrationPoints(
        mIntegrationScheme->GetIntegrationPoints(), GetGeometry());
    return mIntegrationCoefficientsCalculator.Run<>(mIntegrationScheme->GetIntegrationPoints(),
                                                    determinants_of_jacobian, this);
}

std::vector<Matrix> LineInterfaceElement::CalculateConstitutiveMatricesAtIntegrationPoints()
{
    return ::CalculateConstitutiveMatricesAtIntegrationPoints(mConstitutiveLaws, GetProperties());
}

std::vector<Vector> LineInterfaceElement::CalculateRelativeDisplacementsAtIntegrationPoints(
    const std::vector<Matrix>& rLocalBMatrices) const
{
    const Geometry<Node> no_Pw_geometry;
    const auto           dofs = Geo::DofUtilities::ExtractUPwDofsFromNodes(
        GetGeometry(), no_Pw_geometry, GetGeometry().WorkingSpaceDimension());
    auto nodal_displacement_vector = Vector{dofs.size()};
    std::transform(dofs.begin(), dofs.end(), nodal_displacement_vector.begin(),
                   [](auto p_dof) { return p_dof->GetSolutionStepValue(); });

    auto result = std::vector<Vector>{};
    result.reserve(rLocalBMatrices.size());
    auto calculate_relative_displacement_vector = [&nodal_displacement_vector](const auto& rLocalB) {
        return Vector{prod(rLocalB, nodal_displacement_vector)};
    };
    std::transform(rLocalBMatrices.begin(), rLocalBMatrices.end(), std::back_inserter(result),
                   calculate_relative_displacement_vector);

    return result;
}

std::vector<Vector> LineInterfaceElement::CalculateTractionsAtIntegrationPoints(const std::vector<Vector>& rRelativeDisplacements)
{
    // We have to make a copy of each relative displacement vector, since setting it at the
    // constitutive law parameters requires a reference to a _mutable_ object!
    auto calculate_traction = [&properties = GetProperties()](auto RelativeDisplacement, auto& p_law) {
        auto law_parameters = ConstitutiveLaw::Parameters{};
        law_parameters.SetStrainVector(RelativeDisplacement);
        auto result = Vector{};
        law_parameters.SetStressVector(result);
        law_parameters.SetMaterialProperties(properties);
        p_law->CalculateMaterialResponseCauchy(law_parameters);
        return result;
    };
    auto result = std::vector<Vector>{};
    result.reserve(rRelativeDisplacements.size());
    std::transform(rRelativeDisplacements.begin(), rRelativeDisplacements.end(),
                   mConstitutiveLaws.begin(), std::back_inserter(result), calculate_traction);

    return result;
}

} // namespace Kratos
