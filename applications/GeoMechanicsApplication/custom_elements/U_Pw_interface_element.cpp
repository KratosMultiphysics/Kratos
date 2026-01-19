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
//                   Gennady Markelov
//
#include "U_Pw_interface_element.h"

#include "custom_geometries/interface_geometry.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/dof_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/equation_of_motion_utilities.hpp"
#include "custom_utilities/extrapolation_utilities.h"
#include "custom_utilities/generic_utilities.hpp"
#include "custom_utilities/geometry_utilities.h"
#include "custom_utilities/math_utilities.hpp"
#include "geo_aliases.h"
#include "geometries/line_2d_2.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/triangle_3d_3.h"
#include "interface_stress_state.h"
#include "lobatto_integration_scheme.h"
#include "lumped_integration_scheme.h"

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

std::vector<Matrix> CalculateConstitutiveMatricesAtIntegrationPoints(const std::vector<ConstitutiveLaw::Pointer>& rConstitutiveLaws,
                                                                     const Properties& rProperties,
                                                                     const std::vector<Vector>& rRelativeDisplacements,
                                                                     const ProcessInfo& rProcessInfo)
{
    auto get_constitutive_matrix = [&rProperties, &rProcessInfo](const auto& p_constitutive_law,
                                                                 auto rRelativeDisplacement) {
        auto result = Matrix{p_constitutive_law->GetStrainSize(), p_constitutive_law->GetStrainSize()};
        auto law_parameters = ConstitutiveLaw::Parameters{};
        law_parameters.SetMaterialProperties(rProperties);
        law_parameters.SetStrainVector(rRelativeDisplacement);
        law_parameters.SetProcessInfo(rProcessInfo);
        p_constitutive_law->CalculateValue(law_parameters, CONSTITUTIVE_MATRIX, result);
        return result;
    };
    auto result = std::vector<Matrix>{};
    result.reserve(rConstitutiveLaws.size());
    std::transform(rConstitutiveLaws.begin(), rConstitutiveLaws.end(), rRelativeDisplacements.begin(),
                   std::back_inserter(result), get_constitutive_matrix);

    return result;
}

Geo::OptionalGeometryUniquePtr MakeOptionalWaterPressureGeometry(const Geometry<Node>& rDisplacementGeometry,
                                                                 IsDiffOrderElement IsDiffOrder)
{
    // Create a water pressure geometry only if it differs from the displacement geometry
    if (IsDiffOrder == IsDiffOrderElement::No) return std::nullopt;

    KRATOS_DEBUG_ERROR_IF(rDisplacementGeometry.GetGeometryOrderType() != GeometryData::Kratos_Quadratic_Order) << "Only quadratic order interface elements can create a linear order pressure geometry. \n";

    switch (rDisplacementGeometry.GetGeometryFamily()) {
        using enum GeometryData::KratosGeometryFamily;
    case Kratos_Linear: {
        const auto nodes = GeometryUtilities::GetNodesByIndex(rDisplacementGeometry, {0, 1, 3, 4});
        return std::make_unique<InterfaceGeometry<Line2D2<Node>>>(nodes);
    }
    case Kratos_Triangle: {
        const auto nodes = GeometryUtilities::GetNodesByIndex(rDisplacementGeometry, {0, 1, 2, 6, 7, 8});
        return std::make_unique<InterfaceGeometry<Triangle3D3<Node>>>(nodes);
    }
    case Kratos_Quadrilateral: {
        const auto nodes =
            GeometryUtilities::GetNodesByIndex(rDisplacementGeometry, {0, 1, 2, 3, 8, 9, 10, 11});
        return std::make_unique<InterfaceGeometry<Quadrilateral3D4<Node>>>(nodes);
    }
    default:
        break;
    }

    KRATOS_DEBUG_ERROR << "The specified geometry family is not supported for creating a water "
                          "pressure geometry.\n";
}

} // namespace

namespace Kratos
{

UPwInterfaceElement::UPwInterfaceElement(IndexType                          NewId,
                                         const GeometryType::Pointer&       rpGeometry,
                                         const Properties::Pointer&         rpProperties,
                                         std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                         IsDiffOrderElement                 IsDiffOrder)
    : Element(NewId, rpGeometry, rpProperties), mpStressStatePolicy(std::move(pStressStatePolicy))
{
    MakeIntegrationSchemeAndAssignFunction();
    mpOptionalPressureGeometry = MakeOptionalWaterPressureGeometry(GetDisplacementGeometry(), IsDiffOrder);
}

UPwInterfaceElement::UPwInterfaceElement(IndexType                          NewId,
                                         const GeometryType::Pointer&       rpGeometry,
                                         std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                         IsDiffOrderElement                 IsDiffOrder)
    : UPwInterfaceElement(NewId, rpGeometry, nullptr /* no properties */, std::move(pStressStatePolicy), IsDiffOrder)
{
}

void UPwInterfaceElement::MakeIntegrationSchemeAndAssignFunction()
{
    if (GetDisplacementGeometry().GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear) {
        mpIntegrationScheme =
            std::make_unique<LobattoIntegrationScheme>(GetDisplacementMidGeometry().PointsNumber());
        mfpCalculateRotationMatrix = GeometryUtilities::Calculate2DRotationMatrixForLineGeometry;
    } else {
        mpIntegrationScheme =
            std::make_unique<LumpedIntegrationScheme>(GetDisplacementMidGeometry().PointsNumber());
        mfpCalculateRotationMatrix = GeometryUtilities::Calculate3DRotationMatrixForPlaneGeometry;
    }
}

Element::Pointer UPwInterfaceElement::Create(IndexType               NewId,
                                             const NodesArrayType&   rNodes,
                                             PropertiesType::Pointer pProperties) const
{
    return Create(NewId, this->GetDisplacementGeometry().Create(rNodes), pProperties);
}

Element::Pointer UPwInterfaceElement::Create(IndexType               NewId,
                                             GeometryType::Pointer   pGeometry,
                                             PropertiesType::Pointer pProperties) const
{
    const auto is_diff_order = mpOptionalPressureGeometry ? IsDiffOrderElement::Yes : IsDiffOrderElement::No;
    return make_intrusive<UPwInterfaceElement>(NewId, pGeometry, pProperties,
                                               mpStressStatePolicy->Clone(), is_diff_order);
}

void UPwInterfaceElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

void UPwInterfaceElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rProcessInfo)
{
    // Currently, the left-hand side matrix only includes the stiffness matrix. In the future, it
    // will also include water pressure contributions and coupling terms.
    const auto number_of_dofs = GetDofs().size();
    rLeftHandSideMatrix       = ZeroMatrix{number_of_dofs, number_of_dofs};

    const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();
    const auto number_of_u_dofs =
        GetDisplacementGeometry().size() * GetDisplacementGeometry().WorkingSpaceDimension();
    subrange(rLeftHandSideMatrix, 0, number_of_u_dofs, 0, number_of_u_dofs) =
        GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
            local_b_matrices,
            CalculateConstitutiveMatricesAtIntegrationPoints(
                CalculateRelativeDisplacementsAtIntegrationPoints(local_b_matrices), rProcessInfo),
            CalculateIntegrationCoefficients());
}

void UPwInterfaceElement::CalculateRightHandSide(Element::VectorType& rRightHandSideVector,
                                                 const ProcessInfo&   rProcessInfo)
{
    rRightHandSideVector = ZeroVector{GetDofs().size()};

    // Currently, the right-hand side only includes the internal force vector. In the future, it
    // will also include water pressure contributions and coupling terms.
    const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();
    const auto relative_displacements = CalculateRelativeDisplacementsAtIntegrationPoints(local_b_matrices);
    const auto tractions = CalculateTractionsAtIntegrationPoints(relative_displacements, rProcessInfo);
    const auto integration_coefficients = CalculateIntegrationCoefficients();
    const auto number_of_u_dofs =
        GetDisplacementGeometry().size() * GetDisplacementGeometry().WorkingSpaceDimension();
    subrange(rRightHandSideVector, 0, 0 + number_of_u_dofs) =
        -1.0 * GeoEquationOfMotionUtilities::CalculateInternalForceVector(
                   local_b_matrices, tractions, integration_coefficients);
}

void UPwInterfaceElement::CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                                               VectorType&        rRightHandSideVector,
                                               const ProcessInfo& rCurrentProcessInfo)
{
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void UPwInterfaceElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                       std::vector<Vector>&    rOutput,
                                                       const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == GEO_RELATIVE_DISPLACEMENT_VECTOR) {
        const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();
        rOutput = CalculateRelativeDisplacementsAtIntegrationPoints(local_b_matrices);
    } else if (rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
        const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();
        const auto relative_displacements = CalculateRelativeDisplacementsAtIntegrationPoints(local_b_matrices);
        rOutput = CalculateTractionsAtIntegrationPoints(relative_displacements, rCurrentProcessInfo);
    } else {
        Element::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

void UPwInterfaceElement::CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                                       std::vector<ConstitutiveLaw::Pointer>& rOutput,
                                                       const ProcessInfo&)
{
    KRATOS_ERROR_IF_NOT(rVariable == CONSTITUTIVE_LAW)
        << "Cannot calculate on integration points: got unexpected variable " << rVariable.Name() << "\n";

    rOutput = mConstitutiveLaws;
}

void UPwInterfaceElement::Calculate(const Variable<Vector>& rVariable, Vector& rOutput, const ProcessInfo& rProcessInfo)
{
    KRATOS_ERROR_IF_NOT(rVariable == INTERNAL_FORCES_VECTOR || rVariable == EXTERNAL_FORCES_VECTOR)
        << "Variable " << rVariable.Name() << " is unknown for element with Id " << this->GetId() << ".";

    rOutput = ZeroVector{GetDofs().size()};

    // Currently, the right-hand side only includes the internal force vector. In the future, it
    // will also include water pressure contributions and coupling terms.
    const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();
    const auto relative_displacements = CalculateRelativeDisplacementsAtIntegrationPoints(local_b_matrices);
    const auto tractions = CalculateTractionsAtIntegrationPoints(relative_displacements, rProcessInfo);
    const auto integration_coefficients = CalculateIntegrationCoefficients();
    const auto number_of_u_dofs =
        GetDisplacementGeometry().size() * GetDisplacementGeometry().WorkingSpaceDimension();
    if (rVariable == INTERNAL_FORCES_VECTOR) {
        subrange(rOutput, 0, 0 + number_of_u_dofs) = GeoEquationOfMotionUtilities::CalculateInternalForceVector(
            local_b_matrices, tractions, integration_coefficients);
    } else if (rVariable == EXTERNAL_FORCES_VECTOR) {
        subrange(rOutput, 0, 0 + number_of_u_dofs) = Vector{ZeroVector{number_of_u_dofs}};
    }
}

void UPwInterfaceElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

void UPwInterfaceElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    Element::Initialize(rCurrentProcessInfo);

    mConstitutiveLaws.clear();
    for (auto i = std::size_t{0}; i < mpIntegrationScheme->GetNumberOfIntegrationPoints(); ++i) {
        mConstitutiveLaws.push_back(GetProperties()[CONSTITUTIVE_LAW]->Clone());
    }
    // Only interpolate when neighbouring elements that provide nodal stresses were found
    if (this->Has(NEIGHBOUR_ELEMENTS) && this->GetValue(NEIGHBOUR_ELEMENTS).size() > 0) {
        // interface element can at maximum have 2 neighbours
        KRATOS_DEBUG_ERROR_IF(this->GetValue(NEIGHBOUR_ELEMENTS).size() > 2)
            << "Too many neighbour elements for interface element " << this->Id() << std::endl;
        const auto interface_node_ids =
            GenericUtilities::GetIdsFromEntityContents(GetDisplacementGeometry());
        std::vector<std::optional<Vector>> interface_nodal_cauchy_stresses(interface_node_ids.size());
        auto&               r_neighbour_element = this->GetValue(NEIGHBOUR_ELEMENTS).front();
        std::vector<Vector> neighbour_cauchy_stresses;
        // Note that the interface elements don't account for water pressures yet. Consequently,
        // we need to consider the total stresses rather than the effective stresses to calculate
        // the appropriate prestresses to be applied.
        r_neighbour_element.CalculateOnIntegrationPoints(
            TOTAL_STRESS_VECTOR, neighbour_cauchy_stresses, rCurrentProcessInfo);
        interface_nodal_cauchy_stresses = ExtrapolationUtilities::CalculateNodalVectors(
            interface_node_ids, r_neighbour_element, neighbour_cauchy_stresses);
        InterpolateNodalStressesToInitialTractions(interface_nodal_cauchy_stresses);
    }
    const auto shape_function_values_at_integration_points =
        GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(
            mpIntegrationScheme->GetIntegrationPoints(), GetDisplacementGeometry());
    for (auto i = std::size_t{0}; i < mConstitutiveLaws.size(); ++i) {
        mConstitutiveLaws[i]->InitializeMaterial(GetProperties(), GetDisplacementGeometry(),
                                                 shape_function_values_at_integration_points[i]);
    }
}

int UPwInterfaceElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    int error = Element::Check(rCurrentProcessInfo);
    if (error != 0) return error;

    if (this->IsActive()) {
        KRATOS_ERROR_IF(mpIntegrationScheme->GetNumberOfIntegrationPoints() != mConstitutiveLaws.size())
            << "Number of integration points (" << mpIntegrationScheme->GetNumberOfIntegrationPoints()
            << ") and constitutive laws (" << mConstitutiveLaws.size() << ") do not match.\n";

        const auto r_properties  = GetProperties();
        const auto expected_size = mpStressStatePolicy->GetVoigtSize();
        ConstitutiveLawUtilities::CheckStrainSize(r_properties, expected_size, Id());

        error = r_properties[CONSTITUTIVE_LAW]->Check(r_properties, GetDisplacementGeometry(), rCurrentProcessInfo);
        return error;
    }

    return 0;
}

const IntegrationScheme& UPwInterfaceElement::GetIntegrationScheme() const
{
    return *mpIntegrationScheme;
}

const Geometry<Node>& UPwInterfaceElement::GetDisplacementMidGeometry() const
{
    constexpr auto unused_part_index = std::size_t{0};
    return GetDisplacementGeometry().GetGeometryPart(unused_part_index);
}

Element::DofsVectorType UPwInterfaceElement::GetDofs() const
{
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(GetDisplacementGeometry(), GetWaterPressureGeometry(),
                                                      GetDisplacementGeometry().WorkingSpaceDimension());
}

const Element::GeometryType& UPwInterfaceElement::GetDisplacementGeometry() const
{
    return GetGeometry();
}

const Element::GeometryType& UPwInterfaceElement::GetWaterPressureGeometry() const
{
    return mpOptionalPressureGeometry ? *(mpOptionalPressureGeometry.value()) : GetDisplacementGeometry();
}

std::vector<Matrix> UPwInterfaceElement::CalculateLocalBMatricesAtIntegrationPoints() const
{
    const auto& r_integration_points = mpIntegrationScheme->GetIntegrationPoints();
    const auto  shape_function_values_at_integration_points =
        GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(r_integration_points,
                                                                       GetDisplacementGeometry());

    auto result = std::vector<Matrix>{};
    result.reserve(shape_function_values_at_integration_points.size());
    auto calculate_local_b_matrix = [&r_geometry = GetDisplacementGeometry(), this](
                                        const auto& rShapeFunctionValuesAtIntegrationPoint,
                                        const auto& rIntegrationPoint) {
        // For interface elements, the shape function gradients are not used, since these are
        // non-continuum elements. Therefore, we pass an empty matrix.
        const auto dummy_gradients = Matrix{};
        auto       b_matrix        = mpStressStatePolicy->CalculateBMatrix(
            dummy_gradients, rShapeFunctionValuesAtIntegrationPoint, r_geometry);
        ApplyRotationToBMatrix(b_matrix, mfpCalculateRotationMatrix(r_geometry, rIntegrationPoint));
        return b_matrix;
    };
    std::transform(shape_function_values_at_integration_points.begin(),
                   shape_function_values_at_integration_points.end(), r_integration_points.begin(),
                   std::back_inserter(result), calculate_local_b_matrix);

    return result;
}

std::vector<double> UPwInterfaceElement::CalculateIntegrationCoefficients() const
{
    const auto determinants_of_jacobian = CalculateDeterminantsOfJacobiansAtIntegrationPoints(
        mpIntegrationScheme->GetIntegrationPoints(), GetDisplacementGeometry());
    return mIntegrationCoefficientsCalculator.Run<>(mpIntegrationScheme->GetIntegrationPoints(),
                                                    determinants_of_jacobian, this);
}

std::vector<Matrix> UPwInterfaceElement::CalculateConstitutiveMatricesAtIntegrationPoints(
    const std::vector<Vector>& rRelativeDisplacements, const ProcessInfo& rProcessInfo)
{
    return ::CalculateConstitutiveMatricesAtIntegrationPoints(mConstitutiveLaws, GetProperties(),
                                                              rRelativeDisplacements, rProcessInfo);
}

std::vector<Vector> UPwInterfaceElement::CalculateRelativeDisplacementsAtIntegrationPoints(
    const std::vector<Matrix>& rLocalBMatrices) const
{
    // We feel that we'd better have two extra member functions: one that returns the displacement
    // degrees of freedom, and another one that returns the water pressure degrees of freedom. That
    // would make this function body easier to understand.
    const Geometry<Node> no_Pw_geometry;
    const auto           dofs = Geo::DofUtilities::ExtractUPwDofsFromNodes(
        GetDisplacementGeometry(), no_Pw_geometry, GetDisplacementGeometry().WorkingSpaceDimension());
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

std::vector<Vector> UPwInterfaceElement::CalculateTractionsAtIntegrationPoints(const std::vector<Vector>& rRelativeDisplacements,
                                                                               const ProcessInfo& rProcessInfo)
{
    // We have to make a copy of each relative displacement vector, since setting it at the
    // constitutive law parameters requires a reference to a _mutable_ object!
    auto calculate_traction = [&properties = GetProperties(), &rProcessInfo](
                                  auto RelativeDisplacement, auto& p_law) {
        auto law_parameters = ConstitutiveLaw::Parameters{};
        law_parameters.SetStrainVector(RelativeDisplacement);
        auto result = Vector{};
        result.resize(p_law->GetStrainSize());
        law_parameters.SetStressVector(result);
        law_parameters.SetMaterialProperties(properties);
        law_parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
        law_parameters.SetProcessInfo(rProcessInfo);
        p_law->CalculateMaterialResponseCauchy(law_parameters);
        return result;
    };
    auto result = std::vector<Vector>{};
    result.reserve(rRelativeDisplacements.size());
    std::ranges::transform(rRelativeDisplacements, mConstitutiveLaws, std::back_inserter(result), calculate_traction);

    return result;
}

void UPwInterfaceElement::ApplyRotationToBMatrix(Matrix& rBMatrix, const Matrix& rRotationMatrix) const
{
    const auto dim = GetDisplacementGeometry().WorkingSpaceDimension();
    for (auto i = std::size_t{0}; i + dim <= rBMatrix.size2(); i += dim) {
        auto sub_matrix = subrange(rBMatrix, 0, rBMatrix.size1(), i, i + dim);
        sub_matrix.assign(Matrix{prod(sub_matrix, trans(rRotationMatrix))});
    }
}

void UPwInterfaceElement::InterpolateNodalStressesToInitialTractions(
    const std::vector<std::optional<Vector>>& rInterfaceNodalCauchyStresses) const
{
    // interpolate nodal stresses on a chosen side
    const auto number_of_nodes_on_side = GetDisplacementGeometry().PointsNumber() / 2;
    // check which side is shared with the neighbour geometry, by checking for the first nodes of the first side
    auto in_first_side = rInterfaceNodalCauchyStresses[0].has_value();

    const std::size_t   start_index = in_first_side ? 0 : number_of_nodes_on_side;
    std::vector<Vector> nodal_stresses;
    for (auto i = std::size_t{0}; i < number_of_nodes_on_side; ++i) {
        nodal_stresses.push_back(rInterfaceNodalCauchyStresses[start_index + i].value());
    }

    std::size_t integration_point_index = 0;
    for (const auto& r_integration_point : mpIntegrationScheme->GetIntegrationPoints()) {
        const auto integration_point_stress =
            InterpolateNodalStressToIntegrationPoints(r_integration_point, nodal_stresses);

        const auto integration_point_local_stress_tensor =
            RotateStressToLocalCoordinates(r_integration_point, integration_point_stress);
        const auto traction_vector = ConvertLocalStressToTraction(integration_point_local_stress_tensor);

        const auto initial_state =
            make_intrusive<InitialState>(traction_vector, InitialState::InitialImposingType::STRESS_ONLY);
        mConstitutiveLaws[integration_point_index]->SetInitialState(initial_state);
        ++integration_point_index;
    }
}

Vector UPwInterfaceElement::InterpolateNodalStressToIntegrationPoints(const Geo::IntegrationPointType& rIntegrationPoint,
                                                                      const std::vector<Vector>& rNodalStresses) const
{
    auto result                                  = Vector(rNodalStresses[0].size(), 0.0);
    auto integration_point_shape_function_values = Vector{};
    GetDisplacementMidGeometry().ShapeFunctionsValues(integration_point_shape_function_values, rIntegrationPoint);
    for (auto i = std::size_t{0}; i < GetDisplacementMidGeometry().PointsNumber(); ++i) {
        result += integration_point_shape_function_values[i] * rNodalStresses[i];
    }
    return result;
}

Matrix UPwInterfaceElement::RotateStressToLocalCoordinates(const Geo::IntegrationPointType& rIntegrationPoint,
                                                           const Vector& rGlobalStressVector) const
{
    auto rotation_tensor = Matrix(3, 3, 0.0);
    if (GetDisplacementGeometry().GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear) {
        const auto two_d_rotation_tensor =
            mfpCalculateRotationMatrix(GetDisplacementGeometry(), rIntegrationPoint);
        GeoElementUtilities::AddMatrixAtPosition(two_d_rotation_tensor, rotation_tensor, 0, 0);
        rotation_tensor(2, 2) = 1.0;
    } else {
        rotation_tensor = mfpCalculateRotationMatrix(GetDisplacementGeometry(), rIntegrationPoint);
    }
    return GeoMechanicsMathUtilities::RotateSecondOrderTensor(
        MathUtils<>::StressVectorToTensor(rGlobalStressVector), rotation_tensor);
}

Vector UPwInterfaceElement::ConvertLocalStressToTraction(const Matrix& rLocalStress) const
{
    // extract normal and shear components to form initial traction
    auto result = Vector(mConstitutiveLaws[0]->GetStrainSize());
    if (GetDisplacementGeometry().GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear) {
        result[0] = rLocalStress(1, 1);
        result[1] = rLocalStress(0, 1);
    } else {
        result[0] = rLocalStress(2, 2);
        result[1] = rLocalStress(0, 2);
        result[2] = rLocalStress(1, 2);
    }
    return result;
}

// Instances of this class can not be copied but can be moved. Check that at compile time.
static_assert(!std::is_copy_constructible_v<UPwInterfaceElement>);
static_assert(!std::is_copy_assignable_v<UPwInterfaceElement>);
static_assert(std::is_move_constructible_v<UPwInterfaceElement>);
static_assert(std::is_move_assignable_v<UPwInterfaceElement>);
} // namespace Kratos
