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

#include "contribution_calculators/calculation_contribution.h"
#include "contribution_calculators/fluid_body_flow_calculator.hpp"
#include "contribution_calculators/permeability_calculator.hpp"
#include "contribution_calculators/stiffness_calculator.hpp"
#include "custom_geometries/interface_geometry.hpp"
#include "custom_retention/retention_law_factory.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/dof_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/equation_of_motion_utilities.hpp"
#include "custom_utilities/extrapolation_utilities.h"
#include "custom_utilities/generic_utilities.hpp"
#include "custom_utilities/geometry_utilities.h"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
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
    std::ranges::transform(rIntegrationPoints, result.begin(), [&rGeometry](const auto& rIntegrationPoint) {
        return rGeometry.DeterminantOfJacobian(rIntegrationPoint);
    });

    return result;
}

Geo::GeometryUniquePtr MakeOptionalWaterPressureGeometry(const Geometry<Node>& rDisplacementGeometry,
                                                         IsDiffOrderElement IsDiffOrder)
{
    // Create a water pressure geometry only if it differs from the displacement geometry
    if (IsDiffOrder == IsDiffOrderElement::No) return nullptr;

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
    return nullptr; // required for release builds
}

Geo::ProcessInfoGetter CreateProcessInfoGetter(const ProcessInfo& rProcessInfo)
{
    return [&rProcessInfo]() -> const ProcessInfo& { return rProcessInfo; };
}

} // namespace

namespace Kratos
{

UPwInterfaceElement::UPwInterfaceElement(IndexType                          NewId,
                                         const GeometryType::Pointer&       rpGeometry,
                                         const Properties::Pointer&         rpProperties,
                                         std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                         IsDiffOrderElement                 IsDiffOrder,
                                         const std::vector<CalculationContribution>& rContributions)
    : Element(NewId, rpGeometry, rpProperties),
      mpStressStatePolicy(std::move(pStressStatePolicy)),
      mContributions(rContributions)
{
    MakeIntegrationSchemeAndAssignFunction();
    mpOptionalPressureGeometry = MakeOptionalWaterPressureGeometry(GetDisplacementGeometry(), IsDiffOrder);
}

UPwInterfaceElement::UPwInterfaceElement(IndexType                          NewId,
                                         const GeometryType::Pointer&       rpGeometry,
                                         std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                         IsDiffOrderElement                 IsDiffOrder,
                                         const std::vector<CalculationContribution>& rContributions)
    : UPwInterfaceElement(
          NewId, rpGeometry, nullptr /* no properties */, std::move(pStressStatePolicy), IsDiffOrder, rContributions)
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
    return make_intrusive<UPwInterfaceElement>(
        NewId, pGeometry, pProperties, mpStressStatePolicy->Clone(), is_diff_order, mContributions);
}

void UPwInterfaceElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

void UPwInterfaceElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rProcessInfo)
{
    const auto number_of_dofs = GetDofs().size();
    rLeftHandSideMatrix       = ZeroMatrix{number_of_dofs, number_of_dofs};

    // Currently, the left-hand side matrix only includes the stiffness matrix. In the future, it
    // will also include water pressure contributions and coupling terms.
    for (auto contribution : mContributions) {
        switch (contribution) {
        case CalculationContribution::Stiffness:
            CalculateAndAssignStifnessMatrix(rLeftHandSideMatrix, rProcessInfo);
            break;
        case CalculationContribution::Permeability:
            CalculateAndAssignPermeabilityMatrix(rLeftHandSideMatrix);
            break;
        case CalculationContribution::FluidBodyFlow:
            break;
        default:
            KRATOS_ERROR << "This contribution is not supported \n";
        }
    }
}

void UPwInterfaceElement::CalculateAndAssignStifnessMatrix(Element::MatrixType& rLeftHandSideMatrix,
                                                           const ProcessInfo&   rProcessInfo)
{
    switch (NumberOfUDofs()) {
    case 8:
        CalculateAndAssignStiffnessMatrix<8>(rLeftHandSideMatrix, rProcessInfo);
        break;
    case 12:
        CalculateAndAssignStiffnessMatrix<12>(rLeftHandSideMatrix, rProcessInfo);
        break;
    case 18:
        CalculateAndAssignStiffnessMatrix<18>(rLeftHandSideMatrix, rProcessInfo);
        break;
    case 36:
        CalculateAndAssignStiffnessMatrix<36>(rLeftHandSideMatrix, rProcessInfo);
        break;
    case 24:
        CalculateAndAssignStiffnessMatrix<24>(rLeftHandSideMatrix, rProcessInfo);
        break;
    case 48:
        CalculateAndAssignStiffnessMatrix<48>(rLeftHandSideMatrix, rProcessInfo);
        break;
    default:
        KRATOS_ERROR << "This stiffness matrix size is not supported: " << NumberOfUDofs() << "\n";
    }
}

void UPwInterfaceElement::CalculateAndAssignPermeabilityMatrix(Element::MatrixType& rLeftHandSideMatrix)
{
    switch (GetGeometry().PointsNumber()) {
    case 4:
        CalculateAndAssignPermeabilityMatrix<4>(rLeftHandSideMatrix);
        break;
    case 6:
        CalculateAndAssignPermeabilityMatrix<6>(rLeftHandSideMatrix);
        break;
    case 8:
        CalculateAndAssignPermeabilityMatrix<8>(rLeftHandSideMatrix);
        break;
    case 12:
        CalculateAndAssignPermeabilityMatrix<12>(rLeftHandSideMatrix);
        break;
    case 16:
        CalculateAndAssignPermeabilityMatrix<16>(rLeftHandSideMatrix);
        break;
    default:
        KRATOS_ERROR
            << "This number of points is not supported for the permeability matrix calculation: "
            << GetGeometry().PointsNumber() << "\n";
    }
}

void UPwInterfaceElement::CalculateRightHandSide(Element::VectorType& rRightHandSideVector,
                                                 const ProcessInfo&   rProcessInfo)
{
    rRightHandSideVector = ZeroVector{GetDofs().size()};

    // Currently, the right-hand side only includes the internal force vector. In the future, it
    // will also include water pressure contributions and coupling terms.
    for (auto contribution : mContributions) {
        switch (contribution) {
        case CalculationContribution::Stiffness:
            CalculateAndAssembleStifnessForceVector(rRightHandSideVector, rProcessInfo);
            break;
        case CalculationContribution::Permeability:
            CalculateAndAssemblePermeabilityFlowVector(rRightHandSideVector);
            break;
        case CalculationContribution::FluidBodyFlow:
            CalculateAndAssembleFluidBodyFlowVector(rRightHandSideVector);
            break;
        default:
            KRATOS_ERROR << "This contribution is not supported \n";
        }
    }
}

void UPwInterfaceElement::CalculateAndAssembleStifnessForceVector(Element::VectorType& rRightHandSideVector,
                                                                  const ProcessInfo& rProcessInfo)
{
    switch (NumberOfUDofs()) {
    case 8:
        CalculateAndAssembleStiffnesForceVector<8>(rRightHandSideVector, rProcessInfo);
        break;
    case 12:
        CalculateAndAssembleStiffnesForceVector<12>(rRightHandSideVector, rProcessInfo);
        break;
    case 18:
        CalculateAndAssembleStiffnesForceVector<18>(rRightHandSideVector, rProcessInfo);
        break;
    case 36:
        CalculateAndAssembleStiffnesForceVector<36>(rRightHandSideVector, rProcessInfo);
        break;
    case 24:
        CalculateAndAssembleStiffnesForceVector<24>(rRightHandSideVector, rProcessInfo);
        break;
    case 48:
        CalculateAndAssembleStiffnesForceVector<48>(rRightHandSideVector, rProcessInfo);
        break;
    default:
        KRATOS_ERROR << "This stiffness force vector size is not supported: " << NumberOfUDofs() << "\n";
    }
}

void UPwInterfaceElement::CalculateAndAssemblePermeabilityFlowVector(Element::VectorType& rRightHandSideVector)
{
    switch (GetWaterPressureGeometry().PointsNumber()) {
    case 4:
        CalculateAndAssemblePermeabilityFlowVector<4>(rRightHandSideVector);
        break;
    case 6:
        CalculateAndAssemblePermeabilityFlowVector<6>(rRightHandSideVector);
        break;
    case 8:
        CalculateAndAssemblePermeabilityFlowVector<8>(rRightHandSideVector);
        break;
    case 12:
        CalculateAndAssemblePermeabilityFlowVector<12>(rRightHandSideVector);
        break;
    case 16:
        CalculateAndAssemblePermeabilityFlowVector<16>(rRightHandSideVector);
        break;
    default:
        KRATOS_ERROR << "This number of points for the permeability flow vector is not supported: "
                     << GetWaterPressureGeometry().PointsNumber() << "\n";
    }
}

void UPwInterfaceElement::CalculateAndAssembleFluidBodyFlowVector(Element::VectorType& rRightHandSideVector)
{
    switch (GetWaterPressureGeometry().PointsNumber()) {
    case 4:
        CalculateAndAssembleFluidBodyFlowVector<4>(rRightHandSideVector);
        break;
    case 6:
        CalculateAndAssembleFluidBodyFlowVector<6>(rRightHandSideVector);
        break;
    case 8:
        CalculateAndAssembleFluidBodyFlowVector<8>(rRightHandSideVector);
        break;
    case 12:
        CalculateAndAssembleFluidBodyFlowVector<12>(rRightHandSideVector);
        break;
    case 16:
        CalculateAndAssembleFluidBodyFlowVector<16>(rRightHandSideVector);
        break;
    default:
        KRATOS_ERROR << "This number of points for the fluid body flow vector is not supported: "
                     << GetWaterPressureGeometry().PointsNumber() << "\n";
    }
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
        rOutput = StressStrainUtilities::CalculateStressVectorsFromStrainVectors(
            relative_displacements, rCurrentProcessInfo, GetProperties(), mConstitutiveLaws);
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
    const auto tractions = StressStrainUtilities::CalculateStressVectorsFromStrainVectors(
        relative_displacements, rProcessInfo, GetProperties(), mConstitutiveLaws);
    const auto integration_coefficients = CalculateIntegrationCoefficients();
    if (rVariable == INTERNAL_FORCES_VECTOR) {
        GeoElementUtilities::AssignUBlockVector(
            rOutput, GeoEquationOfMotionUtilities::CalculateInternalForceVector(
                         local_b_matrices, tractions, integration_coefficients));
        // Todo: Extend with permeability flow and other p parts
    } else if (rVariable == EXTERNAL_FORCES_VECTOR) {
        GeoElementUtilities::AssignUBlockVector(rOutput, Vector{NumberOfUDofs(), 0.0});
        // Todo: Extend with p parts
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
    mRetentionLaws.clear();
    for (auto i = std::size_t{0}; i < mpIntegrationScheme->GetNumberOfIntegrationPoints(); ++i) {
        mConstitutiveLaws.push_back(GetProperties()[CONSTITUTIVE_LAW]->Clone());
        mRetentionLaws.push_back(RetentionLawFactory::Clone(GetProperties()));
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

const Geometry<Node>& UPwInterfaceElement::GetWaterPressureMidGeometry() const
{
    constexpr auto unused_part_index = std::size_t{0};
    return GetWaterPressureGeometry().GetGeometryPart(unused_part_index);
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
    return mpOptionalPressureGeometry ? *mpOptionalPressureGeometry : GetDisplacementGeometry();
}

std::size_t UPwInterfaceElement::NumberOfUDofs() const
{
    return GetDisplacementGeometry().size() * GetDisplacementGeometry().WorkingSpaceDimension();
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

Matrix UPwInterfaceElement::CalculatePwBMatrix(const Vector& rN, const Geometry<Node>& rGeometry) const
{
    KRATOS_ERROR_IF(rN.empty())
        << "Shape function values are empty. Therefore, the PwB matrix can not be computed.\n";
    KRATOS_ERROR_IF_NOT(rN.size() == rGeometry.size() / 2)
        << "The number of shape functions should be equal to the number of node pairs. Therefore, "
           "the PwB matrix can not be computed.\n";
    auto component_order = ComponentOrder(rGeometry);
    auto result          = Matrix{rGeometry.size(), component_order.size(), 0.0};

    auto number_of_pw_dofs_per_side = result.size1() / 2;
    for (auto i = size_t{0}; i < number_of_pw_dofs_per_side; ++i) {
        for (auto j = size_t{0}; j < component_order.size(); ++j) {
            result(i, component_order[j])                              = -rN[i];
            result(i + number_of_pw_dofs_per_side, component_order[j]) = rN[i];
        }
    }
    return result;
}

std::vector<std::size_t> UPwInterfaceElement::ComponentOrder(const Geometry<Node>& rGeometry) const
{
    return rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear
               ? std::vector<std::size_t>{1, 0}
               : std::vector<std::size_t>{2, 0, 1};
}

Geometry<Node>::ShapeFunctionsGradientsType UPwInterfaceElement::CalculateLocalPwBMatricesAtIntegrationPoints() const
{
    const auto& r_integration_points = mpIntegrationScheme->GetIntegrationPoints();
    const auto  shape_function_values_at_integration_points =
        GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(r_integration_points,
                                                                       GetWaterPressureGeometry());

    auto result = Geometry<Node>::ShapeFunctionsGradientsType{r_integration_points.size()};
    // result.reserve(shape_function_values_at_integration_points.size());
    auto calculate_local_b_matrix = [&r_geometry = GetWaterPressureGeometry(), this](
                                        const auto& rShapeFunctionValuesAtIntegrationPoint,
                                        const auto& rIntegrationPoint) {
        auto b_matrix = CalculatePwBMatrix(rShapeFunctionValuesAtIntegrationPoint, r_geometry);
        ApplyRotationToBMatrix(b_matrix, mfpCalculateRotationMatrix(r_geometry, rIntegrationPoint));
        return b_matrix;
    };
    std::ranges::transform(shape_function_values_at_integration_points, r_integration_points,
                           result.begin(), calculate_local_b_matrix);

    return result;
}

std::vector<double> UPwInterfaceElement::CalculateIntegrationCoefficients() const
{
    const auto determinants_of_jacobian = CalculateDeterminantsOfJacobiansAtIntegrationPoints(
        mpIntegrationScheme->GetIntegrationPoints(), GetDisplacementGeometry());
    return mIntegrationCoefficientsCalculator.Run<>(mpIntegrationScheme->GetIntegrationPoints(),
                                                    determinants_of_jacobian, this);
}

std::vector<Vector> UPwInterfaceElement::CalculateRelativeDisplacementsAtIntegrationPoints(
    const std::vector<Matrix>& rLocalBMatrices) const
{
    const Geometry<Node> no_Pw_geometry;
    const auto           u_dofs = Geo::DofUtilities::ExtractUPwDofsFromNodes(
        GetDisplacementGeometry(), no_Pw_geometry, GetDisplacementGeometry().WorkingSpaceDimension());
    auto nodal_displacement_vector = Vector{u_dofs.size()};
    std::ranges::transform(u_dofs, nodal_displacement_vector.begin(),
                           [](auto p_dof) { return p_dof->GetSolutionStepValue(); });

    auto result = std::vector<Vector>{};
    result.reserve(rLocalBMatrices.size());
    auto calculate_relative_displacement_vector = [&nodal_displacement_vector](const auto& rLocalB) {
        return Vector{prod(rLocalB, nodal_displacement_vector)};
    };
    std::ranges::transform(rLocalBMatrices, std::back_inserter(result), calculate_relative_displacement_vector);

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
    auto rotation_tensor = Matrix{IdentityMatrix(3, 3)};
    if (GetDisplacementGeometry().GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Linear) {
        GeoElementUtilities::AssignMatrixAtPosition(
            rotation_tensor, mfpCalculateRotationMatrix(GetDisplacementGeometry(), rIntegrationPoint), 0, 0);
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

Vector UPwInterfaceElement::GetWaterPressureGeometryNodalVariable(const Variable<double>& rVariable) const
{
    Vector result{this->GetWaterPressureGeometry().size()};
    auto   index = std::size_t{0};
    for (auto& r_node : this->GetWaterPressureGeometry()) {
        result[index] = r_node.FastGetSolutionStepValue(rVariable);
        index++;
    }
    return result;
}

std::vector<double> UPwInterfaceElement::CalculateIntegrationPointFluidPressures() const
{
    Matrix n_container{mpIntegrationScheme->GetIntegrationPoints().size(),
                       GetWaterPressureMidGeometry().PointsNumber()};
    auto   integration_point_index = std::size_t{0};
    for (auto& r_integration_point : mpIntegrationScheme->GetIntegrationPoints()) {
        auto integration_point_shape_function_values = Vector{};
        // water pressure shape function values on integration point ( the integration points are shared with the displacement mid geometry )
        GetWaterPressureMidGeometry().ShapeFunctionsValues(integration_point_shape_function_values,
                                                           r_integration_point);
        row(n_container, integration_point_index) = integration_point_shape_function_values;
        integration_point_index++;
    }
    auto   nodal_water_pressure = GetWaterPressureGeometryNodalVariable(WATER_PRESSURE);
    Vector mid_geometry_water_pressures{nodal_water_pressure.size() / 2};
    for (auto i = std::size_t{0}; i < nodal_water_pressure.size() / 2; ++i) {
        // this choice of averaging should be different for impermeable interface
        mid_geometry_water_pressures[i] =
            (nodal_water_pressure[i] + nodal_water_pressure[i + nodal_water_pressure.size() / 2]) / 2.0;
    }
    return GeoTransportEquationUtilities::CalculateFluidPressures(n_container, mid_geometry_water_pressures);
}

std::vector<Vector> UPwInterfaceElement::CalculateProjectedGravity() const
{
    const auto& r_integration_points = mpIntegrationScheme->GetIntegrationPoints();
    const auto  shape_function_values_at_integration_points =
        GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(
            r_integration_points, GetWaterPressureMidGeometry());
    const auto number_integration_points = r_integration_points.size();
    const auto dimension                 = GetWaterPressureMidGeometry().WorkingSpaceDimension();

    // Get volume acceleration from WaterPressureGeometry and average to WaterPressureMidGeometry
    std::vector<Vector> volume_accelerations;
    volume_accelerations.reserve(GetWaterPressureGeometry().PointsNumber());
    for (const auto& r_node : GetWaterPressureGeometry()) {
        volume_accelerations.emplace_back(r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION));
    }
    auto component_order = ComponentOrder(GetWaterPressureGeometry());

    // average to WaterPressureMidGeometry and sort to element directions
    std::vector<Vector> mid_volume_accelerations;
    const auto          number_of_mid_points = GetWaterPressureMidGeometry().PointsNumber();
    mid_volume_accelerations.reserve(number_of_mid_points);
    for (auto i = std::size_t{0}; i < number_of_mid_points; ++i) {
        auto mean_acceleration = Vector{dimension};
        for (auto idim = std::size_t{0}; idim < dimension; ++idim) {
            mean_acceleration[component_order[idim]] =
                (volume_accelerations[i][idim] + volume_accelerations[i + number_of_mid_points][idim]) / 2.0;
        }
        mid_volume_accelerations.emplace_back(mean_acceleration);
    }

    std::vector<Vector> projected_gravity;
    projected_gravity.reserve(number_integration_points);
    for (auto integration_point_index = std::size_t{0};
         integration_point_index < number_integration_points; ++integration_point_index) {
        auto body_acceleration = Vector{dimension, 0.0};
        for (auto mid_node_index = std::size_t{0}; mid_node_index < number_of_mid_points; ++mid_node_index) {
            body_acceleration +=
                shape_function_values_at_integration_points[integration_point_index][mid_node_index] *
                mid_volume_accelerations[mid_node_index];
        }
        projected_gravity.emplace_back(body_acceleration);
    }

    return projected_gravity;
}

Geo::BMatricesGetter UPwInterfaceElement::CreateBMatricesGetter() const
{
    return [this]() { return this->CalculateLocalBMatricesAtIntegrationPoints(); };
}

Geo::ConstitutiveLawsGetter UPwInterfaceElement::CreateConstitutiveLawsGetter() const
{
    return
        [this]() -> const std::vector<ConstitutiveLaw::Pointer>& { return this->mConstitutiveLaws; };
}

Geo::RetentionLawsGetter UPwInterfaceElement::CreateRetentionLawsGetter() const
{
    return [this]() -> const std::vector<RetentionLaw::Pointer>& { return this->mRetentionLaws; };
}

Geo::MaterialPermeabilityMatrixGetter UPwInterfaceElement::CreateMaterialPermeabilityGetter() const
{
    return [this]() {
        return GeoElementUtilities::FillInterfacePermeabilityMatrix(
            this->GetProperties(), this->GetWaterPressureMidGeometry().WorkingSpaceDimension());
    };
}

Geo::StrainVectorsGetter UPwInterfaceElement::CreateRelativeDisplacementsGetter() const
{
    return [this]() {
        return this->CalculateRelativeDisplacementsAtIntegrationPoints(
            this->CalculateLocalBMatricesAtIntegrationPoints());
    };
}

Geo::IntegrationCoefficientsGetter UPwInterfaceElement::CreateIntegrationCoefficientsGetter() const
{
    return [this]() { return this->CalculateIntegrationCoefficients(); };
}

Geo::NodalValuesGetter UPwInterfaceElement::CreateWaterPressureGeometryNodalVariableGetter() const
{
    return [this](const Variable<double>& rVariable) {
        return this->GetWaterPressureGeometryNodalVariable(rVariable);
    };
}

Geo::PropertiesGetter UPwInterfaceElement::CreatePropertiesGetter() const
{
    return [this]() -> const Properties& { return this->GetProperties(); };
}

Geo::ShapeFunctionGradientsGetter UPwInterfaceElement::CreatePwBMatricesGetter() const
{
    return [this]() { return this->CalculateLocalPwBMatricesAtIntegrationPoints(); };
}

Geo::IntegrationPointValuesGetter UPwInterfaceElement::CreateFluidPressureCalculator() const
{
    return [this]() { return this->CalculateIntegrationPointFluidPressures(); };
}

std::function<std::vector<Vector>()> UPwInterfaceElement::CreateProjectedGravityCalculator() const
{
    return [this]() { return this->CalculateProjectedGravity(); };
}

template <unsigned int MatrixSize>
typename StiffnessCalculator<MatrixSize>::InputProvider UPwInterfaceElement::CreateStiffnessInputProvider(const ProcessInfo& rProcessInfo)
{
    return typename StiffnessCalculator<MatrixSize>::InputProvider(
        CreateBMatricesGetter(), CreateRelativeDisplacementsGetter(), CreateIntegrationCoefficientsGetter(),
        CreatePropertiesGetter(), CreateProcessInfoGetter(rProcessInfo), CreateConstitutiveLawsGetter());
}

template <unsigned int MatrixSize>
auto UPwInterfaceElement::CreateStiffnessCalculator(const ProcessInfo& rProcessInfo)
{
    return StiffnessCalculator<MatrixSize>(CreateStiffnessInputProvider<MatrixSize>(rProcessInfo));
}

template <unsigned int MatrixSize>
void UPwInterfaceElement::CalculateAndAssignStiffnessMatrix(MatrixType&        rLeftHandSideMatrix,
                                                            const ProcessInfo& rProcessInfo)
{
    GeoElementUtilities::AssignUUBlockMatrix(
        rLeftHandSideMatrix, CreateStiffnessCalculator<MatrixSize>(rProcessInfo).LHSContribution().value());
}

template <unsigned int MatrixSize>
void UPwInterfaceElement::CalculateAndAssembleStiffnesForceVector(VectorType& rRightHandSideVector,
                                                                  const ProcessInfo& rProcessInfo)
{
    GeoElementUtilities::AssembleUBlockVector(
        rRightHandSideVector, CreateStiffnessCalculator<MatrixSize>(rProcessInfo).RHSContribution());
}

template <unsigned int TNumNodes>
typename PermeabilityCalculator<TNumNodes>::InputProvider UPwInterfaceElement::CreatePermeabilityInputProvider()
{
    return typename PermeabilityCalculator<TNumNodes>::InputProvider(
        CreatePropertiesGetter(), CreateRetentionLawsGetter(), CreateMaterialPermeabilityGetter(),
        CreateIntegrationCoefficientsGetter(), CreateWaterPressureGeometryNodalVariableGetter(),
        CreatePwBMatricesGetter(), CreateFluidPressureCalculator());
}

template <unsigned int TNumNodes>
auto UPwInterfaceElement::CreatePermeabilityCalculator()
{
    return PermeabilityCalculator<TNumNodes>(CreatePermeabilityInputProvider<TNumNodes>());
}

template <unsigned int TNumNodes>
void UPwInterfaceElement::CalculateAndAssignPermeabilityMatrix(MatrixType& rLeftHandSideMatrix)
{
    GeoElementUtilities::AssignPPBlockMatrix(
        rLeftHandSideMatrix, CreatePermeabilityCalculator<TNumNodes>().LHSContribution().value());
}

template <unsigned int TNumNodes>
void UPwInterfaceElement::CalculateAndAssemblePermeabilityFlowVector(VectorType& rRightHandSideVector)
{
    GeoElementUtilities::AssemblePBlockVector(
        rRightHandSideVector, CreatePermeabilityCalculator<TNumNodes>().RHSContribution());
}

template <unsigned int TNumNodes>
typename FluidBodyFlowCalculator<TNumNodes>::InputProvider UPwInterfaceElement::CreateFluidBodyFlowInputProvider()
{
    return typename FluidBodyFlowCalculator<TNumNodes>::InputProvider(
        CreatePropertiesGetter(), CreateRetentionLawsGetter(), CreateMaterialPermeabilityGetter(),
        CreateIntegrationCoefficientsGetter(), CreateProjectedGravityCalculator(),
        CreatePwBMatricesGetter(), CreateFluidPressureCalculator());
}

template <unsigned int TNumNodes>
auto UPwInterfaceElement::CreateFluidBodyFlowCalculator()
{
    return FluidBodyFlowCalculator<TNumNodes>(CreateFluidBodyFlowInputProvider<TNumNodes>());
}

template <unsigned int TNumNodes>
void UPwInterfaceElement::CalculateAndAssembleFluidBodyFlowVector(VectorType& rRightHandSideVector)
{
    GeoElementUtilities::AssemblePBlockVector(
        rRightHandSideVector, CreateFluidBodyFlowCalculator<TNumNodes>().RHSContribution());
}

// Instances of this class can not be copied but can be moved. Check that at compile time.
static_assert(!std::is_copy_constructible_v<UPwInterfaceElement>);
static_assert(!std::is_copy_assignable_v<UPwInterfaceElement>);
static_assert(std::is_move_constructible_v<UPwInterfaceElement>);
static_assert(std::is_move_assignable_v<UPwInterfaceElement>);
} // namespace Kratos
