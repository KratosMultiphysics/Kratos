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
    mRetentionLawVector.resize(mpIntegrationScheme->GetIntegrationPoints().size());
    auto properties = Properties{};
    if (rpProperties) properties = *rpProperties;
    for (auto& r_retention_law : mRetentionLawVector) {
        r_retention_law = RetentionLawFactory::Clone(properties);
    }
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

    // Currently, the left-hand side matrix includes the stiffness matrix and U Pw coupling terms.
    // In the future, it will also include water pressure contributions.
    for (auto contribution : mContributions) {
        switch (contribution) {
        case CalculationContribution::Stiffness:
            CalculateAndAssignStiffnessMatrix(rLeftHandSideMatrix, rProcessInfo);
            break;
        case CalculationContribution::UPCoupling:
            CalculateAndAssignUPCouplingMatrix(rLeftHandSideMatrix);
            break;
        default:
            KRATOS_ERROR << "This contribution is not supported \n";
        }
    }
}

void UPwInterfaceElement::CalculateAndAssignStiffnessMatrix(Element::MatrixType& rLeftHandSideMatrix,
                                                            const ProcessInfo& rProcessInfo)
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

void UPwInterfaceElement::CalculateAndAssignUPCouplingMatrix(MatrixType& rLeftHandSideMatrix) const
{
    using Key  = std::pair<std::size_t, std::size_t>;
    using Func = void (UPwInterfaceElement::*)(MatrixType&) const;

    static const std::map<Key, Func> dispatch_table = {
        {{8, 4}, &UPwInterfaceElement::CalculateAndAssignUPCouplingMatrix<8, 4>},
        {{12, 6}, &UPwInterfaceElement::CalculateAndAssignUPCouplingMatrix<12, 6>},
        {{12, 4}, &UPwInterfaceElement::CalculateAndAssignUPCouplingMatrix<12, 4>},
        {{18, 6}, &UPwInterfaceElement::CalculateAndAssignUPCouplingMatrix<18, 6>},
        {{36, 12}, &UPwInterfaceElement::CalculateAndAssignUPCouplingMatrix<36, 12>},
        {{36, 6}, &UPwInterfaceElement::CalculateAndAssignUPCouplingMatrix<36, 6>},
        {{24, 8}, &UPwInterfaceElement::CalculateAndAssignUPCouplingMatrix<24, 8>},
        {{48, 16}, &UPwInterfaceElement::CalculateAndAssignUPCouplingMatrix<48, 16>},
        {{48, 8}, &UPwInterfaceElement::CalculateAndAssignUPCouplingMatrix<48, 8>}};

    const auto key = Key{NumberOfUDofs(), GetWaterPressureGeometry().size()};
    KRATOS_ERROR_IF_NOT(dispatch_table.contains(key))
        << "This Coupling matrix size is not supported: " << key.first << "x" << key.second << "\n";
    (this->*dispatch_table.at(key))(rLeftHandSideMatrix);
}

void UPwInterfaceElement::CalculateRightHandSide(Element::VectorType& rRightHandSideVector,
                                                 const ProcessInfo&   rProcessInfo)
{
    rRightHandSideVector = ZeroVector{GetDofs().size()};

    // Currently, the right-hand side includes the internal force vector and U Pw coupling terms. In
    // the future, it will also include water pressure contributions.
    for (auto contribution : mContributions) {
        switch (contribution) {
        case CalculationContribution::Stiffness:
            CalculateAndAssembleStiffnessForceVector(rRightHandSideVector, rProcessInfo);
            break;
        case CalculationContribution::UPCoupling:
            CalculateAndAssembleUPCouplingForceVector(rRightHandSideVector);
            break;
        default:
            KRATOS_ERROR << "This contribution is not supported \n";
        }
    }
}

void UPwInterfaceElement::CalculateAndAssembleStiffnessForceVector(Element::VectorType& rRightHandSideVector,
                                                                   const ProcessInfo& rProcessInfo)
{
    switch (NumberOfUDofs()) {
    case 8:
        CalculateAndAssembleStiffnessForceVector<8>(rRightHandSideVector, rProcessInfo);
        break;
    case 12:
        CalculateAndAssembleStiffnessForceVector<12>(rRightHandSideVector, rProcessInfo);
        break;
    case 18:
        CalculateAndAssembleStiffnessForceVector<18>(rRightHandSideVector, rProcessInfo);
        break;
    case 36:
        CalculateAndAssembleStiffnessForceVector<36>(rRightHandSideVector, rProcessInfo);
        break;
    case 24:
        CalculateAndAssembleStiffnessForceVector<24>(rRightHandSideVector, rProcessInfo);
        break;
    case 48:
        CalculateAndAssembleStiffnessForceVector<48>(rRightHandSideVector, rProcessInfo);
        break;
    default:
        KRATOS_ERROR << "This stiffness force vector size is not supported: " << NumberOfUDofs() << "\n";
    }
}

void UPwInterfaceElement::CalculateAndAssembleUPCouplingForceVector(Element::VectorType& rRightHandSideVector) const
{
    using Key  = std::pair<std::size_t, std::size_t>;
    using Func = void (UPwInterfaceElement::*)(Element::VectorType&) const;

    static const std::map<Key, Func> dispatch_table = {
        {{8, 4}, &UPwInterfaceElement::CalculateAndAssembleUPCouplingForceVector<8, 4>},
        {{12, 6}, &UPwInterfaceElement::CalculateAndAssembleUPCouplingForceVector<12, 6>},
        {{12, 4}, &UPwInterfaceElement::CalculateAndAssembleUPCouplingForceVector<12, 4>},
        {{18, 6}, &UPwInterfaceElement::CalculateAndAssembleUPCouplingForceVector<18, 6>},
        {{36, 12}, &UPwInterfaceElement::CalculateAndAssembleUPCouplingForceVector<36, 12>},
        {{36, 6}, &UPwInterfaceElement::CalculateAndAssembleUPCouplingForceVector<36, 6>},
        {{24, 8}, &UPwInterfaceElement::CalculateAndAssembleUPCouplingForceVector<24, 8>},
        {{48, 16}, &UPwInterfaceElement::CalculateAndAssembleUPCouplingForceVector<48, 16>},
        {{48, 8}, &UPwInterfaceElement::CalculateAndAssembleUPCouplingForceVector<48, 8>}};

    const auto key = Key{NumberOfUDofs(), GetWaterPressureGeometry().size()};
    KRATOS_ERROR_IF_NOT(dispatch_table.contains(key))
        << "This coupling force vector size is not supported: " << key.first << "x" << key.second << "\n";
    (this->*dispatch_table.at(key))(rRightHandSideVector);
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
    } else if (rVariable == EXTERNAL_FORCES_VECTOR) {
        GeoElementUtilities::AssignUBlockVector(rOutput, Vector{NumberOfUDofs(), 0.0});
    }
}

const Geometry<Node>& UPwInterfaceElement::GetWaterPressureMidGeometry() const
{
    constexpr auto unused_part_index = std::size_t{0};
    return GetWaterPressureGeometry().GetGeometryPart(unused_part_index);
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
        r_neighbour_element.CalculateOnIntegrationPoints(
            CAUCHY_STRESS_VECTOR, neighbour_cauchy_stresses, rCurrentProcessInfo);
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

Geo::BMatricesGetter UPwInterfaceElement::CreateBMatricesGetter() const
{
    return [this]() { return this->CalculateLocalBMatricesAtIntegrationPoints(); };
}

Geo::ConstitutiveLawsGetter UPwInterfaceElement::CreateConstitutiveLawsGetter() const
{
    return
        [this]() -> const std::vector<ConstitutiveLaw::Pointer>& { return this->mConstitutiveLaws; };
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

Geo::PropertiesGetter UPwInterfaceElement::CreatePropertiesGetter() const
{
    return [this]() -> const Properties& { return this->GetProperties(); };
}

std::function<const Matrix()> UPwInterfaceElement::CreateNpContainerGetter() const
{
    return [this]() { return this->GetNpContainer(); };
}

std::function<Vector()> UPwInterfaceElement::CreateVoigtVectorGetter() const
{
    return [this]() { return mpStressStatePolicy->GetVoigtVector(); };
}

std::function<std::vector<double>()> UPwInterfaceElement::CreateBiotCoefficientsGetter() const
{
    return [this]() { return this->CalculateBiotCoefficients(); };
}

std::vector<double> UPwInterfaceElement::CalculateBiotCoefficients() const
{
    const auto& r_properties = this->GetProperties();
    const double biot_coefficient = r_properties.Has(BIOT_COEFFICIENT) ? r_properties[BIOT_COEFFICIENT] : 1.0;
    const std::size_t n_points = mpIntegrationScheme->GetNumberOfIntegrationPoints();
    return std::vector<double>(n_points, biot_coefficient);
}

std::function<std::vector<double>()> UPwInterfaceElement::CreateBishopCoefficientsGetter() const
{
    return [this]() { return this->CalculateBishopCoefficients(); };
}

std::vector<double> UPwInterfaceElement::CalculateBishopCoefficients() const
{
    const auto fluid_pressure = CalculateIntegrationPointFluidPressures();
    KRATOS_ERROR_IF_NOT(fluid_pressure.size() == mRetentionLawVector.size());

    auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};

    auto result = std::vector<double>{};
    result.reserve(mRetentionLawVector.size());
    std::transform(mRetentionLawVector.begin(), mRetentionLawVector.end(), fluid_pressure.begin(),
                   std::back_inserter(result),
                   [&retention_law_params](const auto& pRetentionLaw, auto FluidPressure) {
        retention_law_params.SetFluidPressure(FluidPressure);
        return pRetentionLaw->CalculateBishopCoefficient(retention_law_params);
    });
    return result;
}

Matrix UPwInterfaceElement::GetNpContainer() const
{
    const auto total_number_of_nodes    = GetWaterPressureGeometry().size();
    const auto number_of_pressure_nodes = GetWaterPressureMidGeometry().PointsNumber();
    Matrix n_container{mpIntegrationScheme->GetIntegrationPoints().size(), total_number_of_nodes};
    auto   shape_function_values_interface = Vector(total_number_of_nodes, 0.0);

    auto integration_point_index = std::size_t{0};
    for (auto& r_integration_point : mpIntegrationScheme->GetIntegrationPoints()) {
        auto integration_point_shape_function_values = Vector{};
        // water pressure shape function values on integration point ( the integration points are shared with the displacement mid-geometry )
        GetWaterPressureMidGeometry().ShapeFunctionsValues(integration_point_shape_function_values,
                                                           r_integration_point);
        // to interplate nodal values of water pressure on mid-geometry, the shape function is split into two equal contributions
        noalias(subrange(shape_function_values_interface, 0, number_of_pressure_nodes)) =
            0.5 * integration_point_shape_function_values;
        noalias(subrange(shape_function_values_interface, number_of_pressure_nodes, 2 * number_of_pressure_nodes)) =
            subrange(shape_function_values_interface, 0, number_of_pressure_nodes);
        row(n_container, integration_point_index) = shape_function_values_interface;
        integration_point_index++;
    }
    return n_container;
}

Vector UPwInterfaceElement::CalculateIntegrationPointFluidPressures() const
{
    const auto fluid_pressure = GeoTransportEquationUtilities::CalculateFluidPressures(
        GetNpContainer(), GetWaterPressureGeometryNodalVariable());
    auto result = Vector(fluid_pressure.size());
    std::ranges::copy(fluid_pressure, result.begin());
    return result;
}

std::function<Vector()> UPwInterfaceElement::CreateNodalPressuresGetter() const
{
    return [this]() { return this->GetWaterPressureGeometryNodalVariable(); };
}

Vector UPwInterfaceElement::GetWaterPressureGeometryNodalVariable() const
{
    Vector result{this->GetWaterPressureGeometry().size()};
    VariablesUtilities::GetNodalValues(this->GetWaterPressureGeometry(), WATER_PRESSURE, result.begin());
    return result;
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
void UPwInterfaceElement::CalculateAndAssembleStiffnessForceVector(VectorType& rRightHandSideVector,
                                                                   const ProcessInfo& rProcessInfo)
{
    GeoElementUtilities::AssembleUBlockVector(
        rRightHandSideVector, CreateStiffnessCalculator<MatrixSize>(rProcessInfo).RHSContribution());
}

template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
typename UPCouplingCalculator<NumberOfRows, NumberOfColumns>::InputProvider UPwInterfaceElement::CreateUPCouplingInputProvider() const
{
    return typename UPCouplingCalculator<NumberOfRows, NumberOfColumns>::InputProvider(
        CreateNpContainerGetter(), CreateBMatricesGetter(), CreateVoigtVectorGetter(),
        CreateIntegrationCoefficientsGetter(), CreateBiotCoefficientsGetter(),
        CreateBishopCoefficientsGetter(), CreateNodalPressuresGetter());
}

template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
auto UPwInterfaceElement::CreateUPCouplingCalculator() const
{
    return UPCouplingCalculator<NumberOfRows, NumberOfColumns>(
        CreateUPCouplingInputProvider<NumberOfRows, NumberOfColumns>());
}

template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
void UPwInterfaceElement::CalculateAndAssignUPCouplingMatrix(MatrixType& rLeftHandSideMatrix) const
{
    GeoElementUtilities::AssignUPBlockMatrix(
        rLeftHandSideMatrix,
        CreateUPCouplingCalculator<NumberOfRows, NumberOfColumns>().LHSContribution().value());
}

template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
void UPwInterfaceElement::CalculateAndAssembleUPCouplingForceVector(VectorType& rRightHandSideVector) const
{
    GeoElementUtilities::AssembleUBlockVector(
        rRightHandSideVector,
        (-1.0) * CreateUPCouplingCalculator<NumberOfRows, NumberOfColumns>().RHSContribution());
}

// Instances of this class can not be copied but can be moved. Check that at compile time.
static_assert(!std::is_copy_constructible_v<UPwInterfaceElement>);
static_assert(!std::is_copy_assignable_v<UPwInterfaceElement>);
static_assert(std::is_move_constructible_v<UPwInterfaceElement>);
static_assert(std::is_move_assignable_v<UPwInterfaceElement>);
} // namespace Kratos
