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
#include "interface_element.h"

#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/equation_of_motion_utilities.h"
#include "custom_utilities/extrapolation_utilities.h"
#include "custom_utilities/geometry_utilities.h"
#include "custom_utilities/math_utilities.h"
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

} // namespace

namespace Kratos
{

InterfaceElement::InterfaceElement(IndexType                                             NewId,
                                   const Geometry<GeometricalObject::NodeType>::Pointer& rGeometry,
                                   const Properties::Pointer&         rProperties,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy)
    : Element(NewId, rGeometry, rProperties), mpStressStatePolicy(std::move(pStressStatePolicy))
{
    MakeIntegrationSchemeAndAssignFunction();
}

InterfaceElement::InterfaceElement(IndexType                          NewId,
                                   const GeometryType::Pointer&       rGeometry,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy)
    : Element(NewId, rGeometry), mpStressStatePolicy(std::move(pStressStatePolicy))
{
    MakeIntegrationSchemeAndAssignFunction();
}

void InterfaceElement::MakeIntegrationSchemeAndAssignFunction()
{
    if (GetGeometry().LocalSpaceDimension() == 1) {
        mIntegrationScheme = std::make_unique<LobattoIntegrationScheme>(GetGeometry().PointsNumber() / 2);
        mfpCalculateRotationMatrix = GeometryUtilities::Calculate2DRotationMatrixForLineGeometry;
    } else {
        mIntegrationScheme = std::make_unique<LumpedIntegrationScheme>(GetGeometry().PointsNumber() / 2);
        mfpCalculateRotationMatrix = GeometryUtilities::Calculate3DRotationMatrixForPlaneGeometry;
    }
}

Element::Pointer InterfaceElement::Create(IndexType               NewId,
                                          const NodesArrayType&   rNodes,
                                          PropertiesType::Pointer pProperties) const
{
    return Create(NewId, this->GetGeometry().Create(rNodes), pProperties);
}

Element::Pointer InterfaceElement::Create(IndexType               NewId,
                                          GeometryType::Pointer   pGeometry,
                                          PropertiesType::Pointer pProperties) const
{
    return make_intrusive<InterfaceElement>(NewId, pGeometry, pProperties, mpStressStatePolicy->Clone());
}

void InterfaceElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

void InterfaceElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rProcessInfo)
{
    // Currently, the left-hand side matrix only includes the stiffness matrix. In the future, it
    // will also include water pressure contributions and coupling terms.
    const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();
    rLeftHandSideMatrix         = GeoEquationOfMotionUtilities::CalculateStiffnessMatrix(
        local_b_matrices,
        CalculateConstitutiveMatricesAtIntegrationPoints(
            CalculateRelativeDisplacementsAtIntegrationPoints(local_b_matrices), rProcessInfo),
        CalculateIntegrationCoefficients());
}

void InterfaceElement::CalculateRightHandSide(Element::VectorType& rRightHandSideVector, const ProcessInfo& rProcessInfo)
{
    // Currently, the right-hand side only includes the internal force vector. In the future, it
    // will also include water pressure contributions and coupling terms.
    const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();
    const auto relative_displacements = CalculateRelativeDisplacementsAtIntegrationPoints(local_b_matrices);
    const auto tractions = CalculateTractionsAtIntegrationPoints(relative_displacements, rProcessInfo);
    const auto integration_coefficients = CalculateIntegrationCoefficients();
    rRightHandSideVector = -GeoEquationOfMotionUtilities::CalculateInternalForceVector(
        local_b_matrices, tractions, integration_coefficients);
}

void InterfaceElement::CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                                            VectorType&        rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void InterfaceElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                    std::vector<Vector>&    rOutput,
                                                    const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == STRAIN) {
        const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();
        rOutput = CalculateRelativeDisplacementsAtIntegrationPoints(local_b_matrices);
    } else if (rVariable == CAUCHY_STRESS_VECTOR) {
        const auto local_b_matrices = CalculateLocalBMatricesAtIntegrationPoints();
        const auto relative_displacements = CalculateRelativeDisplacementsAtIntegrationPoints(local_b_matrices);
        rOutput = CalculateTractionsAtIntegrationPoints(relative_displacements, rCurrentProcessInfo);
    } else {
        Element::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

void InterfaceElement::CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                                    std::vector<ConstitutiveLaw::Pointer>& rOutput,
                                                    const ProcessInfo&)
{
    KRATOS_ERROR_IF_NOT(rVariable == CONSTITUTIVE_LAW)
        << "Cannot calculate on integration points: got unexpected variable " << rVariable.Name() << "\n";

    rOutput = mConstitutiveLaws;
}

void InterfaceElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

void InterfaceElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    Element::Initialize(rCurrentProcessInfo);

    const auto shape_function_values_at_integration_points =
        GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(
            mIntegrationScheme->GetIntegrationPoints(), GetGeometry());

    mConstitutiveLaws.clear();
    for (auto i = std::size_t{0}; i < mIntegrationScheme->GetNumberOfIntegrationPoints(); ++i) {
        mConstitutiveLaws.push_back(GetProperties()[CONSTITUTIVE_LAW]->Clone());
    }
    // Only interpolate when neighbouring elements that provide nodal stresses were found
    if (this->Has(NEIGHBOUR_ELEMENTS)) {
        // interface element can at maximum have 2 neighbours
        KRATOS_DEBUG_ERROR_IF(this->GetValue(NEIGHBOUR_ELEMENTS).size() > 2)
            << "Too many neighbour elements for interface element " << this->Id() << std::endl;
        const auto interface_node_ids = GeometryUtilities::GetNodeIdsFromGeometry(GetGeometry());
        std::vector<std::optional<Vector>> interface_nodal_cauchy_stresses(interface_node_ids.size());
        for (auto& r_neighbour_element : this->GetValue(NEIGHBOUR_ELEMENTS)) {
            std::vector<Vector> neighbour_cauchy_stresses;
            r_neighbour_element.CalculateOnIntegrationPoints(
                CAUCHY_STRESS_VECTOR, neighbour_cauchy_stresses, rCurrentProcessInfo);
            interface_nodal_cauchy_stresses = ExtrapolationUtilities::CalculateNodalVectors(
                interface_node_ids, r_neighbour_element.GetGeometry(),
                r_neighbour_element.GetIntegrationMethod(), neighbour_cauchy_stresses,
                r_neighbour_element.Id());
            // values of previous iteration completely overwritten. that should probably change
        }
        InterpolateNodalStressesToIntegrationPointTractions(interface_nodal_cauchy_stresses);
    }
    for (auto i = std::size_t{0}; i < mConstitutiveLaws.size(); ++i) {
        mConstitutiveLaws[i]->InitializeMaterial(GetProperties(), GetGeometry(),
                                                 shape_function_values_at_integration_points[i]);
    }
}

int InterfaceElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    int error = Element::Check(rCurrentProcessInfo);
    if (error != 0) return error;

    if (this->IsActive()) {
        KRATOS_ERROR_IF(mIntegrationScheme->GetNumberOfIntegrationPoints() != mConstitutiveLaws.size())
            << "Number of integration points (" << mIntegrationScheme->GetNumberOfIntegrationPoints()
            << ") and constitutive laws (" << mConstitutiveLaws.size() << ") do not match.\n";

        const auto r_properties  = GetProperties();
        const auto expected_size = mpStressStatePolicy->GetVoigtSize();
        ConstitutiveLawUtilities::CheckStrainSize(r_properties, expected_size, Id());

        error = r_properties[CONSTITUTIVE_LAW]->Check(r_properties, GetGeometry(), rCurrentProcessInfo);
        return error;
    }

    return 0;
}

Element::DofsVectorType InterfaceElement::GetDofs() const
{
    const auto no_Pw_geometry_yet = Geometry<Node>{};
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(GetGeometry(), no_Pw_geometry_yet,
                                                      GetGeometry().WorkingSpaceDimension());
}

std::vector<Matrix> InterfaceElement::CalculateLocalBMatricesAtIntegrationPoints() const
{
    const auto& r_integration_points = mIntegrationScheme->GetIntegrationPoints();
    const auto  shape_function_values_at_integration_points =
        GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(r_integration_points, GetGeometry());

    auto result = std::vector<Matrix>{};
    result.reserve(shape_function_values_at_integration_points.size());
    auto calculate_local_b_matrix = [&r_geometry = GetGeometry(), this](const auto& rShapeFunctionValuesAtIntegrationPoint,
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

std::vector<double> InterfaceElement::CalculateIntegrationCoefficients() const
{
    const auto determinants_of_jacobian = CalculateDeterminantsOfJacobiansAtIntegrationPoints(
        mIntegrationScheme->GetIntegrationPoints(), GetGeometry());
    return mIntegrationCoefficientsCalculator.Run<>(mIntegrationScheme->GetIntegrationPoints(),
                                                    determinants_of_jacobian, this);
}

std::vector<Matrix> InterfaceElement::CalculateConstitutiveMatricesAtIntegrationPoints(
    const std::vector<Vector>& rRelativeDisplacements, const ProcessInfo& rProcessInfo)
{
    return ::CalculateConstitutiveMatricesAtIntegrationPoints(mConstitutiveLaws, GetProperties(),
                                                              rRelativeDisplacements, rProcessInfo);
}

std::vector<Vector> InterfaceElement::CalculateRelativeDisplacementsAtIntegrationPoints(const std::vector<Matrix>& rLocalBMatrices) const
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

std::vector<Vector> InterfaceElement::CalculateTractionsAtIntegrationPoints(const std::vector<Vector>& rRelativeDisplacements,
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

void InterfaceElement::ApplyRotationToBMatrix(Matrix& rBMatrix, const Matrix& rRotationMatrix) const
{
    const auto dim = GetGeometry().WorkingSpaceDimension();
    for (auto i = std::size_t{0}; i + dim <= rBMatrix.size2(); i += dim) {
        auto sub_matrix = subrange(rBMatrix, 0, rBMatrix.size1(), i, i + dim);
        sub_matrix.assign(Matrix{prod(sub_matrix, trans(rRotationMatrix))});
    }
}

void InterfaceElement::InterpolateNodalStressesToIntegrationPointTractions(
    const std::vector<std::optional<Vector>>& interface_nodal_cauchy_stresses) const
{
    // interpolate nodal stresses on a chosen side
    auto&      r_interface_geometry = GetGeometry();
    const auto number_of_nodes_on_side =
        GeometryUtilities::GetNodeIdsFromGeometry(r_interface_geometry).size() / 2;
    // check which side is shared with the neighbour geometry, by checking for the first nodes of the first side
    auto in_first_side = interface_nodal_cauchy_stresses[0].has_value();

    constexpr std::size_t start_index = in_first_side ? 0 : number_of_nodes_on_side;
    std::vector<Vector>   nodal_stresses;
    for (auto i = 0; i < number_of_nodes_on_side; ++i) {
        nodal_stresses.push_back(interface_nodal_cauchy_stresses[start_index + i].value());
    }
    KRATOS_DEBUG_ERROR_IF(nodal_stresses.size() < number_of_nodes_on_side)
        << "Less stresses found than needed for 1 side of the interface element with Id."
        << this->Id() << std::endl;

    std::size_t integration_point_index = 0;
    for (const auto& r_integration_point : mIntegrationScheme->GetIntegrationPoints()) {
        // interpolate on interface
        auto integration_point_stress                = Vector{ZeroVector(nodal_stresses[0].size())};
        auto integration_point_shape_function_values = Vector{};
        r_interface_geometry.ShapeFunctionsValues(integration_point_shape_function_values, r_integration_point);
        for (std::size_t i = 0; i < GetGeometry().PointsNumber() / 2; ++i) {
            integration_point_stress += integration_point_shape_function_values[i] * nodal_stresses[i];
        }

        auto rotation_matrix = Matrix(3, 3, 0.0);
        if (r_interface_geometry.LocalSpaceDimension() == 1) {
            auto two_d_rotation_matrix = mfpCalculateRotationMatrix(GetGeometry(), r_integration_point);
            GeoElementUtilities::AssembleUUBlockMatrix(rotation_matrix, two_d_rotation_matrix);
            rotation_matrix(2, 2) = 1.0;
        } else {
            rotation_matrix = mfpCalculateRotationMatrix(GetGeometry(), r_integration_point);
        }

        auto integration_point_stress_tensor = MathUtils<double>::StressVectorToTensor(integration_point_stress);
        auto integration_point_local_stress_tensor = GeoMechanicsMathUtilities::RotateSecondOrderTensor(
            integration_point_stress_tensor, rotation_matrix);

        // extract normal and shear components to form initial traction
        auto traction_vector = Vector(mConstitutiveLaws[0]->GetStrainSize());
        if (r_interface_geometry.LocalSpaceDimension() == 1) {
            traction_vector[0] = integration_point_local_stress_tensor(1, 1);
            traction_vector[1] = integration_point_local_stress_tensor(0, 1);
        } else {
            traction_vector[0] = integration_point_local_stress_tensor(2, 2);
            traction_vector[1] = integration_point_local_stress_tensor(0, 2);
            traction_vector[2] = integration_point_local_stress_tensor(1, 2);
        }
        const auto initial_state =
            make_intrusive<InitialState>(traction_vector, InitialState::InitialImposingType::STRESS_ONLY);
        mConstitutiveLaws[integration_point_index]->SetInitialState(initial_state);
        ++integration_point_index;
    }
}

// Instances of this class can not be copied but can be moved. Check that at compile time.
static_assert(!std::is_copy_constructible_v<InterfaceElement>);
static_assert(!std::is_copy_assignable_v<InterfaceElement>);
static_assert(std::is_move_constructible_v<InterfaceElement>);
static_assert(std::is_move_assignable_v<InterfaceElement>);
} // namespace Kratos
