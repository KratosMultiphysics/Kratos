// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//                   Based on the work of Massimo Petracca and Peter Wilson
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_elements/base_shell_element.h"
#include "custom_utilities/shell_utilities.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/shellt3_corotational_coordinate_transformation.hpp"
#include "custom_utilities/shellq4_corotational_coordinate_transformation.hpp"

namespace Kratos
{

using SizeType = std::size_t;
using IndexType = std::size_t;

template <class TCoordinateTransformation>
BaseShellElement<TCoordinateTransformation>::BaseShellElement(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry),
      mpCoordinateTransformation(Kratos::make_unique<TCoordinateTransformation>(pGeometry))
{
}

template <class TCoordinateTransformation>
BaseShellElement<TCoordinateTransformation>::BaseShellElement(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties),
      mpCoordinateTransformation(Kratos::make_unique<TCoordinateTransformation>(pGeometry))
{
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::EquationIdVector(EquationIdVectorType& rResult,
                                        const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType num_dofs = GetNumberOfDofs();

    if (rResult.size() != num_dofs) {
        rResult.resize(num_dofs, false);
    }

    const auto& r_geom = GetGeometry();

    for (IndexType i = 0; i < r_geom.size(); ++i) {
        const IndexType index = i * 6;
        const NodeType& i_node = r_geom[i];

        rResult[index]     = i_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = i_node.GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = i_node.GetDof(DISPLACEMENT_Z).EquationId();

        rResult[index + 3] = i_node.GetDof(ROTATION_X).EquationId();
        rResult[index + 4] = i_node.GetDof(ROTATION_Y).EquationId();
        rResult[index + 5] = i_node.GetDof(ROTATION_Z).EquationId();
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::GetDofList(DofsVectorType& rElementalDofList,
                                  const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType num_dofs = GetNumberOfDofs();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(num_dofs);

    const auto& r_geom = GetGeometry();

    for (IndexType i = 0; i < r_geom.size(); ++i) {
        const NodeType& i_node = r_geom[i];

        rElementalDofList.push_back(i_node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(i_node.pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(i_node.pGetDof(DISPLACEMENT_Z));

        rElementalDofList.push_back(i_node.pGetDof(ROTATION_X));
        rElementalDofList.push_back(i_node.pGetDof(ROTATION_Y));
        rElementalDofList.push_back(i_node.pGetDof(ROTATION_Z));
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::GetValuesVector(Vector& rValues, int Step) const
{
    const SizeType num_dofs = GetNumberOfDofs();

    if (rValues.size() != num_dofs) {
        rValues.resize(num_dofs, false);
    }

    const auto& r_geom = GetGeometry();

    for (IndexType i = 0; i < r_geom.size(); ++i) {
        const NodeType& i_node = r_geom[i];
        const array_1d<double, 3>& disp = i_node.FastGetSolutionStepValue(DISPLACEMENT, Step);
        const array_1d<double, 3>& rot = i_node.FastGetSolutionStepValue(ROTATION, Step);

        const IndexType index = i * 6;
        rValues[index]     = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];

        rValues[index + 3] = rot[0];
        rValues[index + 4] = rot[1];
        rValues[index + 5] = rot[2];
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    const SizeType num_dofs = GetNumberOfDofs();

    if (rValues.size() != num_dofs) {
        rValues.resize(num_dofs, false);
    }

    const auto& r_geom = GetGeometry();

    for (IndexType i = 0; i < r_geom.size(); ++i) {
        const NodeType& i_node = r_geom[i];
        const array_1d<double, 3>& r_vel = i_node.FastGetSolutionStepValue(VELOCITY, Step);
        const array_1d<double, 3>& r_ang_vel = i_node.FastGetSolutionStepValue(ANGULAR_VELOCITY, Step);

        const IndexType index = i * 6;
        rValues[index]     = r_vel[0];
        rValues[index + 1] = r_vel[1];
        rValues[index + 2] = r_vel[2];

        rValues[index + 3] = r_ang_vel[0];
        rValues[index + 4] = r_ang_vel[1];
        rValues[index + 5] = r_ang_vel[2];
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    const SizeType num_dofs = GetNumberOfDofs();

    if (rValues.size() != num_dofs) {
        rValues.resize(num_dofs, false);
    }

    const auto& r_geom = GetGeometry();

    for (IndexType i = 0; i < r_geom.size(); ++i) {
        const NodeType& i_node = r_geom[i];
        const array_1d<double, 3>& r_acc = i_node.FastGetSolutionStepValue(ACCELERATION, Step);
        const array_1d<double, 3>& r_ang_acc = i_node.FastGetSolutionStepValue(ANGULAR_ACCELERATION, Step);

        const IndexType index = i * 6;
        rValues[index]     = r_acc[0];
        rValues[index + 1] = r_acc[1];
        rValues[index + 2] = r_acc[2];

        rValues[index + 3] = r_ang_acc[0];
        rValues[index + 4] = r_ang_acc[1];
        rValues[index + 5] = r_ang_acc[2];
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::ResetConstitutiveLaw()
{
    KRATOS_TRY

    const auto& r_geom = GetGeometry();
    const Matrix& r_shape_fct_values = r_geom.ShapeFunctionsValues(GetIntegrationMethod());

    const auto& r_props = GetProperties();
    for (IndexType i = 0; i < mSections.size(); ++i) {
        mSections[i]->ResetCrossSection(r_props, r_geom, row(r_shape_fct_values, i));
    }

    KRATOS_CATCH("")
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        const auto& r_geom = GetGeometry();
        const auto& r_props = GetProperties();

        const SizeType num_gps = GetNumberOfGPs();

        if (mSections.size() != num_gps) {
            const Matrix& r_shape_fct_values =
                r_geom.ShapeFunctionsValues(GetIntegrationMethod());

            ShellCrossSection::Pointer p_ref_section;

            if (ShellUtilities::IsOrthotropic(r_props)) {
                // make new instance of shell cross section
                p_ref_section = Kratos::make_shared<ShellCrossSection>();

                // Parse material properties for each layer
                p_ref_section->ParseOrthotropicPropertyMatrix(r_props);
            } else {
                p_ref_section = Kratos::make_shared<ShellCrossSection>();
                const IndexType ply_index = 0;
                const SizeType num_points = 5;
                p_ref_section->BeginStack();
                p_ref_section->AddPly(ply_index, num_points, r_props);
                p_ref_section->EndStack();
            }

            mSections.clear();
            for (SizeType i = 0; i < num_gps; ++i) {
                ShellCrossSection::Pointer p_section_clone = p_ref_section->Clone();
                p_section_clone->SetSectionBehavior(GetSectionBehavior());
                p_section_clone->InitializeCrossSection(r_props, r_geom, row(r_shape_fct_values, i));
                mSections.push_back(p_section_clone);
            }
        }

        if (Has(LOCAL_MATERIAL_AXIS_1)) {
            // calculate the angle between the prescribed direction and the local axis 1
            // this is currently required in teh derived classes TODO refactor

            std::vector<array_1d<double, 3>> local_axes_1;
            std::vector<array_1d<double, 3>> local_axes_2;
            CalculateOnIntegrationPoints(LOCAL_AXIS_1, local_axes_1 , rCurrentProcessInfo);
            CalculateOnIntegrationPoints(LOCAL_AXIS_2, local_axes_2 , rCurrentProcessInfo);

            const array_1d<double, 3> prescribed_direcition = GetValue(LOCAL_MATERIAL_AXIS_1);

            double mat_orientation_angle = MathUtils<double>::VectorsAngle(local_axes_1[0], prescribed_direcition);

            // make sure the angle is positively defined according to right hand rule
            if (inner_prod(local_axes_2[0], prescribed_direcition) < 0.0) {
                // mat_orientation_angle is currently negative, flip to positive definition
                mat_orientation_angle *= -1.0;
            }

            SetValue(MATERIAL_ORIENTATION_ANGLE, mat_orientation_angle);
        }

        mpCoordinateTransformation->Initialize();
        SetupOrientationAngles();
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    mpCoordinateTransformation->InitializeNonLinearIteration();

    const auto& r_geom = GetGeometry();
    const Matrix& r_shape_fct_values = r_geom.ShapeFunctionsValues(GetIntegrationMethod());
    for (IndexType i = 0; i < mSections.size(); ++i) {
        mSections[i]->InitializeNonLinearIteration(GetProperties(), r_geom, row(r_shape_fct_values, i), rCurrentProcessInfo);
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    mpCoordinateTransformation->FinalizeNonLinearIteration();

    const auto& r_geom = GetGeometry();
    const Matrix& r_shape_fct_values = r_geom.ShapeFunctionsValues(GetIntegrationMethod());
    for (IndexType i = 0; i < mSections.size(); ++i) {
        mSections[i]->FinalizeNonLinearIteration(GetProperties(), r_geom, row(r_shape_fct_values, i), rCurrentProcessInfo);
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_props = GetProperties();
    const auto& r_geom = GetGeometry();
    const Matrix& r_shape_fct_values = r_geom.ShapeFunctionsValues(GetIntegrationMethod());

    for (IndexType i = 0; i < mSections.size(); ++i) {
        mSections[i]->InitializeSolutionStep(r_props, r_geom, row(r_shape_fct_values, i), rCurrentProcessInfo);
    }

    mpCoordinateTransformation->InitializeSolutionStep();
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_props = GetProperties();
    const auto& r_geom = GetGeometry();
    const Matrix& r_shape_fct_values = r_geom.ShapeFunctionsValues(GetIntegrationMethod());

    for (IndexType i = 0; i < mSections.size(); ++i) {
        mSections[i]->FinalizeSolutionStep(r_props, r_geom, row(r_shape_fct_values, i), rCurrentProcessInfo);
    }

    mpCoordinateTransformation->FinalizeSolutionStep();
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
{
    // Calculation flags
    const bool calculate_stiffness_matrix_flag = true;
    const bool calculate_residual_vector_flag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                 calculate_stiffness_matrix_flag, calculate_residual_vector_flag);
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo)
{
    // Calculation flags
    const bool calculate_stiffness_matrix_flag = true;
    const bool calculate_residual_vector_flag = true; // TODO check is this can be false => see solids

    Vector dummy;
    CalculateAll(rLeftHandSideMatrix, dummy, rCurrentProcessInfo,
                 calculate_stiffness_matrix_flag, calculate_residual_vector_flag);
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::CalculateRightHandSide(VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
{
    // Calculation flags
    const bool calculate_stiffness_matrix_flag = true; // TODO check is this can be false => see solids
    const bool calculate_residual_vector_flag = true;

    Matrix dummy;
    CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo,
                 calculate_stiffness_matrix_flag, calculate_residual_vector_flag);
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    const bool compute_lumped_mass_matrix = StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(GetProperties(), rCurrentProcessInfo);

    const SizeType num_gps = GetNumberOfGPs();
    const SizeType num_dofs = GetNumberOfDofs();
    const SizeType num_nodes = GetGeometry().PointsNumber();

    if ((rMassMatrix.size1() != num_dofs) || (rMassMatrix.size2() != num_dofs)) {
        rMassMatrix.resize(num_dofs, num_dofs, false);
    }
    noalias(rMassMatrix) = ZeroMatrix(num_dofs, num_dofs);

    // Compute the local coordinate system.
    auto ref_coord_sys(mpCoordinateTransformation->CreateReferenceCoordinateSystem());
    const double ref_area = ref_coord_sys.Area();

    // Average mass per unit area over the whole element
    double av_mass_per_unit_area = 0.0;

    for (SizeType i = 0; i < num_gps; i++) {
        av_mass_per_unit_area += this->mSections[i]->CalculateMassPerUnitArea(GetProperties());
    }
    av_mass_per_unit_area /= static_cast<double>(num_gps);

    if (compute_lumped_mass_matrix) {
        // lumped area
        double lump_area = ref_area / static_cast<double>(num_nodes);
        double nodal_mass = av_mass_per_unit_area * lump_area;

        // loop on nodes
        for (SizeType i=0; i < num_nodes; i++) {
            SizeType index = i * 6;

            // translational mass
            rMassMatrix(index, index) = nodal_mass;
            rMassMatrix(index + 1, index + 1) = nodal_mass;
            rMassMatrix(index + 2, index + 2) = nodal_mass;

            // rotational mass - neglected for the moment...
        }
    } else { // consistent mass matrix
        if (num_nodes == 3) { // triangular element
            // General matrix form as per Felippa plane stress CST eqn 31.27:
            // http://kis.tu.kielce.pl/mo/COLORADO_FEM/colorado/IFEM.Ch31.pdf

            // Density and thickness are averaged over element.

            // Average thickness over the whole element
            double thickness = 0.0;
            for (SizeType i = 0; i < num_gps; i++) {
                thickness += this->mSections[i]->GetThickness(GetProperties());
            }
            thickness /= static_cast<double>(num_gps);

            // Populate mass matrix with integation results
            for (SizeType row = 0; row < num_dofs; row++) {
                if (row % 6 < 3) { // translational entry
                    for (SizeType col = 0; col < 3; col++) {
                        rMassMatrix(row, 6 * col + row % 6) = 1.0;
                    }
                } else { // rotational entry
                    for (SizeType col = 0; col < 3; col++) {
                        rMassMatrix(row, 6 * col + row % 6) = thickness*thickness / 12.0;
                    }
                }

                // Diagonal entry
                rMassMatrix(row, row) *= 2.0;
            }

            rMassMatrix *= av_mass_per_unit_area * ref_area / 12.0;

        } else { // quadrilateral element
            // Get shape function values and setup jacobian
            const GeometryType& geom = GetGeometry();
            const Matrix& shapeFunctions = geom.ShapeFunctionsValues();
            ShellUtilities::JacobianOperator jac_operator;

            // Get integration points
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->mIntegrationMethod);

            // Setup matrix of shape functions
            Matrix N = Matrix(6, 24, 0.0);

            // Gauss loop
            for (SizeType gauss_point = 0; gauss_point < 4; gauss_point++) {
                // Calculate average mass per unit area and thickness at the
                // current GP
                av_mass_per_unit_area =
                    this->mSections[gauss_point]->CalculateMassPerUnitArea(GetProperties());
                const double thickness = this->mSections[gauss_point]->GetThickness(GetProperties());

                // Calc jacobian and weighted dA at current GP
                jac_operator.Calculate(ref_coord_sys, geom.ShapeFunctionLocalGradient(gauss_point));
                const double dA = integration_points[gauss_point].Weight() *
                    jac_operator.Determinant();

                // Assemble shape function matrix over nodes
                for (SizeType node = 0; node < 4; node++) {
                    // translational entries - dofs 1-3
                    for (SizeType dof = 0; dof < 3; dof++) {
                        N(dof, 6 * node + dof) =
                            shapeFunctions(gauss_point, node);
                    }

                    // rotational inertia entries - dofs 4-6
                    for (SizeType dof = 0; dof < 3; dof++) {
                        N(dof + 3, 6 * node + dof + 3) =
                            thickness / std::sqrt(12.0) *
                            shapeFunctions(gauss_point, node);
                    }
                }

                // Add contribution to total mass matrix
                rMassMatrix += prod(trans(N), N)*dA*av_mass_per_unit_area;
            }
        }
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    const std::size_t matrix_size = GetNumberOfDofs();

    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        matrix_size);
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3> >& rVariable,
    std::vector<array_1d<double, 3> >& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == LOCAL_AXIS_1 ||
            rVariable == LOCAL_AXIS_2 ||
            rVariable == LOCAL_AXIS_3) {
        ComputeLocalAxis(rVariable, rOutput);
    } else if (rVariable == LOCAL_MATERIAL_AXIS_1 ||
               rVariable == LOCAL_MATERIAL_AXIS_2 ||
               rVariable == LOCAL_MATERIAL_AXIS_3) {
        ComputeLocalMaterialAxis(rVariable, rOutput);
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::CalculateOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        rValues.clear();
        for (const auto& p_sec : mSections) {
            auto vec_integration_points = p_sec->GetConstitutiveLawsVector(GetProperties());
            rValues.reserve(rValues.size() + vec_integration_points.size());
            for (std::size_t i=0; i<vec_integration_points.size(); ++i) {
                rValues.push_back(vec_integration_points[i]);
            }
        }
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::Calculate(
    const Variable<Matrix>& rVariable, Matrix& Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == LOCAL_ELEMENT_ORIENTATION) {
        Output.resize(3, 3, false);

        // Compute the local coordinate system.
        auto localCoordinateSystem(mpCoordinateTransformation->CreateReferenceCoordinateSystem());
        Output = trans(localCoordinateSystem.Orientation());
    }
}

template <class TCoordinateTransformation>
int BaseShellElement<TCoordinateTransformation>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    Element::Check(rCurrentProcessInfo);

    CheckDofs();
    CheckProperties(rCurrentProcessInfo);

    KRATOS_ERROR_IF(GetGeometry().Area() < std::numeric_limits<double>::epsilon()*1000)
            << "Element #" << Id() << " has an Area of zero!" << std::endl;

    // TODO check ConstLaws

    return 0;

    KRATOS_CATCH("")
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::SetCrossSectionsOnIntegrationPoints(std::vector< ShellCrossSection::Pointer >& crossSections)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(crossSections.size() == GetNumberOfGPs())
            << "The number of cross section is wrong: " << crossSections.size() << std::endl;
    mSections.clear();
    for (IndexType i = 0; i < crossSections.size(); ++i) {
        mSections.push_back(crossSections[i]);
    }
    SetupOrientationAngles();
    KRATOS_CATCH("")
}

template <class TCoordinateTransformation>
const Parameters BaseShellElement<TCoordinateTransformation>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["static","implicit","explicit"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : [],
            "nodal_historical"       : ["DISPLACEMENT","ROTATION","VELOCITY","ACCELERATION"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISPLACEMENT","ROTATION"],
        "required_dofs"              : ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z","ROTATION_X","ROTATION_Y","ROTATION_Z"],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle3D3", "Quadrilateral3D4"],
        "element_integrates_in_time" : false,
        "compatible_constitutive_laws": {
            "type"        : ["PlaneStress"],
            "dimension"   : ["3D"],
            "strain_size" : [3]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   : "Base element for thin/thick shell formulations."
    })");

    return specifications;
}

template <class TCoordinateTransformation>
std::string BaseShellElement<TCoordinateTransformation>::Info() const
{
    std::stringstream buffer;
    buffer << "BaseShellElement #" << Id();
    return buffer.str();
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "BaseShellElement #" << Id();
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::PrintData(std::ostream& rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

template <class TCoordinateTransformation>
SizeType BaseShellElement<TCoordinateTransformation>::GetNumberOfDofs() const
{
    return (6 * GetGeometry().PointsNumber());   // 6 dofs per node
}

template <class TCoordinateTransformation>
SizeType BaseShellElement<TCoordinateTransformation>::GetNumberOfGPs() const
{
    return GetGeometry().IntegrationPoints(mIntegrationMethod).size();
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
)
{
    KRATOS_ERROR << "You have called to the CalculateAll from the base class for shell elements" << std::endl;
}

template <class TCoordinateTransformation>
ShellCrossSection::SectionBehaviorType BaseShellElement<TCoordinateTransformation>::GetSectionBehavior() const
{
    KRATOS_ERROR << "You have called to the GetSectionBehavior from the base class for shell elements" << std::endl;
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::SetupOrientationAngles()
{
    if (Has(MATERIAL_ORIENTATION_ANGLE)) {
        for (auto it = mSections.begin(); it != mSections.end(); ++it) {
            (*it)->SetOrientationAngle(GetValue(MATERIAL_ORIENTATION_ANGLE));
        }
    } else {
        auto lcs(mpCoordinateTransformation->CreateReferenceCoordinateSystem());

        Vector3Type normal;
        noalias(normal) = lcs.Vz();

        Vector3Type dZ;
        dZ(0) = 0.0;
        dZ(1) = 0.0;
        dZ(2) = 1.0; // for the moment let's take this. But the user can specify its own triad! TODO

        Vector3Type dirX;
        MathUtils<double>::CrossProduct(dirX, dZ, normal);

        // try to normalize the x vector. if it is near zero it means that we need
        // to choose a default one.
        double dirX_norm = dirX(0)*dirX(0) + dirX(1)*dirX(1) + dirX(2)*dirX(2);
        if (dirX_norm < 1.0E-12) {
            dirX(0) = 1.0;
            dirX(1) = 0.0;
            dirX(2) = 0.0;
        } else if (dirX_norm != 1.0) {
            dirX_norm = std::sqrt(dirX_norm);
            dirX /= dirX_norm;
        }

        Vector3Type elem_dirX = lcs.Vx();

        // now calculate the angle between the element x direction and the material x direction.
        Vector3Type& a = elem_dirX;
        Vector3Type& b = dirX;
        double a_dot_b = a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
        if (a_dot_b < -1.0) {
            a_dot_b = -1.0;
        }
        if (a_dot_b > 1.0) {
            a_dot_b = 1.0;
        }
        double angle = std::acos(a_dot_b);

        // if they are not counter-clock-wise, let's change the sign of the angle
        if (angle != 0.0) {
            const MatrixType& R = lcs.Orientation();
            if (dirX(0)*R(1, 0) + dirX(1)*R(1, 1) + dirX(2)*R(1, 2) < 0.0) {
                angle = -angle;
            }
        }

        for (auto it = mSections.begin(); it != mSections.end(); ++it) {
            (*it)->SetOrientationAngle(angle);
        }
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::ComputeLocalAxis(
    const Variable<array_1d<double, 3> >& rVariable,
    std::vector<array_1d<double, 3> >& rOutput) const
{
    const SizeType num_gps = GetNumberOfGPs();
    if (rOutput.size() != num_gps) {
        rOutput.resize(num_gps);
    }

    for (IndexType i=1; i<num_gps; ++i) {
        noalias(rOutput[i]) = ZeroVector(3);
    }

    const auto localCoordinateSystem(mpCoordinateTransformation->CreateLocalCoordinateSystem());
    if (rVariable == LOCAL_AXIS_1) {
        noalias(rOutput[0]) = localCoordinateSystem.Vx();
    } else if (rVariable == LOCAL_AXIS_2) {
        noalias(rOutput[0]) = localCoordinateSystem.Vy();
    } else if (rVariable == LOCAL_AXIS_3) {
        noalias(rOutput[0]) = localCoordinateSystem.Vz();
    } else {
        KRATOS_ERROR << "Wrong variable: " << rVariable.Name() << "!" << std::endl;
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::ComputeLocalMaterialAxis(
    const Variable<array_1d<double, 3> >& rVariable,
    std::vector<array_1d<double, 3> >& rOutput) const
{
    const double mat_angle = GetValue(MATERIAL_ORIENTATION_ANGLE);

    const SizeType num_gps = GetNumberOfGPs();
    if (rOutput.size() != num_gps) {
        rOutput.resize(num_gps);
    }

    for (IndexType i=1; i<num_gps; ++i) {
        noalias(rOutput[i]) = ZeroVector(3);
    }

    const auto localCoordinateSystem(mpCoordinateTransformation->CreateLocalCoordinateSystem());

    const auto eZ = localCoordinateSystem.Vz();

    if (rVariable == LOCAL_MATERIAL_AXIS_1) {
        const auto q = QuaternionType::FromAxisAngle(eZ(0), eZ(1), eZ(2), mat_angle);
        q.RotateVector3(localCoordinateSystem.Vx(), rOutput[0]);
    } else if (rVariable == LOCAL_MATERIAL_AXIS_2) {
        const auto q = QuaternionType::FromAxisAngle(eZ(0), eZ(1), eZ(2), mat_angle);
        q.RotateVector3(localCoordinateSystem.Vy(), rOutput[0]);
    } else if (rVariable == LOCAL_MATERIAL_AXIS_3) {
        noalias(rOutput[0]) = eZ;
    } else {
        KRATOS_ERROR << "Wrong variable: " << rVariable.Name() << "!" << std::endl;
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::CheckDofs() const
{
    // verify that the dofs exist
    for (const auto& r_node : GetGeometry().Points()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ROTATION_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ROTATION_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ROTATION_Z, r_node);

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node);

        KRATOS_ERROR_IF(r_node.GetBufferSize() < 2) << "This Element needs "
                << "at least a buffer size = 2" << std::endl;
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::CheckProperties(const ProcessInfo& rCurrentProcessInfo) const
{
    // check properties
    if (pGetProperties() == nullptr) {
        KRATOS_ERROR << "Properties not provided for element " << Id() << std::endl;
    }

    const auto& r_props = GetProperties();

    const auto& r_geom = GetGeometry(); // TODO check if this can be const

    if (r_props.Has(SHELL_ORTHOTROPIC_LAYERS)) {
        CheckSpecificProperties();

        const auto& r_props = GetProperties();

        KRATOS_ERROR_IF(r_props.Has(THICKNESS)) << "Specifying THICKNESS conflicts with the "
                                                << "definition of SHELL_ORTHOTROPIC_LAYERS (where the thickness is also specified)"
                                                << std::endl;

        KRATOS_ERROR_IF(r_props.Has(DENSITY)) << "Specifying DENSITY conflicts with the "
                                              << "definition of SHELL_ORTHOTROPIC_LAYERS (where the density is also specified)"
                                              << std::endl;

        KRATOS_ERROR_IF(r_props.Has(YOUNG_MODULUS)) << "Specifying YOUNG_MODULUS conflicts with the "
                << "definition of SHELL_ORTHOTROPIC_LAYERS (where the youngs-modulus is also specified)"
                << std::endl;

        KRATOS_ERROR_IF(r_props.Has(POISSON_RATIO)) << "Specifying POISSON_RATIO conflicts with the "
                << "definition of SHELL_ORTHOTROPIC_LAYERS (where the poisson-ratio is also specified)"
                << std::endl;

        // perform detailed orthotropic check later in shell_cross_section
    } else { // ... allow the automatic creation of a homogeneous section from a material and a thickness
        CheckSpecificProperties();

        const auto& r_props = GetProperties();

        if (!r_props.Has(THICKNESS)) {
            KRATOS_ERROR << "THICKNESS not provided for element " << Id() << std::endl;
        }
        if (r_props[THICKNESS] <= 0.0) {
            KRATOS_ERROR << "wrong THICKNESS value provided for element " << Id() << std::endl;
        }

        if (!r_props.Has(DENSITY)) {
            KRATOS_ERROR << "DENSITY not provided for element " << Id() << std::endl;
        }
        if (r_props[DENSITY] < 0.0) {
            KRATOS_ERROR << "wrong DENSITY value provided for element " << Id() << std::endl;
        }

        // TODO is this needed???? => it is, the dummy is needed for "Check" => unify!
        ShellCrossSection::Pointer dummySection = ShellCrossSection::Pointer(new ShellCrossSection());
        dummySection->BeginStack();
        dummySection->AddPly(0, 5, GetProperties());
        dummySection->EndStack();
        dummySection->SetSectionBehavior(ShellCrossSection::Thick);
        dummySection->Check(r_props, r_geom, rCurrentProcessInfo);
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::CheckSpecificProperties() const
{
    const auto& r_props = GetProperties();

    if (!r_props.Has(CONSTITUTIVE_LAW)) {
        KRATOS_ERROR << "CONSTITUTIVE_LAW not provided for element " << Id() << std::endl;
    }
    const ConstitutiveLaw::Pointer& claw = r_props[CONSTITUTIVE_LAW];
    if (claw == nullptr) {
        KRATOS_ERROR << "CONSTITUTIVE_LAW not provided for element " << Id() << std::endl;
    }

    ConstitutiveLaw::Features LawFeatures;
    claw->GetLawFeatures(LawFeatures);

    // KRATOS_ERROR_IF(LawFeatures.mOptions.Is(ConstitutiveLaw::ANISOTROPIC) &&
    //                 !Has(MATERIAL_ORIENTATION_ANGLE))
    //     << "Using an Anisotropic Constitutive law requires the specification of "
    //     << "\"MATERIAL_ORIENTATION_ANGLE\" for shell element with Id " << Id() << std::endl;

    if (GetSectionBehavior() == ShellCrossSection::Thick) {
        // Check constitutive law has been verified with Stenberg stabilization
        // applicable for 5-parameter shells only.
        bool stenberg_stabilization_suitable = false;
        claw->GetValue(STENBERG_SHEAR_STABILIZATION_SUITABLE, stenberg_stabilization_suitable);
        KRATOS_WARNING_IF("BaseShellElement", !stenberg_stabilization_suitable)
                << "The current constitutive law has not been checked with Stenberg "
                << "shear stabilization.\nPlease check results carefully." << std::endl;
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::DecimalCorrection(Vector& a)
{
    const double norm = norm_2(a);
    const double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
    for (SizeType i = 0; i < a.size(); i++) {
        if (std::abs(a(i)) < tolerance) {
            a(i) = 0.0;
        }
    }
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    rSerializer.save("Sections", mSections);
    rSerializer.save("CoordinateTransformation", mpCoordinateTransformation);
    rSerializer.save("IntM", (int)mIntegrationMethod);
}

template <class TCoordinateTransformation>
void BaseShellElement<TCoordinateTransformation>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    rSerializer.load("Sections", mSections);
    rSerializer.load("CoordinateTransformation", mpCoordinateTransformation);
    int temp;
    rSerializer.load("IntM", temp);
    mIntegrationMethod = (IntegrationMethod)temp;
}

template class BaseShellElement< ShellT3_CoordinateTransformation >;
template class BaseShellElement< ShellT3_CorotationalCoordinateTransformation  >;
template class BaseShellElement< ShellQ4_CoordinateTransformation >;
template class BaseShellElement< ShellQ4_CorotationalCoordinateTransformation  >;

} // namespace Kratos.
