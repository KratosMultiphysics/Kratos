//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Based on the work of Massimo Petracca and Peter Wilson
//

// System includes

// External includes

// Project includes
#include "custom_elements/base_shell_element.h"
#include "custom_utilities/shell_utilities.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "includes/checks.h"

namespace Kratos
{

using SizeType = std::size_t;
using IndexType = std::size_t;

BaseShellElement::BaseShellElement(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

BaseShellElement::BaseShellElement(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

BaseShellElement::~BaseShellElement()
{
}

void BaseShellElement::EquationIdVector(EquationIdVectorType& rResult,
                                        ProcessInfo& rCurrentProcessInfo)
{
    const SizeType num_dofs = GetNumberOfDofs();

    if (rResult.size() != num_dofs) {
        rResult.resize(num_dofs, false);
    }

    auto& r_geom = GetGeometry();

    for (IndexType i = 0; i < r_geom.size(); ++i) {
        const IndexType index = i * 6;
        NodeType& i_node = r_geom[i];

        rResult[index]     = i_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = i_node.GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = i_node.GetDof(DISPLACEMENT_Z).EquationId();

        rResult[index + 3] = i_node.GetDof(ROTATION_X).EquationId();
        rResult[index + 4] = i_node.GetDof(ROTATION_Y).EquationId();
        rResult[index + 5] = i_node.GetDof(ROTATION_Z).EquationId();
    }
}

void BaseShellElement::GetDofList(DofsVectorType& rElementalDofList,
                                  ProcessInfo& rCurrentProcessInfo)
{
    const SizeType num_dofs = GetNumberOfDofs();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(num_dofs);

    auto& r_geom = GetGeometry();

    for (IndexType i = 0; i < r_geom.size(); ++i) {
        NodeType& i_node = r_geom[i];

        rElementalDofList.push_back(i_node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(i_node.pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(i_node.pGetDof(DISPLACEMENT_Z));

        rElementalDofList.push_back(i_node.pGetDof(ROTATION_X));
        rElementalDofList.push_back(i_node.pGetDof(ROTATION_Y));
        rElementalDofList.push_back(i_node.pGetDof(ROTATION_Z));
    }
}

void BaseShellElement::GetValuesVector(Vector& rValues, int Step)
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

void BaseShellElement::GetFirstDerivativesVector(Vector& rValues, int Step)
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

void BaseShellElement::GetSecondDerivativesVector(Vector& rValues, int Step)
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

void BaseShellElement::ResetConstitutiveLaw()
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

void BaseShellElement::Initialize()
{
    const auto& r_geom = GetGeometry();
    const auto& r_props = GetProperties();

    const SizeType num_gps = GetNumberOfGPs();

    if (mSections.size() != num_gps) {
        const Matrix& r_shape_fct_values =
            r_geom.ShapeFunctionsValues(GetIntegrationMethod());

        ShellCrossSection::Pointer p_ref_section;

        if (r_props.Has(SHELL_CROSS_SECTION)) {
            p_ref_section = r_props[SHELL_CROSS_SECTION];
        } else if (ShellUtilities::IsOrthotropic(r_props)) {
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

    if (this->Has(LOCAL_MATERIAL_AXIS_1)) {
        // calculate the angle between the prescribed direction and the local axis 1
        // this is currently required in teh derived classes TODO refactor

        std::vector<array_1d<double, 3>> local_axes_1;
		std::vector<array_1d<double, 3>> local_axes_2;
        ProcessInfo tmp_process_info; // TODO refactor once Initialize gets ProcessInfo
		this->CalculateOnIntegrationPoints(LOCAL_AXIS_1, local_axes_1 , tmp_process_info);
		this->CalculateOnIntegrationPoints(LOCAL_AXIS_2, local_axes_2 , tmp_process_info);

        const array_1d<double, 3> prescribed_direcition = this->GetValue(LOCAL_MATERIAL_AXIS_1);

        double mat_orientation_angle = MathUtils<double>::VectorsAngle(local_axes_1[0], prescribed_direcition);

        // make sure the angle is positively defined according to right hand rule
        if (inner_prod(local_axes_2[0], prescribed_direcition) < 0.0) {
            // mat_orientation_angle is currently negative, flip to positive definition
            mat_orientation_angle *= -1.0;
        }

        this->SetValue(MATERIAL_ORIENTATION_ANGLE, mat_orientation_angle);
    }
}

void BaseShellElement::BaseInitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geom = this->GetGeometry();
    const Matrix& r_shape_fct_values = r_geom.ShapeFunctionsValues(GetIntegrationMethod());
    for (IndexType i = 0; i < mSections.size(); ++i) {
        mSections[i]->InitializeNonLinearIteration(GetProperties(), r_geom, row(r_shape_fct_values, i), rCurrentProcessInfo);
    }
}

void BaseShellElement::BaseFinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geom = this->GetGeometry();
    const Matrix& r_shape_fct_values = r_geom.ShapeFunctionsValues(GetIntegrationMethod());
    for (IndexType i = 0; i < mSections.size(); ++i) {
        mSections[i]->FinalizeNonLinearIteration(GetProperties(), r_geom, row(r_shape_fct_values, i), rCurrentProcessInfo);
    }
}

void BaseShellElement::BaseInitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_props = GetProperties();
    const auto& r_geom = GetGeometry();
    const Matrix& r_shape_fct_values = r_geom.ShapeFunctionsValues(GetIntegrationMethod());

    for (IndexType i = 0; i < mSections.size(); ++i) {
        mSections[i]->InitializeSolutionStep(r_props, r_geom, row(r_shape_fct_values, i), rCurrentProcessInfo);
    }
}

void BaseShellElement::BaseFinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_props = GetProperties();
    const auto& r_geom = GetGeometry();
    const Matrix& r_shape_fct_values = r_geom.ShapeFunctionsValues(GetIntegrationMethod());

    for (IndexType i = 0; i < mSections.size(); ++i) {
        mSections[i]->FinalizeSolutionStep(r_props, r_geom, row(r_shape_fct_values, i), rCurrentProcessInfo);
    }
}

void BaseShellElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
    // Calculation flags
    const bool calculate_stiffness_matrix_flag = true;
    const bool calculate_residual_vector_flag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                 calculate_stiffness_matrix_flag, calculate_residual_vector_flag);
}


void BaseShellElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo)
{
    // Calculation flags
    const bool calculate_stiffness_matrix_flag = true;
    const bool calculate_residual_vector_flag = true; // TODO check is this can be false => see solids

    Vector dummy;
    CalculateAll(rLeftHandSideMatrix, dummy, rCurrentProcessInfo,
                 calculate_stiffness_matrix_flag, calculate_residual_vector_flag);
}

void BaseShellElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
    // Calculation flags
    const bool calculate_stiffness_matrix_flag = true; // TODO check is this can be false => see solids
    const bool calculate_residual_vector_flag = true;

    Matrix dummy;
    CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo,
                 calculate_stiffness_matrix_flag, calculate_residual_vector_flag);
}

void BaseShellElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    // TODO unify implementation and move it to BaseClass
}

void BaseShellElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
)
{
    const std::size_t matrix_size = GetNumberOfDofs();

    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        matrix_size);
}

int BaseShellElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    Element::Check(rCurrentProcessInfo);

    CheckVariables();
    CheckDofs();
    CheckProperties(rCurrentProcessInfo);

    KRATOS_ERROR_IF(GetGeometry().Area() < std::numeric_limits<double>::epsilon()*1000)
            << "Element #" << Id() << " has an Area of zero!" << std::endl;

    // TODO check ConstLaws

    return 0;

    KRATOS_CATCH("")
}

void BaseShellElement::SetCrossSectionsOnIntegrationPoints(std::vector< ShellCrossSection::Pointer >& crossSections)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(crossSections.size() == GetNumberOfGPs())
            << "The number of cross section is wrong: " << crossSections.size() << std::endl;
    mSections.clear();
    for (IndexType i = 0; i < crossSections.size(); ++i) {
        mSections.push_back(crossSections[i]);
    }
    this->SetupOrientationAngles();
    KRATOS_CATCH("")
}

std::string BaseShellElement::Info() const
{
    std::stringstream buffer;
    buffer << "BaseShellElement #" << Id();
    return buffer.str();
}

void BaseShellElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "BaseShellElement #" << Id();
}

void BaseShellElement::PrintData(std::ostream& rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

SizeType BaseShellElement::GetNumberOfDofs() const
{
    return (6 * GetGeometry().PointsNumber());   // 6 dofs per node
}

SizeType BaseShellElement::GetNumberOfGPs() const
{
    return GetGeometry().IntegrationPoints(mIntegrationMethod).size();
}

void BaseShellElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
)
{
    KRATOS_ERROR << "You have called to the CalculateAll from the base class for shell elements" << std::endl;
}

ShellCrossSection::SectionBehaviorType BaseShellElement::GetSectionBehavior() const
{
    KRATOS_ERROR << "You have called to the GetSectionBehavior from the base class for shell elements" << std::endl;
}

void BaseShellElement::SetupOrientationAngles()
{
    KRATOS_ERROR << "You have called to the SetupOrientationAngles from the base class for shell elements" << std::endl;
}

void BaseShellElement::CheckVariables() const
{
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ROTATION);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ANGULAR_VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(ANGULAR_ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS);
    KRATOS_CHECK_VARIABLE_KEY(SHELL_ORTHOTROPIC_LAYERS);
    KRATOS_CHECK_VARIABLE_KEY(CONSTITUTIVE_LAW);
    KRATOS_CHECK_VARIABLE_KEY(SHELL_CROSS_SECTION);
}

void BaseShellElement::CheckDofs() const
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

void BaseShellElement::CheckProperties(const ProcessInfo& rCurrentProcessInfo) const
{
    // check properties
    if (pGetProperties() == nullptr) {
        KRATOS_ERROR << "Properties not provided for element " << Id() << std::endl;
    }

    const auto& r_props = GetProperties();

    const auto& r_geom = GetGeometry(); // TODO check if this can be const

    if (r_props.Has(SHELL_CROSS_SECTION)) { // if the user specified a cross section ...
        const ShellCrossSection::Pointer& section = r_props[SHELL_CROSS_SECTION];
        if (section == nullptr) {
            KRATOS_ERROR << "SHELL_CROSS_SECTION not provided for element " << Id() << std::endl;
        }

        section->Check(r_props, r_geom, rCurrentProcessInfo);
    } else if (r_props.Has(SHELL_ORTHOTROPIC_LAYERS)) {
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

void BaseShellElement::CheckSpecificProperties() const
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
    //     << "\"MATERIAL_ORIENTATION_ANGLE\" for shell element with Id " << this->Id() << std::endl;

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

void BaseShellElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    rSerializer.save("Sections", mSections);
    rSerializer.save("IntM", (int)mIntegrationMethod);
}

void BaseShellElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    rSerializer.load("Sections", mSections);
    int temp;
    rSerializer.load("IntM", temp);
    mIntegrationMethod = (IntegrationMethod)temp;
}

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, BaseShellElement& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const BaseShellElement& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.
