// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"

// Application includes
#include "structural_mechanics_application_variables.h"

// Include base h
#include "custom_elements/bushing_element.h"

namespace Kratos {

BushingElement::BushingElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
}

BushingElement::BushingElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

BushingElement::BushingElement(BushingElement const& rOther)
    : Element(rOther)
{
}

BushingElement& BushingElement::operator=(BushingElement const& rOther)
{
    Element::operator=(rOther);
    return *this;
}

Element::Pointer BushingElement::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties) const
{
    // NEEDED TO CREATE AN ELEMENT
    return Kratos::make_intrusive<BushingElement>(
        NewId, GetGeometry().Create(rThisNodes), pProperties);
}

Element::Pointer BushingElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    // NEEDED TO CREATE AN ELEMENT
    return Kratos::make_intrusive<BushingElement>(NewId, pGeom, pProperties);
}

Element::Pointer BushingElement::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    // YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
    // ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

    BushingElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

    return Kratos::make_intrusive<BushingElement>(NewElement);
}

BushingElement::~BushingElement()
{
}

void BushingElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    // NEEDED TO DEFINE THE DOFS OF THE ELEMENT

    auto& r_geometry = this->GetGeometry();
    const auto disp_x_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const auto rot_x_pos = r_geometry[0].GetDofPosition(ROTATION_X);

    // Resizing as needed
    if (rElementalDofList.size() != msElementSize) {
        rElementalDofList.resize(msElementSize);
    }

    for (IndexType i = 0; i < msNumNodes; ++i) {
        const IndexType index = i * msLocalSize;
        auto& r_node = r_geometry[i];
        rElementalDofList[index    ] = r_node.pGetDof(DISPLACEMENT_X, disp_x_pos);
        rElementalDofList[index + 1] = r_node.pGetDof(DISPLACEMENT_Y, disp_x_pos + 1);
        rElementalDofList[index + 2] = r_node.pGetDof(DISPLACEMENT_Z, disp_x_pos + 2);
        rElementalDofList[index + 3] = r_node.pGetDof(ROTATION_X, rot_x_pos);
        rElementalDofList[index + 4] = r_node.pGetDof(ROTATION_Y, rot_x_pos + 1);
        rElementalDofList[index + 5] = r_node.pGetDof(ROTATION_Z, rot_x_pos + 2);
    }
}

void BushingElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    // NEEDED TO DEFINE GLOBAL IDS FOR THE CORRECT ASSEMBLY
    if (rResult.size() != msElementSize) {
        rResult.resize(msElementSize, false);
    }

    auto& r_geometry = this->GetGeometry();
    const auto disp_x_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const auto rot_x_pos = r_geometry[0].GetDofPosition(ROTATION_X);

    for (IndexType i = 0; i < msNumNodes; ++i) {
        const IndexType index = i * msLocalSize;
        auto& r_node = r_geometry[i];
        rResult[index] = r_node.GetDof(DISPLACEMENT_X, disp_x_pos).EquationId();
        rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y, disp_x_pos + 1).EquationId();
        rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z, disp_x_pos + 2).EquationId();
        rResult[index + 3] = r_node.GetDof(ROTATION_X, rot_x_pos).EquationId();
        rResult[index + 4] = r_node.GetDof(ROTATION_Y, rot_x_pos + 1).EquationId();
        rResult[index + 5] = r_node.GetDof(ROTATION_Z, rot_x_pos + 2).EquationId();
    }
}

void BushingElement::GetValuesVector(
    Vector& rValues,
    int Step) const
{
    // GIVES THE VECTOR WITH THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT DISPLACEMENTS)
    if (rValues.size() != msElementSize) {
        rValues.resize(msElementSize, false);
    }

    BoundedVector<double, 12> u;
    this->GetValuesVector(u, Step);
    noalias(rValues) = u;
}

void BushingElement::GetFirstDerivativesVector(
    Vector& rValues,
    int Step) const
{
    // GIVES THE VECTOR WITH THE TIME DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT VELOCITIES)
    if (rValues.size() != msElementSize) {
        rValues.resize(msElementSize, false);
    }

    for (IndexType i = 0; i < msNumNodes; ++i) {
        const auto& r_node = GetGeometry()[i];
        const auto& vel = r_node.FastGetSolutionStepValue(VELOCITY, Step);
        const auto& avel = r_node.FastGetSolutionStepValue(ANGULAR_VELOCITY, Step);

        IndexType index = i * msLocalSize;
        rValues[index] = vel[0];
        rValues[index + 1] = vel[1];
        rValues[index + 2] = vel[2];
        rValues[index + 3] = avel[0];
        rValues[index + 4] = avel[1];
        rValues[index + 5] = avel[2];
    }
}

void BushingElement::GetSecondDerivativesVector(
    Vector& rValues,
    int Step) const
{
    // GIVES THE VECTOR WITH THE TIME SECOND DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT ACCELERATIONS)
    if (rValues.size() != msElementSize) {
        rValues.resize(msElementSize, false);
    }

    for (IndexType i = 0; i < GetGeometry().size(); ++i) {
        const auto& r_node = GetGeometry()[i];
        const auto& acc = r_node.FastGetSolutionStepValue(ACCELERATION, Step);
        const auto& aacc = r_node.FastGetSolutionStepValue(ANGULAR_ACCELERATION, Step);
        IndexType index = i * msLocalSize;

        rValues[index] = acc[0];
        rValues[index + 1] = acc[1];
        rValues[index + 2] = acc[2];
        rValues[index + 3] = aacc[0];
        rValues[index + 4] = aacc[1];
        rValues[index + 5] = aacc[2];
    }
}

void BushingElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rLeftHandSideMatrix.size1() != msElementSize || rLeftHandSideMatrix.size2() != msElementSize) {
        rLeftHandSideMatrix.resize(msElementSize, msElementSize, false);
        rRightHandSideVector.resize(msElementSize, false);
    }

    BoundedVector<double, 12> u;
    this->GetValuesVector(u, 0);

    BoundedMatrix<double, 12, 12> lhs;
    this->CalculateStiffnessMatrix(lhs, rCurrentProcessInfo);

    noalias(rLeftHandSideMatrix) = lhs;
    noalias(rRightHandSideVector) = -prod(lhs, u);

    KRATOS_CATCH("");
}

void BushingElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != msElementSize) {
        rRightHandSideVector.resize(msElementSize, false);
    }

    BoundedVector<double, 12> u;
    this->GetValuesVector(u, 0);

    BoundedMatrix<double, 12, 12> lhs;
    this->CalculateStiffnessMatrix(lhs, rCurrentProcessInfo);

    noalias(rRightHandSideVector) = -prod(lhs, u);

    KRATOS_CATCH("");
}

void BushingElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rLeftHandSideMatrix.size1() != msElementSize || rLeftHandSideMatrix.size2() != msElementSize) {
        rLeftHandSideMatrix.resize(msElementSize, msElementSize, false);
    }

    BoundedMatrix<double, 12, 12> lhs;
    this->CalculateStiffnessMatrix(lhs, rCurrentProcessInfo);

    noalias(rLeftHandSideMatrix) = lhs;

    KRATOS_CATCH("");
}

void BushingElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rMassMatrix.size1() != msElementSize || rMassMatrix.size2() != msElementSize) {
        rMassMatrix.resize(msElementSize, msElementSize, false);
    }

    rMassMatrix = ZeroMatrix(msElementSize, msElementSize);

    KRATOS_CATCH("");
}

void BushingElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rDampingMatrix.size1() != msElementSize || rDampingMatrix.size2() != msElementSize) {
        rDampingMatrix.resize(msElementSize, msElementSize, false);
    }

    rDampingMatrix = ZeroMatrix(msElementSize, msElementSize);

    KRATOS_CATCH("");
}

int BushingElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Verify that the dofs exist
    for (const auto& r_node : this->GetGeometry()) {
        // The displacement terms
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node)

        // The rotational terms
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION, r_node)
        KRATOS_CHECK_DOF_IN_NODE(ROTATION_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(ROTATION_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(ROTATION_Z, r_node)
    }

    return 0;

    KRATOS_CATCH("Problem in the Check in the BushingElement")
}

void BushingElement::GetValuesVector(
    BoundedVector<double, 12>& rValues,
    const int Step) const
{
    for (IndexType i = 0; i < msNumNodes; ++i) {
        const auto& r_node = GetGeometry()[i];
        const auto& disp = r_node.FastGetSolutionStepValue(DISPLACEMENT, Step);
        const auto& rot = r_node.FastGetSolutionStepValue(ROTATION, Step);

        const IndexType index = i * msLocalSize;
        rValues[index] = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];
        rValues[index + 3] = rot[0];
        rValues[index + 4] = rot[1];
        rValues[index + 5] = rot[2];
    }
}

void BushingElement::CalculateStiffnessMatrix(
    BoundedMatrix<double, 12, 12>& rStiffnessMatrix,
    const ProcessInfo& rProcessInfo) const
{
    KRATOS_TRY

    noalias(rStiffnessMatrix) = ZeroMatrix(msElementSize, msElementSize); //resetting LHS

    const Matrix& r_shape_functions = this->GetGeometry().ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);

    const array_1d<double, 6> stiffness {
        this->GetProperties().GetValue(NODAL_DISPLACEMENT_STIFFNESS_X,
                                       this->GetGeometry(),
                                       row(r_shape_functions, 0),
                                       rProcessInfo),
        this->GetProperties().GetValue(NODAL_DISPLACEMENT_STIFFNESS_Y,
                                       this->GetGeometry(),
                                       row(r_shape_functions, 0),
                                       rProcessInfo),
        this->GetProperties().GetValue(NODAL_DISPLACEMENT_STIFFNESS_Z,
                                       this->GetGeometry(),
                                       row(r_shape_functions, 0),
                                       rProcessInfo),
        this->GetProperties().GetValue(NODAL_ROTATIONAL_STIFFNESS_X,
                                       this->GetGeometry(),
                                       row(r_shape_functions, 0),
                                       rProcessInfo),
        this->GetProperties().GetValue(NODAL_ROTATIONAL_STIFFNESS_Y,
                                       this->GetGeometry(),
                                       row(r_shape_functions, 0),
                                       rProcessInfo),
        this->GetProperties().GetValue(NODAL_ROTATIONAL_STIFFNESS_Z,
                                       this->GetGeometry(),
                                       row(r_shape_functions, 0),
                                       rProcessInfo)
    };

    const double l = GetGeometry().Length();

    rStiffnessMatrix( 0,  0) = stiffness[0];
    rStiffnessMatrix( 0,  6) = -stiffness[0];
    rStiffnessMatrix( 1,  1) = stiffness[1];
    rStiffnessMatrix( 1,  5) = 0.5 * stiffness[1] * l;
    rStiffnessMatrix( 1,  7) = -stiffness[1];
    rStiffnessMatrix( 1, 11) = 0.5 * stiffness[1] * l;
    rStiffnessMatrix( 2,  2) = stiffness[2];
    rStiffnessMatrix( 2,  4) = -0.5 * stiffness[2] * l;
    rStiffnessMatrix( 2,  8) = -stiffness[2];
    rStiffnessMatrix( 2, 10) = -0.5 * stiffness[2] * l;
    rStiffnessMatrix( 3,  3) = stiffness[3];
    rStiffnessMatrix( 3,  9) = -stiffness[3];
    rStiffnessMatrix( 4,  2) = -0.5 * stiffness[2] * l;
    rStiffnessMatrix( 4,  4) = stiffness[4] + 0.25 * stiffness[2] * l * l;
    rStiffnessMatrix( 4,  8) = 0.5 * stiffness[2] * l;
    rStiffnessMatrix( 4, 10) = -stiffness[4] + 0.25 * stiffness[2] * l * l;
    rStiffnessMatrix( 5,  1) = 0.5 * stiffness[1] * l;
    rStiffnessMatrix( 5,  5) = stiffness[5] + 0.25 * stiffness[1] * l * l;
    rStiffnessMatrix( 5,  7) = -0.5 * stiffness[1] * l;
    rStiffnessMatrix( 5, 11) = -stiffness[5] + 0.25 * stiffness[1] * l * l;
    rStiffnessMatrix( 6,  0) = -stiffness[0];
    rStiffnessMatrix( 6,  6) = stiffness[0];
    rStiffnessMatrix( 7,  1) = -stiffness[1];
    rStiffnessMatrix( 7,  5) = -0.5 * stiffness[1] * l;
    rStiffnessMatrix( 7,  7) = stiffness[1];
    rStiffnessMatrix( 7, 11) = -0.5 * stiffness[1] * l;
    rStiffnessMatrix( 8,  2) = -stiffness[2];
    rStiffnessMatrix( 8,  4) = 0.5 * stiffness[2] * l;
    rStiffnessMatrix( 8,  8) = stiffness[2];
    rStiffnessMatrix( 8, 10) = 0.5 * stiffness[2] * l;
    rStiffnessMatrix( 9,  3) = -stiffness[3];
    rStiffnessMatrix( 9,  9) = stiffness[3];
    rStiffnessMatrix(10,  2) = -0.5 * stiffness[2] * l;
    rStiffnessMatrix(10,  4) = -stiffness[4] + 0.25 * stiffness[2] * l * l;
    rStiffnessMatrix(10,  8) = 0.5 * stiffness[2] * l;
    rStiffnessMatrix(10, 10) = stiffness[4] + 0.25 * stiffness[2] * l * l;
    rStiffnessMatrix(11,  1) = 0.5 * stiffness[1] * l;
    rStiffnessMatrix(11,  5) = -stiffness[5] + 0.25 * stiffness[1] * l * l;
    rStiffnessMatrix(11,  7) = -0.5 * stiffness[1] * l;
    rStiffnessMatrix(11, 11) = stiffness[5] + 0.25 * stiffness[1] * l * l;

    KRATOS_CATCH("");
}

void BushingElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
}
void BushingElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
}

} // Namespace Kratos
