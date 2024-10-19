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
#include <cmath>
#include <limits>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "structural_mechanics_application_variables.h"

// Include base h
#include "custom_elements/bushing_element.h"

namespace Kratos {

double BushingElement::ConstantStiffness::GetValue(
    const double Value,
    const Properties& rProperties) const
{
    KRATOS_TRY

    return rProperties.GetValue(*(this->mpVariable));

    KRATOS_CATCH("");
}

double BushingElement::NonLinearStiffness::GetValue(
    const double Value,
    const Properties& rProperties) const
{
    KRATOS_TRY

    const auto& r_table = rProperties.GetTable(*mpVariableX, *mpVariableY);
    if (std::abs(Value) > std::numeric_limits<double>::epsilon()) {
        return r_table.GetValue(Value) / Value;
    } else {
        return r_table.GetDerivative(Value);
    }

    KRATOS_CATCH("");
}

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

    const double l = this->GetGeometry().Length();

    BoundedVector<double, 12> u;
    this->GetValuesVector(u, 0);

    array_1d<double, 6> stiffness_values;
    this->CalculateStiffnessValues(stiffness_values, u, l);

    BoundedMatrix<double, 12, 12> lhs;
    this->CalculateStiffnessMatrix(lhs, stiffness_values, l);

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

    const double l = this->GetGeometry().Length();

    BoundedVector<double, 12> u;
    this->GetValuesVector(u, 0);

    array_1d<double, 6> stiffness_values;
    this->CalculateStiffnessValues(stiffness_values, u, l);

    BoundedMatrix<double, 12, 12> lhs;
    this->CalculateStiffnessMatrix(lhs, stiffness_values, l);

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

    const double l = this->GetGeometry().Length();

    BoundedVector<double, 12> u;
    this->GetValuesVector(u, 0);

    array_1d<double, 6> stiffness_values;
    this->CalculateStiffnessValues(stiffness_values, u, l);

    BoundedMatrix<double, 12, 12> lhs;
    this->CalculateStiffnessMatrix(lhs, stiffness_values, l);

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

    const auto& r_properties = this->GetProperties();

    const array_1d<bool, 3> is_disp_stiffness_defined{
        r_properties.Has(NODAL_DISPLACEMENT_STIFFNESS_X),
        r_properties.Has(NODAL_DISPLACEMENT_STIFFNESS_Y),
        r_properties.Has(NODAL_DISPLACEMENT_STIFFNESS_Z)
    };

    const array_1d<bool, 3> is_disp_forces_defined{
        r_properties.HasTable(DISPLACEMENT_X, FORCE_X),
        r_properties.HasTable(DISPLACEMENT_Y, FORCE_Y),
        r_properties.HasTable(DISPLACEMENT_Z, FORCE_Z)
    };

    const array_1d<bool, 3> is_rot_stiffness_defined{
        r_properties.Has(NODAL_ROTATIONAL_STIFFNESS_X),
        r_properties.Has(NODAL_ROTATIONAL_STIFFNESS_Y),
        r_properties.Has(NODAL_ROTATIONAL_STIFFNESS_Z)
    };

    const array_1d<bool, 3> is_rot_forces_defined{
        r_properties.HasTable(ROTATION_X, MOMENT_X),
        r_properties.HasTable(ROTATION_Y, MOMENT_Y),
        r_properties.HasTable(ROTATION_Z, MOMENT_Z)
    };

    for (IndexType i = 0; i < 3; ++i) {
        KRATOS_ERROR_IF(is_disp_stiffness_defined[i] && is_disp_forces_defined[i])
            << "Element with id " << this->Id() << " has both NODAL_DISPLACEMENT_STIFFNESS "
            << "and (DISPLACEMENT, FORCE) table defined in direction " << i << " which is not allowed.";

        KRATOS_ERROR_IF_NOT(is_disp_stiffness_defined[i] || is_disp_forces_defined[i])
            << "Element with id " << this->Id() << " has not defined either NODAL_DISPLACEMENT_STIFFNESS "
            << "or (DISPLACEMENT, FORCE) table in direction " << i << ".";

        KRATOS_ERROR_IF(is_rot_stiffness_defined[i] && is_rot_forces_defined[i])
            << "Element with id " << this->Id() << " has both NODAL_ROTATIONAL_STIFFNESS "
            << "and (ROTATION, MOMENT) table defined in direction " << i << " which is not allowed.";

        KRATOS_ERROR_IF_NOT(is_rot_stiffness_defined[i] || is_rot_forces_defined[i])
            << "Element with id " << this->Id() << " has not defined either NODAL_ROTATIONAL_STIFFNESS "
            << "or (ROTATION, MOMENT) table in direction " << i << ".";
    }

    return 0;

    KRATOS_CATCH("Problem in the Check in the BushingElement")
}

void BushingElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_properties = this->GetProperties();

    const array_1d<bool, 3> is_disp_stiffness_defined{
        r_properties.Has(NODAL_DISPLACEMENT_STIFFNESS_X),
        r_properties.Has(NODAL_DISPLACEMENT_STIFFNESS_Y),
        r_properties.Has(NODAL_DISPLACEMENT_STIFFNESS_Z)
    };

    const array_1d<bool, 3> is_rot_stiffness_defined{
        r_properties.Has(NODAL_ROTATIONAL_STIFFNESS_X),
        r_properties.Has(NODAL_ROTATIONAL_STIFFNESS_Y),
        r_properties.Has(NODAL_ROTATIONAL_STIFFNESS_Z)
    };

    if (is_disp_stiffness_defined[0]) {
        mStiffnessGetters[0] = Kratos::make_unique<ConstantStiffness>(NODAL_DISPLACEMENT_STIFFNESS_X);
    } else {
        mStiffnessGetters[0] = Kratos::make_unique<NonLinearStiffness>(DISPLACEMENT_X, FORCE_X);
    }

    if (is_disp_stiffness_defined[1]) {
        mStiffnessGetters[1] = Kratos::make_unique<ConstantStiffness>(NODAL_DISPLACEMENT_STIFFNESS_Y);
    } else {
        mStiffnessGetters[1] = Kratos::make_unique<NonLinearStiffness>(DISPLACEMENT_Y, FORCE_Y);
    }

    if (is_disp_stiffness_defined[2]) {
        mStiffnessGetters[2] = Kratos::make_unique<ConstantStiffness>(NODAL_DISPLACEMENT_STIFFNESS_Z);
    } else {
        mStiffnessGetters[2] = Kratos::make_unique<NonLinearStiffness>(DISPLACEMENT_Z, FORCE_Z);
    }

    if (is_rot_stiffness_defined[0]) {
        mStiffnessGetters[3] = Kratos::make_unique<ConstantStiffness>(NODAL_ROTATIONAL_STIFFNESS_X);
    } else {
        mStiffnessGetters[3] = Kratos::make_unique<NonLinearStiffness>(ROTATION_X, MOMENT_X);
    }

    if (is_rot_stiffness_defined[1]) {
        mStiffnessGetters[4] = Kratos::make_unique<ConstantStiffness>(NODAL_ROTATIONAL_STIFFNESS_Y);
    } else {
        mStiffnessGetters[4] = Kratos::make_unique<NonLinearStiffness>(ROTATION_Y, MOMENT_Y);
    }

    if (is_rot_stiffness_defined[2]) {
        mStiffnessGetters[5] = Kratos::make_unique<ConstantStiffness>(NODAL_ROTATIONAL_STIFFNESS_Z);
    } else {
        mStiffnessGetters[5] = Kratos::make_unique<NonLinearStiffness>(ROTATION_Z, MOMENT_Z);
    }

    KRATOS_CATCH("");
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

void BushingElement::CalculateStiffnessValues(
    array_1d<double, 6>& rStiffnessValues,
    const BoundedVector<double, 12>& rValues,
    const double Length) const
{
    KRATOS_TRY

    const auto& r_properties = this->GetProperties();

    array_1d<double, 6> deflections {
        rValues[6] - rValues[0],
        rValues[7] - rValues[1] - 0.5 * Length * (rValues[5] + rValues[11]),
        rValues[8] - rValues[2] + 0.5 * Length * (rValues[4] + rValues[10]),
        rValues[9] - rValues[3],
        rValues[10] - rValues[4],
        rValues[11] - rValues[5]
    };


    for (IndexType i = 0; i < 6; ++i) {
        rStiffnessValues[i] = mStiffnessGetters[i]->GetValue(deflections[i], r_properties);
    }

    KRATOS_CATCH("");
}

void BushingElement::CalculateStiffnessMatrix(
    BoundedMatrix<double, 12, 12>& rStiffnessMatrix,
    const array_1d<double, 6>& rStiffnessValues,
    const double Length) const
{
    KRATOS_TRY

    noalias(rStiffnessMatrix) = ZeroMatrix(msElementSize, msElementSize); //resetting LHS

    rStiffnessMatrix( 0,  0) = rStiffnessValues[0];
    rStiffnessMatrix( 0,  6) = -rStiffnessValues[0];
    rStiffnessMatrix( 1,  1) = rStiffnessValues[1];
    rStiffnessMatrix( 1,  5) = 0.5 * rStiffnessValues[1] * Length;
    rStiffnessMatrix( 1,  7) = -rStiffnessValues[1];
    rStiffnessMatrix( 1, 11) = 0.5 * rStiffnessValues[1] * Length;
    rStiffnessMatrix( 2,  2) = rStiffnessValues[2];
    rStiffnessMatrix( 2,  4) = -0.5 * rStiffnessValues[2] * Length;
    rStiffnessMatrix( 2,  8) = -rStiffnessValues[2];
    rStiffnessMatrix( 2, 10) = -0.5 * rStiffnessValues[2] * Length;
    rStiffnessMatrix( 3,  3) = rStiffnessValues[3];
    rStiffnessMatrix( 3,  9) = -rStiffnessValues[3];
    rStiffnessMatrix( 4,  2) = -0.5 * rStiffnessValues[2] * Length;
    rStiffnessMatrix( 4,  4) = rStiffnessValues[4] + 0.25 * rStiffnessValues[2] * Length * Length;
    rStiffnessMatrix( 4,  8) = 0.5 * rStiffnessValues[2] * Length;
    rStiffnessMatrix( 4, 10) = -rStiffnessValues[4] + 0.25 * rStiffnessValues[2] * Length * Length;
    rStiffnessMatrix( 5,  1) = 0.5 * rStiffnessValues[1] * Length;
    rStiffnessMatrix( 5,  5) = rStiffnessValues[5] + 0.25 * rStiffnessValues[1] * Length * Length;
    rStiffnessMatrix( 5,  7) = -0.5 * rStiffnessValues[1] * Length;
    rStiffnessMatrix( 5, 11) = -rStiffnessValues[5] + 0.25 * rStiffnessValues[1] * Length * Length;
    rStiffnessMatrix( 6,  0) = -rStiffnessValues[0];
    rStiffnessMatrix( 6,  6) = rStiffnessValues[0];
    rStiffnessMatrix( 7,  1) = -rStiffnessValues[1];
    rStiffnessMatrix( 7,  5) = -0.5 * rStiffnessValues[1] * Length;
    rStiffnessMatrix( 7,  7) = rStiffnessValues[1];
    rStiffnessMatrix( 7, 11) = -0.5 * rStiffnessValues[1] * Length;
    rStiffnessMatrix( 8,  2) = -rStiffnessValues[2];
    rStiffnessMatrix( 8,  4) = 0.5 * rStiffnessValues[2] * Length;
    rStiffnessMatrix( 8,  8) = rStiffnessValues[2];
    rStiffnessMatrix( 8, 10) = 0.5 * rStiffnessValues[2] * Length;
    rStiffnessMatrix( 9,  3) = -rStiffnessValues[3];
    rStiffnessMatrix( 9,  9) = rStiffnessValues[3];
    rStiffnessMatrix(10,  2) = -0.5 * rStiffnessValues[2] * Length;
    rStiffnessMatrix(10,  4) = -rStiffnessValues[4] + 0.25 * rStiffnessValues[2] * Length * Length;
    rStiffnessMatrix(10,  8) = 0.5 * rStiffnessValues[2] * Length;
    rStiffnessMatrix(10, 10) = rStiffnessValues[4] + 0.25 * rStiffnessValues[2] * Length * Length;
    rStiffnessMatrix(11,  1) = 0.5 * rStiffnessValues[1] * Length;
    rStiffnessMatrix(11,  5) = -rStiffnessValues[5] + 0.25 * rStiffnessValues[1] * Length * Length;
    rStiffnessMatrix(11,  7) = -0.5 * rStiffnessValues[1] * Length;
    rStiffnessMatrix(11, 11) = rStiffnessValues[5] + 0.25 * rStiffnessValues[1] * Length * Length;

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
