// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Mahmoud Zidan
//

// System includes
#include<iostream>

// External includes

// Project includes
#include "custom_elements/fiber_beam_column_element_3D2N.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "includes/checks.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Constructor
FiberBeamColumnElement3D2N::FiberBeamColumnElement3D2N(IndexType NewId)
    : Element(NewId) {}

/// Constructor using an array of nodes
FiberBeamColumnElement3D2N::FiberBeamColumnElement3D2N(IndexType NewId, const NodesArrayType& rThisNodes)
    : Element(NewId, rThisNodes) {}

/// Constructor using Geometry
FiberBeamColumnElement3D2N::FiberBeamColumnElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

/// Constructor using Properties
FiberBeamColumnElement3D2N::FiberBeamColumnElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

/// Copy Constructor
FiberBeamColumnElement3D2N::FiberBeamColumnElement3D2N(FiberBeamColumnElement3D2N const& rOther)
    : Element(rOther) {}

/// Destructor
FiberBeamColumnElement3D2N::~FiberBeamColumnElement3D2N() {}

///@}
///@name Operators
///@{

/// Assigment Operator
FiberBeamColumnElement3D2N & FiberBeamColumnElement3D2N::operator=(FiberBeamColumnElement3D2N const& rOther)
{
    BaseType::operator=(rOther);
    Flags::operator=(rOther);
    return *this;
}

///@}
///@name Operations
///@{

Element::Pointer FiberBeamColumnElement3D2N::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<FiberBeamColumnElement3D2N>(NewId, GetGeometry().Create(rThisNodes), pProperties);
    KRATOS_CATCH("")
}

Element::Pointer FiberBeamColumnElement3D2N::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<FiberBeamColumnElement3D2N>(NewId, pGeom, pProperties);
    KRATOS_CATCH("")
}

Element::Pointer FiberBeamColumnElement3D2N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<FiberBeamColumnElement3D2N>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    KRATOS_CATCH("")
}

void FiberBeamColumnElement3D2N::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{

    if (rResult.size() != msGlobalSize) {
        rResult.resize(msGlobalSize, false);
    }

    const unsigned int global_size_per_node = 6;
    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * global_size_per_node;
        rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
        rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
        rResult[index + 5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
    }
}

void FiberBeamColumnElement3D2N::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{

    if (rElementalDofList.size() != msGlobalSize) {
        rElementalDofList.resize(msGlobalSize);
    }

    const unsigned int global_size_per_node = 6;
    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * global_size_per_node;
        rElementalDofList[index] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
        rElementalDofList[index + 3] = GetGeometry()[i].pGetDof(ROTATION_X);
        rElementalDofList[index + 4] = GetGeometry()[i].pGetDof(ROTATION_Y);
        rElementalDofList[index + 5] = GetGeometry()[i].pGetDof(ROTATION_Z);
    }
}

void FiberBeamColumnElement3D2N::Initialize()
{
    KRATOS_TRY

    auto prop = GetProperties();
    const IndexType number_of_sections = prop[NUMBER_OF_SECTIONS];
    mSections.resize(number_of_sections);

    const GeometryData::IntegrationPointsArrayType integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());
    for (SizeType i = 0; i < number_of_sections; ++i) {
        // mSections[i] = FiberBeamColumnSection( i+1, integration_points[i](0), integration_points[i].Weight(), pGetProperties());
        mSections[i] = FiberBeamColumnSection(
            i+1,
            integration_points[i](0),
            integration_points[i].Weight(),
            prop[NUMBER_FIBERS_Y],
            prop[NUMBER_FIBERS_Z],
            prop[BEAM_WIDTH],
            prop[BEAM_HEIGHT]);
    }
    for (SizeType i = 0; i < number_of_sections; ++i) {
        mSections[i].Initialize();
        KRATOS_WATCH(mSections[i])
    }

    KRATOS_CATCH("")
}

GeometryData::IntegrationMethod FiberBeamColumnElement3D2N::GetIntegrationMethod() const
{
    const IndexType number_of_sections = GetProperties()[NUMBER_OF_SECTIONS];
    switch (number_of_sections){
        case 2: return GeometryData::GI_GAUSS_2;
        case 3: return GeometryData::GI_GAUSS_3;
        case 4: return GeometryData::GI_GAUSS_4;
        case 5: return GeometryData::GI_GAUSS_5;
        default:
        {
            KRATOS_ERROR << "NUMBER_OF_SECTIONS can only be in range [2, 5]."
            << "The given NUMBER_OF_SECTIONS is: " << number_of_sections << std::endl;
            break;
        }
    }
}

void FiberBeamColumnElement3D2N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void FiberBeamColumnElement3D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // resizing the matrices + create memory for LHS
    rLeftHandSideMatrix = ZeroMatrix(msGlobalSize, msGlobalSize);
    BoundedMatrix<double,msGlobalSize,msLocalSize> rotation_matrix = CreateTransformationMatrix();
    // creating LHS
    BoundedMatrix<double,msLocalSize,msLocalSize> local_stiffness = CreateElementLocalStiffnessMatrix(rCurrentProcessInfo);
    BoundedMatrix<double,msGlobalSize,msLocalSize> aux_matrix = prod(rotation_matrix, local_stiffness);
    noalias(rLeftHandSideMatrix) += prod(aux_matrix, Matrix(trans(rotation_matrix)));

    KRATOS_CATCH("")
}

BoundedMatrix<double, FiberBeamColumnElement3D2N::msLocalSize, FiberBeamColumnElement3D2N::msLocalSize>
FiberBeamColumnElement3D2N::CreateElementLocalStiffnessMatrix(ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    BoundedMatrix<double,msLocalSize,msLocalSize> local_stiffness_matrix = ZeroMatrix(msLocalSize, msLocalSize);
    // FIXME: calculate local stiffness here
    return local_stiffness_matrix;
    KRATOS_CATCH("")
}

void FiberBeamColumnElement3D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    rRightHandSideVector = ZeroVector(msLocalSize);
    BoundedMatrix<double,msGlobalSize,msLocalSize> rotation_matrix = CreateTransformationMatrix();
    BoundedVector<double,msLocalSize> resisting_forces = CreateElementLocalResistingForces(rCurrentProcessInfo);
    noalias(rRightHandSideVector) -= prod(rotation_matrix, resisting_forces);
    KRATOS_CATCH("")
}

BoundedVector<double, FiberBeamColumnElement3D2N::msLocalSize>
FiberBeamColumnElement3D2N::CreateElementLocalResistingForces(ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    BoundedVector<double,msLocalSize> resisting_forces = ZeroVector(msLocalSize);
    // FIXME: calculate resisting forces here
    return resisting_forces;
    KRATOS_CATCH("")
}

BoundedMatrix<double,FiberBeamColumnElement3D2N::msGlobalSize,FiberBeamColumnElement3D2N::msLocalSize>
FiberBeamColumnElement3D2N::CreateTransformationMatrix() const
{
    KRATOS_TRY
    BoundedMatrix<double,msGlobalSize,msLocalSize> rotation_matrix = ZeroMatrix(msGlobalSize,msLocalSize);
    // FIXME: calculate rotation matrix here
    return rotation_matrix;
    KRATOS_CATCH("")
}

// void FiberBeamColumnElement3D2N::GetValuesVector(
//     Vector& rValues,
//     int Step = 0)
// {}

// int FiberBeamColumnElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
// {
//     KRATOS_TRY

//     KRATOS_ERROR_IF(this->Id() < 1) <<"FiberBeamColumnElement3D2N found with Id 0 or negative" << std::endl;

//     KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0) << "On FiberBeamColumnElement3D2N -> "
//         << this->Id() <<  "; Area cannot be less than or equal to 0" << std::endl;

//     return 0;

//     KRATOS_CATCH("");
// }


std::string FiberBeamColumnElement3D2N::Info() const {
    std::stringstream buffer;
    buffer << "FiberBeamColumnElement3D2N #" << Id();
    return buffer.str();
}

void FiberBeamColumnElement3D2N::PrintInfo(std::ostream& rOStream) const {
    rOStream << "FiberBeamColumnElement3D2N #" << Id();
}

void FiberBeamColumnElement3D2N::PrintData(std::ostream& rOStream) const {
    pGetGeometry()->PrintData(rOStream);
}


void FiberBeamColumnElement3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    rSerializer.save("mSections", mSections);
}
void FiberBeamColumnElement3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    rSerializer.load("mSections", mSections);
}

} // namespace Kratos.
