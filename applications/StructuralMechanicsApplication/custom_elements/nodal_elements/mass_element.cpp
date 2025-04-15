// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes


// External includes


// Project includes
#include "includes/variables.h"
#include "includes/checks.h"
#include "mass_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{

///@name Operations
///@{

Element::Pointer MassElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<MassElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

Element::Pointer MassElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<MassElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

Element::Pointer MassElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<MassElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

void MassElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
    KRATOS_TRY;

    const SizeType system_size = GetGeometry().PointsNumber() * 3;

    if (rResult.size() != system_size) {
        rResult.resize(system_size, false);
    }

    SizeType local_index = 0;
    const SizeType d_pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

    for (const auto& r_node : GetGeometry()) {
        rResult[local_index++] = r_node.GetDof(DISPLACEMENT_X, d_pos).EquationId();
        rResult[local_index++] = r_node.GetDof(DISPLACEMENT_Y, d_pos + 1).EquationId();
        rResult[local_index++] = r_node.GetDof(DISPLACEMENT_Z, d_pos + 2).EquationId();
    }

    KRATOS_CATCH("")
}

void MassElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const
{
    KRATOS_TRY;

    const SizeType system_size = GetGeometry().PointsNumber() * 3;

    if (rElementalDofList.size() != system_size) {
        rElementalDofList.resize(system_size);
    }

    SizeType local_index = 0;

    for (const auto& r_node : GetGeometry()){
        rElementalDofList[local_index++] = r_node.pGetDof(DISPLACEMENT_X);
        rElementalDofList[local_index++] = r_node.pGetDof(DISPLACEMENT_Y);
        rElementalDofList[local_index++] = r_node.pGetDof(DISPLACEMENT_Z);
    }

    KRATOS_CATCH("")
}

void MassElement::GenericGetValuesVector(Vector& rValues, int Step, const ArrayVariableType& rVariable) const
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType system_size = number_of_nodes * 3;

    if (rValues.size() != system_size) {
        rValues.resize(system_size, false);
    }

    for (IndexType i=0; i<number_of_nodes; ++i) {
        int index = i * 3;
        const auto& r_vals = GetGeometry()[i].FastGetSolutionStepValue(rVariable, Step);

        rValues[index]     = r_vals[0];
        rValues[index + 1] = r_vals[1];
        rValues[index + 2] = r_vals[2];
    }

    KRATOS_CATCH("")
}

void MassElement::GetValuesVector(Vector& rValues, int Step) const
{
    GenericGetValuesVector(rValues, Step, DISPLACEMENT);
}

void MassElement::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    GenericGetValuesVector(rValues, Step, VELOCITY);
}

void MassElement::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    GenericGetValuesVector(rValues, Step, ACCELERATION);
}

void MassElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        // currently the mass is calculated using the current configuration
        // hence the mass is saved as a member
        mMass = GetElementMass();
    }
}

void MassElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

void MassElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType system_size = GetGeometry().PointsNumber() * 3;
    if (rLeftHandSideMatrix.size1() != system_size) {
        rLeftHandSideMatrix.resize(system_size, system_size, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(system_size, system_size);
}

void MassElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const SizeType system_size = GetGeometry().PointsNumber() * 3;
    if (rRightHandSideVector.size() != system_size) {
        rRightHandSideVector.resize(system_size, false);
    }

    rRightHandSideVector = ZeroVector(system_size);

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.PointsNumber();

    Vector lump_fact = ZeroVector(number_of_nodes);
    r_geom.LumpingFactors(lump_fact);

    for (SizeType i=0; i<number_of_nodes; ++i) {
        const double temp = lump_fact[i] * mMass;

        for (SizeType j=0; j<3; ++j) {
            const SizeType index = i * 3 + j;
            rRightHandSideVector[index] += temp * r_geom[i].FastGetSolutionStepValue(VOLUME_ACCELERATION)[j];
        }
    }
    KRATOS_CATCH("")
}

void MassElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const auto& r_geom = GetGeometry();

    // lumped mass matrix
    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType system_size = number_of_nodes * 3;

    if (rMassMatrix.size1() != system_size) {
        rMassMatrix.resize(system_size, system_size, false);
    }

    noalias(rMassMatrix) = ZeroMatrix(system_size, system_size);

    Vector lump_fact = ZeroVector(number_of_nodes);
    r_geom.LumpingFactors(lump_fact);

    for (SizeType i=0; i<number_of_nodes; ++i) {
        const double temp = lump_fact[i] * mMass;

        for (SizeType j=0; j<3; ++j) {
            const SizeType index = i * 3 + j;
            rMassMatrix(index, index) = temp;
        }
    }

    KRATOS_CATCH("")
}

void MassElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        GetGeometry().PointsNumber() * 3);
}

double MassElement::GetElementMass() const
{
    // TODO this should be coming from total_mass_process (currently not working bcs not threadsafe!)
    // => we need a better solution for accessing things in the initial condiguration
    KRATOS_TRY

    const auto& r_geom = GetGeometry();
    const SizeType local_dim = GetGeometry().LocalSpaceDimension();

    double total_mass;

    if (local_dim == 1) { // line
        total_mass = GetProperties()[CROSS_AREA] * StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);

    } else if (local_dim == 2) { // tri / quad
        total_mass = GetProperties()[THICKNESS]* r_geom.Area();

    } else {
        KRATOS_ERROR << "Wrong local space dimension!" << std::endl;
    }

    total_mass *= StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);

    return total_mass;

    KRATOS_CATCH("")
}

int MassElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1) << "Element found with Id 0 or negative" << std::endl;

    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType local_dim = GetGeometry().LocalSpaceDimension();

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : GetGeometry()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION, r_node);

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node);
    }

    KRATOS_ERROR_IF_NOT(GetProperties().Has(DENSITY)) << "DENSITY not provided for this element #" << Id() << std::endl;

    // dimension specific checks
    if (local_dim == 1) { // line
        KRATOS_ERROR_IF(number_of_nodes != 2) << "Wrong number of nodes!" << std::endl;

        KRATOS_ERROR_IF(StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this) <= 0) << "On element #" << Id() << "; Length cannot be less than or equal to 0" << std::endl;

        KRATOS_ERROR_IF(GetProperties().Has(CROSS_AREA) == false || GetProperties()[CROSS_AREA] <= numerical_limit) << "CROSS_AREA not provided for this element #" << Id() << std::endl;
    } else if (local_dim == 2) { // tri / quad
        KRATOS_ERROR_IF(number_of_nodes != 3 && number_of_nodes != 4) << "Wrong number of nodes!" << std::endl;

        KRATOS_ERROR_IF(GetGeometry().Area() <= 0) << "On element #" << Id() << "; Area cannot be less than or equal to 0" << std::endl;

        KRATOS_ERROR_IF(GetProperties().Has(THICKNESS) == false || GetProperties()[THICKNESS] <= numerical_limit) << "THICKNESS not provided for this element #" << Id() << std::endl;
    } else {
        KRATOS_ERROR << "Wrong local space dimension, can only be 2 (line) or 3 (surface), got: " << local_dim << " !" << std::endl;
    }

    return 0;

    KRATOS_CATCH("");
}

///@}
///@name Input and output
///@{

/// Turn back information as a string.

std::string MassElement::Info() const
{
    std::stringstream buffer;
    buffer << "MassElement #" << Id();
    return buffer.str();
}

/// Print information about this object.

void MassElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "MassElement #" << Id();
}

/// Print object's data.

void MassElement::PrintData(std::ostream& rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

///@}
///@name Serialization
///@{

void MassElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    rSerializer.save("mass", mMass);
}

void MassElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
    rSerializer.load("mass", mMass);
}

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >> (std::istream& rIStream, MassElement& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const MassElement& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.
