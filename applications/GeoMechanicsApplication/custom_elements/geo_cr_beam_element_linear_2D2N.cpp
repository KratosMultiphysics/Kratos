// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//  license: 	 structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//                Vahid Galavi
//
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/geo_cr_beam_element_linear_2D2N.hpp"
#include "custom_utilities/static_condensation_utility.h"
#include "includes/define.h"
#include "geo_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{
GeoCrBeamElementLinear2D2N::GeoCrBeamElementLinear2D2N(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : GeoCrBeamElement2D2N(NewId, pGeometry) {}

GeoCrBeamElementLinear2D2N::GeoCrBeamElementLinear2D2N(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : GeoCrBeamElement2D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
GeoCrBeamElementLinear2D2N::Create(IndexType NewId,
                                NodesArrayType const& rThisNodes,
                                PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<GeoCrBeamElementLinear2D2N>(
               NewId, rGeom.Create(rThisNodes), pProperties);
}

Element::Pointer
GeoCrBeamElementLinear2D2N::Create(IndexType NewId,
                                GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoCrBeamElementLinear2D2N>(
               NewId, pGeom, pProperties);
}

GeoCrBeamElementLinear2D2N::~GeoCrBeamElementLinear2D2N() {}

void GeoCrBeamElementLinear2D2N::
    CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    // Kt
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    Vector nodal_deformation = ZeroVector(msElementSize);
    GetValuesVector(nodal_deformation);

    rRightHandSideVector = ZeroVector(msElementSize);
    noalias(mInternalGlobalForces) = prod(rLeftHandSideMatrix, nodal_deformation);

    noalias(rRightHandSideVector) -= (mInternalGlobalForces + mInternalGlobalForcesFinalizedPrevious);

    noalias(rRightHandSideVector) += CalculateBodyForces();

    KRATOS_CATCH("")
}

void GeoCrBeamElementLinear2D2N::
    CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    Vector nodal_deformation = ZeroVector(msElementSize);
    GetValuesVector(nodal_deformation);
    rRightHandSideVector = ZeroVector(msElementSize);
    noalias(mInternalGlobalForces) = prod(mK_Master, nodal_deformation);
    noalias(rRightHandSideVector) -= (mInternalGlobalForces + mInternalGlobalForcesFinalizedPrevious);
    noalias(rRightHandSideVector) += CalculateBodyForces();
    KRATOS_CATCH("")
}

void GeoCrBeamElementLinear2D2N::
    CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    rLeftHandSideMatrix = CreateElementStiffnessMatrix_Total();

    //// start static condensation
    if (Has(CONDENSED_DOF_LIST)) {
        Vector dof_list_input = GetValue(CONDENSED_DOF_LIST);
        std::vector<int> dofList(dof_list_input.size());
        for (SizeType i = 0; i < dof_list_input.size(); ++i) {
            dofList[i] = dof_list_input[i];
        }
        GeoStaticCondensationUtility::CondenseLeftHandSide(*this, rLeftHandSideMatrix, dofList);
    }
    //// end static condensation

    GlobalizeMatrix(rLeftHandSideMatrix);
    mK_Master = rLeftHandSideMatrix;
    KRATOS_CATCH("")
}

/////////////////////////////////////////////////
///////////// CUSTOM FUNCTIONS --->>
/////////////////////////////////////////////////

BoundedMatrix<double, GeoCrBeamElement2D2N::msElementSize,
GeoCrBeamElement2D2N::msElementSize>
GeoCrBeamElementLinear2D2N::CreateRotationMatrix()
{
    KRATOS_TRY;
    const double current_element_angle = CalculateInitialElementAngle();
    const double c = std::cos(current_element_angle);
    const double s = std::sin(current_element_angle);

    BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix =
        ZeroMatrix(msElementSize, msElementSize);

    rotation_matrix(0, 0) = c;
    rotation_matrix(0, 1) = -s;
    rotation_matrix(1, 0) = s;
    rotation_matrix(1, 1) = c;
    rotation_matrix(2, 2) = 1.00;

    rotation_matrix(3, 3) = c;
    rotation_matrix(3, 4) = -s;
    rotation_matrix(4, 3) = s;
    rotation_matrix(4, 4) = c;
    rotation_matrix(5, 5) = 1.00;

    return rotation_matrix;
    KRATOS_CATCH("")
}

void GeoCrBeamElementLinear2D2N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    // element with two nodes can only represent results at one node
    const unsigned int& write_points_number =
        GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);
    if (rOutput.size() != write_points_number) {
        rOutput.resize(write_points_number);
    }

    Matrix left_hand_side_matrix = CreateElementStiffnessMatrix_Total();

    Vector nodal_deformation = ZeroVector(msElementSize);
    GetValuesVector(nodal_deformation);

    BoundedMatrix<double, msElementSize, msElementSize> transformation_matrix =
        CreateRotationMatrix();
    // calculate local displacements
    nodal_deformation = prod((trans(transformation_matrix)), nodal_deformation);

    //// start static back condensation
    if (Has(CONDENSED_DOF_LIST)) {
        Vector dof_list_input = GetValue(CONDENSED_DOF_LIST);
        std::vector<int> dofList(dof_list_input.size());
        for (SizeType i = 0; i < dof_list_input.size(); ++i) {
            dofList[i] = dof_list_input[i];
        }
        Vector nodal_deformation_temp = nodal_deformation;
        GeoStaticCondensationUtility::ConvertingCondensation(
            *this, nodal_deformation_temp, nodal_deformation, dofList,
            left_hand_side_matrix);
    }
    //// end static back condensation

    Vector stress = prod(left_hand_side_matrix, nodal_deformation);

    // rOutput[GP 1,2,3][x,y,z]

    if (rVariable == MOMENT) {
        rOutput[0][0] = 0.00;
        rOutput[1][0] = 0.00;
        rOutput[2][0] = 0.00;

        rOutput[0][1] = 0.00;
        rOutput[1][1] = 0.00;
        rOutput[2][1] = 0.00;

        rOutput[0][2] = 1.0 * stress[2] * 0.75 - stress[5] * 0.25;
        rOutput[1][2] = 1.0 * stress[2] * 0.50 - stress[5] * 0.50;
        rOutput[2][2] = 1.0 * stress[2] * 0.25 - stress[5] * 0.75;
    }
    if (rVariable == FORCE) {
        rOutput[0][0] = -1.0 * stress[0] * 0.75 + stress[3] * 0.25;
        rOutput[1][0] = -1.0 * stress[0] * 0.50 + stress[3] * 0.50;
        rOutput[2][0] = -1.0 * stress[0] * 0.25 + stress[3] * 0.75;

        rOutput[0][1] = -1.0 * stress[1] * 0.75 + stress[4] * 0.25;
        rOutput[1][1] = -1.0 * stress[1] * 0.50 + stress[4] * 0.50;
        rOutput[2][1] = -1.0 * stress[1] * 0.25 + stress[4] * 0.75;

        rOutput[0][2] = 0.00;
        rOutput[1][2] = 0.00;
        rOutput[2][2] = 0.00;
    }

    KRATOS_CATCH("")
}

double GeoCrBeamElementLinear2D2N::CalculateLength() const
{
    KRATOS_TRY;
    return GeoStructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
    KRATOS_CATCH("")
}

void GeoCrBeamElementLinear2D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GeoCrBeamElement2D2N);
    rSerializer.save("MasterStiffnessMatrix", mK_Master);
}

void GeoCrBeamElementLinear2D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GeoCrBeamElement2D2N);
    rSerializer.load("MasterStiffnessMatrix", mK_Master);
}

} // namespace Kratos.
