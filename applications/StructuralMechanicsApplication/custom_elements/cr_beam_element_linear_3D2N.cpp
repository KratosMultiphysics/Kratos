// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/cr_beam_element_linear_3D2N.hpp"
#include "custom_utilities/static_condensation_utility.h"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos {

CrBeamElementLinear3D2N::CrBeamElementLinear3D2N(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : CrBeamElement3D2N(NewId, pGeometry) {}

CrBeamElementLinear3D2N::CrBeamElementLinear3D2N(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : CrBeamElement3D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
CrBeamElementLinear3D2N::Create(IndexType NewId,
                                NodesArrayType const &rThisNodes,
                                PropertiesType::Pointer pProperties) const {
  const GeometryType &rGeom = this->GetGeometry();
  return Kratos::make_shared<CrBeamElementLinear3D2N>(
      NewId, rGeom.Create(rThisNodes), pProperties);
}

Element::Pointer
CrBeamElementLinear3D2N::Create(IndexType NewId,
                                 GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const {
  return Kratos::make_shared<CrBeamElementLinear3D2N>(
      NewId, pGeom, pProperties);
}

CrBeamElementLinear3D2N::~CrBeamElementLinear3D2N() {}

void CrBeamElementLinear3D2N::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

  Vector nodal_deformation = ZeroVector(msElementSize);
  this->GetValuesVector(nodal_deformation);
  rRightHandSideVector = ZeroVector(msElementSize);
  rRightHandSideVector -= prod(rLeftHandSideMatrix, nodal_deformation);

  // add bodyforces
  rRightHandSideVector += this->CalculateBodyForces();
  this->IncrementIterationCounter();
  KRATOS_CATCH("")
}

void CrBeamElementLinear3D2N::CalculateRightHandSide(
    VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;
  rRightHandSideVector = ZeroVector(msElementSize);

  Matrix left_hand_side_matrix = ZeroMatrix(msElementSize, msElementSize);
  this->CalculateLeftHandSide(left_hand_side_matrix, rCurrentProcessInfo);
  Vector nodal_deformation = ZeroVector(msElementSize);
  this->GetValuesVector(nodal_deformation);
  rRightHandSideVector = ZeroVector(msElementSize);
  noalias(rRightHandSideVector) -=
      prod(left_hand_side_matrix, nodal_deformation);

  // add bodyforces
  noalias(rRightHandSideVector) += this->CalculateBodyForces();
  KRATOS_CATCH("")
}

void CrBeamElementLinear3D2N::CalculateLeftHandSide(
    MatrixType &rLeftHandSideMatrix, ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY;
  BoundedMatrix<double, msElementSize, msElementSize> transformation_matrix =
      this->CalculateInitialLocalCS();
  rLeftHandSideMatrix = ZeroMatrix(msElementSize, msElementSize);
  noalias(rLeftHandSideMatrix) += this->CreateElementStiffnessMatrix_Material();

  //// start static condensation
  if (this->Has(CONDENSED_DOF_LIST)) {
    Vector dof_list_input = this->GetValue(CONDENSED_DOF_LIST);
    std::vector<int> dofList(dof_list_input.size());
    for (SizeType i = 0; i < dof_list_input.size(); ++i)
      dofList[i] = dof_list_input[i];
    StaticCondensationUtility::CondenseLeftHandSide(*this, rLeftHandSideMatrix,
                                                    dofList);
  }
  //// end static condensation

  BoundedMatrix<double, msElementSize, msElementSize> aux_matrix =
      ZeroMatrix(msElementSize);
  aux_matrix = prod(transformation_matrix, rLeftHandSideMatrix);
  noalias(rLeftHandSideMatrix) =
      prod(aux_matrix, Matrix(trans(transformation_matrix)));

  KRATOS_CATCH("")
}

void CrBeamElementLinear3D2N::CalculateMassMatrix(MatrixType &rMassMatrix,
                                            ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;
  if (rMassMatrix.size1() != msElementSize) {
    rMassMatrix.resize(msElementSize, msElementSize, false);
  }
  rMassMatrix = ZeroMatrix(msElementSize, msElementSize);

  bool use_consistent_mass_matrix = false;

  if (this->GetProperties().Has(USE_CONSISTENT_MASS_MATRIX)) {
    use_consistent_mass_matrix = GetProperties()[USE_CONSISTENT_MASS_MATRIX];
  }

  if (!use_consistent_mass_matrix) {
    this->CalculateLumpedMassMatrix(rMassMatrix, rCurrentProcessInfo);
  } else {
    this->CalculateConsistentMassMatrix(rMassMatrix, rCurrentProcessInfo);
    BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix =
        this->CalculateInitialLocalCS();

    BoundedMatrix<double, msElementSize, msElementSize> aux_matrix =
        prod(rotation_matrix, rMassMatrix);
    rMassMatrix = prod(aux_matrix, Matrix(trans(rotation_matrix)));
  }

  KRATOS_CATCH("")
}


BoundedMatrix<double, CrBeamElement3D2N::msLocalSize,
               CrBeamElement3D2N::msLocalSize>
CrBeamElementLinear3D2N::CalculateDeformationStiffness() {

  KRATOS_TRY
  BoundedMatrix<double, msLocalSize, msLocalSize> Kd =
      ZeroMatrix(msLocalSize, msLocalSize);
  const double E = this->GetProperties()[YOUNG_MODULUS];
  const double G = this->CalculateShearModulus();
  const double A = this->GetProperties()[CROSS_AREA];
  const double L = this->CalculateReferenceLength();

  const double J = this->GetProperties()[TORSIONAL_INERTIA];
  const double Iy = this->GetProperties()[I22];
  const double Iz = this->GetProperties()[I33];

  double Ay = 0.00;
  if (this->GetProperties().Has(AREA_EFFECTIVE_Y)) {
    Ay = GetProperties()[AREA_EFFECTIVE_Y];
  }

  double Az = 0.00;
  if (this->GetProperties().Has(AREA_EFFECTIVE_Z)) {
    Az = GetProperties()[AREA_EFFECTIVE_Z];
  }
  const double Psi_y = this->CalculatePsi(Iy, Az);
  const double Psi_z = this->CalculatePsi(Iz, Ay);

  Kd(0, 0) = G * J / L;
  Kd(1, 1) = E * Iy / L;
  Kd(2, 2) = E * Iz / L;
  Kd(3, 3) = E * A / L;
  Kd(4, 4) = 3.0 * E * Iy * Psi_y / L;
  Kd(5, 5) = 3.0 * E * Iz * Psi_z / L;

  return Kd;
  KRATOS_CATCH("")
}

void CrBeamElementLinear3D2N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>> &rVariable,
    std::vector<array_1d<double, 3>> &rOutput,
    const ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  // element with two nodes can only represent results at one node
  const unsigned int &write_points_number =
      GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);
  if (rOutput.size() != write_points_number) {
    rOutput.resize(write_points_number);
  }

  Matrix left_hand_side_matrix = CreateElementStiffnessMatrix_Material();

  Vector nodal_deformation = ZeroVector(msElementSize);
  this->GetValuesVector(nodal_deformation);

  BoundedMatrix<double, msElementSize, msElementSize> transformation_matrix =
      this->CalculateInitialLocalCS();
  nodal_deformation =
      prod(Matrix(trans(transformation_matrix)), nodal_deformation);

  //// start static back condensation
  if (this->Has(CONDENSED_DOF_LIST)) {
    Vector dof_list_input = this->GetValue(CONDENSED_DOF_LIST);
    std::vector<int> dofList(dof_list_input.size());
    for (SizeType i = 0; i < dof_list_input.size(); ++i)
      dofList[i] = dof_list_input[i];
    Vector nodal_deformation_temp = nodal_deformation;
    StaticCondensationUtility::ConvertingCondensation(
        *this, nodal_deformation_temp, nodal_deformation, dofList,
        left_hand_side_matrix);
  }
  //// end static back condensation

  Vector stress = prod(left_hand_side_matrix, nodal_deformation);

  // rOutput[GP 1,2,3][x,y,z]

  if (rVariable == MOMENT) {
    rOutput[0][0] = -1.0 * stress[3] * 0.75 + stress[9] * 0.25;
    rOutput[1][0] = -1.0 * stress[3] * 0.50 + stress[9] * 0.50;
    rOutput[2][0] = -1.0 * stress[3] * 0.25 + stress[9] * 0.75;

    rOutput[0][1] = -1.0 * stress[4] * 0.75 + stress[10] * 0.25;
    rOutput[1][1] = -1.0 * stress[4] * 0.50 + stress[10] * 0.50;
    rOutput[2][1] = -1.0 * stress[4] * 0.25 + stress[10] * 0.75;

    rOutput[0][2] = 1.0 * stress[5] * 0.75 - stress[11] * 0.25;
    rOutput[1][2] = 1.0 * stress[5] * 0.50 - stress[11] * 0.50;
    rOutput[2][2] = 1.0 * stress[5] * 0.25 - stress[11] * 0.75;
  }
  if (rVariable == FORCE) {
    rOutput[0][0] = -1.0 * stress[0] * 0.75 + stress[6] * 0.25;
    rOutput[1][0] = -1.0 * stress[0] * 0.50 + stress[6] * 0.50;
    rOutput[2][0] = -1.0 * stress[0] * 0.25 + stress[6] * 0.75;

    rOutput[0][1] = -1.0 * stress[1] * 0.75 + stress[7] * 0.25;
    rOutput[1][1] = -1.0 * stress[1] * 0.50 + stress[7] * 0.50;
    rOutput[2][1] = -1.0 * stress[1] * 0.25 + stress[7] * 0.75;

    rOutput[0][2] = -1.0 * stress[2] * 0.75 + stress[8] * 0.25;
    rOutput[1][2] = -1.0 * stress[2] * 0.50 + stress[8] * 0.50;
    rOutput[2][2] = -1.0 * stress[2] * 0.25 + stress[8] * 0.75;
  }

  KRATOS_CATCH("")
}

void CrBeamElementLinear3D2N::CalculateOnIntegrationPoints(
    const Variable<Vector> &rVariable, std::vector<Vector> &rOutput,
    const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  if (rVariable == LOCAL_AXES_VECTOR) {
    BoundedMatrix<double, msElementSize, msElementSize> transformation_matrix;
    transformation_matrix = this->CalculateInitialLocalCS();

    rOutput.resize(3);
    for (int i = 0; i < 3; ++i)
      rOutput[i] = ZeroVector(3);

    for (SizeType i = 0; i < 3; ++i) {
      for (SizeType j = 0; j < 3; ++j) {
        rOutput[i][j] = transformation_matrix(j, i);
      }
    }
  }

  KRATOS_CATCH("");
}

void CrBeamElementLinear3D2N::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, CrBeamElement3D2N);
}

void CrBeamElementLinear3D2N::load(Serializer &rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, CrBeamElement3D2N);
}

} // namespace Kratos.
