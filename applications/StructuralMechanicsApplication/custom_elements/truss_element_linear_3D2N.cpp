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
#include "custom_elements/truss_element_linear_3D2N.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos {
TrussElementLinear3D2N::TrussElementLinear3D2N(IndexType NewId,
                                               GeometryType::Pointer pGeometry)
    : TrussElement3D2N(NewId, pGeometry) {}

TrussElementLinear3D2N::TrussElementLinear3D2N(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : TrussElement3D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
TrussElementLinear3D2N::Create(IndexType NewId,
                               NodesArrayType const &rThisNodes,
                               PropertiesType::Pointer pProperties) const {
  const GeometryType &rGeom = this->GetGeometry();
  return Kratos::make_shared<TrussElementLinear3D2N>(
      NewId, rGeom.Create(rThisNodes), pProperties);
}

Element::Pointer
TrussElementLinear3D2N::Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                               PropertiesType::Pointer pProperties) const {
  return Kratos::make_shared<TrussElementLinear3D2N>(
      NewId, pGeom, pProperties);
}

TrussElementLinear3D2N::~TrussElementLinear3D2N() {}

BoundedMatrix<double, TrussElement3D2N::msLocalSize,
               TrussElement3D2N::msLocalSize>
TrussElementLinear3D2N::CreateElementStiffnessMatrix(
    ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  BoundedMatrix<double, msLocalSize, msLocalSize> LocalStiffnessMatrix =
      ZeroMatrix(msLocalSize, msLocalSize);
  this->CalculateElasticStiffnessMatrix(LocalStiffnessMatrix,
                                        rCurrentProcessInfo);
  return LocalStiffnessMatrix;
  KRATOS_CATCH("")
}

void TrussElementLinear3D2N::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  // resizing the matrices + create memory for LHS
  rLeftHandSideMatrix = ZeroMatrix(msLocalSize, msLocalSize);
  // creating LHS
  rLeftHandSideMatrix = this->CreateElementStiffnessMatrix(rCurrentProcessInfo);

  Vector nodal_deformation = ZeroVector(msLocalSize);
  this->GetValuesVector(nodal_deformation, 0);
  rRightHandSideVector = ZeroVector(msLocalSize);
  noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, nodal_deformation);
  noalias(rRightHandSideVector) += this->CalculateBodyForces();
  this->AddPrestressLinear(rRightHandSideVector);

  KRATOS_CATCH("")
}

void TrussElementLinear3D2N::AddPrestressLinear(
    VectorType &rRightHandSideVector) {
  KRATOS_TRY;
  BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
      ZeroMatrix(msLocalSize, msLocalSize);
  this->CreateTransformationMatrix(transformation_matrix);
  double prestress = 0.00;
  if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
    prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
  }
  const double A = this->GetProperties()[CROSS_AREA];
  const double N = prestress * A;

  // internal force vectors
  BoundedVector<double, msLocalSize> f_local = ZeroVector(msLocalSize);
  f_local[0] = -1.00 * N;
  f_local[3] = 1.00 * N;
  rRightHandSideVector -= prod(transformation_matrix, f_local);
  KRATOS_CATCH("")
}

void TrussElementLinear3D2N::CalculateRightHandSide(
    VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY

  Matrix left_hand_side_matrix = ZeroMatrix(msLocalSize, msLocalSize);
  this->CalculateLeftHandSide(left_hand_side_matrix, rCurrentProcessInfo);
  Vector nodal_deformation = ZeroVector(msLocalSize);
  this->GetValuesVector(nodal_deformation);
  rRightHandSideVector = ZeroVector(msLocalSize);
  noalias(rRightHandSideVector) -=
      prod(left_hand_side_matrix, nodal_deformation);
  this->AddPrestressLinear(rRightHandSideVector);

  // add bodyforces
  noalias(rRightHandSideVector) += this->CalculateBodyForces();
  KRATOS_CATCH("")
}

void TrussElementLinear3D2N::CalculateLeftHandSide(
    MatrixType &rLeftHandSideMatrix, ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  // resizing the matrices + create memory for LHS
  rLeftHandSideMatrix = ZeroMatrix(msLocalSize, msLocalSize);
  // creating LHS
  rLeftHandSideMatrix = this->CreateElementStiffnessMatrix(rCurrentProcessInfo);
  KRATOS_CATCH("")
}

void TrussElementLinear3D2N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>> &rVariable,
    std::vector<array_1d<double, 3>> &rOutput,
    const ProcessInfo &rCurrentProcessInfo) {

  const GeometryType::IntegrationPointsArrayType &integration_points =
      GetGeometry().IntegrationPoints();
  if (rOutput.size() != integration_points.size()) {
    rOutput.resize(integration_points.size());
  }

  if (rVariable == FORCE) {
    BoundedVector<double, msDimension> truss_forces = ZeroVector(msDimension);
    truss_forces[2] = 0.00;
    truss_forces[1] = 0.00;
    const double A = this->GetProperties()[CROSS_AREA];

    double prestress = 0.00;
    if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
      prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
    }

    Matrix left_hand_side_matrix = ZeroMatrix(msLocalSize, msLocalSize);
    ProcessInfo
        dummy_info; // CalculateLeftHandSide does not take const ProcessInfo
    this->CalculateLeftHandSide(left_hand_side_matrix, dummy_info);
    Vector nodal_deformation = ZeroVector(msLocalSize);
    this->GetValuesVector(nodal_deformation);
    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    this->CreateTransformationMatrix(transformation_matrix);
    Vector f_int = prod(left_hand_side_matrix, nodal_deformation);
    f_int = prod(Matrix(trans(transformation_matrix)), f_int);
    truss_forces[0] = f_int[3] + prestress * A;

    rOutput[0] = truss_forces;
  }
}

void TrussElementLinear3D2N::WriteTransformationCoordinates(
    BoundedVector<double, TrussElement3D2N::msLocalSize>
        &rReferenceCoordinates) {
  KRATOS_TRY;
  rReferenceCoordinates = ZeroVector(msLocalSize);
  rReferenceCoordinates[0] = this->GetGeometry()[0].X0();
  rReferenceCoordinates[1] = this->GetGeometry()[0].Y0();
  rReferenceCoordinates[2] = this->GetGeometry()[0].Z0();
  rReferenceCoordinates[3] = this->GetGeometry()[1].X0();
  rReferenceCoordinates[4] = this->GetGeometry()[1].Y0();
  rReferenceCoordinates[5] = this->GetGeometry()[1].Z0();
  KRATOS_CATCH("");
}

void TrussElementLinear3D2N::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TrussElement3D2N);
}
void TrussElementLinear3D2N::load(Serializer &rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TrussElement3D2N);
}

} // namespace Kratos.
