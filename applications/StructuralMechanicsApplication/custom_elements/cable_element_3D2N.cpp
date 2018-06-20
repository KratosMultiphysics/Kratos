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
#include "custom_elements/cable_element_3D2N.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos {
CableElement3D2N::CableElement3D2N(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
    : TrussElement3D2N(NewId, pGeometry) {}

CableElement3D2N::CableElement3D2N(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
    : TrussElement3D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
CableElement3D2N::Create(IndexType NewId, NodesArrayType const &rThisNodes,
                         PropertiesType::Pointer pProperties) const {
  const GeometryType &rGeom = this->GetGeometry();
  return Kratos::make_shared<CableElement3D2N>(NewId, rGeom.Create(rThisNodes),
                                               pProperties);
}

Element::Pointer
CableElement3D2N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                         PropertiesType::Pointer pProperties) const {
  return Kratos::make_shared<CableElement3D2N>(NewId, pGeom,
                                               pProperties);
}

CableElement3D2N::~CableElement3D2N() {}

BoundedMatrix<double, TrussElement3D2N::msLocalSize,
               TrussElement3D2N::msLocalSize>
CableElement3D2N::CreateElementStiffnessMatrix(
    ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  BoundedMatrix<double, msLocalSize, msLocalSize> local_stiffness_matrix =
      ZeroMatrix(msLocalSize, msLocalSize);

  if (this->mIsCompressed) {
    local_stiffness_matrix = ZeroMatrix(msLocalSize, msLocalSize);
  }

  else {
    this->CalculateElasticStiffnessMatrix(local_stiffness_matrix,
                                          rCurrentProcessInfo);
    BoundedMatrix<double, msLocalSize, msLocalSize> K_geo =
        ZeroMatrix(msLocalSize, msLocalSize);
    this->CalculateGeometricStiffnessMatrix(K_geo, rCurrentProcessInfo);

    local_stiffness_matrix += K_geo;
  }

  return local_stiffness_matrix;
  KRATOS_CATCH("")
}

void CableElement3D2N::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                            VectorType &rRightHandSideVector,
                                            ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  // calculate internal forces
  BoundedVector<double, msLocalSize> internal_forces = ZeroVector(msLocalSize);
  this->UpdateInternalForces(internal_forces);
  // resizing the matrices + create memory for LHS
  rLeftHandSideMatrix = ZeroMatrix(msLocalSize, msLocalSize);
  // creating LHS
  noalias(rLeftHandSideMatrix) =
      this->CreateElementStiffnessMatrix(rCurrentProcessInfo);

  if (this->mIsCompressed) {
    rRightHandSideVector = ZeroVector(msLocalSize);
  } else {
    // create+compute RHS
    rRightHandSideVector = ZeroVector(msLocalSize);
    // update Residual
    noalias(rRightHandSideVector) -= internal_forces;
    // add bodyforces
    noalias(rRightHandSideVector) += this->CalculateBodyForces();
  }

  KRATOS_CATCH("")
}

void CableElement3D2N::CalculateRightHandSide(
    VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  BoundedVector<double, msLocalSize> internal_forces = ZeroVector(msLocalSize);
  this->UpdateInternalForces(internal_forces);

  if (this->mIsCompressed) {
    rRightHandSideVector = ZeroVector(msLocalSize);
  } else {
    rRightHandSideVector = ZeroVector(msLocalSize);
    noalias(rRightHandSideVector) -= internal_forces;
    // add bodyforces
    noalias(rRightHandSideVector) += this->CalculateBodyForces();
  }

  KRATOS_CATCH("")
}

void CableElement3D2N::UpdateInternalForces(
    BoundedVector<double, TrussElement3D2N::msLocalSize> &rinternalForces) {

  KRATOS_TRY
  BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
      ZeroMatrix(msLocalSize, msLocalSize);

  this->CreateTransformationMatrix(transformation_matrix);
  const double internal_strain_green_lagrange =
      this->CalculateGreenLagrangeStrain();
  const double l = this->CalculateCurrentLength();
  const double L0 = this->CalculateReferenceLength();
  const double E = this->GetProperties()[YOUNG_MODULUS];
  const double A = this->GetProperties()[CROSS_AREA];

  double prestress = 0.00;
  if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
    prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
  }

  const double normal_force =
      ((E * internal_strain_green_lagrange + prestress) * l * A) / L0;

  if (normal_force < 0.00)
    this->mIsCompressed = true;
  else
    this->mIsCompressed = false;

  // internal force vectors
  BoundedVector<double, msLocalSize> f_local = ZeroVector(msLocalSize);
  f_local[0] = -1.00 * normal_force;
  f_local[3] = 1.00 * normal_force;
  rinternalForces = ZeroVector(msLocalSize);
  noalias(rinternalForces) = prod(transformation_matrix, f_local);
  KRATOS_CATCH("");
}

void CableElement3D2N::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TrussElement3D2N);
  rSerializer.save("mIscompressed", this->mIsCompressed);
}
void CableElement3D2N::load(Serializer &rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TrussElement3D2N);
  rSerializer.load("mIscompressed", this->mIsCompressed);
}
} // namespace Kratos.
