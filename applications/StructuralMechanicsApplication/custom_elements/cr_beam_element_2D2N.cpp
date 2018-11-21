// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//  license: 	 structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/cr_beam_element_2D2N.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos {

CrBeamElement2D2N::CrBeamElement2D2N(IndexType NewId,
                                     GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

CrBeamElement2D2N::CrBeamElement2D2N(IndexType NewId,
                                     GeometryType::Pointer pGeometry,
                                     PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
CrBeamElement2D2N::Create(IndexType NewId, NodesArrayType const &rThisNodes,
                          PropertiesType::Pointer pProperties) const {
  const GeometryType &rGeom = this->GetGeometry();
  return Kratos::make_shared<CrBeamElement2D2N>(NewId, rGeom.Create(rThisNodes),
                                                pProperties);
}

Element::Pointer
CrBeamElement2D2N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                          PropertiesType::Pointer pProperties) const {
  return Kratos::make_shared<CrBeamElement2D2N>(NewId, pGeom,
                                                pProperties);
}

CrBeamElement2D2N::~CrBeamElement2D2N() {}

void CrBeamElement2D2N::EquationIdVector(EquationIdVectorType &rResult,
                                         ProcessInfo &rCurrentProcessInfo) {
  if (rResult.size() != msElementSize)
    rResult.resize(msElementSize);

  for (int i = 0; i < msNumberOfNodes; ++i) {
    int index = i * msLocalSize;
    rResult[index] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
    rResult[index + 1] =
        this->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();

    rResult[index + 2] = this->GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
  }
}

void CrBeamElement2D2N::GetDofList(DofsVectorType &rElementalDofList,
                                   ProcessInfo &rCurrentProcessInfo) {

  if (rElementalDofList.size() != msElementSize) {
    rElementalDofList.resize(msElementSize);
  }

  for (int i = 0; i < msNumberOfNodes; ++i) {
    int index = i * msLocalSize;
    rElementalDofList[index] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
    rElementalDofList[index + 1] =
        this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);

    rElementalDofList[index + 2] = this->GetGeometry()[i].pGetDof(ROTATION_Z);
  }
}

void CrBeamElement2D2N::Initialize() {

  KRATOS_TRY;
  KRATOS_CATCH("")
}

void CrBeamElement2D2N::GetValuesVector(Vector &rValues, int Step) {

  KRATOS_TRY
  if (rValues.size() != msElementSize)
    rValues.resize(msElementSize, false);

  for (int i = 0; i < msNumberOfNodes; ++i) {
    int index = i * msLocalSize;
    rValues[index] =
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
    rValues[index + 1] =
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
    rValues[index + 2] =
        this->GetGeometry()[i].FastGetSolutionStepValue(ROTATION_Z, Step);
  }
  KRATOS_CATCH("")
}

void CrBeamElement2D2N::GetFirstDerivativesVector(Vector &rValues, int Step) {

  KRATOS_TRY
  if (rValues.size() != msElementSize)
    rValues.resize(msElementSize, false);

  for (int i = 0; i < msNumberOfNodes; ++i) {
    int index = i * msLocalSize;
    rValues[index] =
        this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X, Step);
    rValues[index + 1] =
        this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
    rValues[index + 2] = this->GetGeometry()[i].FastGetSolutionStepValue(
        ANGULAR_VELOCITY_Z, Step);
  }

  KRATOS_CATCH("")
}

void CrBeamElement2D2N::GetSecondDerivativesVector(Vector &rValues, int Step) {
  KRATOS_TRY
  if (rValues.size() != msElementSize)
    rValues.resize(msElementSize, false);

  for (int i = 0; i < msNumberOfNodes; ++i) {
    int index = i * msLocalSize;

    rValues[index] =
        this->GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
    rValues[index + 1] =
        this->GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
    rValues[index + 2] = this->GetGeometry()[i].FastGetSolutionStepValue(
        ANGULAR_ACCELERATION_Z, Step);
  }
  KRATOS_CATCH("")
}

void CrBeamElement2D2N::CalculateMassMatrix(MatrixType &rMassMatrix,
                                            ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;
  if (rMassMatrix.size1() != msElementSize) {
    rMassMatrix.resize(msElementSize, msElementSize, false);
  }
  rMassMatrix = ZeroMatrix(msElementSize, msElementSize);

  const double L = this->CalculateLength();
  const double A = this->GetProperties()[CROSS_AREA];
  const double rho = this->GetProperties()[DENSITY];

  const double pre_beam = (rho * A * L) / 420.00;
  const double pre_bar = (rho * A * L) / 6.00;

  // bar part
  rMassMatrix(0, 0) = 2.00 * pre_bar;
  rMassMatrix(0, 3) = 1.00 * pre_bar;
  rMassMatrix(3, 0) = 1.00 * pre_bar;
  rMassMatrix(3, 3) = 2.00 * pre_bar;

  // beam part

  rMassMatrix(1, 1) = pre_beam * 156.00;
  rMassMatrix(1, 2) = pre_beam * 22.00 * L;
  rMassMatrix(1, 4) = pre_beam * 54.00;
  rMassMatrix(1, 5) = pre_beam * (-13.00) * L;

  rMassMatrix(2, 1) = pre_beam * 22.00 * L;
  rMassMatrix(2, 2) = pre_beam * 4.00 * L * L;
  rMassMatrix(2, 4) = pre_beam * 13.00 * L;
  rMassMatrix(2, 5) = pre_beam * (-3.00) * L * L;

  rMassMatrix(4, 1) = pre_beam * 54.00;
  rMassMatrix(4, 2) = pre_beam * 13.00 * L;
  rMassMatrix(4, 4) = pre_beam * 156.00;
  rMassMatrix(4, 5) = pre_beam * (-22.00) * L;

  rMassMatrix(5, 1) = pre_beam * (-13.00) * L;
  rMassMatrix(5, 2) = pre_beam * (-3.00) * L * L;
  rMassMatrix(5, 4) = pre_beam * (-22.00) * L;
  rMassMatrix(5, 5) = pre_beam * (4.00) * L * L;

  this->GlobalizeMatrix(rMassMatrix);

  KRATOS_CATCH("")
}

void CrBeamElement2D2N::CalculateDampingMatrix(
    MatrixType &rDampingMatrix, ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  if (rDampingMatrix.size1() != msElementSize) {
    rDampingMatrix.resize(msElementSize, msElementSize, false);
  }

  rDampingMatrix = ZeroMatrix(msElementSize, msElementSize);

  Matrix stiffness_matrix = ZeroMatrix(msElementSize, msElementSize);

  this->CalculateLeftHandSide(stiffness_matrix, rCurrentProcessInfo);

  Matrix mass_matrix = ZeroMatrix(msElementSize, msElementSize);

  this->CalculateMassMatrix(mass_matrix, rCurrentProcessInfo);

  double alpha = 0.0;
  if (this->GetProperties().Has(RAYLEIGH_ALPHA)) {
    alpha = this->GetProperties()[RAYLEIGH_ALPHA];
  } else if (rCurrentProcessInfo.Has(RAYLEIGH_ALPHA)) {
    alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
  }

  double beta = 0.0;
  if (this->GetProperties().Has(RAYLEIGH_BETA)) {
    beta = this->GetProperties()[RAYLEIGH_BETA];
  } else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA)) {
    beta = rCurrentProcessInfo[RAYLEIGH_BETA];
  }

  noalias(rDampingMatrix) += alpha * mass_matrix;
  noalias(rDampingMatrix) += beta * stiffness_matrix;

  KRATOS_CATCH("")
}

void CrBeamElement2D2N::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                             VectorType &rRightHandSideVector,
                                             ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;
  // t
  this->mDeformationForces = this->CalculateInternalStresses_DeformationModes();

  // qe
  Vector nodal_forces = ZeroVector(msElementSize);
  nodal_forces = this->ReturnElementForces_Local();
  // q
  this->GlobalizeVector(nodal_forces);
  this->mInternalGlobalForces = nodal_forces;
  // Kt
  this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

  // residual >>> r = f_ext - f_int
  rRightHandSideVector = ZeroVector(msElementSize);
  noalias(rRightHandSideVector) -= nodal_forces;

  noalias(rRightHandSideVector) += this->CalculateBodyForces();

  KRATOS_CATCH("")
}

void CrBeamElement2D2N::CalculateRightHandSide(
    VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;
  rRightHandSideVector = ZeroVector(msElementSize);
  noalias(rRightHandSideVector) -= this->mInternalGlobalForces;
  noalias(rRightHandSideVector) += this->CalculateBodyForces();
  KRATOS_CATCH("")
}

void CrBeamElement2D2N::CalculateLeftHandSide(
    MatrixType &rLeftHandSideMatrix, ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;
  rLeftHandSideMatrix = this->CreateElementStiffnessMatrix_Total();
  this->GlobalizeMatrix(rLeftHandSideMatrix);
  KRATOS_CATCH("")
}

/////////////////////////////////////////////////
///////////// CUSTOM FUNCTIONS --->>
/////////////////////////////////////////////////

BoundedVector<double, CrBeamElement2D2N::msElementSize>
CrBeamElement2D2N::CalculateBodyForces() {
    KRATOS_TRY
    // getting shapefunctionvalues for linear SF
    const Matrix &Ncontainer =
        this->GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

    BoundedVector<double, 3> equivalent_line_load = ZeroVector(3);
    BoundedVector<double, msElementSize> body_forces_global =
        ZeroVector(msElementSize);

    const double A = this->GetProperties()[CROSS_AREA];
    const double l = this->CalculateLength();
    const double rho = this->GetProperties()[DENSITY];

    // Calculating equivalent line load
    if (GetProperties().Has( VOLUME_ACCELERATION ))
        noalias(equivalent_line_load) = A * rho * GetProperties()[VOLUME_ACCELERATION];
    else if (this->Has( VOLUME_ACCELERATION ))
        noalias(equivalent_line_load) = A * rho * this->GetValue(VOLUME_ACCELERATION);
    else if( GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION) ) {
        for (int i = 0; i < msNumberOfNodes; ++i) {
            noalias(equivalent_line_load) +=
                A * rho *
                this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION) *
                Ncontainer(0, i);
        }
    }

    // adding the nodal forces
    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msLocalSize;
        for (int j = 0; j < msDimension; ++j) {
        body_forces_global[j + index] =
            equivalent_line_load[j] * Ncontainer(0, i) * l;
        }
    }

    // adding the nodal moments
    this->CalculateAndAddWorkEquivalentNodalForcesLineLoad(equivalent_line_load,
                                                            body_forces_global, l);
    // return the total ForceVector
    return body_forces_global;
    KRATOS_CATCH("")
}

void CrBeamElement2D2N::CalculateAndAddWorkEquivalentNodalForcesLineLoad(
    const BoundedVector<double, 3> ForceInput,
    BoundedVector<double, CrBeamElement2D2N::msElementSize>
        &rRightHandSideVector,
    const double GeometryLength) {
  KRATOS_TRY;
  // calculate orthogonal load vector
  const double numerical_limit = std::numeric_limits<double>::epsilon();
  Vector geometric_orientation = ZeroVector(3);
  geometric_orientation[0] =
      this->GetGeometry()[1].X() - this->GetGeometry()[0].X();
  geometric_orientation[1] =
      this->GetGeometry()[1].Y() - this->GetGeometry()[0].Y();
  geometric_orientation[2] = 0.000;

  const double vector_norm_a = MathUtils<double>::Norm(geometric_orientation);
  if (vector_norm_a > numerical_limit)
    geometric_orientation /= vector_norm_a;

  Vector line_load_direction = ZeroVector(3);
  for (int i = 0; i < msDimension; ++i) {
    line_load_direction[i] = ForceInput[i];
  }

  const double vector_norm_b = MathUtils<double>::Norm(line_load_direction);
  if (vector_norm_b > numerical_limit)
    line_load_direction /= vector_norm_b;

  double cos_angle = 0.00;
  for (int i = 0; i < msDimension; ++i) {
    cos_angle += line_load_direction[i] * geometric_orientation[i];
  }

  const double sin_angle = std::sqrt(1.00 - (cos_angle * cos_angle));
  const double norm_force_vector_orthogonal = sin_angle * vector_norm_b;

  Vector node_a = ZeroVector(3);
  node_a[0] = this->GetGeometry()[0].X();
  node_a[1] = this->GetGeometry()[0].Y();
  node_a[2] = 0.00;

  Vector node_b = ZeroVector(3);
  node_b = node_a + line_load_direction;

  Vector node_c = ZeroVector(3);
  node_c = node_a + (geometric_orientation * cos_angle);

  Vector load_orthogonal_direction = ZeroVector(3);
  load_orthogonal_direction = node_b - node_c;
  const double vector_norm_c =
      MathUtils<double>::Norm(load_orthogonal_direction);
  if (vector_norm_c > numerical_limit)
    load_orthogonal_direction /= vector_norm_c;

  // now caluclate respective work equivilent nodal moments

  const double custom_moment =
      norm_force_vector_orthogonal * GeometryLength * GeometryLength / 12.00;

  Vector moment_node_a = ZeroVector(3);
  moment_node_a = MathUtils<double>::CrossProduct(geometric_orientation,
                                                  load_orthogonal_direction);
  moment_node_a *= custom_moment;

  rRightHandSideVector[msDimension] += moment_node_a[2];
  rRightHandSideVector[(2 * msDimension) + 1] -= moment_node_a[2];

  KRATOS_CATCH("")
}

double CrBeamElement2D2N::CalculateShearModulus() {
  KRATOS_TRY;
  const double nu = this->GetProperties()[POISSON_RATIO];
  const double E = this->GetProperties()[YOUNG_MODULUS];
  const double G = E / (2.0 * (1.0 + nu));
  return G;
  KRATOS_CATCH("")
}

double CrBeamElement2D2N::CalculatePsi(const double I, const double A_eff) {

  KRATOS_TRY;
  const double E = this->GetProperties()[YOUNG_MODULUS];
  const double L = this->CalculateLength();
  const double G = this->CalculateShearModulus();

  const double phi = (12.0 * E * I) / (L * L * G * A_eff);
  double psi;
  // interpret input A_eff == 0 as shearstiff -> psi = 1.0
  if (A_eff == 0.00)
    psi = 1.00;
  else
    psi = 1.0 / (1.0 + phi);

  return psi;
  KRATOS_CATCH("")
}

double CrBeamElement2D2N::CalculateInitialElementAngle() {
  KRATOS_TRY;
  const double numerical_limit = std::numeric_limits<double>::epsilon();

  const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
  const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();

  const double norm = std::sqrt((dx * dx) + (dy * dy));

  double phi;
  if ((dx > numerical_limit) & (std::abs(dy) < numerical_limit))
    phi = 0.00; // dy = 0 and dx > 0
  else if ((dx < -numerical_limit) & (std::abs(dy) < numerical_limit))
    phi = Globals::Pi; // dy = 0 and dx < 0
  else if (std::abs(dx) < numerical_limit) {
    phi = Globals::Pi / 2.00; // dy > 0 and dx = 0
    if (dy < -numerical_limit)
      phi = 1.500 * Globals::Pi; // dy < 0 and dx = 0
  } else {
    phi = (norm - dx) / dy;
    phi = std::atan(phi);
    phi = 2.00 * phi;
  }

  return phi;
  KRATOS_CATCH("")
}

double CrBeamElement2D2N::CalculateDeformedElementAngle() {
  KRATOS_TRY;
  const double numerical_limit = std::numeric_limits<double>::epsilon();

  Vector current_displacement = ZeroVector(msElementSize);
  this->GetValuesVector(current_displacement, 0);

  const double dx = (this->GetGeometry()[1].X0() + current_displacement[3]) -
                    (this->GetGeometry()[0].X0() + current_displacement[0]);
  const double dy = (this->GetGeometry()[1].Y0() + current_displacement[4]) -
                    (this->GetGeometry()[0].Y0() + current_displacement[1]);

  const double norm = std::sqrt((dx * dx) + (dy * dy));

  double phi;
  if ((dx > numerical_limit) & (std::abs(dy) < numerical_limit))
    phi = 0.00; // dy = 0 and dx > 0
  else if ((dx < -numerical_limit) & (std::abs(dy) < numerical_limit))
    phi = Globals::Pi; // dy = 0 and dx < 0
  else if (std::abs(dx) < numerical_limit) {
    phi = Globals::Pi / 2.00; // dy > 0 and dx = 0
    if (dy < -numerical_limit)
      phi = 1.500 * Globals::Pi; // dy < 0 and dx = 0
  } else {
    phi = (norm - dx) / dy;
    phi = std::atan(phi);
    phi = 2.00 * phi;
  }

  return phi;
  KRATOS_CATCH("")
}

double CrBeamElement2D2N::CalculateLength() {
  KRATOS_TRY;
  const double numerical_limit = std::numeric_limits<double>::epsilon();
  const double du =
      this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) -
      this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
  const double dv =
      this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) -
      this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);

  const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
  const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();

  const double l = std::sqrt((du + dx) * (du + dx) + (dv + dy) * (dv + dy));

  KRATOS_ERROR_IF(l < numerical_limit) << "length 0 for element " << this->Id()
                                       << std::endl;
  return l;
  KRATOS_CATCH("")
}

double CrBeamElement2D2N::CalculateReferenceLength() {
  KRATOS_TRY;
  const double numerical_limit = std::numeric_limits<double>::epsilon();
  const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
  const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
  const double L = std::sqrt((dx * dx) + (dy * dy));

  KRATOS_ERROR_IF(L < numerical_limit) << "length 0 for element " << this->Id()
                                       << std::endl;
  return L;
  KRATOS_CATCH("")
}

BoundedMatrix<double, CrBeamElement2D2N::msElementSize,
               CrBeamElement2D2N::msLocalSize>
CrBeamElement2D2N::CalculateTransformationS() {
  KRATOS_TRY;
  const double L = this->CalculateLength();
  BoundedMatrix<double, msElementSize, msLocalSize> S =
      ZeroMatrix(msElementSize, msLocalSize);
  S(0, 0) = -1.00;
  S(1, 2) = 2.00 / L;
  S(2, 1) = -1.00;
  S(2, 2) = 1.00;
  S(3, 0) = 1.00;
  S(4, 2) = -2.00 / L;
  S(5, 1) = 1.00;
  S(5, 2) = 1.00;
  return S;
  KRATOS_CATCH("")
}

BoundedMatrix<double, CrBeamElement2D2N::msLocalSize,
               CrBeamElement2D2N::msLocalSize>
CrBeamElement2D2N::CreateElementStiffnessMatrix_Kd_mat() {
  KRATOS_TRY
  // element properties
  const double E = this->GetProperties()[YOUNG_MODULUS];
  const double A = this->GetProperties()[CROSS_AREA];
  const double L = this->CalculateLength();

  const double Iz = this->GetProperties()[I33];

  double Ay = 0.00;
  if (this->GetProperties().Has(AREA_EFFECTIVE_Y)) {
    Ay = GetProperties()[AREA_EFFECTIVE_Y];
  }

  const double Psi = this->CalculatePsi(Iz, Ay);

  // element material stiffness matrix
  BoundedMatrix<double, msLocalSize, msLocalSize> Kd =
      ZeroMatrix(msLocalSize, msLocalSize);

  Kd(0, 0) = E * A / L;
  Kd(1, 1) = E * Iz / L;
  Kd(2, 2) = 3.00 * Psi * E * Iz / L;
  return Kd;
  KRATOS_CATCH("")
}

BoundedMatrix<double, CrBeamElement2D2N::msLocalSize,
               CrBeamElement2D2N::msLocalSize>
CrBeamElement2D2N::CreateElementStiffnessMatrix_Kd_geo() {
  KRATOS_TRY
  // element properties
  const double L = this->CalculateLength();
  const double N = this->mDeformationForces[0];

  // element material stiffness matrix
  BoundedMatrix<double, msLocalSize, msLocalSize> Kd =
      ZeroMatrix(msLocalSize, msLocalSize);

  Kd(1, 1) = N * L / 12.00;
  Kd(2, 2) = N * L / 20.00;
  return Kd;
  KRATOS_CATCH("")
}

BoundedMatrix<double, CrBeamElement2D2N::msElementSize,
               CrBeamElement2D2N::msElementSize>
CrBeamElement2D2N::CreateElementStiffnessMatrix_Kr() {
  KRATOS_TRY
  // element properties
  const double L = this->CalculateLength();
  const double N = this->mDeformationForces[0];
  const double Q = (-2.00 / L) * this->mDeformationForces[2];

  // element material stiffness matrix
  BoundedMatrix<double, msElementSize, msElementSize> Kr =
      ZeroMatrix(msElementSize, msElementSize);

  Kr(0, 1) = -Q;
  Kr(0, 4) = Q;
  Kr(1, 0) = -Q;
  Kr(1, 1) = N;
  Kr(1, 3) = Q;
  Kr(1, 4) = -N;

  Kr(3, 1) = Q;
  Kr(3, 4) = -Q;
  Kr(4, 0) = Q;
  Kr(4, 1) = -N;
  Kr(4, 3) = -Q;
  Kr(4, 4) = N;
  return Kr;
  KRATOS_CATCH("")
}

BoundedMatrix<double, CrBeamElement2D2N::msElementSize,
               CrBeamElement2D2N::msElementSize>
CrBeamElement2D2N::CreateElementStiffnessMatrix_Total() {
  KRATOS_TRY
  // co-rotating K
  BoundedMatrix<double, msElementSize, msElementSize> K_r =
      this->CreateElementStiffnessMatrix_Kr();

  // element K (mat+geo)
  BoundedMatrix<double, msLocalSize, msLocalSize> K_d_mat =
      this->CreateElementStiffnessMatrix_Kd_mat();
  BoundedMatrix<double, msLocalSize, msLocalSize> K_d_geo =
      this->CreateElementStiffnessMatrix_Kd_geo();
  BoundedMatrix<double, msLocalSize, msLocalSize> K_d = K_d_mat + K_d_geo;

  BoundedMatrix<double, msElementSize, msLocalSize> S =
      this->CalculateTransformationS();
  BoundedMatrix<double, msElementSize, msElementSize> K_d_element =
      prod(K_d, Matrix(trans(S)));
  K_d_element = prod(S, K_d_element);

  // total K
  BoundedMatrix<double, msElementSize, msElementSize> K_total =
      ZeroMatrix(msElementSize, msElementSize);
  K_total += K_r;
  K_total += K_d_element;

  return K_total;
  KRATOS_CATCH("")
}

BoundedVector<double, CrBeamElement2D2N::msLocalSize>
CrBeamElement2D2N::CalculateDeformationParameters() {
  KRATOS_TRY;
  // calculate v

  Vector current_displacement = ZeroVector(msElementSize);
  this->GetValuesVector(current_displacement, 0);

  BoundedVector<double, msLocalSize> deformation_parameters =
      ZeroVector(msLocalSize);
  deformation_parameters[0] =
      this->CalculateLength() - this->CalculateReferenceLength();
  deformation_parameters[1] = current_displacement[5] - current_displacement[2];
  deformation_parameters[2] = current_displacement[5] + current_displacement[2];
  deformation_parameters[2] -= 2.00 * (this->CalculateDeformedElementAngle() -
                                       this->CalculateInitialElementAngle());

  // calculate modulus 2PI for phi_a
  deformation_parameters[2] =
      this->Modulus2Pi(deformation_parameters[2] + Globals::Pi) - Globals::Pi;

  return deformation_parameters;
  KRATOS_CATCH("")
}

BoundedVector<double, CrBeamElement2D2N::msLocalSize>
CrBeamElement2D2N::CalculateInternalStresses_DeformationModes() {
  KRATOS_TRY;
  // calculate t

  BoundedVector<double, msLocalSize> deformation_stresses =
      ZeroVector(msLocalSize);

  BoundedVector<double, msLocalSize> deformation_modes =
      this->CalculateDeformationParameters();

  BoundedMatrix<double, msLocalSize, msLocalSize> K_d_mat =
      this->CreateElementStiffnessMatrix_Kd_mat();
  BoundedMatrix<double, msLocalSize, msLocalSize> K_d_geo =
      this->CreateElementStiffnessMatrix_Kd_geo();
  BoundedMatrix<double, msLocalSize, msLocalSize> K_d = K_d_mat + K_d_geo;

  deformation_stresses = prod(K_d, deformation_modes);

  return deformation_stresses;
  KRATOS_CATCH("")
}

BoundedMatrix<double, CrBeamElement2D2N::msElementSize,
               CrBeamElement2D2N::msElementSize>
CrBeamElement2D2N::CreateRotationMatrix() {
  KRATOS_TRY;
  const double current_element_angle = this->CalculateDeformedElementAngle();
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

void CrBeamElement2D2N::GlobalizeMatrix(Matrix &A) {
  KRATOS_TRY;
  BoundedMatrix<double, msElementSize, msElementSize> R =
      this->CreateRotationMatrix();

  A = prod(A, Matrix(trans(R)));
  A = prod(R, A);
  KRATOS_CATCH("")
}

void CrBeamElement2D2N::GlobalizeVector(Vector &A) {
  KRATOS_TRY;
  BoundedMatrix<double, msElementSize, msElementSize> R =
      this->CreateRotationMatrix();
  A = prod(R, A);
  KRATOS_CATCH("")
}

void CrBeamElement2D2N::CalculateOnIntegrationPoints(
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

  BoundedMatrix<double, msElementSize, msElementSize> transformation_matrix =
      this->CreateRotationMatrix();
  Vector stress = this->mInternalGlobalForces;
  stress = prod(trans(transformation_matrix), stress);

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

void CrBeamElement2D2N::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>> &rVariable,
    std::vector<array_1d<double, 3>> &rOutput,
    const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;
  this->CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
  KRATOS_CATCH("")
}

CrBeamElement2D2N::IntegrationMethod
CrBeamElement2D2N::GetIntegrationMethod() const {
  // do this to have 3GP as an output in GID
  return Kratos::GeometryData::GI_GAUSS_3;
}

BoundedVector<double, CrBeamElement2D2N::msElementSize>
CrBeamElement2D2N::ReturnElementForces_Local() {
  KRATOS_TRY;
  // calculate qe

  BoundedMatrix<double, msElementSize, msLocalSize> S =
      this->CalculateTransformationS();
  BoundedVector<double, msLocalSize> t =
      this->CalculateInternalStresses_DeformationModes();

  BoundedVector<double, msElementSize> qe = prod(S, t);
  return qe;
  KRATOS_CATCH("")
}

double CrBeamElement2D2N::Modulus2Pi(double A) {
  KRATOS_TRY;
  const int B = A / (2.00 * Globals::Pi);
  const double C = A - (B * 2.00 * Globals::Pi);
  return C;
  KRATOS_CATCH("")
}

void CrBeamElement2D2N::AddExplicitContribution(
    const VectorType &rRHSVector, const Variable<VectorType> &rRHSVariable,
    Variable<array_1d<double, 3>> &rDestinationVariable,
    const ProcessInfo &rCurrentProcessInfo) {
  // FORCE- & Moment- Residual is 3D vector
  KRATOS_TRY;

  BoundedVector<double, msElementSize> damping_residual_contribution =
      ZeroVector(msElementSize);
  // calculate damping contribution to residual -->
  if ((this->GetProperties().Has(RAYLEIGH_ALPHA) ||
       this->GetProperties().Has(RAYLEIGH_BETA)) &&
      (rDestinationVariable != NODAL_INERTIA)) {
    Vector current_nodal_velocities = ZeroVector(msElementSize);
    this->GetFirstDerivativesVector(current_nodal_velocities);
    Matrix damping_matrix = ZeroMatrix(msElementSize, msElementSize);
    ProcessInfo temp_process_information; // cant pass const ProcessInfo
    this->CalculateDampingMatrix(damping_matrix, temp_process_information);
    // current residual contribution due to damping
    noalias(damping_residual_contribution) =
        prod(damping_matrix, current_nodal_velocities);
  }

  if (rRHSVariable == RESIDUAL_VECTOR &&
      rDestinationVariable == FORCE_RESIDUAL) {

    for (SizeType i = 0; i < msNumberOfNodes; ++i) {
      SizeType index = msLocalSize * i;

      GetGeometry()[i].SetLock();

      array_1d<double, 3> &r_force_residual =
          GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

      for (SizeType j = 0; j < msDimension; ++j) {
        r_force_residual[j] +=
            rRHSVector[index + j] - damping_residual_contribution[index + j];
      }

      r_force_residual[msDimension] = 0.00;
      GetGeometry()[i].UnSetLock();
    }
  }

  if (rRHSVariable == RESIDUAL_VECTOR &&
      rDestinationVariable == MOMENT_RESIDUAL) {

    for (SizeType i = 0; i < msNumberOfNodes; ++i) {
      SizeType index = (msLocalSize * i) + msDimension;

      GetGeometry()[i].SetLock();

      array_1d<double, 3> &r_moment_residual =
          GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

      for (SizeType j = 0; j < msDimension; ++j) {
        r_moment_residual[j] = 0.00;
      }
      r_moment_residual[msDimension] +=
          rRHSVector[index] - damping_residual_contribution[index];
      GetGeometry()[i].UnSetLock();
    }
  }

  if (rDestinationVariable == NODAL_INERTIA) {
        Matrix element_mass_matrix = ZeroMatrix(msElementSize, msElementSize);
        ProcessInfo temp_info; // Dummy
        this->CalculateMassMatrix(element_mass_matrix, temp_info);

        for (IndexType i = 0; i < msNumberOfNodes; ++i) {
            double aux_nodal_mass = 0.0;
            double aux_nodal_inertia = 0.0;

            const SizeType index = i * msLocalSize;

            for (IndexType j = 0; j < msElementSize; ++j) {
                aux_nodal_mass += element_mass_matrix(index, j);
                aux_nodal_inertia += element_mass_matrix(index + msDimension, j);
            }

            #pragma omp atomic
            GetGeometry()[i].GetValue(NODAL_MASS) += aux_nodal_mass;

            array_1d<double, 3>& r_nodal_inertia = GetGeometry()[i].GetValue(NODAL_INERTIA);
            #pragma omp atomic
            r_nodal_inertia[msDimension] += std::abs(aux_nodal_inertia);
        }
  }

  KRATOS_CATCH("")
}

int CrBeamElement2D2N::Check(const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY
  const double numerical_limit = std::numeric_limits<double>::epsilon();
  if (GetGeometry().WorkingSpaceDimension() != 2 || GetGeometry().size() != 2) {
    KRATOS_ERROR
        << "The beam element works only in 2D and with 2 noded elements"
        << "" << std::endl;
  }
  // verify that the variables are correctly initialized
  if (VELOCITY.Key() == 0) {
    KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is "
                    "correctly registered"
                 << "" << std::endl;
  }
  if (DISPLACEMENT.Key() == 0) {
    KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is "
                    "correctly registered"
                 << "" << std::endl;
  }
  if (ACCELERATION.Key() == 0) {
    KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is "
                    "correctly registered"
                 << "" << std::endl;
  }
  if (DENSITY.Key() == 0) {
    KRATOS_ERROR << "DENSITY has Key zero! (check if the application is "
                    "correctly registered"
                 << "" << std::endl;
  }
  if (CROSS_AREA.Key() == 0) {
    KRATOS_ERROR << "CROSS_AREA has Key zero! (check if the application is "
                    "correctly registered"
                 << "" << std::endl;
  }
  // verify that the dofs exist
  for (unsigned int i = 0; i < this->GetGeometry().size(); ++i) {
    if (this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false) {
      KRATOS_ERROR << "missing variable DISPLACEMENT on node "
                   << this->GetGeometry()[i].Id() << std::endl;
    }
    if (this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false ||
        this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false) {
      KRATOS_ERROR
          << "missing one of the dofs for the variable DISPLACEMENT on node "
          << GetGeometry()[i].Id() << std::endl;
    }
  }

  if (this->GetProperties().Has(CROSS_AREA) == false ||
      this->GetProperties()[CROSS_AREA] <= numerical_limit) {
    KRATOS_ERROR << "CROSS_AREA not provided for this element" << this->Id()
                 << std::endl;
  }

  if (this->GetProperties().Has(YOUNG_MODULUS) == false ||
      this->GetProperties()[YOUNG_MODULUS] <= numerical_limit) {
    KRATOS_ERROR << "YOUNG_MODULUS not provided for this element" << this->Id()
                 << std::endl;
  }
  if (this->GetProperties().Has(DENSITY) == false) {
    KRATOS_ERROR << "DENSITY not provided for this element" << this->Id()
                 << std::endl;
  }

  if (this->GetProperties().Has(POISSON_RATIO) == false) {
    KRATOS_ERROR << "POISSON_RATIO not provided for this element" << this->Id()
                 << std::endl;
  }

  if (this->GetProperties().Has(I33) == false) {
    KRATOS_ERROR << "I33 not provided for this element" << this->Id()
                 << std::endl;
  }
  return 0;

  KRATOS_CATCH("")
}

void CrBeamElement2D2N::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
  rSerializer.save("mDeformationForces", this->mDeformationForces);
  rSerializer.save("GlobalInternalForces", this->mInternalGlobalForces);
}

void CrBeamElement2D2N::load(Serializer &rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
  rSerializer.load("mDeformationForces", this->mDeformationForces);
  rSerializer.load("GlobalInternalForces", this->mInternalGlobalForces);
}
} // namespace Kratos.
