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
#include "custom_elements/truss_element_3D2N.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos {
TrussElement3D2N::TrussElement3D2N(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

TrussElement3D2N::TrussElement3D2N(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
TrussElement3D2N::Create(IndexType NewId, NodesArrayType const &rThisNodes,
                         PropertiesType::Pointer pProperties) const {
  const GeometryType &rGeom = this->GetGeometry();
  return Kratos::make_shared<TrussElement3D2N>(NewId, rGeom.Create(rThisNodes),
                                               pProperties);
}

Element::Pointer
TrussElement3D2N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                         PropertiesType::Pointer pProperties) const {
  return Kratos::make_shared<TrussElement3D2N>(NewId, pGeom,
                                               pProperties);
}

TrussElement3D2N::~TrussElement3D2N() {}

void TrussElement3D2N::EquationIdVector(EquationIdVectorType &rResult,
                                        ProcessInfo &rCurrentProcessInfo) {

  if (rResult.size() != msLocalSize)
    rResult.resize(msLocalSize);

  for (int i = 0; i < msNumberOfNodes; ++i) {
    int index = i * 3;
    rResult[index] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
    rResult[index + 1] =
        this->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[index + 2] =
        this->GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
  }
}
void TrussElement3D2N::GetDofList(DofsVectorType &rElementalDofList,
                                  ProcessInfo &rCurrentProcessInfo) {

  if (rElementalDofList.size() != msLocalSize) {
    rElementalDofList.resize(msLocalSize);
  }

  for (int i = 0; i < msNumberOfNodes; ++i) {
    int index = i * 3;
    rElementalDofList[index] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
    rElementalDofList[index + 1] =
        this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
    rElementalDofList[index + 2] =
        this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
  }
}

void TrussElement3D2N::Initialize() {
  KRATOS_TRY
  if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        this->mConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
    }
  else
    KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
  KRATOS_CATCH("")
}

BoundedMatrix<double, TrussElement3D2N::msLocalSize,
               TrussElement3D2N::msLocalSize>
TrussElement3D2N::CreateElementStiffnessMatrix(
    ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  BoundedMatrix<double, msLocalSize, msLocalSize> local_stiffness_matrix =
      ZeroMatrix(msLocalSize, msLocalSize);

  this->CalculateElasticStiffnessMatrix(local_stiffness_matrix,
                                        rCurrentProcessInfo);
  BoundedMatrix<double, msLocalSize, msLocalSize> K_geo =
      ZeroMatrix(msLocalSize, msLocalSize);
  this->CalculateGeometricStiffnessMatrix(K_geo, rCurrentProcessInfo);

  local_stiffness_matrix += K_geo;

  return local_stiffness_matrix;
  KRATOS_CATCH("")
}

void TrussElement3D2N::CalculateDampingMatrix(
    MatrixType &rDampingMatrix, ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  if (rDampingMatrix.size1() != msLocalSize) {
    rDampingMatrix.resize(msLocalSize, msLocalSize, false);
  }

  rDampingMatrix = ZeroMatrix(msLocalSize, msLocalSize);

  MatrixType stiffness_matrix = ZeroMatrix(msLocalSize, msLocalSize);

  this->CalculateLeftHandSide(stiffness_matrix, rCurrentProcessInfo);

  MatrixType mass_matrix = ZeroMatrix(msLocalSize, msLocalSize);

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

void TrussElement3D2N::CalculateMassMatrix(MatrixType &rMassMatrix,
                                           ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  if (rMassMatrix.size1() != msLocalSize) {
    rMassMatrix.resize(msLocalSize, msLocalSize, false);
  }

  rMassMatrix = ZeroMatrix(msLocalSize, msLocalSize);

  const double A = this->GetProperties()[CROSS_AREA];
  const double L = this->CalculateReferenceLength();
  const double rho = this->GetProperties()[DENSITY];

  const double total_mass = A * L * rho;

  Vector lumping_factor = ZeroVector(msNumberOfNodes);

  lumping_factor = this->GetGeometry().LumpingFactors(lumping_factor);

  for (int i = 0; i < msNumberOfNodes; ++i) {
    double temp = lumping_factor[i] * total_mass;

    for (int j = 0; j < msDimension; ++j) {
      int index = i * msDimension + j;

      rMassMatrix(index, index) = temp;
    }
  }
  KRATOS_CATCH("")
}

BoundedVector<double, TrussElement3D2N::msLocalSize>
TrussElement3D2N::CalculateBodyForces() {

  KRATOS_TRY
  // getting shapefunctionvalues
  const Matrix &Ncontainer =
      this->GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

  // creating necessary values
  const double A = this->GetProperties()[CROSS_AREA];
  const double L = this->CalculateReferenceLength();
  const double rho = this->GetProperties()[DENSITY];

  double total_mass = A * L * rho;
  BoundedVector<double, msDimension> body_forces_node =
      ZeroVector(msDimension);
  BoundedVector<double, msLocalSize> body_forces_global =
      ZeroVector(msLocalSize);

  // assemble global Vector
  for (int i = 0; i < msNumberOfNodes; ++i) {
    body_forces_node =
        total_mass *
        this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION) *
        Ncontainer(0, i);

    for (unsigned int j = 0; j < msDimension; ++j) {
      body_forces_global[(i * msDimension) + j] = body_forces_node[j];
    }
  }

  return body_forces_global;
  KRATOS_CATCH("")
}

void TrussElement3D2N::GetValuesVector(Vector &rValues, int Step) {

  KRATOS_TRY
  if (rValues.size() != msLocalSize)
    rValues.resize(msLocalSize, false);

  for (int i = 0; i < msNumberOfNodes; ++i) {
    int index = i * msDimension;
    const auto &disp =
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);

    rValues[index] = disp[0];
    rValues[index + 1] = disp[1];
    rValues[index + 2] = disp[2];
  }
  KRATOS_CATCH("")
}

void TrussElement3D2N::GetFirstDerivativesVector(Vector &rValues, int Step) {

  KRATOS_TRY
  if (rValues.size() != msLocalSize)
    rValues.resize(msLocalSize, false);

  for (int i = 0; i < msNumberOfNodes; ++i) {
    int index = i * msDimension;
    const auto &vel =
        this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);

    rValues[index] = vel[0];
    rValues[index + 1] = vel[1];
    rValues[index + 2] = vel[2];
  }
  KRATOS_CATCH("")
}

void TrussElement3D2N::GetSecondDerivativesVector(Vector &rValues, int Step) {

  KRATOS_TRY
  if (rValues.size() != msLocalSize)
    rValues.resize(msLocalSize, false);

  for (int i = 0; i < msNumberOfNodes; ++i) {
    int index = i * msDimension;
    const auto &acc =
        this->GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);

    rValues[index] = acc[0];
    rValues[index + 1] = acc[1];
    rValues[index + 2] = acc[2];
  }

  KRATOS_CATCH("")
}

void TrussElement3D2N::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
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

  // create+compute RHS
  rRightHandSideVector = ZeroVector(msLocalSize);
  // update Residual
  noalias(rRightHandSideVector) -= internal_forces;
  // add bodyforces
  noalias(rRightHandSideVector) += this->CalculateBodyForces();

  KRATOS_CATCH("")
}

void TrussElement3D2N::CalculateRightHandSide(
    VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  rRightHandSideVector = ZeroVector(msLocalSize);

  BoundedVector<double, msLocalSize> internal_forces = ZeroVector(msLocalSize);
  this->UpdateInternalForces(internal_forces);
  noalias(rRightHandSideVector) -= internal_forces;

  // add bodyforces
  noalias(rRightHandSideVector) += this->CalculateBodyForces();
  KRATOS_CATCH("")
}

void TrussElement3D2N::CalculateLeftHandSide(MatrixType &rLeftHandSideMatrix,
                                             ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  // resizing the matrices + create memory for LHS
  rLeftHandSideMatrix = ZeroMatrix(msLocalSize, msLocalSize);
  // creating LHS
  noalias(rLeftHandSideMatrix) =
      this->CreateElementStiffnessMatrix(rCurrentProcessInfo);
  KRATOS_CATCH("")
}

void TrussElement3D2N::CalculateOnIntegrationPoints(
    const Variable<double> &rVariable, std::vector<double> &rOutput,
    const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY
  const GeometryType::IntegrationPointsArrayType &integration_points =
      GetGeometry().IntegrationPoints();

  if (rOutput.size() != integration_points.size()) {
    rOutput.resize(integration_points.size());
  }
  if (rVariable == TRUSS_PRESTRESS_PK2) {
    rOutput[0] = 0.00;
    if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
      rOutput[0] = this->GetProperties()[TRUSS_PRESTRESS_PK2];
    }
  }
  KRATOS_CATCH("")
}

void TrussElement3D2N::CalculateOnIntegrationPoints(
    const Variable<Vector> &rVariable, std::vector<Vector> &rOutput,
    const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY
  const GeometryType::IntegrationPointsArrayType &integration_points =
      GetGeometry().IntegrationPoints();
  if (rOutput.size() != integration_points.size()) {
    rOutput.resize(integration_points.size());
  }
  if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
    Vector Strain = ZeroVector(msDimension);
    Strain[0] = this->CalculateGreenLagrangeStrain();
    Strain[1] = 0.00;
    Strain[2] = 0.00;
    rOutput[0] = Strain;
  }
  KRATOS_CATCH("")
}

void TrussElement3D2N::CalculateOnIntegrationPoints(
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

    const double internal_strain_green_lagrange =
        this->CalculateGreenLagrangeStrain();
    const double L0 = this->CalculateReferenceLength();
    const double l = this->CalculateCurrentLength();
    const double E = this->GetProperties()[YOUNG_MODULUS];

    truss_forces[0] =
        ((E * internal_strain_green_lagrange + prestress) * l * A) / L0;

    rOutput[0] = truss_forces;
  }
}

void TrussElement3D2N::GetValueOnIntegrationPoints(
    const Variable<double> &rVariable, std::vector<double> &rValues,
    const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY
  this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
  KRATOS_CATCH("")
}
void TrussElement3D2N::GetValueOnIntegrationPoints(
    const Variable<Vector> &rVariable, std::vector<Vector> &rValues,
    const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY
  this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
  KRATOS_CATCH("")
}

void TrussElement3D2N::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>> &rVariable,
    std::vector<array_1d<double, 3>> &rOutput,
    const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;
  this->CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
  KRATOS_CATCH("")
}

int TrussElement3D2N::Check(const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY
  const double numerical_limit = std::numeric_limits<double>::epsilon();


  this->mConstitutiveLaw->Check(this->GetProperties(),this->GetGeometry(),rCurrentProcessInfo);

  if (this->GetGeometry().WorkingSpaceDimension() != msDimension ||
      this->GetGeometry().PointsNumber() != msNumberOfNodes) {
    KRATOS_ERROR << "The truss element works only in 3D and with 2 noded elements" << std::endl;
      }
  // verify that the variables are correctly initialized
  if (VELOCITY.Key() == 0)
    KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is "
                    "correctly registered"
                 << "" << std::endl;
  if (DISPLACEMENT.Key() == 0)
    KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is "
                    "correctly registered"
                 << "" << std::endl;
  if (ACCELERATION.Key() == 0)
    KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is "
                    "correctly registered"
                 << "" << std::endl;
  if (DENSITY.Key() == 0)
    KRATOS_ERROR << "DENSITY has Key zero! (check if the application is "
                    "correctly registered"
                 << "" << std::endl;
  if (CROSS_AREA.Key() == 0)
    KRATOS_ERROR << "CROSS_AREA has Key zero! (check if the application is "
                    "correctly registered"
                 << "" << std::endl;
  // verify that the dofs exist
  for (unsigned int i = 0; i < this->GetGeometry().PointsNumber(); ++i) {
    if (this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false)
      KRATOS_ERROR << "missing variable DISPLACEMENT on node "
                   << this->GetGeometry()[i].Id() << std::endl;
    if (this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false ||
        this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false ||
        this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false)
      KRATOS_ERROR
          << "missing one of the dofs for the variable DISPLACEMENT on node "
          << GetGeometry()[i].Id() << std::endl;
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

  return 0;

  KRATOS_CATCH("")
}
double TrussElement3D2N::CalculateGreenLagrangeStrain() {

  KRATOS_TRY
  const double l = this->CalculateCurrentLength();
  const double L = this->CalculateReferenceLength();
  const double e = ((l * l - L * L) / (2.00 * L * L));
  return e;
  KRATOS_CATCH("")
}

double TrussElement3D2N::CalculateReferenceLength() {

  KRATOS_TRY;
  const double numerical_limit = std::numeric_limits<double>::epsilon();
  const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
  const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
  const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
  const double L = std::sqrt(dx * dx + dy * dy + dz * dz);

  KRATOS_ERROR_IF(L<=numerical_limit)
   << "Reference Length of element" << this->Id() << "~ 0" << std::endl;
  return L;
  KRATOS_CATCH("")
}
double TrussElement3D2N::CalculateCurrentLength() {

  KRATOS_TRY;
  const double numerical_limit = std::numeric_limits<double>::epsilon();
  const double du =
      this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) -
      this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
  const double dv =
      this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) -
      this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
  const double dw =
      this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z) -
      this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
  const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
  const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
  const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
  const double l = std::sqrt((du + dx) * (du + dx) + (dv + dy) * (dv + dy) +
                             (dw + dz) * (dw + dz));

  KRATOS_ERROR_IF(l<=numerical_limit)
   << "Current Length of element" << this->Id() << "~ 0" << std::endl;
  return l;
  KRATOS_CATCH("")
}
void TrussElement3D2N::UpdateInternalForces(
    BoundedVector<double, TrussElement3D2N::msLocalSize> &rInternalForces) {

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

  // internal force vectors
  BoundedVector<double, msLocalSize> f_local = ZeroVector(msLocalSize);
  f_local[0] = -1.00 * normal_force;
  f_local[3] = 1.00 * normal_force;
  rInternalForces = ZeroVector(msLocalSize);
  noalias(rInternalForces) = prod(transformation_matrix, f_local);
  KRATOS_CATCH("");
}

void TrussElement3D2N::CreateTransformationMatrix(
    BoundedMatrix<double, TrussElement3D2N::msLocalSize,
                   TrussElement3D2N::msLocalSize> &rRotationMatrix) {

  KRATOS_TRY
  // 1st calculate transformation matrix
  typedef BoundedVector<double, msDimension> arraydim;
  typedef BoundedVector<double, msLocalSize> arraylocal;
  arraydim direction_vector_x = ZeroVector(msDimension);
  arraydim direction_vector_y = ZeroVector(msDimension);
  arraydim direction_vector_z = ZeroVector(msDimension);
  arraylocal reference_coordinates = ZeroVector(msLocalSize);
  arraydim global_z_vector = ZeroVector(msDimension);
  global_z_vector[2] = 1.0;

  this->WriteTransformationCoordinates(reference_coordinates);

  for (unsigned int i = 0; i < msDimension; ++i) {
    direction_vector_x[i] =
        (reference_coordinates[i + msDimension] - reference_coordinates[i]);
  }
  // local x-axis (e1_local) is the beam axis  (in GID is e3_local)
  double VectorNorm;
  VectorNorm = MathUtils<double>::Norm(direction_vector_x);
  if (VectorNorm != 0)
    direction_vector_x /= VectorNorm;

  if (direction_vector_x[2] == 1.00) {
    direction_vector_y[1] = 1.0;
    direction_vector_z[0] = -1.0;
  }

  if (direction_vector_x[2] == -1.00) {
    direction_vector_y[1] = 1.0;
    direction_vector_z[0] = 1.0;
  }

  if (fabs(direction_vector_x[2]) != 1.00) {
    MathUtils<double>::UnitCrossProduct(direction_vector_y, direction_vector_x,
                                        global_z_vector);
    MathUtils<double>::UnitCrossProduct(direction_vector_z, direction_vector_y,
                                        direction_vector_x);
  }

  // 2nd fill big rotation matrix
  BoundedMatrix<double, msDimension, msDimension> current_coordinate_system =
      ZeroMatrix(msDimension, msDimension);
  for (unsigned int i = 0; i < msDimension; ++i) {
    current_coordinate_system(i, 0) = direction_vector_x[i];
    current_coordinate_system(i, 1) = direction_vector_y[i];
    current_coordinate_system(i, 2) = direction_vector_z[i];
  }

  rRotationMatrix = ZeroMatrix(msLocalSize, msLocalSize);
  // Building the rotation matrix for the local element matrix
  for (unsigned int kk = 0; kk < msLocalSize; kk += msDimension) {
    for (unsigned int i = 0; i < msDimension; ++i) {
      for (unsigned int j = 0; j < msDimension; ++j) {
        rRotationMatrix(i + kk, j + kk) = current_coordinate_system(i, j);
      }
    }
  }
  KRATOS_CATCH("")
}

void TrussElement3D2N::WriteTransformationCoordinates(
    BoundedVector<double, TrussElement3D2N::msLocalSize>
        &rReferenceCoordinates) {
  KRATOS_TRY;
  rReferenceCoordinates = ZeroVector(msLocalSize);
  Vector current_displacement = ZeroVector(msLocalSize);
  this->GetValuesVector(current_displacement, 0);
  rReferenceCoordinates[0] =
      this->GetGeometry()[0].X0() + current_displacement[0];
  rReferenceCoordinates[1] =
      this->GetGeometry()[0].Y0() + current_displacement[1];
  rReferenceCoordinates[2] =
      this->GetGeometry()[0].Z0() + current_displacement[2];
  rReferenceCoordinates[3] =
      this->GetGeometry()[1].X0() + current_displacement[3];
  rReferenceCoordinates[4] =
      this->GetGeometry()[1].Y0() + current_displacement[4];
  rReferenceCoordinates[5] =
      this->GetGeometry()[1].Z0() + current_displacement[5];
  KRATOS_CATCH("");
}

void TrussElement3D2N::AddExplicitContribution(
    const VectorType &rRHSVector, const Variable<VectorType> &rRHSVariable,
    Variable<array_1d<double, 3>> &rDestinationVariable,
    const ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  BoundedVector<double, msLocalSize> damping_residual_contribution =
      ZeroVector(msLocalSize);
  // calculate damping contribution to residual -->
  if ((this->GetProperties().Has(RAYLEIGH_ALPHA) ||
       this->GetProperties().Has(RAYLEIGH_BETA)) &&
      (rDestinationVariable != NODAL_INERTIA)) {
    Vector current_nodal_velocities = ZeroVector(msLocalSize);
    this->GetFirstDerivativesVector(current_nodal_velocities);
    Matrix damping_matrix = ZeroMatrix(msLocalSize, msLocalSize);
    ProcessInfo temp_process_information; // cant pass const ProcessInfo
    this->CalculateDampingMatrix(damping_matrix, temp_process_information);
    // current residual contribution due to damping
    noalias(damping_residual_contribution) =
        prod(damping_matrix, current_nodal_velocities);
  }

  if (rRHSVariable == RESIDUAL_VECTOR &&
      rDestinationVariable == FORCE_RESIDUAL) {

    for (size_t i = 0; i < msNumberOfNodes; ++i) {
      size_t index = msDimension * i;
      array_1d<double, 3> &r_force_residual =
          GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
      for (size_t j = 0; j < msDimension; ++j) {
#pragma omp atomic
        r_force_residual[j] +=
            rRHSVector[index + j] - damping_residual_contribution[index + j];
      }
    }
  }

  if (rDestinationVariable == NODAL_INERTIA) {

    Matrix element_mass_matrix = ZeroMatrix(msLocalSize, msLocalSize);
    ProcessInfo temp_info; // Dummy
    this->CalculateMassMatrix(element_mass_matrix, temp_info);

    for (int i = 0; i < msNumberOfNodes; ++i) {
      double &r_nodal_mass = GetGeometry()[i].GetValue(NODAL_MASS);
      array_1d<double, msDimension> &r_nodal_inertia =
          GetGeometry()[i].GetValue(NODAL_INERTIA);
      int index = i * msDimension;

      for (SizeType j = 0; j < msLocalSize; ++j) {
#pragma omp atomic
        r_nodal_mass += element_mass_matrix(index, j);
      }
      for (int k = 0; k < msDimension; ++k) {
#pragma omp atomic
        r_nodal_inertia[k] += 0.00;
      }
    }
  }
  KRATOS_CATCH("")
}

void TrussElement3D2N::CalculateGeometricStiffnessMatrix(
    BoundedMatrix<double, TrussElement3D2N::msLocalSize,
                   TrussElement3D2N::msLocalSize> &rGeometricStiffnessMatrix,
    ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  double E = 0.00;
  ConstitutiveLaw::Parameters Values(this->GetGeometry(),this->GetProperties(),rCurrentProcessInfo);
  this->mConstitutiveLaw->CalculateValue(Values,TANGENT_MODULUS,E);

  const double A = this->GetProperties()[CROSS_AREA];

  double prestress = 0.00;
  if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
    prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
  }

  rGeometricStiffnessMatrix = ZeroMatrix(msLocalSize, msLocalSize);

  // du... delta displacement in x-direction
  // dv... delta displacement in y-direction
  // dw... delta displacement in z-direction
  // L... inital member length
  // l... deformed member length
  // e_gl... green_lagrange strain

  double du = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) -
              this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
  double dv = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) -
              this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
  double dw = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z) -
              this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
  const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
  const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
  const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
  const double L = this->CalculateReferenceLength();
  const double l = this->CalculateCurrentLength();
  double e_gL = (l * l - L * L) / (2.00 * L * L);
  const double L3 = L * L * L;

  const double K_sigma = ((E * A * e_gL) / L) + ((prestress * A) / L);
  const double K_uij = (E * A) / L3;

  rGeometricStiffnessMatrix(0, 0) = K_sigma + K_uij * (2 * du * dx + du * du);
  rGeometricStiffnessMatrix(3, 3) = rGeometricStiffnessMatrix(0, 0);

  rGeometricStiffnessMatrix(1, 1) = K_sigma + K_uij * (2 * dv * dy + dv * dv);
  rGeometricStiffnessMatrix(4, 4) = rGeometricStiffnessMatrix(1, 1);

  rGeometricStiffnessMatrix(2, 2) = K_sigma + K_uij * (2 * dw * dz + dw * dw);
  rGeometricStiffnessMatrix(5, 5) = rGeometricStiffnessMatrix(2, 2);

  rGeometricStiffnessMatrix(0, 1) = K_uij * (dx * dv + dy * du + du * dv);
  rGeometricStiffnessMatrix(1, 0) = rGeometricStiffnessMatrix(0, 1);

  rGeometricStiffnessMatrix(0, 2) = K_uij * (dx * dw + dz * du + du * dw);
  rGeometricStiffnessMatrix(2, 0) = rGeometricStiffnessMatrix(0, 2);

  rGeometricStiffnessMatrix(0, 3) = -rGeometricStiffnessMatrix(0, 0);
  rGeometricStiffnessMatrix(3, 0) = rGeometricStiffnessMatrix(0, 3);

  rGeometricStiffnessMatrix(0, 4) = -rGeometricStiffnessMatrix(0, 1);
  rGeometricStiffnessMatrix(4, 0) = rGeometricStiffnessMatrix(0, 4);

  rGeometricStiffnessMatrix(0, 5) = -rGeometricStiffnessMatrix(0, 2);
  rGeometricStiffnessMatrix(5, 0) = rGeometricStiffnessMatrix(0, 5);

  rGeometricStiffnessMatrix(1, 2) = K_uij * (dy * dw + dz * dv + dv * dw);
  rGeometricStiffnessMatrix(2, 1) = rGeometricStiffnessMatrix(1, 2);

  rGeometricStiffnessMatrix(1, 3) = rGeometricStiffnessMatrix(0, 4);
  rGeometricStiffnessMatrix(3, 1) = rGeometricStiffnessMatrix(1, 3);

  rGeometricStiffnessMatrix(1, 4) = -rGeometricStiffnessMatrix(1, 1);
  rGeometricStiffnessMatrix(4, 1) = rGeometricStiffnessMatrix(1, 4);

  rGeometricStiffnessMatrix(1, 5) = -rGeometricStiffnessMatrix(1, 2);
  rGeometricStiffnessMatrix(5, 1) = rGeometricStiffnessMatrix(1, 5);

  rGeometricStiffnessMatrix(2, 3) = -rGeometricStiffnessMatrix(0, 2);
  rGeometricStiffnessMatrix(3, 2) = rGeometricStiffnessMatrix(2, 3);

  rGeometricStiffnessMatrix(2, 4) = -rGeometricStiffnessMatrix(1, 2);
  rGeometricStiffnessMatrix(4, 2) = rGeometricStiffnessMatrix(2, 4);

  rGeometricStiffnessMatrix(2, 5) = -rGeometricStiffnessMatrix(2, 2);
  rGeometricStiffnessMatrix(5, 2) = rGeometricStiffnessMatrix(2, 5);

  rGeometricStiffnessMatrix(3, 4) = rGeometricStiffnessMatrix(0, 1);
  rGeometricStiffnessMatrix(4, 3) = rGeometricStiffnessMatrix(3, 4);

  rGeometricStiffnessMatrix(3, 5) = rGeometricStiffnessMatrix(0, 2);
  rGeometricStiffnessMatrix(5, 3) = rGeometricStiffnessMatrix(3, 5);

  rGeometricStiffnessMatrix(4, 5) = rGeometricStiffnessMatrix(1, 2);
  rGeometricStiffnessMatrix(5, 4) = rGeometricStiffnessMatrix(4, 5);
  KRATOS_CATCH("")
}

void TrussElement3D2N::CalculateElasticStiffnessMatrix(
    BoundedMatrix<double, TrussElement3D2N::msLocalSize,
                   TrussElement3D2N::msLocalSize> &rElasticStiffnessMatrix,
    ProcessInfo &rCurrentProcessInfo) {
  KRATOS_TRY;

  double E = 0.00;
  ConstitutiveLaw::Parameters Values(this->GetGeometry(),this->GetProperties(),rCurrentProcessInfo);
  this->mConstitutiveLaw->CalculateValue(Values,TANGENT_MODULUS,E);

  double A = this->GetProperties()[CROSS_AREA];

  rElasticStiffnessMatrix = ZeroMatrix(msLocalSize, msLocalSize);

  const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
  const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
  const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
  const double L = this->CalculateReferenceLength();
  const double L3 = L * L * L;

  const double EA = E * A;

  rElasticStiffnessMatrix(0, 0) = (EA * dx * dx) / L3;
  rElasticStiffnessMatrix(3, 3) = rElasticStiffnessMatrix(0, 0);

  rElasticStiffnessMatrix(1, 1) = (EA * dy * dy) / L3;
  rElasticStiffnessMatrix(4, 4) = rElasticStiffnessMatrix(1, 1);

  rElasticStiffnessMatrix(2, 2) = (EA * dz * dz) / L3;
  rElasticStiffnessMatrix(5, 5) = rElasticStiffnessMatrix(2, 2);

  rElasticStiffnessMatrix(0, 1) = (EA * dx * dy) / L3;
  rElasticStiffnessMatrix(1, 0) = rElasticStiffnessMatrix(0, 1);

  rElasticStiffnessMatrix(0, 2) = (EA * dx * dz) / L3;
  rElasticStiffnessMatrix(2, 0) = rElasticStiffnessMatrix(0, 2);

  rElasticStiffnessMatrix(0, 3) = -rElasticStiffnessMatrix(0, 0);
  rElasticStiffnessMatrix(3, 0) = rElasticStiffnessMatrix(0, 3);

  rElasticStiffnessMatrix(0, 4) = -rElasticStiffnessMatrix(0, 1);
  rElasticStiffnessMatrix(4, 0) = rElasticStiffnessMatrix(0, 4);

  rElasticStiffnessMatrix(0, 5) = -rElasticStiffnessMatrix(0, 2);
  rElasticStiffnessMatrix(5, 0) = rElasticStiffnessMatrix(0, 5);

  rElasticStiffnessMatrix(1, 2) = (EA * dy * dz) / L3;
  rElasticStiffnessMatrix(2, 1) = rElasticStiffnessMatrix(1, 2);

  rElasticStiffnessMatrix(1, 3) = rElasticStiffnessMatrix(0, 4);
  rElasticStiffnessMatrix(3, 1) = rElasticStiffnessMatrix(1, 3);

  rElasticStiffnessMatrix(1, 4) = -rElasticStiffnessMatrix(1, 1);
  rElasticStiffnessMatrix(4, 1) = rElasticStiffnessMatrix(1, 4);

  rElasticStiffnessMatrix(1, 5) = -rElasticStiffnessMatrix(1, 2);
  rElasticStiffnessMatrix(5, 1) = rElasticStiffnessMatrix(1, 5);

  rElasticStiffnessMatrix(2, 3) = -rElasticStiffnessMatrix(0, 2);
  rElasticStiffnessMatrix(3, 2) = rElasticStiffnessMatrix(2, 3);

  rElasticStiffnessMatrix(2, 4) = -rElasticStiffnessMatrix(1, 2);
  rElasticStiffnessMatrix(4, 2) = rElasticStiffnessMatrix(2, 4);

  rElasticStiffnessMatrix(2, 5) = -rElasticStiffnessMatrix(2, 2);
  rElasticStiffnessMatrix(5, 2) = rElasticStiffnessMatrix(2, 5);

  rElasticStiffnessMatrix(3, 4) = rElasticStiffnessMatrix(0, 1);
  rElasticStiffnessMatrix(4, 3) = rElasticStiffnessMatrix(3, 4);

  rElasticStiffnessMatrix(3, 5) = rElasticStiffnessMatrix(0, 2);
  rElasticStiffnessMatrix(5, 3) = rElasticStiffnessMatrix(3, 5);

  rElasticStiffnessMatrix(4, 5) = rElasticStiffnessMatrix(1, 2);
  rElasticStiffnessMatrix(5, 4) = rElasticStiffnessMatrix(4, 5);
  KRATOS_CATCH("")
}

void TrussElement3D2N::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
  rSerializer.save("mConstitutiveLaw", mConstitutiveLaw);
}
void TrussElement3D2N::load(Serializer &rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
  rSerializer.load("mConstitutiveLaw", mConstitutiveLaw);
}
} // namespace Kratos.
