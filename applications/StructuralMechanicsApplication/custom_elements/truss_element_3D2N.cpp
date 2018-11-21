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
#include "includes/checks.h"

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
        this->mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
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

  noalias(local_stiffness_matrix) += K_geo;

  return local_stiffness_matrix;
  KRATOS_CATCH("")
}

void TrussElement3D2N::CalculateDampingMatrix(
    MatrixType &rDampingMatrix, ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
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

  rDampingMatrix = alpha * mass_matrix;
  noalias(rDampingMatrix) += beta * stiffness_matrix;

  KRATOS_CATCH("")
}

void TrussElement3D2N::CalculateMassMatrix(MatrixType &rMassMatrix,
                                           ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  rMassMatrix = ZeroMatrix(msLocalSize, msLocalSize);

  const double A = this->GetProperties()[CROSS_AREA];
  const double L = this->CalculateReferenceLength();
  const double rho = this->GetProperties()[DENSITY];

  const double total_mass = A * L * rho;

  for (int i = 0; i < msNumberOfNodes; ++i) {

    for (int j = 0; j < msDimension; ++j) {
      int index = i * msDimension + j;

      rMassMatrix(index, index) = total_mass * 0.50;
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

  // Assemble global Vector
  for (int i = 0; i < msNumberOfNodes; ++i) {
    if (GetProperties().Has( VOLUME_ACCELERATION ))
        noalias(body_forces_node) = total_mass * GetProperties()[VOLUME_ACCELERATION];
    else if (this->Has( VOLUME_ACCELERATION ))
        noalias(body_forces_node) = total_mass * this->GetValue(VOLUME_ACCELERATION);
    else if( GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION) ) {
        noalias(body_forces_node) =
            total_mass *
            this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION) *
            Ncontainer(0, i);
    }

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
  // creating LHS
  rLeftHandSideMatrix =
      this->CreateElementStiffnessMatrix(rCurrentProcessInfo);

  // create+compute RHS
  rRightHandSideVector = -internal_forces;
  // add bodyforces
  if (this->HasSelfWeight()) noalias(rRightHandSideVector) += this->CalculateBodyForces();
  KRATOS_CATCH("")
}

void TrussElement3D2N::CalculateRightHandSide(
    VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY
  rRightHandSideVector = ZeroVector(msLocalSize);

  BoundedVector<double,msLocalSize> internal_forces =
    this->GetConstitutiveLawTrialResponse(rCurrentProcessInfo,false);

  BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
      ZeroMatrix(msLocalSize, msLocalSize);
  this->CreateTransformationMatrix(transformation_matrix);

  noalias(rRightHandSideVector) -= prod(transformation_matrix, internal_forces);

  // add bodyforces
  if (this->HasSelfWeight()) noalias(rRightHandSideVector) += this->CalculateBodyForces();
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
    Vector strain = ZeroVector(msDimension);
    strain[0] = this->CalculateGreenLagrangeStrain();
    strain[1] = 0.00;
    strain[2] = 0.00;
    rOutput[0] = strain;
  }
  if (rVariable == PK2_STRESS_VECTOR) {

    array_1d<double, 3 > truss_stresses;
    array_1d<double, msDimension> temp_internal_stresses = ZeroVector(msDimension);
    ProcessInfo temp_process_information;

    ConstitutiveLaw::Parameters Values(this->GetGeometry(),this->GetProperties(),temp_process_information);
    Vector temp_strain = ZeroVector(1);
    temp_strain[0] = this->CalculateGreenLagrangeStrain();
    Values.SetStrainVector(temp_strain);
    this->mpConstitutiveLaw->CalculateValue(Values,FORCE,temp_internal_stresses);

    rOutput[0] = temp_internal_stresses;
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

    const double L0 = this->CalculateReferenceLength();
    const double l = this->CalculateCurrentLength();


    array_1d<double, msDimension> temp_internal_stresses = ZeroVector(msDimension);
    ProcessInfo temp_process_information;
    ConstitutiveLaw::Parameters Values(this->GetGeometry(),this->GetProperties(),temp_process_information);

    Vector temp_strain = ZeroVector(1);
    temp_strain[0] = this->CalculateGreenLagrangeStrain();
    Values.SetStrainVector(temp_strain);
    this->mpConstitutiveLaw->CalculateValue(Values,FORCE,temp_internal_stresses);

    truss_forces[0] =
        ((temp_internal_stresses[0] + prestress) * l * A) / L0;

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
  const SizeType number_of_nodes = this->GetGeometry().size();
  const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

  if(this->mpConstitutiveLaw != nullptr) {
    this->mpConstitutiveLaw->Check(this->GetProperties(),this->GetGeometry(),rCurrentProcessInfo);
  }

  if (dimension != msDimension ||number_of_nodes != msNumberOfNodes) {
    KRATOS_ERROR << "The truss element works only in 3D and with 2 noded elements" << std::endl;
      }
  // verify that the variables are correctly initialized
  KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
  KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
  KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
  KRATOS_CHECK_VARIABLE_KEY(DENSITY);
  KRATOS_CHECK_VARIABLE_KEY(CROSS_AREA);

  // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
  for (IndexType i = 0; i < number_of_nodes; ++i) {
    NodeType &rnode = this->GetGeometry()[i];
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, rnode);

    KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode);
    KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode);
    KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode);
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

  const double l = this->CalculateCurrentLength();
  const double L0 = this->CalculateReferenceLength();
  const double A = this->GetProperties()[CROSS_AREA];

  double prestress = 0.00;
  if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
    prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
  }

  Vector temp_internal_stresses = ZeroVector(msLocalSize);
  ProcessInfo temp_process_information;
  ConstitutiveLaw::Parameters Values(this->GetGeometry(),this->GetProperties(),temp_process_information);


  Vector temp_strain = ZeroVector(1);
  temp_strain[0] = this->CalculateGreenLagrangeStrain();
  Values.SetStrainVector(temp_strain);
  this->mpConstitutiveLaw->CalculateValue(Values,NORMAL_STRESS,temp_internal_stresses);


    const double normal_force =
        ((temp_internal_stresses[3] + prestress) * l * A) / L0;

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
  const double numeric_limit = std::numeric_limits<double>::epsilon();
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
  if (VectorNorm > numeric_limit)
    direction_vector_x /= VectorNorm;

  else KRATOS_ERROR << "length of element" << this->Id() << "~ zero" << std::endl;

  if (std::abs(direction_vector_x[2]-1.00) <= numeric_limit) {
    direction_vector_y[1] = 1.0;
    direction_vector_z[0] = -1.0;
  }

  else if (std::abs(direction_vector_x[2]+1.00) <= numeric_limit) {
    direction_vector_y[1] = 1.0;
    direction_vector_z[0] = 1.0;
  }

  else {
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

    if (rRHSVariable == RESIDUAL_VECTOR &&
        rDestinationVariable == FORCE_RESIDUAL) {

    BoundedVector<double, msLocalSize> damping_residual_contribution =
        ZeroVector(msLocalSize);
      Vector current_nodal_velocities = ZeroVector(msLocalSize);
      this->GetFirstDerivativesVector(current_nodal_velocities);
      Matrix damping_matrix;
      ProcessInfo temp_process_information; // cant pass const ProcessInfo
      this->CalculateDampingMatrix(damping_matrix, temp_process_information);
      // current residual contribution due to damping
      noalias(damping_residual_contribution) =
          prod(damping_matrix, current_nodal_velocities);


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

    else if (rDestinationVariable == NODAL_INERTIA) {

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
  this->mpConstitutiveLaw->CalculateValue(Values,TANGENT_MODULUS,E);

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
  this->mpConstitutiveLaw->CalculateValue(Values,TANGENT_MODULUS,E);

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


void TrussElement3D2N::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY;
  this->GetConstitutiveLawTrialResponse(rCurrentProcessInfo,true);
  KRATOS_CATCH("");
}

void TrussElement3D2N::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY;
  Vector temp_shape_function = ZeroVector(3);
  this->mpConstitutiveLaw->FinalizeNonLinearIteration(this->GetProperties(),
	this->GetGeometry(),temp_shape_function,rCurrentProcessInfo);
  KRATOS_CATCH("");
}


BoundedVector<double,TrussElement3D2N::msLocalSize>
  TrussElement3D2N::GetConstitutiveLawTrialResponse(
   ProcessInfo& rCurrentProcessInfo,const bool& rSaveInternalVariables)
{
    KRATOS_TRY;
    Vector strain_vector = ZeroVector(this->mpConstitutiveLaw->GetStrainSize());
    Vector stress_vector = ZeroVector(this->mpConstitutiveLaw->GetStrainSize());
    strain_vector[0] = this->CalculateGreenLagrangeStrain();

    Matrix temp_matrix;
    Vector temp_vector;

    this->mpConstitutiveLaw->CalculateMaterialResponse(strain_vector,
    temp_matrix,stress_vector,temp_matrix,rCurrentProcessInfo,this->GetProperties(),
    this->GetGeometry(),temp_vector,true,true,rSaveInternalVariables);

    BoundedVector<double,msLocalSize> internal_forces = ZeroVector(msLocalSize);
    const double l = this->CalculateCurrentLength();
    const double L0 = this->CalculateReferenceLength();
    const double A = this->GetProperties()[CROSS_AREA];
    double prestress = 0.00;
    if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
      prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
    }

    const double normal_force =
        ((stress_vector[0] + prestress) * l * A) / L0;

    internal_forces[0] = -1.0 * normal_force;
    internal_forces[3] = +1.0 * normal_force;

    return internal_forces;
    KRATOS_CATCH("");
}


void TrussElement3D2N::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY;
  Vector temp_shape_function = ZeroVector(3);
  this->mpConstitutiveLaw->FinalizeSolutionStep(this->GetProperties(),
	this->GetGeometry(),temp_shape_function,rCurrentProcessInfo);
  KRATOS_CATCH("");
}


void TrussElement3D2N::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
  rSerializer.save("mpConstitutiveLaw", this->mpConstitutiveLaw);
}
void TrussElement3D2N::load(Serializer &rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
  rSerializer.load("mpConstitutiveLaw", this->mpConstitutiveLaw);
}

bool TrussElement3D2N::HasSelfWeight() const
{
  const double norm_self_weight =
   this->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION)[0]*
   this->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION)[0]+
   this->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION)[1]*
   this->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION)[1]+
   this->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION)[2]*
   this->GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION)[2];

  if (norm_self_weight<=std::numeric_limits<double>::epsilon()) return false;
  else return true;
}





} // namespace Kratos.
