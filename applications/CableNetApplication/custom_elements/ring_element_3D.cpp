//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
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
#include "custom_elements/ring_element_3D.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "includes/checks.h"


namespace Kratos {
RingElement3D::RingElement3D(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

RingElement3D::RingElement3D(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
RingElement3D::Create(IndexType NewId, NodesArrayType const &rThisNodes,
                         PropertiesType::Pointer pProperties) const {
  const GeometryType &rGeom = this->GetGeometry();
  return Kratos::make_intrusive<RingElement3D>(NewId, rGeom.Create(rThisNodes),
                                               pProperties);
}

Element::Pointer
RingElement3D::Create(IndexType NewId, GeometryType::Pointer pGeom,
                         PropertiesType::Pointer pProperties) const {
  return Kratos::make_intrusive<RingElement3D>(NewId, pGeom,
                                               pProperties);
}

RingElement3D::~RingElement3D() {}

void RingElement3D::EquationIdVector(EquationIdVectorType &rResult,
                                     const ProcessInfo &rCurrentProcessInfo) const {

  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;


  if (rResult.size() != local_size)
    rResult.resize(local_size);

  for (int i = 0; i < points_number; ++i) {
    int index = i * 3;
    rResult[index] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
    rResult[index + 1] =
        this->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[index + 2] =
        this->GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
  }
}
void RingElement3D::GetDofList(DofsVectorType &rElementalDofList,
                               const ProcessInfo &rCurrentProcessInfo) const {

  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;


  if (rElementalDofList.size() != local_size) {
    rElementalDofList.resize(local_size);
  }

  for (int i = 0; i < points_number; ++i) {
    int index = i * 3;
    rElementalDofList[index] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
    rElementalDofList[index + 1] =
        this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
    rElementalDofList[index + 2] =
        this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
  }
}


void RingElement3D::GetValuesVector(Vector &rValues, int Step) const {

  KRATOS_TRY;
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  if (rValues.size() != local_size)
    rValues.resize(local_size, false);

  for (int i = 0; i < points_number; ++i) {
    int index = i * dimension;
    const auto &disp =
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);

    rValues[index] = disp[0];
    rValues[index + 1] = disp[1];
    rValues[index + 2] = disp[2];
  }
  KRATOS_CATCH("")
}

void RingElement3D::GetFirstDerivativesVector(Vector &rValues, int Step) const {

  KRATOS_TRY;
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  if (rValues.size() != local_size)
    rValues.resize(local_size, false);

  for (int i = 0; i < points_number; ++i) {
    int index = i * dimension;
    const auto &vel =
        this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);

    rValues[index] = vel[0];
    rValues[index + 1] = vel[1];
    rValues[index + 2] = vel[2];
  }
  KRATOS_CATCH("")
}

void RingElement3D::GetSecondDerivativesVector(Vector &rValues, int Step) const {

  KRATOS_TRY;
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  if (rValues.size() != local_size)
    rValues.resize(local_size, false);

  for (int i = 0; i < points_number; ++i) {
    int index = i * dimension;
    const auto &acc =
        this->GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);

    rValues[index] = acc[0];
    rValues[index + 1] = acc[1];
    rValues[index + 2] = acc[2];
  }

  KRATOS_CATCH("")
}

Vector RingElement3D::GetCurrentLengthArray() const
{
  const int points_number = GetGeometry().PointsNumber();
  const int number_of_segments = points_number;

  Vector segment_lengths = ZeroVector(number_of_segments);
  for (int i=0;i<number_of_segments;++i)
  {
    int next_node_id = i+1;
    if (i==points_number-1) next_node_id = 0;

    const double du =
        this->GetGeometry()[next_node_id].FastGetSolutionStepValue(DISPLACEMENT_X) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X);
    const double dv =
        this->GetGeometry()[next_node_id].FastGetSolutionStepValue(DISPLACEMENT_Y) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
    const double dw =
        this->GetGeometry()[next_node_id].FastGetSolutionStepValue(DISPLACEMENT_Z) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Z);
    const double dx = this->GetGeometry()[next_node_id].X0() - this->GetGeometry()[i].X0();
    const double dy = this->GetGeometry()[next_node_id].Y0() - this->GetGeometry()[i].Y0();
    const double dz = this->GetGeometry()[next_node_id].Z0() - this->GetGeometry()[i].Z0();
    segment_lengths[i] = std::sqrt((du + dx) * (du + dx) + (dv + dy) * (dv + dy) +
                             (dw + dz) * (dw + dz));
  }
  return segment_lengths;
}

Vector RingElement3D::GetRefLengthArray() const
{
  const int points_number = GetGeometry().PointsNumber();
  const int number_of_segments = points_number;

  Vector segment_lengths = ZeroVector(number_of_segments);
  for (int i=0;i<number_of_segments;++i)
  {
    int next_node_id = i+1;
    if (i==points_number-1) next_node_id = 0;

    const double dx = this->GetGeometry()[next_node_id].X0() - this->GetGeometry()[i].X0();
    const double dy = this->GetGeometry()[next_node_id].Y0() - this->GetGeometry()[i].Y0();
    const double dz = this->GetGeometry()[next_node_id].Z0() - this->GetGeometry()[i].Z0();
    segment_lengths[i] = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
  }
  return segment_lengths;
}

double RingElement3D::GetCurrentLength() const
{
  const int points_number = GetGeometry().PointsNumber();
  const int number_of_segments = points_number;
  Vector segment_lengths = this->GetCurrentLengthArray();
  double length = 0.0;
  for (int i=0;i<number_of_segments;++i) length += segment_lengths[i];
  return length;
}

double RingElement3D::GetRefLength() const
{
  const int points_number = GetGeometry().PointsNumber();
  const int number_of_segments = points_number;
  Vector segment_lengths = this->GetRefLengthArray();
  double length = 0.0;
  for (int i=0;i<number_of_segments;++i) length += segment_lengths[i];
  return length;
}

Vector RingElement3D::GetDeltaPositions(const int& rDirection) const
{
  const int points_number = GetGeometry().PointsNumber();
  const int number_of_segments = points_number;

  Vector delta_position = ZeroVector(number_of_segments);

  double d_disp = 0.0;
  double d_ref_pos = 0.0;

  for (int i=0;i<number_of_segments;++i)
  {

    int next_node_id = i+1;
    if (i==points_number-1) next_node_id = 0;

    if (rDirection==1)
    {
      d_disp =
        this->GetGeometry()[next_node_id].FastGetSolutionStepValue(DISPLACEMENT_X) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X);
      d_ref_pos = this->GetGeometry()[next_node_id].X0() - this->GetGeometry()[i].X0();
    }

    else if (rDirection==2)
    {
      d_disp=
        this->GetGeometry()[next_node_id].FastGetSolutionStepValue(DISPLACEMENT_Y) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
      d_ref_pos = this->GetGeometry()[next_node_id].Y0() - this->GetGeometry()[i].Y0();
    }

    else if (rDirection==3)
    {
      d_disp =
        this->GetGeometry()[next_node_id].FastGetSolutionStepValue(DISPLACEMENT_Z) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Z);
      d_ref_pos = this->GetGeometry()[next_node_id].Z0() - this->GetGeometry()[i].Z0();
    }

    else KRATOS_ERROR << "maximum 3 dimensions" << std::endl;

    delta_position[i] = d_disp+d_ref_pos;
  }
  return delta_position;
}

double RingElement3D::CalculateGreenLagrangeStrain() const
{
  const double reference_length = this->GetRefLength();
  const double current_length = this->GetCurrentLength();
  return (0.50 * (((current_length*current_length)-(reference_length*reference_length)) / (reference_length*reference_length)));
}

Vector RingElement3D::GetDirectionVectorNt() const
{
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  Vector n_t = ZeroVector(local_size);
  const Vector d_x = this->GetDeltaPositions(1);
  const Vector d_y = this->GetDeltaPositions(2);
  const Vector d_z = this->GetDeltaPositions(3);
  const Vector lengths = this->GetCurrentLengthArray();

  const int last_node = points_number-1;

  n_t[0] = (d_x[last_node]/lengths[last_node])-(d_x[0]/lengths[0]);
  n_t[1] = (d_y[last_node]/lengths[last_node])-(d_y[0]/lengths[0]);
  n_t[2] = (d_z[last_node]/lengths[last_node])-(d_z[0]/lengths[0]);

  n_t[3] = (d_x[0]/lengths[0])-(d_x[1]/lengths[1]);
  n_t[4] = (d_y[0]/lengths[0])-(d_y[1]/lengths[1]);
  n_t[5] = (d_z[0]/lengths[0])-(d_z[1]/lengths[1]);

  n_t[6] = (d_x[1]/lengths[1])-(d_x[2]/lengths[2]);
  n_t[7] = (d_y[1]/lengths[1])-(d_y[2]/lengths[2]);
  n_t[8] = (d_z[1]/lengths[1])-(d_z[2]/lengths[2]);

  if (points_number==4)
  {
    n_t[9]  = (d_x[2]/lengths[2])-(d_x[3]/lengths[3]);
    n_t[10] = (d_y[2]/lengths[2])-(d_y[3]/lengths[3]);
    n_t[11] = (d_z[2]/lengths[2])-(d_z[3]/lengths[3]);
  }

  return n_t;
}

Vector RingElement3D::GetInternalForces() const
{
  const double k_0            = this->LinearStiffness();
  const double strain_gl      = this->CalculateGreenLagrangeStrain();
  const double current_length = this->GetCurrentLength();

  const double total_internal_force = k_0 * strain_gl * current_length;
  const Vector internal_foces = this->GetDirectionVectorNt() * total_internal_force;
  return internal_foces;
}

Matrix RingElement3D::ElasticStiffnessMatrix() const
{
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  Matrix elastic_stiffness_matrix = ZeroMatrix(local_size,local_size);
  const Vector direction_vector = this->GetDirectionVectorNt();
  elastic_stiffness_matrix = outer_prod(direction_vector,direction_vector);
  elastic_stiffness_matrix *= this->LinearStiffness();
  elastic_stiffness_matrix *= (1.0 + (3.0*this->CalculateGreenLagrangeStrain()));
  return elastic_stiffness_matrix;
}

Matrix RingElement3D::GeometricStiffnessMatrix() const
{
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  const double k_0            = this->LinearStiffness();
  const double strain_gl      = this->CalculateGreenLagrangeStrain();
  const double current_length = this->GetCurrentLength();

  const Vector d_x = this->GetDeltaPositions(1);
  const Vector d_y = this->GetDeltaPositions(2);
  const Vector d_z = this->GetDeltaPositions(3);
  const Vector lengths = this->GetCurrentLengthArray();

  Matrix geometric_stiffness_matrix = ZeroMatrix(local_size,local_size);

  for (int i=0;i<points_number;++i)
  {
    const double current_length = lengths[i];
    Matrix sub_stiffness_matrix = ZeroMatrix(dimension,dimension);
    sub_stiffness_matrix(0, 0) = -std::pow(d_x[i],2.0)+std::pow(lengths[i],2.0);
    sub_stiffness_matrix(0, 1) = -d_x[i]*d_y[i];
    sub_stiffness_matrix(0, 2) = -d_x[i]*d_z[i];

    sub_stiffness_matrix(1, 0) = sub_stiffness_matrix(0, 1);
    sub_stiffness_matrix(1, 1) = -std::pow(d_y[i],2.0)+std::pow(lengths[i],2.0);
    sub_stiffness_matrix(1, 2) = -d_y[i]*d_z[i];

    sub_stiffness_matrix(2, 0) = sub_stiffness_matrix(0, 2);
    sub_stiffness_matrix(2, 1) = sub_stiffness_matrix(1, 2);
    sub_stiffness_matrix(2, 2) = -std::pow(d_z[i],2.0)+std::pow(lengths[i],2.0);

    sub_stiffness_matrix /= std::pow(current_length,3.0);


    if (i==points_number-1)
    {
      project(geometric_stiffness_matrix, range(0,3),range(0,3)) += sub_stiffness_matrix;
      project(geometric_stiffness_matrix, range(local_size-dimension,local_size),range(local_size-dimension,local_size)) += sub_stiffness_matrix;

      project(geometric_stiffness_matrix, range(0,3),range(local_size-dimension,local_size)) -= sub_stiffness_matrix;
      project(geometric_stiffness_matrix, range(local_size-dimension,local_size),range(0,3)) -= sub_stiffness_matrix;
    }
    else
    {
      project(geometric_stiffness_matrix, range((i*3),((i+1)*3)),range((i*3),((i+1)*3))) += sub_stiffness_matrix;
      project(geometric_stiffness_matrix, range((i*3)+3,((i+1)*3)+3),range((i*3)+3,((i+1)*3)+3)) += sub_stiffness_matrix;

      project(geometric_stiffness_matrix, range((i*3),((i+1)*3)),range((i*3)+3,((i+1)*3)+3)) -= sub_stiffness_matrix;
      project(geometric_stiffness_matrix, range((i*3)+3,((i+1)*3)+3),range((i*3),((i+1)*3))) -= sub_stiffness_matrix;
    }
  }

  geometric_stiffness_matrix *= k_0 * strain_gl * current_length;

  return geometric_stiffness_matrix;
}

inline Matrix RingElement3D::TotalStiffnessMatrix() const
{
  const Matrix ElasticStiffnessMatrix = this->ElasticStiffnessMatrix();
  const Matrix GeometrixStiffnessMatrix = this->GeometricStiffnessMatrix();
  return (ElasticStiffnessMatrix+GeometrixStiffnessMatrix);
}

void RingElement3D::CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY;
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;
  // resizing the matrices + create memory for LHS
  rLeftHandSideMatrix = ZeroMatrix(local_size, local_size);
  // creating LHS
  noalias(rLeftHandSideMatrix) = this->TotalStiffnessMatrix();

  KRATOS_CATCH("")
}

void RingElement3D::CalculateRightHandSide(
    VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
{
  KRATOS_TRY;
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  rRightHandSideVector = ZeroVector(local_size);
  noalias(rRightHandSideVector) -= this->GetInternalForces();
  if (this->HasSelfWeight()) noalias(rRightHandSideVector) += this->CalculateBodyForces();
  KRATOS_CATCH("")
}

void RingElement3D::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                         VectorType &rRightHandSideVector,
                                         const ProcessInfo &rCurrentProcessInfo)
{
  KRATOS_TRY;
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;


  rLeftHandSideMatrix = ZeroMatrix(local_size, local_size);
  noalias(rLeftHandSideMatrix) = this->TotalStiffnessMatrix();

  rRightHandSideVector = ZeroVector(local_size);
  noalias(rRightHandSideVector) -= this->GetInternalForces();
  if (this->HasSelfWeight()) noalias(rRightHandSideVector) += this->CalculateBodyForces();
  KRATOS_CATCH("")
}

void RingElement3D::CalculateLumpedMassVector(VectorType &rMassVector)
{
    KRATOS_TRY;


    // ATTENTION !!!!
    // this function uses a fictiuous mass for the sliding nodes
    // needs improvement !!!

    const int points_number = GetGeometry().PointsNumber();
    const int dimension = 3;
    const SizeType local_size = dimension*points_number;

    // Clear matrix
    if (rMassVector.size() != local_size)
        rMassVector.resize( local_size );

    const double A = this->GetProperties()[CROSS_AREA];
    const double L = this->GetRefLength();
    const double rho = this->GetProperties()[DENSITY];

    const double total_mass = A * L * rho;

    for (int i = 0; i < points_number; ++i) {
        for (int j = 0; j < dimension; ++j) {
            int index = i * dimension + j;

            rMassVector[index] = total_mass;
        }
    }

    KRATOS_CATCH("")
}

void RingElement3D::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;
    const int points_number = GetGeometry().PointsNumber();
    const int dimension = 3;
    const SizeType local_size = dimension*points_number;

    // Compute lumped mass matrix
    VectorType temp_vector(local_size);
    CalculateLumpedMassVector(temp_vector);

    // Clear matrix
    if (rMassMatrix.size1() != local_size || rMassMatrix.size2() != local_size)
        rMassMatrix.resize( local_size, local_size, false );
    rMassMatrix = ZeroMatrix(local_size, local_size);

    // Fill the matrix
    for (IndexType i = 0; i < local_size; ++i)
        rMassMatrix(i, i) = temp_vector[i];

    KRATOS_CATCH("")
}

void RingElement3D::CalculateDampingMatrix(
    MatrixType &rDampingMatrix, const ProcessInfo &rCurrentProcessInfo) {

  KRATOS_TRY;
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  MatrixType stiffness_matrix = ZeroMatrix(local_size, local_size);

  this->CalculateLeftHandSide(stiffness_matrix, rCurrentProcessInfo);

  MatrixType mass_matrix = ZeroMatrix(local_size, local_size);

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

void RingElement3D::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    const int points_number = GetGeometry().PointsNumber();
    const int dimension = 3;
    const SizeType local_size = dimension*points_number;

    auto& r_geom = GetGeometry();

    if (rDestinationVariable == NODAL_MASS) {
        VectorType element_mass_vector(local_size);
        this->CalculateLumpedMassVector(element_mass_vector);

        for (int i = 0; i < points_number; ++i) {
            double &r_nodal_mass = r_geom[i].GetValue(NODAL_MASS);
            int index = i * dimension;

            #pragma omp atomic
            r_nodal_mass += element_mass_vector(index);
        }
    }
    KRATOS_CATCH("")
}

void RingElement3D::AddExplicitContribution(
    const VectorType &rRHSVector, const Variable<VectorType> &rRHSVariable,
    const Variable<array_1d<double, 3>> &rDestinationVariable,
    const ProcessInfo &rCurrentProcessInfo
    )
{
    KRATOS_TRY;
    const int points_number = GetGeometry().PointsNumber();
    const int dimension = 3;
    const SizeType local_size = dimension*points_number;

    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {

        Vector damping_residual_contribution = ZeroVector(local_size);
        Vector current_nodal_velocities = ZeroVector(local_size);
        this->GetFirstDerivativesVector(current_nodal_velocities);
        Matrix damping_matrix;
        ProcessInfo temp_process_information; // cant pass const ProcessInfo
        this->CalculateDampingMatrix(damping_matrix, temp_process_information);
        // current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);

        for (int i = 0; i < points_number; ++i) {
            size_t index = dimension * i;
            array_1d<double, 3> &r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for (size_t j = 0; j < dimension; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
            }
        }
    } else if (rDestinationVariable == NODAL_INERTIA) {

        // Getting the vector mass
        VectorType mass_vector(local_size);
        CalculateLumpedMassVector(mass_vector);

        for (int i = 0; i < points_number; ++i) {
            double &r_nodal_mass = GetGeometry()[i].GetValue(NODAL_MASS);
            array_1d<double, dimension> &r_nodal_inertia = GetGeometry()[i].GetValue(NODAL_INERTIA);
            int index = i * dimension;

            #pragma omp atomic
            r_nodal_mass += mass_vector[index];

            for (int k = 0; k < dimension; ++k) {
                #pragma omp atomic
                r_nodal_inertia[k] += 0.0;
            }
        }
    }
    KRATOS_CATCH("")
}

int RingElement3D::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF( this->Id() < 1 ) << "Element found with Id " << this->Id() << std::endl;

    const double domain_size = this->GetCurrentLength();
    KRATOS_ERROR_IF( domain_size <= 0.0 ) << "Element " << this->Id() << " has non-positive size " << domain_size << std::endl;


    KRATOS_ERROR_IF_NOT((GetGeometry().PointsNumber()==4) || (GetGeometry().PointsNumber()==3)) << "the ring element " << this->Id() << "does not have 3 || 4 nodes" << std::endl;

    return 0;

    KRATOS_CATCH("")
}

bool RingElement3D::HasSelfWeight() const
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

Vector RingElement3D::CalculateBodyForces() {

    const int points_number = GetGeometry().PointsNumber();
    const int dimension = 3;
    const SizeType local_size = dimension*points_number;

    // creating necessary values
    const double A = this->GetProperties()[CROSS_AREA];
    const Vector l_array = this->GetCurrentLengthArray();
    const double l = this->GetCurrentLength();
    const double L = this->GetRefLength();
    const double rho = this->GetProperties()[DENSITY];

    double total_mass = A * L * rho;
    Vector body_forces_node = ZeroVector(dimension);
    Vector body_forces_global = ZeroVector(local_size);



    // assemble global Vector
    for (int i = 0; i < points_number; ++i) {

      const SizeType node_i = i;
      SizeType node_j = i-1;
      if (i==0) node_j = points_number-1;
      const double weight_fraction = (l_array[node_i]+l_array[node_j])/l;

      body_forces_node = total_mass *
        this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION) *
        weight_fraction * 0.5;

      for (unsigned int j = 0; j < dimension; ++j) {
        body_forces_global[(i * dimension) + j] = body_forces_node[j];
      }
    }

    return body_forces_global;
}

double RingElement3D::LinearStiffness() const
{
    return (this->GetProperties()[CROSS_AREA] * this->GetProperties()[YOUNG_MODULUS] / this->GetRefLength());
}

void RingElement3D::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}
void RingElement3D::load(Serializer &rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}
} // namespace Kratos.