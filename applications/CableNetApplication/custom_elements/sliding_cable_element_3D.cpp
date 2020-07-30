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
#include "custom_elements/sliding_cable_element_3D.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "includes/checks.h"


namespace Kratos {
SlidingCableElement3D::SlidingCableElement3D(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

SlidingCableElement3D::SlidingCableElement3D(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
SlidingCableElement3D::Create(IndexType NewId, NodesArrayType const &rThisNodes,
                         PropertiesType::Pointer pProperties) const {
  const GeometryType &rGeom = this->GetGeometry();
  return Kratos::make_intrusive<SlidingCableElement3D>(NewId, rGeom.Create(rThisNodes),
                                               pProperties);
}

Element::Pointer
SlidingCableElement3D::Create(IndexType NewId, GeometryType::Pointer pGeom,
                         PropertiesType::Pointer pProperties) const {
  return Kratos::make_intrusive<SlidingCableElement3D>(NewId, pGeom,
                                               pProperties);
}

SlidingCableElement3D::~SlidingCableElement3D() {}

void SlidingCableElement3D::EquationIdVector(EquationIdVectorType &rResult,
                                        ProcessInfo &rCurrentProcessInfo) {

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
void SlidingCableElement3D::GetDofList(DofsVectorType &rElementalDofList,
                                  ProcessInfo &rCurrentProcessInfo) {

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

void SlidingCableElement3D::Initialize()
{
    KRATOS_TRY
    mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
    KRATOS_CATCH("")
}

void SlidingCableElement3D::GetValuesVector(Vector &rValues, int Step) {

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

void SlidingCableElement3D::GetFirstDerivativesVector(Vector &rValues, int Step) {

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

void SlidingCableElement3D::GetSecondDerivativesVector(Vector &rValues, int Step) {

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

Vector SlidingCableElement3D::GetCurrentLengthArray(int step) const
{
  const int points_number = GetGeometry().PointsNumber();
  const int number_of_segments = points_number-1;

  Vector segment_lengths = ZeroVector(number_of_segments);
  for (int i=0;i<number_of_segments;++i)
  {
    const double du =
        this->GetGeometry()[i+1].FastGetSolutionStepValue(DISPLACEMENT_X,step) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X,step);
    const double dv =
        this->GetGeometry()[i+1].FastGetSolutionStepValue(DISPLACEMENT_Y,step) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y,step);
    const double dw =
        this->GetGeometry()[i+1].FastGetSolutionStepValue(DISPLACEMENT_Z,step) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Z,step);
    const double dx = this->GetGeometry()[i+1].X0() - this->GetGeometry()[i].X0();
    const double dy = this->GetGeometry()[i+1].Y0() - this->GetGeometry()[i].Y0();
    const double dz = this->GetGeometry()[i+1].Z0() - this->GetGeometry()[i].Z0();
    segment_lengths[i] = std::sqrt((du + dx) * (du + dx) + (dv + dy) * (dv + dy) +
                             (dw + dz) * (dw + dz));
  }
  return segment_lengths;
}

Vector SlidingCableElement3D::GetRefLengthArray() const
{
  const int points_number = GetGeometry().PointsNumber();
  const int number_of_segments = points_number-1;

  Vector segment_lengths = ZeroVector(number_of_segments);
  for (int i=0;i<number_of_segments;++i)
  {
    const double dx = this->GetGeometry()[i+1].X0() - this->GetGeometry()[i].X0();
    const double dy = this->GetGeometry()[i+1].Y0() - this->GetGeometry()[i].Y0();
    const double dz = this->GetGeometry()[i+1].Z0() - this->GetGeometry()[i].Z0();
    segment_lengths[i] = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
  }
  return segment_lengths;
}

double SlidingCableElement3D::GetCurrentLength() const
{
  const int points_number = GetGeometry().PointsNumber();
  const int number_of_segments = points_number-1;
  Vector segment_lengths = this->GetCurrentLengthArray();
  double length = 0.0;
  for (int i=0;i<number_of_segments;++i) length += segment_lengths[i];
  return length;
}

double SlidingCableElement3D::GetRefLength() const
{
  const int points_number = GetGeometry().PointsNumber();
  const int number_of_segments = points_number-1;
  Vector segment_lengths = this->GetRefLengthArray();
  double length = 0.0;
  for (int i=0;i<number_of_segments;++i) length += segment_lengths[i];
  return length;
}

Vector SlidingCableElement3D::GetDeltaPositions(const int& rDirection) const
{
  const int points_number = GetGeometry().PointsNumber();
  const int number_of_segments = points_number-1;

  Vector delta_position = ZeroVector(number_of_segments);

  double d_disp = 0.0;
  double d_ref_pos = 0.0;

  for (int i=0;i<number_of_segments;++i)
  {

    if (rDirection==1)
    {
      d_disp =
        this->GetGeometry()[i+1].FastGetSolutionStepValue(DISPLACEMENT_X) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X);
      d_ref_pos = this->GetGeometry()[i+1].X0() - this->GetGeometry()[i].X0();
    }

    else if (rDirection==2)
    {
      d_disp=
        this->GetGeometry()[i+1].FastGetSolutionStepValue(DISPLACEMENT_Y) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
      d_ref_pos = this->GetGeometry()[i+1].Y0() - this->GetGeometry()[i].Y0();
    }

    else if (rDirection==3)
    {
      d_disp =
        this->GetGeometry()[i+1].FastGetSolutionStepValue(DISPLACEMENT_Z) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Z);
      d_ref_pos = this->GetGeometry()[i+1].Z0() - this->GetGeometry()[i].Z0();
    }

    else KRATOS_ERROR << "maximum 3 dimensions" << std::endl;

    delta_position[i] = d_disp+d_ref_pos;
  }
  return delta_position;
}

double SlidingCableElement3D::CalculateGreenLagrangeStrain() const
{
  const double reference_length = this->GetRefLength();
  const double current_length = this->GetCurrentLength();
  return (0.50 * (((current_length*current_length)-(reference_length*reference_length)) / (reference_length*reference_length)));
}

Vector SlidingCableElement3D::GetDirectionVectorNt() const
{
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  Vector n_t = ZeroVector(local_size);
  const Vector d_x = this->GetDeltaPositions(1);
  const Vector d_y = this->GetDeltaPositions(2);
  const Vector d_z = this->GetDeltaPositions(3);
  const Vector lengths = this->GetCurrentLengthArray();

  n_t[0] = -d_x[0] / lengths[0];
  n_t[1] = -d_y[0] / lengths[0];
  n_t[2] = -d_z[0] / lengths[0];

  for (int i=0;i<points_number-2;++i)
  {
    n_t[(i+1)*3]     = (d_x[i]/lengths[i])-(d_x[i+1]/lengths[i+1]);
    n_t[((i+1)*3)+1] = (d_y[i]/lengths[i])-(d_y[i+1]/lengths[i+1]);
    n_t[((i+1)*3)+2] = (d_z[i]/lengths[i])-(d_z[i+1]/lengths[i+1]);
  }

  n_t[local_size-dimension]   = d_x[points_number-2]/lengths[points_number-2];
  n_t[local_size-dimension+1] = d_y[points_number-2]/lengths[points_number-2];
  n_t[local_size-dimension+2] = d_z[points_number-2]/lengths[points_number-2];

  return n_t;
}


Vector SlidingCableElement3D::GetInternalForces()
{
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const int segments_number = points_number-1;
  const int points_int_number = points_number-2;

  const double current_length = this->GetCurrentLength();
  const double prestress     = this->GetPK2PrestressValue();
  const double area           = this->GetProperties()[CROSS_AREA];
  const double ref_length     = this->GetRefLength();


  Vector temp_internal_stresses = ZeroVector(6);
  ProcessInfo temp_process_information;
  ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),temp_process_information);

  Vector temp_strain = ZeroVector(1);
  temp_strain[0] = CalculateGreenLagrangeStrain();
  Values.SetStrainVector(temp_strain);
  mpConstitutiveLaw->CalculateValue(Values,NORMAL_STRESS,temp_internal_stresses);


  const double total_internal_force = (temp_internal_stresses[3]+prestress) * area * current_length / ref_length;




  Vector internal_forces = total_internal_force*this->GetDirectionVectorNt();

  double friction_coefficient = 0.0;
  if (this->GetProperties().Has(FRICTION_COEFFICIENT))  friction_coefficient = GetProperties()[FRICTION_COEFFICIENT];

  if (friction_coefficient>0.0)
  {
    Vector l_p = this->CalculateProjectionLengths();

    /*   Vector l_0 = this->GetCurrentLengthArray();
      Vector l_1 = this->GetCurrentLengthArray(2);
      Vector d_l = l_0-l_1; */

    Vector l_0 = this->GetRefLengthArray();
    Vector d_l = l_p-l_0;

    Vector e_l = ZeroVector(segments_number);
    for (int i=0;i<segments_number;++i) e_l[i] = d_l[i]/l_0[i];

    //const double total_internal_force = k_0 * strain_gl * current_length;
    Vector deviation_forces = this->GetDirectionVectorNt() * total_internal_force;

    Vector internal_normal_resulting_forces = ZeroVector(segments_number);
    for (int i=0;i<segments_number;++i) internal_normal_resulting_forces[i] = total_internal_force;

    /* starting_segment_id = 0;
    for (int i=0;i<segments_number;++i)
    {
      if (d_l[i]>0.0)
      {
        starting_segment_id = i;
        break;
      }
    } */
    // simplified start at left segment and say N1 = N

    for (int i=0;i<points_int_number;++i)
    {
      int node_i = i+1;
      //double dl_i_0 = d_l[i];
      //double dl_i_1 = d_l[i+1];

      double el_i_0 = e_l[i];
      double el_i_1 = e_l[i+1];

      double deviation_force_i =
      std::sqrt((deviation_forces[node_i*dimension]*deviation_forces[node_i*dimension])+
      (deviation_forces[(node_i*dimension)+1]*deviation_forces[(node_i*dimension)+1])+
      (deviation_forces[(node_i*dimension)+2]*deviation_forces[(node_i*dimension)+2]));

      double friction_force = deviation_force_i*friction_coefficient;

      if ((el_i_0>=0) && (el_i_1>=0)) internal_normal_resulting_forces[i+1] = internal_normal_resulting_forces[i];
      else
      {
        double next_n = 0.0;
        if (el_i_1<el_i_0) next_n = internal_normal_resulting_forces[i] - friction_force;
        else next_n = internal_normal_resulting_forces[i] + friction_force;

        internal_normal_resulting_forces[i+1] = next_n;
      }
    }

    internal_forces = this->GetCustomInternalForceWithFriction(internal_normal_resulting_forces);
  }

  return internal_forces;
}

Vector SlidingCableElement3D::GetCustomInternalForceWithFriction(const Vector& rNormalForces)
{
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  Vector n_t = ZeroVector(local_size);
  const Vector d_x = this->GetDeltaPositions(1);
  const Vector d_y = this->GetDeltaPositions(2);
  const Vector d_z = this->GetDeltaPositions(3);
  const Vector lengths = this->GetCurrentLengthArray();

  n_t[0] = (-d_x[0] / lengths[0])*rNormalForces[0];
  n_t[1] = (-d_y[0] / lengths[0])*rNormalForces[0];
  n_t[2] = (-d_z[0] / lengths[0])*rNormalForces[0];

  for (int i=0;i<points_number-2;++i)
  {
    n_t[(i+1)*3]     = ((d_x[i]/lengths[i])*rNormalForces[i])-((d_x[i+1]/lengths[i+1])*rNormalForces[i+1]);
    n_t[((i+1)*3)+1] = ((d_y[i]/lengths[i])*rNormalForces[i])-((d_y[i+1]/lengths[i+1])*rNormalForces[i+1]);
    n_t[((i+1)*3)+2] = ((d_z[i]/lengths[i])*rNormalForces[i])-((d_z[i+1]/lengths[i+1])*rNormalForces[i+1]);
  }

  n_t[local_size-dimension]   = (d_x[points_number-2]/lengths[points_number-2])*rNormalForces[points_number-2];
  n_t[local_size-dimension+1] = (d_y[points_number-2]/lengths[points_number-2])*rNormalForces[points_number-2];
  n_t[local_size-dimension+2] = (d_z[points_number-2]/lengths[points_number-2])*rNormalForces[points_number-2];

  return n_t;
}


Matrix SlidingCableElement3D::ElasticStiffnessMatrix(const ProcessInfo& rCurrentProcessInfo) const
{
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  const double youngs_mod = ReturnTangentModulus1D(rCurrentProcessInfo);


  const double prestress     = this->GetPK2PrestressValue();

  Matrix elastic_stiffness_matrix = ZeroMatrix(local_size,local_size);
  const Vector direction_vector = this->GetDirectionVectorNt();
  elastic_stiffness_matrix = outer_prod(direction_vector,direction_vector);
  elastic_stiffness_matrix *= this->LinearStiffness();
  elastic_stiffness_matrix *= (1.0 + (3.0*(this->CalculateGreenLagrangeStrain()+(prestress/youngs_mod))));
  return elastic_stiffness_matrix;
}

Matrix SlidingCableElement3D::GeometricStiffnessMatrix(const ProcessInfo& rCurrentProcessInfo) const
{
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  //const double k_0            = this->LinearStiffness();
  const double strain_gl      = this->CalculateGreenLagrangeStrain();
  const double current_length = this->GetCurrentLength();

  const double prestress     = this->GetPK2PrestressValue();
  const double area           = this->GetProperties()[CROSS_AREA];


  const double youngs_mod = ReturnTangentModulus1D(rCurrentProcessInfo);


  const double ref_length     = this->GetRefLength();

  const Vector d_x = this->GetDeltaPositions(1);
  const Vector d_y = this->GetDeltaPositions(2);
  const Vector d_z = this->GetDeltaPositions(3);
  const Vector lengths = this->GetCurrentLengthArray();

  Matrix geometric_stiffness_matrix = ZeroMatrix(local_size,local_size);

  for (int i=0;i<points_number-1;++i)
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

    project(geometric_stiffness_matrix, range((i*3),((i+1)*3)),range((i*3),((i+1)*3))) += sub_stiffness_matrix;
    project(geometric_stiffness_matrix, range((i*3)+3,((i+1)*3)+3),range((i*3)+3,((i+1)*3)+3)) += sub_stiffness_matrix;

    project(geometric_stiffness_matrix, range((i*3),((i+1)*3)),range((i*3)+3,((i+1)*3)+3)) -= sub_stiffness_matrix;
    project(geometric_stiffness_matrix, range((i*3)+3,((i+1)*3)+3),range((i*3),((i+1)*3))) -= sub_stiffness_matrix;
  }

  geometric_stiffness_matrix *= (youngs_mod*strain_gl+prestress) * area * current_length / ref_length;

  return geometric_stiffness_matrix;
}

inline Matrix SlidingCableElement3D::TotalStiffnessMatrix(const ProcessInfo& rCurrentProcessInfo) const
{
  const Matrix ElasticStiffnessMatrix = this->ElasticStiffnessMatrix(rCurrentProcessInfo);
  const Matrix GeometrixStiffnessMatrix = this->GeometricStiffnessMatrix(rCurrentProcessInfo);
  return (ElasticStiffnessMatrix+GeometrixStiffnessMatrix);
}

void SlidingCableElement3D::CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY;
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;
  // resizing the matrices + create memory for LHS
  rLeftHandSideMatrix = ZeroMatrix(local_size, local_size);
  // creating LHS
  noalias(rLeftHandSideMatrix) = this->TotalStiffnessMatrix(rCurrentProcessInfo);

  KRATOS_CATCH("")
}

void SlidingCableElement3D::CalculateRightHandSide(
    VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
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

void SlidingCableElement3D::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                            VectorType &rRightHandSideVector,
                                            ProcessInfo &rCurrentProcessInfo)
{
  KRATOS_TRY;
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;

  rLeftHandSideMatrix = ZeroMatrix(local_size, local_size);
  noalias(rLeftHandSideMatrix) = this->TotalStiffnessMatrix(rCurrentProcessInfo);

  rRightHandSideVector = ZeroVector(local_size);
  noalias(rRightHandSideVector) -= this->GetInternalForces();

  if (this->HasSelfWeight()) noalias(rRightHandSideVector) += this->CalculateBodyForces();
  KRATOS_CATCH("")
}

void SlidingCableElement3D::CalculateLumpedMassVector(VectorType &rMassVector)
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
    const Vector l_array = this->GetCurrentLengthArray();
    const double l = this->GetCurrentLength();

    const double total_mass = A * L * rho;

    for (int i = 0; i < points_number; ++i) {

        double weight_fraction = 0.0;
        if (i==0) weight_fraction = l_array[i]/l;
        else if (i==points_number-1) weight_fraction = l_array[i-1]/l;
        else weight_fraction = (l_array[i]+l_array[i-1])/l;

        double nodal_mass = total_mass * weight_fraction * 0.5;

        for (int j = 0; j < dimension; ++j) {
            int index = i * dimension + j;

            rMassVector[index] = nodal_mass;
            //rMassVector[index] = total_mass;
        }
    }

    KRATOS_CATCH("")
}

void SlidingCableElement3D::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    ProcessInfo &rCurrentProcessInfo)
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

void SlidingCableElement3D::CalculateDampingMatrix(
    MatrixType &rDampingMatrix, ProcessInfo &rCurrentProcessInfo) {

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

void SlidingCableElement3D::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    Variable<double >& rDestinationVariable,
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

void SlidingCableElement3D::AddExplicitContribution(
    const VectorType &rRHSVector, const Variable<VectorType> &rRHSVariable,
    Variable<array_1d<double, 3>> &rDestinationVariable,
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

int SlidingCableElement3D::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF( this->Id() < 1 ) << "Element found with Id " << this->Id() << std::endl;

    const double domain_size = this->GetCurrentLength();
    KRATOS_ERROR_IF( domain_size <= 0.0 ) << "Element " << this->Id() << " has non-positive size " << domain_size << std::endl;


    KRATOS_ERROR_IF_NOT(GetProperties()[CONSTITUTIVE_LAW])  << "A constitutive law needs to be specified for the element with ID " << Id() << std::endl;
    mpConstitutiveLaw->Check(GetProperties(),GetGeometry(),rCurrentProcessInfo);

    return 0;

    KRATOS_CATCH("")
}

bool SlidingCableElement3D::HasSelfWeight() const
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



Vector SlidingCableElement3D::CalculateBodyForces() {

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
      double weight_fraction = 0.0;
      if (i==0) weight_fraction = l_array[i]/l;
      else if (i==points_number-1) weight_fraction = l_array[i-1]/l;
      else weight_fraction = (l_array[i]+l_array[i-1])/l;

      body_forces_node = total_mass *
        this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION) *
        weight_fraction * 0.5;

      for (unsigned int j = 0; j < dimension; ++j) {
        body_forces_global[(i * dimension) + j] = body_forces_node[j];
      }
    }

    return body_forces_global;
}


Vector SlidingCableElement3D::CalculateProjectionLengths()
{
  const int points_number = GetGeometry().PointsNumber();
  const int number_of_segments = points_number-1;

  Vector projection_lengths = ZeroVector(number_of_segments);

  const Vector d_x = this->GetDeltaPositions(1);
  const Vector d_y = this->GetDeltaPositions(2);
  const Vector d_z = this->GetDeltaPositions(3);

  const Vector l_0 = this->GetRefLengthArray();

  for (int i=0;i<number_of_segments;++i)
  {
    double l_p_i = 0.00;

    const double dX = this->GetGeometry()[i+1].X0() - this->GetGeometry()[i].X0();
    const double dY = this->GetGeometry()[i+1].Y0() - this->GetGeometry()[i].Y0();
    const double dZ = this->GetGeometry()[i+1].Z0() - this->GetGeometry()[i].Z0();

    const double l_0_i = l_0[i];

    l_p_i += dX*d_x[i];
    l_p_i += dY*d_y[i];
    l_p_i += dZ*d_z[i];

    l_p_i /= l_0_i;

    projection_lengths[i] = l_p_i;
  }

  return projection_lengths;
}


void SlidingCableElement3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    ConstitutiveLaw::Parameters element_parameters;
    mpConstitutiveLaw->FinalizeMaterialResponse(element_parameters,ConstitutiveLaw::StressMeasure_PK2);
    KRATOS_CATCH("");
}


void SlidingCableElement3D::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    GetConstitutiveLawTrialResponse(rCurrentProcessInfo);
    KRATOS_CATCH("");
}

void SlidingCableElement3D::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    Vector temp_shape_function = ZeroVector(3);
    mpConstitutiveLaw->FinalizeNonLinearIteration(GetProperties(),
            GetGeometry(),temp_shape_function,rCurrentProcessInfo);
    KRATOS_CATCH("");
}


void SlidingCableElement3D::GetConstitutiveLawTrialResponse(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    Vector strain_vector = ZeroVector(mpConstitutiveLaw->GetStrainSize());
    Vector stress_vector = ZeroVector(mpConstitutiveLaw->GetStrainSize());
    strain_vector[0] = CalculateGreenLagrangeStrain();


    ConstitutiveLaw::Parameters element_parameters;
    element_parameters.SetMaterialProperties(GetProperties());
    element_parameters.SetStressVector(stress_vector);
    element_parameters.SetStrainVector(strain_vector);

    mpConstitutiveLaw->CalculateMaterialResponse(element_parameters,ConstitutiveLaw::StressMeasure_PK2);

    KRATOS_CATCH("");
}

double SlidingCableElement3D::ReturnTangentModulus1D(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    double tangent_modulus(0.00);
    Vector strain_vector = ZeroVector(mpConstitutiveLaw->GetStrainSize());
    strain_vector[0] = CalculateGreenLagrangeStrain();

    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
    Values.SetStrainVector(strain_vector);

    mpConstitutiveLaw->CalculateValue(Values,TANGENT_MODULUS,tangent_modulus);
    return tangent_modulus;
    KRATOS_CATCH("");
}


void SlidingCableElement3D::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    rSerializer.save("mpConstitutiveLaw", mpConstitutiveLaw);
}
void SlidingCableElement3D::load(Serializer &rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
  rSerializer.load("mpConstitutiveLaw", mpConstitutiveLaw);
}
} // namespace Kratos.