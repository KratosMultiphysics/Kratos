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
#include "custom_elements/ring_element_ROCCO_3D.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "includes/checks.h"
#include "utilities/atomic_utilities.h"


namespace Kratos {
RingElementROCCO3D::RingElementROCCO3D(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

RingElementROCCO3D::RingElementROCCO3D(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
RingElementROCCO3D::Create(IndexType NewId, NodesArrayType const &rThisNodes,
                         PropertiesType::Pointer pProperties) const {
  const GeometryType &rGeom = this->GetGeometry();
  return Kratos::make_intrusive<RingElementROCCO3D>(NewId, rGeom.Create(rThisNodes),
                                               pProperties);
}

Element::Pointer
RingElementROCCO3D::Create(IndexType NewId, GeometryType::Pointer pGeom,
                         PropertiesType::Pointer pProperties) const {
  return Kratos::make_intrusive<RingElementROCCO3D>(NewId, pGeom,
                                               pProperties);
}

RingElementROCCO3D::~RingElementROCCO3D() {}

void RingElementROCCO3D::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    mDeformedLength = GetProperties()[RING_LIMIT_CIRCUMFERENCE];
    mMaxDeformedLength = GetProperties()[RING_LIMIT_CIRCUMFERENCE];
    mDeformedDiagonal = GetDiagonalLengthArray();
    mMaxDeformedDiagonal = GetDiagonalLengthArray();
    KRATOS_CATCH("")
}

void RingElementROCCO3D::EquationIdVector(EquationIdVectorType &rResult,
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
void RingElementROCCO3D::GetDofList(DofsVectorType &rElementalDofList,
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


void RingElementROCCO3D::GetValuesVector(Vector &rValues, int Step) const {

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

void RingElementROCCO3D::GetFirstDerivativesVector(Vector &rValues, int Step) const {

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

void RingElementROCCO3D::GetSecondDerivativesVector(Vector &rValues, int Step) const {

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


void RingElementROCCO3D::CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY;
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;
  // resizing the matrices + create memory for LHS
  rLeftHandSideMatrix = ZeroMatrix(local_size, local_size); // not implemented yet -> explicit calculations
  KRATOS_CATCH("")
}



void RingElementROCCO3D::InternalForcesCircumference(VectorType &rRightHandSideVector)
{ 
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;

  const double k_t = GetProperties()[RING_TENSILE_STIFFNESS]; 
  const double l_min = GetProperties()[RING_LIMIT_CIRCUMFERENCE]; 
  const double l_current = sum(GetCurrentLengthCircumferenceArray());


  if (l_current > l_min) 
  {
    const double delta_l = l_current - mDeformedLength;
  
    if (delta_l > 0.0)
    {
      // Loading
      if (l_current >= mMaxDeformedLength) mPerimeterForce += delta_l * k_t;
      else mPerimeterForce += delta_l * 3.0 * k_t;
    }
    else 
    {
      // Unloading
      mPerimeterForce += delta_l * 3.0 * k_t;
    }

    // Update
    mDeformedLength = l_current;
    if (l_current >= mMaxDeformedLength) mMaxDeformedLength = l_current;

  }

  

  noalias(rRightHandSideVector) -= GetDirectionVectorCircumference()*mPerimeterForce;

}

void RingElementROCCO3D::InternalForcesDiagonal(VectorType &rRightHandSideVector)
{ 
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;

  const double k_b = GetProperties()[RING_BENDING_STIFFNESS]; 
  const double r_t = GetProperties()[RING_THICKNESS_WIRE]; 
  const double r_n = GetProperties()[RING_NR_WIRES]; 
  const double l_min = GetProperties()[RING_LIMIT_CIRCUMFERENCE]; 
  const double d_min = GetProperties()[RING_DIAMETER]; 


  // only for 4noded element
  
  const Vector diagonals = GetDiagonalLengthArray();


  Vector forces_diagonal = ZeroVector(2);

  for (SizeType i=0;i<2;++i)
  {
      const double delta_diagonal = diagonals[i] - mDeformedDiagonal[i];

      if (diagonals[i] >= d_min)
      {

        if (delta_diagonal>0.0) 
        {
          // Loading
          if (diagonals[i] >= mMaxDeformedDiagonal[i]) mDiagonalForces[i] += delta_diagonal * k_b;
          else mDiagonalForces[i] += delta_diagonal * 100.0 * k_b;
        }
        else 
        {
          // Unloading
          mDiagonalForces[i] += delta_diagonal * 100.0 * k_b;
        }

        // Update
        mDeformedDiagonal[i] = diagonals[i];
        if (diagonals[i] >= mMaxDeformedDiagonal[i]) mMaxDeformedDiagonal[i] = diagonals[i];


      } 
  }

  forces_diagonal = mDiagonalForces;

  Vector global_diagonal_forces = GetDirectionVectorDiagonal();

  for (SizeType i=0;i<3;++i)
  {
    global_diagonal_forces[i]   *= forces_diagonal[0];
    global_diagonal_forces[i+6] *= forces_diagonal[0];

    global_diagonal_forces[i+3] *= forces_diagonal[1];
    global_diagonal_forces[i+9] *= forces_diagonal[1];
  }  


  noalias(rRightHandSideVector) -= global_diagonal_forces;

}

void RingElementROCCO3D::CalculateRightHandSide(
    VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
{
  KRATOS_TRY;
  const int points_number = GetGeometry().PointsNumber();
  const int dimension = 3;
  const SizeType local_size = dimension*points_number;
  rRightHandSideVector = ZeroVector(local_size);

  InternalForcesCircumference(rRightHandSideVector);
  InternalForcesDiagonal(rRightHandSideVector);
  // + body forces

  KRATOS_CATCH("")
}

void RingElementROCCO3D::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                         VectorType &rRightHandSideVector,
                                         const ProcessInfo &rCurrentProcessInfo)
{
  KRATOS_TRY;
  CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
  CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
  KRATOS_CATCH("")
}

Vector RingElementROCCO3D::GetRefLengthCircumferenceArray() const
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


Vector RingElementROCCO3D::GetCurrentLengthCircumferenceArray() const
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

Vector RingElementROCCO3D::DistanceVectorNodes(const int node_a, const int node_b) const
{
    const int dimension = 3;
    Vector distance = ZeroVector(dimension);

    distance[0] = 
        this->GetGeometry()[node_a].FastGetSolutionStepValue(DISPLACEMENT_X) -
        this->GetGeometry()[node_b].FastGetSolutionStepValue(DISPLACEMENT_X) +
        this->GetGeometry()[node_a].X0() -
        this->GetGeometry()[node_b].X0();

    distance[1] = 
        this->GetGeometry()[node_a].FastGetSolutionStepValue(DISPLACEMENT_Y) -
        this->GetGeometry()[node_b].FastGetSolutionStepValue(DISPLACEMENT_Y) +
        this->GetGeometry()[node_a].Y0() -
        this->GetGeometry()[node_b].Y0();

    distance[2] = 
        this->GetGeometry()[node_a].FastGetSolutionStepValue(DISPLACEMENT_Z) -
        this->GetGeometry()[node_b].FastGetSolutionStepValue(DISPLACEMENT_Z) +
        this->GetGeometry()[node_a].Z0() -
        this->GetGeometry()[node_b].Z0();

    return distance;

}

Vector RingElementROCCO3D::GetDirectionVectorDiagonal() const
{
    const int points_number = GetGeometry().PointsNumber();
    const int dimension = 3;
    const SizeType local_size = dimension*points_number;

    Vector n_t = ZeroVector(local_size);

    // only 4 nodes yet

    Vector distance = DistanceVectorNodes(0,2);
    distance /= norm_2(distance);

    for (SizeType i=0;i<dimension; i++)
    {
      n_t[i] = distance[i];
      n_t[i+6] = -distance[i];
    }

    distance = DistanceVectorNodes(1,3);
    distance /= norm_2(distance);

    for (SizeType i=0;i<dimension; i++)
    {
      n_t[i+3] = distance[i];
      n_t[i+9] = -distance[i];
    }

  return n_t;
}


Vector RingElementROCCO3D::GetDirectionVectorCircumference() const
{
    const int points_number = GetGeometry().PointsNumber();
    const int dimension = 3;
    const SizeType local_size = dimension*points_number;

    Vector n_t = ZeroVector(local_size);
    const Vector d_x = this->GetDeltaPositions(1);
    const Vector d_y = this->GetDeltaPositions(2);
    const Vector d_z = this->GetDeltaPositions(3);
    const Vector lengths = this->GetCurrentLengthCircumferenceArray();

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

Vector RingElementROCCO3D::GetDiagonalLengthArray(const int step) const {

  const int points_number = GetGeometry().PointsNumber();

  // only 4 noded so far
  Vector diagonals = ZeroVector(2);
  for (int i=0;i<2;++i)
  {
    int next_node_id = i+2;

    const double du =
        this->GetGeometry()[next_node_id].FastGetSolutionStepValue(DISPLACEMENT_X,step) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X,step);
    const double dv =
        this->GetGeometry()[next_node_id].FastGetSolutionStepValue(DISPLACEMENT_Y,step) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y,step);
    const double dw =
        this->GetGeometry()[next_node_id].FastGetSolutionStepValue(DISPLACEMENT_Z,step) -
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Z,step);
    const double dx = this->GetGeometry()[next_node_id].X0() - this->GetGeometry()[i].X0();
    const double dy = this->GetGeometry()[next_node_id].Y0() - this->GetGeometry()[i].Y0();
    const double dz = this->GetGeometry()[next_node_id].Z0() - this->GetGeometry()[i].Z0();



    diagonals[i] = std::sqrt((du + dx) * (du + dx) + (dv + dy) * (dv + dy) +
                             (dw + dz) * (dw + dz));
  }

  return diagonals;
}


Vector RingElementROCCO3D::GetDeltaPositions(const int& rDirection) const
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


void RingElementROCCO3D::CalculateLumpedMassVector(
    VectorType &rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const int points_number = GetGeometry().PointsNumber();
    const int dimension = 3;
    const SizeType local_size = dimension*points_number;

    // Clear matrix
    if (rLumpedMassVector.size() != local_size) {
      rLumpedMassVector.resize(local_size);
    }
    rLumpedMassVector = ZeroVector(local_size);

    const double total_mass = GetProperties()[RING_MASS];


    const Vector l_array = GetCurrentLengthCircumferenceArray();
    const double l_current = sum(l_array);

    for (IndexType i = 0; i < points_number; ++i)
    {
      int next_node_id = i+1;
      if (i==points_number-1) next_node_id = 0;

      const double dl = l_array[i]/l_current;


      for (IndexType j = 0; j < dimension; ++j)
      {
        rLumpedMassVector[i*dimension + j] += 0.50 * total_mass * dl;
        rLumpedMassVector[next_node_id*dimension + j] += 0.50 * total_mass * dl;
      }
    }

    KRATOS_CATCH("")
}

void RingElementROCCO3D::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;
    const int points_number = GetGeometry().PointsNumber();
    const int dimension = 3;
    const SizeType local_size = dimension*points_number;

    // Compute lumped mass matrix
    VectorType temp_vector(local_size);
    CalculateLumpedMassVector(temp_vector, rCurrentProcessInfo);

    // Clear matrix
    if (rMassMatrix.size1() != local_size || rMassMatrix.size2() != local_size)
        rMassMatrix.resize( local_size, local_size, false );
    rMassMatrix = ZeroMatrix(local_size, local_size);

    // Fill the matrix
    for (IndexType i = 0; i < local_size; ++i)
        rMassMatrix(i, i) = temp_vector[i];

    KRATOS_CATCH("")
}

void RingElementROCCO3D::CalculateDampingMatrix(
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

void RingElementROCCO3D::AddExplicitContribution(
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
        this->CalculateLumpedMassVector(element_mass_vector, rCurrentProcessInfo);

        for (int i = 0; i < points_number; ++i) {
            double &r_nodal_mass = r_geom[i].GetValue(NODAL_MASS);
            int index = i * dimension;

            AtomicAdd(r_nodal_mass, element_mass_vector(index));
        }
    }
    KRATOS_CATCH("")
}

void RingElementROCCO3D::AddExplicitContribution(
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
                AtomicAdd(r_force_residual[j], rRHSVector[index + j] - damping_residual_contribution[index + j] );
            }
        }
    } else if (rDestinationVariable == NODAL_INERTIA) {

        // Getting the vector mass
        VectorType mass_vector(local_size);
        CalculateLumpedMassVector(mass_vector, rCurrentProcessInfo);

        for (int i = 0; i < points_number; ++i) {
            double &r_nodal_mass = GetGeometry()[i].GetValue(NODAL_MASS);
            int index = i * dimension;

            AtomicAdd(r_nodal_mass, mass_vector[index]);
        }
    }
    KRATOS_CATCH("")
}

int RingElementROCCO3D::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF( this->Id() < 1 ) << "Element found with Id " << this->Id() << std::endl;

    KRATOS_ERROR_IF_NOT((GetGeometry().PointsNumber()==4) || (GetGeometry().PointsNumber()==3)) << "the ring element " << this->Id() << "does not have 3 || 4 nodes" << std::endl;

    KRATOS_ERROR_IF((GetGeometry().WorkingSpaceDimension() != 3)) << "The ring elements works only in 3D" << std::endl;

    if (GetProperties().Has(RING_MASS) == false) {
        KRATOS_ERROR << "RING_MASS not provided for this element" << Id()
                     << std::endl;
    }

    if (GetProperties().Has(RING_THICKNESS_WIRE) == false) {
        KRATOS_ERROR << "RING_THICKNESS_WIRE not provided for this element" << Id()
                     << std::endl;
    }

    if (GetProperties().Has(RING_NR_WIRES) == false) {
        KRATOS_ERROR << "RING_THICKNESS_WIRE not provided for this element" << Id()
                     << std::endl;
    }

    if (GetProperties().Has(RING_LIMIT_CIRCUMFERENCE) == false) {
        KRATOS_ERROR << "RING_LIMIT_CIRCUMFERENCE not provided for this element" << Id()
                     << std::endl;
    }

    if (GetProperties().Has(RING_DIAMETER) == false) {
        KRATOS_ERROR << "RING_DIAMETER not provided for this element" << Id()
                     << std::endl;
    }




    if (GetProperties().Has(RING_TENSILE_STIFFNESS) == false) {
        KRATOS_ERROR << "RING_TENSILE_STIFFNESS not provided for this element" << Id()
                     << std::endl;
    }

    if (GetProperties().Has(RING_BENDING_STIFFNESS) == false) {
        KRATOS_ERROR << "RING_BENDING_STIFFNESS not provided for this element" << Id()
                     << std::endl;
    }

    const double l_min = GetProperties()[RING_LIMIT_CIRCUMFERENCE]; 
    const double l_ref = sum(GetRefLengthCircumferenceArray());
    KRATOS_ERROR_IF(l_ref > l_min) << "The ring elements dimensions do not fit the given reference circumference: " << l_min << " < " << l_ref <<  std::endl;

    return 0;

    KRATOS_CATCH("")
}

bool RingElementROCCO3D::HasSelfWeight() const
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

Vector RingElementROCCO3D::CalculateBodyForces() {

    const int points_number = GetGeometry().PointsNumber();
    const int dimension = 3;
    const SizeType local_size = dimension*points_number;
    Vector body_forces_global = ZeroVector(local_size);
    return body_forces_global;
}


void RingElementROCCO3D::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    rOutput.resize(1);

    // temporary variable, delete this later
    if (rVariable == TEMPERATURE) {
        rOutput[0] = sum(GetCurrentLengthCircumferenceArray());
    }
}


void RingElementROCCO3D::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}
void RingElementROCCO3D::load(Serializer &rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}
} // namespace Kratos.