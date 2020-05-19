
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					        Navaneeth K Narayanan
//


#include "rotate_region_process.h"
#include "utilities/math_utils.h"
#include "utilities/quaternion.h"


// Application includes

namespace Kratos
{

/// Constructor.
RotateRegionProcess::RotateRegionProcess(ModelPart &rModelPart, Parameters rParameters)
    : Process(), mrModelPart(rModelPart), mParameters(rParameters)
{

  Parameters default_parameters(R"(
            {
                "model_part_name":"SPECIFY_MODELPART_NAME",
                "torque_model_part_name":"PLEASE_SPECITY",
                "center_of_rotation":[],
                "calculate_torque":false,
                "moment_of_inertia":0.0,
                "rotational_damping":0.0,
                "angular_velocity_radians":0.0,
                "axis_of_rotation":[],
                "is_ale" : false
            }  )");

  mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

  mAngularVelocityRadians =
      mParameters["angular_velocity_radians"].GetDouble();
  mCenterOfRotation = mParameters["center_of_rotation"].GetVector();
  auto axis_of_rotation_raw = mParameters["axis_of_rotation"].GetVector();
  double norm = norm_2(axis_of_rotation_raw);
  KRATOS_ERROR_IF(norm<1e-10)<<"Norm of the Axis of rotation is close to zero. Please check the input!"<<std::endl;
  mAxisOfRotationVector = axis_of_rotation_raw / norm;
  mTheta = 0.0;
  mToCalculateTorque = mParameters["calculate_torque"].GetBool();
  KRATOS_ERROR_IF(mToCalculateTorque && mAngularVelocityRadians != 0.0)
      << "RotateRegionProcess: both \"calculate_torque\" and "
         "\"angular_velocity_radians\" cannot be specified. Please check the "
         "input."
      << std::endl;

  KRATOS_WARNING_IF("RotateRegionProcess",
                    mToCalculateTorque && mParameters["moment_of_inertia"].GetDouble() == 0.0)
      << "No moment_of_inertia specified. So no rotation possible"
      << std::endl;

  if (mToCalculateTorque)
  {
    mpRotationSystem = Kratos::make_shared<RotationSystem>(
        mParameters["moment_of_inertia"].GetDouble(), mParameters["rotational_damping"].GetDouble() );
  }
}

void RotateRegionProcess::SetAngularVelocity(const double NewAngularVelocity)
{
  mAngularVelocityRadians = NewAngularVelocity;
}

void RotateRegionProcess::ExecuteInitializeSolutionStep()
{
  KRATOS_TRY;
  const auto &r_process_info = mrModelPart.GetProcessInfo();
  const int domain_size = r_process_info[DOMAIN_SIZE];
  const double time = r_process_info[TIME];
  if (mTime == time)
    return;

  mTime = time;
  // Does the time integration and calculates the new theta and omega
  CalculateCurrentRotationState();

  const int num_nodes = static_cast<int>(mrModelPart.NumberOfNodes());
  const NodeIteratorType it_node_begin = mrModelPart.NodesBegin();

#pragma omp parallel for schedule(guided, 512)
  for (int i_node = 0; i_node < num_nodes; ++i_node)
  {
    auto it_node = it_node_begin;
    std::advance(it_node, i_node);

    /// Calculating the displacement of the current node
    array_1d<double, 3> transformed_coordinates;
    TransformNode(it_node->GetInitialPosition().Coordinates(),
                  transformed_coordinates, mTheta);

    it_node->X() = transformed_coordinates[0];
    it_node->Y() = transformed_coordinates[1];
    if (domain_size > 2)
      it_node->Z() = transformed_coordinates[2];

    noalias(it_node->FastGetSolutionStepValue(ROTATION_MESH_DISPLACEMENT)) = transformed_coordinates - it_node->GetInitialPosition().Coordinates();

    // Computing the linear velocity at this it_node
    DenseVector<double> radius(3);
    DenseVector<double> linear_velocity(3);
    radius[0] = it_node->X() - mCenterOfRotation[0];
    radius[1] = it_node->Y() - mCenterOfRotation[1];
    radius[2] = it_node->Z() - mCenterOfRotation[2];
    CalculateLinearVelocity(mAxisOfRotationVector, radius, linear_velocity);
    auto& r_rotational_mesh_vel = it_node->FastGetSolutionStepValue(ROTATION_MESH_VELOCITY, 0);

    r_rotational_mesh_vel[0] = mAngularVelocityRadians * linear_velocity[0];
    r_rotational_mesh_vel[1] = mAngularVelocityRadians * linear_velocity[1];
    if (domain_size > 2)
      r_rotational_mesh_vel[2] = mAngularVelocityRadians * linear_velocity[2];

    if (mParameters["is_ale"].GetBool())
    {
      auto& r_mesh_vel = it_node->FastGetSolutionStepValue(MESH_VELOCITY, 0);

      r_mesh_vel[0] = it_node->FastGetSolutionStepValue(ROTATION_MESH_VELOCITY_X, 0);
      r_mesh_vel[1] = it_node->FastGetSolutionStepValue(ROTATION_MESH_VELOCITY_Y, 0);
      if (domain_size > 2)
        r_mesh_vel[2] = it_node->FastGetSolutionStepValue(ROTATION_MESH_VELOCITY_Z, 0);

      auto& r_vel = it_node->FastGetSolutionStepValue(VELOCITY, 0);

      if (it_node->IsFixed(VELOCITY_X))
        r_vel[0] = it_node->FastGetSolutionStepValue(MESH_VELOCITY_X, 0);
      if (it_node->IsFixed(VELOCITY_Y))
        r_vel[1] = it_node->FastGetSolutionStepValue(MESH_VELOCITY_Y, 0);
      if (domain_size > 2)
        if (it_node->IsFixed(VELOCITY_Z))
          r_vel[2] = it_node->FastGetSolutionStepValue(MESH_VELOCITY_Z, 0);
    }
  }

  KRATOS_CATCH("");
}

std::string RotateRegionProcess::Info() const
{
  std::stringstream buffer;
  buffer << "RotateRegionProcess";
  return buffer.str();
}

void RotateRegionProcess::PrintInfo(std::ostream &rOStream) const
{
  rOStream << "RotateRegionProcess";
}

void RotateRegionProcess::PrintData() {}

RotateRegionProcess::RotationSystem::RotationSystem(const double MomentOfInertia,
                                                    const double DampingCoefficient)
    : mMomentOfInertia(MomentOfInertia), mDampingCoeff(DampingCoefficient),
      mTorque(0.0), mTime(0.0)
{
  mBdf2Coeff.resize(3, false);
  mBdf2Coeff[0] = 0.0;
  mBdf2Coeff[1] = 0.0;
  mBdf2Coeff[2] = 0.0;
  mTheta.resize(3, false);
  mTheta[0] = 0.0;
  mTheta[1] = 0.0;
  mTheta[2] = 0.0;
  mOmega.resize(3, false);
  mOmega[0] = 0.0;
  mOmega[1] = 0.0;
  mOmega[2] = 0.0;
}

void RotateRegionProcess::RotationSystem::CloneTimeStep(const double NewTime, const double Dt)
{
  mTime = NewTime;
  // compute BDF coefficients
  mDt = Dt;
  mBdf2Coeff[0] = 1.5 / mDt;
  mBdf2Coeff[1] = -2.0 / mDt;
  mBdf2Coeff[2] = 0.5 / mDt;
  // update theta database
  mTheta[2] = mTheta[1];
  mTheta[1] = mTheta[0];
  // update omega database
  mOmega[2] = mOmega[1];
  mOmega[1] = mOmega[0];
}

void RotateRegionProcess::RotationSystem::ApplyTorque(const double Torque)
{
  mTorque = Torque;
}

double RotateRegionProcess::RotationSystem::CalculateCurrentRotationState()
{
  Predict();
  double LHS = ComputeLHS();
  double RHS = ComputeRHS();
  double d_theta = RHS / LHS;
  Update(d_theta);
  return d_theta;
}

double RotateRegionProcess::RotationSystem::GetCurrentTheta() const
{
  return mTheta[0];
}

double RotateRegionProcess::RotationSystem::GetCurrentOmega() const
{
  return mOmega[0];
}

double RotateRegionProcess::RotationSystem::CalculateInertiaTorque() const
{
  return mMomentOfInertia *
         (mBdf2Coeff[0] * mOmega[0] + mBdf2Coeff[1] * mOmega[1] +
          mBdf2Coeff[2] * mOmega[2]);
}

double RotateRegionProcess::RotationSystem::CalculateDampingTorque() const
{
  return mDampingCoeff * mOmega[0];
}

double RotateRegionProcess::RotationSystem::ComputeLHS() const
{
  return std::pow(mBdf2Coeff[0], 2) * mMomentOfInertia +
         mBdf2Coeff[0] * mDampingCoeff;
}

double RotateRegionProcess::RotationSystem::ComputeRHS() const
{
  double t_inertia = CalculateInertiaTorque();
  double t_damping = CalculateDampingTorque();
  double res = mTorque - t_inertia - t_damping;
  return res;
}

void RotateRegionProcess::RotationSystem::Predict()
{
  double d_theta = mOmega[1] * mDt;
  Update(d_theta);
}

void RotateRegionProcess::RotationSystem::Update(double UpdateTheta)
{
  mTheta[0] += UpdateTheta;
  mOmega[0] = mBdf2Coeff[0] * mTheta[0] + mBdf2Coeff[1] * mTheta[1] +
              mBdf2Coeff[2] * mTheta[2];
}

void RotateRegionProcess::CalculateCurrentRotationState()
{
  const auto &r_process_info = mrModelPart.GetProcessInfo();
  if (mToCalculateTorque)
  {
    const double time = r_process_info[TIME];
    const double delta_t = r_process_info[DELTA_TIME];
    mpRotationSystem->CloneTimeStep(time, delta_t);
    const double torque = CalculateTorque();
    KRATOS_INFO("RotateRegionProcess") << "Current torque             :: " << torque << std::endl;
    mpRotationSystem->ApplyTorque(torque);
    mDTheta = mpRotationSystem->CalculateCurrentRotationState();
    mTheta = mpRotationSystem->GetCurrentTheta();
    mAngularVelocityRadians = mpRotationSystem->GetCurrentOmega();
  }
  else
  {
    const double dt = r_process_info[DELTA_TIME];
    mDTheta = mAngularVelocityRadians * dt;
    mTheta += mDTheta;
  }
  auto &r_model = mrModelPart.GetModel();
  auto &r_torque_model_part = r_model.HasModelPart(mParameters["torque_model_part_name"].GetString()) ? r_model.GetModelPart(mParameters["torque_model_part_name"].GetString()) : mrModelPart;
  KRATOS_INFO("RotateRegionProcess") << "Current angular velocity   :: " << mAngularVelocityRadians << std::endl;
  KRATOS_INFO("RotateRegionProcess") << "Current angle of rotation  :: " << mTheta << std::endl;
  KRATOS_INFO("RotateRegionProcess") << "dCurrent angle of rotation :: " << mDTheta << std::endl;

  r_torque_model_part.SetValue(ROTATIONAL_ANGLE, mTheta);
  r_torque_model_part.SetValue(ROTATIONAL_VELOCITY, mAngularVelocityRadians);
}

void RotateRegionProcess::CalculateLinearVelocity(const DenseVector<double> &rAxisOfRotationVector,
                                                  const DenseVector<double> &rRadius,
                                                  DenseVector<double> &rLinearVelocity)
{
  assert(rAxisOfRotationVector.size() == rRadius.size());
  rLinearVelocity[0] = rAxisOfRotationVector[1] * rRadius[2] -
                       rAxisOfRotationVector[2] * rRadius[1];
  rLinearVelocity[1] = rAxisOfRotationVector[2] * rRadius[0] -
                       rAxisOfRotationVector[0] * rRadius[2];
  rLinearVelocity[2] = rAxisOfRotationVector[0] * rRadius[1] -
                       rAxisOfRotationVector[1] * rRadius[0];
}

void RotateRegionProcess::TransformNode(const array_1d<double, 3> &rCoordinates,
                                        array_1d<double, 3> &rTransformedCoordinates,
                                        double Theta) const
{
  // Changing the origin to the center of rotation
  auto new_coordinates = rCoordinates - mCenterOfRotation;
  Quaternion<double> mQuaternion = Quaternion<double>::FromAxisAngle(
      mAxisOfRotationVector(0), mAxisOfRotationVector(1),
      mAxisOfRotationVector(2), Theta);
  mQuaternion.RotateVector3(new_coordinates, rTransformedCoordinates);
  // Bringing back to the original coordinate system.
  rTransformedCoordinates = rTransformedCoordinates + mCenterOfRotation;
}

double RotateRegionProcess::CalculateTorque() const
{
  double torque = 0.0;
  auto &r_model = mrModelPart.GetModel();
  auto &r_torque_model_part = r_model.HasModelPart(mParameters["torque_model_part_name"].GetString()) ? r_model.GetModelPart(mParameters["torque_model_part_name"].GetString()) : mrModelPart;
  const int num_nodes = static_cast<int>(r_torque_model_part.NumberOfNodes());
  const NodeIteratorType it_node_begin = r_torque_model_part.NodesBegin();

#pragma omp parallel for schedule(guided, 512) reduction(+ \
                                                         : torque)
  for (int i_node = 0; i_node < num_nodes; ++i_node)
  {
    NodeIteratorType it_node = it_node_begin + i_node;

    array_1d<double, 3> torque_vector;
    const array_1d<double, 3> r_vector =
        it_node->Coordinates() - mCenterOfRotation;
    auto reaction = it_node->FastGetSolutionStepValue(REACTION, 0);
    const double rho = it_node->FastGetSolutionStepValue(DENSITY);

    torque_vector[0] = r_vector[1] * -1 * reaction[2] - r_vector[2] * -1 * reaction[1];
    torque_vector[1] = r_vector[2] * -1 * reaction[0] - r_vector[0] * -1 * reaction[2];
    torque_vector[2] = r_vector[0] * -1 * reaction[1] - r_vector[1] * -1 * reaction[0];

    torque += rho * (torque_vector[0] * mAxisOfRotationVector[0] +
                     torque_vector[1] * mAxisOfRotationVector[1] +
                     torque_vector[2] * mAxisOfRotationVector[2]);
  }
  return torque;
}

}; // namespace Kratos.
