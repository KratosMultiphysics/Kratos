//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
// ==============================================================================
//

#ifndef ROTATE_REGION_PROCESS_H
#define ROTATE_REGION_PROCESS_H

// System includes
#include <cmath>
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/variables.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "utilities/quaternion.h"

// Application includes

namespace Kratos {

/**
 * @class RotateRegionProcess
 * @ingroup KratosChimeraApplication
 * @brief This is the to apply rotation to a given modelpart
 * @details The class applies a rotation to a given modelpart given a constant
 * angular velocity or it calculates the torque applied on the modelpart and
 * calculates(by BDF2 time integration) the current theta and omega and rotates
 * the modelpart accordingly.
 * @author Aditya Ghantasala
 */
class KRATOS_API(CHIMERA_APPLICATION) RotateRegionProcess : public Process {
public:
  /// Pointer definition of MoveRotorProcess
  KRATOS_CLASS_POINTER_DEFINITION(RotateRegionProcess);

  typedef ModelPart::NodeIterator NodeIteratorType;
  typedef ProcessInfo ProcessInfoType;
  typedef ProcessInfo::Pointer ProcessInfoPointerType;
  typedef Matrix MatrixType;
  typedef Vector VectorType;

  /// Constructor.
  RotateRegionProcess(ModelPart &rModelPart, Parameters rParameters)
      : Process(Flags()), mrModelPart(rModelPart), mParameters(rParameters) {

    Parameters default_parameters(R"(
            {
                "model_part_name":"SPECIFY_MODELPART_NAME",
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
    mAxisOfRotationVector = axis_of_rotation_raw / norm;
    mTheta = 0.0;
    mToCalculateTorque = mParameters["calculate_torque"].GetBool();
    mMomentOfInertia = mParameters["moment_of_inertia"].GetDouble();
    mRotationalDamping = mParameters["rotational_damping"].GetDouble();
    KRATOS_ERROR_IF(mToCalculateTorque && mAngularVelocityRadians != 0.0)
        << "RotateRegionProcess: both \"calculate_torque\" and "
           "\"angular_velocity_radians\" cannot be specified. Please check the "
           "input."
        << std::endl;

    KRATOS_WARNING_IF("RotateRegionProcess",
                      mToCalculateTorque && mMomentOfInertia == 0.0)
        << "No moment_of_inertia specified. So no rotation possible"
        << std::endl;

    if (mToCalculateTorque) {
      mpRotationSystem = Kratos::make_shared<RotationSystem>(
          mMomentOfInertia, mRotationalDamping);
    }
  }

  /// Destructor.
  virtual ~RotateRegionProcess() {}

  void ExecuteBeforeSolutionLoop() override {}

  void SetAngularVelocity(const double NewAngularVelocity) {
    mAngularVelocityRadians = NewAngularVelocity;
  }

  void SetTorque(const double NewTorque) { mTorque = NewTorque; }

  void ExecuteInitializeSolutionStep() override {
    KRATOS_TRY;
    const auto &r_process_info = mrModelPart.GetProcessInfo();
    const int domain_size = r_process_info[DOMAIN_SIZE];
    // Does the time integration and calculates the new theta and omega
    CalculateCurrentRotationState();

    const int num_nodes = mrModelPart.NumberOfNodes();
    const NodeIteratorType it_node_begin = mrModelPart.NodesBegin();

#pragma omp parallel for schedule(guided, 512)
    for (int i_node = 0; i_node < num_nodes; ++i_node) {
      auto it_node = it_node_begin;
      std::advance(it_node, i_node);

      /// Calculating the displacement of the current node
      array_1d<double, 3> transformed_coordinates;
      TransformNode(it_node->GetInitialPosition().Coordinates(),
                    transformed_coordinates);

      it_node->X() = transformed_coordinates[0];
      it_node->Y() = transformed_coordinates[1];
      if (domain_size > 2)
        it_node->Z() = transformed_coordinates[2];

      // Computing the linear velocity at this it_node
      DenseVector<double> radius(3);
      DenseVector<double> linearVelocity(3);
      radius[0] = it_node->X() - mCenterOfRotation[0];
      radius[1] = it_node->Y() - mCenterOfRotation[1];
      radius[2] = it_node->Z() - mCenterOfRotation[2];
      CalculateLinearVelocity(mAxisOfRotationVector, radius, linearVelocity);
      if (mParameters["is_ale"].GetBool()) {
        it_node->FastGetSolutionStepValue(MESH_VELOCITY_X, 0) =
            mAngularVelocityRadians * linearVelocity[0];
        it_node->FastGetSolutionStepValue(MESH_VELOCITY_Y, 0) =
            mAngularVelocityRadians * linearVelocity[1];
        if (domain_size > 2)
          it_node->FastGetSolutionStepValue(MESH_VELOCITY_Z, 0) =
              mAngularVelocityRadians * linearVelocity[2];

        if (it_node->IsFixed(VELOCITY_X))
          it_node->FastGetSolutionStepValue(VELOCITY_X, 0) =
              it_node->FastGetSolutionStepValue(MESH_VELOCITY_X, 0);

        if (it_node->IsFixed(VELOCITY_Y))
          it_node->FastGetSolutionStepValue(VELOCITY_Y, 0) =
              it_node->FastGetSolutionStepValue(MESH_VELOCITY_Y, 0);

        if (domain_size > 2)
          if (it_node->IsFixed(VELOCITY_Z))
            it_node->FastGetSolutionStepValue(VELOCITY_Z, 0) =
                it_node->FastGetSolutionStepValue(MESH_VELOCITY_Z, 0);
      }
    }

    KRATOS_CATCH("");
  }

  void ExecuteFinalizeSolutionStep() override {}

  void ExecuteAfterOutputStep() override {}

  virtual std::string Info() const override {
    std::stringstream buffer;
    buffer << "RotateRegionProcess";
    return buffer.str();
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream &rOStream) const override {
    rOStream << "RotateRegionProcess";
  }

  /// Print object's data.
  void PrintData() {}

protected:
  ///@name Protected static Member Variables
  ///@{

  ///@}
  ///@name Protected member Variables
  ///@{

  ///@}

private:
  ///@name Private static Member Variables
  ///@{

  ///@}
  ///@name Private member Variables
  ///@{

  /**
   * @class RotationSystem
   * @ingroup KratosChimeraApplication
   * @brief This class does the solve of Euler's equations of rotation around
   * the given axis.
   * @details Uses BDF2 time integration for calculating the current theta and
   * omega
   * @author Aditya Ghantasala
   */
  class RotationSystem {
  public:
    /// Pointer definition of RotationSystem
    KRATOS_CLASS_POINTER_DEFINITION(RotationSystem);
    /*
     * @brief Constructor
     * @param MomentOfInertia The moment of inertia of the system.
     * @param DampingCoefficient The Damping Coefficient of the system.
     */
    RotationSystem(const double MomentOfInertia,
                   const double DampingCoefficient = 0.0)
        : mMomentOfInertia(MomentOfInertia), mDampingCoeff(DampingCoefficient),
          mTorque(0.0), mTime(0.0) {
      mBdf2Coeff.resize(3, false);
      mBdf2Coeff[0] = 0.0; mBdf2Coeff[1] = 0.0; mBdf2Coeff[2] = 0.0;
      mTheta.resize(3, false);
      mTheta[0] = 0.0; mTheta[1] = 0.0; mTheta[2] = 0.0;
      mOmega.resize(3, false);
      mOmega[0] = 0.0; mOmega[1] = 0.0; mOmega[2] = 0.0;
    }

    /*
     * @brief Advances the rotation system in time
     * @param NewTime Time where the system is to be set
     */
    bool CloneTimeStep(const double NewTime, const double Dt) {
      if(mTime == NewTime)
        return false;

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
      return true;
    }

    /*
     * @brief Applies the given torque on the system
     * @param Torque The torque to be applied
     */
    void inline ApplyTorque(const double Torque) { mTorque = Torque; }

    /*
     * @brief Calculates the current state of the system
     */
    void CalculateCurrentRotationState() {
      Predict();
      double LHS = ComputeLHS();
      double RHS = ComputeRHS();
      double d_theta = RHS / LHS;
      Update(d_theta);
    }

    /*
     * @brief Gives the current angle of rotation Theta
     * @out Current angle of rotation Theta
     */
    double GetCurrentTheta() const { return mTheta[0]; }

    /*
     * @brief Gives the current angular velocity Omega
     * @out Current angular velocity Omega
     */
    double GetCurrentOmega() const { return mOmega[0]; }

  private:
    double mDt;
    double mMomentOfInertia;
    double mDampingCoeff;
    double mTorque;
    double mTime;
    DenseVector<double> mBdf2Coeff;
    DenseVector<double> mTheta;
    DenseVector<double> mOmega;

    /*
     * @brief Calculates the Intertial torque acting on the system
     * @out The inertial torque
     */
    double CalculateInertiaTorque() const {
      return mMomentOfInertia *
             (mBdf2Coeff[0] * mOmega[0] + mBdf2Coeff[1] * mOmega[1] +
              mBdf2Coeff[2] * mOmega[2]);
    }

    /*
     * @brief Calculates the Damping torque acting on the system
     * @out The damping torque
     */
    double CalculateDampingTorque() const { return mDampingCoeff * mOmega[0]; }

    /*
     * @brief Computes the LHS of the Euler's equations
     * @out The LHS
     */
    double ComputeLHS() const {
      return std::pow(mBdf2Coeff[0], 2) * mMomentOfInertia +
             mBdf2Coeff[0] * mDampingCoeff;
    }

    /*
     * @brief Computes the RHS of the Euler's equations
     * @out The RHS
     */
    double ComputeRHS() const {
      double t_inertia = CalculateInertiaTorque();
      double t_damping = CalculateDampingTorque();
      double res = mTorque - t_inertia - t_damping;
      return res;
    }

    /*
     * @brief Predicts the next time steps theta
     * @out The predicted value
     */
    void Predict() {
      double d_theta = mOmega[1] * mDt;
      Update(d_theta);
    }

    /*
     * @brief Updates the state of the system with a given update
     * @param UpdateTheta The update of theta
     */
    void Update(double UpdateTheta) {
      mTheta[0] += UpdateTheta;
      mOmega[0] = mBdf2Coeff[0] * mTheta[0] + mBdf2Coeff[1] * mTheta[1] +
                  mBdf2Coeff[2] * mTheta[2];
    }
  };

  ModelPart &mrModelPart;
  Parameters mParameters;
  std::string mSubModelPartName;
  double mAngularVelocityRadians;
  array_1d<double, 3> mAxisOfRotationVector;
  array_1d<double, 3> mCenterOfRotation;
  double mTheta;
  bool mToCalculateTorque;
  double mMomentOfInertia;
  double mRotationalDamping;
  double mTorque;
  RotationSystem::Pointer mpRotationSystem;

  ///@}

  /// Assignment operator.
  RotateRegionProcess &operator=(RotateRegionProcess const &rOther) {
    return *this;
  }

  void CalculateCurrentRotationState() {
    const auto &r_process_info = mrModelPart.GetProcessInfo();
    if (mToCalculateTorque) {
      const double time = r_process_info[TIME];
      const double delta_t = r_process_info[DELTA_TIME];
      const bool is_cloned = mpRotationSystem->CloneTimeStep(time, delta_t);
      if(is_cloned){
        const double torque = CalculateTorque();
        KRATOS_INFO("RotateRegionProcess")<<"Current torque             :: "<<torque<<std::endl;
        mpRotationSystem->ApplyTorque(torque);
        mpRotationSystem->CalculateCurrentRotationState();
      }
      mTheta = mpRotationSystem->GetCurrentTheta();
      mAngularVelocityRadians = mpRotationSystem->GetCurrentOmega();
    } else {
      const double dt = r_process_info[DELTA_TIME];
      mTheta += mAngularVelocityRadians * dt;
    }
    auto& r_model = mrModelPart.GetModel();
    auto& r_torque_model_part = mParameters.Has("torque_model_part_name") ? r_model.GetModelPart(mParameters["torque_model_part_name"].GetString()) : mrModelPart;
    KRATOS_INFO("RotateRegionProcess")<<"Current angular velocity   :: "<<mAngularVelocityRadians<<std::endl;
    KRATOS_INFO("RotateRegionProcess")<<"Current angle of rotation  :: "<<mTheta<<std::endl;

    r_torque_model_part.SetValue(ROTATIONAL_ANGLE, mTheta);
    r_torque_model_part.SetValue(ROTATIONAL_VELOCITY, mAngularVelocityRadians);
  }

  /*
   * @brief Calculates the linear velocity v = r x w
   * @param rAxisOfRotationVector the axis of rotation vector
   * @param rRadius the radius vector of the point for which v is to be
   * calculated.
   * @out   rLinearVelocity the calculated linear velocity.
   */
  void CalculateLinearVelocity(const DenseVector<double> &rAxisOfRotationVector,
                               const DenseVector<double> &rRadius,
                               DenseVector<double> &rLinearVelocity) {
    assert(rAxisOfRotationVector.size() == rRadius.size());
    rLinearVelocity[0] = rAxisOfRotationVector[1] * rRadius[2] -
                         rAxisOfRotationVector[2] * rRadius[1];
    rLinearVelocity[1] = rAxisOfRotationVector[2] * rRadius[0] -
                         rAxisOfRotationVector[0] * rRadius[2];
    rLinearVelocity[2] = rAxisOfRotationVector[0] * rRadius[1] -
                         rAxisOfRotationVector[1] * rRadius[0];
  }

  /*
   * @brief Rotates the given node by mTheta around the mAxisOfRotationVector
   * and mCenterOfRotation This function uses Quaternion for rotations.
   * @param rCoordinates The nodal coordinates which are to be rotated.
   * @out rTransformedCoordinates the rotated nodal coordinates.
   */
  void TransformNode(const array_1d<double, 3> &rCoordinates,
                     array_1d<double, 3> &rTransformedCoordinates) const {
    // Changing the origin to the center of rotation
    auto new_coordinates = rCoordinates - mCenterOfRotation;
    Quaternion<double> mQuaternion = Quaternion<double>::FromAxisAngle(
        mAxisOfRotationVector(0), mAxisOfRotationVector(1),
        mAxisOfRotationVector(2), mTheta);
    mQuaternion.RotateVector3(new_coordinates, rTransformedCoordinates);
    // Bringing back to the original coordinate system.
    rTransformedCoordinates = rTransformedCoordinates + mCenterOfRotation;
  }

  /*
   * @brief Calculates Torque on the given modelpart
   * @out Torque on the modelpart about the axis of rotation.
   */
  double CalculateTorque() const {
    double torque = 0.0;
    auto& r_model = mrModelPart.GetModel();
    auto& r_torque_model_part = mParameters.Has("torque_model_part_name") ? r_model.GetModelPart(mParameters["torque_model_part_name"].GetString()) : mrModelPart;
    const int num_nodes = r_torque_model_part.NumberOfNodes();
    const NodeIteratorType it_node_begin = r_torque_model_part.NodesBegin();

#pragma omp parallel for schedule(guided, 512) reduction(+ : torque)
    for (int i_node = 0; i_node < num_nodes; ++i_node) {
      NodeIteratorType it_node = it_node_begin;
      std::advance(it_node, i_node);

      array_1d<double, 3> torque_vector;
      const array_1d<double, 3> r_vector =
          it_node->Coordinates() - mCenterOfRotation;
      auto reaction = it_node->FastGetSolutionStepValue(REACTION, 0);
      const double rho =  it_node->FastGetSolutionStepValue(DENSITY);

      torque_vector[0] = r_vector[1] * -1*reaction[2] - r_vector[2] * -1*reaction[1];
      torque_vector[1] = r_vector[2] * -1*reaction[0] - r_vector[0] * -1*reaction[2];
      torque_vector[2] = r_vector[0] * -1*reaction[1] - r_vector[1] * -1*reaction[0];

      torque += rho*(torque_vector[0] * mAxisOfRotationVector[0] +
                torque_vector[1] * mAxisOfRotationVector[1] +
                torque_vector[2] * mAxisOfRotationVector[2]);

//       torque += std::sqrt(torque_vector[0]*torque_vector[0] +
//                           torque_vector[1]*torque_vector[1] +
//                           torque_vector[2]*torque_vector[2]);
    }
    return torque;
  }
}; // Class MoveRotorProcess

}; // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
