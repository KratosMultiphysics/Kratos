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

#ifndef ROTATE_REGION_PROCESS_H
#define ROTATE_REGION_PROCESS_H

// System includes
#include <cmath>
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "processes/process.h"

// Application includes
#include "chimera_application_variables.h"

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
  RotateRegionProcess(ModelPart &rModelPart, Parameters rParameters);

  /// Destructor.
  ~RotateRegionProcess() = default;

  void SetAngularVelocity(const double NewAngularVelocity);

  void ExecuteInitializeSolutionStep() override;

  std::string Info() const override;

  /// Print information about this object.
  void PrintInfo(std::ostream &rOStream) const override;

  /// Print object's data.
  void PrintData();

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
                   const double DampingCoefficient = 0.0);

    /*
     * @brief Advances the rotation system in time
     * @param NewTime Time where the system is to be set
     */
    void CloneTimeStep(const double NewTime, const double Dt);

    /*
     * @brief Applies the given torque on the system
     * @param Torque The torque to be applied
     */
    void inline ApplyTorque(const double Torque);

    /*
     * @brief Calculates the current state of the system
     */
    double CalculateCurrentRotationState();
    /*
     * @brief Gives the current angle of rotation Theta
     * @out Current angle of rotation Theta
     */
    double GetCurrentTheta() const;

    /*
     * @brief Gives the current angular velocity Omega
     * @out Current angular velocity Omega
     */
    double GetCurrentOmega() const;

  private:
    double mDt; /// The time step (same as the timestep of the simulation)
    double mMomentOfInertia; /// Moment of inertia of the modelpart (input)
    double mDampingCoeff;    /// Rotational damping of the system (input)
    double mTorque;          /// Torque acting on the modelpart. Calculated and updated every time step.
    double mTime;            /// Current time of simulation.
    DenseVector<double> mBdf2Coeff;  /// The coefficients for BDF2 time integration. Has 3 coefficients fro BDF2.
    DenseVector<double> mTheta;      /// Current The angle of rotation theta.
    DenseVector<double> mOmega;      /// Current angular velocity Omega.

    /*
     * @brief Calculates the Intertial torque acting on the system
     * @out The inertial torque
     */
    double CalculateInertiaTorque() const;

    /*
     * @brief Calculates the Damping torque acting on the system
     * @out The damping torque
     */
    double CalculateDampingTorque() const;

    /*
     * @brief Computes the LHS of the Euler's equations
     * @out The LHS
     */
    double ComputeLHS() const;
    /*
     * @brief Computes the RHS of the Euler's equations
     * @out The RHS
     */
    double ComputeRHS() const;

    /*
     * @brief Predicts the next time steps theta
     * @out The predicted value
     */
    void Predict();

    /*
     * @brief Updates the state of the system with a given update
     * @param UpdateTheta The update of theta
     */
    void Update(double UpdateTheta);
  };

  ModelPart &mrModelPart; /// The modelpart which is to be rotated and on which the torque is to be calculated.
  Parameters mParameters; /// Input Parameters
  double mAngularVelocityRadians; /// When a rotation is applied, the angular velocity in rads/sec
  array_1d<double, 3> mAxisOfRotationVector; /// The axis of rotation around which the rotation has to be taken place.
  array_1d<double, 3> mCenterOfRotation; /// The center of rotation.
  double mTheta; /// Angle of rotation at current time step. This is calculated from beginning beginning.
  double mDTheta; /// Update in the current angle of rotation mTheta for next timestep.
  bool mToCalculateTorque; /// A bool to check if torque is to be automatically calculated on the given modelpart. This also activates the auto rotation of the modelpart. So will not work together with specifying the mAngularVelocityRadians.
  RotationSystem::Pointer mpRotationSystem; /// pointer to the object to the ODE system which does time integration (BDF2) using calculated torque and calculate mTheta and mAngularVelocityRadians.
  double mTime; /// Current Time.

  ///@}

  /// Assignment operator.
  RotateRegionProcess &operator=(RotateRegionProcess const &rOther)
  {
    return *this;
  }

  /*
   * @brief Calculates the current rotation of modelpart.
   */
  void CalculateCurrentRotationState();

  /*
   * @brief Calculates the linear velocity v = r x w
   * @param rAxisOfRotationVector the axis of rotation vector
   * @param rRadius the radius vector of the point for which v is to be
   * calculated.
   * @out   rLinearVelocity the calculated linear velocity.
   */
  void CalculateLinearVelocity(const DenseVector<double> &rAxisOfRotationVector,
                               const DenseVector<double> &rRadius,
                               DenseVector<double> &rLinearVelocity);

  /*
   * @brief Rotates the given node by mTheta around the mAxisOfRotationVector
   * and mCenterOfRotation This function uses Quaternion for rotations.
   * @param rCoordinates The nodal coordinates which are to be rotated.
   * @out rTransformedCoordinates the rotated nodal coordinates.
   */
  void TransformNode(const array_1d<double, 3> &rCoordinates,
                     array_1d<double, 3> &rTransformedCoordinates,
                     double Theta) const;
  /*
   * @brief Calculates Torque on the given modelpart
   * @out Torque on the modelpart about the axis of rotation.
   */
  double CalculateTorque() const;
}; // Class MoveRotorProcess

}; // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
