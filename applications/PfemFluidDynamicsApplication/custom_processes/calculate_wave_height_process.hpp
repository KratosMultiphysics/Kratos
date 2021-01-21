//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics PfemFluidDynamics Application
//
//  License:         BSD License
//    Kratos default license:
//  kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_CALCULATE_WAVE_HEIGHT_PROCESS_H_INCLUDED)
#define KRATOS_CALCULATE_WAVE_HEIGHT_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
   */
class CalculateWaveHeightProcess : public Process
{
public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of CalculateWaveHeightProcess
  KRATOS_CLASS_POINTER_DEFINITION(CalculateWaveHeightProcess);

  ///@}
  ///@name Life Cycle
  ///@{
  CalculateWaveHeightProcess(ModelPart &rModelPart, 
                             int HeightDirection,
                             int PlaneDirection,
                             double PlaneCoordinates = 0.0,
                             double HeightReference = 0.0,
                             double Tolerance = 1.0e-2) : mrModelPart(rModelPart), mHeightDirection(HeightDirection),
                                                                mPlaneDirection(PlaneDirection), mPlaneCoordinates(PlaneCoordinates),
                                                                mHeightReference(HeightReference), mTolerance(Tolerance)
  {
  }

  /// Destructor.
  virtual ~CalculateWaveHeightProcess() {}

  ///@}
  ///@name Operators
  ///@{

  /// This operator is provided to call the process as a function and simply calls the Execute method.
  void operator()()
  {
    Execute();
  }

  ///@}
  ///@name Operations
  ///@{

  /// Execute method is used to execute the CalculateWaveHeightProcess.
  void Execute() override
  {
    KRATOS_TRY;


    KRATOS_CATCH("");
  }

  ///@}
  ///@name Access
  ///@{

  ///@}
  ///@name Inquiry
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  std::string Info() const override
  {
    return "CalculateWaveHeightProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream &rOStream) const override
  {
    rOStream << "CalculateWaveHeightProcess";
  }

  ///@}
  ///@name Friends
  ///@{
  ///@}

protected:
  ///@name Protected static Member Variables
  ///@{
  ///@}
  ///@name Protected member Variables
  ///@{
  ///@}
  ///@name Protected Operators
  ///@{

  /// Copy constructor.
  CalculateWaveHeightProcess(CalculateWaveHeightProcess const &rOther);

  ///@}
  ///@name Protected Operations
  ///@{
  ///@}
  ///@name Protected  Access
  ///@{
  ///@}
  ///@name Protected Inquiry
  ///@{
  ///@}
  ///@name Protected LifeCycle
  ///@{
  ///@}

private:
  ///@name Static Member Variables
  ///@{
  ///@}
  ///@name Member Variables
  ///@{

  ModelPart &mrModelPart;
  int mHeightDirection;
  int mPlaneDirection;

  double mPlaneCoordinates;
  double mHeightReference;
  double mTolerance;

  ///@}
  ///@name Private Operations
  ///@{
  ///@}
  ///@name Private  Access
  ///@{

  /// Assignment operator.
  CalculateWaveHeightProcess &operator=(CalculateWaveHeightProcess const &rOther);

  ///@}
  ///@name Serialization
  ///@{
  ///@name Private Inquiry
  ///@{
  ///@}
  ///@name Un accessible methods
  ///@{
  ///@}

}; // Class CalculateWaveHeightProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                CalculateWaveHeightProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const CalculateWaveHeightProcess &rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_TRANSFER_MODEL_PART_ELEMENTS_PROCESS_H_INCLUDED  defined
