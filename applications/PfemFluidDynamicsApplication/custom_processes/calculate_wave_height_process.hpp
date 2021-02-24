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
#include "processes/process.h"
#include <fstream>
#include <iostream>

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
                               const int HeightDirection,
                               const int PlaneDirection,
                               const double PlaneCoordinates = 0.0,
                               const double HeightReference = 0.0,
                               const double Tolerance = 1.0e-2,
                               const std::string OutputFileName = "WaveHeight",
                               const double TimeInterval = 0.0) : mrModelPart(rModelPart), mHeightDirection(HeightDirection),
                                                                  mPlaneDirection(PlaneDirection), mPlaneCoordinates(PlaneCoordinates),
                                                                  mHeightReference(HeightReference), mTolerance(Tolerance),
                                                                  mOutputFileName(OutputFileName), mTimeInterval(TimeInterval)
    {
      std::ofstream my_file;
      const std::string file_name = mOutputFileName + ".txt";
      my_file.open(file_name, std::ios_base::trunc);
      my_file << "    TIME     Wave Height" << std::endl;
      my_file.close();
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
      KRATOS_TRY
      const double time = mrModelPart.GetProcessInfo()[TIME];
      const int step = mrModelPart.GetProcessInfo()[STEP];

      if (time - mPreviousPlotTime > mTimeInterval || step == 1) {
        // We loop over the nodes...
        const auto it_node_begin = mrModelPart.NodesBegin();
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<double> max_vector(num_threads, -1.0);

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
          auto it_node = it_node_begin + i;

          const int thread_id = OpenMPUtils::ThisThread();
          const auto &r_node_coordinates = it_node->Coordinates();
          if (it_node->IsNot(ISOLATED) &&
              it_node->Is(FREE_SURFACE) &&
              r_node_coordinates(mPlaneDirection) < mPlaneCoordinates + mTolerance &&
              r_node_coordinates(mPlaneDirection) > mPlaneCoordinates - mTolerance)
          {
            const double height = r_node_coordinates(mHeightDirection);

            const double wave_height = std::abs(height - mHeightReference);
            if (wave_height > max_vector[thread_id])
              max_vector[thread_id] = wave_height;
          }
        }
        const double max_height = *std::max_element(max_vector.begin(), max_vector.end());

        // We open the file where we print the wave height values
        if (max_height > -1.0) {
          std::ofstream my_file;
          const std::string file_name = mOutputFileName + ".txt";
          my_file.open(file_name, std::ios_base::app);
          my_file << "  " + std::to_string(time) + "    " + std::to_string(max_height - mHeightReference) << std::endl;
          mPreviousPlotTime = time;
        }

      }
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
    double mTimeInterval;
    double mPreviousPlotTime = 0.0;

    std::string mOutputFileName;

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
