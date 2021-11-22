//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "calculate_wave_height_utility.h"


namespace Kratos
{

CalculateWaveHeightUtility::CalculateWaveHeightUtility(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
) : mrModelPart(rThisModelPart)
{
    Parameters default_parameters(R"({
        "automatic_time_step"   : true,
        "adaptive_time_step"    : true,
        "time_step"             : 1.0,
        "courant_number"        : 1.0,
        "minimum_delta_time"    : 1e-4,
        "maximum_delta_time"    : 1e+6
    })");

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mEstimateDt = ThisParameters["automatic_time_step"].GetBool();
    mAdaptiveDt = ThisParameters["adaptive_time_step"].GetBool();
    mConstantDt = ThisParameters["time_step"].GetDouble();
    mCourant = ThisParameters["courant_number"].GetDouble();
    mMinDt = ThisParameters["minimum_delta_time"].GetDouble();
    mMaxDt = ThisParameters["maximum_delta_time"].GetDouble();

    if (mEstimateDt && !mAdaptiveDt) {
        mConstantDt = EstimateTimeStep();
    }
}

double CalculateWaveHeightUtility::Execute() const
{
      KRATOS_TRY
      const double time = mrModelPart.GetProcessInfo()[TIME];
      const int step = mrModelPart.GetProcessInfo()[STEP];

      if (time - mPreviousPlotTime > mTimeInterval || step == 1)
      {
        // We loop over the nodes...
        const auto it_node_begin = mrModelPart.NodesBegin();
        const int num_threads = ParallelUtilities::GetNumThreads();
        // std::vector<double> max_vector(num_threads, -1.0);

        double counter = 0;
        double heightSum = 0;

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
          auto it_node = it_node_begin + i;

          const int thread_id = OpenMPUtils::ThisThread();
          const auto &r_node_coordinates = it_node->Coordinates();
          if (it_node->IsNot(ISOLATED) && it_node->IsNot(RIGID) &&
              it_node->Is(FREE_SURFACE) &&
              r_node_coordinates(mPlaneDirection) < (mPlaneCoordinates + mTolerance) &&
              r_node_coordinates(mPlaneDirection) > (mPlaneCoordinates - mTolerance))
          {
            const double height = r_node_coordinates(mHeightDirection);
            // const double wave_height = std::abs(height - mHeightReference);
            const double wave_height = height - mHeightReference;
            counter += 1.0;
            heightSum += wave_height;

            // if (wave_height > max_vector[thread_id])
            //   max_vector[thread_id] = wave_height;
          }
        }
        // const double max_height = *std::max_element(max_vector.begin(), max_vector.end());
        const double max_height = heightSum / counter;
        // We open the file where we print the wave height values
        if (max_height > -1.0)
        {
          std::ofstream my_file;
          const std::string file_name = mOutputFileName + ".txt";
          my_file.open(file_name, std::ios_base::app);
          my_file << "  " + std::to_string(time) + "    " + std::to_string(max_height - mHeightReference) << std::endl;
          mPreviousPlotTime = time;
        }
      }
      KRATOS_CATCH("");
}

}  // namespace Kratos.
