#include "includes/define.h"
#include "includes/element.h"
#include "containers/variable.h"

#include "statistics_record.h"
#include "statistics_data.h"

#include "fluid_dynamics_application_variables.h"

namespace Kratos {

void StatisticsRecord::AddResult(StatisticsSampler::Pointer pResult)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mInitialized) << "Trying to add statistical data after Initialization of the internal storage." << std::endl;

    std::size_t result_size = pResult->GetSize();
    pResult->SetOffset(mDataBufferSize);

    mDataBufferSize += result_size;
    mAverageData.push_back(pResult);

    KRATOS_CATCH("")
}

void StatisticsRecord::InitializeStorage()
{
    mUpdateBuffer.resize(mDataBufferSize);
    mInitialized = true;
}

void StatisticsRecord::SampleIntegrationPointResults(ModelPart& rModelPart)
{
    mRecordedSteps++;

    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    std::vector<double> dummy;
    for( auto it_elem = rModelPart.ElementsBegin(); it_elem != rModelPart.ElementsEnd(); ++it_elem)
    {
        it_elem->GetValueOnIntegrationPoints(UPDATE_STATISTICS,dummy,r_process_info);
    }
}

void StatisticsRecord::UpdateStatistics(Element* pElement)
{
    //KRATOS_DEBUG_ERROR_IF(!pElement->Has(TURBULENCE_STATISTICS_DATA)) << "Trying to compute turbulent statistics, but " << pElement->Info() << " does not have TURBULENCE_STATISTICS_DATA defined." << std::endl;
    auto &r_elemental_statistics = pElement->GetValue(TURBULENCE_STATISTICS_DATA);
    r_elemental_statistics.UpdateMeasurement(pElement, mAverageData, mUpdateBuffer, mRecordedSteps);
    //r_elemental_statistics.CalculateUpdateDelta(pElement, mAverageData, mMeasurementBuffer, mUpdateBuffer, mRecordedSteps);
    //r_elemental_statistics.UpdateMeasurement(mMeasurementBuffer);
}

void StatisticsRecord::FinalizeStatistics(ModelPart::ElementsContainerType& rElements)
{
    for (auto it_element = rElements.begin(); it_element != rElements.end(); ++it_element )
    {
        auto& r_statistics = it_element->GetValue(TURBULENCE_STATISTICS_DATA);
        for (std::size_t g = 0; g < r_statistics.NumberOfIntegrationPoints(); g++)
        {
            auto data_iterator = r_statistics.DataIterator(g);
            for (auto it_statistic = mAverageData.begin(); it_statistic != mAverageData.end(); ++it_statistic)
            {
                (*it_statistic)->Finalize(data_iterator, mRecordedSteps);
            }
        }
    }
}


std::vector<double> StatisticsRecord::OutputForTest(ModelPart::ElementsContainerType& rElements)
{
    std::vector<double> result;
    for (auto it_element = rElements.begin(); it_element != rElements.end(); ++it_element )
    {
        auto& r_statistics = it_element->GetValue(TURBULENCE_STATISTICS_DATA);
        for (std::size_t g = 0; g < r_statistics.NumberOfIntegrationPoints(); g++)
        {
            auto data_iterator = r_statistics.DataIterator(g);
            for (auto it = data_iterator.begin(); it != data_iterator.end(); ++it)
            {
                result.push_back(*it);
            }
        }
    }
}

KRATOS_CREATE_VARIABLE( StatisticsRecord::Pointer, STATISTICS_CONTAINER)

//TODO move somewhere else
KRATOS_CREATE_VARIABLE( StatisticsData, TURBULENCE_STATISTICS_DATA)

std::vector<double> StatisticsRecord::mUpdateBuffer = std::vector<double>();

}