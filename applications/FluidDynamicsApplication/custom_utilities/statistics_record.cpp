#include "includes/define.h"
#include "containers/variable.h"

#include "statistics_record.h"
#include "statistics_data.h"

namespace Kratos {

KRATOS_CREATE_VARIABLE( StatisticsRecord::Pointer, STATISTICS_CONTAINER)

//TODO move somewhere else
typedef typename StatisticsData<std::vector<double>>::Pointer StatisticsDataPointerType;
KRATOS_CREATE_VARIABLE( StatisticsDataPointerType, TURBULENCE_STATISTICS_DATA)

std::vector<double> StatisticsRecord::mUpdateBuffer = std::vector<double>();
std::vector<double> StatisticsRecord::mMeasurementBuffer = std::vector<double>();

}