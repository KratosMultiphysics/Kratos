#include "includes/define.h"
#include "containers/variable.h"

#include "statistics_record.h"
#include "statistics_data.h"

namespace Kratos {

KRATOS_CREATE_VARIABLE( StatisticsRecord::Pointer, TURBULENCE_STATISTICS)

//TODO move somewhere else
typedef typename StatisticsData<std::vector<double>>::Pointer StatisticsDataPointerType;
KRATOS_CREATE_VARIABLE( StatisticsDataPointerType, TURBULENCE_STATISTICS_DATA)

}