#include "custom_processes/neighbours_search_process.h"

namespace Kratos
{
    void NeighboursSearchProcess::Execute()
    {
        KRATOS_TRY
        const double h = NeighboursSearchUtilities::ComputeSmoothingLength(mrThisModelPart, mrThisParameters["coefficient"].GetDouble());
        const double radius = 2.0 * h;
        mrThisModelPart.GetProcessInfo().SetValue(SMOOTHING_LENGTH, h);

        NeighboursSearchUtilities::SearchNeighbours(mrThisModelPart, radius);
        KRATOS_CATCH("")
    }

    void NeighboursSearchProcess::ExecuteInitialize()
    {
        this->Execute();
        KRATOS_INFO("NeighboursSearch") << "The neighbours search process was executed" << std::endl;
    }
}