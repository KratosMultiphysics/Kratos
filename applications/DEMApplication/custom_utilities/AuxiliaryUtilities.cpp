//
// Author:
// Miguel Angel Celigueta maceli@cimne.upc.edu
//

#include "AuxiliaryUtilities.h"
#include "DEM_application_variables.h"

// System includes
#include <string>
#include <iostream>

namespace Kratos {

    void AuxiliaryUtilities::UpdateTimeInOneModelPart(ModelPart& r_model_part, const double& time, const double& dt, const bool& is_time_to_print){
        auto& process_info = r_model_part.GetProcessInfo();
        process_info[TIME] = time;
        process_info[DELTA_TIME] = dt;
        process_info[TIME_STEPS] += 1;
        process_info[IS_TIME_TO_PRINT] = is_time_to_print;
    };

} // Namespace Kratos