#ifndef STATIONARITY_CHECKER_H
#define STATIONARITY_CHECKER_H

#include "includes/variables.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "utilities/timer.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) StationarityChecker {

        public:

        KRATOS_CLASS_POINTER_DEFINITION(StationarityChecker);

        StationarityChecker();

        virtual ~StationarityChecker();

        bool CheckIfItsTimeToChangeGravity(ModelPart& rSpheresModelPart,
                                       const double velocity_threshold_for_gravity_change,
                                       const double min_time_between_changes,
                                       const double max_time_between_changes);

        bool CheckIfVariableIsNullInModelPart(ModelPart& rSpheresModelPart,
                                    const Variable<double>& var,
                                    const double tolerance,
                                    const bool ignore_isolated_particles);

        virtual std::string Info() const;

        virtual void PrintInfo(std::ostream& rOStream) const;

        virtual void PrintData(std::ostream& rOStream) const;

        double mPreviousChangeTime;

        protected:

        private:

        StationarityChecker & operator=(StationarityChecker const& rOther);

    }; // Class StationarityChecker
} // namespace Kratos

#endif // STATIONARITY_CHECKER_H