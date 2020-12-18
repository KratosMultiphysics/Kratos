//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Ruben Zorrilla
//
//


// System includes


// External includes


// Project includes
#include "includes/cfd_variables.h"
#include "utilities/element_size_calculator.h"
#include "utilities/geometry_utilities.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "estimate_dt_utilities.h"
#include "fluid_characteristic_numbers_utilities.h"


namespace Kratos
{

    void EstimateDtUtility::SetCFL(const double CFL)
    {
        mCFL = CFL;
    }

    void EstimateDtUtility::SetDtMin(const double DtMin)
    {
        mDtMin = DtMin;
    }

    void EstimateDtUtility::SetDtMax(const double DtMax)
    {
        mDtMax = DtMax;
    }

    double EstimateDtUtility::EstimateDt() const
    {
        KRATOS_TRY;

        // Get the projected element size function according to the corresponding geometry
        // Note that in here it is assumed that all the elements in the model part feature the same geometry
        const auto& r_geom = mrModelPart.ElementsBegin()->GetGeometry();
        ElementSizeFunctionType minimum_h_func = FluidCharacteristicNumbersUtilities::GetMinimumElementSizeFunction(r_geom);

        // Obtain the maximum CFL
        const double current_dt = mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);
        const double current_cfl = block_for_each<MaxReduction<double>>(mrModelPart.Elements(), [&](Element& rElement) -> double {
            return FluidCharacteristicNumbersUtilities::CalculateElementCFL(rElement, minimum_h_func, current_dt);
        });

        // Calculate the new time increment from the maximum local CFL in the mesh
        double new_dt = 0.0;
        if (current_cfl < 1e-10) {
            // Avoid division by 0 when the maximum CFL number is close to 0 (e.g. problem initialization)
            KRATOS_INFO("EstimateDtUtility") << "Setting minimum delta time " << mDtMin << " as current time step." << std::endl;
            new_dt = mDtMin;
        } else {
            // Compute new Dt
            new_dt = mCFL * current_dt / current_cfl;
            // Limit max and min Dt
            if (new_dt > mDtMax) {
                new_dt = mDtMax;
            } else if (new_dt < mDtMin) {
                new_dt = mDtMin;
            }
        }

        // Perform MPI sync if needed
        new_dt = mrModelPart.GetCommunicator().GetDataCommunicator().MinAll(new_dt);

        return new_dt;

        KRATOS_CATCH("")
    }

} // namespace Kratos.