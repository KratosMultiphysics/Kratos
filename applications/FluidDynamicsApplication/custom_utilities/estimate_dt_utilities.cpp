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
    /// Local flags to determine the magnitudes for the Dt estimation
    KRATOS_CREATE_LOCAL_FLAG(EstimateDtUtility, CFL_ESTIMATION, 0);
    KRATOS_CREATE_LOCAL_FLAG(EstimateDtUtility, FOURIER_VISCOSITY_ESTIMATION, 1);
    KRATOS_CREATE_LOCAL_FLAG(EstimateDtUtility, FOURIER_CONDUCTIVITY_ESTIMATION, 2);

    void EstimateDtUtility::SetCFL(const double CFL)
    {
        mCFL = CFL;
    }

    void EstimateDtUtility::SetPecletViscosity(const double PecletViscosity)
    {
        mPecletViscosity = PecletViscosity;
    }

    void EstimateDtUtility::SetPecletConductivity(const double PecletConductivity)
    {
        mPecletConductivity = PecletConductivity;
    }

    void EstimateDtUtility::SetDtMin(const double DtMin)
    {
        mDtMin = DtMin;
    }

    void EstimateDtUtility::SetDtMax(const double DtMax)
    {
        mDtMax = DtMax;
    }

    void EstimateDtUtility::SetDtEstimationMagnitudesFlag()
    {
        mDtEstimationMagnitudesFlags = Flags();
        if (mCFL > 0.0) {
            mDtEstimationMagnitudesFlags.Set(CFL_ESTIMATION);
        }
        if (mPecletViscosity > 0.0) {
            mDtEstimationMagnitudesFlags.Set(FOURIER_VISCOSITY_ESTIMATION);
        }
        if (mPecletConductivity > 0.0) {
            mDtEstimationMagnitudesFlags.Set(FOURIER_CONDUCTIVITY_ESTIMATION);
        }
    }

    template<>
    double EstimateDtUtility::InternalEstimateDt<true,false,false>() const
    {
        KRATOS_TRY;

        // Get the minimum element size function according to the corresponding geometry for the CFL calculation
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

    template<>
    double EstimateDtUtility::InternalEstimateDt<true,true,true,true>() const
    {
        KRATOS_TRY;

        // Get the minimum element size function according to the corresponding geometry for the CFL calculation
        // Note that in here it is assumed that all the elements in the model part feature the same geometry
        const auto& r_geom = mrModelPart.ElementsBegin()->GetGeometry();
        ElementSizeFunctionType minimum_h_func = FluidCharacteristicNumbersUtilities::GetMinimumElementSizeFunction(r_geom);
        ElementSizeFunctionType average_h_func = FluidCharacteristicNumbersUtilities::GetAverageElementSizeFunction(r_geom);

        //TODO: CREATE A FUNCTION CALL TO CHECK IF THE DENSITY IS NODAL AND IF THE ARTIFICIAL MAGNITUDES ARE TO BE ADDED

        // Obtain the maximum CFL and Peclet numbers
        double max_CFL, max_Pe_mu, max_Pe_k;
        const double current_dt = mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);
        typedef CombinedReduction<MaxReduction<double>, MaxReduction<double>, MaxReduction<double>> CombinedMaxReduction;
        std::tie(max_CFL, max_Pe_mu, max_Pe_k) = block_for_each<CombinedMaxReduction>(mrModelPart.Elements(), [&](Element& rElement){
            const double CFL = FluidCharacteristicNumbersUtilities::CalculateElementCFL(rElement, minimum_h_func, current_dt);
            //TODO: IN HERE WE SHOULD CHECK IF THE ARTIFICIAL MAGNITUDES ARE TO BE ADDED OR NOT
            const auto Pe_numbers = FluidCharacteristicNumbersUtilities::CalculateElementPecletNumbers<true,true>(rElement, average_h_func);
            return std::make_tuple(CFL, std::get<0>(Pe_numbers), std::get<1>(Pe_numbers));
        });

        // Calculate the new time increment from the maximum local CFL in the mesh
        double new_dt = 0.0;
        // if (current_cfl < 1e-10) {
        //     // Avoid division by 0 when the maximum CFL number is close to 0 (e.g. problem initialization)
        //     KRATOS_INFO("EstimateDtUtility") << "Setting minimum delta time " << mDtMin << " as current time step." << std::endl;
        //     new_dt = mDtMin;
        // } else {
        //     // Compute new Dt
        //     new_dt = mCFL * current_dt / current_cfl;
        //     // Limit max and min Dt
        //     if (new_dt > mDtMax) {
        //         new_dt = mDtMax;
        //     } else if (new_dt < mDtMin) {
        //         new_dt = mDtMin;
        //     }
        // }

        // Perform MPI sync if needed
        new_dt = mrModelPart.GetCommunicator().GetDataCommunicator().MinAll(new_dt);

        return new_dt;

        KRATOS_CATCH("")
    }

    double EstimateDtUtility::EstimateDt() const
    {
        KRATOS_TRY;

        double new_dt = 0.0;

        if (mDtEstimationMagnitudesFlags.Is(CFL_ESTIMATION) && 
            mDtEstimationMagnitudesFlags.IsNot(FOURIER_VISCOSITY_ESTIMATION) &&
            mDtEstimationMagnitudesFlags.IsNot(FOURIER_CONDUCTIVITY_ESTIMATION)) {
            new_dt = InternalEstimateDt<true, false, false>();
        }

        //USE FOURIER NUMBER
        //FO=ALPHA*DT/dX^2 = CFL/PE

        return new_dt;

        KRATOS_CATCH("")
    }

} // namespace Kratos.