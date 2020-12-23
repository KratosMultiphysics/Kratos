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
    KRATOS_CREATE_LOCAL_FLAG(EstimateDtUtility, VISCOUS_FOURIER_ESTIMATION, 1);
    KRATOS_CREATE_LOCAL_FLAG(EstimateDtUtility, THERMAL_FOURIER_ESTIMATION, 2);

    void EstimateDtUtility::SetCFL(const double CFL)
    {
        mCFL = CFL;
    }

    void EstimateDtUtility::SetViscousFourier(const double ViscousFourier)
    {
        mViscousFourier = ViscousFourier;
    }

    void EstimateDtUtility::SetThermalFourier(const double ThermalFourier)
    {
        mThermalFourier = ThermalFourier;
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
        if (mViscousFourier > 0.0) {
            mDtEstimationMagnitudesFlags.Set(VISCOUS_FOURIER_ESTIMATION);
        }
        if (mThermalFourier > 0.0) {
            mDtEstimationMagnitudesFlags.Set(THERMAL_FOURIER_ESTIMATION);
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
        const double new_dt = CalculateNewDeltaTime(
            current_dt,
            std::make_pair(current_cfl, mCFL));

        return new_dt;

        KRATOS_CATCH("")
    }

    template<>
    double EstimateDtUtility::InternalEstimateDt<true,false,true>() const
    {
        KRATOS_TRY;

        // Get the minimum element size function according to the corresponding geometry for the CFL calculation
        // Note that in here it is assumed that all the elements in the model part feature the same geometry
        const auto& r_geom = mrModelPart.ElementsBegin()->GetGeometry();
        ElementSizeFunctionType minimum_h_func = FluidCharacteristicNumbersUtilities::GetMinimumElementSizeFunction(r_geom);

        // Set the function to calculate the Fourier numbers
        std::function<double(const Element&, const ElementSizeFunctionType&, const double)> thermal_fourier_number_function;
        if (mConsiderArtificialDiffusion) {
            if (mNodalDensityFormulation) {
                thermal_fourier_number_function = FluidCharacteristicNumbersUtilities::CalculateElementThermalFourierNumber<true,true>;
            } else {
                thermal_fourier_number_function = FluidCharacteristicNumbersUtilities::CalculateElementThermalFourierNumber<true,false>;
            }
        } else {
            if (mNodalDensityFormulation) {
                thermal_fourier_number_function = FluidCharacteristicNumbersUtilities::CalculateElementThermalFourierNumber<false,true>;
            } else {
                thermal_fourier_number_function = FluidCharacteristicNumbersUtilities::CalculateElementThermalFourierNumber<false,false>;
            }
        }

        // Obtain the maximum CFL and Peclet numbers
        double max_CFL, max_Fo_k;
        const double current_dt = mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);
        typedef CombinedReduction<MaxReduction<double>, MaxReduction<double>> CombinedMaxReduction;
        std::tie(max_CFL, max_Fo_k) = block_for_each<CombinedMaxReduction>(mrModelPart.Elements(), [&](Element& rElement){
            const double CFL = FluidCharacteristicNumbersUtilities::CalculateElementCFL(rElement, minimum_h_func, current_dt);
            const double thermal_Fo_number = thermal_fourier_number_function(rElement, minimum_h_func, current_dt);
            return std::make_tuple(CFL, thermal_Fo_number);
        });

        // Calculate the new time increment from the maximum characteristic numbers in the mesh
        const double new_dt = CalculateNewDeltaTime(
            current_dt,
            std::make_pair(max_CFL, mCFL),
            std::make_pair(max_Fo_k, mThermalFourier));

        return new_dt;

        KRATOS_CATCH("")
    }

    template<>
    double EstimateDtUtility::InternalEstimateDt<true,true,true>() const
    {
        KRATOS_TRY;

        // Get the minimum element size function according to the corresponding geometry for the CFL calculation
        // Note that in here it is assumed that all the elements in the model part feature the same geometry
        const auto& r_geom = mrModelPart.ElementsBegin()->GetGeometry();
        ElementSizeFunctionType minimum_h_func = FluidCharacteristicNumbersUtilities::GetMinimumElementSizeFunction(r_geom);

        // Set the function to calculate the Fourier numbers
        std::function<std::pair<double,double>(const Element&, const ElementSizeFunctionType&, const double)> fourier_numbers_function;
        if (mConsiderArtificialDiffusion) {
            if (mNodalDensityFormulation) {
                fourier_numbers_function = FluidCharacteristicNumbersUtilities::CalculateElementFourierNumbers<true,true>;
            } else {
                fourier_numbers_function = FluidCharacteristicNumbersUtilities::CalculateElementFourierNumbers<true,false>;
            }
        } else {
            if (mNodalDensityFormulation) {
                fourier_numbers_function = FluidCharacteristicNumbersUtilities::CalculateElementFourierNumbers<false,true>;
            } else {
                fourier_numbers_function = FluidCharacteristicNumbersUtilities::CalculateElementFourierNumbers<false,false>;
            }
        }

        // Obtain the maximum CFL and Peclet numbers
        double max_CFL, max_Fo_mu, max_Fo_k;
        const double current_dt = mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);
        typedef CombinedReduction<MaxReduction<double>, MaxReduction<double>, MaxReduction<double>> CombinedMaxReduction;
        std::tie(max_CFL, max_Fo_mu, max_Fo_k) = block_for_each<CombinedMaxReduction>(mrModelPart.Elements(), [&](Element& rElement){
            const double CFL = FluidCharacteristicNumbersUtilities::CalculateElementCFL(rElement, minimum_h_func, current_dt);
            const auto Fo_numbers = fourier_numbers_function(rElement, minimum_h_func, current_dt);
            return std::make_tuple(CFL, std::get<0>(Fo_numbers), std::get<1>(Fo_numbers));
        });

        // Calculate the new time increment from the maximum characteristic numbers in the mesh
        const double new_dt = CalculateNewDeltaTime(
            current_dt,
            std::make_pair(max_CFL, mCFL),
            std::make_pair(max_Fo_mu, mViscousFourier),
            std::make_pair(max_Fo_k, mThermalFourier));

        return new_dt;

        KRATOS_CATCH("")
    }

    double EstimateDtUtility::EstimateDt() const
    {
        KRATOS_TRY;

        double new_dt = 0.0;

        if (mDtEstimationMagnitudesFlags.Is(CFL_ESTIMATION) && mDtEstimationMagnitudesFlags.IsNot(VISCOUS_FOURIER_ESTIMATION) && mDtEstimationMagnitudesFlags.IsNot(THERMAL_FOURIER_ESTIMATION)) {
            // CFL-based delta time estimation (i.e. incompressible flow)
            new_dt = InternalEstimateDt<true, false, false>();
        } else if (mDtEstimationMagnitudesFlags.Is(CFL_ESTIMATION) && mDtEstimationMagnitudesFlags.IsNot(VISCOUS_FOURIER_ESTIMATION) && mDtEstimationMagnitudesFlags.Is(THERMAL_FOURIER_ESTIMATION)) {
            // CFL and thermal Fourier delta time estimation (i.e. convection-diffusion problems)
            new_dt = InternalEstimateDt<true, false, true>();
        } else if (mDtEstimationMagnitudesFlags.Is(CFL_ESTIMATION) && mDtEstimationMagnitudesFlags.Is(VISCOUS_FOURIER_ESTIMATION) && mDtEstimationMagnitudesFlags.Is(THERMAL_FOURIER_ESTIMATION)) {
            // CFL and both Fourier numbers delta time estimation (i.e. compressible flow)
            new_dt = InternalEstimateDt<true, true, true>();
        } else {
            KRATOS_ERROR << "This option is not supporte yet." << std::endl;
        }

        return new_dt;

        KRATOS_CATCH("")
    }

    void EstimateDtUtility::LimitNewDeltaTime(double& rNewDeltaTime) const
    {
        if (rNewDeltaTime > mDtMax) {
            rNewDeltaTime = mDtMax;
        } else if (rNewDeltaTime < mDtMin) {
            rNewDeltaTime = mDtMin;
        }
    }

} // namespace Kratos.