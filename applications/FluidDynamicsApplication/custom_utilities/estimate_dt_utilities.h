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

#ifndef KRATOS_ESTIMATE_DT_UTILITIES_H
#define	KRATOS_ESTIMATE_DT_UTILITIES_H

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/serializer.h"

// Application includes


namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// Estimate the time step in a fluid problem to obtain a given Courant number.
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) EstimateDtUtility
{
public:

	///@name Type Definitions
	///@{

	/// Pointer definition of EstimateDtUtility
	KRATOS_CLASS_POINTER_DEFINITION(EstimateDtUtility);

    /// Function type for the element size calculator function
    typedef std::function<double(const Geometry<Node<3>>&)> ElementSizeFunctionType;

	///@}
	///@name Life Cycle
	///@{

    /// Constructor for CFD-based time step estimation
    /**
     * @param ModelPart The model part containing the problem mesh
     * @param CFL The user-defined Courant-Friedrichs-Lewy number
     * @param DtMin user-defined minimum time increment allowed
     * @param DtMax user-defined maximum time increment allowed
     */
    EstimateDtUtility(
        ModelPart &ModelPart,
        const double CFL,
        const double DtMin,
        const double DtMax)
        : mrModelPart(ModelPart)
    {
        mCFL = CFL;
        mDtMin = DtMin;
        mDtMax = DtMax;
        mViscousFourier = 0.0;
        mThermalFourier = 0.0;
        mConsiderArtificialDiffusion = false;

        SetDtEstimationMagnitudesFlag();
    }

    /// Complete constructor
    /**
     * @param ModelPart The model part containing the problem mesh
     * @param CFL The user-defined Courant-Friedrichs-Lewy number
     * @param ViscousFourier The user-defined viscosity Peclet number
     * @param ThermalFourier The user-defined thermal conductivity Peclet number
     * @param DtMin user-defined minimum time increment allowed
     * @param DtMax user-defined maximum time increment allowed
     */
    EstimateDtUtility(
        ModelPart &ModelPart,
        const double CFL,
        const double ViscousFourier,
        const double ThermalFourier,
        const bool ConsiderArtificialDiffusion,
        const bool NodalDensityFormulation,
        const double DtMin,
        const double DtMax)
        : mrModelPart(ModelPart)
    {
        mCFL = CFL;
        mDtMin = DtMin;
        mDtMax = DtMax;
        mViscousFourier = ViscousFourier;
        mThermalFourier = ThermalFourier;
        mConsiderArtificialDiffusion = ConsiderArtificialDiffusion;
        mNodalDensityFormulation = NodalDensityFormulation;

        SetDtEstimationMagnitudesFlag();
    }

    /// Constructor with Kratos parameters
    /**
     * @param ModelPart The model part containing the problem mesh
     * @param rParameters Kratos parameters containing the CFL number and max time step
     */
    EstimateDtUtility(
        ModelPart& ModelPart,
        Parameters& rParameters)
        : mrModelPart(ModelPart)
    {
        Parameters defaultParameters(R"({
            "automatic_time_step"           : true,
            "CFL_number"                    : 1.0,
            "Viscous_Fourier_number"        : 0.0,
            "Thermal_Fourier_number"        : 0.0,
            "consider_artificial_diffusion" : false,
            "nodal_density_formulation"     : false,
            "minimum_delta_time"            : 1e-4,
            "maximum_delta_time"            : 0.1
        })");

        rParameters.ValidateAndAssignDefaults(defaultParameters);

        mCFL = rParameters["CFL_number"].GetDouble();
        mViscousFourier = rParameters["Viscous_Fourier_number"].GetDouble();
        mThermalFourier = rParameters["Thermal_Fourier_number"].GetDouble();
        mConsiderArtificialDiffusion = rParameters["consider_artificial_diffusion"].GetBool();
        mNodalDensityFormulation = rParameters["nodal_density_formulation"].GetBool();
        mDtMin = rParameters["minimum_delta_time"].GetDouble();
        mDtMax = rParameters["maximum_delta_time"].GetDouble();

        SetDtEstimationMagnitudesFlag();
    }

    /// Destructor
    ~EstimateDtUtility()
    {}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Set the maximum CFL value allowed
     * This method allows setting the maximum user-defined CFL number
     * @param CFL Tue user-defined maximum CFL number
     */
    void SetCFL(const double CFL);

    /**
     * @brief Set the maximum viscosity Peclet value allowed
     * This method allows setting the maximum user-defined viscosity Peclet number
     * @param ViscousFourier Tue user-defined maximum viscosity Peclet number
     */
    void SetViscousFourier(const double ViscousFourier);

    /**
     * @brief Set the maximum conductivity Peclet value allowed
     * This method allows setting the maximum user-defined thermal conductivity Peclet number
     * @param ThermalFourier Tue user-defined maximum conductivity Peclet number
     */
    void SetThermalFourier(const double ThermalFourier);

    /**
     * @brief Set the minimum time step value allowed
     * This method allows setting the minimum user-defined time increment value
     * @param DtMin The user-defined minimum delta time
     */
    void SetDtMin(const double DtMin);
    
    /**
     * @brief Set the maximum time step value allowed
     * This method allows setting the maximum user-defined time increment value
     * @param DtMax The user-defined maximum delta time
     */
    void SetDtMax(const double DtMax);

    /**
     * @brief Calculate the maximum time step that satisfies the CFL user settings
     * This method calculates the maximum time step that satisfies the Courant-Friedrichs-Lewy (CFL) condition
     * according to the user-defined parameters (CFL and maximum/minimum delta time)
     * @return double A time step value that satisfies the user CFL condition for the current mesh and velocit field
     */
    double EstimateDt() const;

    ///@} // Operators

private:

    ///@name Type definitions
    ///@{

    /// Local flags to determine the magnitudes for the Dt estimation
    KRATOS_DEFINE_LOCAL_FLAG(CFL_ESTIMATION);
    KRATOS_DEFINE_LOCAL_FLAG(VISCOUS_FOURIER_ESTIMATION);
    KRATOS_DEFINE_LOCAL_FLAG(THERMAL_FOURIER_ESTIMATION);

	///@}
    ///@name Member Variables
    ///@{

    double    mCFL;                         // User-defined CFL number
    double    mViscousFourier;              // User-defined viscous Fourier number 
    double    mThermalFourier;              // User-defined thermal Fourier number
    bool      mConsiderArtificialDiffusion; // Speficies if the artificial diffusion values are considered in the Peclet numbers
    bool      mNodalDensityFormulation;     // Specifies if the density is nodally stored (only required for the Peclet number)
    double    mDtMax;                       // User-defined maximum time increment allowed
    double    mDtMin;                       // User-defined minimum time increment allowed
    Flags     mDtEstimationMagnitudesFlags; // Flags indicating the reference magnitudes used in the Dt estimation
    ModelPart &mrModelPart;                 // The problem's model part

	///@}
	///@name Serialization
	///@{


    ///@}
    ///@name Private Operations
    ///@{

    void SetDtEstimationMagnitudesFlag();

    template<const bool ConsiderCFL, const bool ConsiderViscousFourier, const bool ConsiderThermalFourier>
    double InternalEstimateDt() const;

    /**
     * @brief Calculate the new delta time
     * For the provided set of pairs (obtained characteristic number and expected one) this method returns
     * the minimum time increment that fulfils all of them. Note that the minimum delta time is set if the
     * obtained characteristic number is close to zero to avoid the division by zero. Even though this is
     * not the optimal value, it is the safer one.
     * @tparam CharacteristicNumbersPairsType Variadic template argument to specify the obtained and sought characteristic numbers pairs
     * @param CurrentDeltaTime Current delta time
     * @param rCharacteristicNumbersPairs Pairs containing the obtained characteristic number (1st position) and the sought one (2nd position)
     * @return double The minimum delta time among all the provided pairs
     */
    template <class... CharacteristicNumbersPairsType>
    double CalculateNewDeltaTime(
        const double CurrentDeltaTime,
        const CharacteristicNumbersPairsType&... rCharacteristicNumbersPairs) const
    {
        KRATOS_TRY

        // Calculate the corresponding new time increments from the provided pairs
        const double zero_tol = 1.0e-10;
        double new_dt_list[sizeof...(CharacteristicNumbersPairsType)] = {(
            (std::get<0>(rCharacteristicNumbersPairs) > zero_tol) ? CurrentDeltaTime * std::get<1>(rCharacteristicNumbersPairs) / std::get<0>(rCharacteristicNumbersPairs) : mDtMin
        )...};

        // Get the minimum one among all the obtained ones and check the user-defined bounds
        double new_dt = *(std::min_element(new_dt_list, new_dt_list + sizeof...(CharacteristicNumbersPairsType)));
        LimitNewDeltaTime(new_dt);

        // Perform MPI sync if needed
        new_dt = mrModelPart.GetCommunicator().GetDataCommunicator().MinAll(new_dt);

        return new_dt;

        KRATOS_CATCH("");
    }

    /**
     * @brief Limit the new delta time value
     * This method checks if the provided time increment is within the user-defined minimum
     * and maximum bounds. If not, it corrects the provided value accordingly.
     * @param rNewDeltaTime Time increment to be checked
     */
    void LimitNewDeltaTime(double& rNewDeltaTime) const;

    ///@} // Private Operations
};

///@} // Kratos classes

///@} // FluidDynamicsApplication group

} // namespace Kratos.


#endif	/* KRATOS_ESTIMATE_DT_UTILITIES_H */
