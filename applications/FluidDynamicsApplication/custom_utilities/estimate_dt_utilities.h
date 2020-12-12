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
        mPecletViscosity = 0.0;
        mPecletConductivity = 0.0;
        mConsiderArtificialDiffusion = false;

        SetDtEstimationMagnitudesFlag();
    }

    /// Complete constructor
    /**
     * @param ModelPart The model part containing the problem mesh
     * @param CFL The user-defined Courant-Friedrichs-Lewy number
     * @param PecletViscosity The user-defined viscosity Peclet number
     * @param PecletConductivity The user-defined thermal conductivity Peclet number
     * @param DtMin user-defined minimum time increment allowed
     * @param DtMax user-defined maximum time increment allowed
     */
    EstimateDtUtility(
        ModelPart &ModelPart,
        const double CFL,
        const double PecletViscosity,
        const double PecletConductivity,
        const bool ConsiderArtificialDiffusion,
        const bool NodalDensityFormulation,
        const double DtMin,
        const double DtMax)
        : mrModelPart(ModelPart)
    {
        mCFL = CFL;
        mDtMin = DtMin;
        mDtMax = DtMax;
        mPecletViscosity = PecletViscosity;
        mPecletConductivity = PecletConductivity;
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
            "Peclet_number_viscosity"       : 0.0,
            "Peclet_number_conductivity"    : 0.0,
            "consider_artificial_diffusion" : false,
            "nodal_density_formulation"     : false,
            "minimum_delta_time"            : 1e-4,
            "maximum_delta_time"            : 0.1
        })");

        rParameters.ValidateAndAssignDefaults(defaultParameters);

        mCFL = rParameters["CFL_number"].GetDouble();
        mPecletViscosity = rParameters["Peclet_number_viscosity"].GetDouble();
        mPecletConductivity = rParameters["Peclet_number_conductivity"].GetDouble();
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
     * @param PecletViscosity Tue user-defined maximum viscosity Peclet number
     */
    void SetPecletViscosity(const double PecletViscosity);

    /**
     * @brief Set the maximum conductivity Peclet value allowed
     * This method allows setting the maximum user-defined thermal conductivity Peclet number
     * @param PecletConductivity Tue user-defined maximum conductivity Peclet number
     */
    void SetPecletConductivity(const double PecletConductivity);

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
    KRATOS_DEFINE_LOCAL_FLAG(PECLET_VISCOSITY_ESTIMATION);
    KRATOS_DEFINE_LOCAL_FLAG(PECLET_CONDUCTIVITY_ESTIMATION);

	///@}
    ///@name Member Variables
    ///@{

    double    mCFL;                         // User-defined CFL number
    double    mPecletViscosity;             // User-defined viscosity Peclet number 
    double    mPecletConductivity;          // User-defined conductivity Peclet number
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

    template<const bool ConsiderCFL, const bool ConsiderPecletViscosity, const bool ConsiderPecletConductivity, const bool NodalDensityFormulation = false>
    double InternalEstimateDt() const;

    ///@} // Private Operations
};

///@} // Kratos classes

///@} // FluidDynamicsApplication group

} // namespace Kratos.


#endif	/* KRATOS_ESTIMATE_DT_UTILITIES_H */
