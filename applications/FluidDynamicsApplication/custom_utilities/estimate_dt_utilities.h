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
template< unsigned int TDim >
class EstimateDtUtility
{
public:

    ///@name Life Cycle
    ///@{

    /// Constructor
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
            "automatic_time_step"   : true,
            "CFL_number"            : 1.0,
            "minimum_delta_time"    : 1e-4,
            "maximum_delta_time"    : 0.1
        })");

        rParameters.ValidateAndAssignDefaults(defaultParameters);

        mCFL = rParameters["CFL_number"].GetDouble();
        mDtMin = rParameters["minimum_delta_time"].GetDouble();
        mDtMax = rParameters["maximum_delta_time"].GetDouble();
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

    /**
     * @brief Calculate each element's CFL for the current time step and provided model part
     * The elemental CFL is stored in the CFL_NUMBER elemental variable
     * To visualize in the post-process file, remember to print CFL_NUMBER as a Gauss point result
     */
    static void CalculateLocalCFL(ModelPart& rModelPart);

    ///@} // Operators

private:

    ///@name Auxiliary Data types
    ///@{

    /**
     * @brief Geometry data container
     * Auxiliary container to avoid allocations within the parallel loops
     */
    struct GeometryDataContainer
    {
        double Area;
        array_1d<double, TDim+1> N;
        BoundedMatrix<double, TDim+1, TDim> DN_DX;
    };

    ///@}
    ///@name Member Variables
    ///@{

    double    mCFL;         // User-defined CFL number
    double    mDtMax;       // User-defined maximum time increment allowed
    double    mDtMin;       // User-defined minimum time increment allowed
    ModelPart &mrModelPart; // The problem's model part

	///@}
	///@name Serialization
	///@{

	friend class Serializer;

	void save(Serializer& rSerializer) const;

	void load(Serializer& rSerializer);

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Calulate element CFL number
     * For the given element, this method calculates the CFL number
     * @param rElement Element to calculate the CFL number
     * @param rGeometryInfo Auxiliary geometry data container
     * @param Dt Current time increment
     * @return double The element CFL number
     */
    static double CalculateElementCFL(
        Element &rElement,
        GeometryDataContainer& rGeometryInfo,
        double Dt);

    ///@} // Private Operations
};

///@} // Kratos classes

///@} // FluidDynamicsApplication group

} // namespace Kratos.


#endif	/* KRATOS_ESTIMATE_DT_UTILITIES_H */
