//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#ifndef KRATOS_FLUID_CHARACTERISTIC_NUMBERS_UTILITIES_H
#define	KRATOS_FLUID_CHARACTERISTIC_NUMBERS_UTILITIES_H

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

// Application includes


namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidCharacteristicNumbersUtilities
{
public:

	///@name Type Definitions
	///@{

	/// Pointer definition of FluidCharacteristicNumbersUtilities
	KRATOS_CLASS_POINTER_DEFINITION(FluidCharacteristicNumbersUtilities);

    /// Function type for the element size calculator function
    typedef std::function<double(const Geometry<Node<3>>&)> ElementSizeFunctionType;

	///@}
	///@name Life Cycle
	///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Calculate each element's CFL for the current time step and provided model part
     * The elemental CFL is stored in the CFL_NUMBER elemental variable
     * To visualize in the post-process file, remember to print CFL_NUMBER as a Gauss point result
     */
    static void CalculateLocalCFL(ModelPart& rModelPart);

    /**
     * @brief Calulate element CFL number
     * For the given element, this method calculates the CFL number
     * @param rElement Element to calculate the CFL number
     * @param rGeometryInfo Auxiliary geometry data container
     * @param Dt Current time increment
     * @return double The element CFL number
     */
    static double CalculateElementCFL(
        const Element &rElement,
        const ElementSizeFunctionType& rElementSizeCalculator,
        const double Dt);

    /**
     * @brief Get the minimum element size calculation function
     * This method checks the geometry of the provided element and returns the corresponding
     * minimum element size calculation fucntion.
     * @param rGeometry Geoemtry in which the element size is to be computed
     * @return ElementSizeFunctionType Function to calculate the minimum element size
     */
    static ElementSizeFunctionType GetMinimumElementSizeFunction(const Geometry<Node<3>>& rGeometry);
    
    ///@} // Private Operations
};

///@} // Kratos classes

///@} // FluidDynamicsApplication group

} // namespace Kratos.


#endif	/* KRATOS_FLUID_CHARACTERISTIC_NUMBERS_UTILITIES_H */
