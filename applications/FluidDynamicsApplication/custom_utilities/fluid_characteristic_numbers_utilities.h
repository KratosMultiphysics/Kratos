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
#include <utility>

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

/**
 * @brief Fluid characteristic numbers calculation utility
 * This class provides static methods to calculate the common adimensional
 * magnitudes that characterize any fluid flow.
 */
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

    /// Intentionally deleting default constructor
    FluidCharacteristicNumbersUtilities() = delete;

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
     * @brief Calculate the element Prandtl number
     * For the given element, this method calculates the Prandtl number
     * @tparam ConsiderArtificialMagnitudes Template parameter specifying if the artificial values are considered
     * @param rElement Element to calculate the Prandtl number
     * @return double The element Prandtl number
     */
    template<bool ConsiderArtificialMagnitudes>
    static double CalculateElementPrandtlNumber(const Element& rElement);

    /**
     * @brief Calculate the elemental Peclet numbers
     * For the given element, this method calculates the Peclet number considering both the
     * dynamic viscosity and the shear conductivity
     * @tparam ConsiderArtificialMagnitudes Template parameter specifying if the artificial values are considered
     * @tparam DensityIsNodal Template parameter specifying if the density is nodally stored
     * @param rElement Element to calculate the Peclet number
     * @return std::pair<double,double> Pair containing the viscosity (first position) and conductivity (second position) Peclet numbers
     */
    template<bool ConsiderArtificialMagnitudes, bool DensityIsNodal>
    static std::pair<double,double> CalculateElementPecletNumbers(
        const Element& rElement,
        const ElementSizeFunctionType& rElementSizeCalculator);

    /**
     * @brief Calculate the elemental Peclet number
     * For the given element, this method calculates the Peclet number
     * @tparam ConsiderArtificialMagnitudes Template parameter specifying if the artificial values are considered
     * @tparam DensityIsNodal Template parameter specifying if the density is nodally stored
     * @param rElement Element to calculate the Peclet number
     * @return double The element Peclet number
     */
    template<bool ConsiderArtificialMagnitudes, bool DensityIsNodal>
    static double CalculateElementViscousPecletNumber(
        const Element& rElement,
        const ElementSizeFunctionType& rElementSizeCalculator);

    /**
     * @brief Calculate the elemental Peclet number
     * For the given element, this method calculates the Peclet number
     * @tparam ConsiderArtificialMagnitudes Template parameter specifying if the artificial values are considered
     * @tparam DensityIsNodal Template parameter specifying if the density is nodally stored
     * @param rElement Element to calculate the Peclet number
     * @return double The element Peclet number
     */
    template<bool ConsiderArtificialMagnitudes, bool DensityIsNodal>
    static double CalculateElementThermalPecletNumber(
        const Element& rElement,
        const ElementSizeFunctionType& rElementSizeCalculator);

    /**
     * @brief Calculate the elemental Fourier numbers
     * For the given element, this method calculates the Fourier number considering both the
     * dynamic viscosity and the shear conductivity
     * @tparam ConsiderArtificialMagnitudes Template parameter specifying if the artificial values are considered
     * @tparam DensityIsNodal Template parameter specifying if the density is nodally stored
     * @param rElement Element to calculate the Fourier number
     * @return std::pair<double,double> Pair containing the viscosity (first position) and thermal (second position) Fourier numbers
     */
    template<bool ConsiderArtificialMagnitudes, bool DensityIsNodal>
    static std::pair<double,double> CalculateElementFourierNumbers(
        const Element& rElement,
        const ElementSizeFunctionType& rElementSizeCalculator,
        const double Dt);

    /**
     * @brief Calculate the elemental Fourier number
     * For the given element, this method calculates the viscous Fourier number
     * @tparam ConsiderArtificialMagnitudes Template parameter specifying if the artificial values are considered
     * @tparam DensityIsNodal Template parameter specifying if the density is nodally stored
     * @param rElement Element to calculate the Fourier number
     * @return double The element Fourier number
     */
    template<bool ConsiderArtificialMagnitudes, bool DensityIsNodal>
    static double CalculateElementViscousFourierNumber(
        const Element& rElement,
        const ElementSizeFunctionType& rElementSizeCalculator,
        const double Dt);

    /**
     * @brief Calculate the elemental Fourier number
     * For the given element, this method calculates the thermal Fourier number
     * @tparam ConsiderArtificialMagnitudes Template parameter specifying if the artificial values are considered
     * @tparam DensityIsNodal Template parameter specifying if the density is nodally stored
     * @param rElement Element to calculate the Fourier number
     * @return double The element Fourier number
     */
    template<bool ConsiderArtificialMagnitudes, bool DensityIsNodal>
    static double CalculateElementThermalFourierNumber(
        const Element& rElement,
        const ElementSizeFunctionType& rElementSizeCalculator,
        const double Dt);

    /**
     * @brief Calculate the element Mach number
     * For the given element, this method calculates the midpoint Mach number
     * Note that it requires the velocity to be stored in the nodal historical
     * database and the sound velocity to be stored in the non-historical one
     * @param rElement Element to calculate the Mach number
     * @return double The element Mach number
     */
    static double CalculateElementMachNumber(const Element& rElement);

    /**
     * @brief Get the minimum element size calculation function
     * This method checks the geometry of the provided element and returns the corresponding
     * minimum element size calculation fucntion.
     * @param rGeometry Geoemtry in which the element size is to be computed
     * @return ElementSizeFunctionType Function to calculate the minimum element size
     */
    static ElementSizeFunctionType GetMinimumElementSizeFunction(const Geometry<Node<3>>& rGeometry);

    /**
     * @brief Get the average element size calculation function
     * This method checks the geometry of the provided element and returns the corresponding
     * average element size calculation fucntion.
     * @param rGeometry Geoemtry in which the element size is to be computed
     * @return ElementSizeFunctionType Function to calculate the average element size
     */
    static ElementSizeFunctionType GetAverageElementSizeFunction(const Geometry<Node<3>>& rGeometry);

    ///@}

private:

    ///@}
    ///@name Operations
    ///@{

    template<bool IsNodal>
    static double AuxiliaryGetDensity(const Element& rElement);

    template<bool AddArtificialValues>
    static std::pair<double,double> GetDiffusivityValues(const Element& rElement);

    template<bool AddArtificialValues>
    static double GetDynamicViscosityValue(const Element& rElement);

    template<bool AddArtificialValues>
    static double GetConductivityValue(const Element& rElement);

    ///@}
};

///@} // Kratos classes

///@} // FluidDynamicsApplication group

} // namespace Kratos.


#endif	/* KRATOS_FLUID_CHARACTERISTIC_NUMBERS_UTILITIES_H */
