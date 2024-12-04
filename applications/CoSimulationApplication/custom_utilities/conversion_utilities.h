//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//                   Ashish Darekar
//

#if !defined(KRATOS_COSIM_CONVERSION_UTILITIES_H_INCLUDED )
#define  KRATOS_COSIM_CONVERSION_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup CoSimulationApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(CO_SIMULATION_APPLICATION) ConversionUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConversionUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ConversionUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConversionUtilities() = delete;

    /// Assignment operator.
    ConversionUtilities& operator=(ConversionUtilities const& rOther) = delete;

    /// Copy constructor.
    ConversionUtilities(ConversionUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Converts elemental data to nodal data by transpose operator. Nodal value is a sum of elemental values, where each elemental value is divided by the number of nodes. Hence the integrated sum is equal on nodes and elements.
     *
     * This utility function takes elemental data from a specified variable and 
     * converts it to nodal data, storing the result in another specified variable.
     *
     * @tparam TDataType The type of the data to be converted.
     * @param rModelPart The model part containing the elements and nodes.
     * @param rElementalVariable The variable containing the elemental data to be converted.
     * @param rNodalVariable The variable where the converted nodal data will be stored.
     */
    template<class TDataType>
    static void ConvertElementalDataToNodalDataTranspose(
        ModelPart& rModelPart,
        const Variable<TDataType>& rElementalVariable,
        const Variable<TDataType>& rNodalVariable );


    /**
     * @brief Converts elemental data to nodal data directly. Nodal value is a sum of elemental values divided by the number of elements.
     * 
     * This function takes elemental data from a specified variable and converts it to nodal data,
     * storing the result in another specified variable. The conversion is performed directly
     * without any intermediate steps.
     * 
     * @tparam TDataType The type of the data to be converted.
     * @param rModelPart The model part containing the elements and nodes.
     * @param rElementalVariable The variable containing the elemental data to be converted.
     * @param rNodalVariable The variable where the converted nodal data will be stored.
     */
    template<class TDataType>
    static void ConvertElementalDataToNodalDataDirect(
        ModelPart& rModelPart,
        const Variable<TDataType>& rElementalVariable,
        const Variable<TDataType>& rNodalVariable );


    /**
     * @brief Converts nodal data to elemental data. Element value is a sum of nodal values divided by the number of nodes.
     *
     * This utility function transfers data from nodal variables to elemental variables
     * within a given ModelPart. It is templated to work with various data types.
     *
     * @tparam TDataType The type of data to be converted.
     * @param rModelPart The ModelPart containing the nodes and elements.
     * @param rElementalVariable The variable to store the converted elemental data.
     * @param rNodalVariable The variable containing the nodal data to be converted.
     */
    template<class TDataType>
    static void ConvertNodalDataToElementalDataDirect(
        ModelPart& rModelPart,
        const Variable<TDataType>& rElementalVariable,
        const Variable<TDataType>& rNodalVariable );

        /**
     * @brief Converts nodal data to elemental data using transpose operator. Element value is a sum of nodal value contribution, which is split equally to all neighbouring elements.
     *
     * This utility function transfers data from nodal variables to elemental variables
     * within a given ModelPart. It is templated to work with various data types.
     *
     * @tparam TDataType The type of data to be converted.
     * @param rModelPart The ModelPart containing the nodes and elements.
     * @param rElementalVariable The variable to store the converted elemental data.
     * @param rNodalVariable The variable containing the nodal data to be converted.
     */
    template<class TDataType>
    static void ConvertNodalDataToElementalDataTranspose(
        ModelPart& rModelPart,
        const Variable<TDataType>& rElementalVariable,
        const Variable<TDataType>& rNodalVariable );


    ///@}

}; // Class ConversionUtilities

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_COSIM_CONVERSION_UTILITIES_H_INCLUDED defined