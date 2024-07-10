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

    template<class TDataType>
    static void ConvertElementalDataToNodalData(
        ModelPart& rModelPart,
        const Variable<TDataType>& rElementalVariable,
        const Variable<TDataType>& rNodalVariable );


    ///@}

}; // Class ConversionUtilities

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_COSIM_CONVERSION_UTILITIES_H_INCLUDED defined
