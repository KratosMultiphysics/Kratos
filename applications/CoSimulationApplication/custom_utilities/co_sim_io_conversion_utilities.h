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
//

#if !defined(KRATOS_CO_SIM_IO_CONVERSION_UTILITIES_H_INCLUDED )
#define  KRATOS_CO_SIM_IO_CONVERSION_UTILITIES_H_INCLUDED

// System includes

// External includes
#include "custom_external_libraries/CoSimIO/co_sim_io/includes/model_part.hpp"
#include "custom_external_libraries/CoSimIO/co_sim_io/includes/info.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
///@addtogroup CoSimulationApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(CO_SIMULATION_APPLICATION) CoSimIOConversionUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CoSimIOConversionUtilities
    KRATOS_CLASS_POINTER_DEFINITION(CoSimIOConversionUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CoSimIOConversionUtilities() = delete;

    /// Assignment operator.
    CoSimIOConversionUtilities& operator=(CoSimIOConversionUtilities const& rOther) = delete;

    /// Copy constructor.
    CoSimIOConversionUtilities(CoSimIOConversionUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    static void CoSimIOModelPartToKratosModelPart(
        const CoSimIO::ModelPart& rCoSimIOModelPart,
        Kratos::ModelPart& rKratosModelPart,
        const DataCommunicator& rDataComm);

    static void KratosModelPartToCoSimIOModelPart(
        const Kratos::ModelPart& rKratosModelPart,
        CoSimIO::ModelPart& rCoSimIOModelPart);

    static CoSimIO::Info InfoFromParameters(const Parameters rSettings);

    ///@}

}; // Class CoSimIOConversionUtilities

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CO_SIM_IO_CONVERSION_UTILITIES_H_INCLUDED defined
