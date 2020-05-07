//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_MODEL_UTILITIES)
#define KRATOS_MODEL_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

// Some forward declaration to avoid many includes
class Model;
class ModelPart;
class Parameters;

/**
 * @namespace ModelUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for Model operations
 * @author Vicente Mataix Ferrandiz
 */
namespace ModelUtilities
{
    /**
     * @brief Auxiliar method to obtain the corresponding modelpart from a given settings and Model
     * @param rModel The model where the where the simulation is performed
     * @param rParameters The parameters of configuration
     */
    ModelPart& KRATOS_API(KRATOS_CORE) GetModelPartFromModelAndSettings(
        Model& rModel,
        const Parameters& rParameters
        );

}; // namespace ModelUtilities
}  // namespace Kratos
#endif /* KRATOS_MODEL_UTILITIES defined */
