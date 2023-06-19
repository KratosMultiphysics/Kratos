//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SINGLE_IMPORT_MODEL_PART_H)
#define KRATOS_SINGLE_IMPORT_MODEL_PART_H

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_parameters.h"

namespace Kratos 
{

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) SingleImportModelPart
{
public:
    /**
     * @brief This method imports one single ModelPart
     * @param rModelPart The model part where to read the input
     * @param ModelPartImportParameters The model part import parameters
     * @param InputType The input type
     */
    static void Import(
        ModelPart& rModelPart, 
        Parameters ModelPartImportParameters,
        const std::string& InputType
        );
};

///@}

} // namespace Kratos.

#endif // KRATOS_SINGLE_IMPORT_MODEL_PART_H  defined
