//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/serializer.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/**
 * @class FileSerializer
 * @ingroup KratosCore
 * @brief This class provides a simplified interface for serializing data to a file.
 * @details Note that you may not override any load or save method of the Serializer base class, as they are not virtual.
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) FileSerializer 
    : public Serializer 
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Condition
    KRATOS_CLASS_POINTER_DEFINITION(FileSerializer);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor that initializes the FileSerializer.
     * @param Filename The name of the file for serialization.
     * @param rTrace Type of serialization trace to be employed (default: SERIALIZER_NO_TRACE).
     */
    FileSerializer(std::string const& Filename, Serializer::TraceType const& rTrace = SERIALIZER_NO_TRACE);

    /**
     * @brief Destructor
     */
    ~FileSerializer() override {}

    ///@}
private:
    ///@name Private Operators
    ///@{

    /**
     * @brief Deleted assignment operator to prevent unwanted copying.
     * @param rOther Another instance of FileSerializer to assign from.
     * @return FileSerializer& Reference to the assigned FileSerializer object.
     */
    FileSerializer& operator=(FileSerializer const& rOther) = delete;

    ///@}
    ///@name Private Life Cycle
    ///@{

    /**
     * @brief Deleted copy constructor to prevent unwanted copying.
     * @param rOther Another instance of FileSerializer to construct from.
     */
    FileSerializer(FileSerializer const& rOther) = delete;

    ///@}
};
///@}
}
