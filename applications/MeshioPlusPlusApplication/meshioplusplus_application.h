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

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "includes/kratos_application.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class KratosMeshioPlusPlusApplication
 * @ingroup MeshioPlusPlusApplication
 * @brief The Kratos application wrapping the meshio++ library.
 * @details Provides multi-format mesh IO (see @ref MeshioPlusPlusIO) plus the meshio++
 * mesh and data operations (clean, transform, split, refine, partition, quality, ...)
 * exposed as Kratos utilities, modelers and processes.
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION) KratosMeshioPlusPlusApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosMeshioPlusPlusApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosMeshioPlusPlusApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosMeshioPlusPlusApplication();

    /// Copy constructor.
    KratosMeshioPlusPlusApplication(KratosMeshioPlusPlusApplication const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    KratosMeshioPlusPlusApplication& operator=(KratosMeshioPlusPlusApplication const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosMeshioPlusPlusApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

    ///@}
}; // Class KratosMeshioPlusPlusApplication

///@}

} // namespace Kratos.
