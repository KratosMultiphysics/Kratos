//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) TabularDatabaseIO
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(TabularDatabaseIO);

    ///@}
    ///@name Life cycle
    ///@{

    virtual ~TabularDatabaseIO() = default;

    ///@}
    ///@name Input / output
    ///@{

    /**
     * @brief Get the Row Ids which are read.
     *
     * @return std::vector<int>  List of row ids.
     */
    virtual std::vector<int> GetRowIds() const
    {
        KRATOS_ERROR << "Please implement this method in the derived class.";
        return std::vector<int>{};
    }

    /**
     * @brief Reads the value on the row @p RowId and column given by @p rKey.
     */
    virtual void Read(
        bool& rValue,
        const int RowId,
        const std::string& rKey)
    {
        KRATOS_ERROR << "Please implement this method in the derived class.";
    }

    /**
     * @brief Reads the value on the row @p RowId and column given by @p rKey.
     */
    virtual void Read(
        int& rValue,
        const int RowId,
        const std::string& rKey)
    {
        KRATOS_ERROR << "Please implement this method in the derived class.";
    }

    /**
     * @brief Reads the value on the row @p RowId and column given by @p rKey.
     */
    virtual void Read(
        double& rValue,
        const int RowId,
        const std::string& rKey)
    {
        KRATOS_ERROR << "Please implement this method in the derived class.";
    }

    /**
     * @brief Reads the value on the row @p RowId and column given by @p rKey.
     */
    virtual void Read(
        std::string& rValue,
        const int RowId,
        const std::string& rKey)
    {
        KRATOS_ERROR << "Please implement this method in the derived class.";
    }

    /**
     * @brief Write the value to the row @p RowId and column given by @p rKey.
     */
    virtual void Write(
        const bool Value,
        const int RowId,
        const std::string& rKey)
    {
        KRATOS_ERROR << "Please implement this method in the derived class.";
    }

    /**
     * @brief Write the value to the row @p RowId and column given by @p rKey.
     */
    virtual void Write(
        const int Value,
        const int RowId,
        const std::string& rKey)
    {
        KRATOS_ERROR << "Please implement this method in the derived class.";
    }

    /**
     * @brief Write the value to the row @p RowId and column given by @p rKey.
     */
    virtual void Write(
        const double Value,
        const int RowId,
        const std::string& rKey)
    {
        KRATOS_ERROR << "Please implement this method in the derived class.";
    }

    /**
     * @brief Write the value to the row @p RowId and column given by @p rKey.
     */
    virtual void Write(
        const std::string& rValue,
        const int RowId,
        const std::string& rKey)
    {
        KRATOS_ERROR << "Please implement this method in the derived class.";
    }

    /// Turn back information as a string.
    virtual std::string Info() const { return "TabularDatabaseIO"; }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const { rOStream << Info(); }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
};

///@}
///@name Input and output
///@{

inline std::ostream& operator << (std::ostream& rOStream, const TabularDatabaseIO& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos