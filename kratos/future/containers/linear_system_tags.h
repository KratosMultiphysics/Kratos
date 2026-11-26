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

#pragma once

// System includes
#include <array>
#include <string_view>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

class LinearSystemTags final
{
public:

    ///@name Type Definitions
    ///@{

    /// Enum defining the available linear system matrices
    enum class SparseMatrixTag
    {
        LHS = 0,             // Left hand side matrix tag
        MassMatrix = 1,      // Mass matrix tag
        StiffnessMatrix = 2, // Stiffness matrix tag
        DampingMatrix = 3,   // Damping matrix tag
        NumberOfTags = 4     // Sentinel tag with the size
    };

    /// Enum defining the available linear system vectors
    enum class DenseVectorTag
    {
        RHS = 0,         // Right hand side vector tag
        Dx = 1,          // Solution increment vector tag
        Eigvals = 2,     // Eigenvalues vector tag
        NumberOfTags = 3 // Sentinel tag with the size
    };

    /// Enum defining the available linear system dense matrices
    enum class DenseMatrixTag
    {
        RHS = 0,         // Right hand side matrix tag
        Dx = 1,          // Solution increment matrix tag
        Eigvects = 2,    // Eigenvectors matrix tag
        NumberOfTags = 3 // Sentinel tag with the size
    };

    ///@}
    ///@name Life Cycle
    ///@{

    LinearSystemTags() = delete; // Delete default constructor since this class is not meant to be instantiated

    ///@}
    ///@name Static Operations
    ///@{

    /**
     * @brief Converts a string to a SparseMatrixTag enum value
     * @param rTagName The string representation of the tag
     * @return SparseMatrixTag The corresponding enum value
     */
    static SparseMatrixTag SparseMatrixTagFromString(const std::string &rTagName)
    {
        constexpr std::array<std::pair<std::string_view, SparseMatrixTag>, static_cast<std::size_t>(SparseMatrixTag::NumberOfTags)> string_to_tag_array = {{
            {"LHS", SparseMatrixTag::LHS},
            {"MassMatrix", SparseMatrixTag::MassMatrix},
            {"StiffnessMatrix", SparseMatrixTag::StiffnessMatrix},
            {"DampingMatrix", SparseMatrixTag::DampingMatrix}
        }};

        for (const auto &[name, tag] : string_to_tag_array) {
            if (name == rTagName) {
                return tag;
            }
        }

        KRATOS_ERROR << "Invalid SparseMatrixTag name: " << rTagName << std::endl;
    }

    /**
     * @brief Converts a string to a DenseVectorTag enum value
     * @param rTagName The string representation of the tag
     * @return DenseVectorTag The corresponding enum value
     */
    static DenseVectorTag DenseVectorTagFromString(const std::string &rTagName)
    {
        constexpr std::array<std::pair<std::string_view, DenseVectorTag>, static_cast<std::size_t>(DenseVectorTag::NumberOfTags)> string_to_tag_array = {{
            {"RHS", DenseVectorTag::RHS},
            {"Dx", DenseVectorTag::Dx},
            {"Eigvals", DenseVectorTag::Eigvals}
        }};

        for (const auto &[name, tag] : string_to_tag_array) {
            if (name == rTagName) {
                return tag;
            }
        }

        KRATOS_ERROR << "Invalid DenseVectorTag name: " << rTagName << std::endl;
    }

    /**
     * @brief Converts a string to a DenseMatrixTag enum value
     * @param rTagName The string representation of the tag
     * @return DenseMatrixTag The corresponding enum value
     */
    static DenseMatrixTag DenseMatrixTagFromString(const std::string &rTagName)
    {
        constexpr std::array<std::pair<std::string_view, DenseMatrixTag>, static_cast<std::size_t>(DenseMatrixTag::NumberOfTags)> string_to_tag_array = {{
            {"RHS", DenseMatrixTag::RHS},
            {"Dx", DenseMatrixTag::Dx},
            {"Eigvects", DenseMatrixTag::Eigvects}
        }};

        for (const auto &[name, tag] : string_to_tag_array) {
            if (name == rTagName) {
                return tag;
            }
        }

        KRATOS_ERROR << "Invalid DenseMatrixTag name: " << rTagName << std::endl;
    }

    /**
     * @brief Get the index of a tag
     * @param Tag The tag of the linear operator to be retrieved (@see SparseMatrixTag)
     * @return constexpr std::size_t The index of the tag
     */
    static constexpr std::size_t GetTagIndex(SparseMatrixTag Tag)
    {
        return static_cast<std::size_t>(Tag);
    }

    /**
     * @brief Get the Tag Index object
     * @param Tag The tag of the linear operator to be retrieved (@see DenseVectorTag)
     * @return constexpr std::size_t The index of the tag
     */
    static constexpr std::size_t GetTagIndex(DenseVectorTag Tag)
    {
        return static_cast<std::size_t>(Tag);
    }

    /**
     * @brief Get the Tag Index object
     * @param Tag The tag of the linear operator to be retrieved (@see DenseMatrixTag)
     * @return constexpr std::size_t The index of the tag
     */
    static constexpr std::size_t GetTagIndex(DenseMatrixTag Tag)
    {
        return static_cast<std::size_t>(Tag);
    }

    ///@}
}; // Class LinearSystemTags

///@} addtogroup block

} // namespace Kratos