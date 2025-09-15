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
#include <tuple>
#include <utility>

// External includes

// Project includes
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief This provides zero derivatives for given any variable
 *
 * This class can be used to represent the derivative computation
 * when interested residualdoes not have dependence on the derivative
 * variable which results in zero derivatives analytically.
 *
 * @tparam TNumNodes    Number of nodes in the element
 * @tparam TBlockSize   Block size of the dofs
 */
template<unsigned int TNumNodes, unsigned int TBlockSize>
class ZeroDerivatives
{
public:
    ///@name Local type definitions
    ///@{

    using IndexType = std::size_t;

    constexpr static IndexType NumNodes = TNumNodes;

    constexpr static IndexType BlockSize = TBlockSize;

    ///@}
    ///@name Public static operations
    ///@{

    template<IndexType TSize, class... TArgs>
    void inline CalculateGaussPointResidualsDerivativeContributions(
        BoundedVector<double, TSize>& rResidualDerivative,
        TArgs&... rArgs) const
    {
        rResidualDerivative.clear();
    }

    ///@}
};

/**
 * @brief Assembles sub-vectors to element vectors and matrices.
 *
 * This class is used to assemble discretized scalar residual vectors having TNumNodes size for a given equation
 * to the local element residuals. The TSubAssemblyType is supposed to give the sub-vectors of
 * the respective scalar residual equation.
 *
 * @tparam TSubAssemblyType         The contributions calculation class type to be used with the sub-assembly.
 * @tparam TElementDataHolderIndex  The index of the data container used in TSubAssemblyType to compute sub-vectors.
 * @tparam TRowStartingIndex        The row starting index within a block for the respective scalar residual.
 * @tparam TColumnStartingIndex     The column starting index within a block for the respective scalar residual.
 */
template<class TSubAssemblyType, unsigned int TElementDataHolderIndex, unsigned int TRowStartingIndex, unsigned int TColumnStartingIndex = 0>
class SubAssembly : public TSubAssemblyType
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using SubAssemblyType = TSubAssemblyType;

    static constexpr unsigned int ElementDataIndex = TElementDataHolderIndex;

    static constexpr unsigned int RowStartingIndex = TRowStartingIndex;

    static constexpr unsigned int ColumnStartingIndex = TColumnStartingIndex;

    static constexpr unsigned int NumNodes = SubAssemblyType::NumNodes;

    static constexpr unsigned int BlockSize = SubAssemblyType::BlockSize;

    static constexpr unsigned int ResidualSize = NumNodes * BlockSize;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Get the Sub Vector
     *
     * This returns the sub-vector which needs to be filled with the scalar equation information
     * which is assembled to the elemental LHS or RHS.
     *
     * @return BoundedVector<double, ResidualSize>& Returns a reference to the sub-vector
     */
    inline BoundedVector<double, ResidualSize>& GetSubVector()
    {
        return mSubVector;
    }

    /**
     * @brief Get the Element Data Container
     *
     * This returns the element data container with ElementDataIndex from a combined data container (std::tuple)
     *
     * @tparam TCombinedElementDataContainer
     * @param rCombinedDataContainer
     * @return std::tuple_element<ElementDataIndex, TCombinedElementDataContainer>::type&
     */
    template<class TCombinedElementDataContainer>
    inline typename std::tuple_element<ElementDataIndex, TCombinedElementDataContainer>::type& GetElementDataContainer(TCombinedElementDataContainer& rCombinedDataContainer) const
    {
        static_assert(
            ElementDataIndex < std::tuple_size_v<TCombinedElementDataContainer>,
            "Required Element data container index is more than the available element data containers.");
        return std::get<ElementDataIndex>(rCombinedDataContainer);
    }

    /**
     * @brief Assembles sub-vector to the element matrix
     *
     * This assembles the sub-vector to the given matrix rOutput. The row to which this is assembled is
     * determined by the NodeIndex and RowStartingIndex and the TAssemblyRowBlockSize. The sub-vector is
     * distributed in that specific row based on the TAssemblyColumnBlockSize, and ColumnStartingIndex.
     *
     * This assembly is done node wise in an outer loop from this method, hence for rows, NodeIndex is passed
     * The sub-vector holds all the contributions for each node, hence no need to pass the node index for the
     * columns
     *
     * @tparam TAssemblyRowBlockSize        Block size of the rows in the final assembled matrix rOutput
     * @tparam TAssemblyColumnBlockSize     Block size of the columns in the final assembled matrix rOutput
     * @param rOutput                       Assembled matrix
     * @param NodeIndex                     Node index corresponds to which row the sub-vector belongs to.
     */
    template<IndexType TAssemblyRowBlockSize, IndexType TAssemblyColumnBlockSize = TAssemblyRowBlockSize>
    inline void AssembleSubVectorToMatrix(
        Matrix& rOutput,
        const IndexType NodeIndex) const
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(rOutput.size1() != TAssemblyRowBlockSize * NumNodes)
            << "rOutput.size1() does not have the required size. [ rOutput.size1() = "
            << rOutput.size1() << ", required_size = TAssemblyRowBlockSize * NumNodes = "
            << TAssemblyRowBlockSize * NumNodes << " ].\n";

        KRATOS_DEBUG_ERROR_IF(rOutput.size2() != TAssemblyColumnBlockSize * NumNodes)
            << "rOutput.size2() does not have the required size. [ rOutput.size2() = "
            << rOutput.size2() << ", required_size = TAssemblyColumnBlockSize * NumNodes = "
            << TAssemblyColumnBlockSize * NumNodes << " ].\n";

        for (IndexType i = 0; i < NumNodes; ++i) {
            for (IndexType j = 0; j < BlockSize; ++j) {
                rOutput(TAssemblyRowBlockSize * NodeIndex + RowStartingIndex, i * TAssemblyColumnBlockSize + j + ColumnStartingIndex) += mSubVector[i * BlockSize + j];
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Assembles sub-vector to a vector
     *
     * This method assembles mSubVector to rOutput vector (which is the element vector).
     *
     * @tparam TAssemblyRowBlockSize    Block size of the rOutput vector
     * @param rOutput                   Assembled vector
     */
    template<IndexType TAssemblyRowBlockSize>
    inline void AssembleSubVectorToVector(Vector& rOutput) const
    {
        KRATOS_TRY

        static_assert(ColumnStartingIndex == 0);

        KRATOS_DEBUG_ERROR_IF(rOutput.size() != TAssemblyRowBlockSize * NumNodes)
            << "rOutput.size() does not have the required size. [ rOutput.size() = "
            << rOutput.size() << ", required_size = TAssemblyRowBlockSize * NumNodes = "
            << TAssemblyRowBlockSize * NumNodes << " ].\n";

        for (IndexType i = 0; i < NumNodes; ++i) {
            for (IndexType j = 0; j < BlockSize; ++j) {
                rOutput[i * TAssemblyRowBlockSize + j] += mSubVector[i * BlockSize + j];
            }
        }

        KRATOS_CATCH("");
    }

    ///@}
private:
    ///@name Private members
    ///@{

    BoundedVector<double, ResidualSize> mSubVector;

    ///@}
};

/**
 * @brief This is traits clas sto hold combined data containers and calculation containers
 *
 * This class holds combined element data containers (std::tuple<T1, T2, ...>) and required calculation
 * containers (std::tuple<G1, G2, G3, ...>).
 *
 * @tparam TCombinedElementDataContainer       std::tuple of data containers
 * @tparam TCombinedCalculationContainers      std::tuple of calculation containers
 */
template <
    class TCombinedElementDataContainer,
    class TCombinedCalculationContainers
>
class CalculationContainerTraits
{
public:
    ///@name Derivative type definitions
    ///@{

    using CombinedElementDataContainerType = TCombinedElementDataContainer;

    using CombinedCalculationContainersType = TCombinedCalculationContainers;

    ///@}
};
}