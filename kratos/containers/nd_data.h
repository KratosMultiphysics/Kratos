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
#include <atomic>
#include <string>

// External	includes
#include <span/span.hpp>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class NDData
 * @ingroup KratosCore
 * @author Suneth Warnakulasuriya
 * @brief An array with dynamic number of dimensions to interface Kratos C++ objects with numpy easily.
 * @details This class provides possibilities to interface with numpy arrays easily within Kratos C++ interfaces.
 *          This also ensures that, the numpy arrays created from these @ref NDData are sharing
 *          the internal data. So the life time of the internal data will be shared as well, making it possible to
 *          use numpy array pointing to the internal data of an instance of @ref NDData even if
 *          the @ref NDData gets destroyed.
 * @tparam TDataType The type of the data stored in the tensor adaptor.
 */
template<class TDataType>
class KRATOS_API(KRATOS_CORE) NDData {
public:
    ///@name Class definitions
    ///@{

    /**
     * @class PointerWrapper
     * @ingroup KratosCore
     * @author Suneth Warnakulasuriya
     * @brief A class wrapping a pointer to an array.
     * @details This class wraps around a pointer to an array with an @p intrusive_ptr.
     *          It needs to be an @p intrusive_ptr to enable sharing the internal data
     *          represented by @p mpData with numpy, whilst having the possibility to extend
     *          the life time.
     */
    class PointerWrapper
    {
    public:
        ///@name Type definitions
        ///@{

        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(PointerWrapper);

        ///@}
        ///@name Life cycle
        ///@{

        /**
         * @brief Construct a new Pointer Wrapper for a given @p pData pointer.
         * @details This construct a PointerWrapper for the given @p pData pointer.
         *              - @p if IsManaged is true, then the data will be deallocated when the last instance of the
         *                @p intrusive_ptr is destroyed.
         *              - @p if IsManaged is false, then the data will not be managed by this @ref PointerWrapper.
         *
         * @param pData         Pointer to the data
         * @param IsManaged     Whether the data pointed by @p pData is managed by this class or not.
         */
        PointerWrapper(
            TDataType * pData,
            const bool IsManaged)
        : mIsManaged(IsManaged),
          mpData(pData) {}

        ~PointerWrapper() { if (mIsManaged && mpData) delete[] mpData; }

        TDataType * Data() { return mpData; }

        TDataType const * Data() const { return mpData; }

        ///@}

    private:
        ///@name Private member variables
        ///@{

        const bool mIsManaged;

        TDataType * mpData = nullptr;

        ///@}
        ///@name Private operations
        ///@{

        //*********************************************
        //this block is needed for refcounting
        mutable std::atomic<int> mReferenceCounter{0};

        friend void intrusive_ptr_add_ref(const PointerWrapper* x)
        {
            x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
        }

        friend void intrusive_ptr_release(const PointerWrapper* x)
        {
            if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
            std::atomic_thread_fence(std::memory_order_acquire);
            delete x;
            }
        }
        //*********************************************

        ///@}
    };

    ///@}
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(NDData);

    using IndexType = std::size_t;

    using DataType = TDataType;

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Construct a new instance with a provided shape.
     * @details This constructor will only allocate memory for the given shape.
     * @warning The values will be uninitialized.
     *
     * @param rShape Dynamic Dimensional Array shape.
     */
    NDData(
        const DenseVector<unsigned int>& rShape);

    /**
     * @brief Construct a new instance with a given @p rShape, and initialized with the @p Value.
     *
     * @param rShape Array shape.
     * @param Value  Initialization value.
     */
    NDData(
        const DenseVector<unsigned int>& rShape,
        const TDataType Value);

    /**
     * @brief Construct a new instance with given pointer and shape.
     * @details This constructor will construct a dynamic dimensional array with the @p pData pointer
     *          and the given shape.
     *              - If @p Copy = true, then the data in the @p pData will be copied
     *              - If @p Copy = false, then the internal data will point to the @p pData .
     *
     * @warning This may SEGFAULT if a @ref NDData is created with @p Copy = false, and when
     *          any of the @ref ViewData method is called while the @p pData is deallocated.
     *
     * @param pData     Pointer to the data
     * @param rShape    Dynamic Dimensional Array shape.
     * @param Copy      Whether to copy the data @p pData referring to or not.
     */
    NDData(
        TDataType * pData,
        const DenseVector<unsigned int>& rShape,
        const bool Copy = true);

    /**
     * @brief Construct a new instance with @p rShape, copying the array referenced by @p pData.
     *
     * @param pData     Pointer to the data
     * @param rShape    Dynamic Dimensional Array shape.
     */
    NDData(
        TDataType const * pData,
        const DenseVector<unsigned int>& rShape);

    /**
     * @brief Copy constructor.
     * @details Copy constructs having the internal data also copied.
     */
    NDData(
        const NDData& rOther);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Return a view of the internal data structure.
     */
    Kratos::span<const TDataType> ViewData() const;

    /**
     * @brief Return a view of the internal data structure.
     * @throws std::runtime_error If the internal data is already moved.
     */
    Kratos::span<TDataType> ViewData();

    /**
     * @brief Get the underlying pointer which is wrapped around the internal data.
     */
    typename PointerWrapper::Pointer pData() const;

    /**
     * @brief Returns the shape of the tensor.
     * @details This function provides the dimensions of the tensor, where each element in the returned
     *          array corresponds to the size of the tensor in that particular dimension. The first
     *          dimension always represents the number of entities stored in the container.
     *
     * @return A vector containing the size of each dimension of the tensor.
     */
    DenseVector<unsigned int> Shape() const;

    /**
     * @brief Returns the number of elements in the tensor adaptor.
     */
    unsigned int Size() const;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const;

    ///@}

private:
    ///@name private member variables
    ///@{

    const DenseVector<unsigned int> mShape;

    typename PointerWrapper::Pointer mpData;

    ///@}
};

///@}

template<class TDataType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const NDData<TDataType>& rThis)
{
    return rOStream << rThis.Info();
}
} // namespace Kratos