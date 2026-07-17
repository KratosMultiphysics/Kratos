//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file ndarray.hpp
 * @brief `NDArray`: a minimal typed, n-dimensional, contiguous (row-major)
 * array — the storage primitive of `meshioplusplus::Mesh`.
 *
 * `NDArray` is used for points, cell connectivity, and every point/cell/field
 * data array. It either *owns* its buffer (the common case: data produced by
 * a reader) or is a non-owning *view* over externally-owned memory (used to
 * wrap a numpy buffer zero-copy on the write path — see `py_to_mesh` in
 * `bindings/np_conversions.hpp`). The binding layer converts between
 * `NDArray` and numpy at the I/O boundary: owning buffers are moved into a
 * capsule backing a writeable numpy array on read, and numpy buffers are
 * wrapped as views (no copy) on write. `Dtype()` records the element type
 * with an internal `DType` enum rather than a template parameter, so
 * `NDArray` can be stored uniformly (e.g. in `Mesh::mCellData`) regardless of
 * the numpy dtype it came from; `As<T>()` reinterprets the raw buffer as `T`
 * once the caller has determined (or asserted) the appropriate type.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

namespace meshioplusplus {

namespace detail {
/**
 * @brief Allocator that leaves elements *default*-initialized rather than
 * value-initialized.
 *
 * For a trivial type like `std::byte` that means the buffer is left
 * uninitialized instead of zero-filled. `NDArray` uses this (via `ByteBuf`)
 * so a buffer it is about to fully overwrite (reader outputs, reconstruction
 * blocks — see `NDArray::Uninit`) can skip the zero-fill `memset`, which for
 * a fresh large allocation is an entire extra cold pass over just-faulted
 * pages (numpy's `calloc`-backed arrays skip it too, for the same reason).
 * `std::vector` with this allocator stays copyable/movable like a normal
 * vector, unlike a raw `unique_ptr` buffer, so `NDArray` can keep value
 * semantics.
 *
 * @tparam T The element type being allocated (used as `std::byte` here).
 *
 * @note The member names below (`value_type`, `allocate`, `deallocate`,
 * `construct`, `rebind`, `operator==`/`operator!=`) are fixed by the C++
 * standard library's Allocator named requirements and must keep these exact
 * spellings regardless of naming convention.
 */
template <class T>
struct NoInitAllocator {
    using value_type = T;
    NoInitAllocator() = default;
    template <class U>
    NoInitAllocator(const NoInitAllocator<U>&) noexcept {}
    template <class U>
    struct rebind {
        using other = NoInitAllocator<U>;
    };
    T* allocate(std::size_t n) { return std::allocator<T>{}.allocate(n); }
    void deallocate(T* pPtr, std::size_t n) { std::allocator<T>{}.deallocate(pPtr, n); }
    // Default-init (no zeroing) for the no-arg case resize() uses; forward
    // everything else so the vector still behaves normally.
    template <class U>
    void construct(U* pPtr) noexcept(std::is_nothrow_default_constructible_v<U>) {
        ::new (static_cast<void*>(pPtr)) U;
    }
    template <class U, class... Args>
    void construct(U* pPtr, Args&&... args) {
        ::new (static_cast<void*>(pPtr)) U(std::forward<Args>(args)...);
    }
    template <class U>
    bool operator==(const NoInitAllocator<U>&) const noexcept {
        return true;
    }
    template <class U>
    bool operator!=(const NoInitAllocator<U>&) const noexcept {
        return false;
    }
};
}  // namespace detail

/**
 * @brief Scalar element type of an `NDArray`, mirroring the numpy dtypes the
 * binding layer converts to/from.
 */
enum class DType {
    Float32,
    Float64,
    Int8,
    Int16,
    Int32,
    Int64,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
};

/**
 * @brief Size in bytes of one element of the given dtype.
 * @param dt The dtype to query.
 * @return 1, 2, 4, or 8, matching the C++ scalar type `dt` represents.
 */
inline std::size_t dtype_size(DType dt) {
    switch (dt) {
        case DType::Float32:
            return 4;
        case DType::Float64:
            return 8;
        case DType::Int8:
        case DType::UInt8:
            return 1;
        case DType::Int16:
        case DType::UInt16:
            return 2;
        case DType::Int32:
        case DType::UInt32:
            return 4;
        case DType::Int64:
        case DType::UInt64:
            return 8;
    }
    return 0;
}

/**
 * @brief numpy dtype string (kind + itemsize) for a `DType`, e.g. `"f8"`, `"i4"`.
 * @param dt The dtype to convert.
 * @return A numpy-style struct format code understood by `numpy.dtype(...)`.
 */
inline const char* dtype_numpy_str(DType dt) {
    switch (dt) {
        case DType::Float32:
            return "f4";
        case DType::Float64:
            return "f8";
        case DType::Int8:
            return "i1";
        case DType::Int16:
            return "i2";
        case DType::Int32:
            return "i4";
        case DType::Int64:
            return "i8";
        case DType::UInt8:
            return "u1";
        case DType::UInt16:
            return "u2";
        case DType::UInt32:
            return "u4";
        case DType::UInt64:
            return "u8";
    }
    return "f8";
}

/**
 * @brief A minimal typed, n-dimensional, row-major contiguous array.
 *
 * `NDArray` is either *owning* (holds its own `ByteBuf`, freed on
 * destruction) or a non-owning *view* over externally-managed memory
 * (`mView != nullptr`); `IsView()` distinguishes the two, and `Data()`
 * transparently returns whichever buffer is active. Views exist so the
 * write path can wrap a numpy array's memory directly (see
 * `bindings/np_conversions.hpp`'s `py_to_mesh`) without copying it into a
 * C++-owned buffer; `MakeOwned()` is the escape hatch for turning a view
 * into an owning copy when a buffer must outlive the memory it points to.
 * There is no reference counting: a view's caller is responsible for
 * keeping the underlying memory alive for the `NDArray`'s lifetime.
 */
class NDArray {
public:
    NDArray() = default;

    /**
     * @brief Constructs an owning array with a zero-initialized buffer.
     * @param dt Element dtype.
     * @param shape Row-major dimensions; total element count is their product.
     */
    NDArray(DType dt, std::vector<std::size_t> shape) : mDtype(dt), mShape(std::move(shape)) {
        const std::size_t nb = Nbytes();
        mOwned.resize(nb);                  // uninitialised (NoInitAllocator)
        std::memset(mOwned.data(), 0, nb);  // explicit zero-fill
    }

    /**
     * @brief Constructs an owning array whose buffer is left *uninitialized*.
     *
     * Only safe for callers that immediately overwrite every byte — typical
     * uses are reader outputs (the whole buffer is about to be filled from
     * the parsed file) and cell-block reconstruction (e.g.
     * `detail::reconstruct_cells` in `vtk_cells.hpp`). Skips both the extra
     * allocator zero-fill and, more importantly, the cold first-touch page
     * faults a `memset` would otherwise incur on a fresh large allocation —
     * the same optimization numpy applies to its own `calloc`-avoidance path.
     * Prefer the two-argument constructor whenever the buffer might not be
     * fully overwritten.
     *
     * @param dt Element dtype.
     * @param shape Row-major dimensions; total element count is their product.
     * @return A new owning, uninitialized `NDArray`.
     */
    static NDArray Uninit(DType dt, std::vector<std::size_t> shape) {
        NDArray a;
        a.mDtype = dt;
        a.mShape = std::move(shape);
        a.mOwned.resize(a.Nbytes());  // no memset
        return a;
    }

    /**
     * @brief Constructs a non-owning view over externally-owned row-major memory.
     *
     * Used to wrap a numpy array's buffer directly at the write boundary
     * (zero-copy): the C++ writer reads through `pPtr` but never frees it.
     * @param dt Element dtype of the memory at `pPtr`.
     * @param shape Row-major dimensions describing how to interpret `pPtr`.
     * @param pPtr Pointer to caller-owned memory; the caller must keep it
     *            alive for at least the lifetime of the returned `NDArray`
     *            (and of any `NDArray` copies/moves derived from it that
     *            remain a view).
     * @return A new non-owning `NDArray` view.
     */
    static NDArray MakeView(DType dt, std::vector<std::size_t> shape, std::byte* pPtr) {
        NDArray a;
        a.mDtype = dt;
        a.mShape = std::move(shape);
        a.mView = pPtr;
        return a;
    }

    DType Dtype() const { return mDtype; }
    const std::vector<std::size_t>& Shape() const { return mShape; }
    std::size_t Ndim() const { return mShape.size(); }
    /** @brief Whether this array is a non-owning view (vs. owning its buffer). */
    bool IsView() const { return mView != nullptr; }

    /** @brief Total element count (product of `Shape()`), or 0 if `Shape()` is empty. */
    std::size_t Size() const {
        if (mShape.empty())
            return 0;
        return std::accumulate(mShape.begin(), mShape.end(), std::size_t{1},
                               std::multiplies<std::size_t>());
    }
    /** @brief Total buffer size in bytes: `Size() * dtype_size(Dtype())`. */
    std::size_t Nbytes() const { return Size() * dtype_size(mDtype); }

    /** @brief Raw pointer to the active buffer (owned or view), for writing. */
    std::byte* Data() { return mView ? mView : mOwned.data(); }
    /** @brief Raw pointer to the active buffer (owned or view), read-only. */
    const std::byte* Data() const { return mView ? mView : mOwned.data(); }

    /**
     * @brief Changes the logical shape in place without touching the buffer.
     *
     * A no-op if the new shape's element count doesn't match the current
     * one (the mismatched reshape is silently ignored rather than throwing).
     * @param new_shape The desired row-major dimensions.
     */
    void Reshape(std::vector<std::size_t> new_shape) {
        std::size_t n = new_shape.empty()
                            ? 0
                            : std::accumulate(new_shape.begin(), new_shape.end(), std::size_t{1},
                                              std::multiplies<std::size_t>());
        if (n != Size())
            return;  // ignore inconsistent reshape
        mShape = std::move(new_shape);
    }

    /**
     * @brief Turns a view into an owning copy in place; a no-op if already owning.
     *
     * Copies the viewed memory into a freshly-allocated owned buffer and
     * clears the view pointer. Used before handing a buffer's lifetime over
     * to Python via a capsule (`mesh_to_py`), where the destination `NDArray`
     * must actually own the memory it hands off.
     */
    void MakeOwned() {
        if (mView == nullptr)
            return;
        const std::size_t nb = Nbytes();
        ByteBuf buf;
        buf.resize(nb);  // uninitialised; fully overwritten by the memcpy below
        std::memcpy(buf.data(), mView, nb);
        mOwned = std::move(buf);
        mView = nullptr;
    }

    /**
     * @brief Reinterprets the raw buffer as a `T*`. No dtype check is performed
     * — the caller must ensure `T` matches `Dtype()`.
     * @tparam T The scalar type to view the buffer as.
     * @return Pointer to the first element, typed as `T`.
     */
    template <typename T>
    T* As() {
        return reinterpret_cast<T*>(Data());
    }
    /** @brief `const` overload of `As()`. */
    template <typename T>
    const T* As() const {
        return reinterpret_cast<const T*>(Data());
    }

private:
    using ByteBuf = std::vector<std::byte, detail::NoInitAllocator<std::byte>>;
    DType mDtype = DType::Float64;
    std::vector<std::size_t> mShape;
    ByteBuf mOwned;
    std::byte* mView = nullptr;
};

}  // namespace meshioplusplus
