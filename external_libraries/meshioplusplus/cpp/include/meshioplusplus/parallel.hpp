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
 * @file parallel.hpp
 * @brief `parallel_for`/`parallel_for_bw`: a backend-agnostic parallel loop
 * over a compile-time-selected SEQ/STL/OpenMP/TBB implementation.
 *
 * The active backend is chosen at compile time by the `MESHIOPLUSPLUS_PARALLEL_*`
 * preprocessor definitions (set from CMake's `MESHIOPLUSPLUS_PARALLEL_BACKEND` =
 * `AUTO|SEQ|STL|OPENMP|TBB`; `AUTO` prefers OpenMP — portable across
 * manylinux/MSVC/macOS without needing TBB — then falls back to STL(+TBB) if
 * detected, else SEQ). `parallel_backend_name()`/`_core.__parallel_backend__`
 * report which one is active. Iterations passed to `parallel_for` must be
 * independent (no cross-iteration state) since they may run concurrently in
 * any order; the first exception thrown by any iteration is captured and
 * rethrown once the parallel region has joined (via `detail::FirstException`),
 * so callers see ordinary C++ exception semantics rather than `std::terminate`
 * or a lost exception.
 *
 * There are two flavors, distinguished by how many threads they are allowed
 * to use:
 *  - `parallel_for` — uses all available cores (up to `max_threads` if
 *    non-zero). Appropriate for compute-bound loops where per-element work
 *    is real computation, e.g. zlib/base64 encode-decode in
 *    `detail/vtu_binary.hpp` and ASCII value formatting.
 *  - `parallel_for_bw` — caps the thread count to `parallel_bandwidth_threads`
 *    (4). Appropriate for memory-bandwidth-bound loops — byte-swap,
 *    transpose, index gather — which saturate a socket's memory bandwidth
 *    with only a few threads and then *regress* as thread count grows
 *    further (more cache contention and dispatch overhead without more
 *    usable bandwidth), unlike compute-bound loops which keep scaling to all
 *    cores.
 *
 * To add a new backend (e.g. Kokkos, HPX): add one CMake branch that defines
 * a new `MESHIOPLUSPLUS_PARALLEL_<NAME>` macro and links the dependency, then
 * add one `#elif defined(MESHIOPLUSPLUS_PARALLEL_<NAME>)` branch in
 * `detail::parallel_for_impl` below (and extend `parallel_backend_name()`
 * to report it).
 */

// System includes
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <utility>

#if defined(MESHIOPLUSPLUS_PARALLEL_STL)
#include <execution>
#include <thread>
#include <vector>
#endif

// External includes
#if defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
#include <omp.h>
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#endif

namespace meshioplusplus {

/**
 * @brief Default grain size (minimum iterations per dispatched chunk) for
 * `parallel_for`/`parallel_for_bw` when the caller doesn't override it.
 *
 * Below this many total iterations, `parallel_for` runs sequentially rather
 * than paying parallel dispatch overhead (see the `n <= grain` check in
 * `parallel_for` below). Callers with atypically coarse or fine per-iteration
 * work (e.g. one whole zlib block per iteration) pass an explicit smaller
 * `grain` (often `1`) so each iteration dispatches individually.
 */
inline constexpr std::size_t parallel_grain_default = 2048;

/**
 * @brief Thread cap used by `parallel_for_bw` for memory-bandwidth-bound loops.
 *
 * Memory-bandwidth-bound loops (byte-swap, transpose, gather) saturate a
 * socket's bandwidth with only a few threads and then *regress* as thread
 * overhead and cache contention grow — unlike compute-bound loops (zlib,
 * base64) which scale to all cores. Cap the bandwidth-bound loops here.
 */
inline constexpr unsigned parallel_bandwidth_threads = 4;

/**
 * @brief Name of the parallel backend selected at compile time.
 *
 * Reflects whichever of `MESHIOPLUSPLUS_PARALLEL_STL`/`_OPENMP`/`_TBB` was
 * defined (by CMake, based on `MESHIOPLUSPLUS_PARALLEL_BACKEND`); none of
 * them defined means the sequential fallback. Exposed to Python as
 * `_core.__parallel_backend__` so tests/diagnostics can assert which backend
 * actually built.
 * @return One of `"stl"`, `"openmp"`, `"tbb"`, `"seq"`.
 */
constexpr const char* parallel_backend_name() {
#if defined(MESHIOPLUSPLUS_PARALLEL_STL)
    return "stl";
#elif defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
    return "openmp";
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
    return "tbb";
#else
    return "seq";
#endif
}

namespace detail {

/**
 * @brief Captures the first exception thrown by any parallel iteration, to
 * be rethrown by the caller after the parallel region joins.
 *
 * Iterations run on multiple threads cannot let a C++ exception escape
 * across the parallelism boundary (OpenMP/TBB would `std::terminate`), so
 * each backend wraps its per-iteration body in `Run()`, which catches
 * everything and records only the *first* exception (subsequent ones from
 * other threads are discarded — `mRaised` is a one-shot latch via
 * `std::atomic_flag`). After the parallel region has fully joined, the
 * caller calls `RethrowIfAny()` to surface that exception on the calling
 * thread with normal C++ semantics.
 */
class FirstException {
public:
    template <class Body>
    void Run(Body&& body) noexcept {
        try {
            body();
        } catch (...) {
            if (!mRaised.test_and_set(std::memory_order_acq_rel))
                mEptr = std::current_exception();
        }
    }
    void RethrowIfAny() {
        if (mEptr)
            std::rethrow_exception(mEptr);
    }

private:
    std::atomic_flag mRaised = ATOMIC_FLAG_INIT;
    std::exception_ptr mEptr;
};

/**
 * @brief Backend-specific dispatch of `n` independent iterations of `f`.
 *
 * Exactly one `#if`/`#elif` branch compiles, selected by the
 * `MESHIOPLUSPLUS_PARALLEL_*` macro CMake defined:
 *  - **STL**: splits `[0, n)` into up to `hardware_concurrency() * 4` chunks
 *    (fewer if `grain`/`max_threads` constrain it further) and runs them via
 *    `std::for_each(std::execution::par, ...)` over a small chunk table
 *    (iterated explicitly because PSTL algorithms require
 *    `Cpp17ForwardIterator`s, which `iota_view` iterators don't satisfy on
 *    every implementation).
 *  - **OpenMP**: `#pragma omp parallel for schedule(dynamic, chunk)` with
 *    `chunk = max(grain/4, 1)`. Dynamic (not static) scheduling matters on
 *    hybrid P+E-core CPUs, where a static split would leave slow E-cores as
 *    stragglers while fast P-cores idle at the join; `grain/4` keeps
 *    dispatch overhead negligible for fine-grained loops while still
 *    honouring explicitly coarse callers (e.g. VTU zlib blocks pass
 *    `grain=1` because each iteration is already a whole compress, so
 *    per-iteration dispatch is exactly what's wanted — the chunk size must
 *    never be floored above the caller's `grain`).
 *  - **TBB**: `tbb::parallel_for` over a `blocked_range` of grain size
 *    `grain`, optionally under a `tbb::global_control` limiting
 *    `max_allowed_parallelism` to `max_threads`.
 *  - **(none, SEQ)**: a plain sequential loop; `grain`/`max_threads` are
 *    unused (cast to `void` to silence warnings).
 *
 * Every branch funnels per-iteration exceptions through a `FirstException`
 * so exactly one is rethrown after the region joins.
 *
 * @tparam F Callable invoked as `f(std::size_t i)` for each `i` in `[0, n)`.
 * @param n Number of iterations.
 * @param rF The per-iteration body (iterations must be independent).
 * @param grain Minimum unit of work per dispatched chunk/task.
 * @param max_threads Cap on threads used (0 = no cap, use all available).
 */
template <class F>
void parallel_for_impl(std::size_t n, F& rF, std::size_t grain, unsigned max_threads) {
#if defined(MESHIOPLUSPLUS_PARALLEL_STL)
    struct Chunk {
        std::size_t mBegin, mEnd;
    };
    const std::size_t hw = std::max<std::size_t>(1, std::thread::hardware_concurrency());
    std::size_t max_chunks = hw * 4;
    if (max_threads)
        max_chunks = std::min<std::size_t>(max_chunks, max_threads);
    const std::size_t by_grain = (n + grain - 1) / grain;
    const std::size_t nchunks = std::max<std::size_t>(1, std::min(max_chunks, by_grain));
    const std::size_t per = (n + nchunks - 1) / nchunks;
    // PSTL algorithms require Cpp17ForwardIterators (iota_view iterators do
    // not qualify on all implementations), so iterate a small chunk table.
    std::vector<Chunk> chunks;
    chunks.reserve(nchunks);
    for (std::size_t b = 0; b < n; b += per)
        chunks.push_back({b, std::min(b + per, n)});
    FirstException exc;
    std::for_each(std::execution::par, chunks.begin(), chunks.end(), [&](const Chunk& c) {
        exc.Run([&] {
            for (std::size_t i = c.mBegin; i < c.mEnd; ++i)
                rF(i);
        });
    });
    exc.RethrowIfAny();
#elif defined(MESHIOPLUSPLUS_PARALLEL_OPENMP)
    FirstException exc;
    const long long nn = static_cast<long long>(n);
    const int nt = max_threads ? std::min<int>(static_cast<int>(max_threads), omp_get_max_threads())
                               : omp_get_max_threads();
    // Dynamic scheduling: on hybrid CPUs (P + E cores) a static split makes the
    // slow cores stragglers while the fast ones idle at the join; moderately
    // sized dynamic chunks self-balance with negligible dispatch overhead.
    // grain/4 keeps dispatch rare for fine-grained loops while honouring
    // explicitly coarse loops (e.g. the VTU zlib blocks pass grain=1: each
    // iteration is a whole compress, so per-iteration dispatch is ideal).
    const long long chunk = static_cast<long long>(std::max<std::size_t>(grain / 4, 1));
#pragma omp parallel for schedule(dynamic, chunk) num_threads(nt)
    for (long long i = 0; i < nn; ++i) {
        exc.Run([&] { rF(static_cast<std::size_t>(i)); });
    }
    exc.RethrowIfAny();
#elif defined(MESHIOPLUSPLUS_PARALLEL_TBB)
    FirstException exc;
    auto body = [&] {
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, n, grain),
                          [&](const tbb::blocked_range<std::size_t>& r) {
                              exc.Run([&] {
                                  for (std::size_t i = r.begin(); i != r.end(); ++i)
                                      rF(i);
                              });
                          });
    };
    if (max_threads) {
        tbb::global_control gc(tbb::global_control::max_allowed_parallelism, max_threads);
        body();
    } else {
        body();
    }
    exc.RethrowIfAny();
#else  // MESHIOPLUSPLUS_PARALLEL_SEQ (and the safe default)
    (void)grain;
    (void)max_threads;
    for (std::size_t i = 0; i < n; ++i)
        rF(i);
#endif
}

}  // namespace detail

/**
 * @brief Runs `n` independent iterations of `f(i)`, in parallel when it's
 * worthwhile, using the compile-time-selected backend (see
 * `parallel_backend_name()`).
 *
 * If `n <= grain`, runs sequentially in-line — the fixed cost of dispatching
 * a parallel region isn't worth it for small workloads. Otherwise delegates
 * to `detail::parallel_for_impl`. `f` must be safe to invoke concurrently
 * from multiple threads for different `i` (no shared mutable state without
 * external synchronization); the first exception any invocation throws is
 * captured and rethrown on the calling thread after all iterations
 * complete (partial results/side effects from other iterations are not
 * rolled back).
 *
 * @tparam F Callable invoked as `f(std::size_t i)`.
 * @param n Number of iterations; a no-op if `n == 0`.
 * @param f The per-iteration body.
 * @param grain Minimum number of iterations to bother parallelizing, and
 *              (backend-dependent) the target chunk size once it does;
 *              defaults to `parallel_grain_default` (2048). Pass a small
 *              value (e.g. `1`) when each iteration is already coarse work
 *              (a whole zlib block, a whole compress) so dispatch happens
 *              per-iteration rather than being batched further.
 * @param max_threads Cap on threads used; `0` (the default) means "use all
 *                     available". Pass `parallel_bandwidth_threads`
 *                     (or call `parallel_for_bw` instead) for
 *                     memory-bandwidth-bound loops.
 */
template <class F>
void parallel_for(std::size_t n, F&& f, std::size_t grain = parallel_grain_default,
                  unsigned max_threads = 0) {
    if (n == 0)
        return;
    if (n <= grain) {
        for (std::size_t i = 0; i < n; ++i)
            f(i);
        return;
    }
    detail::parallel_for_impl(n, f, grain, max_threads);
}

/**
 * @brief `parallel_for`, thread-capped for memory-bandwidth-bound loops.
 *
 * Convenience wrapper that forwards to `parallel_for` with
 * `max_threads = parallel_bandwidth_threads` (4). Use this for byte-swap,
 * transpose, and index-gather loops: they saturate a socket's memory
 * bandwidth with only a few threads and then *regress* — more threads add
 * cache contention and dispatch overhead without more usable bandwidth —
 * unlike genuinely compute-bound loops (zlib/base64), which should use
 * plain `parallel_for` to scale across all cores.
 *
 * @tparam F Callable invoked as `f(std::size_t i)`.
 * @param n Number of iterations; a no-op if `n == 0`.
 * @param f The per-iteration body.
 * @param grain Minimum iterations per chunk; see `parallel_for`'s `grain`.
 */
template <class F>
void parallel_for_bw(std::size_t n, F&& f, std::size_t grain = parallel_grain_default) {
    parallel_for(n, std::forward<F>(f), grain, parallel_bandwidth_threads);
}

}  // namespace meshioplusplus
