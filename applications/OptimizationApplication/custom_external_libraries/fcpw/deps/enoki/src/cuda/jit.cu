/*
    src/cuda/jit.cpp -- CUDA backend (Tracing JIT compiler)

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyrighe (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include <cuda.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <array>
#include <algorithm>
#include <sstream>
#include <set>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <chrono>
#include <mutex>
#include "common.cuh"

/// Should the implementation use streams to schedule kernels in parallel if possible?
#define ENOKI_CUDA_USE_STREAMS 1

/// Synchronize with the device after each kernel launch (useful for debugging)
#define ENOKI_CUDA_LAUNCH_BLOCKING 0

#if defined(NDEBUG)
#  define ENOKI_CUDA_DEFAULT_LOG_LEVEL 0
#else
#  define ENOKI_CUDA_DEFAULT_LOG_LEVEL 1
#endif

/// Reserved registers for grid-stride loop indexing
#define ENOKI_CUDA_REG_RESERVED 10

NAMESPACE_BEGIN(enoki)

using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

// Forward declarations
ENOKI_EXPORT void cuda_inc_ref_ext(uint32_t);
ENOKI_EXPORT void cuda_inc_ref_int(uint32_t);
ENOKI_EXPORT void cuda_dec_ref_ext(uint32_t);
ENOKI_EXPORT void cuda_dec_ref_int(uint32_t);
ENOKI_EXPORT size_t cuda_register_size(EnokiType type);
ENOKI_EXPORT uint32_t cuda_trace_append(EnokiType type, const char *cmd, uint32_t arg1);

// -----------------------------------------------------------------------
//! @{ \name 'Variable' type that is used to record instruction traces
// -----------------------------------------------------------------------

struct Variable {
    /// Data type of this variable
    EnokiType type;

    /// PTX instruction to compute it
    std::string cmd;

    /// Associated label (mainly for debugging)
    std::string label;

    /// Number of entries
    size_t size = 0;

    /// Pointer to device memory
    void *data = nullptr;

    /// External (i.e. by Enoki) reference count
    uint32_t ref_count_ext = 0;

    /// Internal (i.e. within the PTX instruction stream) reference count
    uint32_t ref_count_int = 0;

    /// Dependencies of this instruction
    std::array<uint32_t, 3> dep = { 0, 0, 0 };

    /// Extra dependency (which is not directly used in arithmetic, e.g. scatter/gather)
    uint32_t extra_dep = 0;

    /// Does the instruction have side effects (e.g. 'scatter')
    bool side_effect = false;

    /// A variable is 'dirty' if there are pending scatter operations to it
    bool dirty = false;

    /// Free 'data' after this variable is no longer referenced?
    bool free = true;

    /// Optimization: is this a direct pointer (rather than an array which stores a pointer?)
    bool direct_pointer = false;

    /// Size of the (heuristic for instruction scheduling)
    uint32_t subtree_size = 0;

    Variable(EnokiType type) : type(type) { }

    ~Variable() { if (free && data != nullptr) cuda_free(data); }

    bool is_collected() const {
        return ref_count_int == 0 && ref_count_ext == 0;
    }
};

ENOKI_EXPORT void cuda_shutdown();

#if ENOKI_CUDA_USE_STREAMS == 1
struct Stream {
    cudaStream_t stream = nullptr;
    cudaEvent_t event = nullptr;

    void init() {
        cuda_check(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
        cuda_check(cudaEventCreateWithFlags(&event, cudaEventDisableTiming));
    }

    void release() {
        cuda_check(cudaEventDestroy(event));
        cuda_check(cudaStreamDestroy(stream));
    }
};
#endif

enum MallocType { Normal, Managed, Host };

struct TaggedSize {
    MallocType type;
    size_t size;

    TaggedSize(MallocType type, size_t size) : type(type), size(size) { }
    bool operator==(const TaggedSize &ts) const { return type == ts.type && size == ts.size; }
    bool operator!=(const TaggedSize &ts) const { return type != ts.type || size != ts.size; }
};

struct TaggedSizeHasher {
    size_t operator()(const TaggedSize &ts) const {
        return std::hash<size_t>()((ts.size << 2) + (size_t) ts.type);
    }
};

struct Context {
    /// Current variable index
    uint32_t ctr = 0;

    /// Enumerates "live" (externally referenced) variables and statements with side effects
    std::set<uint32_t> live;

    /// Enumerates "dirty" variables (targets of 'scatter' operations that have not yet executed)
    std::vector<uint32_t> dirty;

    /// Stores the mapping from variable indices to variables
    std::unordered_map<uint32_t, Variable> variables;

    /// Stores the mapping from pointer addresses to variable indices
    std::unordered_map<const void *, uint32_t> ptr_map;

    /// Current operand array for scatter/gather
    uint32_t scatter_gather_operand = 0;

    /// Current log level (0 == none, 1 == minimal, 2 == moderate, 3 == max.)
    uint32_t log_level = ENOKI_CUDA_DEFAULT_LOG_LEVEL;

    /// Callback that will be invoked before each cuda_eval() call
    std::vector<std::pair<void(*)(void *), void *>> callbacks;

    /// Include printf function declaration in PTX assembly?
    bool include_printf = false;


    /// Hash table of previously compiled kernels
    std::unordered_map<std::string, std::pair<CUmodule, CUfunction>, StringHasher> kernels;

    #if ENOKI_CUDA_USE_STREAMS == 1
        /// Streams for parallel execution
        std::vector<Stream> streams;

        /// Event on default stream
        cudaEvent_t stream_0_event = nullptr;
    #endif

    /// Default thread and block count for kernels
    uint32_t block_count = 0, thread_count = 0;

    /// Map of unused memory regions
    std::unordered_multimap<TaggedSize, void *, TaggedSizeHasher> free_map;

    /// Map of currently used memory regions
    std::unordered_map<void *, TaggedSize> used_map;

    /// Memory usage watermarks
    size_t used = 0, used_managed = 0, used_host = 0;
    size_t watermark = 0, watermark_managed = 0, watermark_host = 0;

    /// Mutex protecting the malloc-related data structures
    std::recursive_mutex malloc_mutex;

    ~Context() { clear(); }

    Variable &operator[](uint32_t i) {
        auto it = variables.find(i);
        if (it == variables.end())
            throw std::runtime_error("CUDABackend: referenced unknown variable " + std::to_string(i));
        return it->second;
    }

    void clear() {
#if !defined(NDEBUG)
        if (log_level >= 1) {
            if (ctr != 0 || !variables.empty())
                std::cerr << "cuda_shutdown()" << std::endl;
            size_t n_live = 0;
            for (auto const &var : variables) {
                if (var.first < ENOKI_CUDA_REG_RESERVED)
                    continue;
                if (n_live < 10) {
                    std::cerr << "cuda_shutdown(): variable " << var.first << " is still live. "<< std::endl;
                    if (n_live == 9)
                        std::cerr << "(skipping remainder)" << std::endl;
                }
                ++n_live;
            }
            if (n_live > 0)
                std::cerr << "cuda_shutdown(): " << n_live
                          << " variables were still live at shutdown." << std::endl;
        }
#endif
        ctr = 0;
        dirty.clear();
        variables.clear();
        live.clear();
        scatter_gather_operand = 0;
        include_printf = false;

        #if ENOKI_CUDA_USE_STREAMS == 1
            for (Stream &stream : streams)
                stream.release();
            streams.clear();
            if (stream_0_event) {
                cuda_check(cudaEventDestroy(stream_0_event));
                stream_0_event = nullptr;
            }
        #endif

        for (auto &kv : kernels)
            cuda_check(cuModuleUnload(kv.second.first));
        kernels.clear();
        cuda_sync();
        cuda_malloc_trim();
    }

    Variable& append(EnokiType type) {
        return variables.emplace(ctr++, type).first->second;
    }
};

static Context *__context = nullptr;
bool installed_shutdown_handler = false;
ENOKI_EXPORT void cuda_init();

inline static Context &context() {
    if (ENOKI_UNLIKELY(__context == nullptr))
        cuda_init();
    return *__context;
}

ENOKI_EXPORT void cuda_init() {
    // initialize CUDA
    cudaFree(0);

    // We don't really use shared memory, so put more into L1 cache.
    cuda_check(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));

    /// Reserve indices for reserved kernel variables
    if (__context)
        delete __context;
    __context = new Context();
    Context &ctx = *__context;
    ctx.clear();
    ctx.append(EnokiType::Invalid);
    ctx.append(EnokiType::UInt64);
    while (ctx.variables.size() != ENOKI_CUDA_REG_RESERVED)
        ctx.append(EnokiType::UInt32);

    #if ENOKI_CUDA_USE_STREAMS == 1
        ctx.streams.resize(5);
        for (size_t i = 0; i < 5; ++i)
            ctx.streams[i].init();
        cuda_check(cudaEventCreateWithFlags(&ctx.stream_0_event, cudaEventDisableTiming));
    #endif

    ctx.kernels.reserve(1000);

    int device, num_sm;
    cudaGetDevice(&device);
    cudaDeviceGetAttribute(&num_sm, cudaDevAttrMultiProcessorCount, device);

    ctx.block_count = next_power_of_two(num_sm) * 2;
    ctx.thread_count = 128;

    ctx.used = 0;
    ctx.used_managed = 0;
    ctx.used_host = 0;
    ctx.watermark = 0;
    ctx.watermark_managed = 0;
    ctx.watermark_host = 0;

    if (!installed_shutdown_handler) {
        installed_shutdown_handler = true;
        atexit(cuda_shutdown);
    }
}

ENOKI_EXPORT void cuda_shutdown() {
    if (__context) {
        __context->clear();
        delete __context;
        __context = nullptr;
    }
}

ENOKI_EXPORT void *cuda_var_ptr(uint32_t index) {
    return context()[index].data;
}

ENOKI_EXPORT size_t cuda_var_size(uint32_t index) {
    return context()[index].size;
}

ENOKI_EXPORT void cuda_var_set_label(uint32_t index, const char *str) {
    Context &ctx = context();
    ctx[index].label = str;
#if !defined(NDEBUG)
    if (ctx.log_level >= 4)
        std::cerr << "cuda_var_set_label(" << index << "): " << str << std::endl;
#endif
}

ENOKI_EXPORT uint32_t cuda_var_set_size(uint32_t index, size_t size, bool copy) {
    Context &ctx = context();
#if !defined(NDEBUG)
    if (ctx.log_level >= 4)
        std::cerr << "cuda_var_set_size(" << index << "): " << size << std::endl;
#endif

    Variable &var = ctx[index];
    if (var.size == size)
        return index;

    if (var.data != nullptr || var.ref_count_int > 0) {
        if (var.size == 1 && copy) {
            uint32_t index_new =
                cuda_trace_append(var.type, "mov.$t1 $r1, $r2", index);
            ctx[index_new].size = size;
            cuda_dec_ref_ext(index);
            return index_new;
        }

        throw std::runtime_error(
            "cuda_var_set_size(): attempted to resize variable " +
            std::to_string(index) +
            " which was already allocated (current size = " +
            std::to_string(var.size) +
            ", requested size = " + std::to_string(size) + ")");
    }
    var.size = size;
    return index;
}

ENOKI_EXPORT uint32_t cuda_var_register(EnokiType type, size_t size,
                                        void *ptr, bool free) {
    Context &ctx = context();
    uint32_t idx = ctx.ctr;
#if !defined(NDEBUG)
    if (ctx.log_level >= 4)
        std::cerr << "cuda_var_register(" << idx << "): " << ptr
                  << ", size=" << size << ", free=" << free << std::endl;
#endif
    if (size == 0)
        throw std::runtime_error("cuda_var_register(): attempted to create a "
                                 "variable of size zero!");
    Variable &v = ctx.append(type);
    v.data = ptr;
    v.size = size;
    v.free = free;
    cuda_inc_ref_ext(idx);
    return idx;
}

ENOKI_EXPORT uint32_t cuda_var_register_ptr(const void *ptr) {
    Context &ctx = context();
    auto it = ctx.ptr_map.find(ptr);
    if (it != ctx.ptr_map.end()) {
        cuda_inc_ref_ext(it->second);
        return it->second;
    }

    uint32_t idx = ctx.ctr;
#if !defined(NDEBUG)
    if (ctx.log_level >= 4)
        std::cerr << "cuda_var_register_ptr(" << idx << "): " << ptr
                  << std::endl;
#endif
    Variable &v = ctx.append(EnokiType::Pointer);
    v.data = (void *) ptr;
    v.size = 1;
    v.free = false;
    v.direct_pointer = true;
    cuda_inc_ref_ext(idx);
    ctx.ptr_map[ptr] = idx;
    return idx;
}

ENOKI_EXPORT uint32_t cuda_var_copy_to_device(EnokiType type, size_t size,
                                              const void *value) {
    size_t total_size = size * cuda_register_size(type);

    void *tmp        = cuda_host_malloc(total_size),
         *device_ptr = cuda_malloc(total_size);

    memcpy(tmp, value, total_size);
    cuda_check(cudaMemcpyAsync(device_ptr, tmp, total_size,
                               cudaMemcpyHostToDevice));

    cuda_host_free(tmp);
    return cuda_var_register(type, size, device_ptr, true);
}

ENOKI_EXPORT void cuda_var_free(uint32_t idx) {
    Context &ctx = context();
    Variable &v = ctx[idx];
#if !defined(NDEBUG)
    if (ctx.log_level >= 5) {
        std::cerr << "cuda_var_free(" << idx << ") = " << v.data;
        if (!v.free)
            std::cerr << " (not deleted)";
        std::cerr << std::endl;
    }
#endif
    if (v.direct_pointer) {
        auto it = ctx.ptr_map.find(v.data);
        assert(it != ctx.ptr_map.end());
        ctx.ptr_map.erase(it);
    }
    for (int i = 0; i < 3; ++i)
        cuda_dec_ref_int(v.dep[i]);
    cuda_dec_ref_ext(v.extra_dep);
    ctx.variables.erase(idx); // invokes Variable destructor + cudaFree().
}

ENOKI_EXPORT void cuda_make_managed(uint32_t idx) {
    if (idx == 0)
        return;

    Context &ctx = context();
    Variable &v = ctx[idx];
#if !defined(NDEBUG)
    if (ctx.log_level >= 4)
        std::cerr << "cuda_make_managed(" << idx << ") = " << v.data << std::endl;
#endif
    size_t total_size = v.size * cuda_register_size(v.type);
    if (total_size == 0)
        return;

    if (v.data == nullptr || v.dirty)
        cuda_eval();

    cudaPointerAttributes attr;
    cuda_check(cudaPointerGetAttributes(&attr, v.data));
    if (attr.type == cudaMemoryTypeManaged)
        return;

    void *ptr = cuda_managed_malloc(total_size);
    cuda_check(cudaMemcpy(ptr, v.data, total_size,
                          cudaMemcpyDeviceToDevice));
    cuda_free(v.data);
    v.data = ptr;
}

ENOKI_EXPORT void cuda_set_scatter_gather_operand(uint32_t idx, bool gather) {
    Context &ctx = context();
    if (idx != 0) {
        Variable &v = ctx[idx];
        if (v.data == nullptr || (gather && v.dirty))
            cuda_eval();
    }
    ctx.scatter_gather_operand = idx;
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Common functionality to distinguish types
// -----------------------------------------------------------------------

ENOKI_EXPORT size_t cuda_register_size(EnokiType type) {
    switch (type) {
        case EnokiType::UInt8:
        case EnokiType::Int8:
        case EnokiType::Bool:    return 1;
        case EnokiType::UInt16:
        case EnokiType::Int16:   return 2;
        case EnokiType::UInt32:
        case EnokiType::Int32:   return 4;
        case EnokiType::Pointer:
        case EnokiType::UInt64:
        case EnokiType::Int64:   return 8;
        case EnokiType::Float32: return 4;
        case EnokiType::Float64: return 8;
        default: return (size_t) -1;
    }
}

ENOKI_EXPORT const char *cuda_register_type(EnokiType type) {
    switch (type) {
        case EnokiType::UInt8: return "u8";
        case EnokiType::Int8: return "s8";
        case EnokiType::UInt16: return "u16";
        case EnokiType::Int16: return "s16";
        case EnokiType::UInt32: return "u32";
        case EnokiType::Int32: return "s32";
        case EnokiType::Pointer:
        case EnokiType::UInt64: return "u64";
        case EnokiType::Int64: return "s64";
        case EnokiType::Float16: return "f16";
        case EnokiType::Float32: return "f32";
        case EnokiType::Float64: return "f64";
        case EnokiType::Bool: return "pred";
        default: return nullptr;
    }
}

ENOKI_EXPORT const char *cuda_register_type_bin(EnokiType type) {
    switch (type) {
        case EnokiType::UInt8:
        case EnokiType::Int8: return "b8";
        case EnokiType::UInt16:
        case EnokiType::Float16:
        case EnokiType::Int16: return "b16";
        case EnokiType::Float32:
        case EnokiType::UInt32:
        case EnokiType::Int32: return "b32";
        case EnokiType::Pointer:
        case EnokiType::Float64:
        case EnokiType::UInt64:
        case EnokiType::Int64: return "b64";
        case EnokiType::Bool: return "pred";
        default: return nullptr;
    }
}

ENOKI_EXPORT const char *cuda_register_name(EnokiType type) {
    switch (type) {
        case EnokiType::UInt8:
        case EnokiType::Int8:   return "%b";
        case EnokiType::UInt16:
        case EnokiType::Int16:   return "%w";
        case EnokiType::UInt32:
        case EnokiType::Int32:   return "%r";
        case EnokiType::Pointer:
        case EnokiType::UInt64:
        case EnokiType::Int64:   return "%rd";
        case EnokiType::Float32: return "%f";
        case EnokiType::Float64: return "%d";
        case EnokiType::Bool:    return "%p";
        default: return nullptr;
    }
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Reference counting (internal means: dependency within
//           JIT trace, external means: referenced by an Enoki array)
// -----------------------------------------------------------------------

ENOKI_EXPORT void cuda_inc_ref_ext(uint32_t index) {
    if (index < ENOKI_CUDA_REG_RESERVED)
        return;
    Context &ctx = context();
    Variable &v = ctx[index];
    v.ref_count_ext++;

#if !defined(NDEBUG)
    if (ctx.log_level >= 5)
        std::cerr << "cuda_inc_ref_ext(" << index << ") -> "
                  << v.ref_count_ext << std::endl;
#endif
}

ENOKI_EXPORT void cuda_inc_ref_int(uint32_t index) {
    if (index < ENOKI_CUDA_REG_RESERVED)
        return;
    Context &ctx = context();
    Variable &v = ctx[index];
    v.ref_count_int++;

#if !defined(NDEBUG)
    if (ctx.log_level >= 5)
        std::cerr << "cuda_inc_ref_int(" << index << ") -> "
                  << v.ref_count_int << std::endl;
#endif
}

ENOKI_EXPORT void cuda_dec_ref_ext(uint32_t index) {
    Context &ctx = context();
    if (index < ENOKI_CUDA_REG_RESERVED || ctx.variables.empty())
        return;
    Variable &v = ctx[index];

    if (ENOKI_UNLIKELY(v.ref_count_ext == 0)) {
        fprintf(stderr, "cuda_dec_ref_ext(): Node %u has no external references!\n", index);
        exit(EXIT_FAILURE);
    }

#if !defined(NDEBUG)
    if (ctx.log_level >= 5)
        std::cerr << "cuda_dec_ref_ext(" << index << ") -> "
                  << (v.ref_count_ext - 1) << std::endl;
#endif

    v.ref_count_ext--;

    if (v.ref_count_ext == 0 && !v.side_effect)
        ctx.live.erase(index);

    if (v.is_collected())
        cuda_var_free(index);
}

ENOKI_EXPORT void cuda_dec_ref_int(uint32_t index) {
    if (index < ENOKI_CUDA_REG_RESERVED)
        return;
    Context &ctx = context();
    Variable &v = ctx[index];

    if (ENOKI_UNLIKELY(v.ref_count_int == 0)) {
        fprintf(stderr, "cuda_dec_ref_int(): Node %u has no internal references!\n", index);
        exit(EXIT_FAILURE);
    }

#if !defined(NDEBUG)
    if (ctx.log_level >= 5)
        std::cerr << "cuda_dec_ref_int(" << index << ") -> "
                  << (v.ref_count_int - 1) << std::endl;
#endif

    v.ref_count_int--;

    if (v.is_collected())
        cuda_var_free(index);
}

ENOKI_EXPORT void cuda_var_mark_side_effect(uint32_t index) {
    Context &ctx = context();
#if !defined(NDEBUG)
    if (ctx.log_level >= 4)
        std::cerr << "cuda_var_mark_side_effect(" << index << ")" << std::endl;
#endif

    assert(index >= ENOKI_CUDA_REG_RESERVED);
    ctx[index].side_effect = true;
}

ENOKI_EXPORT void cuda_var_mark_dirty(uint32_t index) {
    Context &ctx = context();
#if !defined(NDEBUG)
    if (ctx.log_level >= 4)
        std::cerr << "cuda_var_mark_dirty(" << index << ")" << std::endl;
#endif

    assert(index >= ENOKI_CUDA_REG_RESERVED);
    ctx[index].dirty = true;
    ctx.dirty.push_back(index);
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name JIT trace append routines
// -----------------------------------------------------------------------

static void strip_ftz(Variable &v) {
    if (v.type != EnokiType::Float32) {
        size_t offset = v.cmd.find(".ftz");
        if (offset != std::string::npos)
            v.cmd.replace(offset, 4, "");
    }
}

ENOKI_EXPORT uint32_t cuda_trace_append(EnokiType type,
                                        const char *cmd) {
    Context &ctx = context();
    uint32_t idx = ctx.ctr;
#if !defined(NDEBUG)
    if (ctx.log_level >= 4)
        std::cerr << "cuda_trace_append(" << idx << "): " << cmd << std::endl;
#endif
    Variable &v = ctx.append(type);
    v.cmd = cmd;
    v.subtree_size = 1;
    v.size = 1;
    cuda_inc_ref_ext(idx);
    ctx.live.insert(idx);
    strip_ftz(v);
    return idx;
}

ENOKI_EXPORT uint32_t cuda_trace_append(EnokiType type,
                                        const char *cmd,
                                        uint32_t arg1) {
    if (ENOKI_UNLIKELY(arg1 == 0))
        throw std::runtime_error(
            "cuda_trace_append(): arithmetic involving "
                        "uninitialized variable!");
    Context &ctx = context();
    const Variable &v1 = ctx[arg1];

    if (v1.dirty)
        cuda_eval();

    uint32_t idx = ctx.ctr;

#if !defined(NDEBUG)
    if (ctx.log_level >= 4)
        std::cerr << "cuda_trace_append(" << idx << " <- " << arg1 << "): " << cmd
                  << std::endl;
#endif

    Variable &v = ctx.append(type);
    v.size = v1.size;
    v.dep[0] = arg1;
    v.cmd = cmd;
    v.subtree_size = v1.subtree_size + 1;
    cuda_inc_ref_int(arg1);
    cuda_inc_ref_ext(idx);
    ctx.live.insert(idx);
    strip_ftz(v);
    return idx;
}

ENOKI_EXPORT uint32_t cuda_trace_append(EnokiType type,
                                        const char *cmd,
                                        uint32_t arg1,
                                        uint32_t arg2) {
    if (ENOKI_UNLIKELY(arg1 == 0 || arg2 == 0))
        throw std::runtime_error(
            "cuda_trace_append(): arithmetic involving "
                        "uninitialized variable!");

    Context &ctx = context();
    const Variable &v1 = ctx[arg1],
                   &v2 = ctx[arg2];

    if (v1.dirty || v2.dirty)
        cuda_eval();

    uint32_t idx = ctx.ctr;

#if !defined(NDEBUG)
    if (ctx.log_level >= 4)
        std::cerr << "cuda_trace_append(" << idx << " <- " << arg1 << ", " << arg2
                  << "): " << cmd << std::endl;
#endif

    size_t size = std::max(v1.size, v2.size);
    if (ENOKI_UNLIKELY((v1.size != 1 && v1.size != size) ||
                       (v2.size != 1 && v2.size != size)))
        throw std::runtime_error("cuda_trace_append(): arithmetic involving "
                                 "arrays of incompatible size (" +
                                 std::to_string(v1.size) + " and " + std::to_string(v2.size) +
                                 "). The instruction was \"" + cmd + "\".");

    Variable &v = ctx.append(type);
    v.size = size;
    v.dep = { arg1, arg2, 0 };
    v.cmd = cmd;
    v.subtree_size = v1.subtree_size + v2.subtree_size + 1;
    cuda_inc_ref_int(arg1);
    cuda_inc_ref_int(arg2);
    cuda_inc_ref_ext(idx);
    ctx.live.insert(idx);

    if (v.cmd.find("ld.global") != std::string::npos) {
        v.extra_dep = ctx.scatter_gather_operand;
        cuda_inc_ref_ext(v.extra_dep);
    } else {
        strip_ftz(v);
    }

    return idx;
}

ENOKI_EXPORT uint32_t cuda_trace_append(EnokiType type,
                                        const char *cmd,
                                        uint32_t arg1,
                                        uint32_t arg2,
                                        uint32_t arg3) {
    if (ENOKI_UNLIKELY(arg1 == 0 || arg2 == 0 || arg3 == 0))
        throw std::runtime_error("cuda_trace_append(): arithmetic involving "
                                 "uninitialized variable!");

    Context &ctx = context();
    const Variable &v1 = ctx[arg1],
                   &v2 = ctx[arg2],
                   &v3 = ctx[arg3];

    if (v1.dirty || v2.dirty || v3.dirty)
        cuda_eval();

    uint32_t idx = ctx.ctr;

#if !defined(NDEBUG)
    if (ctx.log_level >= 4)
        std::cerr << "cuda_trace_append(" << idx << " <- " << arg1 << ", " << arg2
                  << ", " << arg3 << "): " << cmd << std::endl;
#endif

    size_t size = std::max(std::max(v1.size, v2.size), v3.size);
    if (ENOKI_UNLIKELY((v1.size != 1 && v1.size != size) ||
                       (v2.size != 1 && v2.size != size) ||
                       (v3.size != 1 && v3.size != size)))
        throw std::runtime_error("cuda_trace_append(): arithmetic involving "
                                 "arrays of incompatible size (" +
                                 std::to_string(v1.size) + ", " + std::to_string(v2.size) +
                                 " and " + std::to_string(v3.size) + "). The instruction was \"" +
                                 cmd + "\".");

    Variable &v = ctx.append(type);
    v.size = size;
    v.dep = { arg1, arg2, arg3 };
    v.cmd = cmd;
    v.subtree_size = v1.subtree_size +
                     v2.subtree_size +
                     v3.subtree_size + 1;
    cuda_inc_ref_int(arg1);
    cuda_inc_ref_int(arg2);
    cuda_inc_ref_int(arg3);
    cuda_inc_ref_ext(idx);
    ctx.live.insert(idx);

    if (v.cmd.find("st.global") != std::string::npos ||
        v.cmd.find("atom.global.add") != std::string::npos) {
        v.extra_dep = ctx.scatter_gather_operand;
        cuda_inc_ref_ext(v.extra_dep);
    } else {
        strip_ftz(v);
    }

    return idx;
}

ENOKI_EXPORT void cuda_trace_printf(const char *fmt, uint32_t narg, uint32_t* arg) {
    auto &ctx = context();
    std::ostringstream oss;
    oss << "{" << std::endl
        << "        .global .align 1 .b8 fmt[] = {";
    for (int i = 0;; ++i) {
        oss << (unsigned) fmt[i];
        if (fmt[i] == '\0')
            break;
        oss << ", ";
    }
    oss << "};" << std::endl
        << "        .local .align 8 .b8 buf[" << 8*narg << "];" << std::endl
        << "        .reg.b64 %fmt_r, %buf_r;" << std::endl;
    for (uint32_t i = 0; i < narg; ++i) {
        switch (ctx[arg[i]].type) {
            case EnokiType::Float32:
                oss << "        cvt.f64.f32 %d0, $r" << i + 2 << ";" << std::endl
                    << "        st.local.f64 [buf + " << i * 8 << "], %d0;" << std::endl;
                break;

            default:
                oss << "        st.local.$t" << i + 2 << " [buf + " << i * 8
                    << "], $r" << i + 2 << ";" << std::endl;
                break;
        }
    }
    oss << "        cvta.global.u64 %fmt_r, fmt;" << std::endl
        << "        cvta.local.u64 %buf_r, buf;" << std::endl
        << "        {" << std::endl
        << "            .param .b64 fmt_p;" << std::endl
        << "            .param .b64 buf_p;" << std::endl
        << "            .param .b32 rv_p;" << std::endl
        << "            st.param.b64 [fmt_p], %fmt_r;" << std::endl
        << "            st.param.b64 [buf_p], %buf_r;" << std::endl
        << "            call.uni (rv_p), vprintf, (fmt_p, buf_p);" << std::endl
        << "        }" << std::endl
        << "    }" << std::endl;

    uint32_t idx = 0;
    if (narg == 0)
        idx = cuda_trace_append(EnokiType::UInt32, oss.str().c_str());
    else if (narg == 1)
        idx = cuda_trace_append(EnokiType::UInt32, oss.str().c_str(), arg[0]);
    else if (narg == 2)
        idx = cuda_trace_append(EnokiType::UInt32, oss.str().c_str(), arg[0], arg[1]);
    else if (narg == 3)
        idx = cuda_trace_append(EnokiType::UInt32, oss.str().c_str(), arg[0], arg[1], arg[2]);
    else
        throw std::runtime_error("cuda_trace_print(): at most 3 variables can be printed at once!");

    cuda_var_mark_side_effect(idx);
    ctx.include_printf = true;
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name JIT trace generation
// -----------------------------------------------------------------------

static void cuda_render_cmd(std::ostringstream &oss,
                            Context &ctx,
                            const std::unordered_map<uint32_t, uint32_t> &reg_map,
                            uint32_t index) {
    const Variable &var = ctx[index];
    const std::string &cmd = var.cmd;

    oss << "    ";
    for (size_t i = 0; i < cmd.length(); ++i) {
        if (cmd[i] != '$' || i + 2 >= cmd.length()) {
            oss << cmd[i];
            continue;
        } else {
            uint8_t type = cmd[i + 1],
                    dep_offset = cmd[i + 2] - '0';

            if (type != 't' && type != 'b' && type != 'r')
                throw std::runtime_error("cuda_render_cmd: invalid '$' template!");

            if (dep_offset < 1 && dep_offset > 4)
                throw std::runtime_error("cuda_render_cmd: out of bounds!");

            uint32_t dep =
                dep_offset == 1 ? index : var.dep[dep_offset - 2];
            EnokiType dep_type = ctx[dep].type;
            const char *result = nullptr;

            if (type == 't')
                result = cuda_register_type(dep_type);
            else if (type == 'b')
                result = cuda_register_type_bin(dep_type);
            else if (type == 'r')
                result = cuda_register_name(dep_type);

            if (result == nullptr)
                throw std::runtime_error(
                    "CUDABackend: internal error -- variable " +
                    std::to_string(index) + " references " + std::to_string(dep) +
                    " with unsupported type: " + std::string(cmd));

            oss << result;

            if (type == 'r') {
                auto it = reg_map.find(dep);
                if (it == reg_map.end())
                    throw std::runtime_error(
                        "CUDABackend: internal error -- variable not found!");
                oss << it->second;
            }

            i += 2;
        }
    }

    if (!cmd.empty() && cmd[cmd.length() - 1] != '\n')
        oss << ";" << std::endl;
}

static std::pair<std::string, std::vector<void *>>
cuda_jit_assemble(size_t size, const std::vector<uint32_t> &sweep, bool include_printf) {
    Context &ctx = context();
    std::ostringstream oss;
    std::vector<void *> ptrs;
    size_t n_in = 0, n_out = 0, n_arith = 0,
           n_total = 0;

    oss << ".version 6.3" << std::endl
        << ".target sm_" << ENOKI_CUDA_COMPUTE_CAPABILITY << std::endl
        << ".address_size 64" << std::endl
        << std::endl;

#if !defined(NDEBUG)
    if (ctx.log_level >= 4)
        std::cerr << "Register map:" << std::endl;
#endif

    uint32_t n_vars = ENOKI_CUDA_REG_RESERVED;
    std::unordered_map<uint32_t, uint32_t> reg_map;
    for (uint32_t index : sweep) {
        const Variable &v = ctx[index];
#if !defined(NDEBUG)
        if (ctx.log_level >= 4) {
            std::cerr << "    " << cuda_register_name(v.type) << n_vars << " -> " << index;
            const std::string &label = v.label;
            if (!label.empty())
                std::cerr << " \"" << label << "\"";
            if (v.size == 1)
                std::cerr << " [scalar]";
            if (v.data != nullptr)
                std::cerr << " [in]";
            else if (v.side_effect)
                std::cerr << " [se]";
            else if (v.ref_count_ext > 0)
                std::cerr << " [out]";
            std::cerr << std::endl;
        }
#endif
        reg_map[index] = n_vars++;

        if (v.data || v.direct_pointer) {
            n_total++;
        } else {
            if (!v.side_effect &&
                 v.ref_count_ext > 0 &&
                 v.size == size)
            n_total++;
        }
    }
    reg_map[2] = 2;


    /**
     * rd0: ptr
     * r1: N
     * r2: index
     * r3: stride
     * r4: threadIdx
     * r5: blockIdx
     * r6: blockDim
     * r7: gridDim
     * rd8, rd9: address I/O
     */

    if (include_printf) {
        oss << ".extern .func  (.param .b32 rv) vprintf (" << std::endl
            << "    .param .b64 fmt," << std::endl
            << "    .param .b64 buf" << std::endl
            << ");" << std::endl
            << std::endl;
    }

    bool parameter_direct = n_total < 128;

    auto param_ref = [parameter_direct](size_t idx) -> std::string {
        if (parameter_direct)
            return "[arg" + std::to_string(idx) + "]";
        else
            return "[%rd0 + " + std::to_string(idx * 8) + "]";
    };

    if (parameter_direct) {
        oss << ".visible .entry enoki_@@@@@@@@(.param .u32 size," << std::endl;
        for (int i = 0; i < n_total; ++i) {
            oss << "                               .param .u64 arg" << i
                << ((i + 1 < n_total) ? "," : ") {") << std::endl;
        }
    } else {
        oss << ".visible .entry enoki_@@@@@@@@(.param .u32 size," << std::endl
            << "                               .param .u64 ptr) {" << std::endl;
    }


    oss << "    .reg.b8 %b<" << n_vars << ">;" << std::endl
        << "    .reg.b16 %w<" << n_vars << ">;" << std::endl
        << "    .reg.b32 %r<" << n_vars << ">;" << std::endl
        << "    .reg.b64 %rd<" << n_vars << ">;" << std::endl
        << "    .reg.f32 %f<" << n_vars << ">;" << std::endl
        << "    .reg.f64 %d<" << n_vars << ">;" << std::endl
        << "    .reg.pred %p<" << n_vars << ">;" << std::endl << std::endl
        << std::endl
        << "    // Grid-stride loop setup" << std::endl;

    if (!parameter_direct)
        oss << "    ld.param.u64 %rd0, [ptr];" << std::endl;

    oss << "    ld.param.u32 %r1, [size];" << std::endl
        << "    mov.u32 %r4, %tid.x;" << std::endl
        << "    mov.u32 %r5, %ctaid.x;" << std::endl
        << "    mov.u32 %r6, %ntid.x;" << std::endl
        << "    mad.lo.u32 %r2, %r5, %r6, %r4;" << std::endl
        << "    setp.ge.u32 %p0, %r2, %r1;" << std::endl
        << "    @%p0 bra L0;" << std::endl
        << std::endl
        << "    mov.u32 %r7, %nctaid.x;" << std::endl
        << "    mul.lo.u32 %r3, %r6, %r7;" << std::endl
        << std::endl
        << "L1:" << std::endl
        << "    // Loop body" << std::endl;

    for (uint32_t index : sweep) {
        Variable &var = ctx[index];

        if (var.is_collected() || (var.cmd.empty() && var.data == nullptr && !var.direct_pointer))
            throw std::runtime_error(
                "CUDABackend: found invalid/expired variable " + std::to_string(index) + " in schedule! ");

        if (var.size != 1 && var.size != size)
            throw std::runtime_error(
                "CUDABackend: encountered arrays of incompatible size! (" +
                std::to_string(size) + " vs " + std::to_string(var.size) + ")");

        oss << std::endl;
        if (var.data || var.direct_pointer) {
            size_t idx = ptrs.size();
            ptrs.push_back(var.data);

            oss << std::endl
                << "    // Load register " << cuda_register_name(var.type) << reg_map[index];
            if (!var.label.empty())
                oss << ": " << var.label;
            oss << std::endl;

            if (!var.direct_pointer) {
                oss << "    " << (parameter_direct ? "ld.param.u64 %rd8, " : "ld.global.u64 %rd8, ")
                    << param_ref(idx) << ";" << std::endl;

                const char *load_instr = "ldu";
                if (var.size != 1) {
                    oss << "    mul.wide.u32 %rd9, %r2, " << cuda_register_size(var.type) << ";" << std::endl
                        << "    add.u64 %rd8, %rd8, %rd9;" << std::endl;
                    load_instr = "ld";
                }
                if (var.type != EnokiType::Bool) {
                    oss << "    " << load_instr << ".global." << cuda_register_type(var.type) << " "
                        << cuda_register_name(var.type) << reg_map[index] << ", [%rd8]"
                        << ";" << std::endl;
                } else {
                    oss << "    " << load_instr << ".global.u8 %w1, [%rd8];" << std::endl
                        << "    setp.ne.u16 " << cuda_register_name(var.type) << reg_map[index] << ", %w1, 0;";
                }
            } else {
                oss << "    " << (parameter_direct ? "ld.param.u64 " : "ldu.global.u64 ")
                    << cuda_register_name(var.type)
                    << reg_map[index] << ", " << param_ref(idx) << ";";
            }

            n_in++;
        } else {
            if (!var.label.empty())
                oss << "    // Compute register "
                    << cuda_register_name(var.type) << reg_map[index] << ": "
                    << var.label << std::endl;
            cuda_render_cmd(oss, ctx, reg_map, index);
            n_arith++;

            if (var.side_effect) {
                n_out++;
                continue;
            }

            if (var.ref_count_ext == 0)
                continue;

            if (var.size != size)
                continue;

            size_t size_in_bytes =
                cuda_var_size(index) * cuda_register_size(var.type);

            var.data = cuda_malloc(size_in_bytes);
            var.subtree_size = 1;
#if !defined(NDEBUG)
            if (ctx.log_level >= 4)
                std::cerr << "cuda_eval(): allocated variable " << index
                          << " -> " << var.data << " (" << size_in_bytes
                          << " bytes)" << std::endl;
#endif
            size_t idx = ptrs.size();
            ptrs.push_back(var.data);
            n_out++;

            oss << std::endl
                << "    // Store register " << cuda_register_name(var.type) << reg_map[index];
            if (!var.label.empty())
                oss << ": " << var.label;
            oss << std::endl
                << "    " << (parameter_direct ? "ld.param.u64 %rd8, " : "ld.global.u64 %rd8, ")
                << param_ref(idx) << ";" << std::endl;
            if (var.size != 1) {
                oss << "    mul.wide.u32 %rd9, %r2, " << cuda_register_size(var.type) << ";" << std::endl
                    << "    add.u64 %rd8, %rd8, %rd9;" << std::endl;
            }
            if (var.type != EnokiType::Bool) {
                oss << "    st.global." << cuda_register_type(var.type) << " [%rd8], "
                    << cuda_register_name(var.type) << reg_map[index] << ";"
                    << std::endl;
            } else {
                oss << "    selp.u16 %w1, 1, 0, " << cuda_register_name(var.type)
                    << reg_map[index] << ";" << std::endl;
                oss << "    st.global.u8" << " [%rd8], %w1;" << std::endl;
            }
        }
    }

    oss << std::endl
        << "    add.u32     %r2, %r2, %r3;" << std::endl
        << "    setp.ge.u32 %p0, %r2, %r1;" << std::endl
        << "    @!%p0 bra L1;" << std::endl;

    oss << std::endl
        << "L0:" << std::endl
        << "    ret;" << std::endl
        << "}" << std::endl;

    if (ctx.log_level >= 1)
        std::cerr << "cuda_eval(): launching kernel (n=" << size << ", in="
                  << n_in << ", out=" << n_out << ", ops=" << n_arith
                  << ")" << std::endl;

    assert(ptrs.size() == n_total);

    return { oss.str(), ptrs };
}

ENOKI_EXPORT void cuda_jit_run(Context &ctx,
                               std::string &&source_,
                               const std::vector<void *> &ptrs,
                               uint32_t size,
                               uint32_t stream_idx,
                               TimePoint start,
                               TimePoint mid) {
    if (source_.empty())
        return;

    uint32_t hash = (uint32_t) StringHasher()(source_);
    char kernel_name[9];
    snprintf(kernel_name, 9, "%08x", hash);
    char *id = strchr((char *) source_.c_str(), '@');
    memcpy(id, kernel_name, 8);

    auto hash_entry = ctx.kernels.emplace(
        std::move(source_), std::pair<CUmodule, CUfunction>{ nullptr, nullptr });
    const std::string &source = hash_entry.first->first;
    CUmodule &module = hash_entry.first->second.first;
    CUfunction &kernel = hash_entry.first->second.second;

    if (ctx.log_level >= 3)
        std::cout << source << std::endl;

    size_t duration_1 = std::chrono::duration_cast<
            std::chrono::microseconds>(mid - start).count();

    if (hash_entry.second) {
        CUjit_option arg[5];
        void *argv[5];
        char error_log[8192], info_log[8192];
        uintptr_t log_size = 8192;
        arg[0] = CU_JIT_INFO_LOG_BUFFER;
        argv[0] = (void *) info_log;
        arg[1] = CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES;
        argv[1] = (void *) log_size;
        arg[2] = CU_JIT_ERROR_LOG_BUFFER;
        argv[2] = (void *) error_log;
        arg[3] = CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES;
        argv[3] = (void *) log_size;
        arg[4] = CU_JIT_LOG_VERBOSE;
        argv[4] = (void *) 1;

        CUlinkState link_state;
        cuda_check(cuLinkCreate(5, arg, argv, &link_state));

        int rt = cuLinkAddData(link_state, CU_JIT_INPUT_PTX, (void *) source.c_str(),
                               source.length(), nullptr, 0, nullptr, nullptr);
        if (rt != CUDA_SUCCESS) {
            std::cerr << "PTX linker error:" << std::endl << error_log << std::endl;
            exit(EXIT_FAILURE);
        }

        void *link_output = nullptr;
        size_t link_output_size = 0;
        cuda_check(cuLinkComplete(link_state, &link_output, &link_output_size));
        if (rt != CUDA_SUCCESS) {
            std::cerr << "PTX linker error:" << std::endl << error_log << std::endl;
            exit(EXIT_FAILURE);
        }

        TimePoint end = std::chrono::high_resolution_clock::now();
        size_t duration_2 = std::chrono::duration_cast<
                std::chrono::microseconds>(end - mid).count();

        if (ctx.log_level >= 2) {
            char *ptax_details = strstr(info_log, "ptxas info");
            char *details = strstr(info_log, "\ninfo    : used");
            if (details) {
                details += 16;
                char *details_len = strstr(details, "registers,");
                if (details_len)
                    details_len[9] = '\0';
                std::cerr << "cuda_jit_run(): "
                    << ((ptax_details == nullptr) ? "cache hit, " : "cache miss, ")
                    << "jit: " << time_string(duration_1)
                    << ", ptx compilation: " << time_string(duration_2)
                    << ", " << details << std::endl;
            }
            if (ctx.log_level >= 3)
                std::cerr << "Detailed linker output:" << std::endl
                          << info_log << std::endl;
        }

        CUresult ret = cuModuleLoadData(&module, link_output);
        if (ret == CUDA_ERROR_OUT_OF_MEMORY) {
            cuda_malloc_trim();
            ret = cuModuleLoadData(&module, link_output);
        }
        cuda_check(ret);

        // Locate the kernel entry point
        cuda_check(cuModuleGetFunction(&kernel, module, (std::string("enoki_") + kernel_name).c_str()));

        // Destroy the linker invocation
        cuda_check(cuLinkDestroy(link_state));
    } else {
        if (ctx.log_level >= 2) {
            std::cerr << "cuda_jit_run(): cache hit, jit: "
                      << time_string(duration_1) << std::endl;
        }
    }

    cudaStream_t cuda_stream = nullptr;
    #if ENOKI_CUDA_USE_STREAMS == 1
        cuda_stream = ctx.streams[stream_idx].stream;
    #endif

    bool parameter_direct = ptrs.size() < 128;
    uint8_t arg_buffer[4096];
    size_t arg_buffer_size = parameter_direct ? (sizeof(uint64_t) * (ptrs.size() + 1)) :
                                                (sizeof(uint64_t) * 2);

    void *config[] = {
        CU_LAUNCH_PARAM_BUFFER_POINTER, arg_buffer,
        CU_LAUNCH_PARAM_BUFFER_SIZE,    &arg_buffer_size,
        CU_LAUNCH_PARAM_END
    };

    memcpy(arg_buffer, &size, sizeof(uint32_t));

    void *host_args = nullptr, *device_args = nullptr;
    if (parameter_direct) {
        memcpy(arg_buffer + sizeof(uint64_t), ptrs.data(), sizeof(uint64_t) * ptrs.size());
    } else {
        size_t total_size = ptrs.size() * sizeof(void*);

        host_args   = cuda_host_malloc(total_size);
        device_args = cuda_malloc(total_size);
        memcpy(host_args, ptrs.data(), total_size);
        cuda_check(cudaMemcpyAsync(device_args, host_args, total_size,
                                   cudaMemcpyHostToDevice, cuda_stream));

        memcpy(arg_buffer + sizeof(uint64_t), &device_args, sizeof(uint64_t));
    }

    uint32_t thread_count = ctx.thread_count,
             block_count = ctx.block_count;

    if (size == 1)
        thread_count = block_count = 1;

    cuda_check(cuLaunchKernel(kernel, block_count, 1, 1, thread_count,
                              1, 1, 0, cuda_stream, nullptr, config));

    if (!parameter_direct) {
        cuda_host_free(host_args, cuda_stream);
        cuda_free(device_args, cuda_stream);
    }

    #if ENOKI_CUDA_LAUNCH_BLOCKING == 1
        cuda_check(cudaStreamSynchronize(cuda_stream));
    #endif
}

static void sweep_recursive(Context &ctx,
                            std::unordered_set<uint32_t> &visited,
                            std::vector<uint32_t> &sweep,
                            uint32_t idx) {
    if (visited.find(idx) != visited.end())
        return;
    visited.insert(idx);

    std::array<uint32_t, 3> deps = ctx[idx].dep;

    auto prio = [&](uint32_t i) -> uint32_t {
        uint32_t k = deps[i];
        if (k >= ENOKI_CUDA_REG_RESERVED)
            return ctx[k].subtree_size;
        else
            return 0;
    };

    if (prio(1) < prio(2))
        std::swap(deps[1], deps[2]);
    if (prio(0) < prio(2))
        std::swap(deps[0], deps[2]);
    if (prio(0) < prio(1))
        std::swap(deps[0], deps[1]);

    for (uint32_t k : deps) {
        if (k >= ENOKI_CUDA_REG_RESERVED)
            sweep_recursive(ctx, visited, sweep, k);
    }

    sweep.push_back(idx);
}

ENOKI_EXPORT void cuda_eval(bool log_assembly) {
    Context &ctx = context();

    for (auto callback: ctx.callbacks)
        callback.first(callback.second);

    std::map<size_t, std::pair<std::unordered_set<uint32_t>,
                               std::vector<uint32_t>>> sweeps;
    for (uint32_t idx : ctx.live) {
        auto &sweep = sweeps[ctx[idx].size];
        sweep_recursive(ctx, std::get<0>(sweep), std::get<1>(sweep), idx);
    }
    for (uint32_t idx : ctx.dirty)
        ctx[idx].dirty = false;

    ctx.live.clear();
    ctx.dirty.clear();

    if (ctx.log_level >= 2 && sweeps.size() > 1)
        std::cerr << "cuda_eval(): begin parallel group" << std::endl;

    #if ENOKI_CUDA_USE_STREAMS == 1
        if (ctx.streams.size() < sweeps.size()) {
            size_t cur = ctx.streams.size();
            ctx.streams.resize(sweeps.size());
            for (size_t i = cur; i < ctx.streams.size(); ++i)
                ctx.streams[i].init();
        }

        cuda_check(cudaEventRecord(ctx.stream_0_event, nullptr));
    #endif

    size_t stream_idx = 0;
    for (auto it = sweeps.rbegin(); it != sweeps.rend(); ++it) {
        size_t size = std::get<0>(*it);
        const std::vector<uint32_t> &schedule = std::get<1>(std::get<1>(*it));

        #if ENOKI_CUDA_USE_STREAMS == 1
            Stream &stream = ctx.streams[stream_idx];
            cuda_check(cudaStreamWaitEvent(stream.stream, ctx.stream_0_event, 0));
        #endif

        TimePoint start = std::chrono::high_resolution_clock::now();
        auto result = cuda_jit_assemble(size, schedule, ctx.include_printf);
        if (std::get<0>(result).empty())
            continue;
        TimePoint mid = std::chrono::high_resolution_clock::now();

        cuda_jit_run(ctx,
                     std::move(std::get<0>(result)),
                     std::get<1>(result),
                     (uint32_t) size, (uint32_t) stream_idx,
                     start, mid);

        #if ENOKI_CUDA_USE_STREAMS == 1
            cuda_check(cudaEventRecord(stream.event, stream.stream));
            cuda_check(cudaStreamWaitEvent(nullptr, stream.event, 0));
        #endif

        stream_idx++;
    }
    ctx.include_printf = false;

    if (ctx.log_level >= 2 && sweeps.size() > 1)
        std::cerr << "cuda_eval(): end parallel group" << std::endl;

    for (auto const &sweep : sweeps) {
        const std::vector<uint32_t> &schedule =
            std::get<1>(std::get<1>(sweep));

        for (uint32_t idx : schedule) {
            auto it = ctx.variables.find(idx);
            if (it == ctx.variables.end())
                continue;

            Variable &v = it->second;

            if (v.data != nullptr && !v.cmd.empty()) {
                for (int j = 0; j < 3; ++j) {
                    cuda_dec_ref_int(v.dep[j]);
                    v.dep[j] = 0;
                }
                cuda_dec_ref_ext(v.extra_dep);
                v.extra_dep = 0;
            }

            if (v.side_effect)
                cuda_dec_ref_ext(idx);
        }
    }
}

ENOKI_EXPORT void cuda_eval_var(uint32_t index, bool log_assembly) {
    Variable &var = context()[index];
    if (var.data == nullptr || var.dirty)
        cuda_eval(log_assembly);
    assert(!var.dirty);
}

//! @}
// -----------------------------------------------------------------------

ENOKI_EXPORT void cuda_fetch_element(void *dst, uint32_t src, size_t offset, size_t size) {
    Variable &var = context()[src];

    if (var.data == nullptr || var.dirty)
        cuda_eval();

    if (var.dirty)
        throw std::runtime_error("cuda_fetch_element(): element is still "
                                 "marked as 'dirty' even after cuda_eval()!");
    else if (var.data == nullptr)
        throw std::runtime_error(
            "cuda_fetch_element(): tried to read from invalid/uninitialized CUDA array!");

    if (var.size == 1)
        offset = 0;

    cuda_check(cudaMemcpy(dst, (uint8_t *) var.data + size * offset,
                          size, cudaMemcpyDeviceToHost));
}

ENOKI_EXPORT void cuda_set_log_level(uint32_t level) {
#if defined(NDEBUG)
    if (level >= 4)
        throw std::runtime_error("cuda_set_log_level(): log levels >= 4 are only supported when Enoki is compiled in debug mode!");
#endif
    context().log_level = level;
}

ENOKI_EXPORT uint32_t cuda_log_level() {
    return context().log_level;
}

ENOKI_EXPORT void cuda_register_callback(void (*callback)(void *), void *payload) {
    context().callbacks.emplace_back(callback, payload);
}

ENOKI_EXPORT void cuda_unregister_callback(void (*callback)(void *), void *payload) {
    auto &cb = context().callbacks;
    auto it = std::find(cb.begin(), cb.end(), std::make_pair(callback, payload));
    if (it == cb.end())
        throw std::runtime_error("cuda_unregister_callback(): entry not found!");
    cb.erase(it);
}

ENOKI_EXPORT char *cuda_whos() {
    std::ostringstream oss;
    oss << std::endl
        << "  ID        Type   E/I Refs   Size        Memory     Ready    Label" << std::endl
        << "  =================================================================" << std::endl;
    auto &ctx = context();

    std::vector<uint32_t> indices;
    indices.reserve(ctx.variables.size());
    for (const auto& it : ctx.variables)
        indices.push_back(it.first);
    std::sort(indices.begin(), indices.end());

    size_t mem_size_scheduled = 0,
           mem_size_ready = 0,
           mem_size_arith = 0;
    for (uint32_t id : indices) {
        if (id < ENOKI_CUDA_REG_RESERVED)
            continue;
        const Variable &v = ctx[id];
        oss << "  " << std::left << std::setw(9) << id << " ";
        switch (v.type) {
            case EnokiType::Int8:    oss << "i8 "; break;
            case EnokiType::UInt8:   oss << "u8 "; break;
            case EnokiType::Int16:   oss << "i16"; break;
            case EnokiType::UInt16:  oss << "u16"; break;
            case EnokiType::Int32:   oss << "i32"; break;
            case EnokiType::UInt32:  oss << "u32"; break;
            case EnokiType::Int64:   oss << "i64"; break;
            case EnokiType::UInt64:  oss << "u64"; break;
            case EnokiType::Float16: oss << "f16"; break;
            case EnokiType::Float32: oss << "f32"; break;
            case EnokiType::Float64: oss << "f64"; break;
            case EnokiType::Bool:    oss << "msk"; break;
            case EnokiType::Pointer: oss << "ptr"; break;
            default: throw std::runtime_error("Invalid array type!");
        }
        size_t mem_size = v.size * cuda_register_size(v.type);
        oss << "    ";
        oss << std::left << std::setw(10) << (std::to_string(v.ref_count_ext) + " / " + std::to_string(v.ref_count_int)) << " ";
        oss << std::left << std::setw(12) << v.size;
        oss << std::left << std::setw(12) << mem_string(mem_size);
        oss << (v.data ? "[x]" : "[ ]") << "     ";
        oss << v.label;
        oss << std::endl;
        if (v.data) {
            mem_size_ready += mem_size;
        } else {
            if (v.ref_count_ext == 0)
                mem_size_arith += mem_size;
            else
                mem_size_scheduled += mem_size;
        }
    }

    oss << "  =================================================================" << std::endl << std::endl
        << "  Memory usage (ready)     : " << mem_string(mem_size_ready)     << std::endl
        << "  Memory usage (scheduled) : " << mem_string(mem_size_ready) << " + "
        << mem_string(mem_size_scheduled) << " = " << mem_string(mem_size_ready + mem_size_scheduled) << std::endl
        << "  Memory savings           : " << mem_string(mem_size_arith)     << std::endl << std::endl
        << "  cuda_malloc() usage: "
        << mem_string(ctx.used) << " device, "
        << mem_string(ctx.used_managed) << " managed, "
        << mem_string(ctx.used_host) << " host." << std::endl
        << "           max. usage: "
        << mem_string(ctx.watermark) << " device, "
        << mem_string(ctx.watermark_managed) << " managed, "
        << mem_string(ctx.watermark_host) << " host." << std::endl;

    return strdup(oss.str().c_str());
}

ENOKI_EXPORT void cuda_malloc_trim() {
    std::unordered_multimap<TaggedSize, void *, TaggedSizeHasher> free_map;

    Context &ctx = context();
    /* Critical section */ {
        std::lock_guard<std::recursive_mutex> guard(ctx.malloc_mutex);
        free_map.swap(ctx.free_map);
    }

    size_t freed_normal = 0, freed_managed = 0, freed_host = 0,
           count = free_map.size();

    for (auto kv : free_map) {
        switch (kv.first.type) {
            case Normal:
                freed_normal += kv.first.size;
                cuda_check(cudaFree(kv.second));
                break;

            case Managed:
                freed_managed += kv.first.size;
                cuda_check(cudaFree(kv.second));
                break;

            case Host:
                freed_host += kv.first.size;
                cuda_check(cudaFreeHost(kv.second));
                break;

            default:
                throw std::runtime_error("cuda_malloc_trim(): internal error!");
        }
    }
    if (ctx.log_level >= 4 && count > 0)
        std::cerr << "cuda_malloc_trim(): freed " << count << " arrays ("
                  << mem_string(freed_normal) << " device memory, "
                  << mem_string(freed_managed) << " unified memory, and "
                  << mem_string(freed_host) << " host memory)." << std::endl;
}

ENOKI_EXPORT void cuda_sync() {
    Context &ctx = context();
    if (ctx.log_level >= 4)
        std::cerr << "cuda_sync()." << std::endl;
    cuda_check(cudaDeviceSynchronize());
}

size_t malloc_round(size_t x) {
    /* Round to next higher power of two */
    x -= 1;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    x += 1;
    return x;
}

ENOKI_EXPORT void* cuda_malloc(size_t size) {
    if (size == 0)
        return nullptr;

    size = malloc_round(size);
    TaggedSize ts(Normal, size);
    void *ptr = nullptr;

    Context &ctx = context();
    /* Critical section */ {
        std::lock_guard<std::recursive_mutex> guard(ctx.malloc_mutex);
        auto it = ctx.free_map.find(ts);

        if (it != ctx.free_map.end()) {
            ptr = it->second;
            ctx.free_map.erase(it);
        }
    }

    if (ptr == nullptr) {
        cudaError_t ret = cudaMalloc(&ptr, size);
        if (ret != cudaSuccess) {
            cuda_sync();
            cuda_malloc_trim();
            cudaError_t ret = cudaMalloc(&ptr, size);
            if (ret != cudaSuccess)
                throw std::runtime_error("cuda_malloc(): out of memory!");
        }
    }

    /* Critical section */ {
        std::lock_guard<std::recursive_mutex> guard(ctx.malloc_mutex);
        auto result = ctx.used_map.insert(std::make_pair(ptr, ts));
        if (!result.second) {
            fprintf(stderr, "cuda_malloc(): internal error!\n");
            exit(EXIT_FAILURE);
        }

        ctx.used += size;
        ctx.watermark = std::max(ctx.watermark, ctx.used);
    }

    return ptr;
}

ENOKI_EXPORT void* cuda_managed_malloc(size_t size) {
    if (size == 0)
        return nullptr;

    size = malloc_round(size);
    TaggedSize ts(Managed, size);
    void *ptr = nullptr;

    Context &ctx = context();
    /* Critical section */ {
        std::lock_guard<std::recursive_mutex> guard(ctx.malloc_mutex);
        auto it = ctx.free_map.find(ts);

        if (it != ctx.free_map.end()) {
            ptr = it->second;
            ctx.free_map.erase(it);
        }
    }

    if (ptr == nullptr) {
        cudaError_t ret = cudaMallocManaged(&ptr, size);
        if (ret != cudaSuccess) {
            cuda_sync();
            cuda_malloc_trim();
            cudaError_t ret = cudaMallocManaged(&ptr, size);
            if (ret != cudaSuccess)
                throw std::runtime_error("cuda_managed_malloc(): out of memory!");
        }
    }

    /* Critical section */ {
        std::lock_guard<std::recursive_mutex> guard(ctx.malloc_mutex);
        auto result = ctx.used_map.insert(std::make_pair(ptr, ts));
        if (!result.second) {
            fprintf(stderr, "cuda_managed_malloc(): internal error!\n");
            exit(EXIT_FAILURE);
        }

        ctx.used_managed += size;
        ctx.watermark_managed = std::max(ctx.watermark_managed, ctx.used_managed);
    }

    return ptr;
}

ENOKI_EXPORT void* cuda_host_malloc(size_t size) {
    if (size == 0)
        return nullptr;

    size = malloc_round(size);
    TaggedSize ts(Host, size);
    void *ptr = nullptr;

    Context &ctx = context();
    /* Critical section */ {
        std::lock_guard<std::recursive_mutex> guard(ctx.malloc_mutex);
        auto it = ctx.free_map.find(ts);

        if (it != ctx.free_map.end()) {
            ptr = it->second;
            ctx.free_map.erase(it);
        }
    }

    if (ptr == nullptr) {
        cudaError_t ret = cudaMallocHost(&ptr, size);
        if (ret != cudaSuccess) {
            cuda_sync();
            cuda_malloc_trim();
            cudaError_t ret = cudaMallocHost(&ptr, size);
            if (ret != cudaSuccess)
                throw std::runtime_error("cuda_host_malloc(): out of memory!");
        }
    }

    /* Critical section */ {
        std::lock_guard<std::recursive_mutex> guard(ctx.malloc_mutex);
        auto result = ctx.used_map.insert(std::make_pair(ptr, ts));
        if (!result.second) {
            fprintf(stderr, "cuda_host_malloc(): internal error!\n");
            exit(EXIT_FAILURE);
        }

        ctx.used_host += size;
        ctx.watermark_host = std::max(ctx.watermark_host, ctx.used_host);
    }

    return ptr;
}

ENOKI_EXPORT void cuda_free(void *ptr, cudaStream_t stream) {
    if (ptr == nullptr)
        return;

    cudaStreamAddCallback(
        stream, [](cudaStream_t stream, cudaError_t status, void *data) {
            Context &ctx = context();
            std::lock_guard<std::recursive_mutex> guard(ctx.malloc_mutex);
            auto it = ctx.used_map.find(data);

            if (it == ctx.used_map.end()) {
                fprintf(stderr, "cuda_host_free(): unknown/unregistered pointer!\n");
                exit(EXIT_FAILURE);
            }

            if (it->second.type == Normal) {
                ctx.used -= it->second.size;
            } else if (it->second.type == Managed) {
                ctx.used_managed -= it->second.size;
            } else {
                fprintf(stderr, "cuda_host_free(): tried to free a host pointer!");
                exit(EXIT_FAILURE);
            }

            ctx.free_map.insert(std::make_pair(it->second, data));
            ctx.used_map.erase(it);
        },
        ptr, 0
    );
}

ENOKI_EXPORT void cuda_free(void *ptr) {
    cuda_free(ptr, nullptr);
}

ENOKI_EXPORT void cuda_host_free(void *ptr, cudaStream_t stream) {
    if (ptr == nullptr)
        return;

    cudaStreamAddCallback(
        stream, [](cudaStream_t stream, cudaError_t status, void *data) {
            Context &ctx = context();
            std::lock_guard<std::recursive_mutex> guard(ctx.malloc_mutex);
            auto it = ctx.used_map.find(data);

            if (it == ctx.used_map.end()) {
                fprintf(stderr, "cuda_host_free(): unknown/unregistered pointer!\n");
                exit(EXIT_FAILURE);
            }

            if (it->second.type != Host) {
                fprintf(stderr, "cuda_host_free(): tried to free a device pointer!");
                exit(EXIT_FAILURE);
            }

            ctx.free_map.insert(std::make_pair(it->second, data));
            ctx.used_host -= it->second.size;
            ctx.used_map.erase(it);
        },
        ptr, 0
    );
}

ENOKI_EXPORT void cuda_host_free(void *ptr) {
    cuda_host_free(ptr, nullptr);
}

NAMESPACE_END(enoki)
