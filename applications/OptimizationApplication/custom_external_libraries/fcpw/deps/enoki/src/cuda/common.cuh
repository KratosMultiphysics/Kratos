#if !defined(NAMESPACE_BEGIN)
#  define NAMESPACE_BEGIN(name) namespace name {
#endif

#if !defined(NAMESPACE_END)
#  define NAMESPACE_END(name) }
#endif
#
#if defined(_MSC_VER)
#  define ENOKI_EXPORT                 __declspec(dllexport)
#  define ENOKI_LIKELY(x)              x
#  define ENOKI_UNLIKELY(x)            x
#else
#  define ENOKI_EXPORT                 __attribute__ ((visibility("default")))
#  define ENOKI_LIKELY(x)              __builtin_expect(!!(x), 1)
#  define ENOKI_UNLIKELY(x)            __builtin_expect(!!(x), 0)
#endif

#if defined(__SSE4_2__)
#  include <x86intrin.h>
#endif

NAMESPACE_BEGIN(enoki)

enum EnokiType { Invalid = 0, Int8, UInt8, Int16, UInt16,
                 Int32, UInt32, Int64, UInt64, Float16,
                 Float32, Float64, Bool, Pointer };

#define cuda_check(err) cuda_check_impl(err, __FILE__, __LINE__)
ENOKI_EXPORT extern void cuda_check_impl(CUresult errval, const char *file, const int line);
ENOKI_EXPORT extern void cuda_check_impl(cudaError_t errval, const char *file, const int line);

inline uint32_t next_power_of_two(uint32_t n) {
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

ENOKI_EXPORT extern void* cuda_malloc(size_t size);
ENOKI_EXPORT extern void* cuda_managed_malloc(size_t size);
ENOKI_EXPORT extern void* cuda_host_malloc(size_t size);
ENOKI_EXPORT extern void cuda_free(void *p, cudaStream_t stream);
ENOKI_EXPORT extern void cuda_free(void *p);
ENOKI_EXPORT extern void cuda_host_free(void *p, cudaStream_t stream);
ENOKI_EXPORT extern void cuda_host_free(void *p);
ENOKI_EXPORT extern void cuda_malloc_trim();
ENOKI_EXPORT extern void cuda_sync();
ENOKI_EXPORT void cuda_eval(bool log_assembly = false);

extern std::string mem_string(size_t size);
extern std::string time_string(size_t size);

struct StringHasher {
    size_t operator()(const std::string& k) const {
#if defined(__SSE4_2__)
        const char *ptr = k.c_str(),
                   *end = ptr + k.length();

        uint64_t state64 = 0;
        while (ptr + 8 < end) {
            state64 = _mm_crc32_u64(state64, *((uint64_t *) ptr));
            ptr += 8;
        }
        uint32_t state32 = (uint32_t) state64;
        while (ptr < end)
            state32 = _mm_crc32_u8(state32, *ptr++);
        return (size_t) state32;
#else
        return std::hash<std::string>()(k);
#endif
    }
};

#define cuda_check_maybe_redo(expr)                                            \
    for (int i = 0; i < 2; ++i) {                                              \
        cudaError_t rv = expr;                                                 \
        if (rv == cudaErrorMemoryAllocation && i == 0) {                       \
            cuda_malloc_trim();                                                \
        } else {                                                               \
            cuda_check(rv);                                                    \
            break;                                                             \
        }                                                                      \
    }                                                                          \


NAMESPACE_END(enoki)
