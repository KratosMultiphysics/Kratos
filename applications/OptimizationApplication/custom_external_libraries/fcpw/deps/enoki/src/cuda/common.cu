/*
    src/cuda/common.cpp -- CUDA backend (wrapper routines)

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyrighe (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include <cuda.h>
#include <stdio.h>
#include <iostream>
#include "common.cuh"

NAMESPACE_BEGIN(enoki)

std::string mem_string(size_t size) {
    const char *orders[] = {
        "B", "KiB", "MiB", "GiB",
        "TiB", "PiB", "EiB"
    };
    float value = (float) size;

    int i = 0;
    for (i = 0; i < 6 && value >= 1024.f; ++i)
        value /= 1024.f;

    char buf[32];
    snprintf(buf, 32, "%.5g %s", value, orders[i]);

    return buf;
}

std::string time_string(size_t value_) {
    struct Order { float factor; const char* suffix; };
    const Order orders[] = { { 0, "us" },   { 1000, "ms" },
                             { 1000, "s" }, { 60, "m" },
                             { 60, "h" },   { 24, "d" },
                             { 7, "w" },    { (float) 52.1429, "y" } };

    int i = 0;
    float value = (float) value_;
    for (i = 0; i < 7 && value > orders[i+1].factor; ++i)
        value /= orders[i+1].factor;

    char buf[32];
    snprintf(buf, 32, "%.5g %s", value, orders[i].suffix);

    return buf;
}

template <typename T> __global__ void fill(T *out, T value, size_t n) {
    for (size_t i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
         i += blockDim.x * gridDim.x)
        out[i] = value;
}

template <typename T> __global__ void set_value(T *ptr, size_t idx, T value) {
    ptr[idx] = value;
}

ENOKI_EXPORT void cuda_fill(uint8_t *ptr, uint8_t value, size_t size) {
    cuda_check(cudaMemsetAsync(ptr, value, size));
}

ENOKI_EXPORT void cuda_fill(uint16_t *ptr, uint16_t value, size_t size) {
    fill<<<256, 256>>>(ptr, value, size);
}

ENOKI_EXPORT void cuda_fill(uint32_t *ptr, uint32_t value, size_t size) {
    fill<<<256, 256>>>(ptr, value, size);
}

ENOKI_EXPORT void cuda_fill(uint64_t *ptr, uint64_t value, size_t size) {
    fill<<<256, 256>>>(ptr, value, size);
}

template <typename T> __global__ void reverse(T *out, const T *in, size_t n) {
    for (size_t i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
         i += blockDim.x * gridDim.x)
        out[i] = in[n - 1 - i];
}

ENOKI_EXPORT void cuda_reverse(uint8_t *out, const uint8_t *in, size_t size) {
    reverse<<<256, 256>>>(out, in, size);
}

ENOKI_EXPORT void cuda_reverse(uint16_t *out, const uint16_t *in, size_t size) {
    reverse<<<256, 256>>>(out, in, size);
}

ENOKI_EXPORT void cuda_reverse(uint32_t *out, const uint32_t *in, size_t size) {
    reverse<<<256, 256>>>(out, in, size);
}

ENOKI_EXPORT void cuda_reverse(uint64_t *out, const uint64_t *in, size_t size) {
    reverse<<<256, 256>>>(out, in, size);
}

ENOKI_EXPORT void cuda_memcpy_to_device(void *dst, const void *src, size_t size) {
    cuda_check(cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice));
}

ENOKI_EXPORT void cuda_memcpy_from_device(void *dst, const void *src, size_t size) {
    cuda_check(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToHost));
}

ENOKI_EXPORT void cuda_memcpy_to_device_async(void *dst, const void *src, size_t size) {
    cuda_check(cudaMemcpyAsync(dst, src, size, cudaMemcpyHostToDevice));
}

ENOKI_EXPORT void cuda_memcpy_from_device_async(void *dst, const void *src, size_t size) {
    cuda_check(cudaMemcpyAsync(dst, src, size, cudaMemcpyDeviceToHost));
}

ENOKI_EXPORT void cuda_mem_get_info(size_t *free, size_t *total) {
    cuda_check(cudaMemGetInfo(free, total));
}

struct CUDAErrorList {
    CUresult id;
    const char *value;
};

static CUDAErrorList __cuda_error_list[] = {
    { CUDA_SUCCESS,
     "CUDA_SUCCESS"},
    { CUDA_ERROR_INVALID_VALUE,
     "CUDA_ERROR_INVALID_VALUE"},
    { CUDA_ERROR_OUT_OF_MEMORY,
     "CUDA_ERROR_OUT_OF_MEMORY"},
    { CUDA_ERROR_NOT_INITIALIZED,
     "CUDA_ERROR_NOT_INITIALIZED"},
    { CUDA_ERROR_DEINITIALIZED,
     "CUDA_ERROR_DEINITIALIZED"},
    { CUDA_ERROR_PROFILER_DISABLED,
     "CUDA_ERROR_PROFILER_DISABLED"},
    { CUDA_ERROR_PROFILER_NOT_INITIALIZED,
     "CUDA_ERROR_PROFILER_NOT_INITIALIZED"},
    { CUDA_ERROR_PROFILER_ALREADY_STARTED,
     "CUDA_ERROR_PROFILER_ALREADY_STARTED"},
    { CUDA_ERROR_PROFILER_ALREADY_STOPPED,
     "CUDA_ERROR_PROFILER_ALREADY_STOPPED"},
    { CUDA_ERROR_NO_DEVICE,
     "CUDA_ERROR_NO_DEVICE"},
    { CUDA_ERROR_INVALID_DEVICE,
     "CUDA_ERROR_INVALID_DEVICE"},
    { CUDA_ERROR_INVALID_IMAGE,
     "CUDA_ERROR_INVALID_IMAGE"},
    { CUDA_ERROR_INVALID_CONTEXT,
     "CUDA_ERROR_INVALID_CONTEXT"},
    { CUDA_ERROR_CONTEXT_ALREADY_CURRENT,
     "CUDA_ERROR_CONTEXT_ALREADY_CURRENT"},
    { CUDA_ERROR_MAP_FAILED,
     "CUDA_ERROR_MAP_FAILED"},
    { CUDA_ERROR_UNMAP_FAILED,
     "CUDA_ERROR_UNMAP_FAILED"},
    { CUDA_ERROR_ARRAY_IS_MAPPED,
     "CUDA_ERROR_ARRAY_IS_MAPPED"},
    { CUDA_ERROR_ALREADY_MAPPED,
     "CUDA_ERROR_ALREADY_MAPPED"},
    { CUDA_ERROR_NO_BINARY_FOR_GPU,
     "CUDA_ERROR_NO_BINARY_FOR_GPU"},
    { CUDA_ERROR_ALREADY_ACQUIRED,
     "CUDA_ERROR_ALREADY_ACQUIRED"},
    { CUDA_ERROR_NOT_MAPPED,
     "CUDA_ERROR_NOT_MAPPED"},
    { CUDA_ERROR_NOT_MAPPED_AS_ARRAY,
     "CUDA_ERROR_NOT_MAPPED_AS_ARRAY"},
    { CUDA_ERROR_NOT_MAPPED_AS_POINTER,
     "CUDA_ERROR_NOT_MAPPED_AS_POINTER"},
    { CUDA_ERROR_ECC_UNCORRECTABLE,
     "CUDA_ERROR_ECC_UNCORRECTABLE"},
    { CUDA_ERROR_UNSUPPORTED_LIMIT,
     "CUDA_ERROR_UNSUPPORTED_LIMIT"},
    { CUDA_ERROR_CONTEXT_ALREADY_IN_USE,
     "CUDA_ERROR_CONTEXT_ALREADY_IN_USE"},
    { CUDA_ERROR_PEER_ACCESS_UNSUPPORTED,
     "CUDA_ERROR_PEER_ACCESS_UNSUPPORTED"},
    { CUDA_ERROR_INVALID_PTX,
     "CUDA_ERROR_INVALID_PTX"},
    { CUDA_ERROR_INVALID_GRAPHICS_CONTEXT,
     "CUDA_ERROR_INVALID_GRAPHICS_CONTEXT"},
    { CUDA_ERROR_NVLINK_UNCORRECTABLE,
     "CUDA_ERROR_NVLINK_UNCORRECTABLE"},
    { CUDA_ERROR_JIT_COMPILER_NOT_FOUND,
     "CUDA_ERROR_JIT_COMPILER_NOT_FOUND"},
    { CUDA_ERROR_INVALID_SOURCE,
     "CUDA_ERROR_INVALID_SOURCE"},
    { CUDA_ERROR_FILE_NOT_FOUND,
     "CUDA_ERROR_FILE_NOT_FOUND"},
    { CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND,
     "CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND"},
    { CUDA_ERROR_SHARED_OBJECT_INIT_FAILED,
     "CUDA_ERROR_SHARED_OBJECT_INIT_FAILED"},
    { CUDA_ERROR_OPERATING_SYSTEM,
     "CUDA_ERROR_OPERATING_SYSTEM"},
    { CUDA_ERROR_INVALID_HANDLE,
     "CUDA_ERROR_INVALID_HANDLE"},
    { CUDA_ERROR_NOT_FOUND,
     "CUDA_ERROR_NOT_FOUND"},
    { CUDA_ERROR_NOT_READY,
     "CUDA_ERROR_NOT_READY"},
    { CUDA_ERROR_ILLEGAL_ADDRESS,
     "CUDA_ERROR_ILLEGAL_ADDRESS"},
    { CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES,
     "CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES"},
    { CUDA_ERROR_LAUNCH_TIMEOUT,
     "CUDA_ERROR_LAUNCH_TIMEOUT"},
    { CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING,
     "CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING"},
    { CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED,
     "CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED"},
    { CUDA_ERROR_PEER_ACCESS_NOT_ENABLED,
     "CUDA_ERROR_PEER_ACCESS_NOT_ENABLED"},
    { CUDA_ERROR_PRIMARY_CONTEXT_ACTIVE,
     "CUDA_ERROR_PRIMARY_CONTEXT_ACTIVE"},
    { CUDA_ERROR_CONTEXT_IS_DESTROYED,
     "CUDA_ERROR_CONTEXT_IS_DESTROYED"},
    { CUDA_ERROR_ASSERT,
     "CUDA_ERROR_ASSERT"},
    { CUDA_ERROR_TOO_MANY_PEERS,
     "CUDA_ERROR_TOO_MANY_PEERS"},
    { CUDA_ERROR_HOST_MEMORY_ALREADY_REGISTERED,
     "CUDA_ERROR_HOST_MEMORY_ALREADY_REGISTERED"},
    { CUDA_ERROR_HOST_MEMORY_NOT_REGISTERED,
     "CUDA_ERROR_HOST_MEMORY_NOT_REGISTERED"},
    { CUDA_ERROR_HARDWARE_STACK_ERROR,
     "CUDA_ERROR_HARDWARE_STACK_ERROR"},
    { CUDA_ERROR_ILLEGAL_INSTRUCTION,
     "CUDA_ERROR_ILLEGAL_INSTRUCTION"},
    { CUDA_ERROR_MISALIGNED_ADDRESS,
     "CUDA_ERROR_MISALIGNED_ADDRESS"},
    { CUDA_ERROR_INVALID_ADDRESS_SPACE,
     "CUDA_ERROR_INVALID_ADDRESS_SPACE"},
    { CUDA_ERROR_INVALID_PC,
     "CUDA_ERROR_INVALID_PC"},
    { CUDA_ERROR_LAUNCH_FAILED,
     "CUDA_ERROR_LAUNCH_FAILED"},
    { CUDA_ERROR_COOPERATIVE_LAUNCH_TOO_LARGE,
     "CUDA_ERROR_COOPERATIVE_LAUNCH_TOO_LARGE"},
    { CUDA_ERROR_NOT_PERMITTED,
     "CUDA_ERROR_NOT_PERMITTED"},
    { CUDA_ERROR_NOT_SUPPORTED,
     "CUDA_ERROR_NOT_SUPPORTED"},
    { CUDA_ERROR_UNKNOWN,
     "CUDA_ERROR_UNKNOWN"},
    { (CUresult) -1, nullptr }
};

ENOKI_EXPORT const char *cuda_error_string(CUresult id) {
    int index = 0;

    while (__cuda_error_list[index].id != id &&
           __cuda_error_list[index].id != (CUresult) -1)
        index++;

    if (__cuda_error_list[index].id == id)
        return __cuda_error_list[index].value;
    else
        return "Invalid CUDA error status!";
}

ENOKI_EXPORT void cuda_check_impl(CUresult errval, const char *file, const int line) {
    if (errval != CUDA_SUCCESS && errval != CUDA_ERROR_DEINITIALIZED) {
        const char *err_msg = cuda_error_string(errval);
        fprintf(stderr,
                "cuda_check(): driver API error = %04d \"%s\" in "
                "%s:%i.\n", (int) errval, err_msg, file, line);
        exit(EXIT_FAILURE);
    }
}

ENOKI_EXPORT void cuda_check_impl(cudaError_t errval, const char *file, const int line) {
    if (errval != cudaSuccess && errval != cudaErrorCudartUnloading) {
        const char *err_msg = cudaGetErrorName(errval);
        fprintf(stderr,
                "cuda_check(): runtime API error = %04d \"%s\" in "
                "%s:%i.\n", (int) errval, err_msg, file, line);
        exit(EXIT_FAILURE);
    }
}

NAMESPACE_END(enoki)
