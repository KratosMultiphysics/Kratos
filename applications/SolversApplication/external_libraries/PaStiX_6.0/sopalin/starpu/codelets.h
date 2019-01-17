/**
 *
 * @file codelets.h
 *
 * @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 **/

#ifndef _codelets_h_
#define _codelets_h_

#include "common.h"

#ifdef STARPU_CUDA_ASYNC
#define CODELET_CUDA_FLAGS(flags) .cuda_flags = {(flags)},
#else
#define CODELET_CUDA_FLAGS(flags)
#endif

#define CODELETS_ALL( _unit_, _name_, _nbuffers_, _cpu_func_name_, _cuda_func_name_, _original_location_, _cuda_flags_) \
    struct starpu_codelet cl_##_name_##_##_unit_ = {                    \
        .where     = (_original_location_),                             \
        .cpu_funcs[0] = (_cpu_func_name_),                              \
        CODELET_CUDA_FLAGS(_cuda_flags_)                                \
        .cuda_funcs[0] = (_cuda_func_name_),                            \
        .nbuffers  = (_nbuffers_),                                      \
        .name      = #_name_ ,                                          \
        .model     = &(starpu_##_name_##_model)                         \
    };

#if defined(PASTIX_STARPU_SIMULATION)

#define CODELETS_CPU(_name_, _nbuffers_ )                                  \
    CODELETS_ALL( cpu, _name_, _nbuffers_, (starpu_cpu_func_t) 1, NULL, STARPU_CPU, 0 )

#define CODELETS_GPU(_name_, _nbuffers_, _cuda_flags_)                       \
    CODELETS_ALL( gpu, _name_, _nbuffers_, (starpu_cpu_func_t) 1, (starpu_cuda_func_t) 1, STARPU_CPU | STARPU_CUDA, _cuda_flags_ )

#else

#define CODELETS_CPU(_name_, _nbuffers_ )                                  \
    CODELETS_ALL( cpu, _name_, _nbuffers_, fct_##_name_##_cpu, NULL, STARPU_CPU, 0 )

#define CODELETS_GPU(_name_, _nbuffers_, _cuda_flags_)                       \
    CODELETS_ALL( gpu, _name_, _nbuffers_, fct_##_name_##_cpu, fct_##_name_##_gpu, STARPU_CPU | STARPU_CUDA, _cuda_flags_ )

#endif

#if !defined(PASTIX_WITH_CUDA)
#undef CODELETS_GPU
#define CODELETS_GPU(_name_, _nbuffers_, _cuda_flags_) 
#endif

#endif /* _codelets_h_ */
