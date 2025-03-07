#pragma once

#include <enoki/random.h>
#include "common.h"
#include "docstr.h"

#define D(...) DOC(__VA_ARGS__)

template <typename PCG32>
void bind_pcg32(py::module &m, py::module s, const char *name) {
    using UInt64 = typename PCG32::UInt64;
    using Mask = mask_t<typename PCG32::Float32>;

    py::class_<PCG32>(s, name, D(PCG32))
        .def(py::init<UInt64, UInt64>(),
             "initstate"_a = PCG32_DEFAULT_STATE,
             "initseq"_a = PCG32_DEFAULT_STREAM,
             D(PCG32, PCG32))
        .def("seed", &PCG32::seed,
             "initstate"_a = PCG32_DEFAULT_STATE,
             "initseq"_a = PCG32_DEFAULT_STREAM,
             D(PCG32, seed))
        .def(py::self - py::self, D(PCG32, operator, sub))
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("advance", &PCG32::advance, D(PCG32, advance))
        .def("next_uint32",
             [](PCG32 &pcg) { return pcg.next_uint32(); },
             D(PCG32, next_uint32))
        .def("next_uint32",
             [](PCG32 &pcg, const Mask &mask) {
                 return pcg.next_uint32(mask);
             },
             "mask"_a, D(PCG32, next_uint32, 2))
        .def("next_uint64",
             [](PCG32 &pcg) { return pcg.next_uint64(); },
             D(PCG32, next_uint64))
        .def("next_uint64",
             [](PCG32 &pcg, const Mask &mask) {
                 return pcg.next_uint64(mask);
             },
             "mask"_a, D(PCG32, next_uint64, 2))
        .def("next_uint32_bounded",
             [](PCG32 &pcg, uint32_t bound) {
                 return pcg.next_uint32_bounded(bound);
             },
             "bound"_a, D(PCG32, next_uint32_bounded))
        .def("next_uint32_bounded",
             [](PCG32 &pcg, uint32_t bound, const Mask &mask) {
                 return pcg.next_uint32_bounded(bound, mask);
             },
             "bound"_a, "mask"_a)
        .def("next_uint64_bounded",
             [](PCG32 &pcg, uint64_t bound) {
                 return pcg.next_uint64_bounded(bound);
             },
             "bound"_a, D(PCG32, next_uint64_bounded))
        .def("next_uint64_bounded",
             [](PCG32 &pcg, uint64_t bound, const Mask &mask) {
                 return pcg.next_uint64_bounded(bound, mask);
             },
             "bound"_a, "mask"_a)
        .def("next_float32",
             [](PCG32 &pcg) { return pcg.next_float32(); },
             D(PCG32, next_float32))
        .def("next_float32",
             [](PCG32 &pcg, const Mask &mask) {
                 return pcg.next_float32(mask);
             },
             "mask"_a, D(PCG32, next_float32, 2))
        .def_readwrite("state", &PCG32::state)
        .def_readwrite("inc", &PCG32::inc)
        .def("__repr__", [](const PCG32 &pcg) {
            std::ostringstream oss;
            oss << "PCG32[state=" << pcg.state << ", inc=" << pcg.inc << "]";
            return oss.str();
        });
}
