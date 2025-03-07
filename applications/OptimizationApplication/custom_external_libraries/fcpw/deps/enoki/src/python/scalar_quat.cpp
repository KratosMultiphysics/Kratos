#include "quat.h"

void bind_scalar_quaternion(py::module& m, py::module& s) {
    bind_quaternion<Quaternion4f>(m, s, "Quaternion4f");
    bind_quaternion<Quaternion4d>(m, s, "Quaternion4d");
}
