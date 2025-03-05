#include "matrix.h"
#include <enoki/transform.h>

void bind_scalar_matrix(py::module& m, py::module& s) {
    bind_matrix_mask<Matrix2m>(m, s, "Matrix2m");
    bind_matrix_mask<Matrix3m>(m, s, "Matrix3m");
    bind_matrix_mask<Matrix4m>(m, s, "Matrix4m");
    bind_matrix_mask<Matrix44m>(m, s, "Matrix44m");

    bind_matrix<Matrix2f>(m, s, "Matrix2f");
    bind_matrix<Matrix3f>(m, s, "Matrix3f");
    bind_matrix<Matrix4f>(m, s, "Matrix4f");
    bind_matrix<Matrix44f>(m, s, "Matrix44f");

    bind_matrix<Matrix2d>(m, s, "Matrix2d");
    bind_matrix<Matrix3d>(m, s, "Matrix3d");
    bind_matrix<Matrix4d>(m, s, "Matrix4d");
    bind_matrix<Matrix44d>(m, s, "Matrix44d");

    bind_matrix_mask<Matrix41m>(m, s, "Matrix41m");
    bind_matrix<Matrix41f>(m, s, "Matrix41f");
    bind_matrix<Matrix41d>(m, s, "Matrix41d");

    m.def("transform_decompose", [](const Matrix4f &m) { return transform_decompose(m); });
    m.def("transform_compose",
          [](const Matrix3f &m, const Quaternion<Float32> &q,
             const Array<Float32, 3> &v) { return transform_compose(m, q, v); });
}
