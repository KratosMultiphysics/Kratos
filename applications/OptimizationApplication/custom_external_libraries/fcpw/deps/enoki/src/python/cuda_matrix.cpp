#include "matrix.h"

void bind_cuda_matrix(py::module& m, py::module& s) {
    bind_matrix_mask<Matrix2mC>(m, s, "Matrix2m");
    bind_matrix_mask<Matrix3mC>(m, s, "Matrix3m");
    bind_matrix_mask<Matrix4mC>(m, s, "Matrix4m");
    bind_matrix_mask<Matrix44mC>(m, s, "Matrix44m");

    bind_matrix<Matrix2fC>(m, s, "Matrix2f");
    bind_matrix<Matrix3fC>(m, s, "Matrix3f");
    bind_matrix<Matrix4fC>(m, s, "Matrix4f");
    bind_matrix<Matrix44fC>(m, s, "Matrix44f");

    bind_matrix<Matrix2dC>(m, s, "Matrix2d");
    bind_matrix<Matrix3dC>(m, s, "Matrix3d");
    bind_matrix<Matrix4dC>(m, s, "Matrix4d");
    bind_matrix<Matrix44dC>(m, s, "Matrix44d");
}
