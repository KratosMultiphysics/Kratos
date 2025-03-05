#include "matrix.h"

void bind_dynamic_matrix(py::module& m, py::module& s) {
    bind_matrix_mask<Matrix2mX>(m, s, "Matrix2m");
    bind_matrix_mask<Matrix3mX>(m, s, "Matrix3m");
    bind_matrix_mask<Matrix4mX>(m, s, "Matrix4m");
    bind_matrix_mask<Matrix44mX>(m, s, "Matrix44m");

    bind_matrix<Matrix2fX>(m, s, "Matrix2f");
    bind_matrix<Matrix3fX>(m, s, "Matrix3f");
    bind_matrix<Matrix4fX>(m, s, "Matrix4f");
    bind_matrix<Matrix44fX>(m, s, "Matrix44f");

    bind_matrix<Matrix2dX>(m, s, "Matrix2d");
    bind_matrix<Matrix3dX>(m, s, "Matrix3d");
    bind_matrix<Matrix4dX>(m, s, "Matrix4d");
    bind_matrix<Matrix44dX>(m, s, "Matrix44d");
}
