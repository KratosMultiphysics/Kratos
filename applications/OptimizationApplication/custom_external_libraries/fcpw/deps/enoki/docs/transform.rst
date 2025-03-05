.. cpp:namespace:: enoki

Homogeneous transformations
===========================

Enoki provides a number of convenience functions to construct 3D homogeneous
coordinate transformations (rotations, translations, scales, perspective
transformation matrices, etc.). To use them, include the following header file:

.. code-block:: cpp

    #include <enoki/transform.h>

Reference
---------

.. cpp:function:: template <typename Matrix, typename Vector3> Matrix translate(Vector3 v)

    Constructs a homogeneous coordinate transformation, which translates points by ``v``.

.. cpp:function:: template <typename Matrix, typename Vector3> Matrix scale(Vector3 v)

    Constructs a homogeneous coordinate transformation, which scales points by ``v``.

.. cpp:function:: template <typename Matrix, typename Vector3, typename Float> Matrix rotate(Vector3 v, Float angle)

    Constructs a homogeneous coordinate transformation, which rotates by ``angle`` radians
    around the axis ``v``. The function requires ``v`` to be normalized.

.. cpp:function:: template <typename Matrix> auto transform_decompose(Matrix m)

    Performs a polar decomposition of a non-perspective 4x4 homogeneous
    coordinate matrix and returns a tuple of

    1. A positive definite 3x3 matrix containing an inhomogeneous scaling operation

    2. A rotation quaternion

    3. A 3D translation vector

    This representation is helpful when animating keyframe animations.

    The function also handles singular inputs ``m``, in which case the rotation
    component is set to the identity quaternion and the scaling part simply
    copies the input matrix.

.. cpp:function:: template <typename Matrix3, typename Quaternion, typename Vector3> auto transform_compose(Matrix3 scale, Quaternion rotation, Vector3 translate)

    This function composes a 4x4 homogeneous coordinate transformation from the
    given scale, rotation, and translation. It performs the reverse of
    ``transform_decompose``.

.. cpp:function:: template <typename Matrix3, typename Quaternion, typename Vector3> auto transform_compose_inverse(Matrix3 scale, Quaternion rotation, Vector3 translate)

    This function composes a 4x4 homogeneous *inverse* coordinate
    transformation from the given scale, rotation, and translation. It is the
    equivalent to (but more efficient than) the expression
    ``inverse(transform_compose(...))``.

.. cpp:function:: template <typename Matrix, typename Point3, typename Vector3> Matrix look_at(Point3 origin, Point3, target, Vector3 up)

    Constructs a homogeneous coordinate transformation, which translates to
    :math:`\mathrm{origin}`, maps the negative :math:`z` axis to
    :math:`\mathrm{target}-\mathrm{origin}` (normalized) and the positive
    :math:`y` axis to :math:`\mathrm{up}` (if orthogonal to
    :math:`\mathrm{target}-\mathrm{origin}`). The algorithm performs
    Gram-Schmidt orthogonalization to ensure that the returned matrix is
    orthonormal.

.. cpp:function:: template <typename Matrix, typename Float> Matrix perspective(Float fov, Float near, Float far)

    Constructs an OpenGL-compatible perspective projection matrix with the
    specified field of view (in radians) and near and far clip planes. The
    returned matrix performs the transformation

    .. math::

        \begin{pmatrix}
        x\\y\\z\end{pmatrix}
        \mapsto
        \begin{pmatrix}
        -c\,x/z\\ -c\,x/z\\
        \frac{2\,\mathrm{far}\,\mathrm{near}\,+\,z\,(\mathrm{far}+\mathrm{near})}{z\, (\mathrm{far}-\mathrm{near})}
        \end{pmatrix},

    where

    .. math::

        c = \mathrm{cot}\!\left(0.5\, \textrm{fov}\right),

    which maps :math:`(0, 0, -\mathrm{near})^T` to :math:`(0, 0, -1)^T` and
    :math:`(0, 0, -\mathrm{far})^T` to :math:`(0, 0, 1)^T`.

.. cpp:function:: template <typename Matrix, typename Float> Matrix frustum(Float left, Float right, Float bottom, Float top, Float near, Float far)

    Constructs an OpenGL-compatible perspective projection matrix. The provided
    parameters specify the intersection of the camera frustum with the near
    clipping plane. Specifically, the returned transformation maps
    :math:`(\mathrm{left}, \mathrm{bottom}, -\mathrm{near})` to :math:`(-1, -1,
    -1)` and :math:`(\mathrm{right}, \mathrm{top}, -\mathrm{near})` to
    :math:`(1, 1, -1)`.

.. cpp:function:: template <typename Matrix, typename Float> Matrix ortho(Float left, Float right, Float bottom, Float top, Float near, Float far)

    Constructs an OpenGL-compatible orthographic projection matrix. The
    provided parameters specify the intersection of the camera frustum with the
    near clipping plane. Specifically, the returned transformation maps
    :math:`(\mathrm{left}, \mathrm{bottom}, -\mathrm{near})` to :math:`(-1, -1,
    -1)` and :math:`(\mathrm{right}, \mathrm{top}, -\mathrm{near})` to
    :math:`(1, 1, -1)`.
