#include "matrix.h"

namespace AMatrix {

template <typename TDataType>
using Matrix11 = Matrix<TDataType, 1, 1>;

template <typename TDataType>
using Matrix22 = Matrix<TDataType, 2, 2>;

template <typename TDataType>
using Matrix33 = Matrix<TDataType, 3, 3>;

template <typename TDataType>
using Matrix44 = Matrix<TDataType, 4, 4>;

template <typename TDataType>
using Matrix55 = Matrix<TDataType, 5, 5>;

template <typename TDataType>
using Matrix66 = Matrix<TDataType, 6, 6>;

template <typename TDataType>
using Matrix77 = Matrix<TDataType, 7, 7>;

template <typename TDataType>
using Matrix88 = Matrix<TDataType, 8, 8>;

template <typename TDataType>
using Matrix99 = Matrix<TDataType, 9, 9>;

template <typename TDataType, std::size_t TSize>
using Vector = Matrix<TDataType, TSize, 1>;

template <typename TDataType>
using Vector1 = Matrix<TDataType, 1, 1>;

template <typename TDataType>
using Vector2 = Matrix<TDataType, 2, 1>;

template <typename TDataType>
using Vector3 = Matrix<TDataType, 3, 1>;

template <typename TDataType>
using Vector4 = Matrix<TDataType, 4, 1>;

template <typename TDataType>
using Vector5 = Matrix<TDataType, 5, 1>;

template <typename TDataType>
using Vector6 = Matrix<TDataType, 6, 1>;

template <typename TDataType>
using Vector7 = Matrix<TDataType, 7, 1>;

template <typename TDataType>
using Vector8 = Matrix<TDataType, 8, 1>;

template <typename TDataType>
using Vector9 = Matrix<TDataType, 9, 1>;

template <typename TDataType, std::size_t TSize>
using RowVector = Matrix<TDataType, 1, TSize>;

template <typename TDataType>
using RowVector1 = Matrix<TDataType, 1, 1>;

template <typename TDataType>
using RowVector2 = Matrix<TDataType, 1, 2>;

template <typename TDataType>
using RowVector3 = Matrix<TDataType, 1, 3>;

template <typename TDataType>
using RowVector4 = Matrix<TDataType, 1, 4>;

template <typename TDataType>
using RowVector5 = Matrix<TDataType, 1, 5>;

template <typename TDataType>
using RowVector6 = Matrix<TDataType, 1, 6>;

template <typename TDataType>
using RowVector7 = Matrix<TDataType, 1, 7>;

template <typename TDataType>
using RowVector8 = Matrix<TDataType, 1, 8>;

template <typename TDataType>
using RowVector9 = Matrix<TDataType, 1, 9>;

template <typename TDataType>
using ZeroVector = ZeroMatrix<TDataType>;
}