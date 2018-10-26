#pragma once

#include <iostream>

namespace AMatrix {
constexpr std::size_t dynamic = 0;
constexpr std::size_t row_major_access = 1;
constexpr std::size_t column_major_access = 2;
constexpr std::size_t unordered_access = 3;

template <std::size_t TCategory1, std::size_t TCategory2>
class AccessTrait {
   public:
    static constexpr std::size_t category = unordered_access;
};

template <>
class AccessTrait<row_major_access, row_major_access> {
   public:
    static constexpr std::size_t category = row_major_access;
};

template <typename TExpressionType, std::size_t TCategory = unordered_access>
class MatrixExpression {
   public:
    static constexpr std::size_t category = TCategory;

    // using value_type = TExpressionType::value_type;
    MatrixExpression() {}

    MatrixExpression(MatrixExpression const& Other) = default;

    MatrixExpression(MatrixExpression&& Other) = default;

    MatrixExpression& operator=(MatrixExpression const& Other) = default;

    MatrixExpression& operator=(MatrixExpression&& Other) = default;

    TExpressionType& expression() {
        return *static_cast<TExpressionType*>(this);
    }
    TExpressionType const& expression() const {
        return *static_cast<const TExpressionType* const>(this);
    }

    TExpressionType& noalias() { return expression(); }
};

template <typename TExpressionType>
class TransposeMatrix
    : public MatrixExpression<TransposeMatrix<TExpressionType>> {
    TExpressionType const& _original_expression;

   public:
    using data_type = typename TExpressionType::data_type;
    TransposeMatrix() = delete;

    TransposeMatrix(TExpressionType const& Original)
        : _original_expression(Original) {}

    inline data_type operator()(std::size_t i, std::size_t j) const {
        return _original_expression(j, i);
    }

    inline std::size_t size1() const { return _original_expression.size2(); }
    inline std::size_t size2() const { return _original_expression.size1(); }
};

template <typename TExpressionType>
class MatrixRow : public MatrixExpression<MatrixRow<TExpressionType>> {
    TExpressionType& _original_expression;
    std::size_t _row_index;

   public:
    using data_type = typename TExpressionType::data_type;
    MatrixRow() = delete;

    MatrixRow(TExpressionType& Original, std::size_t RowIndex)
        : _original_expression(Original), _row_index(RowIndex) {}

    template <typename TOtherExpressionType, std::size_t TCategory>
    MatrixRow& operator=(
        MatrixExpression<TOtherExpressionType, TCategory> const& Other) {
        auto& other_expression = Other.expression();
        for (std::size_t j = 0; j < size2(); j++)
            _original_expression(_row_index, j) = Other.expression()(0, j);
        return *this;
    }

    template <typename TOtherExpressionType>
    MatrixRow& operator=(
        MatrixExpression<TOtherExpressionType, row_major_access> const& Other) {
        auto& other_expression = Other.expression();
        for (std::size_t j = 0; j < size2(); j++)
            _original_expression(_row_index, j) = Other.expression()[j];
        return *this;
    }

    inline data_type const& operator()(std::size_t i, std::size_t j) const {
        //assert(i == 0)
        return _original_expression(_row_index, j);
    }

    inline data_type& operator()(std::size_t i, std::size_t j) {
        //assert(i == 0)
        return _original_expression(_row_index, j);
    }

    inline data_type const& operator[](std::size_t i) const {
        return _original_expression(_row_index, i);
    }

    inline data_type& operator[](std::size_t i) {
        return _original_expression(_row_index, i);
    }

    inline std::size_t size() const { return _original_expression.size2(); }
    inline std::size_t size1() const { return 1; }
    inline std::size_t size2() const { return _original_expression.size2(); }
};

template <typename TExpressionType>
class MatrixColumn : public MatrixExpression<MatrixColumn<TExpressionType>> {
    TExpressionType& _original_expression;
    std::size_t _column_index;

   public:
    using data_type = typename TExpressionType::data_type;
    MatrixColumn() = delete;

    MatrixColumn(TExpressionType& Original, std::size_t ColumnIndex)
        : _original_expression(Original), _column_index(ColumnIndex) {}

    template <typename TOtherExpressionType, std::size_t TCategory>
    MatrixColumn& operator=(
        MatrixExpression<TOtherExpressionType, TCategory> const& Other) {
        auto& other_expression = Other.expression();
        for (std::size_t i = 0; i < size1(); i++)
            _original_expression(i, _column_index) = Other.expression()(i, 0);
        return *this;
    }

    template <typename TOtherExpressionType>
    MatrixColumn& operator=(
        MatrixExpression<TOtherExpressionType, row_major_access> const& Other) {
        auto& other_expression = Other.expression();
        for (std::size_t i = 0; i < size1(); i++)
            _original_expression(i, _column_index) = Other.expression()[i];
        return *this;
    }

    inline data_type const& operator()(std::size_t i, std::size_t j) const {
        return _original_expression(i, _column_index);
    }

    inline data_type& operator()(std::size_t i, std::size_t j) {
		//assert(j == 0)
        return _original_expression(i, _column_index);
    }

    inline data_type const& operator[](std::size_t i) const {
        //assert(j == 0)
        return _original_expression(i, _column_index);
    }

    inline data_type& operator[](std::size_t i) {
        return _original_expression(i, _column_index);
    }

    inline std::size_t size() const { return _original_expression.size1(); }
    inline std::size_t size1() const { return _original_expression.size1(); }
    inline std::size_t size2() const { return 1; }
};

template <typename TExpressionType>
class SubMatrix : public MatrixExpression<SubMatrix<TExpressionType>> {
    TExpressionType& _original_expression;
    std::size_t _origin_index1;
    std::size_t _origin_index2;
    std::size_t _size1;
    std::size_t _size2;

   public:
    using data_type = typename TExpressionType::data_type;
    SubMatrix() = delete;

    SubMatrix(TExpressionType& Original, std::size_t OriginIndex1,
        std::size_t OriginIndex2, std::size_t TheSize1, std::size_t TheSize2)
        : _original_expression(Original),
          _origin_index1(OriginIndex1),
          _origin_index2(OriginIndex2),
          _size1(TheSize1),
          _size2(TheSize2) {}

    template <typename TOtherExpressionType, std::size_t TCategory>
    SubMatrix& operator=(
        MatrixExpression<TOtherExpressionType, TCategory> const& Other) {
        auto& other_expression = Other.expression();
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                _original_expression(i + _origin_index1, j + _origin_index2) =
                    Other.expression()(i, j);
        return *this;
    }

    template <typename TOtherExpressionType>
    SubMatrix& operator=(
        MatrixExpression<TOtherExpressionType, row_major_access> const& Other) {
        auto& other_expression = Other.expression();
        std::size_t k = 0;
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                _original_expression(i + _origin_index1, j + _origin_index2) =
                    Other.expression()[k++];
        return *this;
    }

    inline data_type const& operator()(std::size_t i, std::size_t j) const {
        return _original_expression(i + _origin_index1, j + _origin_index2);
    }

    inline data_type& operator()(std::size_t i, std::size_t j) {
        return _original_expression(i + _origin_index1, j + _origin_index2);
    }

    inline std::size_t size() const { return _size1 * _size2; }
    inline std::size_t size1() const { return _size1; }
    inline std::size_t size2() const { return _size2; }
};

template <typename TExpressionType>
class SubVector
    : public MatrixExpression<SubVector<TExpressionType>, row_major_access> {
    TExpressionType& _original_expression;
    std::size_t _origin_index;
    std::size_t _size;

   public:
    using data_type = typename TExpressionType::data_type;
    SubVector() = delete;

    SubVector(
        TExpressionType& Original, std::size_t OriginIndex, std::size_t TheSize)
        : _original_expression(Original),
          _origin_index(OriginIndex),
          _size(TheSize) {}

    template <typename TOtherExpressionType>
    SubVector& operator=(TOtherExpressionType const& Other) {
        for (std::size_t i = 0; i < _size; i++)
            _original_expression[i + _origin_index] = Other[i];

        return *this;
    }

    template <typename TOtherExpressionType>
    SubVector& operator+=(TOtherExpressionType const& Other) {
        for (std::size_t i = 0; i < _size; i++)
            _original_expression[i + _origin_index] += Other[i];

        return *this;
    }

    template <typename TOtherExpressionType>
    SubVector& operator-=(TOtherExpressionType const& Other) {
        for (std::size_t i = 0; i < _size; i++)
            _original_expression[i + _origin_index] -= Other[i];

        return *this;
    }

    inline data_type operator()(std::size_t i, std::size_t j) const {
		//assert(i == 0)
        return _original_expression(0, j + _origin_index);
    }

    inline data_type& operator()(std::size_t i, std::size_t j) {
        //assert(i == 0)
        return _original_expression(0, j + _origin_index);
    }

    inline data_type const& operator[](std::size_t i) const {
        return _original_expression[i + _origin_index];
    }

    inline data_type& operator[](std::size_t i) {
        return _original_expression[i + _origin_index];
    }

    inline std::size_t size() const { return _size; }

    inline std::size_t size1() const { return _size; }
    inline std::size_t size2() const { return 1; }

    data_type* data() { return &_original_expression[_origin_index]; }

    data_type const* data() const { return &_original_expression[_origin_index]; }
};

template <typename TDataType>
class ZeroMatrix
    : public MatrixExpression<ZeroMatrix<TDataType>, row_major_access> {
    std::size_t _size1;
    std::size_t _size2;

   public:
    using data_type = TDataType;

    ZeroMatrix() = delete;

    ZeroMatrix(std::size_t Size) : _size1(Size), _size2(1) {}

    ZeroMatrix(std::size_t Size1, std::size_t Size2)
        : _size1(Size1), _size2(Size2) {}

    inline TDataType operator()(std::size_t i, std::size_t j) const {
        return TDataType();
    }

    inline TDataType operator[](std::size_t i) const { return TDataType(); }

    inline std::size_t size1() const { return _size1; }
    inline std::size_t size2() const { return _size2; }

    inline std::size_t size() const { return _size1 * _size2; }
};

template <typename TDataType>
class IdentityMatrix
    : public MatrixExpression<IdentityMatrix<TDataType>, unordered_access> {
    std::size_t _size;

   public:
    using data_type = TDataType;

    IdentityMatrix() = delete;

    IdentityMatrix(std::size_t Size) : _size(Size) {}

    inline TDataType operator()(std::size_t i, std::size_t j) const {
        return 1.00 * (i == j);
    }

    inline std::size_t size1() const { return _size; }
    inline std::size_t size2() const { return _size; }
    inline std::size_t size() const { return _size * _size; }
};

template <typename TExpression1Type, typename TExpression2Type>
class MatrixSumExpression
    : public MatrixExpression<
          MatrixSumExpression<TExpression1Type, TExpression2Type>,
          AccessTrait<TExpression1Type::category,
              TExpression2Type::category>::category> {
    TExpression1Type const& _first;
    TExpression2Type const& _second;

   public:
    MatrixSumExpression(
        TExpression1Type const& First, TExpression2Type const& Second)
        : _first(First), _second(Second) {}
    using data_type = typename TExpression1Type::data_type;

    std::size_t size1() const { return _first.size1(); }

    std::size_t size2() const { return _first.size2(); }

    std::size_t size() const { return _first.size(); }

    inline data_type operator()(std::size_t i, std::size_t j) const {
        return _first(i, j) + _second(i, j);
    }

    inline data_type operator[](std::size_t i) const {
        return _first[i] + _second[i];
    }
};

template <typename TExpression1Type, typename TExpression2Type,
    std::size_t TCategory1, std::size_t TCategory2>
MatrixSumExpression<TExpression1Type, TExpression2Type> operator+(
    MatrixExpression<TExpression1Type, TCategory1> const& First,
    MatrixExpression<TExpression2Type, TCategory2> const& Second) {
    return MatrixSumExpression<TExpression1Type, TExpression2Type>(
        First.expression(), Second.expression());
}

template <typename TExpression1Type, typename TExpression2Type>
class MatrixMinusExpression
    : public MatrixExpression<
          MatrixMinusExpression<TExpression1Type, TExpression2Type>,
          AccessTrait<TExpression1Type::category,
              TExpression2Type::category>::category> {
    TExpression1Type const& _first;
    TExpression2Type const& _second;

   public:
    MatrixMinusExpression(
        TExpression1Type const& First, TExpression2Type const& Second)
        : _first(First), _second(Second) {}
    using data_type = typename TExpression1Type::data_type;

    std::size_t size1() const { return _first.size1(); }

    std::size_t size2() const { return _first.size2(); }

    std::size_t size() const { return _first.size(); }

    inline data_type operator()(std::size_t i, std::size_t j) const {
        return _first(i, j) - _second(i, j);
    }

    inline data_type operator[](std::size_t i) const {
        return _first[i] - _second[i];
    }
};

template <typename TExpression1Type, typename TExpression2Type,
    std::size_t TCategory1, std::size_t TCategory2>
MatrixMinusExpression<TExpression1Type, TExpression2Type> operator-(
    MatrixExpression<TExpression1Type, TCategory1> const& First,
    MatrixExpression<TExpression2Type, TCategory2> const& Second) {
    return MatrixMinusExpression<TExpression1Type, TExpression2Type>(
        First.expression(), Second.expression());
}

template <typename TExpressionType>
class MatrixUnaryMinusExpression
    : public MatrixExpression<MatrixUnaryMinusExpression<TExpressionType>,
          AccessTrait<TExpressionType::category, row_major_access>::category> {
    TExpressionType const& _original_expression;

   public:
    using data_type = typename TExpressionType::data_type;

    MatrixUnaryMinusExpression(TExpressionType const& TheExpression)
        : _original_expression(TheExpression) {}
    std::size_t size1() const { return _original_expression.size1(); }

    std::size_t size2() const { return _original_expression.size2(); }

    std::size_t size() const { return _original_expression.size(); }

    inline data_type operator()(std::size_t i, std::size_t j) const {
        return -_original_expression(i, j);
    }

    inline data_type operator[](std::size_t i) const {
        return -_original_expression[i];
    }
};

template <typename TExpressionType>
class MatrixScalarProductExpression
    : public MatrixExpression<MatrixScalarProductExpression<TExpressionType>,
          AccessTrait<TExpressionType::category, row_major_access>::category> {
    typename TExpressionType::data_type const& _first;
    TExpressionType const& _second;

   public:
    using data_type = typename TExpressionType::data_type;

    MatrixScalarProductExpression(
        data_type const& First, TExpressionType const& Second)
        : _first(First), _second(Second) {}
    std::size_t size1() const { return _second.size1(); }

    std::size_t size2() const { return _second.size2(); }

    std::size_t size() const { return _second.size(); }

    inline data_type operator()(std::size_t i, std::size_t j) const {
        return _first * _second(i, j);
    }

    inline data_type operator[](std::size_t i) const {
        return _first * _second[i];
    }
};

template <typename TExpressionType, std::size_t TCategory>
MatrixScalarProductExpression<TExpressionType> operator*(
    typename TExpressionType::data_type const& First,
    MatrixExpression<TExpressionType, TCategory> const& Second) {
    return MatrixScalarProductExpression<TExpressionType>(
        First, Second.expression());
}

template <typename TExpressionType, std::size_t TCategory>
MatrixScalarProductExpression<TExpressionType> operator*(
    MatrixExpression<TExpressionType, TCategory> const& First,
    typename TExpressionType::data_type const& Second) {
    return MatrixScalarProductExpression<TExpressionType>(
        Second, First.expression());
}

template <typename TExpressionType>
class MatrixScalarDivisionExpression
    : public MatrixExpression<MatrixScalarDivisionExpression<TExpressionType>,
          AccessTrait<TExpressionType::category, row_major_access>::category> {
    TExpressionType const& _first;
    typename TExpressionType::data_type const _inverse_of_second;

   public:
    using data_type = typename TExpressionType::data_type;

    MatrixScalarDivisionExpression(
        TExpressionType const& First, data_type const& Second)
        : _first(First), _inverse_of_second(data_type(1) / Second) {}
    std::size_t size1() const { return _first.size1(); }

    std::size_t size2() const { return _first.size2(); }

    std::size_t size() const { return _first.size(); }

    inline data_type operator()(std::size_t i, std::size_t j) const {
        return _first(i, j) * _inverse_of_second;
    }

    inline data_type operator[](std::size_t i) const {
        return _first[i] * _inverse_of_second;
    }
};

template <typename TExpressionType, std::size_t TCategory>
MatrixScalarDivisionExpression<TExpressionType> operator/(
    MatrixExpression<TExpressionType, TCategory> const& First,
    typename TExpressionType::data_type const& Second) {
    return MatrixScalarDivisionExpression<TExpressionType>(
        First.expression(), Second);
}

template <typename TExpression1Type, typename TExpression2Type>
class MatrixProductExpression
    : public MatrixExpression<
          MatrixProductExpression<TExpression1Type, TExpression2Type>,
          unordered_access> {
    TExpression1Type const& _first;
    TExpression2Type const& _second;

   public:
    MatrixProductExpression(
        TExpression1Type const& First, TExpression2Type const& Second)
        : _first(First), _second(Second) {}
    using data_type = typename TExpression1Type::data_type;

    std::size_t size1() const { return _first.size1(); }

    std::size_t size2() const { return _second.size2(); }

    std::size_t size() const { return size1() * size2(); }

    inline data_type operator()(std::size_t i, std::size_t j) const {
        data_type result = data_type();
        for (std::size_t k = 0; k < _first.size2(); k++)
            result += _first(i, k) * _second(k, j);
        return result;
    }
};

template <typename TExpression1Type, typename TExpression2Type,
    std::size_t TCategory1, std::size_t TCategory2>
MatrixProductExpression<TExpression1Type, TExpression2Type> operator*(
    MatrixExpression<TExpression1Type, TCategory1> const& First,
    MatrixExpression<TExpression2Type, TCategory2> const& Second) {
    return MatrixProductExpression<TExpression1Type, TExpression2Type>(
        First.expression(), Second.expression());
}

template <typename TExpression1Type, typename TExpression2Type>
class VectorOuterProductExpression
    : public MatrixExpression<
          VectorOuterProductExpression<TExpression1Type, TExpression2Type>,
          unordered_access> {
    TExpression1Type const& _first;
    TExpression2Type const& _second;

   public:
    VectorOuterProductExpression(
        TExpression1Type const& First, TExpression2Type const& Second)
        : _first(First), _second(Second) {}
    using data_type = typename TExpression1Type::data_type;

    std::size_t size1() const { return _first.size(); }

    std::size_t size2() const { return _second.size(); }

    std::size_t size() const { return size() * size(); }

    inline data_type operator()(std::size_t i, std::size_t j) const {
        return _first[i] * _second[j];
    }
};

template <typename TExpression1Type, typename TExpression2Type,
    std::size_t TCategory1, std::size_t TCategory2>
VectorOuterProductExpression<TExpression1Type, TExpression2Type> OuterProduct(
    MatrixExpression<TExpression1Type, TCategory1> const& First,
    MatrixExpression<TExpression2Type, TCategory2> const& Second) {
    return VectorOuterProductExpression<TExpression1Type, TExpression2Type>(
        First.expression(), Second.expression());
}

template <typename TMatrixType, typename TPermutationVectorType>
class LUFactorization
    : public MatrixExpression<
          LUFactorization<TMatrixType, TPermutationVectorType>> {
    TMatrixType& _matrix;
    TPermutationVectorType _permutation_vector;
    std::size_t number_of_pivoting;

   public:
    using data_type = typename TMatrixType::data_type;
    LUFactorization() = delete;

    LUFactorization(TMatrixType& Original) : _matrix(Original) {
        perform_lu();
    }

    inline data_type const& operator()(std::size_t i, std::size_t j) const {
        return _matrix(_permutation_vector[i], j);
    }

    inline std::size_t size1() const { return _matrix.size1(); }
    inline std::size_t size2() const { return _matrix.size2(); }

    /// The algorithm is based on wikipedia implemenation which
    /// can be found in https://en.wikipedia.org/wiki/LU_decomposition
    double determinant() {
        const std::size_t size = size1();
        double result = _matrix(_permutation_vector[0], 0);

        for (std::size_t i = 1; i < size; i++)
            result *= _matrix(_permutation_vector[i], i);

        if ((number_of_pivoting) % 2 == 0)
            return result;
        else
            return -result;
    }

    /// The algorithm is based on wikipedia implemenation which
    /// can be found in https://en.wikipedia.org/wiki/LU_decomposition
    TMatrixType inverse() {
        const std::size_t size = size1();
        TMatrixType result(size, size);

        for (std::size_t j = 0; j < size; j++) {
            for (std::size_t i = 0; i < size; i++) {
                if (_permutation_vector[i] == j)
                    result(i, j) = 1.0;
                else
                    result(i, j) = 0.0;

                for (std::size_t k = 0; k < i; k++)
                    result(i, j) -=
                        _matrix(_permutation_vector[i], k) * result(k, j);
            }

            for (int i = size - 1; i >= 0; i--) {
                for (std::size_t k = i + 1; k < size; k++)
                    result(i, j) -=
                        _matrix(_permutation_vector[i], k) * result(k, j);

                result(i, j) /= _matrix(_permutation_vector[i], i);
            }
        }

		return result;
    }

    /// The algorithm is based on wikipedia implemenation which
    /// can be found in https://en.wikipedia.org/wiki/LU_decomposition
    template <typename TVectorType>
    TVectorType solve(TVectorType const& RHS) {
        const std::size_t size = size1();
        TVectorType result(size);

        for (std::size_t i = 0; i < size; i++) {
            result[i] = RHS[_permutation_vector[i]];

            for (std::size_t k = 0; k < i; k++)
                result[i] -= _matrix(_permutation_vector[i], k) * result[k];
        }

         for (int i = size - 1; i >= 0; i--) {
            for (std::size_t k = i + 1; k < size; k++)
                 result[i] -= _matrix(_permutation_vector[i], k) * result[k];

            result[i] /= _matrix(_permutation_vector[i], i);
        }

        return result;
    }


   private:
    /// The algorithm is based on wikipedia implemenation which
    /// can be found in https://en.wikipedia.org/wiki/LU_decomposition
    int perform_lu() {
        constexpr double tolerance = std::numeric_limits<double>::epsilon();
        std::size_t size1 = _matrix.size1();
        number_of_pivoting = 0;

        initialize_permutation_vector();

        for (std::size_t i = 0; i < size1; i++) {
            double max_pivot = 0.0;
            std::size_t i_max = i;
            double abs_max_pivot = 0.0;

            for (std::size_t k = i; k < size1; k++)
                if ((abs_max_pivot = fabs(_matrix(_permutation_vector[k], i))) >
                    max_pivot) {
                    max_pivot = abs_max_pivot;
                    i_max = k;
                }

            if (max_pivot < tolerance)
                return 0;  // Failed

            if (i_max != i) {
                //pivoting
                std::size_t j = _permutation_vector[i];
                _permutation_vector[i] = _permutation_vector[i_max];
                _permutation_vector[i_max] = j;

                //counting pivots to be used for determinant
                number_of_pivoting++;
            }

            for (std::size_t j = i + 1; j < size1; j++) {
                _matrix(_permutation_vector[j], i) /=
                    _matrix(_permutation_vector[i], i);

                for (std::size_t k = i + 1; k < size1; k++)
                    _matrix(_permutation_vector[j], k) -=
                        _matrix(_permutation_vector[j], i) *
                        _matrix(_permutation_vector[i], k);
            }
        }

        return 1;  // Success
    }

    void initialize_permutation_vector() {
        std::size_t size = _matrix.size1();
        _permutation_vector.resize(size);

        for (std::size_t i = 0; i < size; i++)
            _permutation_vector[i] = i;
    }
};

}  // namespace AMatrix