#pragma once

#include "dense_storage.h"

namespace AMatrix {

template <typename TDataType, std::size_t TSize1, std::size_t TSize2>
class MatrixStorage : public DenseStorage<TDataType, TSize1 * TSize2> {
   public:
    using base_type = DenseStorage<TDataType, TSize1 * TSize2>;
    using base_type::at;
    using base_type::data;
    using base_type::size;

    MatrixStorage() {}

    explicit MatrixStorage(std::size_t TheSize1, std::size_t TheSize2) {}

    MatrixStorage(MatrixStorage const& Other) : base_type(Other) {}

    MatrixStorage(MatrixStorage&& Other) = default;

    template <typename TExpressionType, std::size_t TCategory>
    explicit MatrixStorage(
        MatrixExpression<TExpressionType, TCategory> const& Other)
        : base_type(Other) {}

    template <typename TOtherMatrixType>
    explicit MatrixStorage(TOtherMatrixType const& Other) : base_type(Other) {}

    explicit MatrixStorage(std::initializer_list<TDataType> InitialValues)
        : base_type(InitialValues) {}

    template <typename TExpressionType, std::size_t TCategory>
    MatrixStorage& operator=(
        MatrixExpression<TExpressionType, TCategory> const& Other) {
        base_type::operator=(Other);
        return *this;
    }

    template <typename TOtherMatrixType>
    MatrixStorage& operator=(TOtherMatrixType const& Other) {
        base_type::operator=(Other);
        return *this;
    }

    MatrixStorage& operator=(MatrixStorage const& Other) {
        base_type::operator=(Other);
        return *this;
    }

    MatrixStorage& operator=(MatrixStorage&& Other) = default;

    TDataType& operator()(std::size_t i, std::size_t j) { return at(i, j); }

    TDataType const& operator()(std::size_t i, std::size_t j) const {
        return at(i, j);
    }

    TDataType& at(std::size_t i, std::size_t j) { return at(i * TSize2 + j); }

    TDataType const& at(std::size_t i, std::size_t j) const {
        return at(i * TSize2 + j);
    }

    static constexpr std::size_t size1() { return TSize1; }

    static constexpr std::size_t size2() { return TSize2; }

    void resize(std::size_t NewSize) {
    }
};

template <typename TDataType>
class MatrixStorage<TDataType, dynamic, dynamic> {
    std::size_t _size1;
    std::size_t _size2;
    TDataType* _data;

   public:
    MatrixStorage() = delete;

    MatrixStorage(std::size_t TheSize1, std::size_t TheSize2)
        : _size1(TheSize1), _size2(TheSize2) {
        _data = new TDataType[size()];
    }

    explicit MatrixStorage(std::size_t TheSize1, std::size_t TheSize2,
        TDataType const& InitialValue)
        : _size1(TheSize1), _size2(TheSize2) {
        _data = new TDataType[size()];
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = InitialValue;
    }

    MatrixStorage(MatrixStorage const& Other)
        : _size1(Other.size1()), _size2(Other.size2()) {
        _data = new TDataType[size()];
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other._data[i];
    }

    MatrixStorage(MatrixStorage&& Other)
        : _size1(Other.size1()), _size2(Other.size2()), _data(Other._data) {
        Other._data = nullptr;
    }

    explicit MatrixStorage(std::initializer_list<TDataType> InitialValues)
        : _size1(1), _size2(InitialValues.size()), _data(nullptr) {
        if (_size2 == 0)
            return;
        _data = new TDataType[size()];
        std::size_t position = 0;
        for (auto& i : InitialValues) {
            _data[position++] = i;
        }
    }

    virtual ~MatrixStorage() {
        if (_data)
            delete[] _data;
    }

    template <typename TExpressionType, std::size_t TCategory>
    explicit MatrixStorage(
        MatrixExpression<TExpressionType, TCategory> const& Other)
        : MatrixStorage(Other.expression()) {}

    template <typename TExpressionType>
    explicit MatrixStorage(
        MatrixExpression<TExpressionType, row_major_access> const& Other)
        : _size1(Other.expression().size1()),
          _size2(Other.expression().size2()) {
        _data = new TDataType[size()];
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other.expression()[i];
    }

    template <typename TOtherMatrixType>
    explicit MatrixStorage(TOtherMatrixType const& Other)
        : _size1(Other.size1()), _size2(Other.size2()) {
        _data = new TDataType[size()];
        auto i_data = _data;
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                *(i_data++) = Other(i, j);
    }

    template <typename TExpressionType, std::size_t TCategory>
    MatrixStorage& operator=(
        MatrixExpression<TExpressionType, TCategory> const& Other) {
        auto& other_expression = Other.expression();
        resize(other_expression.size1(), other_expression.size2());
        auto i_data = _data;
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                *(i_data++) = other_expression(i, j);
        return *this;
    }

    template <typename TExpressionType>
    MatrixStorage& operator=(
        MatrixExpression<TExpressionType, row_major_access> const& Other) {
        auto the_expression = Other.expression();
        resize(the_expression.size1(), the_expression.size2());
        auto i_data = _data;
        for (std::size_t i = 0; i < size(); i++)
            *(i_data++) = the_expression[i];
        return *this;
    }

    MatrixStorage& operator=(MatrixStorage const& Other) {
        resize(Other.size1(), Other.size2());
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other._data[i];
        return *this;
    }

    MatrixStorage& operator=(MatrixStorage&& Other) {
        if (_data)
            delete[] _data;

        _size1 = Other.size1();
        _size2 = Other.size2();
        _data = Other._data;
        Other._data = nullptr;

        return *this;
    }

    TDataType& operator()(std::size_t i, std::size_t j) { return at(i, j); }

    TDataType const& operator()(std::size_t i, std::size_t j) const {
        return at(i, j);
    }

    TDataType& at(std::size_t i, std::size_t j) {
        return _data[i * _size2 + j];
    }

    TDataType const& at(std::size_t i, std::size_t j) const {
        return _data[i * _size2 + j];
    }

    TDataType& operator[](std::size_t i) { return at(i); }

    TDataType const& operator[](std::size_t i) const { return at(i); }

    TDataType& at(std::size_t i) { return _data[i]; }

    TDataType const& at(std::size_t i) const { return _data[i]; }

    std::size_t size1() const { return _size1; }

    std::size_t size2() const { return _size2; }

    std::size_t size() const { return _size1 * _size2; }

    void resize(std::size_t NewSize1, std::size_t NewSize2) {
        std::size_t new_size = NewSize1 * NewSize2;
        if (size() != new_size) {
            delete[] _data;
            _data = new TDataType[new_size];
        }
        _size1 = NewSize1;
        _size2 = NewSize2;
    }

    void swap(MatrixStorage& Other){
        TDataType* p_temp = _data;
        _data = Other._data;
        Other._data = p_temp;
    }

    TDataType* data() { return _data; }

    TDataType const* data() const { return _data; }
};

template <typename TDataType, std::size_t TSize1>
class MatrixStorage<TDataType, TSize1, dynamic> {
    std::size_t _size2;
    TDataType* _data;

   public:
    MatrixStorage() = delete;

    MatrixStorage(std::size_t TheSize1, std::size_t TheSize2)
        : _size2(TheSize2) {
        _data = new TDataType[size()];
    }

    explicit MatrixStorage(std::size_t TheSize1, std::size_t TheSize2,
        TDataType const& InitialValue)
        : _size2(TheSize2) {
        _data = new TDataType[size()];
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = InitialValue;
    }

    MatrixStorage(MatrixStorage const& Other) : _size2(Other.size2()) {
        _data = new TDataType[size()];
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other._data[i];
    }

    MatrixStorage(MatrixStorage&& Other)
        : _size2(Other.size2()), _data(Other._data) {
        Other._data = nullptr;
    }

    explicit MatrixStorage(std::initializer_list<TDataType> InitialValues)
        : _size2(InitialValues.size() / TSize1), _data(nullptr) {
        if (_size2 == 0)
            return;
        _data = new TDataType[size()];
        std::size_t position = 0;
        for (auto& i : InitialValues) {
            _data[position++] = i;
        }
    }

    virtual ~MatrixStorage() {
        if (_data)
            delete[] _data;
    }

    template <typename TExpressionType, std::size_t TCategory>
    explicit MatrixStorage(
        MatrixExpression<TExpressionType, TCategory> const& Other)
        : MatrixStorage(Other.expression()) {}

    template <typename TExpressionType>
    explicit MatrixStorage(
        MatrixExpression<TExpressionType, row_major_access> const& Other)
        : _size2(Other.size2()) {
        _data = new TDataType[size()];
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other[i];
    }

    template <typename TOtherMatrixType>
    explicit MatrixStorage(TOtherMatrixType const& Other)
        : _size2(Other.size2()) {
        _data = new TDataType[size()];
        auto i_data = _data;
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                *(i_data++) = Other(i, j);
    }

    template <typename TExpressionType, std::size_t TCategory>
    MatrixStorage& operator=(
        MatrixExpression<TExpressionType, TCategory> const& Other) {
        auto& other_expression = Other.expression();
        std::size_t new_size =
            other_expression.size1() * other_expression.size2();
        if (size() != new_size) {
            delete[] _data;
            _data = new TDataType[new_size];
        }
        _size2 = other_expression.size2();

        auto i_data = _data;
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                *(i_data++) = Other.expression()(i, j);
        return *this;
    }

    template <typename TOtherMatrixType>
    MatrixStorage& operator=(TOtherMatrixType const& Other) {
        std::size_t new_size = Other.size1() * Other.size2();
        if (size() != new_size) {
            delete[] _data;
            _data = new TDataType[new_size];
        }
        _size2 = Other.size2();

        auto i_data = _data;
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                *(i_data++) = Other(i, j);
        return *this;
    }

    MatrixStorage& operator=(MatrixStorage const& Other) {
        std::size_t new_size = Other.size1() * Other.size2();
        if (size() != new_size) {
            delete[] _data;
            _data = new TDataType[new_size];
        }
        _size2 = Other.size2();

        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other._data[i];
        return *this;
    }

    MatrixStorage& operator=(MatrixStorage&& Other) {
        if (_data)
            delete[] _data;

        _size2 = Other.size2();
        _data = Other._data;
        Other._data = nullptr;

        return *this;
    }

    TDataType& operator()(std::size_t i, std::size_t j) { return at(i, j); }

    TDataType const& operator()(std::size_t i, std::size_t j) const {
        return at(i, j);
    }

    TDataType& at(std::size_t i, std::size_t j) {
        return _data[i * _size2 + j];
    }

    TDataType const& at(std::size_t i, std::size_t j) const {
        return _data[i * _size2 + j];
    }

    TDataType& operator[](std::size_t i) { return at(i); }

    TDataType const& operator[](std::size_t i) const { return at(i); }

    TDataType& at(std::size_t i) { return _data[i]; }

    TDataType const& at(std::size_t i) const { return _data[i]; }

    constexpr std::size_t size1() const { return TSize1; }

    std::size_t size2() const { return _size2; }

    std::size_t size() const { return TSize1 * _size2; }

    void resize(std::size_t NewSize) {
        if (size2() != NewSize) {
            delete[] _data;
            _data = new TDataType[TSize1 * NewSize];
        }
        _size2 = NewSize;
    }

    void swap(MatrixStorage& Other){
        TDataType* p_temp = _data;
        _data = Other._data;
        Other._data = p_temp;
    }

    TDataType* data() { return _data; }

    TDataType const* data() const { return _data; }
};

template <typename TDataType, std::size_t TSize2>
class MatrixStorage<TDataType, dynamic, TSize2> {
    std::size_t _size1;
    TDataType* _data;

   public:
    MatrixStorage() : _size1(0), _data(nullptr) {}

    MatrixStorage(std::size_t TheSize1, std::size_t TheSize2 = TSize2)
        : _size1(TheSize1) {
        _data = new TDataType[size()];
    }

    explicit MatrixStorage(std::size_t TheSize1, std::size_t TheSize2,
        TDataType const& InitialValue)
        : _size1(TheSize1) {
        _data = new TDataType[size()];
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = InitialValue;
    }

    MatrixStorage(MatrixStorage const& Other) : _size1(Other.size1()) {
        _data = new TDataType[size()];
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other._data[i];
    }

    MatrixStorage(MatrixStorage&& Other)
        : _size1(Other.size1()), _data(Other._data) {
        Other._data = nullptr;
    }

    explicit MatrixStorage(std::initializer_list<TDataType> InitialValues)
        : _size1(InitialValues.size() / TSize2), _data(nullptr) {
        if (_size1 == 0)
            return;
        _data = new TDataType[size()];
        std::size_t position = 0;
        for (auto& i : InitialValues) {
            _data[position++] = i;
        }
    }

    virtual ~MatrixStorage() {
        if (_data)
            delete[] _data;
    }

    template <typename TExpressionType, std::size_t TCategory>
    explicit MatrixStorage(
        MatrixExpression<TExpressionType, TCategory> const& Other)
        : MatrixStorage(Other.expression()) {}

    template <typename TExpressionType>
    explicit MatrixStorage(
        MatrixExpression<TExpressionType, row_major_access> const& Other)
        : _size1(Other.expression().size1()) {
        _data = new TDataType[size()];
        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other.expression()[i];
    }

    template <typename TOtherMatrixType>
    explicit MatrixStorage(TOtherMatrixType const& Other)
        : _size1(Other.size1()) {
        _data = new TDataType[size()];
        auto i_data = _data;
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                *(i_data++) = Other(i, j);
    }

    template <typename TExpressionType, std::size_t TCategory>
    MatrixStorage& operator=(
        MatrixExpression<TExpressionType, TCategory> const& Other) {
        auto& other_expression = Other.expression();
        std::size_t new_size =
            other_expression.size1() * other_expression.size2();
        if (size() != new_size) {
            delete[] _data;
            _data = new TDataType[new_size];
        }
        _size1 = other_expression.size1();

        auto i_data = _data;
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                *(i_data++) = Other.expression()(i, j);
        return *this;
    }

    template <typename TOtherMatrixType>
    MatrixStorage& operator=(TOtherMatrixType const& Other) {
        std::size_t new_size = Other.size();
        if (size() != new_size) {
            delete[] _data;
            _data = new TDataType[new_size];
        }
        _size1 = Other.size1();

        auto i_data = _data;
        for (std::size_t i = 0; i < size1(); i++)
            for (std::size_t j = 0; j < size2(); j++)
                *(i_data++) = Other(i, j);
        return *this;
    }

    MatrixStorage& operator=(MatrixStorage const& Other) {
        std::size_t new_size = Other.size1() * Other.size2();
        if (size() != new_size) {
            delete[] _data;
            _data = new TDataType[new_size];
        }
        _size1 = Other.size1();

        for (std::size_t i = 0; i < size(); i++)
            _data[i] = Other._data[i];
        return *this;
    }

    MatrixStorage& operator=(MatrixStorage&& Other) {
        if (_data)
            delete[] _data;

        _size1 = Other.size1();
        _data = Other._data;
        Other._data = nullptr;

        return *this;
    }

    TDataType& operator()(std::size_t i, std::size_t j) { return at(i, j); }

    TDataType const& operator()(std::size_t i, std::size_t j) const {
        return at(i, j);
    }

    TDataType& at(std::size_t i, std::size_t j) {
        return _data[i * TSize2 + j];
    }

    TDataType const& at(std::size_t i, std::size_t j) const {
        return _data[i * TSize2 + j];
    }

    TDataType& operator[](std::size_t i) { return at(i); }

    TDataType const& operator[](std::size_t i) const { return at(i); }

    TDataType& at(std::size_t i) { return _data[i]; }

    TDataType const& at(std::size_t i) const { return _data[i]; }

    std::size_t size1() const { return _size1; }

    constexpr std::size_t size2() const { return TSize2; }

    std::size_t size() const { return _size1 * TSize2; }

    void resize(std::size_t NewSize) {
        if (size1() != NewSize) {
            delete[] _data;
            _data = new TDataType[NewSize * TSize2];
        }
        _size1 = NewSize;
    }

    void swap(MatrixStorage& Other){
        TDataType* p_temp = _data;
        _data = Other._data;
        Other._data = p_temp;
    }

    TDataType* data() { return _data; }

    TDataType const* data() const { return _data; }
};

}  // namespace AMatrix