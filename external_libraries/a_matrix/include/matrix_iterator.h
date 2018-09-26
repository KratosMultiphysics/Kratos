#include <iterator>

namespace AMatrix {

template <typename TDataType>
class RandomAccessIterator
    : public std::iterator<std::random_access_iterator_tag, TDataType> {
    TDataType* _p_data;

   public:
    using base_type = std::iterator<std::random_access_iterator_tag, TDataType>;
    using typename base_type::difference_type;
    RandomAccessIterator() = default;
    RandomAccessIterator(RandomAccessIterator const& Other) = default;
    RandomAccessIterator(RandomAccessIterator&& Other) = default;
    RandomAccessIterator(TDataType* pData) : _p_data(pData) {}

    template <typename TOtherDataType>
    RandomAccessIterator(RandomAccessIterator<TOtherDataType> const& Other):
        _p_data(Other.operator->()) {}

    RandomAccessIterator& operator=(
        RandomAccessIterator const& Other) = default;

    bool operator==(const RandomAccessIterator& Other) const {
        return _p_data == Other._p_data;
    }

    bool operator!=(const RandomAccessIterator& Other) const {
        return _p_data != Other._p_data;
    }

    TDataType const& operator*() const { return *_p_data; }

    TDataType& operator*() { return *_p_data; }

    TDataType* const operator->() const { return _p_data; }

    TDataType* operator->() { return _p_data; }

    RandomAccessIterator& operator++() {
        ++_p_data;
        return *this;
    }

    RandomAccessIterator operator++(int) {
        RandomAccessIterator temp(*this);
        ++_p_data;
        return temp;
    }

    RandomAccessIterator& operator--() {
        --_p_data;
        return *this;
    }

    RandomAccessIterator operator--(int) {
        RandomAccessIterator temp(*this);
        --_p_data;
        return temp;
    }

    RandomAccessIterator& operator+=(difference_type Offset) {
        _p_data += Offset;
        return *this;
    }

    RandomAccessIterator& operator-=(difference_type Offset) {
        _p_data -= Offset;
        return *this;
    }

    friend RandomAccessIterator operator+(
        RandomAccessIterator First, difference_type Second) {
        First += Second;
        return First;
    }

    friend RandomAccessIterator operator+(
        difference_type First, RandomAccessIterator Second) {
        Second += First;
        return Second;
    }

    friend RandomAccessIterator operator-(
        RandomAccessIterator First, difference_type Second) {
        First -= Second;
        return First;
    }

    friend difference_type operator-(
        RandomAccessIterator const& First, RandomAccessIterator const& Second) {
        return First._p_data - Second._p_data;
    }

    friend bool operator<(
        RandomAccessIterator const& First, RandomAccessIterator const& Second) {
        return First._p_data < Second._p_data;
    }

    friend bool operator>(
        RandomAccessIterator const& First, RandomAccessIterator const& Second) {
        return Second < First;
    }
    friend bool operator<=(
        RandomAccessIterator const& First, RandomAccessIterator const& Second) {
        return !(First > Second);
    }
    friend bool operator>=(
        RandomAccessIterator const& First, RandomAccessIterator const& Second) {
        return !(First < Second);
    }
};

}  // namespace AMatrix