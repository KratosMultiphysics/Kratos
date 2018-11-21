#pragma once

#include "Math.h"

#include <vector>

namespace ANurbs {

template <typename TData>
class Grid
{
public:
    using DataType = TData;

private:
    int m_nbRows;
    int m_nbCols;
    std::vector<DataType> m_values;

public:
    Grid(
        const int nbRows,
        const int nbCols)
        : m_nbRows(nbRows)
        , m_nbCols(nbCols)
        , m_values(nbRows * nbCols)
    {
    }

    Grid(
        const int nbRows,
        const int nbCols,
        const std::vector<DataType> values)
        : m_nbRows(nbRows)
        , m_nbCols(nbCols)
        , m_values(values)
    {
    }

    int
    NbRows() const
    {
        return m_nbRows;
    }

    int
    NbCols() const
    {
        return m_nbCols;
    }

    int
    NbValues() const
    {
        return NbRows() * NbCols();
    }

    const DataType&
    operator[](
        const int index) const
    {
        return m_values[index];
    }

    DataType&
    operator[](
        const int index)
    {
        return m_values[index];
    }

    const DataType&
    operator()(
        const int row,
        const int col) const
    {
        int index = Math::MatrixIndex(NbRows(), NbCols(), row, col);

        return m_values[index];
    }

    void
    SetValue(
        const int index,
        const DataType value)
    {
        m_values[index] = value;
    }

    void
    SetValue(
        const int row,
        const int col,
        const DataType& value)
    {
        int index = Math::MatrixIndex(NbRows(), NbCols(), row, col);

        m_values[index] = value;
    }

    void
    Resize(
        const int nbRows,
        const int nbCols)
    {
        m_nbRows = nbRows;
        m_nbCols = nbCols;
        m_values.resize(nbRows * nbCols);
    }

    void
    Clear()
    {
        m_values.clear();
        m_values.resize(NbValues());
    }
};

} // namespace ANurbs
