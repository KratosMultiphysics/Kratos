#ifndef AMGCL_ADAPTER_EPETRA_HPP
#define AMGCL_ADAPTER_EPETRA_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>
Copyright (c) 2014, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
\file    amgcl/adapter/epetra.hpp
\author  Denis Demidov <dennis.demidov@gmail.com>
\brief   Adapt Epetra_CrsMatrix from Trilinos.
\ingroup adapters
*/

#include <vector>

#include <Epetra_CrsMatrix.h>
#include <Epetra_IntVector.h>
#include <Epetra_Import.h>

#include <amgcl/backend/interface.hpp>

namespace amgcl {
namespace adapter {

/// Adapts Epetra_CrsMatrix
class epetra_map {
    public:
        epetra_map(const Epetra_CrsMatrix &A)
            : A(A), order(A.ColMap())
        {
            const Epetra_Map& row_map = A.RowMap();
            const Epetra_Map& col_map = A.ColMap();

            int entries_before;
            int local_entries = row_map.NumMyElements();
            A.Comm().ScanSum(&local_entries, &entries_before, 1);
            entries_before -= local_entries;

            Epetra_IntVector perm(row_map);
            for(int i = 0, j = entries_before; i < local_entries; ++i, ++j)
                perm[i] = j;

            Epetra_Import importer = Epetra_Import(col_map, row_map);

            order.Import(perm, importer, Insert);
        }

        size_t rows() const {
            return A.NumMyRows();
        }

        size_t cols() const {
            return A.NumGlobalCols();
        }

        size_t nonzeros() const {
            return A.NumMyNonzeros();
        }

        class row_iterator {
            public:
                typedef int    col_type;
                typedef double val_type;

                row_iterator(
                        const Epetra_CrsMatrix &A,
                        const Epetra_IntVector &order,
                        int row
                        ) : order(order)
                {
                    int nnz;
                    A.ExtractMyRowView(row, nnz, m_val, m_col);
                    m_end = m_col + nnz;
                }

                operator bool() const {
                    return m_col != m_end;
                }

                row_iterator& operator++() {
                    ++m_col;
                    ++m_val;
                    return *this;
                }

                col_type col() const {
                    return order[*m_col];
                }

                val_type value() const {
                    return *m_val;
                }

            private:
                const Epetra_IntVector &order;

                col_type * m_col;
                col_type * m_end;
                val_type * m_val;
        };

        row_iterator row_begin(int row) const {
            return row_iterator(A, order, row);
        }
    private:
        const Epetra_CrsMatrix &A;
        Epetra_IntVector order;
};

/// Adapts Epetra_CrsMatrix
inline epetra_map map(const Epetra_CrsMatrix &A) {
    return epetra_map(A);
}

} // namespace adapter

//---------------------------------------------------------------------------
// Specialization of matrix interface
//---------------------------------------------------------------------------
namespace backend {

template <>
struct value_type < adapter::epetra_map > {
    typedef double type;
};

template <>
struct rows_impl < adapter::epetra_map > {
    static size_t get(const adapter::epetra_map &A) {
        return A.rows();
    }
};

template <>
struct cols_impl < adapter::epetra_map > {
    static size_t get(const adapter::epetra_map &A) {
        return A.cols();
    }
};

template <>
struct nonzeros_impl < adapter::epetra_map > {
    static size_t get(const adapter::epetra_map &A) {
        return A.nonzeros();
    }
};

template <>
struct row_iterator < adapter::epetra_map > {
    typedef adapter::epetra_map::row_iterator type;
};

template <>
struct row_begin_impl< adapter::epetra_map >
{
    static typename row_iterator< adapter::epetra_map >::type
    get(const adapter::epetra_map &A, size_t row) {
        return A.row_begin(row);
    }
};

} // namespace backend
} // namespace amgcl


#endif
