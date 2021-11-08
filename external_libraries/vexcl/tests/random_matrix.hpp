#ifndef TESTS_RANDOM_MATRIX_HPP
#define TESTS_RANDOM_MATRIX_HPP

#include <vector>
#include <set>
#include "random_vector.hpp"

template <typename RT, typename CT, typename VT>
void random_matrix(size_t n, size_t m, size_t nnz_per_row,
        std::vector<RT> &row,
        std::vector<CT> &col,
        std::vector<VT> &val
        )
{
    row.clear();
    col.clear();

    row.reserve(n + 1);
    col.reserve(nnz_per_row * n);

    std::default_random_engine rng( std::rand() );
    std::uniform_int_distribution<size_t> random_width(0, nnz_per_row - 1);
    std::uniform_int_distribution<size_t> random_column(0, m - 1);

    row.push_back(0);
    for(size_t k = 0; k < n; k++) {
        size_t width = random_width(rng);

        std::set<CT> cs;
        while(cs.size() < width)
            cs.insert(static_cast<CT>(random_column(rng)));

        for(auto c = cs.begin(); c != cs.end(); c++)
            col.push_back(*c);

        row.push_back(static_cast<RT>(col.size()));
    }

    random_vector<VT>( col.size() ).swap(val);
}

#endif
