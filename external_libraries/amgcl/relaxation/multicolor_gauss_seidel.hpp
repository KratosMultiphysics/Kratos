#ifndef AMGCL_RELAXATION_MULTICOLOR_GAUSS_SEIDEL_HPP
#define AMGCL_RELAXATION_MULTICOLOR_GAUSS_SEIDEL_HPP

/*
The MIT License

Copyright (c) 2012-2017 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/relaxation/multicolor_gauss_seidel.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Multicolor Gauss-Seidel relaxation scheme.
 */

#include <algorithm>
#include <numeric>

#include <boost/graph/sequential_vertex_coloring.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <amgcl/backend/interface.hpp>
#include <amgcl/util.hpp>

namespace amgcl {

namespace detail {

//---------------------------------------------------------------------------
// Adapt backend::crs<> for use with Boost Graph Library.
//---------------------------------------------------------------------------
template <class Matrix>
struct graph {
    const Matrix &A;

    typedef typename Matrix::col_type vertex_descriptor;
    typedef typename Matrix::ptr_type edge_descriptor;
    typedef boost::directed_tag directed_category;
    typedef boost::allow_parallel_edge_tag edge_parallel_category;

    class traversal_category:
        public boost::incidence_graph_tag,
        public boost::adjacency_graph_tag,
        public boost::vertex_list_graph_tag,
        public boost::edge_list_graph_tag
    {};

    typedef boost::counting_iterator<vertex_descriptor> vertex_iterator;
    typedef vertex_descriptor vertices_size_type;

    typedef const vertex_descriptor* adjacency_iterator;

    static vertex_descriptor null_vertex() { return -1; }

    graph(const Matrix &A) : A(A) {}
};

template <class Matrix>
std::pair<
    typename graph<Matrix>::vertex_iterator,
    typename graph<Matrix>::vertex_iterator
    >
vertices(const graph<Matrix> &G) {
    return std::make_pair(
            boost::counting_iterator<typename graph<Matrix>::vertex_descriptor>(0),
            boost::counting_iterator<typename graph<Matrix>::vertex_descriptor>(backend::rows(G.A))
            );
}

template <class Matrix>
std::pair<
    typename graph<Matrix>::adjacency_iterator,
    typename graph<Matrix>::adjacency_iterator
    >
adjacent_vertices(ptrdiff_t v, const graph<Matrix> &G) {
    typename Matrix::ptr_type row_beg = G.A.ptr[v];
    typename Matrix::ptr_type row_end = G.A.ptr[v + 1];

    return std::make_pair(
            G.A.col + row_beg,
            G.A.col + row_end
            );
}

template <class Matrix>
typename graph<Matrix>::vertices_size_type
num_vertices(const graph<Matrix> &G) {
    return backend::rows(G.A);
}

template <class Matrix>
graph<Matrix> as_graph(const Matrix &A) {
    return graph<Matrix>(A);
}

} // namespace detail

namespace relaxation {

template <class Backend>
struct multicolor_gauss_seidel {
    typedef amgcl::detail::empty_params params;

    template <class Matrix>
    multicolor_gauss_seidel(
            const Matrix &A,
            const params&,
            const typename Backend::params&
            ) : order( backend::rows(A) )
    {
        const size_t n = backend::rows(A);

        std::vector<int> color(n);
        num_colors = boost::sequential_vertex_coloring(amgcl::detail::as_graph(A), &color[0]);

        ptr.resize(num_colors + 1, 0);

        for(size_t i = 0; i < n; ++i) {
            ++ptr[ color[i] + 1];
            order[i] = i;
        }

        std::stable_sort(order.begin(), order.end(), order_by(color));
        std::partial_sum(ptr.begin(), ptr.end(), ptr.begin());
    }

    template <class Matrix, class VecRHS, class VecX, class VecTMP>
    void apply_pre(const Matrix &A, const VecRHS &rhs, VecX &x, VecTMP&, const params&) const
    {
        for(int c = 0; c < num_colors; ++c) iterate(A, rhs, x, c);
    }

    template <class Matrix, class VecRHS, class VecX, class VecTMP>
    void apply_post(const Matrix &A, const VecRHS &rhs, VecX &x, VecTMP&, const params&) const
    {
        for(int c = num_colors; c --> 0; ) iterate(A, rhs, x, c);
    }

    template <class Matrix, class VecRHS, class VecX>
    void apply(const Matrix &A, const VecRHS &rhs, VecX &x, const params&) const
    {
        backend::clear(x);
        for(int c = 0; c < num_colors; ++c) iterate(A, rhs, x, c);
        for(int c = num_colors; c --> 0; )  iterate(A, rhs, x, c);
    }

    private:
        struct order_by {
            const std::vector<int> &ref;

            order_by(const std::vector<int> &ref) : ref(ref) {}

            bool operator()(ptrdiff_t i, ptrdiff_t j) const {
                return ref[i] < ref[j];
            }
        };

        int num_colors;
        std::vector<ptrdiff_t> order, ptr;

        template <class Matrix, class VectorRHS, class VectorX>
        void iterate(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, int c
                ) const
        {
            typedef typename backend::row_iterator<Matrix>::type row_iterator;
            typedef typename backend::value_type<Matrix>::type    val_type;
            typedef typename backend::value_type<VectorRHS>::type rhs_type;

#pragma omp parallel for
            for(ptrdiff_t j = ptr[c]; j < ptr[c+1]; ++j) {
                ptrdiff_t i = order[j];

                rhs_type temp = rhs[i];
                val_type diag = math::identity<val_type>();
                for (row_iterator a = backend::row_begin(A, i); a; ++a) {
                    if (a.col() == i)
                        diag = a.value();
                    else
                        temp -= a.value() * x[a.col()];
                }

                x[i] = math::inverse(diag) * temp;
            }
        }
};

} // namespace relaxation

namespace backend {

template <class Backend>
struct relaxation_is_supported<
    Backend,
    relaxation::multicolor_gauss_seidel,
    typename boost::disable_if<
            typename Backend::provides_row_iterator
        >::type
    > : boost::false_type
{};

} // namespace backend
} // namespace amgcl


#endif
