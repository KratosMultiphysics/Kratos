#ifndef AMGCL_BACKEND_HPX_HPP
#define AMGCL_BACKEND_HPX_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/backend/hpx.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  HPX backend.
 */

#include <vector>

#include <hpx/hpx.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/parallel_for_each.hpp>
#include <hpx/include/parallel_transform_reduce.hpp>

#include <boost/range/irange.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/iterator_range.hpp>

#include <amgcl/util.hpp>
#include <amgcl/backend/builtin.hpp>

namespace amgcl {
namespace backend {

/// The matrix is a thin wrapper on top of amgcl::builtin::crs<>.
template <typename real>
class hpx_matrix {
    public:
        typedef real      value_type;
        typedef ptrdiff_t index_type;

        typedef crs<value_type, index_type> Base;
        typedef typename Base::row_iterator row_iterator;

        // For each of the output segments y[i] in y = A * x it stores a range of
        // segments in x that y[i] depends on.
        std::vector<std::tuple<index_type, index_type>> xrange;

        // And yrange stores the inverted dependencies: for each segment in x
        // yrange stores range of segments in y that depend on x. This is required
        // to determine when a segment of x in [y = A * x] is safe to update.
        std::vector<std::tuple<index_type, index_type>> yrange;

        // Creates the matrix from builtin datatype, sets up xrange.
        hpx_matrix(std::shared_ptr<Base> A, int grain_size) : base(A)
        {
            index_type n = backend::rows(*A);
            index_type m = backend::cols(*A);

            index_type nseg = (n + grain_size - 1) / grain_size;
            index_type mseg = (m + grain_size - 1) / grain_size;

            xrange.resize(nseg);
            yrange.resize(mseg, std::make_tuple(m, 0));

            auto range = boost::irange<index_type>(0, nseg);

            hpx::parallel::for_each(
                    hpx::parallel::par,
                    boost::begin(range), boost::end(range),
                    [this, grain_size, n, A](index_type seg) {
                        index_type i = seg * grain_size;
                        index_type beg = A->ptr[i];
                        index_type end = A->ptr[std::min<index_type>(i + grain_size, n)];

                        auto mm = std::minmax_element(A->col + beg, A->col + end);

                        index_type xbeg = *std::get<0>(mm) / grain_size;
                        index_type xend = *std::get<1>(mm) / grain_size + 1;

                        xrange[seg] = std::make_tuple(xbeg, xend);

                        for(index_type i = xbeg; i < xend; ++i) {
                            std::get<0>(yrange[i]) = std::min(seg,   std::get<0>(yrange[i]));
                            std::get<1>(yrange[i]) = std::max(seg+1, std::get<1>(yrange[i]));
                        }
                    });
        }

        size_t rows()     const { return backend::rows(*base);     }
        size_t cols()     const { return backend::cols(*base);     }
        size_t nonzeros() const { return backend::nonzeros(*base); }

        row_iterator row_begin(size_t row) const {
            return base->row_begin(row);
        }
    private:
        // Base matrix is stored in shared_ptr<> to reduce the overhead
        // of data transfer from builtin datatypes (used for AMG setup) to the
        // backend datatypes.
        std::shared_ptr<Base> base;

};

/// The vector type to be used with HPX backend.
/**
 * hpx_vector is a thin wrapper on top of std::vector.
 * The vector consists of continuous segments of fixed size (grain_size) except
 * may be the last one that is allowed to be shorter.
 * A vector of shared_futures corresponding to each of the segments is stored
 * along the data vector to facilitate construction of HPX dependency graph.
 */
template < typename real >
class hpx_vector {
    public:
        typedef real                          value_type;
        typedef std::vector<real>             Base;
        typedef typename Base::iterator       iterator;
        typedef typename Base::const_iterator const_iterator;

        int nseg;        // Number of segments in the vector
        int grain_size;  // Segment size.

        // Futures associated with each segment:
        mutable std::vector<hpx::shared_future<void>> safe_to_read;
        mutable std::vector<hpx::shared_future<void>> safe_to_write;

        hpx_vector(size_t n, int grain_size)
            : nseg( (n + grain_size - 1) / grain_size ),
              grain_size( grain_size ),
              buf( std::make_shared<Base>(n) )
        {
            precondition(grain_size > 0, "grain size should be positive");
            init_futures();
        }

        template <class Other>
        hpx_vector(std::shared_ptr<Other> o, int grain_size)
            : nseg( (o->size() + grain_size - 1) / grain_size ),
              grain_size( grain_size ),
              buf(std::make_shared<Base>(o->data(), o->data() + o->size()))
        {
            precondition(grain_size > 0, "grain size should be positive");
            init_futures();
        }

        size_t size() const { return buf->size(); }

        const real & operator[](size_t i) const { return (*buf)[i]; }
        real & operator[](size_t i) { return (*buf)[i]; }

        const real* data() const { return buf->data(); }
        real*       data()       { return buf->data(); }

        iterator begin() { return buf->begin(); }
        iterator end()   { return buf->end();   }

        const_iterator begin() const { return buf->cbegin(); }
        const_iterator end()   const { return buf->cend();   }

        const_iterator cbegin() const { return buf->cbegin(); }
        const_iterator cend()   const { return buf->cend();   }

        template <class IdxTuple>
        boost::iterator_range<
            typename std::vector< hpx::shared_future<void> >::iterator
            >
        safe_range(IdxTuple idx) const {
            return boost::make_iterator_range(
                    safe_to_read.begin() + std::get<0>(idx),
                    safe_to_read.begin() + std::get<1>(idx)
                    );
        }
    private:
        // Segments stored in a continuous array.
        // The base vector is stored with shared_ptr for the same reason as with
        // hpx_matrix above: to reduce the overhead of data transfer.
        std::shared_ptr<Base> buf;

        void init_futures() {
            safe_to_read.reserve(nseg);
            safe_to_write.reserve(nseg);
            for(ptrdiff_t i = 0; i < nseg; ++i) {
                safe_to_read.push_back(hpx::make_ready_future());
                safe_to_write.push_back(hpx::make_ready_future());
            }
        }
};

/// HPX backend
/**
 * This is a backend that is based on HPX -- a general purpose C++ runtime
 * system for parallel and distributed applications of any scale
 * http://stellar-group.org/libraries/hpx.
 */
template <typename real>
struct HPX {
    typedef real      value_type;
    typedef ptrdiff_t index_type;

    struct provides_row_iterator : std::false_type {};

    struct params {
        /// Number of vector elements in a single segment.
        int grain_size;

        params() : grain_size(4096) {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, grain_size)
        {
            check_params(p, {"grain_size"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, grain_size);
        }
#endif
    };

    typedef hpx_matrix<value_type>         matrix;
    typedef hpx_vector<value_type>         vector;
    typedef hpx_vector<value_type>         matrix_diagonal;

    struct direct_solver : public solver::skyline_lu<value_type> {
        typedef solver::skyline_lu<value_type> Base;
        typedef typename Base::params params;

        template <class Matrix>
        direct_solver(const Matrix &A, const params &prm = params())
            : Base(A, prm)
        {}

        struct call_base {
            const Base *base;
            const real *fptr;
            real       *xptr;

            template <class... T>
            void operator()(T&&...) const {
                (*base)(fptr, xptr);
            }
        };

        void operator()(const vector &rhs, vector &x) const {
            const real *fptr = &rhs[0];
            real       *xptr = &x[0];

            using hpx::dataflow;

            hpx::shared_future<void> solve = dataflow(
                    hpx::launch::async,
                    call_base{this, fptr, xptr},
                    rhs.safe_to_read,
                    x.safe_to_write
                    );

            boost::fill(x.safe_to_read, solve);
        }
    };

    static std::string name() { return "HPX"; }

    /// Copy matrix.
    static std::shared_ptr<matrix>
    copy_matrix(std::shared_ptr<typename matrix::Base> A, const params &p)
    {
        return std::make_shared<matrix>(A, p.grain_size);
    }

    /// Copy vector to builtin backend.
    static std::shared_ptr<vector>
    copy_vector(const typename vector::Base &x, const params &p)
    {
        return std::make_shared<vector>(
                std::make_shared<typename vector::Base>(x), p.grain_size
                );
    }

    /// Copy vector to builtin backend.
    template <typename Other>
    static std::shared_ptr<hpx_vector<typename Other::value_type>>
    copy_vector(std::shared_ptr<Other> x, const params &p)
    {
        return std::make_shared<hpx_vector<typename Other::value_type>>(x, p.grain_size);
    }

    /// Create vector of the specified size.
    static std::shared_ptr<vector>
    create_vector(size_t size, const params &p)
    {
        return std::make_shared<vector>(size, p.grain_size);
    }

    /// Create direct solver for coarse level
    static std::shared_ptr<direct_solver>
    create_solver(std::shared_ptr<typename matrix::Base> A, const params&) {
        return std::make_shared<direct_solver>(*A);
    }
};

//---------------------------------------------------------------------------
// Backend interface implementation
//---------------------------------------------------------------------------
template < typename real >
struct rows_impl< hpx_matrix<real> > {
    static size_t get(const hpx_matrix<real> &A) {
        return A.rows();
    }
};

template < typename real >
struct cols_impl< hpx_matrix<real> > {
    static size_t get(const hpx_matrix<real> &A) {
        return A.cols();
    }
};

template < typename real >
struct nonzeros_impl< hpx_matrix<real> > {
    static size_t get(const hpx_matrix<real> &A) {
        return A.nonzeros();
    }
};

template < typename Alpha, typename Beta, typename real >
struct spmv_impl<
    Alpha, hpx_matrix<real>, hpx_vector<real>,
    Beta,  hpx_vector<real>
    >
{
    typedef hpx_matrix<real> matrix;
    typedef hpx_vector<real> vector;

    struct process_ab {
        Alpha                   alpha;
        const hpx_matrix<real> &A;
        const real             *xptr;
        Beta                    beta;
        real                   *yptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        void operator()(T&&...) const {
            for(ptrdiff_t i = beg; i < end; ++i) {
                real sum = 0;
                for(auto a = A.row_begin(i); a; ++a)
                    sum += a.value() * xptr[a.col()];
                yptr[i] = alpha * sum + beta * yptr[i];
            }
        }
    };

    struct process_a {
        Alpha                   alpha;
        const hpx_matrix<real> &A;
        const real             *xptr;
        real                   *yptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        void operator()(T&&...) const {
            for(ptrdiff_t i = beg; i < end; ++i) {
                real sum = 0;
                for(auto a = A.row_begin(i); a; ++a)
                    sum += a.value() * xptr[a.col()];
                yptr[i] = alpha * sum;
            }
        }
    };

    struct wait_for_it {
        template <class T>
        void operator()(T&&) const {}
    };

    static void apply(Alpha alpha, const matrix &A, const vector &x,
            Beta beta, vector &y)
    {
        const real *xptr = &x[0];
        real       *yptr = &y[0];

        using hpx::dataflow;

        auto range = boost::irange(0, y.nseg);

        if (beta) {
            // y = alpha * A * x + beta * y
            hpx::parallel::for_each(
                    hpx::parallel::par,
                    boost::begin(range), boost::end(range),
                    [alpha, &A, &x, beta, &y, xptr, yptr](ptrdiff_t seg) {
                        ptrdiff_t beg = seg * y.grain_size;
                        ptrdiff_t end = std::min<ptrdiff_t>(beg + y.grain_size, y.size());

                        y.safe_to_read[seg] = dataflow(
                                hpx::launch::async,
                                process_ab{alpha, A, xptr, beta, yptr, beg, end},
                                y.safe_to_read[seg],
                                y.safe_to_write[seg],
                                x.safe_range(A.xrange[seg])
                                );
                    });
        } else {
            // y = alpha * A * x
            hpx::parallel::for_each(
                    hpx::parallel::par,
                    boost::begin(range), boost::end(range),
                    [alpha, &A, &x, &y, xptr, yptr](ptrdiff_t seg) {
                        ptrdiff_t beg = seg * y.grain_size;
                        ptrdiff_t end = std::min<ptrdiff_t>(beg + y.grain_size, y.size());

                        y.safe_to_read[seg] = dataflow(hpx::launch::async,
                                process_a{alpha, A, xptr, yptr, beg, end},
                                y.safe_to_write[seg],
                                x.safe_range(A.xrange[seg])
                                );
                    });
        }

        // Do not update x until y is ready.
        range = boost::irange(0, x.nseg);
        hpx::parallel::for_each(
                hpx::parallel::par,
                boost::begin(range), boost::end(range),
                [&A, &x, &y](ptrdiff_t seg) {
                    x.safe_to_write[seg] = dataflow(hpx::launch::async,
                            wait_for_it(),
                            y.safe_range(A.yrange[seg])
                            );
                });
    }
};

template < typename real >
struct residual_impl<
    hpx_matrix<real>,
    hpx_vector<real>,
    hpx_vector<real>,
    hpx_vector<real>
    >
{
    typedef hpx_matrix<real> matrix;
    typedef hpx_vector<real> vector;

    struct process {
        const real             *fptr;
        const hpx_matrix<real> &A;
        const real             *xptr;
        real                   *rptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        void operator()(T&&...) const {
            for(ptrdiff_t i = beg; i < end; ++i) {
                real sum = fptr[i];
                for(auto a = A.row_begin(i); a; ++a)
                    sum -= a.value() * xptr[a.col()];
                rptr[i] = sum;
            }
        }
    };

    struct wait_for_it {
        template <class T>
        void operator()(T&&) const {}
    };

    static void apply(const vector &f, const matrix &A, const vector &x,
            vector &r)
    {
        const real *xptr = &x[0];
        const real *fptr = &f[0];
        real       *rptr = &r[0];

        using hpx::dataflow;

        auto range = boost::irange(0, f.nseg);
        hpx::parallel::for_each(
                hpx::parallel::par,
                boost::begin(range), boost::end(range),
                [&f, &A, &x, &r, xptr, fptr, rptr](ptrdiff_t seg) {
                    ptrdiff_t beg = seg * f.grain_size;
                    ptrdiff_t end = std::min<ptrdiff_t>(beg + f.grain_size, f.size());

                    r.safe_to_read[seg] = dataflow(hpx::launch::async,
                            process{fptr, A, xptr, rptr, beg, end},
                            f.safe_to_read[seg],
                            r.safe_to_write[seg],
                            x.safe_range(A.xrange[seg])
                            );
                });

        // Do not update x until r is ready.
        range = boost::irange(0, x.nseg);
        hpx::parallel::for_each(
                hpx::parallel::par,
                boost::begin(range), boost::end(range),
                [&A, &x, &r](ptrdiff_t seg) {
                    x.safe_to_write[seg] = dataflow(hpx::launch::async,
                            wait_for_it(),
                            r.safe_range(A.yrange[seg])
                            );
                });
    }
};

template < typename real >
struct clear_impl<
    hpx_vector<real>
    >
{
    typedef hpx_vector<real> vector;

    struct process {
        real *xptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        void operator()(T&&...) const {
            for(ptrdiff_t i = beg; i < end; ++i)
                xptr[i] = 0;
        }
    };

    static void apply(vector &x) {
        real *xptr = &x[0];

        using hpx::dataflow;

        auto range = boost::irange(0, x.nseg);
        hpx::parallel::for_each(
                hpx::parallel::par,
                boost::begin(range), boost::end(range),
                [&x, xptr](ptrdiff_t seg) {
                    ptrdiff_t beg = seg * x.grain_size;
                    ptrdiff_t end = std::min<ptrdiff_t>(beg + x.grain_size, x.size());

                    x.safe_to_read[seg] = dataflow(hpx::launch::async,
                            process{xptr, beg, end},
                            x.safe_to_write[seg]
                            );
                });
    }
};

template < typename real >
struct copy_impl<
    hpx_vector<real>,
    hpx_vector<real>
    >
{
    typedef hpx_vector<real> vector;

    struct process {
        const real *xptr;
        real       *yptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        void operator()(T&&...) const {
            for(ptrdiff_t i = beg; i < end; ++i)
                yptr[i] = xptr[i];
        }
    };

    static void apply(const vector &x, vector &y)
    {
        const real *xptr = &x[0];
        real       *yptr = &y[0];

        using hpx::dataflow;

        auto range = boost::irange(0, x.nseg);
        hpx::parallel::for_each(
                hpx::parallel::par,
                boost::begin(range), boost::end(range),
                [&x, &y, xptr, yptr](ptrdiff_t seg) {
                    ptrdiff_t beg = seg * x.grain_size;
                    ptrdiff_t end = std::min<ptrdiff_t>(beg + x.grain_size, x.size());

                    y.safe_to_read[seg] = dataflow(hpx::launch::async,
                            process{xptr, yptr, beg, end},
                            x.safe_to_read[seg],
                            y.safe_to_write[seg]
                            );
                });
    }
};

template < typename real >
struct copy<
    std::vector<real>,
    hpx_vector<real>
    >
{
    typedef hpx_vector<real> vector;

    struct process {
        const real *xptr;
        real       *yptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        void operator()(T&&...) const {
            for(ptrdiff_t i = beg; i < end; ++i)
                yptr[i] = xptr[i];
        }
    };

    static void apply(const std::vector<real> &x, vector &y)
    {
        const real *xptr = &x[0];
        real       *yptr = &y[0];

        using hpx::dataflow;

        auto range = boost::irange(0, y.nseg);
        hpx::parallel::for_each(
                hpx::parallel::par,
                boost::begin(range), boost::end(range),
                [&y, xptr, yptr](ptrdiff_t seg) {
                    ptrdiff_t beg = seg * y.grain_size;
                    ptrdiff_t end = std::min<ptrdiff_t>(beg + y.grain_size, y.size());

                    y.safe_to_read[seg] = dataflow(hpx::launch::async,
                            process{xptr, yptr, beg, end},
                            y.safe_to_write[seg]
                            );
                });
    }
};

template < typename real >
struct inner_product_impl<
    hpx_vector<real>,
    hpx_vector<real>
    >
{
    typedef hpx_vector<real> vector;

    struct process {
        const real *xptr;
        const real *yptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        double operator()(T&&...) const {
            real sum = 0;

            for(ptrdiff_t i = beg; i < end; ++i)
                sum += xptr[i] * yptr[i];

            return sum;
        }
    };

    static real get(const vector &x, const vector &y)
    {
        const real *xptr = &x[0];
        const real *yptr = &y[0];

        using hpx::dataflow;

        auto range = boost::irange(0, x.nseg);
        return hpx::parallel::transform_reduce(
                hpx::parallel::par,
                boost::begin(range), boost::end(range),
                math::zero<real>(), std::plus<real>(),
                [&x, &y, xptr, yptr](ptrdiff_t seg) {
                    ptrdiff_t beg = seg * x.grain_size;
                    ptrdiff_t end = std::min<ptrdiff_t>(beg + x.grain_size, x.size());

                    return dataflow(hpx::launch::async,
                            process{xptr, yptr, beg, end},
                            x.safe_to_read[seg],
                            y.safe_to_read[seg]
                            ).get();
                }
                );
    }
};

template < typename A, typename B, typename real >
struct axpby_impl<
    A, hpx_vector<real>,
    B, hpx_vector<real>
    >
{
    typedef hpx_vector<real> vector;

    struct process_ab {
        typedef void result_type;

        A           a;
        const real *xptr;
        B           b;
        real       *yptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        void operator()(T&&...) const {
            for(ptrdiff_t i = beg; i < end; ++i)
                yptr[i] = a * xptr[i] + b * yptr[i];
        }
    };

    struct process_a {
        typedef void result_type;

        A           a;
        const real *xptr;
        real       *yptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        void operator()(T&&...) const {
            for(ptrdiff_t i = beg; i < end; ++i)
                yptr[i] = a * xptr[i];
        }
    };

    static void apply(A a, const vector &x, B b, vector &y)
    {
        const real *xptr = &x[0];
        real       *yptr = &y[0];

        using hpx::dataflow;

        auto range = boost::irange(0, x.nseg);
        if (b) {
            // y = a * x + b * y;
            hpx::parallel::for_each(
                    hpx::parallel::par,
                    boost::begin(range), boost::end(range),
                    [a, &x, b, &y, xptr, yptr](ptrdiff_t seg) {
                        ptrdiff_t beg = seg * x.grain_size;
                        ptrdiff_t end = std::min<ptrdiff_t>(beg + x.grain_size, x.size());

                        y.safe_to_read[seg] = dataflow(hpx::launch::async,
                                process_ab{a, xptr, b, yptr, beg, end},
                                x.safe_to_read[seg],
                                y.safe_to_read[seg],
                                y.safe_to_write[seg]
                                );
                    });
        } else {
            // y = a * x;
            hpx::parallel::for_each(
                    hpx::parallel::par,
                    boost::begin(range), boost::end(range),
                    [a, &x, &y, xptr, yptr](ptrdiff_t seg) {
                        ptrdiff_t beg = seg * x.grain_size;
                        ptrdiff_t end = std::min<ptrdiff_t>(beg + x.grain_size, x.size());

                        y.safe_to_read[seg] = dataflow(hpx::launch::async,
                                process_a{a, xptr, yptr, beg, end},
                                x.safe_to_read[seg],
                                y.safe_to_write[seg]
                                );
                    });
        }
    }
};

template < typename A, typename B, typename C, typename real >
struct axpbypcz_impl<
    A, hpx_vector<real>,
    B, hpx_vector<real>,
    C, hpx_vector<real>
    >
{
    typedef hpx_vector<real> vector;

    struct process_abc {
        A           a;
        const real *xptr;
        B           b;
        const real *yptr;
        C           c;
        real       *zptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        void operator()(T&&...) const {
            for(ptrdiff_t i = beg; i < end; ++i)
                zptr[i] = a * xptr[i] + b * yptr[i] + c * zptr[i];
        }
    };

    struct process_ab {
        A           a;
        const real *xptr;
        B           b;
        const real *yptr;
        real       *zptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        void operator()(T&&...) const {
            for(ptrdiff_t i = beg; i < end; ++i)
                zptr[i] = a * xptr[i] + b * yptr[i];
        }
    };

    static void apply(
            A a, const vector &x,
            B b, const vector &y,
            C c,       vector &z
            )
    {
        const real *xptr = &x[0];
        const real *yptr = &y[0];
        real       *zptr = &z[0];

        using hpx::dataflow;

        auto range = boost::irange(0, x.nseg);
        if (c) {
            //z = a * x + b * y + c * z;
            hpx::parallel::for_each(
                    hpx::parallel::par,
                    boost::begin(range), boost::end(range),
                    [a, &x, b, &y, c, &z, xptr, yptr, zptr](ptrdiff_t seg) {
                        ptrdiff_t beg = seg * x.grain_size;
                        ptrdiff_t end = std::min<ptrdiff_t>(beg + x.grain_size, x.size());

                        z.safe_to_read[seg] = dataflow(hpx::launch::async,
                                process_abc{a, xptr, b, yptr, c, zptr, beg, end},
                                x.safe_to_read[seg],
                                y.safe_to_read[seg],
                                z.safe_to_read[seg],
                                z.safe_to_write[seg]
                                );
                    });
        } else {
            //z = a * x + b * y;
            hpx::parallel::for_each(
                    hpx::parallel::par,
                    boost::begin(range), boost::end(range),
                    [a, &x, b, &y, &z, xptr, yptr, zptr](ptrdiff_t seg) {
                        ptrdiff_t beg = seg * x.grain_size;
                        ptrdiff_t end = std::min<ptrdiff_t>(beg + x.grain_size, x.size());

                        z.safe_to_read[seg] = dataflow(hpx::launch::async,
                                process_ab{a, xptr, b, yptr, zptr, beg, end},
                                x.safe_to_read[seg],
                                y.safe_to_read[seg],
                                z.safe_to_write[seg]
                                );
                    });
        }
    }
};

template < typename A, typename B, typename real >
struct vmul_impl<
    A, hpx_vector<real>, hpx_vector<real>,
    B, hpx_vector<real>
    >
{
    typedef hpx_vector<real> vector;

    struct process_ab {
        A           a;
        const real *xptr;
        const real *yptr;
        B           b;
        real       *zptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        void operator()(T&&...) const {
            for(ptrdiff_t i = beg; i < end; ++i)
                zptr[i] = a * xptr[i] * yptr[i] + b * zptr[i];
        }
    };

    struct process_a {
        A           a;
        const real *xptr;
        const real *yptr;
        real       *zptr;

        ptrdiff_t beg;
        ptrdiff_t end;

        template <class... T>
        void operator()(T&&...) const {
            for(ptrdiff_t i = beg; i < end; ++i)
                zptr[i] = a * xptr[i] * yptr[i];
        }
    };

    static void apply(A a, const vector &x, const vector &y, B b, vector &z)
    {
        const real *xptr = &x[0];
        const real *yptr = &y[0];
        real       *zptr = &z[0];

        using hpx::dataflow;

        auto range = boost::irange(0, x.nseg);
        if (b) {
            //z = a * x * y + b * z;
            hpx::parallel::for_each(
                    hpx::parallel::par,
                    boost::begin(range), boost::end(range),
                    [a, &x, &y, b, &z, xptr, yptr, zptr](ptrdiff_t seg) {
                        ptrdiff_t beg = seg * x.grain_size;
                        ptrdiff_t end = std::min<ptrdiff_t>(beg + x.grain_size, x.size());

                        z.safe_to_read[seg] = dataflow(hpx::launch::async,
                                process_ab{a, xptr, yptr, b, zptr, beg, end},
                                x.safe_to_read[seg],
                                y.safe_to_read[seg],
                                z.safe_to_read[seg],
                                z.safe_to_write[seg]
                                );
                    });
        } else {
            //z = a * x * y;
            hpx::parallel::for_each(
                    hpx::parallel::par,
                    boost::begin(range), boost::end(range),
                    [a, &x, &y, &z, xptr, yptr, zptr](ptrdiff_t seg) {
                        ptrdiff_t beg = seg * x.grain_size;
                        ptrdiff_t end = std::min<ptrdiff_t>(beg + x.grain_size, x.size());

                        z.safe_to_read[seg] = dataflow(hpx::launch::async,
                                process_a{a, xptr, yptr, zptr, beg, end},
                                x.safe_to_read[seg],
                                y.safe_to_read[seg],
                                z.safe_to_write[seg]
                                );
                    });
        }
    }
};

} // namespace backend
} // namespace amgcl

#endif
