#ifndef AMGCL_ADAPTER_ZERO_COPY_HPP
#define AMGCL_ADAPTER_ZERO_COPY_HPP

/*
The MIT License

Copyright (c) 2012-2016 Denis Demidov <dennis.demidov@gmail.com>

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
\file    amgcl/adapter/zero_copy.hpp
\author  Denis Demidov <dennis.demidov@gmail.com>
\brief   Zero-copy adapter for input matrix in CRS format.
\ingroup adapters
*/

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include <amgcl/util.hpp>
#include <amgcl/backend/builtin.hpp>

namespace amgcl {
namespace adapter {

template <typename value_type>
struct zero_copy_adapter
    : public backend::crs<value_type, ptrdiff_t, ptrdiff_t>
{
    typedef backend::crs<value_type, ptrdiff_t, ptrdiff_t> Base;

    typedef ptrdiff_t  ptr_type;
    typedef ptrdiff_t  col_type;
    typedef value_type val_type;

    const ptr_type *ptr;
    const col_type *col;
    const val_type *val;

    template <typename Ptr, typename Col>
    zero_copy_adapter(size_t n, const Ptr *ptr, const Col *col, const value_type *val)
        : ptr(reinterpret_cast<const ptr_type*>(ptr)),
          col(reinterpret_cast<const col_type*>(col)),
          val(val)
    {
        // Check that Ptr and Col types are binary-compatible with ptrdiff_t:
        BOOST_STATIC_ASSERT(boost::is_integral<Ptr>::value);
        BOOST_STATIC_ASSERT(boost::is_integral<Col>::value);
        BOOST_STATIC_ASSERT(sizeof(Ptr) == sizeof(ptr_type));
        BOOST_STATIC_ASSERT(sizeof(Col) == sizeof(col_type));

        this->nrows = this->ncols = n;
    }

    virtual const ptr_type* ptr_data() const { return ptr; }
    virtual const col_type* col_data() const { return col; }
    virtual const val_type* val_data() const { return val; }

    virtual ptr_type* ptr_data() { precondition(false, "The matrix is read-only"); return 0; }
    virtual col_type* col_data() { precondition(false, "The matrix is read-only"); return 0; }
    virtual val_type* val_data() { precondition(false, "The matrix is read-only"); return 0; }
};

template <typename Ptr, typename Col, typename Val>
boost::shared_ptr< backend::crs<Val, ptrdiff_t, ptrdiff_t> >
zero_copy(size_t n, const Ptr *ptr, const Col *col, const Val *val) {
    return boost::make_shared< zero_copy_adapter<Val> >(n, ptr, col, val);
}

} // namespace adapter
} // namespace amgcl

#endif
