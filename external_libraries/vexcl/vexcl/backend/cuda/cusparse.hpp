#ifndef VEXCL_BACKEND_CUDA_CUSPARSE_HPP
#define VEXCL_BACKEND_CUDA_CUSPARSE_HPP

#include <map>
#include <memory>

#include <cusparse_v2.h>

#include <vexcl/vector.hpp>
#include <vexcl/cache.hpp>
#include <vexcl/backend/cuda/error.hpp>
#include <vexcl/backend/cuda/context.hpp>
#include <vexcl/detail/backtrace.hpp>

namespace vex {
namespace backend {
namespace cuda {

/// Send human-readable representation of CUresult to the output stream.
inline std::ostream& operator<<(std::ostream &os, cusparseStatus_t rc) {
    os << "CUSPARSE Error (";
#define VEXCL_CUDA_ERR2TXT(e) case e: os << static_cast<int>(e) << " - " << #e; break
    switch(rc) {
        VEXCL_CUDA_ERR2TXT(CUSPARSE_STATUS_SUCCESS);
        VEXCL_CUDA_ERR2TXT(CUSPARSE_STATUS_NOT_INITIALIZED);
        VEXCL_CUDA_ERR2TXT(CUSPARSE_STATUS_ALLOC_FAILED);
        VEXCL_CUDA_ERR2TXT(CUSPARSE_STATUS_INVALID_VALUE);
        VEXCL_CUDA_ERR2TXT(CUSPARSE_STATUS_ARCH_MISMATCH);
        VEXCL_CUDA_ERR2TXT(CUSPARSE_STATUS_MAPPING_ERROR);
        VEXCL_CUDA_ERR2TXT(CUSPARSE_STATUS_EXECUTION_FAILED);
        VEXCL_CUDA_ERR2TXT(CUSPARSE_STATUS_INTERNAL_ERROR);
        VEXCL_CUDA_ERR2TXT(CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED);
        default:
            os << "Unknown error";
    }
#undef VEXCL_CUDA_ERR2TXT
    return os << ")";
}

inline void check(cusparseStatus_t rc, const char *file, int line) {
    if (rc != CUSPARSE_STATUS_SUCCESS) {
        vex::detail::print_backtrace();
        throw error(rc, file, line);
    }
}

namespace detail {

template <>
struct deleter_impl<cusparseHandle_t> {
    static void dispose(cusparseHandle_t handle) {
        cuda_check( cusparseDestroy(handle) );
    }
};

template <>
struct deleter_impl<cusparseMatDescr_t> {
    static void dispose(cusparseMatDescr_t handle) {
        cuda_check( cusparseDestroyMatDescr(handle) );
    }
};

template <>
struct deleter_impl<cusparseHybMat_t> {
    static void dispose(cusparseHybMat_t handle) {
        cuda_check( cusparseDestroyHybMat(handle) );
    }
};

} // namespace detail

inline cusparseHandle_t cusparse_handle(const command_queue &q) {
    typedef std::shared_ptr<std::remove_pointer<cusparseHandle_t>::type> smart_handle;
    typedef vex::detail::object_cache<vex::detail::index_by_context, smart_handle> cache_type;

    static cache_type cache;

    auto h = cache.find(q);

    if (h == cache.end()) {
        select_context(q);
        cusparseHandle_t handle;
        cuda_check( cusparseCreate(&handle) );
        cuda_check( cusparseSetStream(handle, q.raw()) );

        h = cache.insert(q, smart_handle(handle, detail::deleter(q.context().raw())));
    }

    return h->second.get();
}

template <typename val_t>
class spmat_hyb {
    static_assert(
            std::is_same<val_t, float>::value ||
            std::is_same<val_t, double>::value,
            "Unsupported value type for spmat_cusparse"
            );

    public:
        template <typename row_t, typename col_t>
        spmat_hyb(
                const command_queue &queue,
                int n, int m,
                const row_t *row_begin,
                const col_t *col_begin,
                const val_t *val_begin
                )
            : handle( cusparse_handle(queue) ),
              desc  ( create_description(), detail::deleter(queue.context().raw()) ),
              mat   ( create_matrix(),      detail::deleter(queue.context().raw()) )
        {
            cuda_check( cusparseSetMatType(desc.get(), CUSPARSE_MATRIX_TYPE_GENERAL) );
            cuda_check( cusparseSetMatIndexBase(desc.get(), CUSPARSE_INDEX_BASE_ZERO) );

            fill_matrix(queue, n, m, row_begin, col_begin, val_begin);
        }

        void apply(const vex::vector<val_t> &x, vex::vector<val_t> &y,
                 val_t alpha = 1, bool append = false) const
        {
            precondition(x.nparts() == 1 && y.nparts() == 1,
                    "Incompatible vectors");

            mul(x(0), y(0), alpha, append);
        }

        void mul(const device_vector<float> &x, device_vector<float> &y,
                 float alpha = 1, bool append = false) const
        {
            float beta = append ? 1.0f : 0.0f;

            cuda_check(
                    cusparseShybmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha, desc.get(), mat.get(),
                        x.raw_ptr(), &beta, y.raw_ptr()
                        )
                    );
        }

        void mul(const device_vector<double> &x, device_vector<double> &y,
                 double alpha = 1, bool append = false) const
        {
            double beta = append ? 1.0 : 0.0;

            cuda_check(
                    cusparseDhybmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha, desc.get(), mat.get(),
                        x.raw_ptr(), &beta, y.raw_ptr()
                        )
                    );
        }
    private:
        cusparseHandle_t handle;

        std::shared_ptr<std::remove_pointer<cusparseMatDescr_t>::type> desc;
        std::shared_ptr<std::remove_pointer<cusparseHybMat_t>::type>   mat;

        static cusparseMatDescr_t create_description() {
            cusparseMatDescr_t desc;
            cuda_check( cusparseCreateMatDescr(&desc) );
            return desc;
        }

        static cusparseHybMat_t create_matrix() {
            cusparseHybMat_t mat;
            cuda_check( cusparseCreateHybMat(&mat) );
            return mat;
        }

        template <typename row_t, typename col_t>
        void fill_matrix(const command_queue &q,
                int n, int m, const row_t *row, const col_t *col, const float *val)
        {
            device_vector<int>   r(q, n + 1,  row);
            device_vector<int>   c(q, row[n], col + row[0]);
            device_vector<float> v(q, row[n], val + row[0]);

            if (row[0] != 0) vector<int>(q, r) -= row[0];

            cuda_check(
                    cusparseScsr2hyb(handle, n, m, desc.get(),
                        v.raw_ptr(), r.raw_ptr(), c.raw_ptr(), mat.get(), 0,
                        CUSPARSE_HYB_PARTITION_AUTO
                        )
                    );
        }

        template <typename row_t, typename col_t>
        void fill_matrix(const command_queue &q,
                int n, int m, const row_t *row, const col_t *col, const double *val)
        {
            device_vector<int>    r(q, n + 1,  row);
            device_vector<int>    c(q, row[n], col + row[0]);
            device_vector<double> v(q, row[n], val + row[0]);

            if (row[0] != 0) vector<int>(q, r) -= row[0];

            cuda_check(
                    cusparseDcsr2hyb(handle, n, m, desc.get(),
                        v.raw_ptr(), r.raw_ptr(), c.raw_ptr(), mat.get(), 0,
                        CUSPARSE_HYB_PARTITION_AUTO
                        )
                    );
        }

};

template <typename T>
additive_operator< spmat_hyb<T>, vector<T> >
operator*(const spmat_hyb<T> &A, const vector<T> &x) {
    return additive_operator< spmat_hyb<T>, vector<T> >(A, x);
}

template <typename val_t>
class spmat_crs {
    static_assert(
            std::is_same<val_t, float>::value ||
            std::is_same<val_t, double>::value,
            "Unsupported value type for spmat_cusparse"
            );

    public:
        template <typename row_t, typename col_t>
        spmat_crs(
                const command_queue &queue,
                int n, int m,
                const row_t *row_begin,
                const col_t *col_begin,
                const val_t *val_begin
                )
            : n(n), m(m), nnz(static_cast<unsigned>(row_begin[n] - row_begin[0])),
              handle( cusparse_handle(queue) ),
              desc  ( create_description(), detail::deleter(queue.context().raw()) ),
              row(queue, n+1, row_begin),
              col(queue, nnz, col_begin + row_begin[0]),
              val(queue, nnz, val_begin + row_begin[0])
        {
            if (row_begin[0] != 0)
                vector<int>(queue, row) -= row_begin[0];

            cuda_check( cusparseSetMatType(desc.get(), CUSPARSE_MATRIX_TYPE_GENERAL) );
            cuda_check( cusparseSetMatIndexBase(desc.get(), CUSPARSE_INDEX_BASE_ZERO) );
        }

        void apply(const vex::vector<val_t> &x, vex::vector<val_t> &y,
                 val_t alpha = 1, bool append = false) const
        {
            precondition(x.nparts() == 1 && y.nparts() == 1,
                    "Incompatible vectors");

            mul(x(0), y(0), alpha, append);
        }

        void mul(const device_vector<float> &x, device_vector<float> &y,
                 float alpha = 1, bool append = false) const
        {
            float beta = append ? 1.0f : 0.0f;

            cuda_check(
                    cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        n, m, nnz, &alpha, desc.get(),
                        val.raw_ptr(), row.raw_ptr(), col.raw_ptr(),
                        x.raw_ptr(), &beta, y.raw_ptr()
                        )
                 );
        }

        void mul(const device_vector<double> &x, device_vector<double> &y,
                 double alpha = 1, bool append = false) const
        {
            double beta = append ? 1.0 : 0.0;

            cuda_check(
                    cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                        n, m, nnz, &alpha, desc.get(),
                        val.raw_ptr(), row.raw_ptr(), col.raw_ptr(),
                        x.raw_ptr(), &beta, y.raw_ptr()
                        )
                 );
        }
    private:
        unsigned n, m, nnz;

        cusparseHandle_t handle;

        std::shared_ptr<std::remove_pointer<cusparseMatDescr_t>::type> desc;

        device_vector<int>   row;
        device_vector<int>   col;
        device_vector<val_t> val;

        static cusparseMatDescr_t create_description() {
            cusparseMatDescr_t desc;
            cuda_check( cusparseCreateMatDescr(&desc) );
            return desc;
        }
};

template <typename T>
additive_operator< spmat_crs<T>, vector<T> >
operator*(const spmat_crs<T> &A, const vector<T> &x) {
    return additive_operator< spmat_crs<T>, vector<T> >(A, x);
}

} // namespace cuda
} // namespace backend
} // namespace vex

#endif
