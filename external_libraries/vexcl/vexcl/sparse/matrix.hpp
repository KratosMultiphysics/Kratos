#ifndef VEXCL_SPARSE_MATRIX_HPP
#define VEXCL_SPARSE_MATRIX_HPP

#include <vexcl/sparse/ell.hpp>
#include <vexcl/sparse/csr.hpp>

namespace vex {
namespace sparse {

template <typename Val, typename Col = int, typename Ptr = Col>
class matrix {
    public:
        typedef Val value_type;

        typedef Val val_type;
        typedef Col col_type;
        typedef Ptr ptr_type;

        template <class PtrRange, class ColRange, class ValRange>
        matrix(
                const std::vector<backend::command_queue> &q,
                size_t nrows, size_t ncols,
                const PtrRange &ptr,
                const ColRange &col,
                const ValRange &val,
                bool fast_setup = true
           ) : q(q[0])
        {
            if (is_cpu(q[0])) {
                Acpu = std::make_shared<Csr>(q, nrows, ncols, ptr, col, val);
            } else {
                Agpu = std::make_shared<Ell>(q, nrows, ncols, ptr, col, val, fast_setup);
            }
        }

        // Dummy matrix
        matrix() {}
        
        // Dummy matrix; used internally to pass empty parameters to kernels.
        matrix(const backend::command_queue &q) : q(q) {}

        template <class Expr>
        friend
        typename std::enable_if<
            boost::proto::matches<
                typename boost::proto::result_of::as_expr<Expr>::type,
                vector_expr_grammar
            >::value,
            matrix_vector_product<matrix, Expr>
        >::type
        operator*(const matrix &A, const Expr &x) {
            return matrix_vector_product<matrix, Expr>(A, x);
        }

        template <class Vector>
        static void terminal_preamble(const Vector &x, backend::source_generator &src,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
        {
            if (is_cpu(q)) {
                Csr::terminal_preamble(x, src, q, prm_name, state);
            } else {
                Ell::terminal_preamble(x, src, q, prm_name, state);
            }
        }

        template <class Vector>
        static void local_terminal_init(const Vector &x, backend::source_generator &src,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
        {
            if (is_cpu(q)) {
                Csr::local_terminal_init(x, src, q, prm_name, state);
            } else {
                Ell::local_terminal_init(x, src, q, prm_name, state);
            }
        }

        template <class Vector>
        static void kernel_param_declaration(const Vector &x, backend::source_generator &src,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
        {
            if (is_cpu(q)) {
                Csr::kernel_param_declaration(x, src, q, prm_name, state);
            } else {
                Ell::kernel_param_declaration(x, src, q, prm_name, state);
            }
        }

        template <class Vector>
        static void partial_vector_expr(const Vector &x, backend::source_generator &src,
            const backend::command_queue &q, const std::string &prm_name,
            detail::kernel_generator_state_ptr state)
        {
            if (is_cpu(q)) {
                Csr::partial_vector_expr(x, src, q, prm_name, state);
            } else {
                Ell::partial_vector_expr(x, src, q, prm_name, state);
            }
        }

        template <class Vector>
        void kernel_arg_setter(const Vector &x,
            backend::kernel &kernel, unsigned part, size_t index_offset,
            detail::kernel_generator_state_ptr state) const
        {
            if (is_cpu(q)) {
                if (Acpu) {
                    Acpu->kernel_arg_setter(x, kernel, part, index_offset, state);
                } else {
                    Csr dummy_A(q);
                    dummy_A.kernel_arg_setter(x, kernel, part, index_offset, state);
                }
            } else {
                if (Agpu) {
                    Agpu->kernel_arg_setter(x, kernel, part, index_offset, state);
                } else {
                    Ell dummy_A(q);
                    dummy_A.kernel_arg_setter(x, kernel, part, index_offset, state);
                }
            }
        }

        template <class Vector>
        void expression_properties(const Vector &x,
            std::vector<backend::command_queue> &queue_list,
            std::vector<size_t> &partition,
            size_t &size) const
        {
            if (Acpu) {
                Acpu->expression_properties(x, queue_list, partition, size);
            } else if (Agpu) {
                Agpu->expression_properties(x, queue_list, partition, size);
            }
        }

        size_t rows()     const { return Acpu ? Acpu->rows()     : Agpu->rows();     }
        size_t cols()     const { return Acpu ? Acpu->cols()     : Agpu->cols();     }
        size_t nonzeros() const { return Acpu ? Acpu->nonzeros() : Agpu->nonzeros(); }
    private:
        typedef ell<Val, Col, Ptr> Ell;
        typedef csr<Val, Col, Ptr> Csr;

        backend::command_queue q;

        std::shared_ptr<Ell> Agpu;
        std::shared_ptr<Csr> Acpu;

};

} // namespace sparse
} // namespace vex

#endif
