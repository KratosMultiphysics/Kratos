#include <iostream>
#include <vector>

#include <vexcl/vexcl.hpp>
#include <vexcl/external/viennacl.hpp>

#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/bicgstab.hpp>

typedef double real;

template <class Tag>
void do_solve(const vex::SpMat<real> &A, vex::vector<real> f, const Tag &tag)
{
    vex::vector<real> x = viennacl::linalg::solve(A, f, tag);

    std::cout << "  Iterations: " << tag.iters() << std::endl
              << "  Error:      " << tag.error() << std::endl;

    // Test for convergence.
    f -= A * x;

    static vex::Reductor<real, vex::MAX> max(vex::current_context().queue());
    std::cout << "  max(residual) = " << max(fabs(f)) << std::endl;
}

int main() {
    try {
        size_t n = 1024;
        real h2i = (n - 1) * (n - 1);

        std::vector<size_t> row;
        std::vector<size_t> col;
        std::vector<real>   val;
        std::vector<real>   rhs;

        // Prepare problem (1D Poisson equation).
        row.reserve(n + 1);
        col.reserve(2 + (n - 2) * 3);
        val.reserve(2 + (n - 2) * 3);
        rhs.reserve(n);

        row.push_back(0);
        for(size_t i = 0; i < n; i++) {
            if (i == 0 || i == n - 1) {
                col.push_back(i);
                val.push_back(1);
                rhs.push_back(0);
                row.push_back(row.back() + 1);
            } else {
                col.push_back(i-1);
                val.push_back(-h2i);

                col.push_back(i);
                val.push_back(2 * h2i);

                col.push_back(i+1);
                val.push_back(-h2i);

                rhs.push_back(2);
                row.push_back(row.back() + 3);
            }
        }

        // Move data to GPU(s).
        vex::Context ctx(vex::Filter::Env && vex::Filter::DoublePrecision);
        if (!ctx) throw std::runtime_error("No GPUs with double precision found");
        std::cout << ctx << std::endl;

        vex::SpMat <real> A(ctx, n, n, row.data(), col.data(), val.data());
        vex::vector<real> f(ctx, rhs);

        vex::profiler<> prof(ctx);

        // Solve problem with ViennaCL's solvers:
        std::cout << "CG" << std::endl;
        prof.tic_cl("CG");
        do_solve(A, f, viennacl::linalg::cg_tag(1e-8, n));
        prof.toc("CG");

        std::cout << "BiCGStab" << std::endl;
        prof.tic_cl("BiCGStab");
        do_solve(A, f, viennacl::linalg::bicgstab_tag(1e-8, n));
        prof.toc("BiCGStab");

        std::cout << prof << std::endl;

    } catch(const std::exception &err) {
        std::cerr << "Error: " << err.what() << std::endl;
    }
}

// vim: et
