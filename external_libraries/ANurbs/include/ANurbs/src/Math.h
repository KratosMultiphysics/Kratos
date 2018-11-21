#pragma once

#include <vector>

namespace ANurbs {
namespace Math {

static constexpr inline int
Binom(
    const int n,
    const int k) noexcept
{
    // clang-format off
    return (k > n               ) ? 0           :  // out of range
           (k == 0 || k == n    ) ? 1           :  // edge
           (k == 1 || k == n - 1) ? n           :  // first
           (k + k < n           ) ?                // recursive:
           (Binom(n - 1, k - 1  ) * n) / k      :  //   path to k = 1     faster
           (Binom(n - 1, k      ) * n) / (n - k);  //   path to k = n - 1 faster
    // clang-format on
}

static constexpr inline int
MatrixIndex(
    const int rows,
    const int cols,
    const int row,
    const int col) noexcept
{
    return row * cols + col;
}

template <typename TScalar, typename TFunction>
TScalar
Romberg(
    TFunction f,
    TScalar a,
    TScalar b,
    int maxSteps,
    TScalar tolerance)
{
    std::vector<TScalar> r1(maxSteps);
    std::vector<TScalar> r2(maxSteps);

    TScalar* Rp = &r1[0]; // current row
    TScalar* Rc = &r2[0]; // previous row

    TScalar h = (b - a); // step size

    Rp[0] = (f(a) + f(b)) * h * 0.5; // first trapezoidal step

    for (int i = 1; i < maxSteps; ++i) {
        h /= 2.0;

        TScalar c = 0;

        int ep = 1 << (i - 1); // 2^(n-1)

        for (int j = 1; j <= ep; ++j) {
            c += f(a + (2 * j - 1) * h);
        }

        Rc[0] = h * c + 0.5 * Rp[0]; // R(i,0)

        for (int j = 1; j <= i; ++j) {
            TScalar n_k = pow(4, j);
            Rc[j] = (n_k * Rc[j - 1] - Rp[j - 1]) / (n_k - 1); // compute R(i,j)
        }

        if (i > 1 && fabs(Rp[i - 1] - Rc[i]) < tolerance) {
            return Rc[i - 1];
        }

        // swap Rn and Rc as we only need the last row
        TScalar* rt = Rp;

        Rp = Rc;
        Rc = rt;
    }

    return Rp[maxSteps - 1];
}

} // namespace Math
} // namespace ANurbs
