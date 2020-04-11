#ifndef AMGCL_COARSENING_RIGID_BODY_MODES_HPP
#define AMGCL_COARSENING_RIGID_BODY_MODES_HPP

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
 * \file   amgcl/coarsening/rigid_body_modes.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Create rigid body modes from coordinates.
 */

#include <amgcl/util.hpp>

namespace amgcl {
namespace coarsening {

// Create rigid body modes from coordinate vector.
// To be used as near-nullspace vectors with aggregation coarsening
// for 2D or 3D elasticity problems.
template <class Vector>
int rigid_body_modes(int ndim, const Vector &coo, std::vector<double> &B) {
    precondition(ndim == 2 || ndim == 3, "Only 2D or 3D problems are supported");
    precondition(coo.size() % ndim == 0, "Coordinate vector size should be divisible by ndim");

    size_t n = coo.size();
    int nmodes = (ndim == 2 ? 3 : 6);
    B.resize(n * nmodes, 0.0);

    double sn = 1 / sqrt(n);

    if (ndim == 2) {
        for(size_t i = 0; i < n; ++i) {
            size_t nod = i / ndim;
            size_t dim = i % ndim;

            double x = coo[nod * 2 + 0];
            double y = coo[nod * 2 + 1];

            // Translation
            B[i * nmodes + dim] = sn;

            // Rotation
            switch(dim) {
                case 0:
                    B[i * nmodes + 2] = -y;
                    break;
                case 1:
                    B[i * nmodes + 2] = x;
                    break;
            }
        }
    } else if (ndim == 3) {
        for(size_t i = 0; i < n; ++i) {
            size_t nod = i / ndim;
            size_t dim = i % ndim;

            double x = coo[nod * 3 + 0];
            double y = coo[nod * 3 + 1];
            double z = coo[nod * 3 + 2];

            // Translation
            B[i * nmodes + dim] = sn;

            // Rotation
            switch(dim) {
                case 0:
                    B[i * nmodes + 3] = y;
                    B[i * nmodes + 5] = z;
                    break;
                case 1:
                    B[i * nmodes + 3] = -x;
                    B[i * nmodes + 4] = -z;
                    break;
                case 2:
                    B[i * nmodes + 4] =  y;
                    B[i * nmodes + 5] = -x;
                    break;
            }
        }
    }

    // Orthonormalization
    std::array<double, 6> dot;
    for(int i = ndim; i < nmodes; ++i) {
        std::fill(dot.begin(), dot.end(), 0.0);
        for(size_t j = 0; j < n; ++j) {
            for(int k = 0; k < i; ++k)
                dot[k] += B[j * nmodes + k] * B[j * nmodes + i];
        }
        double s = 0.0;
        for(size_t j = 0; j < n; ++j) {
            for(int k = 0; k < i; ++k)
                B[j * nmodes + i] -= dot[k] * B[j * nmodes + k];
            s += B[j * nmodes + i] * B[j * nmodes + i];
        }
        s = sqrt(s);
        for(size_t j = 0; j < n; ++j)
            B[j * nmodes + i] /= s;
    }

    return nmodes;
}

} // namespace coarsening
} // namespace amgcl

#endif
