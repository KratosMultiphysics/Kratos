// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef NNLS_INCLUDE_H
#define NNLS_INCLUDE_H

//// STL includes
#include <vector>

namespace queso {
namespace nnls {

typedef std::vector<std::vector<double>> MatrixType;
typedef std::vector<double> VectorType;

// Wrapper for nnls solver
double nnls(MatrixType& A, const VectorType& b, VectorType& x);

} // End Namespace nnls
} // End namespace queso

#endif // NNLS_INCLUDE_H