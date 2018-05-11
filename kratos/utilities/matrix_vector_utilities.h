//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//                    
//


#if !defined(KRATOS_MATRIX_VECTOR_UTILITIES_INCLUDED )
#define  KRATOS_MATRIX_VECTOR_UTILITIES_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
class MatrixVectorUtils
{
public:
    static inline void InitializeMatrix(Matrix& rMat, std::size_t Dim1, std::size_t Dim2, bool Zero=false)
    {
        if (rMat.size1() != Dim1 || rMat.size2() != Dim2)
            rMat.resize(Dim1, Dim2, false);
        if (Zero)
            noalias(rMat) = ZeroMatrix(Dim1, Dim2);
    }

    static inline void InitializeVector(Vector& rVec, std::size_t Dim, bool Zero=false)
    {
        if (rVec.size() != Dim)
            rVec.resize(Dim, false);
        if (Zero)
            noalias(rVec) = ZeroVector(Dim);
    }
};

}  // namespace Kratos.

#endif // KRATOS_MATRIX_VECTOR_UTILITIES_INCLUDED  defined 


