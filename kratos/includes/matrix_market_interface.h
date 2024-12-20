//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes
#include <stdio.h>


// External includes

// To avoid linking problems
extern "C"
{
#include "includes/mmio.h"
}

#include <boost/numeric/ublas/matrix_sparse.hpp>


// Project includes

namespace Kratos
{

// Type checks
template<typename T>
constexpr bool KRATOS_API(KRATOS_CORE) IsCorrectType(MM_typecode& mm_code);

// Matrix I/O routines
bool KRATOS_API(KRATOS_CORE) ReadMatrixMarketMatrixEntry(FILE *f, int& I, int& J, double& V);
bool KRATOS_API(KRATOS_CORE) ReadMatrixMarketMatrixEntry(FILE *f, int& I, int& J, std::complex<double>& V);

template <typename CompressedMatrixType> 
bool KRATOS_API(KRATOS_CORE) ReadMatrixMarketMatrix(const char *FileName, CompressedMatrixType &M);


void KRATOS_API(KRATOS_CORE) SetMatrixMarketValueTypeCode(MM_typecode& mm_code, const double& value);
void KRATOS_API(KRATOS_CORE) SetMatrixMarketValueTypeCode(MM_typecode& mm_code, const std::complex<double>& value);

int KRATOS_API(KRATOS_CORE) WriteMatrixMarketMatrixEntry(FILE *f, int I, int J, const double& entry);
int KRATOS_API(KRATOS_CORE) WriteMatrixMarketMatrixEntry(FILE *f, int I, int J, const std::complex<double>& entry);

template <typename CompressedMatrixType> 
bool KRATOS_API(KRATOS_CORE) WriteMatrixMarketMatrix(const char *FileName, CompressedMatrixType &M, bool Symmetric);

// Vector I/O routines
bool KRATOS_API(KRATOS_CORE) ReadMatrixMarketVectorEntry(FILE *f, double& entry);
bool KRATOS_API(KRATOS_CORE) ReadMatrixMarketVectorEntry(FILE *f, std::complex<double>& entry);

template <typename VectorType> 
bool KRATOS_API(KRATOS_CORE) ReadMatrixMarketVector(const char *FileName, VectorType &V);

int KRATOS_API(KRATOS_CORE) WriteMatrixMarketVectorEntry(FILE *f, const double& entry);
int KRATOS_API(KRATOS_CORE) WriteMatrixMarketVectorEntry(FILE *f, const std::complex<double>& entry);

template <typename VectorType> 
bool KRATOS_API(KRATOS_CORE) WriteMatrixMarketVector(const char *FileName, VectorType &V);

} // namespace Kratos