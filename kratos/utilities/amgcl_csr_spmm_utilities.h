
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_CSR_SPMM_UTILITIES_H_INCLUDED)
#define  KRATOS_CSR_SPMM_UTILITIES_H_INCLUDED

// System includes

// External includes
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/zero_copy.hpp>

// Project includes
#include "containers/csr_matrix.h"
#include "utilities/amgcl_csr_conversion_utilities.h"

namespace Kratos
{

/**
Utilities to convert the distributed_csr matrix to other libraries
 */
class AmgclCSRSpMMUtilities
{

public:

    /**
     This function returns a shared pointer to an Amgcl distributed_matrix
     */
    template< class TDataType=double, class TIndexType=std::size_t >
    static typename CsrMatrix<TDataType, TIndexType>::Pointer SparseMultiply(
        const CsrMatrix<TDataType, TIndexType>& rA,
        const CsrMatrix<TDataType, TIndexType>& rB
    )
    {
        auto pAamgcl = AmgclCSRConversionUtilities::ConvertToAmgcl<TDataType,TIndexType>(rA);
        auto pBamgcl = AmgclCSRConversionUtilities::ConvertToAmgcl<TDataType,TIndexType>(rB);

        auto Camgcl = amgcl::backend::product(*pAamgcl, *pBamgcl);
        amgcl::backend::sort_rows(*Camgcl);

        return AmgclCSRConversionUtilities::ConvertToCsrMatrix<TDataType,TIndexType>(*Camgcl);
    }



};

}

#endif // KRATOS_CSR_SPMM_UTILITIES_H_INCLUDED
