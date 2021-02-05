
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

#include "containers/csr_matrix.h"
#include "utilities/amgcl_csr_conversion_utilities.h"

#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/zero_copy.hpp>

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
    static void SparseMultiply(
        const CsrMatrix<TDataType, TIndexType>& rA,
        const CsrMatrix<TDataType, TIndexType>& rB,
        CsrMatrix<TDataType, TIndexType>& rOutput
        )
	{
        //bool move_to_backend=false;
        auto pAamgcl = AmgclCSRConversionUtilities::ConvertToAmgcl<TDataType,IndexType>(rA);
        auto pBamgcl = AmgclCSRConversionUtilities::ConvertToAmgcl<TDataType,IndexType>(rB); 

        auto Camgcl = amgcl::backend::product(*pAamgcl, *pBamgcl);
        amgcl::backend::sort_rows(*Camgcl);

        AmgclCSRConversionUtilities::ConvertToCsrMatrix<TDataType,IndexType>(*Camgcl,rOutput);
	}



};

}

#endif // KRATOS_CSR_SPMM_UTILITIES_H_INCLUDED
