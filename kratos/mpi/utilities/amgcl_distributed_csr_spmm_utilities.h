
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
#include "mpi/utilities/amgcl_distributed_csr_conversion_utilities.h"

#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/zero_copy.hpp>

namespace Kratos
{

/**
Utilities to convert the distributed_csr matrix to other libraries
 */
class DistributedAmgclCSRSpMMUtilities
{

public:

    /**
     This function returns a shared pointer to an Amgcl distributed_matrix
     */
    template< class TDataType=double, class TIndexType=std::size_t >
    static DistributedCsrMatrix<TDataType, TIndexType> SparseMultiply(
        const DistributedCsrMatrix<TDataType, TIndexType>& rA,
        const DistributedCsrMatrix<TDataType, TIndexType>& rB
        )
	{
                
        bool move_to_backend=false;

        auto Aoffdiag_global_index2 = rA.GetOffDiagonalIndex2DataInGlobalNumbering();

    	auto pAamgcl = AmgclDistributedCSRConversionUtilities::ConvertToAmgcl<double,IndexType>
            (rA,Aoffdiag_global_index2,move_to_backend);

        auto Boffdiag_global_index2 = rB.GetOffDiagonalIndex2DataInGlobalNumbering();
    	auto pBamgcl = AmgclDistributedCSRConversionUtilities::ConvertToAmgcl<double,IndexType>
            (rB,Boffdiag_global_index2,move_to_backend);

        auto Camgcl = product(*pAamgcl, *pBamgcl);

        sort_rows(*Camgcl);

        return std::move(AmgclDistributedCSRConversionUtilities::ConvertToCsrMatrix<TDataType,IndexType>(*Camgcl));
	}



};

}

#endif // KRATOS_CSR_SPMM_UTILITIES_H_INCLUDED
