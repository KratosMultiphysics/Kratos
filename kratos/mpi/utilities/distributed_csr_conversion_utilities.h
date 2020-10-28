
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

#if !defined(KRATOS_DISTRIBUTED_CSR_CONVERSION_UTILITIES_H_INCLUDED)
#define  KRATOS_DISTRIBUTED_CSR_CONVERSION_UTILITIES_H_INCLUDED

#include "containers/distributed_csr_matrix.h"
#include "../../../external_libraries/amgcl/backend/builtin.hpp"
#include "../../../external_libraries/amgcl/backend/interface.hpp"
#include "../../../external_libraries/amgcl/mpi/distributed_matrix.hpp"
#include "../../../external_libraries/amgcl/adapter/zero_copy.hpp"

namespace Kratos
{

/**
Utilities to convert the distributed_csr matrix to other libraries
 */
class DistributedCSRConversionUtilities
{

public:

    /**
     This function returns a shared pointer to an Amgcl distributed_matrix
     */
	template< class TDataType, class TIndexType >
	static Kratos::shared_ptr<amgcl::mpi::distributed_matrix<amgcl::backend::builtin<double>>> ConvertToAmgcl(
			const DistributedCsrMatrix<TDataType, TIndexType>& rA,
			DenseVector<TIndexType>& global_index2)
	{
		IndexType chunk = rA.local_size1();
    	auto loc_a = amgcl::adapter::zero_copy(chunk, 
			&rA.GetDiagBlock().index1_data()[0],
			&rA.GetDiagBlock().index2_data()[0],
			&rA.GetDiagBlock().value_data()[0]
			);

		auto rem_a = amgcl::adapter::zero_copy(chunk, 
			&rA.GetOffDiagBlock().index1_data()[0],
			&global_index2[0],
			&rA.GetOffDiagBlock().value_data()[0]
		);

		amgcl::mpi::communicator comm(MPI_COMM_WORLD);
		auto pAmgcl = Kratos::make_shared<amgcl::mpi::distributed_matrix<amgcl::backend::builtin<double>>>(comm, loc_a, rem_a);
		pAmgcl->move_to_backend();
		return pAmgcl;
	}

};

}

#endif // KRATOS_DISTRIBUTED_CSR_CONVERSION_UTILITIES_H_INCLUDED
