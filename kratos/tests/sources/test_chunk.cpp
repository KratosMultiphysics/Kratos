//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   
//
	           
// System includes
#include <set>

// External includes 


// Project includes
#include "testing/testing.h"
#include "includes/chunk.h"




namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(ChunkInitialize, KratosCoreFastSuite)
		{
#ifdef _OPENMP
				int max_threads = omp_get_max_threads();
#else
				int max_threads = 1;
#endif
			std::size_t block_size_in_bytes = 5; // the aligned block size is 8
			std::size_t chunk_size_in_bytes = 5+2*max_threads* sizeof(Chunk::SizeType);
			Chunk too_small_chunk(block_size_in_bytes, chunk_size_in_bytes);
			KRATOS_CHECK_EQUAL(too_small_chunk.GetNumberOfAvailableBlocks(), 0) << " Available block :" << too_small_chunk.GetNumberOfAvailableBlocks();

			chunk_size_in_bytes = 1024;
			Chunk chunk(block_size_in_bytes, chunk_size_in_bytes);
			std::size_t block_size_after_alignment = 8;
			std::size_t available_blocks_should_be = (chunk_size_in_bytes - 2* max_threads * sizeof(Chunk::SizeType)) / block_size_after_alignment;
			KRATOS_CHECK_EQUAL(chunk.GetNumberOfAvailableBlocks(), available_blocks_should_be) << " Available block :" << chunk.GetNumberOfAvailableBlocks() << " vs " << available_blocks_should_be;

		}

		KRATOS_TEST_CASE_IN_SUITE(ChunkAllocateDeallocate, KratosCoreFastSuite)
		{
			std::size_t block_size_in_bytes = 5; 
			std::size_t chunk_size_in_bytes = 1024;
			Chunk chunk(block_size_in_bytes, chunk_size_in_bytes);

			chunk.Release();
			chunk.Initialize();

			std::set<void *> pointer_set;
			while (chunk.HasAvailableBlock())
			{
				void* p = chunk.Allocate();
				KRATOS_ERROR_IF(pointer_set.find(p) != pointer_set.end()) << "The allocated position #" << pointer_set.size() + 1 << " was allocated before";
				pointer_set.insert(p);
			}
			KRATOS_CHECK_EQUAL(pointer_set.size(), chunk.GetNumberOfBlocks()) << " (pointer set size : " << pointer_set.size() << ")" << std::endl;
			
			for (auto i_pointer = pointer_set.begin(); i_pointer != pointer_set.end(); i_pointer++) {
				chunk.Deallocate(*i_pointer);
			}
			KRATOS_CHECK_EQUAL(chunk.GetNumberOfAvailableBlocks(), chunk.GetNumberOfBlocks());
		}

		KRATOS_TEST_CASE_IN_SUITE(ChunkParallelAllocate, KratosCoreFastSuite)
		{
			std::size_t block_size = 5;
			std::size_t chunk_size_in_bytes = 1024;

			auto repeat_number = 10;
			std::stringstream buffer;
#pragma omp parallel for 
			for (auto i_repeat = 0; i_repeat < repeat_number; i_repeat++)
			{
				try {
					Chunk chunk(block_size, chunk_size_in_bytes);
					std::size_t size = 100;

					for (std::size_t i = 0; i < size; i++)
					{
						if (chunk.HasAvailableBlock()) {
							void* p = chunk.Allocate();
							KRATOS_CHECK_NOT_EQUAL(p, nullptr);
						}
					}
				}
				catch (Exception& e) {
#pragma omp critical
					buffer << e;
					break;
				}
				catch (...) {
#pragma omp critical
					buffer << "UknownError";
					break;
				}
			}
			if (!buffer.str().empty())
				KRATOS_ERROR << buffer.str();

		}


		KRATOS_TEST_CASE_IN_SUITE(ChunkParallelAllocateDeallocate, EXCLUDED_KratosCoreFastSuite)
		{
			std::size_t block_size = 5;
			std::size_t chunk_size_in_bytes = 1024;

			auto repeat_number = 10;
			std::stringstream buffer;

#pragma omp parallel for 
			for (auto i_repeat = 0; i_repeat < repeat_number; i_repeat++)
			{
				try {
					Chunk chunk(block_size, chunk_size_in_bytes);
					std::size_t size = 100;
					std::vector<void*> the_vector(size);


					for (std::size_t i = 0; i < size; i++)
					{
						the_vector[i] = chunk.Allocate();
						KRATOS_CHECK_NOT_EQUAL(the_vector[i], nullptr);
					}

					for (std::size_t i = 0; i < size; i++)
						chunk.Deallocate(the_vector[i]);
				}
				catch (Exception& e) {
#pragma omp critical
					buffer << e;
					break;
				}
				catch (...) {
#pragma omp critical
					buffer << "UknownError";
					break;
				}
			}
			if (!buffer.str().empty())
				KRATOS_ERROR << buffer.str();

		}



	}
}  // namespace Kratos.


