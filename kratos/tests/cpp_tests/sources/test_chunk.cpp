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
			std::size_t block_size_in_bytes = 5; // the aligned block size is 8
			std::size_t chunk_size_in_bytes = 1024;

			Chunk chunk(block_size_in_bytes, chunk_size_in_bytes);
			chunk.Initialize();

			std::size_t block_size_after_alignment = 8;
			std::size_t available_blocks_should_be = (chunk_size_in_bytes) / block_size_after_alignment;
			KRATOS_EXPECT_EQ(chunk.GetNumberOfAvailableBlocks(), available_blocks_should_be) << " Available block :" << chunk.GetNumberOfAvailableBlocks() << " vs " << available_blocks_should_be;
		}

		KRATOS_TEST_CASE_IN_SUITE(ChunkInitializeSmallBlock, KratosCoreFastSuite)
		{
			std::size_t block_size_in_bytes = 5; // the aligned block size is 8
			std::size_t chunk_size_in_bytes = 5;

			Chunk too_small_chunk(block_size_in_bytes, chunk_size_in_bytes);
			too_small_chunk.Initialize();

			KRATOS_EXPECT_EQ(too_small_chunk.GetNumberOfAvailableBlocks(), 0) << " Available block :" << too_small_chunk.GetNumberOfAvailableBlocks();
		}


		KRATOS_TEST_CASE_IN_SUITE(ChunkHasPointer, KratosCoreFastSuite)
		{
			std::size_t block_size_in_bytes = 5; // the aligned block size is 8
			std::size_t chunk_size_in_bytes = 1024;

			Chunk chunk(block_size_in_bytes, chunk_size_in_bytes);
			chunk.Initialize();

			std::size_t block_size_after_alignment = 8;
			std::size_t available_blocks_should_be = (chunk_size_in_bytes) / block_size_after_alignment;
			KRATOS_EXPECT_EQ(chunk.GetNumberOfAvailableBlocks(), available_blocks_should_be) << " Available block :" << chunk.GetNumberOfAvailableBlocks() << " vs " << available_blocks_should_be;

			for(std::size_t i = 0 ; i < available_blocks_should_be; i++) {
				const void* pointer = chunk.pGetData() + i*block_size_after_alignment / sizeof(Chunk::BlockType);
				KRATOS_EXPECT_EQ(chunk.Has(pointer), true) << chunk << " should have the pointer " << pointer << std::endl;
			}
			for(std::size_t i = available_blocks_should_be ; i < available_blocks_should_be + 5; i++) {
				const void* pointer = chunk.pGetData() + i*block_size_after_alignment / sizeof(Chunk::BlockType);
				KRATOS_EXPECT_EQ(chunk.Has(pointer), false) << chunk << " should Not have the pointer " << pointer << std::endl;
			}
		}

		KRATOS_TEST_CASE_IN_SUITE(ChunkParallelInitialize, KratosCoreFastSuite)
		{
			std::size_t block_size_in_bytes = 5; // the aligned block size is 8
			std::size_t chunk_size_in_bytes = 1024;

			std::size_t block_size_after_alignment = 8;
			std::size_t available_blocks_should_be = (chunk_size_in_bytes) / block_size_after_alignment;

			auto repeat_number = 10;
			#pragma omp parallel for
			for (auto i_repeat = 0; i_repeat < repeat_number; i_repeat++)
			{
				Chunk chunk(block_size_in_bytes, chunk_size_in_bytes);
				chunk.Initialize();

				KRATOS_EXPECT_EQ(chunk.GetNumberOfAvailableBlocks(), available_blocks_should_be) << " Available block :" << chunk.GetNumberOfAvailableBlocks() << " vs " << available_blocks_should_be;
			}
		}

		KRATOS_TEST_CASE_IN_SUITE(ChunkParallelInitializeSmallBlock, KratosCoreFastSuite)
		{
			std::size_t block_size_in_bytes = 5; // the aligned block size is 8
			std::size_t chunk_size_in_bytes = 5;

			auto repeat_number = 10;
			#pragma omp parallel for
			for (auto i_repeat = 0; i_repeat < repeat_number; i_repeat++)
			{
				Chunk too_small_chunk(block_size_in_bytes, chunk_size_in_bytes);
				too_small_chunk.Initialize();

				KRATOS_EXPECT_EQ(too_small_chunk.GetNumberOfAvailableBlocks(), 0) << " Available block :" << too_small_chunk.GetNumberOfAvailableBlocks();
			}
		}

		KRATOS_TEST_CASE_IN_SUITE(ChunkAllocateDeallocate, KratosCoreFastSuite)
		{
			std::size_t block_size_in_bytes = 5;
		    std::size_t chunk_size_in_bytes =  1024;
			Chunk chunk(block_size_in_bytes, chunk_size_in_bytes);

			chunk.Initialize();

			std::set<void *> pointer_set;
			while (chunk.HasAvailableBlock())
			{
				void* p = chunk.Allocate();
				KRATOS_ERROR_IF(pointer_set.find(p) != pointer_set.end()) << "The allocated position #" << pointer_set.size() + 1 << " was allocated before";
				pointer_set.insert(p);
			}
			KRATOS_EXPECT_EQ(pointer_set.size(), chunk.GetNumberOfBlocks()) << " (pointer set size : " << pointer_set.size() << ")" << std::endl;

			for (auto i_pointer = pointer_set.begin(); i_pointer != pointer_set.end(); i_pointer++) {
				chunk.Deallocate(*i_pointer);
			}
			KRATOS_EXPECT_EQ(chunk.GetNumberOfAvailableBlocks(), chunk.GetNumberOfBlocks());
		}

		KRATOS_TEST_CASE_IN_SUITE(ChunkParallelAllocate, KratosCoreFastSuite)
		{
			std::size_t block_size_in_bytes = 5;
		  	std::size_t chunk_size_in_bytes =  1024;

			auto repeat_number = 10;
			std::stringstream buffer;
			
			#pragma omp parallel for
			for (auto i_repeat = 0; i_repeat < repeat_number; i_repeat++)
			{
				try {
					Chunk chunk(block_size_in_bytes, chunk_size_in_bytes);
					chunk.Initialize();
					std::size_t size = 100;

					for (std::size_t i = 0; i < size; i++)
					{
						if (chunk.HasAvailableBlock()) {
							void* p = chunk.Allocate();
							KRATOS_EXPECT_NE(p, nullptr);
						}
					}
				}
				catch (Exception& e) {

					buffer << e;
				}
				catch (...) {

					buffer << "UknownError";
				}
			}
			if (!buffer.str().empty())
				KRATOS_ERROR << buffer.str();

		}


		KRATOS_TEST_CASE_IN_SUITE(ChunkParallelAllocateDeallocate, EXCLUDED_KratosCoreFastSuite)
		{
			std::size_t block_size_in_bytes = 5;
		  	std::size_t chunk_size_in_bytes =  1024;

			auto repeat_number = 100;
			std::stringstream buffer;

			#pragma omp parallel for
			for (auto i_repeat = 0; i_repeat < repeat_number; i_repeat++)
			{
				try {
					Chunk chunk(block_size_in_bytes, chunk_size_in_bytes);
					chunk.Initialize();
					std::size_t size = 100;
					std::vector<void*> the_vector(size);


					for (std::size_t i = 0; i < size; i++)
					{
							the_vector[i] = chunk.Allocate();
							KRATOS_EXPECT_NE(the_vector[i], nullptr) << " for i = " << i << " and repeat = " << i_repeat;

					}

					for (std::size_t i = 0; i < size; i++)
						chunk.Deallocate(the_vector[i]);
				}
				catch (Exception& e) {
#pragma omp critical
					buffer << e;

				}
				catch (...) {
#pragma omp critical
					buffer << "UknownError";

				}
			}
			if (!buffer.str().empty())
				KRATOS_ERROR << buffer.str();

		}



	}
}  // namespace Kratos.
