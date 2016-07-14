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
			std::size_t block_size = 5;
			unsigned char number_of_blocks = 2;
			Chunk modified_chunk(block_size, number_of_blocks);
			KRATOS_CHECK_EQUAL(modified_chunk.GetNumberOfAvailableBlocks(), number_of_blocks);
		}

		KRATOS_TEST_CASE_IN_SUITE(ChunkAllocateDeallocate, KratosCoreFastSuite)
		{
			std::size_t block_size = 5; 
			unsigned char number_of_blocks = 127;
			Chunk chunk(block_size, number_of_blocks);

			chunk.Release();

			number_of_blocks = 255;
			chunk.Initialize(block_size, number_of_blocks);


			std::set<void *> pointer_set;
			while (chunk.HasAvailableBlock())
			{
				void* p = chunk.Allocate(block_size);
				KRATOS_ERROR_IF(pointer_set.find(p) != pointer_set.end()) << "The allocated position #" << pointer_set.size() + 1 << " was allocated before";
				pointer_set.insert(p);
			}
			KRATOS_CHECK_EQUAL(pointer_set.size(), number_of_blocks) << " (pointer set size : " << pointer_set.size() << ")" << std::endl;
			
			for (auto i_pointer = pointer_set.begin(); i_pointer != pointer_set.end(); i_pointer++) {
				chunk.Deallocate(*i_pointer, block_size);
			}
			KRATOS_CHECK_EQUAL(chunk.GetNumberOfAvailableBlocks(), number_of_blocks);
		}

	}
}  // namespace Kratos.


