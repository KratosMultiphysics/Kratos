---
type: C++ Method
title: Chunk
resource: includes/chunk.h#L72-L80
generated:
  by: okf-rs/0.7.0
---

# Signature

`Chunk(std::size_t BlockSizeInBytes, SizeType SizeInBytes) noexcept : mpData(nullptr) , mpEnd(nullptr) , mpUninitializedMemory(nullptr) , mSize(SizeInBytes) , mBlockSizeInBytes(BlockSizeInBytes) , mNumberOfAvailableBlocks(0) // to be initialized in initialize , mFirstAvailableBlock(mpData) , mBlockSizeAfterAlignment((BlockSizeInBytes + sizeof(BlockType) - 1) / sizeof(BlockType))`