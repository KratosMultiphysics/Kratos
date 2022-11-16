import dislib as ds
from dislib.data.array import Array

def load_blocks_array(blocks, shape, block_size):
    if shape[0] < block_size[0] or  shape[1] < block_size[1]:
        raise ValueError("The block size is greater than the ds-array")
    return Array(blocks, shape=shape, top_left_shape=block_size,
                     reg_shape=block_size, sparse=False)

def load_blocks_rechunk(blocks, shape, block_size, new_block_size):
    if shape[0] < new_block_size[0] or  shape[1] < new_block_size[1]:
        raise ValueError("The block size requested for rechunk"
                         "is greater than the ds-array")
    final_blocks = [[]]
    # Este bucle lo puse por si los Future objects se guardan en una lista, en caso de que la forma de guardarlos cambie, también cambiará un poco este bucle.
    # Si blocks se pasa ya como (p. ej) [[Future_object, Future_object]] no hace falta.
    for block in blocks:
        final_blocks[0].append(block)
    arr = load_blocks_array(final_blocks, shape, block_size)
    return arr.rechunk(new_block_size)
