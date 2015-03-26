template <typename T>
void AllocateGrid(T ** grid,
    const uint &X, const uint &Y, const uint &Z) {

  *grid = (T *)malloc(sizeof(T) * 
    (X+BW) * (Y+BW) * (Z+BW));
}

template <typename T>
void AllocateGrid(T * grid,
    const uint &X, const uint &Y, const uint &Z, const uint &align) {

  uint size = (X+BW) * (Y+BW) * (Z+BW);

  while((sizeof(T) * size) % align) size++;

  *grid = (T *)memalign(align,sizeof(T) * size);
}

template <typename T>
void ReleaseGrid(T * grid) {

  free(*grid);
}