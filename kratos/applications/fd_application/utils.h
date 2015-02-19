template <typename T>
void AllocateGrid(T ** grid,
    const uint &X, const uint &Y, const uint &Z) {

  *grid = (T *)malloc(sizeof(T) * 
    (X+BW) * (Y+BW) * (Z+BW));
}

template <typename T>
void AllocateGrid(T ** grid,
    const uint &X, const uint &Y, const uint &Z, const uint &align) {

  uint size = (X+BW) * (Y+BW) * (Z+BW);

  while((sizeof(T) * size) % align) size++;

  *grid = (T *)memalign(align,sizeof(T) * size);
}

template <typename T>
void ReleaseGrid(T ** grid) {

  free(*grid);
}

template <typename T>
void WriteGrid(T * grid, 
    const uint &X, const uint &Y, const uint &Z,
    const char * fileName) {

  std::ofstream outputFile(fileName,std::ofstream::app);

  outputFile << std::fixed;
  outputFile << std::setprecision(2);

  for(uint k = 0; k < Z + BW; k++) {
    for(uint j = 0; j < Y + BW; j++) {
      for(uint i = 0; i < X + BW; i++) {
        outputFile << grid[k*(Y+BW)*(X+BW)+j*(X+BW)+i] << " ";
      } 
      outputFile << std::endl;
    }
    outputFile << std::endl;
  }
  outputFile << "%" << std::endl;
  
}

template <typename T>
void WriteGridWipe(T * grid, 
    const uint &X, const uint &Y, const uint &Z,
    const char * fileName) {

  std::ofstream outputFile(fileName);

  outputFile << std::fixed;
  outputFile << std::setprecision(2);

  for(uint k = 0; k < Z + BW; k++) {
    for(uint j = 0; j < Y + BW; j++) {
      for(uint i = 0; i < X + BW; i++) {
        outputFile << grid[k*(Y+BW)*(X+BW)+j*(X+BW)+i] << " ";
      } 
      outputFile << std::endl;
    }
    outputFile << std::endl;
  }
  outputFile << "%" << std::endl;
  
}