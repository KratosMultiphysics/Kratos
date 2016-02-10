#include "custom_utilities/grid_printer.h"

GridPrinter::GridPrinter(const double &dx, std::size_t * numCells, std::size_t * borderWidth)
: mDx(dx), mNumCells(numCells), mBorderWidth(borderWidth) {};

GridPrinter::~GridPrinter() {
  GiD_ClosePostResultFile();
};

void GridPrinter::Initialize(const char * name, const std::size_t &N) {
  SetMeshName(name,N);
  setPostName(name,N);

  name_raw << name;

  GiD_OpenPostMeshFile(name_mesh.str().c_str(), GiD_PostAscii);
  GiD_OpenPostResultFile(name_post.str().c_str(), GiD_PostAscii);
}

/** Writes the mesh in GiD format.
 * @X:        X-Size of the grid
 * @Y:        Y-Size of the grid
 * @Z:        Z-Size of the grid
 **/
void GridPrinter::WriteGidMeshWithSkinBinary() {

  int elemi[8];

  GiD_BeginMesh(name_raw.str().c_str(), GiD_3D, GiD_Hexahedra, 8);

  std::size_t gridLength[3] = {
    mNumCells[0] + mBorderWidth[0] * 2,
    mNumCells[1] + mBorderWidth[1] * 2,
    mNumCells[2] + mBorderWidth[2] * 2,
  };

  GiD_BeginCoordinates();
  for(std::size_t k = 0; k < gridLength[2]; k++) {
    for(std::size_t j = 0; j < gridLength[1]; j++) {
      for(std::size_t i = 0; i < gridLength[0]; i++) {
        std::size_t cell = k*(gridLength[2])*(gridLength[1])+j*(gridLength[1])+mBorderWidth[0]+i;
        GiD_WriteCoordinates(
          (int)cell,
          (double)1.0f-i*mDx,
          (double)1.0f-j*mDx,
          (double)1.0f-k*mDx
        );
      }
    }
  }
  GiD_EndCoordinates();

  GiD_BeginElements();
  for(std::size_t k = 0; k < gridLength[2] - 1; k++) {
    for(std::size_t j = 0; j < gridLength[1] - 1; j++) {
      for(std::size_t i = 0; i < gridLength[0] - 1; i++) {
        std::size_t cell = k*(gridLength[2])*(gridLength[1])+j*(gridLength[1])+1+i;
        elemi[0] = (int)(cell);
        elemi[1] = (int)(cell+1);
        elemi[2] = (int)(cell+1+(gridLength[1]));
        elemi[3] = (int)(cell+(gridLength[1]));
        elemi[4] = (int)(cell+(gridLength[2])*(gridLength[1]));
        elemi[5] = (int)(cell+1+(gridLength[2])*(gridLength[1]));
        elemi[6] = (int)(cell+1+(gridLength[2])*(gridLength[1])+(gridLength[1]));
        elemi[7] = (int)(cell+(gridLength[2])*(gridLength[1])+(gridLength[1]));

        GiD_WriteElement(
          (int)cell,
          elemi);
      }
    }
  }

  GiD_EndElements();
  GiD_EndMesh();
}


/** Writes the mesh in GiD format.
 * @X:        X-Size of the grid
 * @Y:        Y-Size of the grid
 * @Z:        Z-Size of the grid
 **/
void GridPrinter::WriteGidMeshBinary() {

  int elemi[8];

  std::size_t gridLength[3] = {
    mNumCells[0] + mBorderWidth[0] * 2,
    mNumCells[1] + mBorderWidth[1] * 2,
    mNumCells[2] + mBorderWidth[2] * 2,
  };

  GiD_BeginMesh(name_raw.str().c_str(), GiD_3D, GiD_Hexahedra, 8);

  GiD_BeginCoordinates();
  for(std::size_t k = 0; k < gridLength[2]; k++) {
    for(std::size_t j = 0; j < gridLength[1]; j++) {
      for(std::size_t i = 0; i < gridLength[0]; i++) {
        std::size_t cell = k*(gridLength[2])*(gridLength[1])+j*(gridLength[1])+mBorderWidth[0]+i;
        GiD_WriteCoordinates(
          (int)cell,
          (double)1.0f-i*mDx,
          (double)1.0f-j*mDx,
          (double)1.0f-k*mDx
        );
      }
    }
  }
  GiD_EndCoordinates();

  GiD_BeginElements();
  for(std::size_t k = mBorderWidth[2]; k < gridLength[2] - mBorderWidth[2] - 1; k++) {
    for(std::size_t j = mBorderWidth[1]; j < gridLength[1] - mBorderWidth[1] - 1; j++) {
      for(std::size_t i = mBorderWidth[0]; i < gridLength[0] - mBorderWidth[0] - 1; i++) {
        std::size_t cell = k*(gridLength[2])*(gridLength[2])+j*(gridLength[1])+mBorderWidth[0]+i;
        elemi[0] = (int)(cell);
        elemi[1] = (int)(cell+1);
        elemi[2] = (int)(cell+1+(gridLength[1]));
        elemi[3] = (int)(cell+(gridLength[1]));
        elemi[4] = (int)(cell+(gridLength[2])*(gridLength[1]));
        elemi[5] = (int)(cell+1+(gridLength[2])*(gridLength[2]));
        elemi[6] = (int)(cell+1+(gridLength[2])*(gridLength[1])+(gridLength[1]));
        elemi[7] = (int)(cell+(gridLength[2])*(gridLength[1])+(gridLength[1]));

        GiD_WriteElement(
          (int)cell,
          elemi);
      }
    }
  }

  GiD_EndElements();
  GiD_EndMesh();
}

/** Writes the results in GiD format.
 * @grid:     Value of the grid in Local or Global coordinatr system
 * @X:        X-Size of the grid
 * @Y:        Y-Size of the grid
 * @Z:        Z-Size of the grid
 * @step:     Step of the result
 **/
template<typename _Tp>
void GridPrinter::WriteGidResultsBinary1D(
    _Tp * grid,
    const int step,
    const char * name) {

  std::size_t gridLength[3] = {
    mNumCells[0] + mBorderWidth[0] * 2,
    mNumCells[1] + mBorderWidth[1] * 2,
    mNumCells[2] + mBorderWidth[2] * 2,
  };

  GiD_BeginResult(name,"Static",step,GiD_Scalar,GiD_OnNodes,NULL,NULL,0,NULL);

  for(std::size_t k = 0; k < gridLength[2]; k++) {
    for(std::size_t j = 0; j < gridLength[1]; j++) {
      for(std::size_t i = 0; i < gridLength[0]; i++) {

        std::size_t celln = k*(gridLength[2])*(gridLength[1])+j*(gridLength[1])+i;
        std::size_t cell = celln; //interleave64(i,j,k);

        GiD_WriteScalar(
          (int)(celln+1),
          (_Tp)grid[cell]);
      }
    }
  }
  GiD_EndResult();
  GiD_FlushPostFile();
}

/** Writes the results in GiD format.
 * @grid:     Value of the grid in Local or Global coordinatr system
 * @X:        X-Size of the grid
 * @Y:        Y-Size of the grid
 * @Z:        Z-Size of the grid
 * @step:     Step of the result
 **/
template<typename _Tp>
void GridPrinter::WriteGidResultsBinary3D(
    _Tp * grid,
    const int step,
    const char * name) {

  std::size_t gridLength[3] = {
    mNumCells[0] + mBorderWidth[0] * 2,
    mNumCells[1] + mBorderWidth[1] * 2,
    mNumCells[2] + mBorderWidth[2] * 2,
  };

  GiD_BeginResult(name,"Static",step,GiD_Vector,GiD_OnNodes,NULL,NULL,0,NULL);

  for(std::size_t k = 0; k < gridLength[2]; k++) {
    for(std::size_t j = 0; j < gridLength[1]; j++) {
      for(std::size_t i = 0; i < gridLength[0]; i++) {

        std::size_t celln = k*(gridLength[2])*(gridLength[1])+j*(gridLength[1])+i;
        std::size_t cell = celln; //interleave64(i,j,k);

        GiD_WriteVector(
          (int)(celln+1),
          (_Tp)grid[cell*3+0],
          (_Tp)grid[cell*3+1],
          (_Tp)grid[cell*3+2]);
      }
    }
  }
  GiD_EndResult();
  GiD_FlushPostFile();
}

/** Sets the name for the mesh file
 * @name:     name of the mesh file
 **/
void GridPrinter::SetMeshName(const char * mesh, const std::size_t &N) {
  name_mesh << mesh << N << ".post.msh";
}

/** Sets the name for the post file
 * @name:     name of the post file
 **/
void GridPrinter::setPostName(const char * post, const std::size_t &N) {
  name_post << post << N << ".post.res";
}

// Function instantiation
template void GridPrinter::WriteGidResultsBinary1D(int *, const int, const char *);
template void GridPrinter::WriteGidResultsBinary3D(int *, const int, const char *);
template void GridPrinter::WriteGidResultsBinary1D(float *, const int, const char *);
template void GridPrinter::WriteGidResultsBinary3D(float *, const int, const char *);
template void GridPrinter::WriteGidResultsBinary1D(double *, const int, const char *);
template void GridPrinter::WriteGidResultsBinary3D(double *, const int, const char *);
