#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <memory>
#include <list>

// GiD IO
#include "gidpost/source/gidpost.h"

class GridPrinter {
public:
  // Creator & destructor

  GridPrinter(const double &dx, std::size_t * numCells, std::size_t * borderWidth);
  ~GridPrinter();

  void Initialize(const char * name, const std::size_t &N);

  void WriteGidMeshWithSkinBinary();
  void WriteGidMeshBinary();

  template<typename _Tp>
  void WriteGidResultsBinary1D(_Tp * grid, const int step, const char * name);

  template<typename _Tp>
  void WriteGidResultsBinary3D(_Tp * grid, const int step, const char * name);

private:

  void SetMeshName(const char * mesh, const std::size_t &N);
  void setPostName(const char * post, const std::size_t &N);

  std::stringstream name_mesh;
  std::stringstream name_post;
  std::stringstream name_raw;

  // double mDt;
  double mDx;
  // double mIdx;

  std::size_t * mNumCells;
  std::size_t * mBorderWidth;
};
