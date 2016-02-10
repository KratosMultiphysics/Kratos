#if !defined(KRATOS_BFECC_CONVECTER)
#define KRATOS_BFECC_CONVECTER

#include <cstddef>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "custom_utilities/bfecc_utils.h"

constexpr std::size_t Dimension = 3;

class BfeccConvecter {
public:

  BfeccConvecter();
  BfeccConvecter(int * flags, std::size_t * numCells, std::size_t * borderWidth, double dt, double dx, std::size_t dim);

  ~BfeccConvecter();

  void Convect(double * const input, double * aux, double * velocity, double * output);

private:
  void ApplyBack(double * output, double * input, double * velocity, const std::size_t &i, const std::size_t &j, const std::size_t &k);
  void ApplyForth(double * output, double * convectedBack, double * original, double * velocity, const std::size_t &i, const std::size_t &j, const std::size_t &k);
  void ApplyEcc(double * output, double * input, double * velocity, const std::size_t &i, const std::size_t &j, const std::size_t &k);

  std::size_t GetIndex(const std::size_t & i, const std::size_t & j, const std::size_t & k);

  void Interpolate(const double & Idx,double * input,double * output, double * Coords, const std::size_t &Dim);

  double mDt;
  double mDx;
  double mIdx;

  std::size_t mDim;

  std::size_t * mNumCells;
  std::size_t * mBorderWidth;
  int * mFlags;

};

#endif // KRATOS_BFECC_UTILS defined
