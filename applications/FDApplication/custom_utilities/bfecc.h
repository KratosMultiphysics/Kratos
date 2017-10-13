#if !defined(KRATOS_BFECC_CONVECTER)
#define KRATOS_BFECC_CONVECTER

#include <stdlib.h>

#include <cstddef>
#include <string>
#include <vector>
#include <iostream>
#include <functional>

#include "custom_utilities/bfecc_utils.h"

constexpr std::size_t Dimension = 3;

class BfeccConvecter {
public:
  BfeccConvecter();
  BfeccConvecter(int * flags, std::vector<std::size_t> numCells, std::vector<std::size_t> borderWidth, double dt, double dx, std::size_t dim);

  ~BfeccConvecter();

  void Convect(double * const input, double * aux, double * velocity, double * output);

private:

  template <typename func_t>
  inline void OmpParallelLoop(func_t StepOperator) {
    #pragma omp parallel for
    for(std::size_t k = mBorderWidth[2]; k < mNumCells[2] + mBorderWidth[2]; k++) {
      for(std::size_t j = mBorderWidth[1]; j < mNumCells[1] + mBorderWidth[1]; j++) {
        for(std::size_t i = mBorderWidth[0]; i < mNumCells[0] + mBorderWidth[0]; i++) {
          StepOperator(i,j,k);
        }
      }
    }
  };

  template <typename func_t>
  inline void OmpParallelBlockLoop(func_t StepOperator, std::vector<std::size_t> NumBlocks) {

    #define BOT(_A,_D) std::max(mBorderWidth[_D],(_A * mCellsPerBlock[_D]))
    #define TOP(_A,_D) mBorderWidth[_D] + std::min(mCellsPerBlock[_D]*NumBlocks[_D]-mBorderWidth[_D],((_A+1) * mCellsPerBlock[_D]))

    std::size_t mCellsPerBlock[3] {
      mNumCells[0] / NumBlocks[0],
      mNumCells[1] / NumBlocks[1],
      mNumCells[2] / NumBlocks[2]
    };

    #pragma omp parallel for
    for(std::size_t kk = 0; kk < NumBlocks[2]; kk++ ) {
      for(std::size_t jj = 0; jj < NumBlocks[1]; jj++ ) {
        for(std::size_t ii = 0; ii < NumBlocks[0]; ii++ ) {
          for(std::size_t k = BOT(kk,2); k < TOP(kk,2); k++) {
            for(std::size_t j = BOT(jj,1); j < TOP(jj,1); j++) {
              for(std::size_t i = BOT(ii,0); i < TOP(ii,0); i++) {
                StepOperator(i,j,k);
              }
            }
          }
        }
      }
    }

    #undef BOT
    #undef TOP
  };

  void ApplyBack( double * output, double * input, double * velocity, const std::size_t &i, const std::size_t &j, const std::size_t &k);
  void ApplyForth(double * output, double * convectedBack, double * original, double * velocity, const std::size_t &i, const std::size_t &j, const std::size_t &k);
  void ApplyEcc(  double * output, double * input, double * velocity, const std::size_t &i, const std::size_t &j, const std::size_t &k);

  std::size_t GetIndex(const std::size_t & i, const std::size_t & j, const std::size_t & k);

  void Interpolate(const double & Idx,double * input,double * output, double * Coords, const std::size_t &Dim);

  double mDt;
  double mDx;
  double mIdx;

  std::size_t mDim;

  std::vector<std::size_t> mNumCells;
  std::vector<std::size_t> mBorderWidth;
  int * mFlags;

};

#endif // KRATOS_BFECC_UTILS defined
