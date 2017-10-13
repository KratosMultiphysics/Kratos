#include "custom_utilities/bfecc.h"

BfeccConvecter::BfeccConvecter() {}

BfeccConvecter::BfeccConvecter(int * flags, std::vector<std::size_t> numCells, std::vector<std::size_t> borderWidth, double dt, double dx, std::size_t dim)
  : mDt(dt), mDx(dx), mIdx(1.0f/dx), mDim(dim), mNumCells(numCells), mBorderWidth(borderWidth), mFlags(flags) {
}

BfeccConvecter::~BfeccConvecter() {
}

inline std::size_t BfeccConvecter::GetIndex(const std::size_t & i, const std::size_t & j, const std::size_t & k) {
  return k * (mNumCells[2] + mBorderWidth[2] * 2) * (mNumCells[1] + mBorderWidth[1] * 2) + j * (mNumCells[1] + mBorderWidth[1] * 2) + i;
}

void BfeccConvecter::Interpolate(const double & Idx,double * input,double * output, double * Coords, const std::size_t &Dim) {
  uint lo_i,lo_j,lo_k,hi_i,hi_j,hi_k;

  // This may cause errors
  for(std::size_t i = 0; i < Dimension; i++) {
    Coords[i] = Coords[i] < 0.0f ? 0.0f : Coords[i] > 1.0f ? 1.0f : Coords[i];
  }

  BfeccUtils::GlobalToLocal(Coords,mIdx,Dimension);

  lo_i = (uint)(Coords[0]); hi_i = lo_i+1;
  lo_j = (uint)(Coords[1]); hi_j = lo_j+1;
  lo_k = (uint)(Coords[2]); hi_k = lo_k+1;

  double Nx, Ny, Nz;

  Nx = 1.0f-(Coords[0] - (double)lo_i);
  Ny = 1.0f-(Coords[1] - (double)lo_j);
  Nz = 1.0f-(Coords[2] - (double)lo_k);

  double aa = (       Nx) * (       Ny);
  double ba = (1.0f - Nx) * (       Ny);
  double ab = (       Nx) * (1.0f - Ny);
  double bb = (1.0f - Nx) * (1.0f - Ny);

  auto padd = mNumCells[2] + mBorderWidth[2] * 2;

  auto cell_A = GetIndex(lo_i,lo_j,lo_k);
  auto cell_B = cell_A + padd;
  for(size_t d = 0; d < Dim; d++) {
    double &A = input[cell_A*Dim+d];
    double &B = input[cell_B*Dim+d];
    *(output+d)  = aa*(Nz*(A-B)+B);
  }

  auto cell_C = GetIndex(hi_i,lo_j,lo_k);
  auto cell_D = cell_C + padd;
  for(size_t d = 0; d < Dim; d++) {
    double &C = input[cell_C*Dim+d];
    double &D = input[cell_D*Dim+d];
    *(output+d) += ba*(Nz*(C-D)+D);
  }

  auto cell_E = GetIndex(lo_i,hi_j,lo_k);
  auto cell_F = cell_E + padd;
  for(size_t d = 0; d < Dim; d++) {
    double &E = input[cell_E*Dim+d];
    double &F = input[cell_F*Dim+d];
    *(output+d) += ab*(Nz*(E-F)+F);
  }

  auto cell_G = GetIndex(hi_i,hi_j,lo_k);
  auto cell_H = cell_G + padd;
  for(size_t d = 0; d < Dim; d++) {
    double &G = input[cell_H*Dim+d];
    double &H = input[cell_H*Dim+d];
    *(output+d) += bb*(Nz*(G-H)+H);
  }
}

void BfeccConvecter::Convect(double * const input, double * aux, double * velocity, double * output) {

  auto ApplyBK = std::bind(&BfeccConvecter::ApplyBack,  this, output, input, velocity, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
  auto ApplyFW = std::bind(&BfeccConvecter::ApplyForth, this, aux, output, input, velocity, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
  auto ApplyEC = std::bind(&BfeccConvecter::ApplyEcc,   this, output, aux, velocity, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

  OmpParallelLoop(ApplyBK);
  OmpParallelLoop(ApplyFW);
  OmpParallelLoop(ApplyEC);

  // OmpParallelBlockLoop(ApplyBK, {4, 4, 4});
  // OmpParallelBlockLoop(ApplyFW, {4, 4, 4});
  // OmpParallelBlockLoop(ApplyEC, {4, 4, 4});
}

void BfeccConvecter::ApplyBack(
    double * output, double * input,
    double * velocity,
    const std::size_t &i, const std::size_t &j, const std::size_t &k) {

  auto cell = GetIndex(i,j,k);

  double convectedBack[Dimension];
  double origin[Dimension] = {i * mDx, j * mDx, k * mDx}; // May be source of errors while casting
  double displacement[Dimension];

  if( k > mBorderWidth[2] && k < mNumCells[2] + mBorderWidth[2] - 1 &&
      j > mBorderWidth[1] && j < mNumCells[1] + mBorderWidth[1] - 1 &&
      i > mBorderWidth[0] && i < mNumCells[0] + mBorderWidth[0] - 1 ) {

    for(std::size_t d = 0; d < 3; d++) {
      displacement[d] = origin[d] - velocity[cell*mDim+d] * mDt;
      if(displacement[d] < 0 || displacement[d] > mNumCells[d]) {
        mFlags[cell] |= OUT_OF_BOUNDS;
      }
    }

    if(!(mFlags[cell] & OUT_OF_BOUNDS)) {
      Interpolate(mIdx,input,(double *)convectedBack,displacement,mDim);

      for(std::size_t d = 0; d < 3; d++) {
        if(!(mFlags[cell]))
          output[cell*mDim+d] = convectedBack[d];
      }
    }
  } else {
    for(std::size_t d = 0; d < 3; d++) {
      output[cell*mDim+d] = input[cell*mDim+d];
    }
  }
}

void BfeccConvecter::ApplyForth(
      double * output, double * convectedBack, double * original,
      double * velocity,
      const std::size_t &i, const std::size_t &j, const std::size_t &k) {

  auto cell = GetIndex(i,j,k);

  double convectedForth[Dimension];
  double origin[Dimension] = {i * mDx, j * mDx, k * mDx}; // May be source of errors while casting
  double displacement[Dimension];

  if( k > mBorderWidth[2] && k < mNumCells[2] + mBorderWidth[2] - 1 &&
      j > mBorderWidth[1] && j < mNumCells[1] + mBorderWidth[1] - 1 &&
      i > mBorderWidth[0] && i < mNumCells[0] + mBorderWidth[0] - 1 ) {

    for(std::size_t d = 0; d < 3; d++) {
      displacement[d] = origin[d] + velocity[cell*mDim+d] * mDt;
      if(displacement[d] < 0 || displacement[d] > mNumCells[d]) {
        mFlags[cell] |= OUT_OF_BOUNDS;
      }
    }

    if(!(mFlags[cell] & OUT_OF_BOUNDS)) {
      Interpolate(mIdx,convectedBack,(double*)convectedForth,displacement,mDim);

      for(std::size_t d = 0; d < 3; d++) {
        if(!(mFlags[cell]))
          output[cell*mDim+d] = 1.5f * original[cell*mDim+d] - 0.5f * convectedForth[d];
      }
    }
  } else {
    for(std::size_t d = 0; d < 3; d++) {
      output[cell*mDim+d] = original[cell*mDim+d];
    }
  }
}

void BfeccConvecter::ApplyEcc(
    double * output, double * input,
    double * velocity,
    const std::size_t &i, const std::size_t &j, const std::size_t &k) {

  auto cell = GetIndex(i,j,k);

  double corrected[Dimension];
  double origin[Dimension] = {i * mDx, j * mDx, k * mDx}; // May be source of errors while casting
  double displacement[Dimension];

  if( k > mBorderWidth[2] && k < mNumCells[2] + mBorderWidth[2] - 1 &&
      j > mBorderWidth[1] && j < mNumCells[1] + mBorderWidth[1] - 1 &&
      i > mBorderWidth[0] && i < mNumCells[0] + mBorderWidth[0] - 1 ) {

    for(std::size_t d = 0; d < 3; d++) {
      displacement[d] = origin[d] - velocity[cell*mDim+d] * mDt;
      if(displacement[d] < 0 || displacement[d] > mNumCells[d]) {
        mFlags[cell] |= OUT_OF_BOUNDS;
      }
    }

    if(!(mFlags[cell] & OUT_OF_BOUNDS)) {
      Interpolate(mIdx,input,(double*)corrected,displacement,mDim);

      for(std::size_t d = 0; d < 3; d++) {
        if(!(mFlags[cell]))
          output[cell*mDim+d] = corrected[d];
      }
    }

    for(std::size_t d = 0; d < 3; d++) {
      mFlags[cell] &= ~OUT_OF_BOUNDS;
    }
  } else {
    for(std::size_t d = 0; d < 3; d++) {
      output[cell*mDim+d] = input[cell*mDim+d];
    }
  }
}
