#include "custom_utilities/bfecc.h"

BfeccConvecter::BfeccConvecter() {}

BfeccConvecter::BfeccConvecter(int * flags, std::size_t * numCells, std::size_t * borderWidth, double dt, double dx, std::size_t dim)
  : mDt(dt), mDx(dx), mIdx(1.0f/dx), mDim(dim), mNumCells(numCells), mBorderWidth(borderWidth), mFlags(flags) {
}

BfeccConvecter::~BfeccConvecter() {
}

inline std::size_t BfeccConvecter::GetIndex(const std::size_t & i, const std::size_t & j, const std::size_t & k) {
  return k * (mNumCells[2] + mBorderWidth[2] * 2) * (mNumCells[1] + mBorderWidth[1] * 2) + j * (mNumCells[1] + mBorderWidth[1] * 2) + i;
}

void BfeccConvecter::Interpolate(const double & Idx,double * input,double * output, double * Coords, const std::size_t &Dim) {
  uint pi,pj,pk,ni,nj,nk;

  // This may cause errors
  for(std::size_t i = 0; i < Dimension; i++) {
    Coords[i] = Coords[i] < 0.0f ? 0.0f : Coords[i] > 1.0f ? 1.0f : Coords[i];
  }

  BfeccUtils::GlobalToLocal(Coords,mIdx,Dimension);

  pi = (uint)(Coords[0]); ni = pi+1;
  pj = (uint)(Coords[1]); nj = pj+1;
  pk = (uint)(Coords[2]); nk = pk+1;

  double Nx, Ny, Nz;

  Nx = 1.0f-(Coords[0] - (double)pi);
  Ny = 1.0f-(Coords[1] - (double)pj);
  Nz = 1.0f-(Coords[2] - (double)pk);

  for(size_t d = 0; d < Dim; d++) {
    *(output+d) = (
      input[GetIndex(pi,pj,pk)*Dim+d] * (       Nx) * (       Ny) * (       Nz) +
      input[GetIndex(ni,pj,pk)*Dim+d] * (1.0f - Nx) * (       Ny) * (       Nz) +
      input[GetIndex(pi,nj,pk)*Dim+d] * (       Nx) * (1.0f - Ny) * (       Nz) +
      input[GetIndex(ni,nj,pk)*Dim+d] * (1.0f - Nx) * (1.0f - Ny) * (       Nz) +
      input[GetIndex(pi,pj,nk)*Dim+d] * (       Nx) * (       Ny) * (1.0f - Nz) +
      input[GetIndex(ni,pj,nk)*Dim+d] * (1.0f - Nx) * (       Ny) * (1.0f - Nz) +
      input[GetIndex(pi,nj,nk)*Dim+d] * (       Nx) * (1.0f - Ny) * (1.0f - Nz) +
      input[GetIndex(ni,nj,nk)*Dim+d] * (1.0f - Nx) * (1.0f - Ny) * (1.0f - Nz)
    );
  }
}

void BfeccConvecter::Convect(double * const input, double * aux, double * velocity, double * output) {

  #pragma omp parallel for
  for(std::size_t k = mBorderWidth[2]; k < mNumCells[2] + mBorderWidth[2]; k++) {
    for(std::size_t j = mBorderWidth[1]; j < mNumCells[1] + mBorderWidth[1]; j++) {
      for(std::size_t i = mBorderWidth[0]; i < mNumCells[0] + mBorderWidth[0]; i++) {
        ApplyBack(output,input,velocity,i,j,k);
      }
    }
  }

  #pragma omp parallel for
  for(std::size_t k = mBorderWidth[2]; k < mNumCells[2] + mBorderWidth[2]; k++) {
    for(std::size_t j = mBorderWidth[1]; j < mNumCells[1] + mBorderWidth[1]; j++) {
      for(std::size_t i = mBorderWidth[0]; i < mNumCells[0] + mBorderWidth[0]; i++) {
        ApplyForth(aux,output,input,velocity,i,j,k);
      }
    }
  }

  #pragma omp parallel for
  for(std::size_t k = mBorderWidth[2]; k < mNumCells[2] + mBorderWidth[2]; k++) {
    for(std::size_t j = mBorderWidth[1]; j < mNumCells[1] + mBorderWidth[1]; j++) {
      for(std::size_t i = mBorderWidth[0]; i < mNumCells[0] + mBorderWidth[0]; i++) {
        ApplyEcc(output,aux,velocity,i,j,k);
      }
    }
  }
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
