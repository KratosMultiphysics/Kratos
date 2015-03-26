#include <iostream>
#include <iomanip>
#include <fstream>

#include <iostream>     // std::cout, std::right, std::endl
#include <iomanip>      // std::setw
#include <sstream>

#include <sys/time.h>
#include <sys/types.h>
#include <omp.h>
#include <math.h>
#include <malloc.h>

// Vectorial
#include <emmintrin.h>
#include <immintrin.h>

// Solver
#include "defines.h"
#include "stencils.h"
#include "file_io.h"

#include "utils.h"

const double PI = 3.14159265;

uint N = 0;
uint OutputStep = 0;

double dx = 0;
double dt = 0.1;
double h = 16;
double omega = 1;
double maxv = 0.0;
double CFL = 2;

typedef double Triple[3];

inline double angle(double deltaX, double deltaY) {
  return atan2(deltaY, deltaX);
}

template <typename T>
void Initialize(T * gridA, T * gridB, 
    const uint &X, const uint &Y, const uint &Z) {

  for(uint k = BWP; k < Z - BWP; k++) {
    for(uint j = BWP; j < Y - BWP; j++) {
      for(uint i = BWP; i < X - BWP; i++ ) {
        gridA[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i] = 0.0;
        gridB[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i] = 0.0;
      }
    }
  }

}

void InitializeVelocity(Triple * field,
    const uint &X, const uint &Y, const uint &Z) {

  for(uint k = 0; k < Z + BW; k++) {
    for(uint j = 0; j < Y + BW; j++) {
      for(uint i = 0; i < X + BW; i++) {
        field[k*(Y+BW)*(X+BW)+j*(X+BW)+i][0] = 0.0;
        field[k*(Y+BW)*(X+BW)+j*(X+BW)+i][1] = 0.0;
        field[k*(Y+BW)*(X+BW)+j*(X+BW)+i][2] = 0.0;
      } 
    }
  }

  for(uint k = BWP; k < Z + BWP; k++) {
    for(uint j = BWP; j < Y + BWP; j++) {
      for(uint i = BWP; i < X + BWP; i++ ) {
        field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][0] = -omega * (double)(j-(Y+1.0)/2.0);
        field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][1] =  omega * (double)(i-(X+1.0)/2.0);
        field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][2] =  0.0;

        maxv = std::max((double)abs(field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][0]),maxv);
        maxv = std::max((double)abs(field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][1]),maxv);
      }
    }
  }

}

template <typename T>
void WriteHeatFocus(T * gridA,
    const uint &X, const uint &Y, const uint &Z) {

  uint Xc, Yc, Zc;

  Xc = 2.0/5.0*(X);
  Yc = 2.0/5.5*(Y);
  Zc = 1.0/2.0*(Z);

  for(uint k = 0; k < Z + BW; k++) {
    for(uint j = 0; j < Y + BW; j++) {
      for(uint i = 0; i < Z + BW; i++) {

        double d2 = pow((Xc - (double)(i)),2) + pow((Yc - (double)(j)),2) + pow((Zc - (double)(k)),2); 
        double rr = pow(X/6.0,2);  
        
        if(d2 < rr)
          gridA[k*(Y+BW)*(X+BW)+j*(X+BW)+i] = 1.0 - d2/rr;
      }
    }
  }

}

void bfceeInterpolationSteep(Triple prevDelta, Triple prevVeloc, Triple * field,
    const uint &X, const uint &Y, const uint &Z) {

  uint pi,pj,pk,ni,nj,nk;

  pi = floor(prevDelta[0]); ni = pi+1;
  pj = floor(prevDelta[1]); nj = pj+1;
  pk = floor(prevDelta[2]); nk = pk+1;

  double Nx, Ny, Nz;

  Nx = 1-(prevDelta[0] - floor(prevDelta[0]));
  Ny = 1-(prevDelta[1] - floor(prevDelta[1]));
  Nz = 1-(prevDelta[2] - floor(prevDelta[2]));

  for(int d = 0; d < 3; d++) {
    prevVeloc[d] = (
      field[pk*(Z+BW)*(Y+BW)+pj*(Y+BW)+pi][d] * (    Nx) * (    Ny) * (    Nz) +
      field[pk*(Z+BW)*(Y+BW)+pj*(Y+BW)+ni][d] * (1 - Nx) * (    Ny) * (    Nz) +
      field[pk*(Z+BW)*(Y+BW)+nj*(Y+BW)+pi][d] * (    Nx) * (1 - Ny) * (    Nz) +
      field[pk*(Z+BW)*(Y+BW)+nj*(Y+BW)+ni][d] * (1 - Nx) * (1 - Ny) * (    Nz) +
      field[nk*(Z+BW)*(Y+BW)+pj*(Y+BW)+pi][d] * (    Nx) * (    Ny) * (1 - Nz) +
      field[nk*(Z+BW)*(Y+BW)+pj*(Y+BW)+ni][d] * (1 - Nx) * (    Ny) * (1 - Nz) +
      field[nk*(Z+BW)*(Y+BW)+nj*(Y+BW)+pi][d] * (    Nx) * (1 - Ny) * (1 - Nz) +
      field[nk*(Z+BW)*(Y+BW)+nj*(Y+BW)+ni][d] * (1 - Nx) * (1 - Ny) * (1 - Nz)
    );
  }
}

template <typename T>
double bfceeInterpolationSteep(Triple prevDelta, T * gridA,
    const uint &X, const uint &Y, const uint &Z) {

  uint pi,pj,pk,ni,nj,nk;

  pi = floor(prevDelta[0]); ni = pi+1;
  pj = floor(prevDelta[1]); nj = pj+1;
  pk = floor(prevDelta[2]); nk = pk+1;

  double Nx, Ny, Nz;

  Nx = 1-(prevDelta[0] - floor(prevDelta[0]));
  Ny = 1-(prevDelta[1] - floor(prevDelta[1]));
  Nz = 1-(prevDelta[2] - floor(prevDelta[2]));

  return (
    gridA[pk*(Z+BW)*(Y+BW)+pj*(Y+BW)+pi] * (    Nx) * (    Ny) * (    Nz) +
    gridA[pk*(Z+BW)*(Y+BW)+pj*(Y+BW)+ni] * (1 - Nx) * (    Ny) * (    Nz) +
    gridA[pk*(Z+BW)*(Y+BW)+nj*(Y+BW)+pi] * (    Nx) * (1 - Ny) * (    Nz) +
    gridA[pk*(Z+BW)*(Y+BW)+nj*(Y+BW)+ni] * (1 - Nx) * (1 - Ny) * (    Nz) +
    gridA[nk*(Z+BW)*(Y+BW)+pj*(Y+BW)+pi] * (    Nx) * (    Ny) * (1 - Nz) +
    gridA[nk*(Z+BW)*(Y+BW)+pj*(Y+BW)+ni] * (1 - Nx) * (    Ny) * (1 - Nz) +
    gridA[nk*(Z+BW)*(Y+BW)+nj*(Y+BW)+pi] * (    Nx) * (1 - Ny) * (1 - Nz) +
    gridA[nk*(Z+BW)*(Y+BW)+nj*(Y+BW)+ni] * (1 - Nx) * (1 - Ny) * (1 - Nz)
  );
}

template <typename T, typename U>
void advection(T * gridA, T * gridB, U * fieldA, U * fieldB, Triple * dCorrection,
    const uint &X, const uint &Y, const uint &Z) {

  // First loop. Estimate var(n+1)
  for(uint k = BWP + omp_get_thread_num(); k < Z + BWP; k+=omp_get_num_threads()) {
    for(uint j = BWP; j < Y + BWP; j++) {
      for(uint i = BWP; i < X + BWP; i++) {
        uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+i;

        Triple origin;
        Triple backward;
        Triple backwardVel;

        origin[0] = i;
        origin[1] = j;
        origin[2] = k;

        for(int d = 0; d < 3; d++) {
          backward[d] = origin[d]*h-fieldA[cell][d]*dt;
        }

        bfceeInterpolationSteep(backward,backwardVel,fieldA,N,N,N);

        for(int d = 0; d < 3; d++) {
          backward[d] = origin[d]*h-backwardVel[d]*dt;
        }

        gridB[cell] = bfceeInterpolationSteep(backward,gridA,N,N,N);
      }
    }
  }

  // First loop. Estimate var(n+1)
  for(uint k = BWP + omp_get_thread_num(); k < Z + BWP; k+=omp_get_num_threads()) {
    for(uint j = BWP; j < Y + BWP; j++) {
      for(uint i = BWP; i < X + BWP; i++) {
        uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+i;

        Triple origin;
        Triple forward;
        Triple forwardVel;

        origin[0] = i;
        origin[1] = j;
        origin[2] = k;

        for(int d = 0; d < 3; d++) {
          forward[d] = origin[d]*h+fieldA[cell][d]*dt;
        }

        bfceeInterpolationSteep(forward,forwardVel,fieldA,N,N,N);

        for(int d = 0; d < 3; d++) {
          forward[d] = origin[d]*h-forwardVel[d]*dt;
        }
 
        gridA[cell] = 1.5 * gridB[cell] - 0.5 * bfceeInterpolationSteep(forward,gridB,N,N,N);
      }
    } 
  }
}

template <typename T>
void difussion(T * &gridA, T * &gridB,
    const uint &X, const uint &Y, const uint &Z) {

  for(uint k = BWP + omp_get_thread_num(); k < Z + BWP; k+=omp_get_num_threads()) {
    for(uint j = BWP; j < Y + BWP; j++) {
      uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+BWP;
      for(uint i = BWP; i < X + BWP; i++) {
        stencilCross(gridA,gridB,cell++,X,Y,Z);
      } 
    }
  }
}

int main(int argc, char *argv[]) {

  N  = atoi(argv[1]);
  int steeps = atoi(argv[2]);
  h = 1;
  dx = h/N;

  FileIO<PrecisionType> io("grid",N);

  PrecisionType * step0 = NULL;
  PrecisionType * step1 = NULL;

  Triple * velf0 = NULL;
  Triple * velf1 = NULL;
  Triple * form0 = NULL;

  // Temperature
  AllocateGrid(&step0,N,N,N);
  AllocateGrid(&step1,N,N,N);

  // Velocity
  AllocateGrid(&velf0,N,N,N);
  AllocateGrid(&velf1,N,N,N);

  // Difusion correction
  AllocateGrid(&form0,N,N,N);

  printf("Allocation correct\n");

  printf("Initialize...\n");
  Initialize(step0,step1,N,N,N);
  WriteHeatFocus(step0,N,N,N);
  InitializeVelocity(velf0,N,N,N);

  dt = 0.25 * CFL*h/maxv;

  io.WriteGidMesh(step0,N,N,N);
  printf("Begin...\n");

  for(int i = 0; i < steeps; i++) {
      advection(step0,step1,velf0,velf1,form0,N,N,N);
          
      if(OutputStep == 0) {
        io.WriteGidResults(step0,N,N,N,i);
        OutputStep = 10;
      }
      //difussion(step1,step2,N,N,N);
      //std::swap(step0,step1);

      OutputStep--;
  }
  //WriteGridXYZ(velf0,N,N,N,"velc.dat");

  ReleaseGrid(&step0);
  ReleaseGrid(&step1);
  ReleaseGrid(&velf0);
  ReleaseGrid(&velf1);
  ReleaseGrid(&form0);

  printf("De-Allocation correct\n");
}