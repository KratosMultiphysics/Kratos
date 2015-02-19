#include <iostream>
#include <iomanip>
#include <fstream>

#include <sys/time.h>
#include <sys/types.h>
#include <omp.h>
#include <malloc.h>

// Vectorial
#include <emmintrin.h>
#include <immintrin.h>

// Solver
#include "defines.h"
#include "stencils.h"
#include "simd.h"
#include "hacks.h"
#include "utils.h"

const uint NSTEPS = 100000;    //  Steps

uint N;
uint NB;
uint NE;

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

  for(uint k = 0; k < BWP; k++) {
    for(uint j = 0; j < (Y+BW); j++) {
      for(uint i = 0; i < (X+BW); i++ ) {
        gridA[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i] = 9.0;
        gridB[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i] = 9.0;
      }
    }
  }
}

template <typename T>
void Solve(T * &gridA, T * &gridB,
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

void CalculateMortonIndex(uint * IndexMatrix,
    const uint &X, const uint &Y, const uint &Z) {

  for(uint k = 0 + omp_get_thread_num(); k < Z; k+=omp_get_num_threads()) {
    for(uint j = 0; j < Y; j++) {
      for(uint i = 0; i < X; i++) {
        uint morton_index = interleave64(i,j,k);

        uint mk = morton_index / (X * Y) + BWP;
        uint mj = (( morton_index % (X * Y) ) / Y ) + BWP;
        uint mi = (( morton_index % (X * Y) ) % Y ) + BWP;

        IndexMatrix[(X*Y)*k+(Y)*j+i] = mk*(Z+BW)*(Y+BW)+mj*(Y+BW)+mi;
      }
    }
  }
}

template <typename T>
void SolveMorton(T * &gridA, T * &gridB, uint * IndexMatrix,
    const uint &X, const uint &Y, const uint &Z) {

  for(uint kk = 0 + omp_get_thread_num(); kk < NB; kk+=omp_get_num_threads()) {
    for(uint jj = 0; jj < NB; jj++) {
      for(uint ii = 0; ii < NB; ii++) {
        for(uint k = (kk * NE); k < ((kk+1) * NE); k++) {
          for(uint j = (jj * NE); j < ((jj+1) * NE); j++) {
            for(uint i = (ii * NE); i < ((ii+1) * NE); i++) {
              uint cell = IndexMatrix[(X*Y)*k+(Y)*j+i];
              stencilCross(gridA,gridB,cell,X,Y,Z);
            }
          }
        }
      }
    }
  }

}

void Solve(VectorType * &gridA, VectorType * &gridB, 
    const uint &X, const uint &Y, const uint &Z) {

  uint cell, pcellb, pcelle;

  PrecisionType __attribute__((aligned(ALING))) tmpa[VP];
  PrecisionType __attribute__((aligned(ALING))) tmpb[VP];

  for(uint kk = 0 + omp_get_thread_num(); kk < NB; kk+=omp_get_num_threads()) {
    for(uint jj = 0; jj < NB; jj++) {
      for(uint k = BWP + (kk * NE); k < BWP + ((kk+1) * NE); k++) {
        for(uint j = BWP + (jj * NE); j < BWP + ((jj+1) * NE); j++) {

          pcellb = k*(X+BW)*(X+BW)/VP+j*(X+BW)/VP+(BWP/VP);
          pcelle = k*(X+BW)*(X+BW)/VP+j*(X+BW)/VP+((X/VP) + BWP/VP - 1);

          VectorType * left  = &gridA[pcellb-1];
          VectorType * right = &gridA[pcellb+1];
          VectorType * bott  = &gridA[pcellb-(X+BW)/VP];
          VectorType * top   = &gridA[pcellb+(X+BW)/VP];
          VectorType * front = &gridA[pcellb-(X+BW)*(X+BW)/VP];
          VectorType * back  = &gridA[pcellb+(X+BW)*(X+BW)/VP];

          // Prefix
          cell = pcellb;
          VSTORE(tmpa,gridA[pcelle]);
          tmpb[0] = 0.0;
          for(uint p = 1; p < VP; p++) {
            tmpb[p] = tmpa[p-1];
          }

          left++;
          gridB[cell++] = VSTENCIL(VLOAD(tmpb),*right++,*top++,*bott++,*front++,*back++);

          // Body
          for(uint i = BWP/VP + 1; i < (X/VP) + BWP/VP - 1; i++) {
            gridB[cell++] = VSTENCIL(*left++,*right++,*top++,*bott++,*front++,*back++);
          }

          // Sufix
          cell = pcelle;
          VSTORE(tmpa,gridA[pcellb]);
          for(uint p = 1; p < VP; p++) {
            tmpb[p-1] = tmpa[p];
          }
          tmpb[VP-1] = 0.0;

          gridB[cell] = VSTENCIL(*left++,VLOAD(tmpb),*top++,*bott++,*front++,*back++);
        }
      }
    }
  }

}

void SolveBlock(VectorType * &gridA, VectorType * &gridB, 
    const uint &X, const uint &Y, const uint &Z) {

  uint cell, pcellb, pcelle;
  uint NVE = NE/VP;

  PrecisionType __attribute__((aligned(ALING))) tmpa[VP];
  PrecisionType __attribute__((aligned(ALING))) tmpb[VP];

  for(uint kk = 0 + omp_get_thread_num(); kk < NB; kk+=omp_get_num_threads()) {
    for(uint jj = 0; jj < NB; jj++) {
      for(uint ii = 0; ii < NB; ii++) {
        for(uint k = BWP + (kk * NE); k < BWP + ((kk+1) * NE); k++) {
          for(uint j = BWP + (jj * NE); j < BWP + ((jj+1) * NE); j++) {

            pcellb = k*(Y+BW)*(X+BW)/VP+j*(X+BW)/VP+(BWP/VP)+(ii*NVE);
            pcelle = k*(Y+BW)*(X+BW)/VP+j*(X+BW)/VP+(((ii+1)*NVE) + BWP/VP - 1);

            VectorType * left  = &gridA[pcellb-1];
            VectorType * right = &gridA[pcellb+1];
            VectorType * bott  = &gridA[pcellb-(X+BW)/VP];
            VectorType * top   = &gridA[pcellb+(X+BW)/VP];
            VectorType * front = &gridA[pcellb-(Y+BW)*(X+BW)/VP];
            VectorType * back  = &gridA[pcellb+(Y+BW)*(X+BW)/VP];

            cell = pcellb;
            if(ii == 0) {
              // Prefix
              VSTORE(tmpa,gridA[pcelle]);
              tmpb[0] = 0.0;
              for(uint p = 1; p < VP; p++) {
                tmpb[p] = tmpa[p-1];
              }

              left++;
              gridB[cell++] = VSTENCIL(VLOAD(tmpb),*right++,*top++,*bott++,*front++,*back++);
            } else {
              gridB[cell++] = VSTENCIL(*left++,*right++,*top++,*bott++,*front++,*back++);
            }

            // Body
            for(uint i = BWP/VP + (ii*NVE) + 1; i < BWP/VP + ((ii+1)*NVE) - 1; i++) {
              gridB[cell++] = VSTENCIL(*left++,*right++,*top++,*bott++,*front++,*back++);
            }

            cell = pcelle;
            if(ii == NB-1) {
              // Sufix
              VSTORE(tmpa,gridA[pcellb]);
              for(uint p = 1; p < VP; p++) {
                tmpb[p-1] = tmpa[p];
              }
              tmpb[VP-1] = 0.0;

              gridB[cell] = VSTENCIL(*left++,VLOAD(tmpb),*top++,*bott++,*front++,*back++);
            } else {
              gridB[cell] = VSTENCIL(*left++,*right++,*top++,*bott++,*front++,*back++);
            }
          }
        }
      }
    }
  }

}


int main(int argc, char *argv[]) {

  double douration_a = 0.0;
  double douration_b = 0.0;
  double douration_c = 0.0;

  std::cout << "Reading problem data..." << std::endl;

  N  = atoi(argv[1]);
  NB = atoi(argv[2]);
  NE = N / NB;

  std::cout << "N: \t" << N  << std::endl;
  std::cout << "NB:\t" << NB << std::endl;
  std::cout << "NE:\t" << NE << std::endl;

  PrecisionType * gridA = NULL;
  PrecisionType * gridB = NULL;

  AllocateGrid(&gridA,N,N,N,ALING);
  AllocateGrid(&gridB,N,N,N,ALING);

  uint * IndexMatrix = (uint *)malloc(sizeof(uint)*N*N*N);

  CalculateMortonIndex(IndexMatrix,N,N,N);

  struct timeval start, end;

  /***************************************************/

  // Initialize(gridA,gridB,N,N,N);

  // gettimeofday(&start, NULL);

  // #pragma omp parallel
  // for(uint step = 0; step < NSTEPS; step++) {
  //   Solve(gridA,gridB,N,N,N);
  //   #pragma omp single
  //     std::swap(gridA,gridB);
  // }

  // gettimeofday(&end, NULL);
  // douration_a = FETCHTIME

  /***************************************************/

  // Initialize(gridA,gridB,N,N,N);

  // gettimeofday(&start, NULL);

  // #pragma omp parallel
  // for(uint step = 0; step < NSTEPS; step++) {
  //   SolveMorton(gridA,gridB,IndexMatrix,N,N,N);
  //   #pragma omp single
  //     std::swap(gridA,gridB);
  // }

  // gettimeofday(&end, NULL);
  // douration_b = FETCHTIME

  /***************************************************/

  Initialize(gridA,gridB,N,N,N);

  gettimeofday(&start, NULL);

  VectorType* vectorGridA = (VectorType*)gridA;
  VectorType* vectorGridB = (VectorType*)gridB;

  #pragma omp parallel
  for(uint step = 0; step < NSTEPS/2; step++) {
    Solve(vectorGridA,vectorGridB,N,N,N);
    Solve(vectorGridB,vectorGridA,N,N,N);
  }

  gettimeofday(&end, NULL);
  douration_c = FETCHTIME

  /***************************************************/
  
  WriteGrid(gridA,N,N,N,"grid.dat");

  ReleaseGrid(&gridA);
  ReleaseGrid(&gridB); 

  std::cout << std::endl;

  std::cout << "OMP  time:\t" << douration_a << std::endl;
  std::cout << "FUUN time:\t" << douration_b << std::endl;
  std::cout << "SIMD time:\t" << douration_c << std::endl;

  std::cout << "SpeedUP:  \t" << douration_a/douration_c << std::endl;

  return 0;
}