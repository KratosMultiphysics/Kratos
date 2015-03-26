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
// #include "simd.h"
// #include "hacks.h"

#include "utils.h"

const double PI = 3.14159265;

uint N = 0;
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
        field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][1] = omega * (double)(i-(X+1.0)/2.0);
        field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][2] = 0.0;

        maxv = std::max((double)abs(field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][0]),maxv);
        maxv = std::max((double)abs(field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][1]),maxv);
      }
    }
  }

}

void InitializeCorrectionField(Triple * field,
    const uint &X, const uint &Y, const uint &Z) {

  for(uint k = 0; k < Z + BW; k++) {
    for(uint j = 0; j < Y + BW; j++) {
      for(uint i = 0; i < X + BW; i++) {
        field[k*(Y+BW)*(X+BW)+j*(X+BW)+i][0] = 0.5;
        field[k*(Y+BW)*(X+BW)+j*(X+BW)+i][1] = 0.5;
        field[k*(Y+BW)*(X+BW)+j*(X+BW)+i][2] = 0.5;
      } 
    }
  }

}

template <typename T>
void WriteHeatFocus(T * gridA,
    const uint &X, const uint &Y, const uint &Z) {

  uint Xc, Yc, Zc;

  Xc = 1.0/3.0*(X);
  Yc = 1.0/3.0*(Y);
  Zc = 1.0/2.0*(Z);

  for(uint k = 0; k < Z + BW; k++) {
    for(uint j = 0; j < Y + BW; j++) {
      for(uint i = 0; i < Z + BW; i++) {

        double d2 = pow((Xc - (double)(i)),2) + pow((Yc - (double)(j)),2) + pow((Zc - (double)(k)),2); 
        double rr = pow(X/6.0,2);  
        
        if(d2 < rr)
          gridA[k*(Y+BW)*(X+BW)+j*(X+BW)+i] = 1.0;// - d2/rr;
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

void calculateFormFactor(Triple prevDelta, Triple * correction, uint cell) {

  correction[cell][0] = 1-(prevDelta[0] - floor(prevDelta[0]));
  correction[cell][1] = 1-(prevDelta[1] - floor(prevDelta[1]));
  correction[cell][2] = 1-(prevDelta[2] - floor(prevDelta[2]));
}

template <typename T>
double bfceeInterpolationSteepWithCorrection(Triple prevDelta, T * gridA, Triple * dc,
    const uint &X, const uint &Y, const uint &Z, uint cell) {

  uint pi,pj,pk,ni,nj,nk;

  pi = floor(prevDelta[0]); ni = pi+1;
  pj = floor(prevDelta[1]); nj = pj+1;
  pk = floor(prevDelta[2]); nk = pk+1;

  double Nx, Ny, Nz;

  Nx = 1-(prevDelta[0] - floor(prevDelta[0]));
  Ny = 1-(prevDelta[1] - floor(prevDelta[1]));
  Nz = 1-(prevDelta[2] - floor(prevDelta[2]));

  uint index_a = pk*(Z+BW)*(Y+BW)+pj*(Y+BW)+pi;
  uint index_b = pk*(Z+BW)*(Y+BW)+pj*(Y+BW)+ni;
  uint index_c = pk*(Z+BW)*(Y+BW)+nj*(Y+BW)+pi;
  uint index_d = pk*(Z+BW)*(Y+BW)+nj*(Y+BW)+ni;
  uint index_e = nk*(Z+BW)*(Y+BW)+pj*(Y+BW)+pi;
  uint index_f = nk*(Z+BW)*(Y+BW)+pj*(Y+BW)+ni;
  uint index_g = nk*(Z+BW)*(Y+BW)+nj*(Y+BW)+pi;
  uint index_h = nk*(Z+BW)*(Y+BW)+nj*(Y+BW)+ni;

  double a_fact = 0.5;
  double b_fact = 0.5;

  return (
    gridA[index_a] * (    Nx * a_fact + dc[cell][0] * b_fact) * (    Ny * a_fact + dc[cell][1] * b_fact) * (    Nz * a_fact + dc[cell][2] * b_fact) +
    gridA[index_b] * (1 - Nx * a_fact - dc[cell][0] * b_fact) * (    Ny * a_fact + dc[cell][1] * b_fact) * (    Nz * a_fact + dc[cell][2] * b_fact) +
    gridA[index_c] * (    Nx * a_fact + dc[cell][0] * b_fact) * (1 - Ny * a_fact - dc[cell][1] * b_fact) * (    Nz * a_fact + dc[cell][2] * b_fact) +
    gridA[index_d] * (1 - Nx * a_fact - dc[cell][0] * b_fact) * (1 - Ny * a_fact - dc[cell][1] * b_fact) * (    Nz * a_fact + dc[cell][2] * b_fact) +
    gridA[index_e] * (    Nx * a_fact + dc[cell][0] * b_fact) * (    Ny * a_fact + dc[cell][1] * b_fact) * (1 - Nz * a_fact - dc[cell][2] * b_fact) +
    gridA[index_f] * (1 - Nx * a_fact - dc[cell][0] * b_fact) * (    Ny * a_fact + dc[cell][1] * b_fact) * (1 - Nz * a_fact - dc[cell][2] * b_fact) +
    gridA[index_g] * (    Nx * a_fact + dc[cell][0] * b_fact) * (1 - Ny * a_fact - dc[cell][1] * b_fact) * (1 - Nz * a_fact - dc[cell][2] * b_fact) +
    gridA[index_h] * (1 - Nx * a_fact - dc[cell][0] * b_fact) * (1 - Ny * a_fact - dc[cell][1] * b_fact) * (1 - Nz * a_fact - dc[cell][2] * b_fact)
  );

  // return (
  //   gridA[pk*(Z+BW)*(Y+BW)+pj*(Y+BW)+pi] * (    Nx) * (    Ny) * (    Nz) +
  //   gridA[pk*(Z+BW)*(Y+BW)+pj*(Y+BW)+ni] * (1 - Nx) * (    Ny) * (    Nz) +
  //   gridA[pk*(Z+BW)*(Y+BW)+nj*(Y+BW)+pi] * (    Nx) * (1 - Ny) * (    Nz) +
  //   gridA[pk*(Z+BW)*(Y+BW)+nj*(Y+BW)+ni] * (1 - Nx) * (1 - Ny) * (    Nz) +
  //   gridA[nk*(Z+BW)*(Y+BW)+pj*(Y+BW)+pi] * (    Nx) * (    Ny) * (1 - Nz) +
  //   gridA[nk*(Z+BW)*(Y+BW)+pj*(Y+BW)+ni] * (1 - Nx) * (    Ny) * (1 - Nz) +
  //   gridA[nk*(Z+BW)*(Y+BW)+nj*(Y+BW)+pi] * (    Nx) * (1 - Ny) * (1 - Nz) +
  //   gridA[nk*(Z+BW)*(Y+BW)+nj*(Y+BW)+ni] * (1 - Nx) * (1 - Ny) * (1 - Nz)
  // );

}

template <typename T, typename U>
void advection(T * gridA, T * gridB, U * fieldA, U * fieldB, Triple * dCorrection,
    const uint &X, const uint &Y, const uint &Z) {

  for(uint k = BWP + omp_get_thread_num(); k < Z + BWP; k+=omp_get_num_threads()) {
    for(uint j = BWP; j < Y + BWP; j++) {
      for(uint i = BWP; i < X + BWP; i++) {
        uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+i;

        Triple origin;
        Triple backward;

        origin[0] = i; origin[1] = j; origin[2] = k;

        for(int d = 0; d < 3; d++) {
          backward[d] = origin[d]*h-fieldA[cell][d]*dt;
        }

        static int firstime = 1;
        if(firstime) {
          firstime = 0;
          calculateFormFactor(backward,dCorrection,cell);
        }

        gridB[cell] = bfceeInterpolationSteepWithCorrection(backward,gridA,dCorrection,N,N,N,cell);

        calculateFormFactor(backward,dCorrection,cell);
      }
    }
  }

  for(uint k = BWP + omp_get_thread_num(); k < Z + BWP; k+=omp_get_num_threads()) {
    for(uint j = BWP; j < Y + BWP; j++) {
      for(uint i = BWP; i < X + BWP; i++) {
        uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+i;

        Triple origin;

        Triple backward;
        Triple backwardVel;
        Triple forward;

        origin[0] = i; origin[1] = j; origin[2] = k;

        for(int d = 0; d < 3; d++) {
          backward[d] = origin[d]*h-fieldA[cell][d]*dt;
        }

        bfceeInterpolationSteep(backward,backwardVel,fieldA,N,N,N);

        for(int d = 0; d < 3; d++) {
          forward[d] = backward[d] + backwardVel[d] * dt;
        }

        gridA[cell] = 3.0/2.0 * gridB[cell] - 0.5 * bfceeInterpolationSteep(forward,gridB,N,N,N);
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

  FileIO<PrecisionType> io;

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
  InitializeCorrectionField(form0,N,N,N);

  dt = CFL*h/maxv;

  std::stringstream name_mesh;
  std::stringstream name_post;

  name_mesh << "grid" << N << ".post.msh";
  name_post << "grid" << N << ".post.res";

  std::ofstream results(name_post.str().c_str());
  results << "GiD Post Results File 1.0" << std::endl << std::endl;

  io.WriteGidMesh(step0,N,N,N,name_mesh.str().c_str());
  printf("Begin...\n");
  for(int i = 0; i < steeps; i++) {
      advection(step0,step1,velf0,velf1,form0,N,N,N);
          
      io.WriteGidResults(step0,N,N,N,results,i);
      //difussion(step1,step2,N,N,N);
      //std::swap(step0,step1);
  }
  //WriteGridXYZ(velf0,N,N,N,"velc.dat");

  ReleaseGrid(&step0);
  ReleaseGrid(&step1);

  ReleaseGrid(&velf0);
  ReleaseGrid(&velf1);

  ReleaseGrid(&form0);

  printf("De-Allocation correct\n");
}