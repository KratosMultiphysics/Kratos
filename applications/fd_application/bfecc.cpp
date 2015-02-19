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
        field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][0] = -omega * (double)(j-(Y+1)/2.0);
        field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][1] = omega * (double)(i-(X+1)/2.0);
        field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][2] = 0.0;

    		maxv = std::max((double)abs(field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][0]),maxv);
    		maxv = std::max((double)abs(field[k*(Z+BW)*(Y+BW)+j*(Y+BW)+i][1]),maxv);
      }
    }
  }

}

template <typename T>
void WriteHeatFocus(T * gridA,
    const uint &X, const uint &Y, const uint &Z) {

  for(uint k = 4; k < 5; k++) {
    for(uint i = 1 ; i < (X+1); i++ ) {
      gridA[(Z+1)/2*(Z+BW)*(Y+BW)+(Y+1)/2*(Y+BW)+i] = 8.8;
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
void advection(T * gridA, T * gridB, U * fieldA, U * fieldB,
    const uint &X, const uint &Y, const uint &Z) {

  for(uint k = BWP + omp_get_thread_num(); k < Z + BWP; k+=omp_get_num_threads()) {
    for(uint j = BWP; j < Y + BWP; j++) {
      for(uint i = BWP; i < X + BWP; i++) {
        uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+i;

        Triple origin;

        Triple backward;
        Triple backwardVel;

        Triple forward;
        Triple error;
        Triple corrected;

        origin[0] = i; origin[1] = j; origin[2] = k;

        for(int d = 0; d < 3; d++) {
          backward[d] = origin[d]*h-fieldA[cell][d]*dt;
        }

        bfceeInterpolationSteep(backward,backwardVel,fieldA,N,N,N);

        for(int d = 0; d < 3; d++) {
          forward[d] = backward[d] + backwardVel[d] * dt;
          error[d] = -1/2 * (origin[d]-forward[d]);
          corrected[d] = origin[d]-error[d]; 
          backward[d] = corrected[d] * h - fieldA[cell][d] * dt;
        }

        gridB[cell] = bfceeInterpolationSteep(backward,gridA,N,N,N);
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

template <typename T>
void WriteGidMesh(T * grid, 
    const uint &X, const uint &Y, const uint &Z,
    const char * fileName) {

  std::ofstream outputFile("grid.post.msh");

  outputFile << "MESH \"Grid\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
  outputFile << "# color 96 96 96" << std::endl;
  outputFile << "Coordinates" << std::endl;
  outputFile << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;

  for(uint k = 0; k < Z + BW; k++) {
    for(uint j = 0; j < Y + BW; j++) {
      uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+BWP;
      for(uint i = 0; i < X + BW; i++) {
        outputFile << cell++ << "  " << i << "  " << j << "  " << k << std::endl;
      }
    }
  }

  outputFile << "end coordinates" << std::endl;
  outputFile << "Elements" << std::endl;
  outputFile << "# Element node_1 node_2 node_3 node_4 node_5 node_6 node_7 node_8" << std::endl;

  for(uint k = BWP; k < Z + BWP-1; k++) {
    for(uint j = BWP; j < Y + BWP-1; j++) {
      uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+BWP;
      for(uint i = BWP; i < X + BWP; i++) {
        outputFile << cell++ << " ";

        outputFile << cell                        << " " << cell+1                    << "  ";
        outputFile << cell+1+(Y+BW)               << " " << cell+(Y+BW)               << "  ";
        outputFile << cell+(Z+BW)*(Y+BW)          << " " << cell+1+(Z+BW)*(Y+BW)      << "  ";
        outputFile << cell+1+(Z+BW)*(Y+BW)+(Y+BW) << " " << cell+(Z+BW)*(Y+BW)+(Y+BW) << "  ";

        outputFile << std::endl;
      }
    }
  }

  outputFile << "end Elements" << std::endl;
}

template <typename T>
void WriteGidResults(T * grid, 
    const uint &X, const uint &Y, const uint &Z,
    const char * fileName) {

  std::ofstream results(fileName);

  results << "GiD Post Results File 1.0" << std::endl << std::endl;
  results << "Result \"Temperature\" \"Kratos\" 1 Scalar OnNodes" << std::endl;
  results << "Values" << std::endl;

  for(uint k = 0; k < Z + BW; k++) {
    for(uint j = 0; j < Y + BW; j++) {
      uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+BWP;
      for(uint i = 0; i < X + BW; i++) {
        results << cell << "  " << grid[cell-1] << std::endl; cell++;
      }
    }
  }

  results << "End Values" << std::endl;
}

int main(int argc, char *argv[]) {

  N  = atoi(argv[1]);
  int steeps = atoi(argv[2]);
  h = 1;
  dx = h/N;

  double * step0 = NULL;
  double * step1 = NULL;

  Triple * velf0 = NULL;
  Triple * velf1 = NULL;

  // Temperature
  AllocateGrid(&step0,N,N,N);
  AllocateGrid(&step1,N,N,N);

  // Velocity
  AllocateGrid(&velf0,N,N,N);
  AllocateGrid(&velf1,N,N,N);

  printf("Allocation correct\n");

  printf("Initialize...\n");
  Initialize(step0,step1,N,N,N);
  WriteHeatFocus(step0,N,N,N);
  InitializeVelocity(velf0,N,N,N);

  dt = CFL*h/maxv;

  WriteGidMesh(step0,N,N,N,"grid.dat");
  printf("Begin...\n");
  for(int i = 0; i < steeps; i++) {
      advection(step0,step1,velf0,velf1,N,N,N);
      std::stringstream name;
      name << "grid.post.res." << i;
      std::cout << name.str() << std::endl;
      WriteGidResults(step0,N,N,N,name.str().c_str());
      //difussion(step1,step2,N,N,N);
      std::swap(step0,step1);
  }
  //WriteGridXYZ(velf0,N,N,N,"velc.dat");

  ReleaseGrid(&step0);
  ReleaseGrid(&step1);

  ReleaseGrid(&velf0);
  ReleaseGrid(&velf1);

  printf("De-Allocation correct\n");
}
