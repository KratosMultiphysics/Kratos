#include "fractional_step_explicit_solver.h"

void FractionalStepExplicitSolver::calculateAcceleration(
    double * initlVel,
    double * finalVel,
    double * acceleration,
    const std::size_t &Dim) {

  for(std::size_t kk = mBorderWidth[2]; kk < rZ + mBorderWidth[2]; kk += mSs) {
    #pragma omp task                                        \
      depend(in: initlVel[(kk)*mSlice3D:mSs*mSlice3D])      \
      depend(in: finalVel[(kk)*mSlice3D:mSs*mSlice3D])      \
      depend(out:acceleration[(kk)*mSlice3D:mSs*mSlice3D])  \
      shared(Dim)
    for(std::size_t k = kk; k < kk+mSs; k++) {
      for(std::size_t j = mBorderWidth[1]; j < rY + mBorderWidth[1]; j++) {
        for(std::size_t i = mBorderWidth[0]; i < rX + mBorderWidth[0]; i++) {
          std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP+i;
          for (std::size_t d = 0; d < Dim; d++) {
            acceleration[cell*Dim+d] = (finalVel[cell*Dim+d] - initlVel[cell*Dim+d]) * rIdt;
          }
        }
      }
    }
  }
}

void FractionalStepExplicitSolver::calculateLapplacian(
    double * input,
    double * output,
    const std::size_t &Dim) {

  double stencilSize = 0.0f;

  for(std::size_t kk = mBorderWidth[2]; kk < rZ + mBorderWidth[2]; kk+=mSs) {
    #pragma omp task \
      depend(in: input[(kk)*mSlice3D:mSs*mSlice3D]) \
      depend(out:output[(kk)*mSlice3D:mSs*mSlice3D]) \
      shared(Dim)
    for(std::size_t k = kk; k < kk+mSs; k++) {
      for(std::size_t j = mBorderWidth[1]; j < rY + mBorderWidth[1]; j++) {
        std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
        for(std::size_t i = mBorderWidth[0]; i < rX + mBorderWidth[0]; i++) {
          for (std::size_t d = 0; d < Dim; d++) {

            output[cell*Dim+d] = 0.0f;
            stencilSize        = 0.0f;

            if(i != rBWP)          {stencilSize += 1.0f; output[cell*Dim+d] += input[(cell - 1)*Dim+d];}
            if(i != rX + rBWP - 1) {stencilSize += 1.0f; output[cell*Dim+d] += input[(cell + 1)*Dim+d];}
            if(j != rBWP)          {stencilSize += 1.0f; output[cell*Dim+d] += input[(cell - (rX+BW))*Dim+d];}
            if(j != rY + rBWP - 1) {stencilSize += 1.0f; output[cell*Dim+d] += input[(cell + (rX+BW))*Dim+d];}
            if(k != rBWP)          {stencilSize += 1.0f; output[cell*Dim+d] += input[(cell - (rY+BW)*(rX+BW))*Dim+d];}
            if(k != rZ + rBWP - 1) {stencilSize += 1.0f; output[cell*Dim+d] += input[(cell - (rY+BW)*(rX+BW))*Dim+d];}

            output[cell*Dim+d] -= (stencilSize * input[(cell)*Dim+d]);
            output[cell*Dim+d] *= (rIdx * rIdx);
          }

          cell++;
        }
      }
    }
  }
}

void FractionalStepExplicitSolver::calculateGradient(
    double * press,
    double * gridB ) {

  double pressGrad[3];

  std::size_t offset[3] = {1,(rX+BW),(rY+BW)*(rX+BW)};
  std::size_t limit[3]  = {rX,rY,rZ};

  // Apply the pressure gradient
  for(std::size_t kk = mBorderWidth[2]; kk < rZ + mBorderWidth[2]; kk+=mSs) {
    #pragma omp task \
      depend(in: press[(kk)*mSlice:mSs*mSlice]) \
      depend(out:gridB[(kk)*mSlice3D:mSs*mSlice3D])
    for(std::size_t k = kk; k < kk+mSs; k++) {
      for(std::size_t j = mBorderWidth[1]; j < rY + mBorderWidth[1]; j++) {
        std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
        for(std::size_t i = mBorderWidth[0]; i < rX + mBorderWidth[0]; i++) {

          std::size_t cellIndex[3] = {i,j,k};

          for(std::size_t d = 0; d < rDim; d++) {
            if(cellIndex[d] == rBWP) {
              pressGrad[d] = (
                press[(cell + offset[d])] -
                press[(cell           )]) * 2.0f;
            } else if (cellIndex[d] == limit[d] + rBWP - 1) {
              pressGrad[d] = (
                press[(cell            )] -
                press[(cell - offset[d])]) * 2.0f;
            } else {
              pressGrad[d] = (
                press[(cell + offset[d])] -
                press[(cell - offset[d])]);
            }

            gridB[cell*rDim+d] = pressGrad[d] * 0.5f * rIdx;
          }

          cell++;
        }
      }
    }
  }
}

void FractionalStepExplicitSolver::calculateDivergence(
    double * gridA,
    double * gridB ) {

  std::size_t offset[3] = {1,(rX+BW),(rY+BW)*(rX+BW)};
  std::size_t limit[3]  = {rX,rY,rZ};

  // Combine it all together and store it back in A
  for(std::size_t kk = mBorderWidth[2]; kk < rZ + mBorderWidth[2]; kk+=mSs) {
    #pragma omp task \
      depend(in: gridA[(kk)*mSlice3D:mSs*mSlice3D]) \
      depend(out:gridB[(kk)*mSlice:mSs*mSlice])
    for(std::size_t k = kk; k < kk+mSs; k++) {
      for(std::size_t j = mBorderWidth[1]; j < rY + mBorderWidth[1]; j++) {
        std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
        for(std::size_t i = mBorderWidth[0]; i < rX + mBorderWidth[0]; i++) {
          gridB[cell] = 0.0f;

          std::size_t cellIndex[3] = {i,j,k};

          for(std::size_t d = 0; d < rDim; d++) {
            if(cellIndex[d] == rBWP) {
              gridB[cell] += (
                gridA[(cell + offset[d]) * rDim + d] -
                gridA[(cell            ) * rDim + d]) * 2.0f;
            } else if (cellIndex[d] == limit[d] + rBWP - 1) {
              gridB[cell] += (
                gridA[(cell            ) * rDim + d] -
                gridA[(cell - offset[d]) * rDim + d]) * 2.0f;
            } else {
              gridB[cell] += (
                gridA[(cell + offset[d]) * rDim + d] -
                gridA[(cell - offset[d]) * rDim + d]);
            }
          }

          gridB[cell] = gridB[cell] * 0.5f * rIdx;
          cell++;
        }
      }
    }
  }
}

void FractionalStepExplicitSolver::calculateVelocity(
    double * velLapp,
    double * pressGrad,
    double * force,
    double * acc,
    double * initVel) {

  for(std::size_t kk = mBorderWidth[2]; kk < rZ + mBorderWidth[2]; kk+=mSs) {
    #pragma omp task \
      depend(in:velLapp[(kk)*mSlice3D:mSs*mSlice3D]) \
      depend(in:pressGrad[(kk)*mSlice3D:mSs*mSlice3D]) \
      depend(in:acc[(kk)*mSlice3D:mSs*mSlice3D]) \
      depend(out:initVel[(kk)*mSlice3D:mSs*mSlice3D])
    for(std::size_t k = kk; k < kk+mSs; k++) {
      for(std::size_t j = mBorderWidth[1]; j < rY + mBorderWidth[1]; j++) {
        std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
        for(std::size_t i = mBorderWidth[0]; i < rX + mBorderWidth[0]; i++) {
          if(!(pFlags[cell] & FIXED_VELOCITY_X))
            initVel[cell*rDim+0] += ((rMu * velLapp[cell*rDim+0] / rRo) - (pressGrad[cell*rDim+0] / rRo) + (force[0] / rRo) - acc[cell*rDim+0]) * rDt;
          if(!(pFlags[cell] & FIXED_VELOCITY_Y))
            initVel[cell*rDim+1] += ((rMu * velLapp[cell*rDim+1] / rRo) - (pressGrad[cell*rDim+1] / rRo) + (force[1] / rRo) - acc[cell*rDim+1]) * rDt;
          if(!(pFlags[cell] & FIXED_VELOCITY_Z))
            initVel[cell*rDim+2] += ((rMu * velLapp[cell*rDim+2] / rRo) - (pressGrad[cell*rDim+2] / rRo) + (force[2] / rRo) - acc[cell*rDim+2]) * rDt;
          cell++;
        }
      }
    }
  }
}

void FractionalStepExplicitSolver::calculatePressureDiff(
    double * velDiv,
    double * initVel,
    double * pressDiff) {

  for(std::size_t kk = mBorderWidth[2]; kk < rZ + mBorderWidth[2]; kk+=mSs) {
    #pragma omp task \
      depend(in:initVel[(kk)*mSlice3D:mSs*mSlice3D]) \
      depend(out:pressDiff[(kk)*mSlice:mSs*mSlice])
    for(std::size_t k = kk; k < kk+mSs; k++) {
      for(std::size_t j = mBorderWidth[1]; j < rY + mBorderWidth[1]; j++) {
        std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
        for(std::size_t i = mBorderWidth[0]; i < rX + mBorderWidth[0]; i++) {
          pressDiff[cell] = -rRo*rCC2*rDt * velDiv[cell];
          cell++;
        }
      }
    }
  }
}

void FractionalStepExplicitSolver::calculatePressure(
    double * pressDiff,
    double * press) {

  // Combine it all together and store it back in A
  for(std::size_t kk = mBorderWidth[2]; kk < rZ + mBorderWidth[2]; kk+=mSs) {
    #pragma omp task \
      depend(in:pressDiff[(kk)*mSlice:mSs*mSlice]) \
      depend(out:press[(kk)*mSlice:mSs*mSlice])
    for(std::size_t k = kk; k < kk+mSs; k++) {
      for(std::size_t j = mBorderWidth[1]; j < rY + mBorderWidth[1]; j++) {
        std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
        for(std::size_t i = mBorderWidth[0]; i < rX + mBorderWidth[0]; i++) {
          if(!(pFlags[cell] & FIXED_PRESSURE))
            press[cell] += pressDiff[cell];// + lapp_fact * pressLapp[cell];
          cell++;
        }
      }
    }
  }
}

void FractionalStepExplicitSolver::calculateSmooth(
    double * gridA,
    double * gridB,
    const std::size_t &Dim) {

  double stencilSize = 0.0f;

  for(std::size_t kk = mBorderWidth[2]; kk < rZ + mBorderWidth[2]; kk+=mSs) {
    #pragma omp task \
      depend(in: gridA[(kk)*mSlice3D:mSs*mSlice3D]) \
      depend(out:gridB[(kk)*mSlice3D:mSs*mSlice3D]) \
      shared(Dim)
    for(std::size_t k = kk; k < kk+mSs; k++) {
      for(std::size_t j = mBorderWidth[1]; j < rY + mBorderWidth[1]; j++) {
        std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
        for(std::size_t i = mBorderWidth[0]; i < rX + mBorderWidth[0]; i++) {
          for (std::size_t d = 0; d < Dim; d++) {

            gridB[cell*Dim+d] = 0.0f;
            stencilSize       = 0.0f;

            if(i != rBWP)          {stencilSize += 1.0f; gridB[cell*Dim+d] += gridA[(cell - 1)*Dim+d];}
            if(i != rX + rBWP - 1) {stencilSize += 1.0f; gridB[cell*Dim+d] += gridA[(cell + 1)*Dim+d];}
            if(j != rBWP)          {stencilSize += 1.0f; gridB[cell*Dim+d] += gridA[(cell - (rX+BW))*Dim+d];}
            if(j != rY + rBWP - 1) {stencilSize += 1.0f; gridB[cell*Dim+d] += gridA[(cell + (rX+BW))*Dim+d];}
            if(k != rBWP)          {stencilSize += 1.0f; gridB[cell*Dim+d] += gridA[(cell - (rY+BW)*(rX+BW))*Dim+d];}
            if(k != rZ + rBWP - 1) {stencilSize += 1.0f; gridB[cell*Dim+d] += gridA[(cell - (rY+BW)*(rX+BW))*Dim+d];}

            gridB[cell*Dim+d] = (gridB[cell*Dim+d] + gridA[cell*Dim+d]) / (stencilSize+1.0f);
          }

          cell++;
        }
      }
    }
  }
}

FractionalStepExplicitSolver::FractionalStepExplicitSolver::(
      Block * block,
      const double& Dt,
      const double& Pdt) :
      Solver(block,Dt,Pdt) {
  }

FractionalStepExplicitSolver::~FractionalStepExplicitSolver() {
  // NYI
}

void FractionalStepExplicitSolver::Prepare_impl() {
  // NYI
}

void FractionalStepExplicitSolver::Finish_impl() {
  // NYI
}

void FractionalStepExplicitSolver::Execute_impl() {

  // Alias for the buffers
  double * initVel   = pBuffers[VELOCITY];
  double * vel       = pBuffers[AUX_3D_1];
  double * acc       = pBuffers[AUX_3D_3];
  double * press     = pBuffers[PRESSURE];

  double * pressGrad = pBuffers[AUX_3D_4];
  double * velLapp   = pBuffers[AUX_3D_2];

  double * velDiv    = pBuffers[AUX_3D_5];
  double * pressDiff = pBuffers[AUX_3D_6];
  double * pressLapp = pBuffers[AUX_3D_7];

  double force[3]    = {0.0f, 0.0f, -9.8f};
  // double force[3]    = {0.0f, 0.0f, 0.0f};

  std::size_t listL[rX*rX];
  std::size_t listR[rX*rX];
  std::size_t listF[rX*rX];
  std::size_t listB[rX*rX];
  std::size_t listT[rX*rX];
  std::size_t listD[rX*rX];

  int normalL[3] = {0,-1,0};
  int normalR[3] = {0,1,0};
  int normalF[3] = {-1,0,0};
  int normalB[3] = {1,0,0};
  int normalD[3] = {0,0,-1};
  int normalT[3] = {0,0,1};

  int normal2L[3] = {0,-2,0};
  int normal2R[3] = {0,2,0};
  int normal2F[3] = {-2,0,0};
  int normal2B[3] = {2,0,0};
  int normal2D[3] = {0,0,-2};
  int normal2T[3] = {0,0,2};

  uint counter = 0;

  for(uint a = rBWP; a < rZ + rBWP; a++) {
    for(uint b = rBWP; b < rY + rBWP; b++) {

      listL[counter] = a*(rZ+rBW)*(rY+rBW)+2*(rZ+rBW)+b;
      listR[counter] = a*(rZ+rBW)*(rY+rBW)+(rY-1)*(rZ+rBW)+b;
      listF[counter] = a*(rZ+rBW)*(rY+rBW)+b*(rZ+rBW)+2;
      listB[counter] = a*(rZ+rBW)*(rY+rBW)+b*(rZ+rBW)+(rX-1);
      listD[counter] = 2*(rZ+rBW)*(rY+rBW)+a*(rZ+rBW)+b;
      listT[counter] = (rZ-1)*(rZ+rBW)*(rY+rBW)+a*(rZ+rBW)+b;

      counter++;
    }
  }

  ///////////////////////////////////////////////////////////////////////////

  // // Calculate acceleration
  // #pragma omp parallel for
  // for(std::size_t k = rBWP; k < rZ + rBWP; k++) {
  //   for(std::size_t j = rBWP; j < rY + rBWP; j++) {
  //     std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
  //     for(std::size_t i = rBWP; i < rX + rBWP; i++) {
  //       calculateAcceleration(initVel,vel,acc,cell++,3);
  //     }
  //   }
  // }
  //
  // // Apply the pressure gradient
  // #pragma omp parallel for
  // for(std::size_t k = rBWP; k < rZ + rBWP; k++) {
  //   for(std::size_t j = rBWP; j < rY + rBWP; j++) {
  //     std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
  //     for(std::size_t i = rBWP; i < rX + rBWP; i++) {
  //       gradient(press,pressGrad,cell++);
  //     }
  //   }
  // }
  //
  // // divergence of the gradient of the velocity
  // #pragma omp parallel for
  // for(std::size_t k = rBWP; k < rZ + rBWP; k++) {
  //   for(std::size_t j = rBWP; j < rY + rBWP; j++) {
  //     std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
  //     for(std::size_t i = rBWP; i < rX + rBWP; i++) {
  //       lapplacian(vel,velLapp,cell++,3);
  //     }
  //   }
  // }
  //
  // // Combine it all together and store it back in A
  // #pragma omp parallel for
  // for(std::size_t k = rBWP; k < rZ + rBWP; k++) {
  //   for(std::size_t j = rBWP; j < rY + rBWP; j++) {
  //     std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
  //     for(std::size_t i = rBWP; i < rX + rBWP; i++) {
  //       if(!(pFlags[cell] & FIXED_VELOCITY_X))
  //         initVel[cell*rDim+0] += (rMu * velLapp[cell*rDim+0] - pressGrad[cell*rDim+0] + force[0] / rRo - acc[cell*rDim+0]) * rDt;
  //       if(!(pFlags[cell] & FIXED_VELOCITY_Y))
  //         initVel[cell*rDim+1] += (rMu * velLapp[cell*rDim+1] - pressGrad[cell*rDim+1] + force[1] / rRo - acc[cell*rDim+1]) * rDt;
  //       if(!(pFlags[cell] & FIXED_VELOCITY_Z))
  //         initVel[cell*rDim+2] += (rMu * velLapp[cell*rDim+2] - pressGrad[cell*rDim+2] + force[2] / rRo - acc[cell*rDim+2]) * rDt;
  //       cell++;
  //     }
  //   }
  // }
  //
  // // Combine it all together and store it back in A
  // #pragma omp parallel for
  // for(std::size_t k = rBWP; k < rZ + rBWP; k++) {
  //   for(std::size_t j = rBWP; j < rY + rBWP; j++) {
  //     std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
  //     for(std::size_t i = rBWP; i < rX + rBWP; i++) {
  //       divergence(initVel,velDiv,cell);
  //       pressDiff[cell] = -rRo*rCC2*rDt * velDiv[cell];
  //       cell++;
  //     }
  //   }
  // }
  //
  // // Combine it all together and store it back in A
  // #pragma omp parallel for
  // for(std::size_t k = rBWP; k < rZ + rBWP; k++) {
  //   for(std::size_t j = rBWP; j < rY + rBWP; j++) {
  //     std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
  //     for(std::size_t i = rBWP; i < rX + rBWP; i++) {
  //       lapplacian(pressDiff,pressLapp,cell,1);
  //       pressLapp[cell] *= (rDx * rDx);
  //       cell++;
  //     }
  //   }
  // }
  //
  //
  // // Combine it all together and store it back in A
  // #pragma omp parallel for
  // for(std::size_t k = rBWP; k < rZ + rBWP; k++) {
  //   for(std::size_t j = rBWP; j < rY + rBWP; j++) {
  //     std::size_t cell = k*(rZ+rBW)*(rY+rBW)+j*(rY+BW)+rBWP;
  //     for(std::size_t i = rBWP; i < rX + rBWP; i++) {
  //       if(!(pFlags[cell] & FIXED_PRESSURE))
  //         press[cell] += pressDiff[cell] + pressLapp[cell] * ( rDt / rRo );
  //       cell++;
  //     }
  //   }
  // }
  //
  // applyBc(press,listL,rX*rX,normalL,1,1);
  // applyBc(press,listR,rX*rX,normalR,1,1);
  // applyBc(press,listF,rX*rX,normalF,1,1);
  // applyBc(press,listB,rX*rX,normalB,1,1);
  // applyBc(press,listD,rX*rX,normalD,1,1);
  //
  // applyBc(press,listL,rX*rX,normal2L,1,1);
  // applyBc(press,listR,rX*rX,normal2R,1,1);
  // applyBc(press,listF,rX*rX,normal2F,1,1);
  // applyBc(press,listB,rX*rX,normal2B,1,1);
  // applyBc(press,listD,rX*rX,normal2D,1,1);
  //
  //
  // applyBc(press,listT,rX*rX,normalT,1,1);
  // // applyBc(press,listT,rX*rX,normal2T,1,1);

}

void FractionalStepExplicitSolver::ExecuteBlock_impl() {}

void FractionalStepExplicitSolver::ExecuteTask_impl() {

  // Alias for the buffers
  double * initVel   = pBuffers[VELOCITY];
  double * vel       = pBuffers[AUX_3D_1];
  double * acc       = pBuffers[AUX_3D_3];
  double * press     = pBuffers[PRESSURE];

  double * pressGrad = pBuffers[AUX_3D_4];
  double * velLapp   = pBuffers[AUX_3D_2];

  double * velDiv    = pBuffers[AUX_3D_5];
  double * pressDiff = pBuffers[AUX_3D_6];
  double * pressLapp = pBuffers[AUX_3D_7];

  double force[3]    = {0.0f, 0.0f, -9.8f};
  // double force[3]    = {0.0f, 0.0f, 0.0f};

  std::size_t listL[rX*rX];
  std::size_t listR[rX*rX];
  std::size_t listF[rX*rX];
  std::size_t listB[rX*rX];
  std::size_t listT[rX*rX];
  std::size_t listD[rX*rX];

  int normalL[3] = {0,-1,0};
  int normalR[3] = {0,1,0};
  int normalF[3] = {-1,0,0};
  int normalB[3] = {1,0,0};
  int normalT[3] = {0,0,-1};
  int normalD[3] = {0,0,1};

  uint counter = 0;

  for(uint a = rBWP; a < rZ + rBWP; a++) {
    for(uint b = rBWP; b < rY + rBWP; b++) {

      listL[counter] = a*(rZ+rBW)*(rY+rBW)+1*(rZ+rBW)+b;
      listR[counter] = a*(rZ+rBW)*(rY+rBW)+(rY)*(rZ+rBW)+b;
      listF[counter] = a*(rZ+rBW)*(rY+rBW)+b*(rZ+rBW)+1;
      listB[counter] = a*(rZ+rBW)*(rY+rBW)+b*(rZ+rBW)+(rX);
      listT[counter] = 1*(rZ+rBW)*(rY+rBW)+a*(rZ+rBW)+b;
      listD[counter] = (rZ)*(rZ+rBW)*(rY+rBW)+a*(rZ+rBW)+b;

      counter++;
    }
  }

  ////////////////////////////////////////////////////////////////////////////

  mBt       = (rX + rBW);
  mSs       = 1;
  mSlice    = mBt * mBt;
  mSlice3D  = mSlice * 3;

  #pragma omp parallel
  #pragma omp single
  {
    #pragma omp taskwait
    calculateAcceleration(vel,initVel,acc,3);

    #pragma omp taskwait
    calculateGradient(press,pressGrad);

    #pragma omp taskwait
    calculateLapplacian(initVel,velLapp,3);

    #pragma omp taskwait
    calculateVelocity(velLapp,pressGrad,force,vel,initVel);

    #pragma omp taskwait
    calculateDivergence(initVel,velDiv);

    #pragma omp taskwait
    calculatePressureDiff(velDiv,initVel,pressDiff);

    #pragma omp taskwait
    calculatePressure(pressDiff,press);

    // #pragma omp taskwait
    // calculateSmooth(pressLapp,press);
  }
}

void FractionalStepExplicitSolver::ExecuteVector_impl() {

  // TODO: Implement this with the new arrays

  // uint cell, pcellb, pcelle;
  //
  // double __attribute__((aligned(ALIGN))) tmpa[VP];
  // double __attribute__((aligned(ALIGN))) tmpb[VP];
  //
  // for(uint kk = 0; kk < rNB; kk++) {
  //   for(uint jj = 0; jj < rNB; jj++) {
  //     for(uint k = rBWP + (kk * rNE); k < rBWP + ((kk+1) * rNE); k++) {
  //       for(uint j = rBWP + (jj * rNE); j < rBWP + ((jj+1) * rNE); j++) {
  //
  //         pcellb = k*(rX+rBW)*(rX+rBW)/VP+j*(rX+rBW)/VP+(rBWP/VP);
  //         pcelle = k*(rX+rBW)*(rX+rBW)/VP+j*(rX+rBW)/VP+((rX/VP) + rBWP/VP - 1);
  //
  //         VectorType * left  = &pPhiA[pcellb-1];
  //         VectorType * right = &pPhiA[pcellb+1];
  //         VectorType * bott  = &pPhiA[pcellb-(rX+rBW)/VP];
  //         VectorType * top   = &pPhiA[pcellb+(rX+rBW)/VP];
  //         VectorType * front = &pPhiA[pcellb-(rX+rBW)*(rX+rBW)/VP];
  //         VectorType * back  = &pPhiA[pcellb+(rX+rBW)*(rX+rBW)/VP];
  //
  //         // Prefix
  //         cell = pcellb;
  //         VSTORE(tmpa,pPhiA[pcelle]);
  //         tmpb[0] = 0.0;
  //         for(uint p = 1; p < VP; p++) {
  //           tmpb[p] = tmpa[p-1];
  //         }
  //
  //         left++;
  //         pPhiB[cell++] = VSTENCIL(VLOAD(tmpb),*right++,*top++,*bott++,*front++,*back++);
  //
  //         // Body
  //         for(uint i = rBWP/VP + 1; i < (rX/VP) + BWP/VP - 1; i++) {
  //           pPhiB[cell++] = VSTENCIL(*left++,*right++,*top++,*bott++,*front++,*back++);
  //         }
  //
  //         // Sufix
  //         cell = pcelle;
  //         VSTORE(tmpa,pPhiA[pcellb]);
  //         for(uint p = 1; p < VP; p++) {
  //           tmpb[p-1] = tmpa[p];
  //         }
  //         tmpb[VP-1] = 0.0;
  //
  //         pPhiB[cell] = VSTENCIL(*left++,VLOAD(tmpb),*top++,*bott++,*front++,*back++);
  //       }
  //     }
  //   }
  // }

}

void FractionalStepExplicitSolver::SetDiffTerm(double diffTerm) {
  mDiffTerm = diffTerm;
}
