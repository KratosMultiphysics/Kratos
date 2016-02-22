#include "solver.hpp"
#include "simd.hpp"

class FractionalStepExplicitSolver : public Solver<StencilSolver> {
private:

  void calculateAcceleration(
      double * gridA,
      double * gridB,
      double * gridC,
      const size_t &Dim);

  void calculateLapplacian(
      double * gridA,
      double * gridB,
      const size_t &Dim);

  void calculateGradient(
      double * press,
      double * gridB );

  void calculateDivergence(
      double * gridA,
      double * gridB );

  void calculateVelocity(
      double * velLapp,
      double * pressGrad,
      double * force,
      double * acc,
      double * initVel);

  void calculatePressureDiff(
      double * velDiv,
      double * initVel,
      double * pressDiff);

  void calculatePressure(
      double * pressDiff,
      double * press);

  void calculateSmooth(
      double * gridA,
      double * gridB,
      const size_t &Dim);

public:

  FractionalStepExplicitSolver(
      Block * block,
      const double& Dt,
      const double& Pdt);

  ~FractionalStepExplicitSolver();

  void Prepare_impl();

  void Finish_impl();

  void Execute_impl();

  void ExecuteBlock_impl();

  void ExecuteTask_impl();

  void ExecuteVector_impl();

  void SetDiffTerm(double diffTerm);

private:

  double mDiffTerm;

  int mBt;
  int mSs;
  int mSlice;
  int mSlice3D;

  ////////////////////////////////////////////////////////////////////////////

  double mDt;
  double mDx;
  double mIdx;

  std::size_t mDim;

  std::size_t * mNumCells;
  std::size_t * mBorderWidth;
  int * mFlags;

};
