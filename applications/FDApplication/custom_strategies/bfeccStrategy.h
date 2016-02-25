#if !defined(KRATOS_BFECC_SOLVER_STRATEGY)
#define KRATOS_BFECC_SOLVER_STRATEGY

// External includes

// System includes
#include <unordered_map>
#include <memory>
#include <limits>
#include <cmath>

// Kratos includes
#include "includes/define.h"
#include "spaces/ublas_space.h"

// Project includes
#include "fd_application.h"
#include "custom_utilities/bfecc.h"
#include "custom_utilities/bfecc_utils.h"
#include "custom_utilities/grid_printer.h"

namespace Kratos {

template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class BfeccSolverStrategy: public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> {

public:
  /// Pointer definition of BfeccSolverStrategy
  KRATOS_CLASS_POINTER_DEFINITION(BfeccSolverStrategy);

  BfeccSolverStrategy(ModelPart& model_part) :
      SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, true),
      mCellSize(std::vector<double>(3)),     // Size of the cells
      mNumCells(std::vector<std::size_t>(3)),     // Number of cells
      mBorderWidth(std::vector<std::size_t>(3)),  // Size of the border
      mConditionsModelPart(model_part) {
  }

  ~BfeccSolverStrategy(){
    std::cout << "Deleting grids..." << std::endl;
    for(auto grid: mGrids) {
      std::cout << "\tDeleting " << grid.first << " grid..." << std::endl;
      free(grid.second);
    }

    std::cout << "Deleting falgs..." << std::endl;
    free(mFlags);

    std::cout << "Deleting process finished" << std::endl;
    free(mPrinter);
  }

  virtual double Solve() {
    return 0.0f;
  }

  virtual void Clear() {
  }

  virtual bool IsConverged() {
    return true;
  }

  virtual void CalculateOutputData() {
  }

  virtual void Predict() {
  }

  virtual void Initialize() {

    std::cout << "Defining grids..." << std::endl;

    std::vector<std::string> gridNames = {
      "VARIABLE",
      "VARIABLE_BUFFER_1",
      "VELOCITY",
      "VELOCITY_BUFFER_1",
      "ACCELERATION",
      "PRESSURE",
      "PRESSURE_BUFFER_1"
    };


    // mDt = 1.0f/32.0f;

    for(std::size_t i = 0; i < 3; i++) {
      // mNumCells[i] = 32;
      mCellSize[i] = 1.0f/mNumCells[i];
      std::cout << "CellSize for " << i << ": " << mCellSize[i] << std::endl;
      // mBorderWidth[i] = 1;
    }

    std::cout << "Creating grids..." << std::endl;

    std::size_t numberOfCells =
      (mNumCells[0] + mBorderWidth[0] * 2) *
      (mNumCells[1] + mBorderWidth[1] * 2) *
      (mNumCells[2] + mBorderWidth[2] * 2);

    for(auto name: gridNames) {
      mGrids.emplace(std::pair<std::string,double*>(
        name,
        std::move((double*)malloc(sizeof(double) * numberOfCells * 3))
      ));
    }

    for(auto grid: mGrids) {
      std::cout << "\tCreated " << grid.first << " grid" << std::endl;
      std::cout << "\t  Number of cells: " << numberOfCells << std::endl;
      std::cout << "\t  X: " << mNumCells[0] << " cells of size " << mCellSize[0] << std::endl;
      std::cout << "\t  Y: " << mNumCells[1] << " cells of size " << mCellSize[1] << std::endl;
      std::cout << "\t  Z: " << mNumCells[2] << " cells of size " << mCellSize[2] << std::endl;
    }

    std::cout << "Creating flags..." << std::endl;
    mFlags = (int *)malloc(sizeof(int) * numberOfCells);

    for(std::size_t i = 0 ; i < numberOfCells; i++) {
      mFlags[i] = 0;
    }

    std::cout << "Applying BC from the model..." << std::endl;
    double bbMin[3] = {DBL_MAX};
    double bbMax[3] = {DBL_MIN};
    double bbSize[3];

    // Parse the ConditionModelpart setting the size and conditions
    for(auto elem: mConditionsModelPart.Conditions()) {
      for(auto point: elem.GetGeometry().Points()) {
        for(std::size_t i = 0; i < 3 ; i++) {
          bbMin[i] = point[i] < bbMin[i] ? point[i] : bbMin[i];
          bbMax[i] = point[i] > bbMax[i] ? point[i] : bbMax[i];
        }
      }
    }

    for(std::size_t i = 0; i < 3; i++) {
      bbSize[i]  = bbMax[i] - bbMin[i];
    }

    std::cout << "BB: " << std::endl;
    std::cout << "\t" << bbMin[0] << " " << bbMin[1] << " " << bbMin[2] << std::endl;
    std::cout << "\t" << bbMax[0] << " " << bbMax[1] << " " << bbMax[2] << std::endl;

    // Map the BC from the model to our cell
    for(auto elem: mConditionsModelPart.Conditions()) {

      double lcMin[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
      double lcMax[3] = {DBL_MIN, DBL_MIN, DBL_MIN};

      bool fixedPressure = false;
      std::vector<bool> fixedVelocity = {false, false, false};

      // std::cout << "Condition: " << elem.Id() << std::endl;

      for(auto point: elem.GetGeometry().Points()) {
        for(std::size_t i = 0; i < 3 ; i++) {
          fixedPressure    = fixedPressure    | point.IsFixed(PRESSURE);
          fixedVelocity[0] = fixedVelocity[0] | point.IsFixed(VELOCITY_X);
          fixedVelocity[1] = fixedVelocity[1] | point.IsFixed(VELOCITY_Y);
          fixedVelocity[2] = fixedVelocity[2] | point.IsFixed(VELOCITY_Z);

          lcMin[i] = point[i] < lcMin[i] ? point[i] : lcMin[i];
          lcMax[i] = point[i] > lcMax[i] ? point[i] : lcMax[i];
        }
      }

      if(fixedPressure    || fixedVelocity[0] ||
         fixedVelocity[1] || fixedVelocity[2]) {

        double condBegin[3] = {
          (lcMin[0] - bbMin[0]) / bbSize[0],
          (lcMin[1] - bbMin[1]) / bbSize[1],
          (lcMin[2] - bbMin[2]) / bbSize[2]
        };

        double condEnd[3] = {
          (lcMax[0] - bbMin[0]) / bbSize[0],
          (lcMax[1] - bbMin[1]) / bbSize[1],
          (lcMax[2] - bbMin[2]) / bbSize[2]
        };

        std::size_t cellBegin[3] = {
          (std::size_t)std::floor(condBegin[0] / mCellSize[0]) + mBorderWidth[0],
          (std::size_t)std::floor(condBegin[1] / mCellSize[1]) + mBorderWidth[1],
          (std::size_t)std::floor(condBegin[2] / mCellSize[2]) + mBorderWidth[2]
        };

        std::size_t cellEnd[3] = {
          (std::size_t)std::ceil(condEnd[0] / mCellSize[0]) + mBorderWidth[0],
          (std::size_t)std::ceil(condEnd[1] / mCellSize[1]) + mBorderWidth[1],
          (std::size_t)std::ceil(condEnd[2] / mCellSize[2]) + mBorderWidth[2]
        };

        // Correct parallel planes
        for(std::size_t d = 0; d < 3; d++) {
          if(cellBegin[d] == cellEnd[d]) {
            if(cellBegin[d] > mBorderWidth[d]) {
              cellBegin[d]--;
            } else {
              cellEnd[d]++;
            }
          }
        }

        std::cout << "Condition marking for:" << std::endl;
        std::cout << "\t" << cellBegin[0] << " " << cellBegin[1] << " " << cellBegin[2] << std::endl;
        std::cout << "\t" << cellEnd[0] << " " << cellEnd[1] << " " << cellEnd[2] << std::endl;
        std::cout << "\t" << fixedPressure << " " << fixedVelocity[0] << " " << fixedVelocity[1] << " " << fixedVelocity[2] << std::endl;

        for(std::size_t i = cellBegin[0]; i < cellEnd[0]; i++) {
          for(std::size_t j = cellBegin[1]; j < cellEnd[1]; j++) {
            for(std::size_t k = cellBegin[2]; k < cellEnd[2]; k++) {
              std::size_t cellIndex = BfeccUtils::Index(i, j, k, mNumCells, mBorderWidth);
              // std::cout << cellIndex << " " << fixedPressure << " " << fixedVelocity[0] << " " << fixedVelocity[1] << " " << fixedVelocity[2] << std::endl;
              mFlags[cellIndex] |= FIXED_PRESSURE   * fixedPressure;
              mFlags[cellIndex] |= FIXED_VELOCITY_X * fixedVelocity[0];
              mFlags[cellIndex] |= FIXED_VELOCITY_Y * fixedVelocity[1];
              mFlags[cellIndex] |= FIXED_VELOCITY_Z * fixedVelocity[2];
            }
          }
        }
      }
    }

    std::cout << "Creating process finished" << std::endl;

    mPrinter = new GridPrinter(mCellSize[0], mNumCells, mBorderWidth);
    mPrinter->Initialize("Grid",32);
    mPrinter->WriteGidMeshBinary();

    /// Init ///

    #pragma omp parallel for
    for(std::size_t k = 0; k < mNumCells[2] + mBorderWidth[2] * 2; k++) {
      for(std::size_t j = 0; j < mNumCells[1] + mBorderWidth[1] * 2; j++) {
        for(std::size_t i = 0; i < mNumCells[0] + mBorderWidth[0] * 2; i++ ) {
          mGrids["VELOCITY"][BfeccUtils::Index(i,j,k,mNumCells,mBorderWidth)*3+0] = -1.0f * (double)(j-(mNumCells[1]+1.0f)/2.0f) * mCellSize[1];
          mGrids["VELOCITY"][BfeccUtils::Index(i,j,k,mNumCells,mBorderWidth)*3+1] =  1.0f * (double)(i-(mNumCells[0]+1.0f)/2.0f) * mCellSize[0];
          mGrids["VELOCITY"][BfeccUtils::Index(i,j,k,mNumCells,mBorderWidth)*3+2] =  0.0f;
        }
      }
    }

    WriteHeatFocus(mGrids["VARIABLE"]);

    mPrinter->WriteGidResultsBinary3D(mGrids["VELOCITY"],0,"VELOCITY");
    mPrinter->WriteGidResultsBinary3D(mGrids["VARIABLE"],0,"VARIABLE");
    mPrinter->WriteGidResultsBinary3D(mGrids["VARIABLE_BUFFER_1"],0,"BFECC_CONVECTED_VAR");
    mPrinter->WriteGidResultsBinary3D(mGrids["PRESSURE"],0,"PRESSURE");
    mPrinter->WriteGidResultsBinary1D(mFlags,0,"FLAGS");
  }

  void WriteHeatFocus(double * buffer) {

    std::size_t Xc, Yc, Zc;

    Xc = (std::size_t)(2.0f / 7.0f * (double)(mNumCells[0]));
  	Yc = (std::size_t)(2.0f / 7.5f * (double)(mNumCells[1]));
  	Zc = (std::size_t)(1.0f / 2.0f * (double)(mNumCells[2]));

    #pragma omp parallel for
    for(std::size_t k = 0; k < mNumCells[2] + mBorderWidth[2] * 2; k++) {
      for(std::size_t j = 0; j < mNumCells[1] + mBorderWidth[1] * 2; j++) {
        for(std::size_t i = 0; i < mNumCells[0] + mBorderWidth[0] * 2; i++ ) {

          double d2 =
            pow(((double)Xc - (double)(i)),2.0f) +
            pow(((double)Yc - (double)(j)),2.0f) +
            pow(((double)Zc - (double)(k)),2.0f);

          double rr =
            pow((double)mNumCells[0]/8.0,2.0f);

          for(std::size_t d = 0; d < 3; d++) {
            buffer[BfeccUtils::Index(i,j,k,mNumCells,mBorderWidth)*3+d] = 0.0f;
          }

          if(d2 < rr) {
            for(std::size_t d = 0; d < 3; d++) {
              buffer[BfeccUtils::Index(i,j,k,mNumCells,mBorderWidth)*3+d] = 1.0f-d2/rr;
            }
          }
        }
      }
    }
  }

  virtual void InitializeSolutionStep() {
  }

  virtual void FinalizeSolutionStep() {
  }

  virtual bool SolveSolutionStep() {

    auto & variable = mGrids["VARIABLE"];
    auto & convVar  = mGrids["VARIABLE_BUFFER_1"];
    auto & initVel  = mGrids["VELOCITY"];
    auto & auxgrid  = mGrids["VELOCITY_BUFFER_1"];

    BfeccConvecter bfecc(mFlags, mNumCells, mBorderWidth, mDt, mCellSize[0], 3);
    // FractionalStepExplicitSolver fractional(mFlags, mNumCells, mBorderWidth, mDt, mGrids, mCellSize[0]);

    bfecc.Convect(variable, auxgrid, initVel, convVar);
    // fractional.Solve();

    // Swap
    auto tmp = variable;
    variable = convVar;
    convVar  = tmp;

    return true;
  }

  virtual void WriteResults(double timeStep) {

    std::vector<std::string> resultList = {
      "VELOCITY",
      "VARIABLE",
    };

    for(auto res: resultList) {
      mPrinter->WriteGidResultsBinary3D(mGrids[res],timeStep,res);
    }
    mPrinter->WriteGidResultsBinary1D(mFlags,timeStep,"FLAGS");
  }

  double GetDt() {
    return mDt;
  }

  void SetDt(double dt) {
    mDt = dt;
  }

  std::vector<std::size_t> GetNumCells() {
    return mNumCells;
  }

  void SetNumCells(std::vector<std::size_t> numCells) {
    mNumCells = numCells;
  }

  std::vector<std::size_t> GetBorderWidth() {
    return mBorderWidth;
  }

  void SetBorderWidth(std::vector<std::size_t> borderWidth) {
    mBorderWidth = borderWidth;
  }

private:

  double      mDt;
  int  *      mFlags;
  std::vector<double> mCellSize;     // Size of the cells
  std::vector<std::size_t> mNumCells;     // Number of cells
  std::vector<std::size_t> mBorderWidth;  // Size of the border

  std::unordered_map<std::string,double*> mGrids; // Buffers

  GridPrinter * mPrinter;

  ModelPart & mConditionsModelPart;

};

}  // namespace Kratos.

#endif // KRATOS_BFECC_SOLVER_STRATEGY defined
