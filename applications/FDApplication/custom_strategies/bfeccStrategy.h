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
      "VELOCITY",
      "PRESSURE",
      "AUXGRID0",
      "TESTVARI",
      "SWAPAUXM"
    };

    mDt = 1.0f/8.0f;

    for(std::size_t i = 0; i < 3; i++) {
      mNumCells[i] = 8;
      mCellSize[i] = 1.0f/mNumCells[i];
      mBorderWidth[i] = 1;
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

      double lcMin[3] = {DBL_MAX};
      double lcMax[3] = {DBL_MIN};

      bool fixedPress = true;

      for(auto point: elem.GetGeometry().Points()) {
        for(std::size_t i = 0; i < 3 ; i++) {
          fixedPress &= point.IsFixed(PRESSURE);
          lcMin[i] = point[i] < lcMin[i] ? point[i] : lcMin[i];
          lcMax[i] = point[i] > lcMax[i] ? point[i] : lcMax[i];
        }
      }

      if(fixedPress) {
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

        for(std::size_t i = cellBegin[0]; i < cellEnd[0]; i++) {
          for(std::size_t j = cellBegin[1]; j < cellEnd[1]; j++) {
            for(std::size_t k = cellBegin[2]; k < cellEnd[2]; k++) {
              std::size_t cellIndex = BfeccUtils::Index(i, j, k, mNumCells, mBorderWidth);
              mFlags[cellIndex] |= FIXED_PRESSURE;
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

    WriteHeatFocus(mGrids["TESTVARI"]);

    mPrinter->WriteGidResultsBinary3D(mGrids["VELOCITY"],0,"VELOCITY");
    mPrinter->WriteGidResultsBinary3D(mGrids["AUXGRID0"],0,"AUXGRID0");
    mPrinter->WriteGidResultsBinary3D(mGrids["SWAPAUXM"],0,"SWAPAUXM");
    mPrinter->WriteGidResultsBinary3D(mGrids["TESTVARI"],0,"TESTVARI");
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

    // auto & variable = mGrids["VELOCITY"];
    auto & velocity = mGrids["VELOCITY"];
    // auto & pressure = mGrids["PRESSURE"];
    auto & auxgrida = mGrids["AUXGRID0"];
    auto & testvari = mGrids["TESTVARI"];
    auto & swapauxm = mGrids["SWAPAUXM"];

    BfeccConvecter bfecc(mFlags, mNumCells, mBorderWidth, mDt, mCellSize[0], 3);

    bfecc.Convect(testvari, auxgrida, velocity, swapauxm);

    static int a = 1;

    auto tmp = testvari;
    testvari = swapauxm;
    swapauxm = tmp;

    if(!(a%100)) {
      mPrinter->WriteGidResultsBinary3D(velocity,a,"VELOCITY");
      mPrinter->WriteGidResultsBinary3D(auxgrida,a,"AUXGRID0");
      mPrinter->WriteGidResultsBinary3D(swapauxm,a,"SWAPAUXM");
      mPrinter->WriteGidResultsBinary3D(testvari,a,"TESTVARI");
      mPrinter->WriteGidResultsBinary1D(mFlags,a,"FLAGS");
    }

    a++;

    return true;
  }

private:

  double      mDt;
  int  *      mFlags;
  double      mCellSize[3];     // Size of the cells
  std::size_t mNumCells[3];   // Number of cells
  std::size_t mBorderWidth[3];  // Size of the border

  std::unordered_map<std::string,double*> mGrids; // Buffers

  GridPrinter * mPrinter;

  ModelPart & mConditionsModelPart;

};

}  // namespace Kratos.

#endif // KRATOS_BFECC_SOLVER_STRATEGY defined
