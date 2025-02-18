//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
// 
#include "rve_wall_boundary_thermal_2d.h"

namespace Kratos
{
    //------------------------------------------------------------------------------------------------------------
    void RVEWallBoundaryThermal2D::PreProcessGlobalResults(void) {
        RVEWallBoundary2D::PreProcessGlobalResults();
        mConductivityTensor      = ZeroMatrix(mDim,mDim);
        mConductivityTensorInner = ZeroMatrix(mDim,mDim);
    }

    //------------------------------------------------------------------------------------------------------------
    void RVEWallBoundaryThermal2D::ProcessGlobalResults(void) {
        
    }

    //------------------------------------------------------------------------------------------------------------
    bool RVEWallBoundaryThermal2D::Homogenize(void) {
        if (!RVEWallBoundary2D::Homogenize()) return false;
        HomogenizeConductivity();
        return true;
    }

    //------------------------------------------------------------------------------------------------------------
    void RVEWallBoundaryThermal2D::HomogenizeConductivity(void) {
        for (unsigned int i = 0; i < mDim; i++) {
            for (unsigned int j = 0; j < mDim; j++) {
                mConductivityTensor(i,j)      /= mVolTotal;
                mConductivityTensorInner(i,j) /= mVolInner;
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------
    void RVEWallBoundaryThermal2D::WriteFileHeadersTensorResults(void) {
        mFileTensorResults << "R = ROW, C = COLUMN, F = FABRIC TENSOR, S = STRESS TENSOR, K = CONDUCTIVITY TENSOR" << std::endl;
        mFileTensorResults << "Ri -> C1: STEP | C2: TIME" << std::endl;
        mFileTensorResults << "Rj -> C1-C4: [F11 F12 F21 F22] | C5: FABRIC_INDEX | C6: ANISOTROPY (ALL CONTACTS)" << std::endl;
        mFileTensorResults << "Rk -> C1-C4: [F11 F12 F21 F22] | C5: FABRIC_INDEX | C6: ANISOTROPY (INN CONTACTS)" << std::endl;
        mFileTensorResults << "Rl -> C1-C4: [S11 S12 S21 S22] | C5: VOL STRESS | C6: DEV STRESS  | C7: WALL STRESS (ALL CONTACTS)" << std::endl;
        mFileTensorResults << "Rm -> C1-C4: [S11 S12 S21 S22] | C5: VOL STRESS | C6: DEV STRESS (INN CONTACTS)" << std::endl;
        mFileTensorResults << "Rn -> C1-C4: [K11 K12 K21 K22] (ALL CONTACTS)" << std::endl;
        mFileTensorResults << "Ro -> C1-C4: [K11 K12 K21 K22] (INN CONTACTS)" << std::endl;
    }

    //------------------------------------------------------------------------------------------------------------
    void RVEWallBoundaryThermal2D::WriteResultFilesTensorResults(void) {
        ProcessInfo& r_process_info = mDemModelPart->GetProcessInfo();

        mFileTensorResults << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME_STEPS] << " "
                           << std::setw(WIDTH_DEFAULT) << std::defaultfloat << r_process_info[TIME]       << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensor(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensor(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensor(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensor(1,1) << " "
                           << "] "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mFidx              << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mAnisotropy        << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensorInner(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensorInner(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensorInner(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mFabricTensorInner(1,1) << " "
                           << "] "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mFidxInner              << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mAnisotropyInner        << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensor(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensor(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensor(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensor(1,1) << " "
                           << "] "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mEffStress         << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mDevStress         << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mWallStress        << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensorInner(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensorInner(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensorInner(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mStressTensorInner(1,1) << " "
                           << "] "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mEffStressInner         << " "
                           << std::setw(WIDTH_FLOAT06) << std::fixed << std::setprecision(6)  << mDevStressInner         << " "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensor(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensor(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensor(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensor(1,1) << " "
                           << "] "
                           << std::endl;

        mFileTensorResults << "[ " 
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensorInner(0,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensorInner(0,1) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensorInner(1,0) << " "
                           << std::setw(WIDTH_FLOAT10) << std::fixed << std::setprecision(10) << mConductivityTensorInner(1,1) << " "
                           << "] "
                           << std::endl;

        mFileTensorResults << std::endl;
    }
}
