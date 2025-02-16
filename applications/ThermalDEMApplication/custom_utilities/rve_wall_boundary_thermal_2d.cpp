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
        RVEWallBoundary2D::ProcessGlobalResults();
        // TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    void RVEWallBoundaryThermal2D::WriteFileHeaders(void) {
        RVEWallBoundary2D::WriteFileHeaders();
        // TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    //------------------------------------------------------------------------------------------------------------
    void RVEWallBoundaryThermal2D::WriteResultFiles(void) {
        RVEWallBoundary2D::WriteResultFiles();
        // TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
}
