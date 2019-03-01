//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//               Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//

// System includes

// External includes 

// Project includes
#include "brep_boundary_loop.h"


namespace Kratos
{
    const std::vector<BrepTrimmingCurve>& BrepBoundaryLoop::GetTrimmingCurves() const
    {
        return mBrepTrimmingCurves;
    }


    bool BrepBoundaryLoop::IsOuterLoop()
    {
        return mIsOuterLoop;
    }

    /// Constructor.
    BrepBoundaryLoop::BrepBoundaryLoop(
        std::vector<BrepTrimmingCurve>& rBrepTrimmingCurves,
        bool rIsOuterLoop)
        : mBrepTrimmingCurves(rBrepTrimmingCurves),
          mIsOuterLoop(rIsOuterLoop)
    {
    }
}  // namespace Kratos.

