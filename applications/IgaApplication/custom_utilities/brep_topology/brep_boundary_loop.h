//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//


#if !defined(KRATOS_BREP_BOUNDARY_LOOP_H_INCLUDED )
#define  KRATOS_BREP_BOUNDARY_LOOP_H_INCLUDED


// System includes

// Project includes
#include "brep_trimming_curve.h"

#include "iga_application_variables.h"

namespace Kratos
{

    ///@name Kratos Classes
    ///@{
    /// Collection of all trimming curves which comply one full boundary loop.
    /** This class contains a collection of trimming curves that imply a boundary
    *   loop of the underlying surface. The loop can either be an inner or outer loop.
    */
    class BrepBoundaryLoop
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KratosNurbsBrepApplication
        KRATOS_CLASS_POINTER_DEFINITION(BrepBoundaryLoop);

        ///@}
        ///@name Life Cycle
        ///@{

        const std::vector<BrepTrimmingCurve>& GetTrimmingCurves() const;

        bool IsOuterLoop();

        //important function, to be implemented!!
        //std::vector<array_1d<double, 2>> GetBoundaryPolygon(unsigned int number_polygon_points);

        /// Constructor.
        BrepBoundaryLoop(
            std::vector<BrepTrimmingCurve>& rBrepTrimmingCurves,
            bool rIsOuterLoop);

        /// Destructor.
        virtual ~BrepBoundaryLoop() {};

        ///@}
    private:
        ///@name Member Variables
        ///@{
        std::vector<BrepTrimmingCurve> mBrepTrimmingCurves;
        bool mIsOuterLoop;
        ///@}
    }; // Class BrepBoundaryLoop

} // namespace Kratos.

#endif // KRATOS_BREP_BOUNDARY_LOOP_H_INCLUDED  defined