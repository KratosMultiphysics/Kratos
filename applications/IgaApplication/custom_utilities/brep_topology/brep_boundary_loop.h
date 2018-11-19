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
    loop of the underlying surface. The loop can either be an inner or outer loop.
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

        std::vector<BrepTrimmingCurve> GetTrimmingCurves();
        bool& IsOuterLoop();
        //important function, to be implemented!!
        //std::vector<array_1d<double, 2>> GetBoundaryPolygon(unsigned int number_polygon_points);

        /// Constructor.
        BrepBoundaryLoop(
            std::vector<BrepTrimmingCurve>& rBrepTrimmingCurves,
            bool rIsOuterLoop);

        /// Destructor.
        virtual ~BrepBoundaryLoop() {};

        ///@}
    protected:

    private:
        ///@name Member Variables
        ///@{
        std::vector<BrepTrimmingCurve> m_brep_trimming_curves;
        bool m_is_outer_loop;
        ///@}
    }; // Class BrepBoundaryLoop

} // namespace Kratos.

#endif // KRATOS_BREP_BOUNDARY_LOOP_H_INCLUDED  defined