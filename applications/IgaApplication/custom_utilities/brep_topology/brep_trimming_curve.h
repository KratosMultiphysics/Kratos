#if !defined(KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED )
#define  KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED


// System includes

// External includes
#include <ANurbs/src/Curve.h>

// Project includes
#include "iga_application_variables.h"


namespace Kratos
{
    ///@name Kratos Classes
    ///@{
    class BrepTrimmingCurve
    {


    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KratosNurbsBrepApplication
        KRATOS_CLASS_POINTER_DEFINITION(BrepTrimmingCurve);

        ///@}
        ///@name Life Cycle
        ///@{
        // Utilities
        /* Returns trimming curve index
        */
        int& GetTrimIndex();

        std::shared_ptr<ANurbs::Curve2D> GetCurve2D();

        /// Constructor.
        BrepTrimmingCurve(
            int& rTrimIndex,
            Vector& rKnotVector,
            int& rDegree,
            std::vector<BoundedVector<double, 4>>& rControlPoints,
            bool rCurveDirection,
            bool rIsRational,
            Vector& rActiveRange);

        /// Destructor.
        virtual ~BrepTrimmingCurve() {};

        ///@}
    protected:

    private:
        ///@name Member Variables
        ///@{
        int m_trim_index;
        bool m_curve_direction;

        std::shared_ptr<ANurbs::CurveGeometry2D> m_geometry;
        std::shared_ptr<ANurbs::Curve2D> m_curve;
        ///@}
        ///@name Un accessible methods
        ///@{


        ///@}

    }; // Class BrepTrimmingCurve

}  // namespace Kratos.

#endif // KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED  defined