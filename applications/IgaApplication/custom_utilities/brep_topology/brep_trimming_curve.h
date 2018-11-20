#if !defined(KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED )
#define  KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED


// System includes

// External includes
//#pragma once
//#include <ANurbs/src/Curve.h>
//#include <ANurbs/src/CurveGeometry.h>
//#include <ANurbs/Integration>
#include "custom_utilities/anurbs.h"
//#include "../node_surface_geometry_3d.h"
//#include <ANurbs/Integration>
//#include "anurbs.h"

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

        const Kratos::shared_ptr<Curve<2>> GetCurve2D() const;

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

        std::shared_ptr<CurveGeometry<2>> m_geometry;
        std::shared_ptr<Curve<2>> m_curve;
        ///@}
        ///@name Un accessible methods
        ///@{


        ///@}

    }; // Class BrepTrimmingCurve

}  // namespace Kratos.

#endif // KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED  defined