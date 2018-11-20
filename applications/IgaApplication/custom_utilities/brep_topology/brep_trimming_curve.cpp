//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//               Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//

//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//

// System includes

// External includes 

// Project includes
#include "brep_trimming_curve.h"


namespace Kratos
{
    int& BrepTrimmingCurve::GetTrimIndex()
    {
        return m_trim_index;
    }

    const Kratos::shared_ptr<ANurbs::Curve<Kratos::array_1d<double, 2>>> BrepTrimmingCurve::GetCurve2D() const
    {
        return m_curve;
    }

    ///Constructor
    BrepTrimmingCurve::BrepTrimmingCurve(
        int& rTrimIndex,
        Vector& rKnotVector,
        int& rDegree,
        std::vector<BoundedVector<double, 4>>& rControlPoints,
        bool rCurveDirection,
        bool rIsRational,
        Vector& rActiveRange)
        : m_curve_direction(rCurveDirection),
          m_trim_index(rTrimIndex)
    {
        int number_poles = rControlPoints.size();

        m_geometry = Kratos::make_shared<ANurbs::CurveGeometry2D>(
            rDegree, number_poles, rIsRational);

        //for (int i = 0; i < rKnotVector.size() - 2; ++i)
        //{
        //    m_geometry->SetKnot(i, rKnotVector(i + 1));
        //}

        //for (int i = 0; i < number_poles; ++i)
        //{
            //Kratos::array_1d<double, 2> cp;
            //cp[0] = rControlPoints[i][0];
            //cp[1] = rControlPoints[i][1];
            //m_geometry->SetPole(i, cp);
            //if (rIsRational)
            //{
            //    m_geometry->SetWeight(i, rControlPoints[i][3]);
            //}
        //}

        m_curve = Kratos::make_unique<ANurbs::Curve<Kratos::array_1d<double, 2>>>(m_geometry, rActiveRange[0], rActiveRange[1]);
    }

}  // namespace Kratos.

