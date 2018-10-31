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
#inclue "anurbs.h"

// Project includes
#include "brep_trimming_curve.h"
#include "iga_application.h"
#include "iga_application_variables.h"


namespace Kratos
{
    unsigned int& BrepTrimmingCurve::GetTrimIndex()
    {
        return m_trim_index;
    }

    //Constructor
    BrepTrimmingCurve::BrepTrimmingCurve(unsigned int trim_index, bool curve_direction, Vector& knot_vector_u,
        unsigned int p, ControlPointVector& control_points,
        Vector& active_range)
        : m_knot_vector_u(knot_vector_u),
        m_curve_direction(curve_direction),
        m_p(p),
        m_control_points(control_points),
        m_active_range(active_range),
        m_trim_index(trim_index)
    {
    }

}  // namespace Kratos.

