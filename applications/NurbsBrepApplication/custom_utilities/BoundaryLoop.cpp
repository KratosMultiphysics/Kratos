//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//

// System includes


// External includes 


// Project includes
#include "BoundaryLoop.h"
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"


namespace Kratos
{
// --------------------------------------------------------------------------


  std::vector<array_1d<double, 2>> BoundaryLoop::GetBoundaryPolygon()
  {
    std::vector<array_1d<double, 2>> boundary_polygon;
    unsigned int number_polygon_points = 500;
    for (unsigned int curve_i = 0; curve_i < m_trimming_curves.size(); curve_i++)
    {
      std::vector<array_1d<double, 2>> boundary_polygon_edge = m_trimming_curves[curve_i].CreatePolygon(number_polygon_points);
      unsigned int old_length = boundary_polygon.size();
      boundary_polygon.resize(boundary_polygon.size() + number_polygon_points);
      for (unsigned int polygon_i = 0; polygon_i < boundary_polygon_edge.size(); polygon_i++)
      {
        boundary_polygon[old_length + polygon_i] = boundary_polygon_edge[polygon_i];
      }
    }
    return boundary_polygon;
  }

  //TrimmingCurveVector& BoundaryLoop::GetTrimmingCurves()
  //{
  //  return m_trimming_curves;
  //}
  bool& BoundaryLoop::IsOuterLoop()
  {
    return m_is_outer_loop;
  }
//Constructor
BoundaryLoop::BoundaryLoop(TrimmingCurveVector& trimming_curves, bool is_outer_loop)
  : m_trimming_curves(trimming_curves),
    m_is_outer_loop(is_outer_loop)
{
}
//Destructor
BoundaryLoop::~BoundaryLoop()
{}

}  // namespace Kratos.

