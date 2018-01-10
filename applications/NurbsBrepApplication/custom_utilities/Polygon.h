#if !defined(KRATOS_POLYGON_H_INCLUDED )
#define  KRATOS_POLYGON_H_INCLUDED


// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <vector>
#include <list>

#include <deque>

#include <boost/geometry.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>
#include <boost/geometry/algorithms/reverse.hpp> 

//#include <boost/geometry/geometries/cartesian2d.hpp>
//#include <boost/geometry/geometries/adapted/c_array_cartesian.hpp>

#include <boost/foreach.hpp>
// ==============================================================================
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"

#include "brep/BrepBoundaryLoop.h"

namespace Kratos
{
  //dynamic programming state for minimum-weight triangulation
  struct DPState {
    bool visible;
    double weight;
    long bestvertex;
  };
  struct Diagonal {
    long index1;
    long index2;
  };
  ///@name Kratos Globals
  ///@{ 
  ///@} 
  ///@name Type Definitions
  ///@{ 
  ///@}
  ///@name  Enum's
  ///@{
  ///@}
  ///@name  Functions 
  ///@{
  ///@}
  ///@name Kratos Classes
  ///@{
  /// Short class definition.
  /** Detail class definition.
  */
  class Polygon
  {
  public:
    ///@name Type Definitions
    ///@{
    //using namespace boost::geometry;

    typedef boost::geometry::model::d2::point_xy<double> PointXYType;
    typedef boost::geometry::model::polygon<PointXYType> PolygonType;
    typedef std::vector<PolygonType> PolygonVectorType;

    /// Pointer definition of KratosNurbsTestcaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(Polygon);
    

    std::vector<Matrix> Triangulate();
    ///@}
    ///@name Life Cycle 
    ///@{
    Polygon clipByKnotSpan(const Vector& parameter_span_u, const Vector& parameter_span_v);

    bool IsInside(const double& u, const double& v);

    double GetArea();

    bool IsFullKnotSpan();
    bool m_is_full_knot_span = false;


    //TODO: you need to give reading access to your internals through the Calculate function
    /// Constructor.
    Polygon(PolygonType polygon);
    Polygon(PolygonVectorType polygon);
    Polygon(std::vector<array_1d<double, 2>> polygon);
    Polygon(std::vector<BrepBoundaryLoop>& boundary_loops);


    /// Destructor.
    virtual ~Polygon();

    /// Copy constructor.
    //Polygon(Polygon const& rOther);

    /// Assignment operator.
    //Polygon& operator=(Polygon const& rOther);
    ///@} 
  protected:

  private:


    ///@name Private methods
    ///@{ 
    bool Triangulate_OPT(const PolygonType& polygon,
      std::vector<Matrix>& triangles);

    void clipPolygon(
      const PolygonType& polygon_1,
      const PolygonType& polygon_2,
      PolygonVectorType& polygon_vector);

    void Reverse(const int& index);

    double Distance(array_1d<double, 2> point_1, array_1d<double, 2> point_2);
    double Distance(PointXYType point_1, PointXYType point_2);

    bool Intersects(array_1d<double, 2> &p11, array_1d<double, 2> &p12,
      array_1d<double, 2> &p21, array_1d<double, 2> &p22);
    bool Intersects(PointXYType &p11, PointXYType &p12,
      PointXYType &p21, PointXYType &p22);

    bool InCone(array_1d<double, 2> &p1, array_1d<double, 2> &p2,
      array_1d<double, 2> &p3, array_1d<double, 2> &p);
    bool InCone(PointXYType &p1, PointXYType &p2,
      PointXYType &p3, PointXYType &p);

    bool IsConvex(
      const array_1d<double, 2>& p1, const array_1d<double, 2>& p2,
      const array_1d<double, 2>& p3);
    bool IsConvex(
      const PointXYType& p1, const PointXYType& p2,
      const PointXYType& p3);
    bool GetOrientation();
    void Invert();
    double GetAreaOfTriangle(const Matrix& triangle);
    ///@} 
    ///@name Member Variables
    ///@{ 

    PolygonVectorType m_polygon_list;

    ///@}    

  }; // Class Polygon 

}  // namespace Kratos.

#endif // KRATOS_POLYGON_H_INCLUDED  defined