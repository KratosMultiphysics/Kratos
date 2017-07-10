#if !defined(KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED )
#define  KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED


// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"

#include "../nurbs_utilities.h"
#include "../knot_span/KnotSpan1d.h"


namespace Kratos
{

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
  class BrepTrimmingCurve
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<array_1d<double, 4>> ControlPointVector;
    
    /// Pointer definition of KratosNurbsBrepApplication
    KRATOS_CLASS_POINTER_DEFINITION(BrepTrimmingCurve);

    ///@}
    ///@name Life Cycle 
    ///@{ 
    std::vector<array_1d<double, 2>> CreatePolygon(unsigned int number_polygon_points);
    std::vector<array_1d<double, 3>> CreatePolygonWithParameter(unsigned int number_polygon_points);
    void EvaluateCurvePoint(Point<3>& rCurvePoint, double parameter_u);
    unsigned int& GetIndex();
    void PrintData();


    std::vector<array_1d<double, 3>> GetQuadraturePoints(std::vector<double> span, double polynomial_order_p);
    std::vector<double> FindIntersections(const int& p, const int& q, const Vector& knot_vector_u, const Vector& knot_vector_v);
    std::vector<double> FindIntersectionsWithPoints(std::vector<Point<2>> intersection_points);
    double EvaluateIntersection(double initial_u, int intersection_base, const double& coordinate_base);
    array_1d<double, 2> GetBaseVector(const int& u);
    void GetCurveDerivatives(std::vector<Vector>& derivatives, const int& order, const double& u);
    void EvaluateCurveDerivatives(Matrix& DN_De, const int& order, const double& u);
    void EvaluateLocalParameter(double& parameter, bool& converged, int baseVec,
    double baseComp, double uinit, double itmax, double iteps);
    bool GetClosestPoint(const Point<2>& closest_point, double& parameter);
    bool ProjectionNewtonRaphson(double& parameter, const Point<2>& closest_point);
    bool ProjectionBisection(double& parameter, const Point<2>& closest_point);
    bool ProjectionBisection(double& parameter, const Point<2>& closest_point, double parameter_min, double parameter_max);
    bool GetClosestPointBisection(double& parameter, const Point<2>& closest_point);
    bool GetKnotIntersection(double& parameter, const int& direction, const double& knot);
    bool Bisection(double& parameter, const double& parameter_1, const double& parameter_2,
      const double& direction, const double& knot);

    void EvaluateRationalCurveDerivatives(std::vector<Vector>& DN_De, const int& order, const double& u);


    //Shape Functions:
    void EvaluateRationalCurveDerivativesPrint(std::vector<Vector>& DN_De, const int& order, const double& u);
    void EvaluateCurveDerivativesPrint(Matrix& DN_De, const int& order, const double& u);


    //TODO: you need to give reading access to your internals through the Calculate function
    /// Constructor.
    //TODO: pass by reference not by value
    //TODO: why control points have size 4??? pass them as kratos nodes
    BrepTrimmingCurve(unsigned int trim_index, bool curve_direction, Vector& knot_vector_u,
      unsigned int p, ControlPointVector& control_points,
      Vector& active_range);

    /// Destructor.
    virtual ~BrepTrimmingCurve();

    /// Copy constructor.
    //BrepTrimmingCurve(BrepTrimmingCurve const& rOther);

    /// Assignment operator.
    //BrepTrimmingCurve& operator=(BrepTrimmingCurve const& rOther);
    ///@} 
  protected:

  private:
    ///@name Member Variables
    ///@{ 
    unsigned int m_trim_index;
    bool m_curve_direction;
    Vector m_knot_vector_u;
    unsigned int m_p;
    ControlPointVector m_control_points;
    Vector m_active_range;
    ///@}    
     
    ///@name Un accessible methods 
    ///@{ 



    ///@}    

  }; // Class BrepTrimmingCurve 

}  // namespace Kratos.

#endif // KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED  defined