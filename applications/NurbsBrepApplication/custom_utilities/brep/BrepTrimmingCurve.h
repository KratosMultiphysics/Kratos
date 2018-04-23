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
  ///@name Kratos Classes
  ///@{
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
	// Discretization function
    std::vector<array_1d<double, 2>> CreatePolygon(unsigned int number_polygon_points);
    std::vector<array_1d<double, 3>> CreatePolygonWithParameter(unsigned int number_polygon_points);

	// Integration functions
    std::vector<array_1d<double, 3>> GetQuadraturePoints(std::vector<double> span, double polynomial_order_p);
	std::vector<array_1d<double, 4>> GetIntegrationPoints(const std::vector<double>& rSpan, const double& rDegree);

	// Knot intersections for entire curve:
    std::vector<double> GetKnotIntersections(const int& p, const int& q, const Vector& knot_vector_u, const Vector& knot_vector_v, const int& rNumberOfPolygonPoints);

	// Closest Point projections:
	bool GetClosestPoint(const array_1d<double, 2>& closest_point, double& parameter);
	bool GetClosestPoint(
		const array_1d<double, 2>& rClosestPoint,
		const double& rKnotPercentage,
		double& rParameter);
	std::vector<double> GetClosestPoints(std::vector<array_1d<double, 2>> points, const int& rNumberOfPolygonPoints);

	// Projections towards curve
    bool ProjectionNewtonRaphson(double& parameter, const array_1d<double, 2>& closest_point);
    bool ProjectionBisection(double& parameter, const array_1d<double, 2>& closest_point);
    bool ProjectionBisection(double& parameter, const array_1d<double, 2>& closest_point, double parameter_min, double parameter_max);

	// Knot intersections for specific knot
    bool GetKnotIntersection(double& parameter, const int& direction, const double& knot);
    bool GetKnotIntersectionBisection(double& parameter, const double& parameter_1, const double& parameter_2,
      const double& direction, const double& knot);

	// Geoemetrical functions
	void EvaluateCurvePoint(Point& rCurvePoint, double parameter_u);
	void EvaluateCurveDerivatives(std::vector<Vector>& derivatives, const int& order, const double& u);

	// Shape Functions
    void EvaluateRationalCurveDerivatives(std::vector<Vector>& DN_De, const int& order, const double& u);

	void RefineKnotVector(const Vector& rRu);

	// Utilities
	unsigned int& GetIndex();

    //TODO: you need to give reading access to your internals through the Calculate function
    /// Constructor.
    //TODO: pass by reference not by value
    BrepTrimmingCurve(unsigned int trim_index, bool curve_direction, Vector& knot_vector_u,
      unsigned int p, ControlPointVector& control_points,
      Vector& active_range);

    /// Destructor.
    virtual ~BrepTrimmingCurve();

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