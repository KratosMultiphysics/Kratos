#if !defined(KRATOS_BREP_FACE_H_INCLUDED )
#define  KRATOS_BREP_FACE_H_INCLUDED



// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "BrepTrimmingCurve.h"
#include "BrepBoundaryLoop.h"
#include "../knot_span/KnotSpan2d.h"
#include "../knot_span/KnotSpan2dNIntegrate.h"
#include "../nurbs_utilities.h"
#include "../../kratos/includes/node.h"
#include "../Polygon.h"

#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"

// ==============================================================================

namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 
  ///@} 
  ///@name Type Definitions
  ///@{ 
  typedef std::vector<int> IntVector;
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
  class BrepFace : public IndexedObject, public Flags
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<int> IntVector;
    typedef std::vector<BrepBoundaryLoop> TrimmingLoopVector;
    typedef std::vector<BrepTrimmingCurve> TrimmingCurveVector;

    /// Pointer definition of KratosNurbsBrepApplication
    KRATOS_CLASS_POINTER_DEFINITION(BrepFace);

    ///@}
    ///@name Life Cycle 
    ///@{ 

    //Get functions
    Vector& GetUKnotVector();
    Vector& GetVKnotVector();

    IntVector GetIntegerVector(const Vector& vector, const int& tolerance);
    //IntVector GetIntegerVKnotVector(const int& tolerance);
    //Closest Point Functions
    void MapNodeNewtonRaphson(const Node<3>::Pointer& node,
      Node<3>::Pointer& node_on_geometry);
    void GetProjectPoint(Node<3>::Pointer& node_on_geometry, const Node<3>::Pointer& node_location, const int& shapefunction_order);
    std::vector<Node<3>::Pointer> GetQuadraturePoints(const int& shapefunction_order);
    std::vector<Node<3>::Pointer> GetQuadraturePointsTrimmed(const int& shapefunction_order);
    std::vector<Node<3>::Pointer> GetQuadraturePointsOfTrimmingCurve(const int& shapefunction_order, const int& trim_index);
    std::vector<Node<3>::Pointer> GetQuadraturePointsOfTrimmingCurveWithPoints(
      const int& shapefunction_order, const int& trim_index, std::vector<Point<3>> intersection_points);
    bool CheckIfPointIsInside(Vector node_parameters);
    void EvaluateSurfacePoint(Point<3>& rSurfacePoint, const double& u, const double& v);
    void EvaluateShapeFunctionsSlaveNode(double& const u, double& const v, const int& shapefunction_order, Node<3>::Pointer node);
    Node<3>::Pointer EvaluateNode(double u, double v, const int& shapefunction_order);
    //void GetLocalParameterOfPoint(const Point<3>& point, double& u, double& v);
    void GetLocalParameterOfPointOnTrimmingCurve(const Point<3>& point, const BrepTrimmingCurve& trimming_curve, double& u, double& v);
    void GetClosestPoint(const Point<3>& point, double& u, double& v);
    bool NewtonRaphson(const Point<3>& point, double& u, double& v);
    std::vector<Node<3>::Pointer> EnhanceShapeFunctions(std::vector<array_1d<double, 3>>& points, const int& shapefunction_order);
    void EnhanceNode(Node<3>::Pointer& node, const double& u, const double& v, const int& shapefunction_order);
    void EnhanceShapeFunctionsSlave(
      std::vector<Node<3>::Pointer>& nodes, const int& trim_index, const int& shapefunction_order);
    BrepTrimmingCurve GetTrimmingCurve(const int& trim_index);
    std::vector<Point<3>> GetIntersectionPoints(const int& trim_index);
    //TODO: you need to give reading access to your internals through the Calculate function
    /// Constructor.
    BrepFace(unsigned int brep_id,
      TrimmingLoopVector& trimming_loops,
      Vector& knot_vector_u, Vector& knot_vector_v,
      unsigned int& p, unsigned int& q, IntVector& control_point_ids,
      ModelPart& model_part);

    /// Destructor.
    virtual ~BrepFace();

    /// Copy constructor.
    //BrepFace(BrepFace const& rOther);

    /// Assignment operator.
    //BrepFace& operator=(BrepFace const& rOther);
    ///@} 
  protected:

  private:

    ///@name Private methods
    ///@{ 
    //Geometry functions
    void EvaluateGradientsForClosestPointSearch(Vector QminP, Matrix& Hessian, Vector& Gradient, double& u, double& v);
    void EvaluateNURBSFunctions(int span_u, int span_v, double _u, double _v, Matrix& R);
    void EvaluateNURBSFunctionsDerivatives(int span_u, int span_v, double _u, double _v, Matrix& _dR, Matrix& _ddR);
    void EvaluateNURBSFunctionsAndDerivative(int span_u, int span_v, double _u, double _v, Matrix& R, std::vector<Matrix>& dR);

    ///@} 
    ///@name Member Variables
    ///@{ 

    TrimmingLoopVector m_trimming_loops;
    Vector m_knot_vector_u;
    Vector m_knot_vector_v;
    unsigned int m_p;
    unsigned int m_q;
    IntVector m_control_points_ids;
    ModelPart& m_model_part;

    ///@}    

  }; // Class BrepFace 

}  // namespace Kratos.

#endif // KRATOS_BREP_FACE_H_INCLUDED  defined