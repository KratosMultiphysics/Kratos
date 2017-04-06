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
#include "../nurbs_utilities.h"
#include "../../kratos/includes/node.h"

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
    IntVector GetIntegerUKnotVector(const int& tolerance);
    IntVector GetIntegerVKnotVector(const int& tolerance);

    //Closest Point Functions
    void MapNodeNewtonRaphson(const Node<3>::Pointer& node,
      Node<3>::Pointer& node_on_geometry);
    std::vector<Node<3>::Pointer> GetQuadraturePoints(const int& shapefunction_order);
    bool CheckIfPointIsInside(Vector node_parameters);
    void EvaluateSurfacePoint(Node<3>::Pointer& rSurfacePoint, double u, double v);
    Node<3>::Pointer EvaluateNode(double u, double v, const int& shapefunction_order);
    std::vector<Node<3>::Pointer> EnhanceShapeFunctions(std::vector<array_1d<double, 3>>& points, const int& shapefunction_order);


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