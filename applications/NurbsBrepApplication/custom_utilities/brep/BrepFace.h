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


#include "includes/define.h"
// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include "boost/make_shared.hpp"
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


#include "includes/model_part.h"

#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"

// ==============================================================================

namespace Kratos
{
  ///@name Type Definitions
  ///@{
  typedef std::vector<int> IntVector;
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
	/* Geometry Refinement Parameters are used to pass the full refinement information
	*  needed for face-refinement.
	*/
	struct GeometryRefinementParameters
	{
		/// knot variables
		Vector knot_insertions_u;
		Vector knot_insertions_v;
		int multiply_knots_u;
		int multiply_knots_v;
		double max_element_size_u;
		double max_element_size_v;
		/// degree variables
		int order_elevation_p;
		int order_elevation_q;
		int min_order_p;
		int min_order_q;

		/* Constructor */
		GeometryRefinementParameters() {
			knot_insertions_u = ZeroVector(0);
			knot_insertions_v = ZeroVector(0);

			multiply_knots_u = 0;
			multiply_knots_v = 0;

			max_element_size_u = 0.0;
			max_element_size_v = 0.0;

			order_elevation_p = 0;
			order_elevation_q = 0;

			min_order_p = 0;
			min_order_p = 0;
		}
	};

	struct EmbeddedPoint{
		int trim_index;
		Vector local_coordinates;

		EmbeddedPoint(const int& rTrimIndex, const Vector& rLocalCoordinates)
		{
			trim_index = rTrimIndex;
			local_coordinates = rLocalCoordinates;
		}
	};

	typedef std::vector<int>               IntVector;
	typedef std::vector<BrepBoundaryLoop>  TrimmingLoopVector;
	typedef std::vector<BrepTrimmingCurve> TrimmingCurveVector;

	/// Pointer definition of KratosNurbsBrepApplication
	KRATOS_CLASS_POINTER_DEFINITION(BrepFace);
    ///@}
    ///@name Life Cycle
    ///@{


	//Refinement functions
	void ApplyGeometryRefinement(const GeometryRefinementParameters& rGeometryRefinementParameters);

	// Mathematical projection functions
    bool ProjectionNewtonRaphson(const Point& point, double& u, double& v,
		const double& rAccuracy, const int& rMaxIterations);

	// Projection functions
	bool GetClosestPoint(const Point& rPoint, double& u, double& v,
		const double& rAccuracy, const int& rMaxIterations);
	std::vector<array_1d<double, 2>> GetClosestPointsTrimmingCurve(
		const std::vector<Point>& rPoints,
		const int& rTrimIndex,
		const double& rAccuracy, const int& rMaxIterations);
	std::vector<array_1d<double, 2>> GetClosestNodesTrimmingCurve(
		const std::vector<Node<3>::Pointer>& rNodes,
		const int& rTrimIndex,
		const double& rAccuracy, const double& rModelTolerance, const int& rMaxIterations);

	// Finds location of node with initial guess
	bool GetClosestIntegrationNode(
		Node<3>::Pointer& rClosestNode,
		const Node<3>::Pointer& rSpaceNode,
		const int& rShapefunctionOrder,
		const double& rAccuracy, const int& rMaxIterations);

	void GetClosestPolygonPointParameters(
		const std::vector<array_1d<double, 2>>& rPolygon,
		const Point& rPoint,
		double& rU, double& rV);


	// Add Variables
	//Node<3>::Pointer EvaluateNode(double u, double v, const int& shapefunction_order);
	// Enhance nicht ausreichend
	//void EnhanceNode(Node<3>::Pointer& node, const double& u, const double& v, const int& shapefunction_order);
	// gleiches Spiel - hier wird noch mehr gemacht wie beispielweise die Projektion auf die Oberfl�che
	// das muss auf jeden Fall generalisiert werden, wie beispielweise mit einer Funkion
	// Project node
    //void EnhanceShapeFunctionsSlave(
    //  std::vector<Node<3>::Pointer>& nodes, const int& trim_index, const int& shapefunction_order);
	// name is wrong - hier werden nodes aus points erstellt --
	//std::vector<Node<3>::Pointer> EnhanceShapeFunctions(std::vector<array_1d<double, 3>>& points, const int& shapefunction_order);
	// wird das �berhaupt gebraucht
	void EvaluateSurfacePoint(Point& rSurfacePoint, const double& u, const double& v);
	//void EvaluateShapeFunctionsSlaveNode(const double& u, const double& v, const int& shapefunction_order, Node<3>::Pointer node);


	// Integration domain
	//// not needed - deutlich mehr generalisieren
	//std::vector<Node<3>::Pointer> GetQuadraturePoints(const int& shapefunction_order);
	//std::vector<Node<3>::Pointer> GetQuadraturePointsTrimmed(const int& shapefunction_order);
	//std::vector<Node<3>::Pointer> GetQuadraturePointsEmbedded(const int& shapefunction_order);
	//std::vector<Node<3>::Pointer> GetQuadraturePointsReversed(const int& shapefunction_order);
	//// Curves - gleicher Name eventuell nur �berladen vor allem erste Funktion nochmal verwenden wenn m�glich
	//std::vector<Node<3>::Pointer> GetQuadraturePointsOfTrimmingCurve(const int& shapefunction_order, const int& trim_index);
	//std::vector<Node<3>::Pointer> GetQuadraturePointsOfTrimmingCurveWithPoints(
	//	const int& shapefunction_order, const int& trim_index, std::vector<Point> intersection_points);

	// Integration domain surfaces
	std::vector<Node<3>::Pointer> GetIntegrationNodesSurface(
		const int& rShapefunctionOrder, const int& rPolygonDiscretization);
	std::vector<Node<3>::Pointer> GetIntegrationNodesSurface(
		Polygon& rBoundaryPolygon,
		const int& rShapefunctionOrder);
	std::vector<Node<3>::Pointer> GetIntegrationNodesEmbedded(
		const int& rShapefunctionOrder, const int& rPolygonDiscretization);
	std::vector<Node<3>::Pointer> GetIntegrationNodesSurfaceReversed(
		const int& rShapefunctionOrder, const int& rPolygonDiscretization);

	// Integration domain trimming curves
	std::vector<Node<3>::Pointer> GetIntegrationNodesTrimmingCurve(
		const int& rShapefunctionOrder, const int& rTrimIndex, const double& rAccuracy);
	std::vector<Node<3>::Pointer> GetIntegrationNodesTrimmingCurveMaster(
		const std::vector<Point>& rPoints,
		const int& rShapefunctionOrder, const int& rTrimIndex, const int& rPQSlave,
		const double& rAccuracy, const double& rModelTolerance, const int& rMaxIterations);

	// Integration domain points
	Node<3>::Pointer GetIntegrationNodePoint(
		const int& rTrimIndex,
		const int& rShapefunctionOrder);

	void EvaluateIntegrationNodesTrimmingCurveSlave(
		std::vector<Node<3>::Pointer>& rNodes, const int& rShapefunctionOrder,
		const int& rTrimIndex,
		const double& rAccuracy, const double& rModelTolerance, const int& rMaxIterations);

	std::vector<Node<3>::Pointer> GetIntegrationNodes(
		std::vector<array_1d<double, 3>>& rPoints,
		const int& rShapefunctionOrder);

	std::vector<Node<3>::Pointer> GetIntegrationNodes(
		std::vector<array_1d<double, 4>>& rPoints, const int& rShapefunctionOrder);

    bool GetIntegrationNodeUpdated(
        Node<3>::Pointer& rClosestNode,
        const Node<3>::Pointer& rSpaceNode,
        const int& rShapefunctionOrder,
        const double& rAccuracy, const int& rMaxIterations);

	void EvaluateIntegrationNode(
		const double& rU, const double& rV,
		const int& rShapefunctionOrder,
		Node<3>::Pointer pNode);
	void EvaluateIntegrationNodeSlave(
		const double& rU, const double& rV,
		const int& rShapefunctionOrder,
		Node<3>::Pointer pNode);

	// Compute knot intersections
    std::vector<Point> GetKnotIntersections(const int& rTrimIndex);

	// Utilities
	bool CheckIfPointIsInside(const Vector& rNodeParameters);
	// zu irgendeiner utility Klasse?
	IntVector GetIntegerVector(const Vector& vector, const int& tolerance);

    //TODO: you need to give reading access to your internals through the Calculate function
	BrepTrimmingCurve GetTrimmingCurve(const int& trim_index);
	//Get functions
	Vector& GetUKnotVector();
	Vector& GetVKnotVector();
	int GetP();
	int GetQ();

    /// Constructor.
    BrepFace(unsigned int brep_id,
        bool is_trimmed, bool is_rational,
        TrimmingLoopVector& trimming_loops,
        TrimmingLoopVector& embedded_loops,
        std::vector<EmbeddedPoint>& embedded_points,
        Vector& knot_vector_u, Vector& knot_vector_v,
        unsigned int& p, unsigned int& q,
        IntVector& control_point_ids,
        Kratos::shared_ptr<ModelPart> model_part);

    /// Destructor.
    virtual ~BrepFace();

    ///@}
  protected:

private:

    ///@name Private methods
    ///@{
	//Refinement functions
	void RefineKnotVector(const Vector& rRu, const Vector& rRv);
	void DegreeElevate(const int& tp, const int& tq);


    // Geometry functions
    void EvaluateGradientsForClosestPointSearch(Vector QminP, Matrix& Hessian, Vector& Gradient, double& u, double& v);
    // All in one function!!
	void EvaluateNURBSFunctions(int span_u, int span_v, double _u, double _v, Matrix& R);
    void EvaluateNURBSFunctionsDerivatives(int span_u, int span_v, double _u, double _v, Matrix& _dR, Matrix& _ddR);
    void EvaluateNURBSFunctionsAndDerivative(int span_u, int span_v, double _u, double _v, Matrix& R, std::vector<Matrix>& dR);

    ///@}
    ///@name Member Variables
    ///@{
	bool m_is_trimmed;
	bool m_is_rational;
    TrimmingLoopVector m_trimming_loops;
    TrimmingLoopVector m_embedded_loops;
	std::vector<EmbeddedPoint> m_embedded_points;
    Vector m_knot_vector_u;
    Vector m_knot_vector_v;
    unsigned int m_p;
    unsigned int m_q;
    IntVector m_control_points_ids;
    Kratos::shared_ptr<ModelPart> mp_model_part;
    ///@}

  }; // Class BrepFace

}  // namespace Kratos.

#endif // KRATOS_BREP_FACE_H_INCLUDED  defined