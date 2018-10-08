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

// Project includes
#include "BrepFace.h"

//#include "define.h"

namespace Kratos
{
	// --------------------------------------------------------------------------
	Vector& BrepFace::GetUKnotVector()
	{
	return m_knot_vector_u;
	}
	Vector& BrepFace::GetVKnotVector()
	{
	return m_knot_vector_v;
	}
	int BrepFace::GetP()
	{
		return (int) m_p;
	}
	int BrepFace::GetQ()
	{
		return (int) m_q;
	}

	BrepTrimmingCurve BrepFace::GetTrimmingCurve(const int& trim_index)
	{
		for (unsigned int i = 0; i < m_trimming_loops.size(); i++)
		{
			std::vector<BrepTrimmingCurve> trimming_curves = m_trimming_loops[i].GetTrimmingCurves();
			for (unsigned int j = 0; j < trimming_curves.size(); j++)
			{
				if (trimming_curves[j].GetIndex() == trim_index)
				{
					return trimming_curves[j];
				}
			}
		}
		KRATOS_ERROR << "Brep Trimming Curve with index " << trim_index << " was not found." << std::endl;
	}

	IntVector BrepFace::GetIntegerVector(const Vector& vector, const int& tolerance)
	{
		IntVector new_vector;
		new_vector.resize(vector.size());

		for (unsigned int i = 0; i < vector.size(); i++)
		{
			new_vector[i] = (int)(vector[i] * tolerance);
		}

		return new_vector;
	}

	// --------------------------------------------------------------------------
	// to be deleted
 // std::vector<Node<3>::Pointer> BrepFace::GetQuadraturePoints(const int& shapefunction_order)
	//{
	//	int tolerance = 10e7;
	//
	//	std::vector<Node<3>::Pointer> NodeVector;
	//
	//	IntVector knot_vector_u = GetIntegerVector(m_knot_vector_u, tolerance);
	//	IntVector knot_vector_v = GetIntegerVector(m_knot_vector_v, tolerance);
	//
	//	Vector parameter_span_u = ZeroVector(2);
	//	Vector parameter_span_v = ZeroVector(2);
	//
	//	for (unsigned int i = 0; i < knot_vector_u.size() - 1; i++)
	//	{
	//		if (abs(knot_vector_u[i + 1] - knot_vector_u[i]) > 1)
	//		{
	//			parameter_span_u[0] = m_knot_vector_u[i];
	//			parameter_span_u[1] = m_knot_vector_u[i + 1];
	//
	//			for (unsigned int j = 0; j < knot_vector_v.size() - 1; j++)
	//			{
	//				if (abs(knot_vector_v[j + 1] - knot_vector_v[j]) > 1)
	//				{
	//					parameter_span_v[0] = m_knot_vector_v[j];
	//					parameter_span_v[1] = m_knot_vector_v[j + 1];
	//
	//					KnotSpan2d knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v);
	//					std::vector<array_1d<double, 3>> points = knot_span.getIntegrationPointsInParameterDomain();
	//
	//					std::vector<Node<3>::Pointer> NodeVectorElement = GetIntegrationNodes(points, shapefunction_order);
	//					for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
	//					{
	//						NodeVector.push_back(NodeVectorElement[k]);
	//					}
	//				}
	//			}
	//		}
	//	}
	//	return NodeVector;
	//}
	// to be deleted
	//std::vector<Node<3>::Pointer> BrepFace::GetQuadraturePointsTrimmed(const int& shapefunction_order)
	//{
	//	Polygon boundaries(m_trimming_loops);
	//
	//	std::vector<Node<3>::Pointer> NodeVector;
	//
	//	int tolerance = 10e6;
	//
	//	IntVector knot_vector_u = GetIntegerVector(m_knot_vector_u, tolerance);
	//	IntVector knot_vector_v = GetIntegerVector(m_knot_vector_v, tolerance);
	//
	//	Vector parameter_span_u = ZeroVector(2);
	//	Vector parameter_span_v = ZeroVector(2);
	//
	//	for (unsigned int i = 0; i < knot_vector_u.size() - 1; i++)
	//	{
	//		if (abs(knot_vector_u[i + 1] - knot_vector_u[i]) > 1)
	//		{
	//			parameter_span_u[0] = m_knot_vector_u[i];
	//			parameter_span_u[1] = m_knot_vector_u[i + 1];
	//
	//			for (unsigned int j = 0; j < knot_vector_v.size() - 1; j++)
	//			{
	//				if (abs(knot_vector_v[j + 1] - knot_vector_v[j]) > 1)
	//				{
	//					parameter_span_v[0] = m_knot_vector_v[j];
	//					parameter_span_v[1] = m_knot_vector_v[j + 1];
	//
	//					Polygon boundary_polygon = boundaries.clipByKnotSpan(parameter_span_u, parameter_span_v);
	//
	//					std::vector<array_1d<double, 3>> points;
	//					if (!boundary_polygon.IsFullKnotSpan())
	//					{
	//						KnotSpan2dNIntegrate knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v, boundary_polygon);
	//						points = knot_span.getIntegrationPointsInParameterDomain();
	//					}
	//					else
	//					{
	//						KnotSpan2d knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v);
	//						points = knot_span.getIntegrationPointsInParameterDomain();
	//					}
	//
	//					std::vector<Node<3>::Pointer> NodeVectorElement = GetIntegrationNodes(points, shapefunction_order);
	//					for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
	//					{
	//						NodeVector.push_back(NodeVectorElement[k]);
	//					}
	//				}
	//			}
	//		}
	//	}
	//	return NodeVector;
	//}
	// to be deleted
	//std::vector<Node<3>::Pointer> BrepFace::GetQuadraturePointsEmbedded(const int& shapefunction_order)
	//{
	//	Polygon boundaries(m_embedded_loops);
	//
	//	std::vector<Node<3>::Pointer> NodeVector;
	//
	//	int tolerance = 10e6;
	//
	//	IntVector knot_vector_u = GetIntegerVector(m_knot_vector_u, tolerance);
	//	IntVector knot_vector_v = GetIntegerVector(m_knot_vector_v, tolerance);
	//
	//	Vector parameter_span_u = ZeroVector(2);
	//	Vector parameter_span_v = ZeroVector(2);
	//
	//	for (unsigned int i = 0; i < knot_vector_u.size() - 1; i++)
	//	{
	//		if (abs(knot_vector_u[i + 1] - knot_vector_u[i]) > 1)
	//		{
	//			parameter_span_u[0] = m_knot_vector_u[i];
	//			parameter_span_u[1] = m_knot_vector_u[i + 1];
	//
	//			for (unsigned int j = 0; j < knot_vector_v.size() - 1; j++)
	//			{
	//				if (abs(knot_vector_v[j + 1] - knot_vector_v[j]) > 1)
	//				{
	//					parameter_span_v[0] = m_knot_vector_v[j];
	//					parameter_span_v[1] = m_knot_vector_v[j + 1];
	//
	//					Polygon boundary_polygon = boundaries.clipByKnotSpan(parameter_span_u, parameter_span_v);
	//
	//					std::vector<array_1d<double, 3>> points;
	//					if (!boundary_polygon.IsFullKnotSpan())
	//					{
	//						KnotSpan2dNIntegrate knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v, boundary_polygon);
	//						points = knot_span.getIntegrationPointsInParameterDomain();
	//					}
	//					else
	//					{
	//						KnotSpan2d knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v);
	//						points = knot_span.getIntegrationPointsInParameterDomain();
	//					}
	//					std::vector<Node<3>::Pointer> NodeVectorElement = GetIntegrationNodes(points, shapefunction_order);
	//					for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
	//					{
	//						NodeVector.push_back(NodeVectorElement[k]);
	//					}
	//				}
	//			}
	//		}
	//	}
	//	return NodeVector;
	//}
//to be deleted
  //std::vector<Node<3>::Pointer> BrepFace::GetQuadraturePointsReversed(const int& shapefunction_order)
  //{
	 // Polygon boundaries(m_embedded_loops);
	 // Polygon full_boundaries(m_trimming_loops);
	 // Polygon difference = full_boundaries.GetDifference(boundaries);
	 //
	 // std::vector<Node<3>::Pointer> NodeVector;
	 //
	 // int tolerance = 10e6;
	 //
	 // IntVector knot_vector_u = GetIntegerVector(m_knot_vector_u, tolerance);
	 // IntVector knot_vector_v = GetIntegerVector(m_knot_vector_v, tolerance);
	 //
	 // Vector parameter_span_u = ZeroVector(2);
	 // Vector parameter_span_v = ZeroVector(2);
	 //
	 // for (unsigned int i = 0; i < knot_vector_u.size() - 1; i++)
	 // {
		//  if (abs(knot_vector_u[i + 1] - knot_vector_u[i]) > 1)
		//  {
		//	  parameter_span_u[0] = m_knot_vector_u[i];
		//	  parameter_span_u[1] = m_knot_vector_u[i + 1];
		//
		//	  for (unsigned int j = 0; j < knot_vector_v.size() - 1; j++)
		//	  {
		//		  if (abs(knot_vector_v[j + 1] - knot_vector_v[j]) > 1)
		//		  {
		//			  parameter_span_v[0] = m_knot_vector_v[j];
		//			  parameter_span_v[1] = m_knot_vector_v[j + 1];
		//
		//			  Polygon boundary_polygon = difference.clipByKnotSpan(parameter_span_u, parameter_span_v);
		//
		//			  std::vector<array_1d<double, 3>> points;
		//			  if (!boundary_polygon.IsFullKnotSpan())
		//			  {
		//				  KnotSpan2dNIntegrate knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v, boundary_polygon);
		//				  points = knot_span.getIntegrationPointsInParameterDomain();
		//			  }
		//			  else
		//			  {
		//				  KnotSpan2d knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v);
		//				  points = knot_span.getIntegrationPointsInParameterDomain();
		//			  }
		//			  std::vector<Node<3>::Pointer> NodeVectorElement = GetIntegrationNodes(points, shapefunction_order);
		//			  for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
		//			  {
		//				  NodeVector.push_back(NodeVectorElement[k]);
		//			  }
		//		  }
		//	  }
		//  }
	 // }
	 // return NodeVector;
  //}
  // to be deleted
  //std::vector<Node<3>::Pointer> BrepFace::GetQuadraturePointsOfTrimmingCurve(const int& shapefunction_order, const int& trim_index)
  //{
  //  BrepTrimmingCurve trimming_curve = GetTrimmingCurve(trim_index);
  //
  //  std::vector<double> intersections = trimming_curve.GetKnotIntersections(m_p, m_q, m_knot_vector_u, m_knot_vector_v, 200);
  //
  //  int degree = std::max(m_p, m_q);
  //
  //  std::vector<array_1d<double, 3>> quadrature_points = trimming_curve.GetQuadraturePoints(intersections, degree);
  //
  //  std::vector<Node<3>::Pointer> NodeVectorElement = GetIntegrationNodes(quadrature_points, shapefunction_order);
  //  for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
  //  {
		//std::vector<Vector> derivatives;
  //    trimming_curve.EvaluateCurveDerivatives(derivatives, 2, quadrature_points[k][2]);
  //    Vector tangents_basis_vector(2);
  //    tangents_basis_vector[0] = derivatives[1][0];
  //    tangents_basis_vector[1] = derivatives[1][1];
  //    NodeVectorElement[k]->SetValue(TANGENTS_BASIS_VECTOR, tangents_basis_vector);
  //  }
  //  return NodeVectorElement;
  //}

	/* Compues the surface integration points. Consideres the possibility of untrimmed surfaces. */
	std::vector<Node<3>::Pointer> BrepFace::GetIntegrationNodesSurface(const int& rShapefunctionOrder, const int& rPolygonDiscretization)
	{
		Polygon boundaries(m_trimming_loops, rPolygonDiscretization);

		std::vector<Node<3>::Pointer> NodeVector;

		int tolerance = 10e6;

		IntVector knot_vector_u = GetIntegerVector(m_knot_vector_u, tolerance);
		IntVector knot_vector_v = GetIntegerVector(m_knot_vector_v, tolerance);

		Vector parameter_span_u = ZeroVector(2);
		Vector parameter_span_v = ZeroVector(2);

		for (unsigned int i = 0; i < knot_vector_u.size() - 1; i++)
		{
			if (std::abs(knot_vector_u[i + 1] - knot_vector_u[i]) > 1)
			{
				parameter_span_u[0] = m_knot_vector_u[i];
				parameter_span_u[1] = m_knot_vector_u[i + 1];

				for (unsigned int j = 0; j < knot_vector_v.size() - 1; j++)
				{
					if (std::abs(knot_vector_v[j + 1] - knot_vector_v[j]) > 1)
					{
						parameter_span_v[0] = m_knot_vector_v[j];
						parameter_span_v[1] = m_knot_vector_v[j + 1];

						if (!m_is_trimmed)
						{
							KnotSpan2d knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v);
							std::vector<array_1d<double, 3>> points = knot_span.getIntegrationPointsInParameterDomain();

							std::vector<Node<3>::Pointer> NodeVectorElement = GetIntegrationNodes(points, rShapefunctionOrder);
							for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
							{
								NodeVector.push_back(NodeVectorElement[k]);
							}
						}
						else
						{
							Polygon boundary_polygon = boundaries.clipByKnotSpan(parameter_span_u, parameter_span_v);

							std::vector<array_1d<double, 3>> points;
							if (!boundary_polygon.IsFullKnotSpan())
							{
								KnotSpan2dNIntegrate knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v, boundary_polygon);
								points = knot_span.getIntegrationPointsInParameterDomain();
							}
							else
							{
								KnotSpan2d knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v);
								points = knot_span.getIntegrationPointsInParameterDomain();
							}
							std::vector<Node<3>::Pointer> NodeVectorElement = GetIntegrationNodes(points, rShapefunctionOrder);
							for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
							{
								NodeVector.push_back(NodeVectorElement[k]);
							}
						}
					}
				}
			}
		}
		return NodeVector;
	}

	/* Compues the surface integration points.*/
	std::vector<Node<3>::Pointer> BrepFace::GetIntegrationNodesSurface(Polygon& rBoundaryPolygon, const int& rShapefunctionOrder)
	{
		std::vector<Node<3>::Pointer> NodeVector;

		int tolerance = 10e6;

		IntVector knot_vector_u = GetIntegerVector(m_knot_vector_u, tolerance);
		IntVector knot_vector_v = GetIntegerVector(m_knot_vector_v, tolerance);

		Vector parameter_span_u = ZeroVector(2);
		Vector parameter_span_v = ZeroVector(2);

		for (unsigned int i = 0; i < knot_vector_u.size() - 1; i++)
		{
			if (std::abs(knot_vector_u[i + 1] - knot_vector_u[i]) > 1)
			{
				parameter_span_u[0] = m_knot_vector_u[i];
				parameter_span_u[1] = m_knot_vector_u[i + 1];

				for (unsigned int j = 0; j < knot_vector_v.size() - 1; j++)
				{
					if (std::abs(knot_vector_v[j + 1] - knot_vector_v[j]) > 1)
					{
						parameter_span_v[0] = m_knot_vector_v[j];
						parameter_span_v[1] = m_knot_vector_v[j + 1];

						Polygon boundary_polygon = rBoundaryPolygon.clipByKnotSpan(parameter_span_u, parameter_span_v);

						std::vector<array_1d<double, 3>> points;
						if (!boundary_polygon.IsFullKnotSpan())
						{
							KnotSpan2dNIntegrate knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v, boundary_polygon);
							points = knot_span.getIntegrationPointsInParameterDomain();
						}
						else
						{
							KnotSpan2d knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v);
							points = knot_span.getIntegrationPointsInParameterDomain();
						}
						std::vector<Node<3>::Pointer> NodeVectorElement = GetIntegrationNodes(points, rShapefunctionOrder);
						for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
						{
							NodeVector.push_back(NodeVectorElement[k]);
						}
					}
				}
			}
		}
		return NodeVector;
	}



	/* Compues the surface integration points which lay inside the embedded loop. */
	std::vector<Node<3>::Pointer> BrepFace::GetIntegrationNodesEmbedded(const int& rShapefunctionOrder, const int& rPolygonDiscretization)
	{
		Polygon boundaries(m_embedded_loops, rPolygonDiscretization);

		return GetIntegrationNodesSurface(boundaries, rShapefunctionOrder);
	}

	/* Compues the surface integration points which lays outside of the the trimmed domain
	*  but still inside the domain of the patch. */
	std::vector<Node<3>::Pointer> BrepFace::GetIntegrationNodesSurfaceReversed(const int& rShapefunctionOrder, const int& rPolygonDiscretization)
	{
		Polygon boundaries(m_embedded_loops, rPolygonDiscretization);
		Polygon full_boundaries(m_trimming_loops, rPolygonDiscretization);
		Polygon difference = full_boundaries.GetDifference(boundaries);

		return GetIntegrationNodesSurface(difference, rShapefunctionOrder);
	}
	/* Compues the edge integration points which lay on the borders of the trimming edge. */
	std::vector<Node<3>::Pointer> BrepFace::GetIntegrationNodesTrimmingCurve(
		const int& rShapefunctionOrder,
		const int& rTrimIndex,
		const double& rAccuracy)
	{
		BrepTrimmingCurve trimming_curve = GetTrimmingCurve(rTrimIndex);

		int number_of_knots = (m_knot_vector_u.size() + m_knot_vector_v.size());
		std::vector<double> intersections = trimming_curve.GetKnotIntersections(m_p, m_q, m_knot_vector_u, m_knot_vector_v, number_of_knots*10);

		int degree = m_p + m_q + 1;

		std::vector<array_1d<double, 3>> quadrature_points = trimming_curve.GetQuadraturePoints(intersections, degree);

		std::vector<Node<3>::Pointer> NodeVectorElement = GetIntegrationNodes(quadrature_points, rShapefunctionOrder);
		for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
		{
			std::vector<Vector> derivatives;
			trimming_curve.EvaluateCurveDerivatives(derivatives, 2, quadrature_points[k][2]);
			NodeVectorElement[k]->SetValue(TANGENTS_BASIS_VECTOR, derivatives[1]);
		}
		return NodeVectorElement;
	}

	std::vector<Node<3>::Pointer> BrepFace::GetIntegrationNodesTrimmingCurveMaster(
		const std::vector<Point>& rPoints,
		const int& rShapefunctionOrder,
		const int& rTrimIndex,
		const int& rPQSlave,
		const double& rAccuracy, const double& rModelTolerance, const int& rMaxIterations)
	{
		std::vector<array_1d<double, 2>> slave_knot_intersection_points = GetClosestPointsTrimmingCurve(rPoints, rTrimIndex, rAccuracy, rMaxIterations);

		int number_of_knots = (m_knot_vector_u.size() + m_knot_vector_v.size());
		BrepTrimmingCurve trimming_curve = GetTrimmingCurve(rTrimIndex);
		std::vector<double> slave_knot_intersection_parameters = trimming_curve.GetClosestPoints(slave_knot_intersection_points, number_of_knots*10);

		std::vector<double> master_knot_intersection_parameters = trimming_curve.GetKnotIntersections(m_p, m_q, m_knot_vector_u, m_knot_vector_v, number_of_knots*10);
		for (unsigned int i = 0; i < master_knot_intersection_parameters.size(); i++)
		{
			slave_knot_intersection_parameters.push_back(master_knot_intersection_parameters[i]);
		}
		std::sort(slave_knot_intersection_parameters.begin(), slave_knot_intersection_parameters.end());
		std::vector<double> intersections;
		intersections.push_back(slave_knot_intersection_parameters[0]);
		for (unsigned int i = 1; i < slave_knot_intersection_parameters.size(); i++)
		{
			if (std::abs(slave_knot_intersection_parameters[i - 1] - slave_knot_intersection_parameters[i]) > rAccuracy)
			{
				intersections.push_back(slave_knot_intersection_parameters[i]);
			}
		}
		int highest_polynomial_degree = std::max(rPQSlave, (int) (m_p + m_q + 1));
		std::vector<array_1d<double, 4>> integration_points = trimming_curve.GetIntegrationPoints(intersections, highest_polynomial_degree);

		std::vector<Node<3>::Pointer> NodeVectorElement = GetIntegrationNodes(integration_points, rShapefunctionOrder);
		for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
		{
			std::vector<Vector> derivatives;
			trimming_curve.EvaluateCurveDerivatives(derivatives, 2, integration_points[k][2]);
			NodeVectorElement[k]->SetValue(TANGENTS_BASIS_VECTOR, derivatives[1]);
		}

		return NodeVectorElement;
	}

	Node<3>::Pointer BrepFace::GetIntegrationNodePoint(
		const int& rTrimIndex,
		const int& rShapefunctionOrder)
	{
		Node<3>::Pointer node = Node<3>::Pointer(new Node<3>(1, 0.0, 0.0, 0.0));
		for (int i = 0; i < m_embedded_points.size(); i++)
		{
			if (m_embedded_points[i].trim_index == rTrimIndex)
			{
				double u = m_embedded_points[i].local_coordinates[0];
				double v = m_embedded_points[i].local_coordinates[1];

				EvaluateIntegrationNode(u, v, rShapefunctionOrder, node);
			}
		}
		return node;
	}

	void BrepFace::EvaluateIntegrationNodesTrimmingCurveSlave(
		std::vector<Node<3>::Pointer>& rNodes,
		const int& rShapefunctionOrder,
		const int& rTrimIndex,
		const double& rAccuracy, const double& rModelTolerance, const int& rMaxIterations)
	{
		BrepTrimmingCurve trimming_curve = GetTrimmingCurve(rTrimIndex);
		std::vector<array_1d<double,2>> closest_points = GetClosestNodesTrimmingCurve(rNodes, rTrimIndex, rAccuracy, rModelTolerance, rMaxIterations);
		//int number_of_knots = (m_knot_vector_u.size() + m_knot_vector_v.size());
		//std::vector<double> closest_parameters = trimming_curve.GetClosestPoints(closest_points, number_of_knots * 25);
		for (unsigned int i = 0; i < rNodes.size(); i++)
		{
			double parameter = 0;

			//// obtain location of closest point in curve parameters
			double percentage = 0.5;
			if (rNodes[i]->Has(CURVE_PARAMETER_KNOT_LOCATION_PERCENTAGE))
				double percentage = rNodes[i]->GetValue(CURVE_PARAMETER_KNOT_LOCATION_PERCENTAGE);

			trimming_curve.GetClosestPoint(closest_points[i], percentage, parameter);


			Point point3d;
			EvaluateSurfacePoint(point3d, closest_points[i][0], closest_points[i][1]);

			std::vector<Vector> derivatives;
			trimming_curve.EvaluateCurveDerivatives(derivatives, 2, parameter);
			double error_distance = sqrt((point3d[0] - derivatives[0][0]) * (point3d[0] - derivatives[0][0])
				+ (point3d[1] - derivatives[0][1]) * (point3d[1] - derivatives[0][1]));
			//if (error_distance > rAccuracy*100)
			//	KRATOS_ERROR <<  "Bad projection: " << error_distance << std::endl;

			EvaluateIntegrationNodeSlave(closest_points[i][0], closest_points[i][1], rShapefunctionOrder, rNodes[i]);
			rNodes[i]->SetValue(TANGENTS_BASIS_VECTOR_SLAVE, derivatives[1]);
		}
	}



  //To be deleted
 // std::vector<Node<3>::Pointer> BrepFace::GetQuadraturePointsOfTrimmingCurveWithPoints(const int& shapefunction_order, const int& trim_index, std::vector<Point> intersection_points)
 // {
 //   BrepTrimmingCurve trimming_curve = GetTrimmingCurve(trim_index);
 //   //trimming_curve.PrintData();
 //
 //   std::vector<array_1d<double, 2>> intersection_points_2d;
 //   for (unsigned int i = 0; i < intersection_points.size(); i++)
 //   {
	//	std::cout << "Intersection point: X: " << intersection_points[i].X() << ", Y: " << intersection_points[i].Y() << ", Z: " << intersection_points[i].Z() << std::endl;
 //     //GetLocalParameterOfPointOnTrimmingCurve(intersection_points[i], trimming_curve, )
 //     double u = 0;
 //     double v = 0;
 //     bool success = ProjectionNewtonRaphson(intersection_points[i], u, v, 1e-7,20);
	//  std::cout << "After first Newton Raphson: u: " << u << ", v: " << v << std::endl;
 //     //GetClosestPoint(intersection_points[i], u, v);
 //     if (!success)
 //     {
 //       std::vector<array_1d<double, 2>> trimming_curve_loop = trimming_curve.CreatePolygon(40);
 //       double distance = 1e10;
 //       array_1d<double, 2> initial_guess;
 //       for (auto point_itr = trimming_curve_loop.begin(); point_itr != trimming_curve_loop.end(); ++point_itr)
 //       {
 //         Point point_global;
 //         EvaluateSurfacePoint(point_global, (*point_itr)[0], (*point_itr)[1]);
 //         double new_distance = sqrt((intersection_points[i][0] - point_global[0]) * (intersection_points[i][0] - point_global[0])
 //           + (intersection_points[i][1] - point_global[1]) * (intersection_points[i][1] - point_global[1])
 //           + (intersection_points[i][2] - point_global[2]) * (intersection_points[i][2] - point_global[2]));
 //         if (new_distance < distance)
 //         {
 //           distance = new_distance;
 //           initial_guess = (*point_itr);
 //         }
 //       }
 //       u = initial_guess[0];
 //       v = initial_guess[1];
	//	std::cout << "Initial guesses: u: " << u << ", v: " << v << std::endl;
 //       bool success = ProjectionNewtonRaphson(intersection_points[i], u, v, 1e-7,20);
 //       if (!success) {
 //         std::cout << "100 iteration points needed." << std::endl;
 //         trimming_curve_loop = trimming_curve.CreatePolygon(100);
 //         distance = 1e10;
 //         for (auto point_itr = trimming_curve_loop.begin(); point_itr != trimming_curve_loop.end(); ++point_itr)
 //         {
 //           Point point_global;
 //           EvaluateSurfacePoint(point_global, (*point_itr)[0], (*point_itr)[1]);
 //           double new_distance = sqrt((intersection_points[i][0] - point_global[0]) * (intersection_points[i][0] - point_global[0])
 //             + (intersection_points[i][1] - point_global[1]) * (intersection_points[i][1] - point_global[1])
 //             + (intersection_points[i][2] - point_global[2]) * (intersection_points[i][2] - point_global[2]));
 //           if (new_distance < distance)
 //           {
 //             distance = new_distance;
 //             initial_guess = (*point_itr);
 //           }
 //         }
 //         u = initial_guess[0];
 //         v = initial_guess[1];
 //         std::cout << "Initial: u=" << u << ", v=" << v << std::endl;
 //         bool success = ProjectionNewtonRaphson(intersection_points[i], u, v, 1e-7,20);
 //         if (!success)
 //         {
 //           u = initial_guess[0];
 //           v = initial_guess[1];
 //           //KRATOS_ERROR << "Point not found after 100 pts." << std::endl;
 //         }
 //       }
 //     }
 //     array_1d<double, 2> point;
	//  point[0] = u;
	//  point[1] = v;
 //     intersection_points_2d.push_back(point);
 //   }
 //
	//for (auto no = intersection_points_2d.begin(); no != intersection_points_2d.end(); ++no)
	//{/*
	// NodeVectorElement2.push_back(no->get());*/
	//	KRATOS_WATCH(*no)
	//}
	//
 //   std::vector<double> intersections_slave = trimming_curve.GetClosestPoints(intersection_points_2d, 100);
 //   std::cout << "intersections slave: ";
 //   for (unsigned int i = 0; i < intersections_slave.size(); i++)
 //   {
 //     std::cout << intersections_slave[i] << ", ";
 //   }
 //   std::cout << std::endl;
 //
 //   std::vector<double> intersections_master = trimming_curve.GetKnotIntersections(m_p, m_q, m_knot_vector_u, m_knot_vector_v, 200);
 //   std::cout << "intersections master: ";
 //   for (unsigned int i = 0; i < intersections_master.size(); i++)
 //   {
 //     std::cout << intersections_master[i] << ", ";
 //   }
 //   std::cout << std::endl;
 //   for (unsigned int i = 0; i < intersections_master.size(); i++)
 //   {
 //     intersections_slave.push_back(intersections_master[i]);
 //   }
 //   std::sort(intersections_slave.begin(), intersections_slave.end());
 //   std::vector<double> intersections;
 //   intersections.push_back(intersections_slave[0]);
 //   for (unsigned int i = 1; i < intersections_slave.size(); i++)
 //   {
 //     if (intersections_slave[i - 1] != intersections_slave[i])
 //       intersections.push_back(intersections_slave[i]);
 //   }
 //
 //   int highest_polynomial_order = std::max(m_p, m_q);
 //   std::vector<array_1d<double, 3>> quadrature_points = trimming_curve.GetQuadraturePoints(intersections, highest_polynomial_order);
 //
 //   std::vector<Node<3>::Pointer> NodeVectorElement = GetIntegrationNodes(quadrature_points, shapefunction_order);
 //   for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
 //   {
	//	std::vector<Vector> derivatives;
	//	trimming_curve.EvaluateCurveDerivatives(derivatives, 2, quadrature_points[k][2]);
 //     //array_1d<double, 2> tangents_basis = trimming_curve.GetBaseVector(quadrature_points[k][2]);
 //     Vector tangents_basis_vector(2);
 //     tangents_basis_vector[0] = derivatives[1][0];
 //     tangents_basis_vector[1] = derivatives[1][1];
 //     NodeVectorElement[k]->SetValue(TANGENTS_BASIS_VECTOR, tangents_basis_vector);
 //     //KRATOS_WATCH(tangents_basis_vector)
 //   }
 //
	////std::vector<Node<3>::Pointer> NodeVectorElement2;
	////std::cout << "here 1.1: " << NodeVectorElement .size() << std::endl;
	////for (auto no = NodeVectorElement.begin(); no != NodeVectorElement.end(); ++no)
	////{/*
	////	NodeVectorElement2.push_back(no->get());*/
	////	KRATOS_WATCH(no->get()->X())
	////		KRATOS_WATCH(no->get()->Y())
	////		KRATOS_WATCH(no->get()->Z())
	////}
	////for (unsigned int i = 0; i < NodeVectorElement.size(); i++)
	////{
	////	nodevector.push_back(NodeVectorElement[i]);
	////}
	////for (auto no = nodevector.begin(); no != nodevector.end(); ++no)
	////{/*
	//// NodeVectorElement2.push_back(no->get());*/
	////	KRATOS_WATCH(no->get()->X())
	////	KRATOS_WATCH(no->get()->Y())
	////	KRATOS_WATCH(no->get()->Z())
	////}
 //   return NodeVectorElement;
 // }

	/* Computes all intersections of knots with specific trimming curves.
	*  @param trim_index index of trimming curve
	*/
	std::vector<Point> BrepFace::GetKnotIntersections(const int& rTrimIndex)
	{
		BrepTrimmingCurve trimming_curve = GetTrimmingCurve(rTrimIndex);
		std::vector<double> intersections = trimming_curve.GetKnotIntersections(m_p, m_q, m_knot_vector_u, m_knot_vector_v, 200);
		std::vector<Point> points;
		for (unsigned int i = 0; i < intersections.size(); i++)
		{
			Point point_parameter, point_global;
			trimming_curve.EvaluateCurvePoint(point_parameter, intersections[i]);

			EvaluateSurfacePoint(point_global, point_parameter[0], point_parameter[1]);
			points.push_back(point_global);
		}
		return points;
	}

	/* Constructs Kratos nodes out of Points. The nodes are enhanced with rShapefunctionOrder of
	*  shape functions, the global location, the integration weight, the id of the patch and the
	*  local parameter location on the patch. */
	std::vector<Node<3>::Pointer> BrepFace::GetIntegrationNodes(
		std::vector<array_1d<double, 3>>& rPoints, const int& rShapefunctionOrder)
	{
		std::vector<Node<3>::Pointer> NodeVector;
		for (unsigned int i = 0; i < rPoints.size(); i++)
		{
			Node<3>::Pointer node = Node<3>::Pointer( new Node<3>(1, 0.0, 0.0, 0.0) ); // = Kratos::make_shared< Node<3> >(1,0.0,0.0,0.0);
			EvaluateIntegrationNode(rPoints[i][0], rPoints[i][1], rShapefunctionOrder, node);
			node->SetValue(INTEGRATION_WEIGHT, rPoints[i][2]);
			NodeVector.push_back(node);
		}
		return NodeVector;
	}

	/* Constructs Kratos nodes out of Points. The nodes are enhanced with rShapefunctionOrder of
	*  shape functions, the global location, the integration weight, the id of the patch and the
	*  local parameter location on the patch. */
	std::vector<Node<3>::Pointer> BrepFace::GetIntegrationNodes(
		std::vector<array_1d<double, 4>>& rPoints, const int& rShapefunctionOrder)
	{
		std::vector<Node<3>::Pointer> NodeVector;
		for (unsigned int i = 0; i < rPoints.size(); i++)
		{
			Node<3>::Pointer node = Node<3>::Pointer(new Node<3>(1, 0.0, 0.0, 0.0)); // = Kratos::make_shared< Node<3> >(1,0.0,0.0,0.0);
			EvaluateIntegrationNode(rPoints[i][0], rPoints[i][1], rShapefunctionOrder, node);
			node->SetValue(INTEGRATION_WEIGHT, rPoints[i][2]);
			node->SetValue(CURVE_PARAMETER_KNOT_LOCATION_PERCENTAGE, rPoints[i][3]);
			NodeVector.push_back(node);
		}
		return NodeVector;
	}


    void BrepFace::EvaluateIntegrationNode(
        const double& rU, const double& rV,
        const int& rShapefunctionOrder,
        Node<3>::Pointer pNode)
    {
        Vector new_point = ZeroVector(3);

        int span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, rU);
        int span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, rV);

        Vector N = ZeroVector((m_q + 1)*(m_p + 1));

        Vector ControlPointIDs = ZeroVector((m_q + 1)*(m_p + 1));
        Vector local_parameter(2);
        local_parameter[0] = rU;
        local_parameter[1] = rV;

        Matrix R;
        EvaluateNURBSFunctions(span_u, span_v, rU, rV, R);
        int k = 0;
        for (int c = 0; c <= m_q; c++)
        {
            for (int b = 0; b <= m_p; b++)
            {
                // the control point vector is filled up by first going over u, then over v
                int ui = span_u - m_p + b;
                int vi = span_v - m_q + c;
                int m_n_u = m_knot_vector_u.size() - m_p - 1;
                int control_point_index = vi * m_n_u + ui;

                if (rShapefunctionOrder > -1)
                    N(k) = R(b, c);

                ControlPointIDs(k) = m_control_points_ids[control_point_index];

                new_point[0] += R(b, c) * mp_model_part->pGetNode(m_control_points_ids[control_point_index])->X();
                new_point[1] += R(b, c) * mp_model_part->pGetNode(m_control_points_ids[control_point_index])->Y();
                new_point[2] += R(b, c) * mp_model_part->pGetNode(m_control_points_ids[control_point_index])->Z();

                k++;
            }
        }
        pNode->X() = new_point(0);
        pNode->Y() = new_point(1);
        pNode->Z() = new_point(2);

        if (rShapefunctionOrder > -1)
            pNode->SetValue(NURBS_SHAPE_FUNCTIONS, N);

        if (rShapefunctionOrder > 0)
        {
            Matrix DN_De;
            Matrix DDN_DDe;
            EvaluateNURBSFunctionsDerivatives(span_u, span_v, rU, rV, DN_De, DDN_DDe);
            pNode->SetValue(NURBS_SHAPE_FUNCTION_DERIVATIVES, DN_De);

            if (rShapefunctionOrder > 1)
                pNode->SetValue(NURBS_SHAPE_FUNCTION_SECOND_DERIVATIVES, DDN_DDe);
        }
        pNode->SetValue(FACE_BREP_ID, this->Id());
        pNode->SetValue(LOCAL_PARAMETERS, local_parameter);
        pNode->SetValue(CONTROL_POINT_IDS, ControlPointIDs);
    }

	void BrepFace::EvaluateIntegrationNodeSlave(const double& rU, const double& rV,
		const int& rShapefunctionOrder,
		Node<3>::Pointer pNode)
	{
		int span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, rU);
		int span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, rV);

		Vector new_point = ZeroVector(3);

		Vector N = ZeroVector((m_q + 1)*(m_p + 1));

		Vector ControlPointIDs = ZeroVector((m_q + 1)*(m_p + 1));
		Vector local_parameter(2);
		local_parameter[0] = rU;
		local_parameter[1] = rV;

		Matrix R;
		EvaluateNURBSFunctions(span_u, span_v, rU, rV, R);

		int k = 0;
		for (int c = 0; c <= m_q; c++)
		{
			for (int b = 0; b <= m_p; b++)
			{
				// the control point vector is filled up by first going over u, then over v
				int ui = span_u - m_p + b;
				int vi = span_v - m_q + c;
				int m_n_u = m_knot_vector_u.size() - m_p - 1;
				int control_point_index = vi*m_n_u + ui;

				if (rShapefunctionOrder > -1)
					N(k) = R(b, c);

				ControlPointIDs(k) = m_control_points_ids[control_point_index];

				new_point[0] += R(b, c) * mp_model_part->pGetNode(m_control_points_ids[control_point_index])->X();
				new_point[1] += R(b, c) * mp_model_part->pGetNode(m_control_points_ids[control_point_index])->Y();
				new_point[2] += R(b, c) * mp_model_part->pGetNode(m_control_points_ids[control_point_index])->Z();

				k++;
			}
		}
		pNode->SetValue(LOCATION_SLAVE, new_point);

		double error_distance = sqrt(pow(pNode->X() - new_point[0], 2) + pow(pNode->Y() - new_point[1], 2) + pow(pNode->Z() - new_point[2], 2));
		if (error_distance > 0.1)
		{
			std::cout << "Error in master-slave projection: " << error_distance << std::endl;
			std::cout << "Punkt X: " << pNode->X() << std::endl;
			std::cout << "Punkt Y: " << pNode->Y() << std::endl;
			std::cout << "Punkt Z: " << pNode->Z() << std::endl;
            KRATOS_WATCH(new_point);
		}

		if (rShapefunctionOrder > -1)
			pNode->SetValue(NURBS_SHAPE_FUNCTIONS_SLAVE, N);

		if (rShapefunctionOrder > 0)
		{
			Matrix DN_De;
			Matrix DDN_DDe;
			EvaluateNURBSFunctionsDerivatives(span_u, span_v, rU, rV, DN_De, DDN_DDe);
			pNode->SetValue(NURBS_SHAPE_FUNCTION_DERIVATIVES_SLAVE, DN_De);

			if (rShapefunctionOrder > 1)
				pNode->SetValue(NURBS_SHAPE_FUNCTION_SECOND_DERIVATIVES_SLAVE, DDN_DDe);
		}
		pNode->SetValue(FACE_BREP_ID_SLAVE, this->Id());
		pNode->SetValue(LOCAL_PARAMETERS_SLAVE, local_parameter);
		pNode->SetValue(CONTROL_POINT_IDS_SLAVE, ControlPointIDs);
	}

	/* Uses Newton Raphson scheme to project a certain point (rPoint) towards its
	*  closest point on the surface. u and v are used as initial guess.
	*/
	bool BrepFace::ProjectionNewtonRaphson(const Point& rPoint, double& u, double& v,
		const double& rAccuracy, const int& rMaxIterations)
	{
		double norm_delta_u = 1e10;
        //KRATOS_WATCH(rMaxIterations)
        for (int i = 0; i < rMaxIterations; ++i)
        {
            // newton_raphson_point is evaluated
            Point newton_raphson_point;
            EvaluateSurfacePoint(newton_raphson_point, u, v);

            Vector difference = ZeroVector(3);
            // Distance between current Q_k and P
            // The distance between Q (on the CAD surface) and P (on the FE-mesh) is evaluated
            difference(0) = newton_raphson_point[0] - rPoint[0];
            difference(1) = newton_raphson_point[1] - rPoint[1];
            difference(2) = newton_raphson_point[2] - rPoint[2];

            Matrix hessian = ZeroMatrix(2, 2);
            Vector gradient = ZeroVector(2);
            // The distance is used to compute Hessian and gradient
            EvaluateGradientsForClosestPointSearch(difference, hessian, gradient, u, v);

            //std::cout << "4" << std::endl;
            double det_H = 0;
            Matrix inv_H = ZeroMatrix(2, 2);

            // u and v are updated
            MathUtils<double>::InvertMatrix(hessian, inv_H, det_H);
            Vector delta_u = prod(inv_H, gradient);
            u -= delta_u(0);
            v -= delta_u(1);
            //std::cout << "u: " << u << ", v: " << v << std::endl;
            //KRATOS_WATCH(delta_u)

            //if (u >= (std::max(m_knot_vector_u[m_knot_vector_u.size() - 1], m_knot_vector_u[0])+1e-6))
            //{
            //    //std::cout << "u: " << u
            //    //    << ", m_knot_vector_u[m_knot_vector_u.size() - 1]: " << m_knot_vector_u[m_knot_vector_u.size() - 1]
            //    //    << ", m_knot_vector_v[0]: " << m_knot_vector_u[0] << std::endl;
            //    u = std::max(m_knot_vector_u[m_knot_vector_u.size() - 1], m_knot_vector_u[0]) - 1e-7;
            //    //std::cout << "out of borders" << std::endl;
            //    //std::cout << "u: " << u << std::endl;
            //    //return false;
            //}
            //if (u <= (std::min(m_knot_vector_u[m_knot_vector_u.size() - 1], m_knot_vector_u[0]) - 1e-6))
            //{
            //    //std::cout << "u: " << u
            //    //    << ", m_knot_vector_u[m_knot_vector_u.size() - 1]: " << m_knot_vector_u[m_knot_vector_u.size() - 1]
            //    //    << ", m_knot_vector_u[0]: " << m_knot_vector_u[0] << std::endl;
            //    u = std::min(m_knot_vector_u[m_knot_vector_u.size() - 1], m_knot_vector_u[0]) + 1e-7;
            //    //std::cout << "out of borders" << std::endl;
            //    //std::cout << "u: " << u << std::endl;
            //    //return false;
            //}

            //if (v >= (std::max(m_knot_vector_v[m_knot_vector_v.size() - 1], m_knot_vector_v[0]) + 1e-6))
            //{
            //    v = std::max(m_knot_vector_v[m_knot_vector_v.size() - 1], m_knot_vector_v[0]) - 1e-7;
            //    //std::cout << "v: " << v
            //    //    << ", m_knot_vector_v[m_knot_vector_v.size() - 1]: " << m_knot_vector_v[m_knot_vector_v.size() - 1]
            //    //    << ", m_knot_vector_v[0]: " << m_knot_vector_v[0] << std::endl;
            //    //std::cout << "out of borders" << std::endl;
            //    //std::cout << "v: " << v << std::endl;
            //    //return false;
            //}
            //if (v <= (std::min(m_knot_vector_v[m_knot_vector_v.size() - 1], m_knot_vector_v[0]) - 1e-6))
            //{
            //    v = std::min(m_knot_vector_v[m_knot_vector_v.size() - 1], m_knot_vector_v[0]) + 1e-7;
            //    //std::cout << "v: " << v 
            //    //    << ", m_knot_vector_v[m_knot_vector_v.size() - 1]: " << m_knot_vector_v[m_knot_vector_v.size() - 1] 
            //    //    << ", m_knot_vector_v[0]: " << m_knot_vector_v[0] << std::endl;
            //    //std::cout << "out of borders" << std::endl;
            //    //std::cout << "v: " << v << std::endl;
            //    //return false;
            //}

            EvaluateSurfacePoint(newton_raphson_point, u, v);
            //KRATOS_WATCH(newton_raphson_point)
            difference(0) = newton_raphson_point[0] - rPoint[0];
            difference(1) = newton_raphson_point[1] - rPoint[1];
            difference(2) = newton_raphson_point[2] - rPoint[2];

            //std::cout << "5" << std::endl;
            norm_delta_u = norm_2(difference);
            
            if (norm_2(delta_u) < rAccuracy)
            {
                return true;
                //std::cout << "projection completed with orthogonal acceptance" << std::endl;
            }
            if (norm_delta_u < rAccuracy)
            {
                return true;
                //std::cout << "projection completed" << std::endl;
            }
        }
        std::cout << "projection not completed with norm delta u: " << norm_delta_u << std::endl;
		return false;
	}
	/* Uses Newton-Raphson projection to project point to closest point on surface.
	*  ...for global naming issues */
	bool BrepFace::GetClosestPoint(const Point& rPoint, double& u, double& v,
		const double& rAccuracy, const int& rMaxIterations)
	{
		return ProjectionNewtonRaphson(rPoint, u, v,
			rAccuracy, rMaxIterations);
	}

	std::vector<array_1d<double, 2>> BrepFace::GetClosestPointsTrimmingCurve(const std::vector<Point>& rPoints, const int& rTrimIndex,
		const double& rAccuracy, const int& rMaxIterations)
	{
		BrepTrimmingCurve trimming_curve = GetTrimmingCurve(rTrimIndex);
		std::vector<array_1d<double, 2>> trimming_curve_polygon;

		int number_of_knots = (m_knot_vector_u.size() + m_knot_vector_v.size());

		std::vector<array_1d<double, 2>> closest_points;
		for (unsigned int i = 0; i < rPoints.size(); i++)
		{

			double u = 0;
			double v = 0;
			bool success = ProjectionNewtonRaphson(rPoints[i], u, v, rAccuracy, rMaxIterations);
			if (!success)
			{
				if (trimming_curve_polygon.size()<(number_of_knots * 3))
					trimming_curve_polygon = trimming_curve.CreatePolygon(number_of_knots * 3);
				double distance = 1e10;
				array_1d<double, 2> initial_guess;

				for (auto point_itr = trimming_curve_polygon.begin(); point_itr != trimming_curve_polygon.end(); ++point_itr)
				{
					Point point_global;
					EvaluateSurfacePoint(point_global, (*point_itr)[0], (*point_itr)[1]);
					double new_distance = sqrt(
						pow((rPoints[i][0] - point_global[0]),2)
						+ pow((rPoints[i][1] - point_global[1]),2)
						+ pow((rPoints[i][2] - point_global[2]),2));

					if (new_distance < distance)
					{
						distance = new_distance;
						initial_guess = (*point_itr);
					}
				}
				u = initial_guess[0];
				v = initial_guess[1];
				bool success = ProjectionNewtonRaphson(rPoints[i], u, v, rAccuracy, rMaxIterations);
				if (!success) {
					if (trimming_curve_polygon.size()<(number_of_knots * 25))
						trimming_curve_polygon = trimming_curve.CreatePolygon(number_of_knots * 25);
					distance = 1e10;
					for (auto point_itr = trimming_curve_polygon.begin(); point_itr != trimming_curve_polygon.end(); ++point_itr)
					{
						Point point_global;
						EvaluateSurfacePoint(point_global, (*point_itr)[0], (*point_itr)[1]);
						double new_distance = sqrt((rPoints[i][0] - point_global[0]) * (rPoints[i][0] - point_global[0])
							+ (rPoints[i][1] - point_global[1]) * (rPoints[i][1] - point_global[1])
							+ (rPoints[i][2] - point_global[2]) * (rPoints[i][2] - point_global[2]));
						if (new_distance < distance)
						{
							distance = new_distance;
							initial_guess = (*point_itr);
						}
					}
					u = initial_guess[0];
					v = initial_guess[1];
					bool success = ProjectionNewtonRaphson(rPoints[i], u, v, rAccuracy, rMaxIterations);

					if (!success)
					{
						//u = initial_guess[0];
						//v = initial_guess[1];

						std::cout << "Intersection point: X: " << rPoints[i][0]
							<< ", Y: " << rPoints[i][1] << ", Z: " << rPoints[i][2] << std::endl;
						std::cout << "WARNING: Point not found! With " << number_of_knots * 25 << " needed iteration points." << std::endl;
					}
				}
			}
			array_1d<double, 2> point;
			point[0] = u;
			point[1] = v;
			closest_points.push_back(point);
		}
		return closest_points;
	}

	std::vector<array_1d<double, 2>> BrepFace::GetClosestNodesTrimmingCurve(
		const std::vector<Node<3>::Pointer>& rNodes,
		const int& rTrimIndex,
		const double& rAccuracy, const double& rModelTolerance, const int& rMaxIterations)
	{
		BrepTrimmingCurve trimming_curve = GetTrimmingCurve(rTrimIndex);
		std::vector<array_1d<double, 2>> trimming_curve_polygon;

		int number_of_knots = (m_knot_vector_u.size() + m_knot_vector_v.size());

		std::vector<array_1d<double, 2>> closest_points;
		for (unsigned int i = 0; i < rNodes.size(); i++)
		{
			Point point(rNodes[i]->X(), rNodes[i]->Y(), rNodes[i]->Z());
			double u = 0;
			double v = 0;
			bool success = ProjectionNewtonRaphson(point, u, v, rAccuracy, rMaxIterations);

			if (!success)
			{
				if (trimming_curve_polygon.size()<(number_of_knots * 3))
					trimming_curve_polygon = trimming_curve.CreatePolygon(number_of_knots * 3);

				GetClosestPolygonPointParameters(trimming_curve_polygon, point, u, v);

				bool success = ProjectionNewtonRaphson(point, u, v, rAccuracy, rMaxIterations);

				if (!success) {
					if (trimming_curve_polygon.size()<(number_of_knots * 25))
						trimming_curve_polygon = trimming_curve.CreatePolygon(number_of_knots * 25);

					GetClosestPolygonPointParameters(trimming_curve_polygon, point, u, v);

					bool success = ProjectionNewtonRaphson(point, u, v, rAccuracy, rMaxIterations);

					if (!success)
					{
						Point difference_point(0, 0, 0);
						EvaluateSurfacePoint(difference_point, u, v);

						double error_distance = sqrt(pow(difference_point[0] - point[0], 2)
							+ pow(difference_point[1] - point[1], 2) + pow(difference_point[2] - point[2], 2));

						if (error_distance < rModelTolerance)
						{
							success = true;
						}
						else
						{
							std::cout << "WARNING: Point not found! With " << number_of_knots * 25 << " iteration points." << std::endl;
							std::cout << "Error distance: " << error_distance << std::endl;
							std::cout << "Closest point: X: " << point[0] << ", Y: " << point[1] << ", Z: " << point[2] << std::endl;
						}
					}
				}
			}
			array_1d<double, 2> closest_point;
			closest_point[0] = u;
			closest_point[1] = v;
			closest_points.push_back(closest_point);
		}

		return closest_points;
	}

	void BrepFace::GetClosestPolygonPointParameters(const std::vector<array_1d<double, 2>>& rPolygon, const Point& rPoint, double& rU, double& rV)
	{
		array_1d<double, 2> initial_guess;
		double distance = 1e10;
		for (auto point_itr = rPolygon.begin(); point_itr != rPolygon.end(); ++point_itr)
		{
			Point point_global;
			EvaluateSurfacePoint(point_global, (*point_itr)[0], (*point_itr)[1]);
			double new_distance = sqrt(
				pow((rPoint[0] - point_global[0]), 2)
				+ pow((rPoint[1] - point_global[1]), 2)
				+ pow((rPoint[2] - point_global[2]), 2));
			if (new_distance < distance)
			{
				distance = new_distance;
				initial_guess = (*point_itr);
			}
		}
		rU = initial_guess[0];
		rV = initial_guess[1];
	}

  //void BrepFace::GetClosestPoint(const Point& point, double& u, double& v)
  //{
  //  double norm_delta_u = 100000000;
  //  unsigned int k = 0;
  //  unsigned int max_itr = 40;
  //
  //  while (norm_delta_u > 1e-5)
  //  {
  //    // newton_raphson_point is evaluated
  //    Point newton_raphson_point;
  //    EvaluateSurfacePoint(newton_raphson_point, u, v);
  //
  //    Vector difference = ZeroVector(3); // Distance between current Q_k and P
  //    // The distance between Q (on the CAD surface) and P (on the FE-mesh) is evaluated
  //    difference(0) = newton_raphson_point[0] - point[0];
  //    difference(1) = newton_raphson_point[1] - point[1];
  //    difference(2) = newton_raphson_point[2] - point[2];
  //
  //    Matrix hessian = ZeroMatrix(2, 2);
  //    Vector gradient = ZeroVector(2);
  //    // The distance is used to compute Hessian and gradient
  //    EvaluateGradientsForClosestPointSearch(difference, hessian, gradient, u, v);
  //
  //    double det_H = 0;
  //    Matrix inv_H = ZeroMatrix(2, 2);
  //
  //    // u and v are updated
  //    MathUtils<double>::InvertMatrix(hessian, inv_H, det_H);
  //    Vector delta_u = prod(inv_H, gradient);
  //    u -= delta_u(0);
  //    v -= delta_u(1);
  //
  //    norm_delta_u = norm_2(delta_u);
  //
  //    k++;
  //    if (k>max_itr)
  //      KRATOS_THROW_ERROR(std::runtime_error, "Newton-Raphson to find closest point did not converge in the following number of iterations: ", k - 1);
  //  }
  //}

	/* Obtains the closest surface point including the specified number of
	*  shape function and derivatives, the local parameter location, the
	*/
    bool BrepFace::GetClosestIntegrationNode(
        Node<3>::Pointer& rClosestNode,
        const Node<3>::Pointer& rSpaceNode,
        const int& rShapefunctionOrder,
        const double& rAccuracy, const int& rMaxIterations)
    {
        //std::cout << "shit debugging" << std::endl;
        Point point(rSpaceNode->X(), rSpaceNode->Y(), rSpaceNode->Z());
        //std::cout << rSpaceNode->X() << "  " << rSpaceNode->Y() << "  " << rSpaceNode->Z() << std::endl;
        //std::cout << "shit debugging 1" << std::endl;

        Vector local_parameter = rClosestNode->GetValue(LOCAL_PARAMETERS);
        //std::cout << "shit debugging 2" << std::endl;
        //KRATOS_WATCH(local_parameter)
        double u = local_parameter(0);
        double v = local_parameter(1);
        //KRATOS_WATCH(u)
            //KRATOS_WATCH(v)
        double success = ProjectionNewtonRaphson(point, u, v, rAccuracy, rMaxIterations);

        if (success == false)
        {
            std::cout << "no success" << std::endl;
            return success;
        }
        //std::cout << "shit debugging 3" << std::endl;
    //KRATOS_WATCH(u)
    //KRATOS_WATCH(v)
        EvaluateIntegrationNode(u, v, rShapefunctionOrder, rClosestNode);
        return success;
        //KRATOS_WATCH(rClosestNode->GetValue(CONTROL_POINT_IDS))
    }


    /* Obtains the closest surface point including the specified number of
    *  shape function and derivatives, the local parameter location, the
    */
    bool BrepFace::GetIntegrationNodeUpdated(
        Node<3>::Pointer& rClosestNode,
        const Node<3>::Pointer& rSpaceNode,
        const int& rShapefunctionOrder,
        const double& rAccuracy, const int& rMaxIterations)
    {
        //std::cout << "shit debugging" << std::endl;
        Point point(rSpaceNode->X(), rSpaceNode->Y(), rSpaceNode->Z());
        //std::cout << rSpaceNode->X() << "  " << rSpaceNode->Y() << "  " << rSpaceNode->Z() << std::endl;
        //std::cout << "shit debugging 1" << std::endl;

        Vector local_parameter = rClosestNode->GetValue(LOCAL_PARAMETERS);
        //std::cout << "shit debugging 2" << std::endl;
        //KRATOS_WATCH(local_parameter)
        double u = local_parameter(0);
        double v = local_parameter(1);
        //KRATOS_WATCH(u)
        //KRATOS_WATCH(v)
        //double success = ProjectionNewtonRaphson(point, u, v, rAccuracy, rMaxIterations);

        //if (success == false)
        //{
        //    std::cout << "no success" << std::endl;
        //    return success;
        //}
        //std::cout << "shit debugging 3" << std::endl;
        //KRATOS_WATCH(u)
        //KRATOS_WATCH(v)
        EvaluateIntegrationNode(u, v, rShapefunctionOrder, rClosestNode);
        return true;
        //KRATOS_WATCH(rClosestNode->GetValue(CONTROL_POINT_IDS))
    }

  //// --------------------------------------------------------------------------
  //void BrepFace::MapNodeNewtonRaphson(const Node<3>::Pointer& node, Node<3>::Pointer& node_on_geometry)
  //{
  //  std::cout << "test hier" << std::endl;
  //  // Initialize P: point on the mesh
  //  Vector P = ZeroVector(3);
  //  P(0) = node->X();
  //  P(1) = node->Y();
  //  P(2) = node->Z();
  //  // Initialize Q_k: point on the CAD surface
  //  Vector Q_k = ZeroVector(3);
  //  Q_k(0) = node_on_geometry->X();
  //  Q_k(1) = node_on_geometry->Y();
  //  Q_k(2) = node_on_geometry->Z();
  //  // Initialize what's needed in the Newton-Raphson iteration
  //  Vector Q_minus_P = ZeroVector(3); // Distance between current Q_k and P
  //  Matrix myHessian = ZeroMatrix(2, 2);
  //  Vector myGradient = ZeroVector(2);
  //  double det_H = 0;
  //  Matrix InvH = ZeroMatrix(2, 2);
  //  Vector local_parameter = node_on_geometry->GetValue(LOCAL_PARAMETERS);
  //  double u_k = local_parameter(0);
  //  double v_k = local_parameter(1);
  //  //Node<3>::Pointer newtonRaphsonPoint;
  //
  //  double norm_delta_u = 100000000;
  //  unsigned int k = 0;
  //  unsigned int max_itr = 20;
  //  while (norm_delta_u > 1e-8)
  //  {
  //    // The distance between Q (on the CAD surface) and P (on the FE-mesh) is evaluated
  //    Q_minus_P(0) = Q_k(0) - P(0);
  //    Q_minus_P(1) = Q_k(1) - P(1);
  //    Q_minus_P(2) = Q_k(2) - P(2);
  //
  //    // The distance is used to compute Hessian and gradient
  //    EvaluateGradientsForClosestPointSearch(Q_minus_P, myHessian, myGradient, u_k, v_k);
  //
  //    // u_k and v_k are updated
  //    MathUtils<double>::InvertMatrix(myHessian, InvH, det_H);
  //    Vector delta_u = prod(InvH, myGradient);
  //    u_k -= delta_u(0);
  //    v_k -= delta_u(1);
  //
  //    // Q is updated
  //    Point point;
  //    EvaluateSurfacePoint(point, u_k, v_k);
  //    Q_k(0) = point[0];
  //    Q_k(1) = point[1];
  //    Q_k(2) = point[2];
  //
  //    KRATOS_WATCH(Q_k)
  //
  //      //Q_k(0) = newtonRaphsonPoint[0];
  //      //Q_k(1) = newtonRaphsonPoint[1];
  //      //Q_k(2) = newtonRaphsonPoint[2];
  //      norm_delta_u = norm_2(delta_u);
  //
  //    k++;
  //
  //    if (k>max_itr)
  //      KRATOS_THROW_ERROR(std::runtime_error, "Newton-Raphson to find closest point did not converge in the following number of iterations: ", k - 1);
  //  }
  //  node_on_geometry->X() = Q_k(0);
  //  node_on_geometry->Y() = Q_k(1);
  //  node_on_geometry->Z() = Q_k(2);
  //}
  //

	/* Check if point lays inside of the trimmed domain.
	*/
	bool BrepFace::CheckIfPointIsInside(const Vector& rNodeParameters)
	{
		Polygon polygon(m_trimming_loops, 10);
		return polygon.IsInside(rNodeParameters[0], rNodeParameters[1]);
	}



	//REFINEMENT OPERATIONS
	/* applies all geometry refinement operations
	*  @param[in] rGeometryRefinementParameters
	*/
	void BrepFace::ApplyGeometryRefinement(
		const GeometryRefinementParameters& rGeometryRefinementParameters)
	{
		/// Degree elevation
		// newDegree = Max(minP or P+deltaP)
		int this_order_elevation_p = rGeometryRefinementParameters.order_elevation_p;
		if (rGeometryRefinementParameters.min_order_p - m_p > 0)
			this_order_elevation_p = std::max((int)(rGeometryRefinementParameters.min_order_p - m_p),
				this_order_elevation_p);

		int this_order_elevation_q = rGeometryRefinementParameters.order_elevation_q;
		if (rGeometryRefinementParameters.min_order_q - m_q > 0)
			this_order_elevation_q = std::max((int)(rGeometryRefinementParameters.min_order_q - m_q),
				this_order_elevation_q);

		//std::cout << "check here degree elevation this_order_elevation_p: " << this_order_elevation_p << ", this_order_elevation_q: " << this_order_elevation_q << std::endl;
		if (this_order_elevation_p > 0 ||
			this_order_elevation_q > 0)
			DegreeElevate(this_order_elevation_p, this_order_elevation_q);

		/// Knot insertions
		//check whether increasing/decreasing

		// 1. insert defined knots
		Vector this_knot_insertions_u = ZeroVector(0);
		if (rGeometryRefinementParameters.knot_insertions_u.size() > 0)
			this_knot_insertions_u = rGeometryRefinementParameters.knot_insertions_u;

		Vector this_knot_insertions_v = ZeroVector(0);
		if (rGeometryRefinementParameters.knot_insertions_v.size() > 0)
			this_knot_insertions_v = rGeometryRefinementParameters.knot_insertions_v;

		if (this_knot_insertions_u.size() > 0 ||
			this_knot_insertions_v.size() > 0)
			RefineKnotVector(this_knot_insertions_u, this_knot_insertions_v);

		// 2. insert knots if higher number of elements is defined
		this_knot_insertions_u = ZeroVector(0);
		this_knot_insertions_v = ZeroVector(0);
		if (rGeometryRefinementParameters.multiply_knots_u > 1)
		{
			for (size_t i = 0; i<m_knot_vector_u.size() - 1; i++)
			{
				if (m_knot_vector_u[i] != m_knot_vector_u[i + 1])
				{
					for (int j = 0; j < rGeometryRefinementParameters.multiply_knots_u - 1; j++)
					{
						this_knot_insertions_u.resize(this_knot_insertions_u.size() + 1, true);
						this_knot_insertions_u[this_knot_insertions_u.size() - 1] = ((1.0 + j) / rGeometryRefinementParameters.multiply_knots_u*(m_knot_vector_u[i + 1] - m_knot_vector_u[i]) + m_knot_vector_u[i]);
					}
				}
			}
			//double new_knot_size = (m_knot_vector_u[m_knot_vector_u.size() - 1] - m_knot_vector_u[0]) / rGeometryRefinementParameters.multiply_knots_u;
			//double knot_insertion = m_knot_vector_u[0] + new_knot_size;
			//for (int i = 0; i<rGeometryRefinementParameters.multiply_knots_u - 1; i++)
			//{
			//	this_knot_insertions_u.resize(this_knot_insertions_u.size() + 1);
			//	this_knot_insertions_u[this_knot_insertions_u.size() - 1] = knot_insertion;
			//	knot_insertion += new_knot_size;
			//}
		}
		if (rGeometryRefinementParameters.multiply_knots_v > 1)
		{
			for (size_t i = 0; i<m_knot_vector_v.size() - 1; i++)
			{
				if (m_knot_vector_v[i] != m_knot_vector_v[i + 1])
				{
					for (int j = 0; j < rGeometryRefinementParameters.multiply_knots_v - 1; j++)
					{
						this_knot_insertions_v.resize(this_knot_insertions_v.size() + 1, true);
						this_knot_insertions_v[this_knot_insertions_v.size() - 1] = ((1.0 + j) / rGeometryRefinementParameters.multiply_knots_v*(m_knot_vector_v[i + 1] - m_knot_vector_v[i]) + m_knot_vector_v[i]);
					}
				}
			}
			//double new_knot_size = (m_knot_vector_u[m_knot_vector_u.size() - 1] - m_knot_vector_u[0]) / rGeometryRefinementParameters.multiply_knots_u;
			//double knot_insertion = m_knot_vector_u[0] + new_knot_size;
			//for (int i = 0; i<rGeometryRefinementParameters.multiply_knots_u - 1; i++)
			//{
			//	this_knot_insertions_u.resize(this_knot_insertions_u.size() + 1);
			//	this_knot_insertions_u[this_knot_insertions_u.size() - 1] = knot_insertion;
			//	knot_insertion += new_knot_size;
			//}
		}

		if (this_knot_insertions_u.size() > 0 ||
			this_knot_insertions_v.size() > 0)
			RefineKnotVector(this_knot_insertions_u, this_knot_insertions_v);


	}


	/* refines the the knot vectors with the additional knots Ru and Rv
	*  Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
	*  Algorithm A5.5
	*  @param[in] rRu knots in u-direction
	*  @param[in] rRv knots in v-direction
	*/
	void BrepFace::RefineKnotVector(const Vector& rRu, const Vector& rRv)
	{
		std::vector<double> Ru, Rv; //reduced vectors for knots that exceed the borders
		double tolerance = 10e-8;
		for (int i = 0; i < rRu.size(); i++) //
		{
			if ((m_knot_vector_u[0] - rRu[i]) < tolerance && (rRu[i] - m_knot_vector_u[m_knot_vector_u.size() - 1]) < tolerance)
			{
				double iu = rRu[i];
				Ru.push_back(iu);
			}
		}
		for (int i = 0; i < rRv.size(); i++)
		{
			if ((m_knot_vector_v[0] - rRv[i]) <tolerance && (rRv[i] - m_knot_vector_v[m_knot_vector_v.size() - 1]) < tolerance)
				Rv.push_back(rRv[i]);
		}
		int number_control_points_u = m_knot_vector_u.size() - m_p - 1;
		int number_control_points_v = m_knot_vector_v.size() - m_q - 1;

		matrix<array_1d<double, 4>> Pw(number_control_points_u, number_control_points_v); // Reference control points
		/* *---------*
		   |  v->    |
		   |   u|    |
		   |    v    |
		   *---------* */
		for (int i = 0; i < number_control_points_u; i++)
		{
			for (int j = 0; j < number_control_points_v; j++)
			{
				array_1d<double, 4> control_point;
				int control_point_index = j*number_control_points_u + i;
				control_point[0] = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->X();
				control_point[1] = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->Y();
				control_point[2] = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->Z();
				control_point[3] = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->GetValue(CONTROL_POINT_WEIGHT);

				mp_model_part->RemoveNodeFromAllLevels(m_control_points_ids[control_point_index]);

				Pw(i,j) = control_point;
			}
		}
		matrix<array_1d<double, 4>> Qw(number_control_points_u + Ru.size(), number_control_points_v); // New control points

		if (Ru.size() > 0)
		{
			int mu = number_control_points_u + m_p + 1;
			int a  = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, Ru[0]);
			int b  = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, Ru[Ru.size() - 1]) + 1;
			for (int col = 0; col < number_control_points_v; col++)
			{
				for (int j = 0; j <= a - m_p; j++)
				{
					Qw(j, col) = Pw(j, col);
				}
				for (int j = b - 1; j < number_control_points_u; j++)
				{
					Qw(j + Ru.size(), col) = Pw(j, col);
				}
			}

			Vector knot_vector_u = ZeroVector(m_knot_vector_u.size() + Ru.size());
			for (int j = 0; j <= a; j++) knot_vector_u[j] = m_knot_vector_u[j];
			for (int j = b + m_p; j <= mu; j++) knot_vector_u[j + Ru.size()-1] = m_knot_vector_u[j-1];

			int index_i = b + m_p - 1;
			int index_k = index_i + Ru.size();

			for (int j = Ru.size() - 1; j >= 0; j--)
			{
				while (Ru[j] <= m_knot_vector_u[index_i] && index_i > a)
				{
					for (int col = 0; col < number_control_points_v; col++)
					{
						Qw(index_k - m_p - 1, col) = Pw(index_i - m_p - 1, col);
					}
					knot_vector_u[index_k] = m_knot_vector_u[index_i];
					index_k--;
					index_i--;
				}
				for (int col = 0; col < number_control_points_v; col++)
				{
					Qw(index_k - m_p - 1, col) = Qw(index_k - m_p, col);
				}
				for (int l = 1; l <= m_p; l++)
				{
					int ind = index_k - m_p + l;
					double alfa = (Ru[j] - knot_vector_u[index_k + l]) / (m_knot_vector_u[index_i - m_p + l] - knot_vector_u[index_k + l]);
					for (int col = 0; col < number_control_points_v; col++)
					{
						Qw(ind - 1, col) = alfa*Qw(ind - 1, col) + (1 - alfa)*Qw(ind, col);
					}
				}
				knot_vector_u[index_k] = Ru[j];
				index_k = index_k - 1;
			}
			number_control_points_u += Ru.size();

			m_knot_vector_u.resize(m_knot_vector_u.size(), true);
			m_knot_vector_u = knot_vector_u;
			KRATOS_WATCH(Pw)
			Pw = Qw;
			KRATOS_WATCH(Pw)
		}
		Qw.resize(number_control_points_u, number_control_points_v + Rv.size(), true); // New control points
		std::cout << "here" << std::endl;
		if (Rv.size() > 0)
		{
			int mv = number_control_points_v + m_q + 1;
			int a = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, Rv[0]);
			int b = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, Rv[Rv.size() - 1]) + 1;

			for (int row = 0; row < number_control_points_u; row++)
			{
				for (int j = 0; j <= a - m_q; j++)
				{
					Qw(row, j) = Pw(row, j);
					KRATOS_ERROR_IF(Qw.size1() <= row) << "here";
					KRATOS_ERROR_IF(Qw.size2() <= j) << "here";
				}
				for (int j = b - 1; j < number_control_points_v; j++)
				{
					Qw(row, j + Rv.size()) = Pw(row, j);
					KRATOS_ERROR_IF(Qw.size1() <= row) << "here";
					KRATOS_ERROR_IF(Qw.size2() <= (j + Rv.size())) << "here";
				}
			}
			Vector knot_vector_v = ZeroVector(m_knot_vector_v.size() + Rv.size());
			for (int j = 0; j <= a; j++)
			{
				knot_vector_v[j] = m_knot_vector_v[j];
				KRATOS_ERROR_IF(knot_vector_v.size() <= j) << "here";
			}
			for (int j = b + m_q; j <= mv; j++) {
				//KRATOS_ERROR_IF(knot_vector_v.size() <= j + Rv.size()) << "here guess" << j + Rv.size() << " size: " << knot_vector_v.size();
				knot_vector_v[j + Rv.size()-1] = m_knot_vector_v[j-1];
				//KRATOS_ERROR_IF(m_knot_vector_v.size() == j ) << "here guess" << j;
			}

			KRATOS_WATCH(knot_vector_v)

			int index_i = b + m_q - 1;
			int index_k = index_i + Rv.size();

			for (int j = Rv.size() - 1; j >= 0; j--)
			{
				while ((Rv[j] <= m_knot_vector_v[index_i]) && (index_i > a))
				{
					for (int row = 0; row < number_control_points_u; row++)
					{
						Qw(row, index_k - m_q - 1) = Pw(row, index_i - m_q - 1);
						KRATOS_ERROR_IF(Qw.size1() <= row) << "here";
						KRATOS_ERROR_IF(Qw.size2() <= index_k - m_q - 1) << "here";
						KRATOS_ERROR_IF(Qw.size2() <= index_i - m_q -1) << "here";
					}
					knot_vector_v[index_k] = m_knot_vector_v[index_i];
					KRATOS_ERROR_IF(knot_vector_v.size() <= index_k) << "here";
					index_k--;
					index_i--;
				}
				for (int row = 0; row < number_control_points_u; row++)
				{
					Qw(row, index_k - m_q - 1) = Qw(row, index_k - m_q);
					KRATOS_ERROR_IF(Qw.size1() <= row) << "here";
					KRATOS_ERROR_IF(Qw.size2() <= index_k - m_q - 1) << "here";
					KRATOS_ERROR_IF(Qw.size2() <= index_k - m_q) << "here";
				}
				for (int l = 1; l <= m_q; l++)
				{
					int ind = index_k - m_q + l;
					double alfa = (Rv[j] - knot_vector_v[index_k + l]) / (m_knot_vector_v[index_i - m_q + l] - knot_vector_v[index_k + l]);
					KRATOS_ERROR_IF(knot_vector_v.size() <= index_k + l) << "here";
					for (int row = 0; row < number_control_points_u; row++)
					{
						Qw(row, ind - 1) = alfa*Qw(row, ind - 1) + (1 - alfa)*Qw(row, ind);
						KRATOS_ERROR_IF(Qw.size1() <= row) << "here";
						KRATOS_ERROR_IF(Qw.size2() <= ind - 1) << "here";
						KRATOS_ERROR_IF(Qw.size2() <= ind) << "here";
					}
				}
				knot_vector_v[index_k] = Rv[j];
				KRATOS_ERROR_IF(knot_vector_v.size() <= index_k) << "here";
				KRATOS_ERROR_IF(Rv.size() <= j) << "here";
				index_k = index_k - 1;
			}
			m_knot_vector_v.resize(knot_vector_v.size());
			m_knot_vector_v = knot_vector_v;
		}
		//KRATOS_WATCH(m_knot_vector_u)
		//KRATOS_WATCH(m_knot_vector_v)

		//m_knot_vector_u = knot_vector_u;
		//m_knot_vector_v = knot_vector_v;
		KRATOS_WATCH(m_knot_vector_u)
		KRATOS_WATCH(m_knot_vector_v)
		m_control_points_ids.resize(Qw.size1()*Qw.size2(), true);
		std::cout << "control point ids:";
		for (int i = 0; i < Qw.size1(); i++)
		{
			for (int j = 0; j < Qw.size2(); j++)
			{
				int new_id = 1;
				if (mp_model_part->Nodes().size() > 0)
					new_id = mp_model_part->GetRootModelPart().Nodes().back().Id() + 1;
				int control_point_index = j*Qw.size1() + i;
				m_control_points_ids[control_point_index] = new_id;

				std::cout << " " << new_id;

				mp_model_part->CreateNewNode(new_id, Qw(i, j)[0], Qw(i, j)[1], Qw(i, j)[2])->SetValue(CONTROL_POINT_WEIGHT, Qw(i, j)[3]);
			}
		}
		std::cout << std::endl;
		KRATOS_WATCH(Qw)
		KRATOS_WATCH(m_control_points_ids.size())
		KRATOS_WATCH(mp_model_part->NumberOfNodes())
  }

	/* elevates polinomial degrees p and q by tp and tq
	*  Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
	*  Algorithm A5.10
	*  @param[in] tp order elevation of polynomial degree p
	*  @param[in] tq order elevation of polynomial degree q
	*/
	void BrepFace::DegreeElevate(const int& tp, const int& tq)
	{
		int number_control_points_u = m_knot_vector_u.size() - m_p - 1;
		int number_control_points_v = m_knot_vector_v.size() - m_q - 1;

		matrix<array_1d<double, 4>> Pw(number_control_points_v, number_control_points_u); // Reference control points
		/*---------*
		 |  v-->   |
		 |u|       |
		 | v       |
		**---------**/
		for (int i = 0; i < number_control_points_v; i++)
		{
			for (int j = 0; j < number_control_points_u; j++)
			{
				array_1d<double, 4> control_point;
				int control_point_index = j*number_control_points_v + i;

				control_point[0] = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->X();
				control_point[1] = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->Y();
				control_point[2] = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->Z();
				control_point[3] = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->GetValue(CONTROL_POINT_WEIGHT);

				mp_model_part->RemoveNodeFromAllLevels(m_control_points_ids[control_point_index]);

				Pw(i,j) = control_point;
			}
		}
		matrix<array_1d<double, 4>> Qw(number_control_points_u + tp, number_control_points_v); // New control points

		// DEGREE ELEVATION p
		if (tp > 0)
		{
			Matrix bezalfs = ZeroMatrix(m_p + tp + 1, m_p + 1); // coefficients for degree elevating the Bézier segments
			matrix<array_1d<double, 4>> bpts(m_p + 1, number_control_points_v); //pth-degree Bézier control points of the current segment
			matrix<array_1d<double, 4>> ebpts(m_p + tp + 1, number_control_points_v); //((p+t)th-degree Bézier control points of the current segment
			matrix<array_1d<double, 4>> Nextbpts(m_p - 1, number_control_points_v); // leftmost control points of the next Bézier segment

			int m = number_control_points_v + m_p;
			int ph = m_p + tp;
			int ku = 1;
			//preallocate u_ele
			for (int i = 0; i<m-1; i++)
				if (m_knot_vector_u[i] != m_knot_vector_u[i + 1]) ku++;
			Vector u_ele = ZeroVector(m + 1 + tp*ku);

			bezalfs(0, 0) = 1.0;
			bezalfs(ph, m_p) = 1.0;

			for (int i = 1; i <= (ph / 2); i++) // bezalfs are symmetric to row ph/2
			{
				double inv = 1.0 / NurbsUtilities::binom(ph, i);
				double mpi = std::min((int) m_p, i);
				for (int j = std::max(0, i - tp); j <= mpi; j++)
				{
					bezalfs(i, j) = inv*NurbsUtilities::binom(m_p, j)*NurbsUtilities::binom(tp, i - j);
					bezalfs(ph - i, m_p - j) = bezalfs(i, j); //symmetry
				}
			}
			int mh = ph;
			int kind = ph + 1;
			int r = -1;
			int a = m_p;
			int b = m_p + 1;
			int cind = 1;
			double ua = m_knot_vector_u[a];
			for (int h = 0; h < number_control_points_v; h++)
				Qw(0, h) = Pw(0, h);
			for (int i = 0; i < ph + 1; i++) // left end of Ur
				u_ele[i] = ua;
			for (int i = 0; i < m_p + 1; i++) // initialize first Bezier seg
			{
				for (int h = 0; h < number_control_points_v; h++)
					bpts(i, h) = Pw(i, h);
			}
			while (b<m) // big loop through knot vector
			{
				int index_i = b;
				while (b<m && m_knot_vector_u[b] == m_knot_vector_u[b + 1])
					b++; // make b the rightmost ocurrence of ub
				int mult = b - index_i + 1; // multiplicity of ub
				mh = mh + mult + tp;
				double ub = m_knot_vector_u[b];
				int oldr = r; // r from last segment
				r = m_p - mult; // ub to be inserted r times

				int lbz = 0;
				int rbz = 0;

				if (oldr>0) lbz = ((oldr + 2) / 2);
				else lbz = 1;
				if (r>0)
					rbz = ph - ((r + 1) / 2);
				else
					rbz = ph;

				Vector alfs = ZeroVector(m_p - mult + 1);
				if (r>0)
				{
					double numer = ub - ua;
					for (int k = m_p; k>mult; k--)
						alfs[k - mult - 1] = numer / (m_knot_vector_u[a + k] - ua); // alfa for knot insertion
					for (int j = 1; j <= r; j++) // r times knot insertion
					{
						int save = r - j;
						int s = mult + j;
						for (int k = m_p; k >= s; k--) // new CP due to knot insertion
						{
							for (int h = 0; h < number_control_points_v; h++)
								bpts(k, h) = alfs[k - s] * bpts(k, h) + (1.0 - alfs[k - s])*bpts(k - 1, h);
						}
						for (int h = 0; h < number_control_points_v; h++)
							Nextbpts(save, h) = bpts(m_p, h);
					}
				} // end of inserting knots
				// degree elevate Bezier
				for (int i = lbz; i <= ph; i++)
				{
					// only points lbz,..,pr are used below
					for (int h = 0; h < number_control_points_v; h++)
					{
						ebpts(i, h)[0] = 0.0;
						ebpts(i, h)[1] = 0.0;
						ebpts(i, h)[2] = 0.0;
						ebpts(i, h)[3] = 0.0;
					}
					int mpi = std::min((int)m_p, i);
					for (int j = std::max(0, i - tp); j <= mpi; j++) // new CP due to degree elevation
					{
						for (int h = 0; h < number_control_points_v; h++)
							ebpts(i, h) = ebpts(i, h) + bezalfs(i, j)*bpts(j, h);
					}
				} // end of degree elevating Bezier
				if (oldr>1) // knot removal ua oldr times
				{
					int first = kind - 2;
					int last = kind;
					double den = ub - ua;
					double bet = (ub - u_ele[kind - 1]) / den;
					for (int tr = 1; tr<oldr; tr++)
					{
						index_i = first;
						int j = last;
						int kj = j - kind + 1;
						while ((j - index_i)>tr) // loop and compute the new CP for one removal step
						{
							if (index_i<cind)
							{
								double alf = (ub - u_ele[index_i]) / (ua - u_ele[index_i]);
								for (int h = 0; h < number_control_points_v; h++)
									Qw(index_i, h) = (alf*Qw(index_i, h) + (1.0 - alf)*Qw(index_i - 1, h));
							}
							if (j >= lbz)
							{
								if ((j - tr) <= (kind - ph + oldr))
								{
									double gam = (ub - u_ele[j - tr]) / den;
									for (int h = 0; h < number_control_points_v; h++)
										ebpts(kj, h) = gam*ebpts(kj, h) + (1.0 - gam)*ebpts(kj + 1, h);
								}
								else
								{
									for (int h = 0; h < number_control_points_v; h++)
										ebpts(kj, h) = bet*ebpts(kj, h) + (1.0 - bet)*ebpts(kj + 1, h);
								}
							}
							index_i++;
							j = j - 1;
							kj = kj - 1;
						}
						first--;
						last++;
					}
				} // end of removing knot
				if (a != m_p) // load the knot ua
				{
					for (int i = 0; i<(ph - oldr); i++)
					{
						u_ele[kind] = ua;
						kind = kind + 1;
					}
				}
				for (int j = lbz; j <= rbz; j++) // load CPs into Qw
				{
					if ((int)Qw.size1() == cind)
					{
						int y = Qw.size1();
						int x = Qw.size2();
						y++;
						Qw.resize(y, x, true);
					}
					for (int h = 0; h < number_control_points_v; h++)
						Qw(cind, h) = ebpts(j, h);
					cind = cind + 1;
				}
				if (b<m) // setup for the next pass through loop
				{
					for (int j = 0; j<r; j++)
						for (int h = 0; h < number_control_points_v; h++)
							bpts(j, h) = Nextbpts(j, h);
					for (int j = r; j <= m_p; j++)
					{
						for (int h = 0; h < number_control_points_v; h++)
							bpts(j, h) = Pw(b - m_p + j, h);
					}
					a = b;
					b = b + 1;
					ua = ub;
				}
				else // end knots
				{
					for (int i = 0; i <= ph; i++)
					{
						u_ele[kind + i] = ub;
					}
				}
			}
			m_p = ph;
			m_knot_vector_u.resize(u_ele.size(), true);
			m_knot_vector_u = u_ele;
			number_control_points_u = number_control_points_u + tp;

			Pw.resize(Qw.size1(),Qw.size2(), true);
			Pw = Qw;
		} // end of big loop through knot vector

		if (tq) {
			Qw.resize(number_control_points_u, number_control_points_v + tq, true);
			Matrix bezalfs = ZeroMatrix(m_q + tq + 1, m_q + 1); // coefficients for degree elevating the Bézier segments
			matrix<array_1d<double, 4>> bpts(number_control_points_u, m_q + 1); //pth-degree Bézier control points of the current segment
			matrix<array_1d<double, 4>> ebpts(number_control_points_u, m_q + tq + 1); //((p+t)th-degree Bézier control points of the current segment
			matrix<array_1d<double, 4>> Nextbpts(number_control_points_u, m_q - 1); // leftmost control points of the next Bézier segment

			int m = number_control_points_v + m_q;
			int qh = m_q + tq;
			int kv = 1;
			for (int i = 0; i < m - 1; i++) // preallocate v_ele
			{
				if (m_knot_vector_v[i] != m_knot_vector_v[i + 1]) kv++;
			}
			Vector v_ele = ZeroVector(m + 1 + tq*kv);

			// compute Bezier degree elevation coefficients
			bezalfs(0, 0) = 1.0;
			bezalfs(qh, m_q) = 1.0;

			for (int i = 1; i <= (qh / 2); i++) // bezalfs are symmetric to row qr/2
			{
				double inv = 1.0 / NurbsUtilities::binom(qh, i);
				int mpi = std::min((int)m_q, i);
				for (int j = std::max(0, i - tq); j <= mpi; j++) // alfas for degree elevation
				{
					bezalfs(i, j) = inv*NurbsUtilities::binom(m_q, j)*NurbsUtilities::binom(tq, i - j);
					bezalfs(qh - i, m_q - j) = bezalfs(i, j); // symmetry !
				}
			}
			int mh = qh;
			int kind = qh + 1;
			int r = -1;
			int a = m_q;
			int b = m_q + 1;
			int cind = 1;
			double ua = m_knot_vector_v[a];
			for (int h = 0; h < number_control_points_u; h++)
				Qw(h, 0) = Pw(h, 0);
			for (int i = 0; i < qh + 1; i++)  // left end of v_ele
			{
				v_ele[i] = ua;
			}
			// initialize first Bezier seg
			for (int i = 0; i < m_q + 1; i++)
			{
				for (int h = 0; h < number_control_points_u; h++)
				{
					bpts(h, i) = Pw(h, i);
				}
			}
			while (b < m) // big loop through knot vector
			{
				int index_i = b;
				while (b < m && m_knot_vector_v[b] == m_knot_vector_v[b + 1]) // make b the rightmost ocurrence of ub
				{
					b = b + 1;
				}
				int mult = b - index_i + 1; // multiplicity of ub
				mh = mh + mult + tq;
				double ub = m_knot_vector_v[b];
				int oldr = r; // r from last segment
				r = m_q - mult; // ub to be inserted r times

				int lbz = 0;
				int rbz = 0;

				if (oldr > 0)
					lbz = ((oldr + 2) / 2);
				else
					lbz = 1;

				if (r > 0) // insert knots to get Bezier segment
					rbz = qh - ((r + 1) / 2);
				else
					rbz = qh;
				Vector alfs = ZeroVector(m_q - mult + 1);

				if (r > 0)
				{
					double numer = ub - ua;
					for (int k = m_q; k > mult; k--)
					{
						alfs[k - mult - 1] = numer / (m_knot_vector_v[a + k] - ua); // alfa for knot insertion
					}
					for (int j = 1; j <= r; j++) // r times knot insertion
					{
						int save = r - j;
						int s = mult + j;
						for (int k = m_q; k >= s; k--) // new CP due to knot insertion
						{
							for (int h = 0; h < number_control_points_u; h++)
								bpts(h, k) = alfs[k - s] * bpts(h, k) + (1.0 - alfs[k - s])*bpts(h, k - 1);
						}
						for (int h = 0; h < number_control_points_u; h++)
							Nextbpts(h, save) = bpts(h, m_q);
					}
				} // end of inserting knots
				// degree elevate Bezier
				for (int i = lbz; i <= qh; i++)
				{
					// degree elevate Bezier
					for (int h = 0; h < number_control_points_u; h++)
					{
						ebpts(h, i)[0] = 0.0;
						ebpts(h, i)[1] = 0.0;
						ebpts(h, i)[2] = 0.0;
						ebpts(h, i)[3] = 0.0;
					}
					int mpi = std::min((int)m_q, i);
					for (int j = std::max(0, i - tq); j <= mpi; j++) // new CP due to degree elevation
					{
						for (int h = 0; h < number_control_points_u; h++)
							ebpts(h, i) = ebpts(h, i) + bezalfs(i, j)*bpts(h, j);
					}
				} // end of degree elevating Bezier
				if (oldr > 1) // knot removal ua oldr times
				{
					int first = kind - 2;
					int last = kind;
					double den = ub - ua;
					double bet = (ub - v_ele[kind - 1]) / den;
					for (int tr = 1; tr < oldr; tr++)
					{
						index_i = first;
						int j = last;
						int kj = j - kind + 1;
						while ((j - index_i) > tr) // loop and compute the new CP for one removal step
						{
							if (index_i < cind)
							{
								double alf = (ub - v_ele[index_i]) / (ua - v_ele[index_i]);
								for (int h = 0; h < number_control_points_u; h++)
									Qw(h, index_i) = (alf*Qw(h, index_i) + (1.0 - alf)*Qw(h, index_i - 1));
							}
							if (j >= lbz)
							{
								if ((j - tr) <= (kind - qh + oldr))
								{
									double gam = (ub - v_ele[j - tr]) / den;
									for (int h = 0; h < number_control_points_u; h++)
										ebpts(h, kj) = gam*ebpts(h, kj) + (1.0 - gam)*ebpts(h, kj + 1);
								}
								else
									for (int h = 0; h < number_control_points_u; h++)
										ebpts(h, kj) = bet*ebpts(h, kj) + (1.0 - bet)*ebpts(h, kj + 1);
							}
							index_i++;
							j--;
							kj--;
						}
						first = first - 1;
						last = last + 1;
					}
				} // end of removing knot

				if (a != m_q) // load the knot ua
				{
					KRATOS_WATCH(qh)
					KRATOS_WATCH(oldr)
					KRATOS_WATCH(ua)
					for (int i = 0; i < qh - oldr; i++)
					{
						v_ele[kind] = ua;
						kind = kind + 1;
					}
				}
				for (int j = lbz; j <= rbz; j++) // load CPs into Qw
				{
					if ((int)Qw.size2() == cind)
					{
						int y = Qw.size1();
						int x = Qw.size2();
						x++;
						Qw.resize(y, x, true);
					}
					for (int h = 0; h < number_control_points_u; h++)
						Qw(h, cind) = ebpts(h, j);
					cind = cind + 1;

				}
				if (b < m) // setup for the next pass through loop
				{
					for (int j = 0; j < r; j++)
					{
						for (int h = 0; h < number_control_points_u; h++)
							bpts(h, j) = Nextbpts(h, j);
					}
					for (int j = r; j <= m_q; j++)
					{
						for (int h = 0; h < number_control_points_u; h++)
						{
							bpts(h, j) = Pw(h, b - m_q + j);
						}
					}
					a = b;
					b = b + 1;
					ua = ub;
				}
				else // end knots
				{
					KRATOS_WATCH(qh)
						KRATOS_WATCH(kind)
					for (int i = 0; i <= qh; i++)
					{
						v_ele[kind + i] = ub;
					}
				}
			}// end of big loop through knot vector

			m_q = qh;
			m_knot_vector_v.resize(v_ele.size(), true);
			m_knot_vector_v = v_ele;
		}

		m_control_points_ids.resize(Qw.size1()*Qw.size2(), true);
		std::cout << "Control ids: ";
		for (int i = 0; i < Qw.size1(); i++)
		{
			for (int j = 0; j < Qw.size2(); j++)
			{
				int new_id = 1;
				if (mp_model_part->Nodes().size() > 0)
					new_id = mp_model_part->GetRootModelPart().Nodes().back().Id() + 1;
				int control_point_index = j*Qw.size1() + i;
				m_control_points_ids[control_point_index] = new_id;

				std::cout << " " << new_id;

				mp_model_part->CreateNewNode(new_id, Qw(i, j)[0], Qw(i, j)[1], Qw(i, j)[2])->SetValue(CONTROL_POINT_WEIGHT, Qw(i, j)[3]);
			}
		}
		std::cout << std::endl;

		KRATOS_WATCH(Qw)

		KRATOS_WATCH(m_p)
		KRATOS_WATCH(m_q)

		KRATOS_WATCH(m_knot_vector_u)
		KRATOS_WATCH(m_knot_vector_v)
	}


  //GEOMETRY FUNCTIONS:
  /* @Author Daniel Baumgaertner
  *  @date   December, 2016
  *  returns the cartesian coordinates (global) for a specific point
  *  located on the NURBS surface S(u=fixed and v=fixed)
  *  Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
  *  Algorithm A4.3
  *
  *  @param[in]  rSurfacePoint  evaluated point
  *  @param[in]  u  local parameter in u-direction
  *  @param[in]  v  local parameter in v-direction
  */
  void BrepFace::EvaluateSurfacePoint(Point& rSurfacePoint, const double& u, const double& v)
  {
    rSurfacePoint[0] = 0;
    rSurfacePoint[1] = 0;
    rSurfacePoint[2] = 0;

    int span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, u);
    int span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, v);

    //Vector ShapeFunctionsN = ZeroVector((m_q + 1)*(m_p + 1));
    Matrix R;
    EvaluateNURBSFunctions(span_u, span_v, u, v, R);

    for (int c = 0; c <= m_q; c++)
    {
      for (int b = 0; b <= m_p; b++)
      {
        // the control point vector is filled up by first going over u, then over v
        int ui = span_u - m_p + b;
        int vi = span_v - m_q + c;
        int m_n_u = m_knot_vector_u.size() - m_p - 1;
        int control_point_index = vi*m_n_u + ui;

		KRATOS_ERROR_IF(control_point_index > m_control_points_ids.size() - 1) << "There is a bug" << std::endl;
        rSurfacePoint[0] += R(b, c) * mp_model_part->pGetNode(m_control_points_ids[control_point_index])->X();
        rSurfacePoint[1] += R(b, c) * mp_model_part->pGetNode(m_control_points_ids[control_point_index])->Y();
        rSurfacePoint[2] += R(b, c) * mp_model_part->pGetNode(m_control_points_ids[control_point_index])->Z();
      }
    }
  }
  // #######################################################################################
  //
  //  \details    evaluate Hessian and Gradient modifying the input objects
  //
  // ======================================================================================
  //  \param[in]  QminP    	 	Distance Vector
  //  \param[in]  H		     	Hessian reference
  //  \param[in]  Gradient    	Gradient reference
  //  \param[in]  v    			parameter
  //  \param[in]  u 				parameter
  //
  // ======================================================================================
  //  \author     Giovanni Filomeno (1/2017) && Massimo Sferza (1/2017)
  //
  //########################################################################################
  void BrepFace::EvaluateGradientsForClosestPointSearch(Vector QminP, Matrix& Hessian, Vector& Gradient, double& u, double& v)
  {
    // The derivatives of the basis functions are evaluated
    Matrix dR;
    Matrix ddR;
    EvaluateNURBSFunctionsDerivatives(-1, -1, u, v, dR, ddR);

    // The derivatives of Q(u,v) are evaluated
    Vector dQdu = ZeroVector(3);
    Vector dQdv = ZeroVector(3);
    Vector dQdudu = ZeroVector(3);
    Vector dQdvdv = ZeroVector(3);
    Vector dQdudv = ZeroVector(3);

    int span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, u);
    int span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, v);

    int k = 0;
    for (int c = 0; c <= m_q; c++)
    {
      for (int b = 0; b <= m_p; b++)
      {
        int ui = span_u - m_p + b;
        int vi = span_v - m_q + c;
        int m_n_u = m_knot_vector_u.size() - m_p - 1;
        int control_point_index = vi*m_n_u + ui;

        double cp_x = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->X();
        double cp_y = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->Y();
        double cp_z = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->Z();

        dQdu(0) += dR(k, 0) * cp_x;
        dQdu(1) += dR(k, 0) * cp_y;
        dQdu(2) += dR(k, 0) * cp_z;

        dQdv(0) += dR(k, 1) * cp_x;
        dQdv(1) += dR(k, 1) * cp_y;
        dQdv(2) += dR(k, 1) * cp_z;

        dQdudu(0) += ddR(k, 0) * cp_x;
        dQdudu(1) += ddR(k, 0) * cp_y;
        dQdudu(2) += ddR(k, 0) * cp_z;

        dQdvdv(0) += ddR(k, 1) * cp_x;
        dQdvdv(1) += ddR(k, 1) * cp_y;
        dQdvdv(2) += ddR(k, 1) * cp_z;

        dQdudv(0) += ddR(k, 2) * cp_x;
        dQdudv(1) += ddR(k, 2) * cp_y;
        dQdudv(2) += ddR(k, 2) * cp_z;

        k++;
      }
    }
    // Hessian and gradient are evaluated
    Hessian(0, 0) = 2 * (inner_prod(dQdudu, QminP) + inner_prod(dQdu, dQdu));
    Hessian(0, 1) = 2 * (inner_prod(dQdudv, QminP) + inner_prod(dQdu, dQdv));
    Hessian(1, 0) = 2 * (inner_prod(dQdudv, QminP) + inner_prod(dQdu, dQdv));
    Hessian(1, 1) = 2 * (inner_prod(dQdvdv, QminP) + inner_prod(dQdv, dQdv));

    Gradient(0) = 2 * inner_prod(dQdu, QminP);
    Gradient(1) = 2 * inner_prod(dQdv, QminP);
  }
  //  #####################################################################################
  // #######################################################################################
  //
  //  \details    returns the basis functions of NURBS basis function w.r.t. u,v
  //              span_u, span_v are the knot span indices. if unknown, insert 0!
  //
  // ======================================================================================
  //  \param[in]  span_u     knotspan index in u-direction
  //  \param[in]  span_v     knotspan index in v-direction
  //  \param[in]  _u         local parameter in u-direction
  //  \param[in]  _v         local parameter in v-direction
  //  \param[out] R         basis func
  //
  // ======================================================================================
  //  \author     Daniel Baumgärtner (12/2016)
  //
  //########################################################################################
  void BrepFace::EvaluateNURBSFunctions(int span_u, int span_v, double _u, double _v, Matrix& R)
  {
    if (span_u == -1) span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, _u);
    if (span_v == -1) span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, _v);

    Vector N;
    Vector M;

    R.resize(m_p + 1, m_q + 1, true);
    noalias(R) = ZeroMatrix(m_p + 1, m_q + 1);

    // Evaluate basis functions with derivatives
    NurbsUtilities::eval_nonzero_basis_function(N, m_knot_vector_u, _u, span_u, m_p);
    NurbsUtilities::eval_nonzero_basis_function(M, m_knot_vector_v, _v, span_v, m_q);

    double sum = 0.0;

    for (int c = 0; c <= m_q; c++)
    {
      for (int b = 0; b <= m_p; b++)
      {
        // the control point vector is filled up by first going over u, then over v
        int ui = span_u - m_p + b;
        int vi = span_v - m_q + c;
        int n_u = m_knot_vector_u.size() - m_p - 1;
		int control_point_index = vi * n_u + ui;
        // Evaluate basis function
		//if (control_point_index > m_control_points_ids.size() - 1)
		//{
		//	KRATOS_WATCH(_u)
		//	KRATOS_WATCH(_v)
		//	KRATOS_WATCH(m_control_points_ids.size())
		//	KRATOS_WATCH(span_u)
		//	KRATOS_WATCH(span_v)
		//	KRATOS_WATCH(m_knot_vector_u)
		//	KRATOS_WATCH(m_knot_vector_v)
		//	KRATOS_WATCH(n_u)
		//	KRATOS_WATCH(m_p)
		//	KRATOS_WATCH(m_q)
		//	KRATOS_WATCH(m_control_points_ids)
		//	KRATOS_WATCH(control_point_index)
		//	KRATOS_ERROR << "There is a bug" << std::endl;
		//}
		
        R(b, c) = N(b)*M(c)*mp_model_part->pGetNode(m_control_points_ids[control_point_index])->GetValue(CONTROL_POINT_WEIGHT);
        sum += R(b, c);
      }
    }

    // divide by sum only required in terms of rational basis functions
    //if (std::abs(sum-weight)> cepsilon) //Breitenberger 18.06.2014
    double inv_sum = 1 / sum;
    // divide through by sum
    for (int c = 0; c <= m_q; c++)
      for (int b = 0; b <= m_p; b++)
        R(b, c) = inv_sum*R(b, c);
  }

  /**
  * @Author M.Breitenberger in Carat (12/2009)
  * @date   March, 2017
  * @brief   rreturns the first and second derivative of NURBS basis function w.r.t. u,v
  *              span_u,span_v are the knot span indices. if unknown, insert -1!
  *
  * @param[in]  span_u    knotspan index in u-direction
  * @param[in]  u     	  local parameter in u-direction
  * @param[in]  span_v 	  knotspan index in v-direction
  * @param[in]  v         local parameter in v-direction
  * @param[out] DN_De     1st derivatives
  * @param[out] DDN_DDe   2nd derivatives
  */
  void BrepFace::EvaluateNURBSFunctionsDerivatives(int span_u, int span_v, double u, double v,
    Matrix& DN_De, Matrix& DDN_DDe)
  {
    if (span_u == -1) span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, u);
    if (span_v == -1) span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, v);

    int number_of_control_points = (m_p + 1)*(m_q + 1); // Control Points per element
    Matrix N;              // Basisfunc at _u
    Matrix M;              // Basisfunc at _v
    NurbsUtilities::eval_nonzero_basis_function_with_derivatives(N, m_knot_vector_u, u, span_u, m_p, 2);
    NurbsUtilities::eval_nonzero_basis_function_with_derivatives(M, m_knot_vector_v, v, span_v, m_q, 2);

    vector<double> r(number_of_control_points);
    r.clear();
    DN_De.resize(number_of_control_points, 2, true);
    DDN_DDe.resize(number_of_control_points, 3, true);

    double sum = 0.0;
    Vector dsum = ZeroVector(2);
    Vector ddsum = ZeroVector(3);
    double weight;

    int k = 0;
    for (int c = 0; c <= m_q; c++)
    {
      for (int b = 0; b <= m_p; b++)
      {
        // the control point vector is filled up by first going over u, then over v
        int ui = span_u - m_p + b;
        int vi = span_v - m_q + c;
        int m_n_u = m_knot_vector_u.size() - m_p - 1;
        int control_point_index = vi*m_n_u + ui;

		KRATOS_ERROR_IF(control_point_index > m_control_points_ids.size() - 1) << "There is a bug" << std::endl;
        // Evaluate basis function
        weight = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->GetValue(CONTROL_POINT_WEIGHT);

        r[k] = N(0, b)*M(0, c)*weight;
        sum += r[k];
        //First derivatives
        DN_De(k, 0) = N(1, b)*M(0, c)*weight;
        dsum[0] += DN_De(k, 0);
        DN_De(k, 1) = N(0, b)*M(1, c)*weight;
        dsum(1) += DN_De(k, 1);
        //Second derivatives  1-du^2, 2-dv^2, 3-dudv
        DDN_DDe(k, 0) = N(2, b)*M(0, c)*weight;
        ddsum(0) = ddsum(0) + DDN_DDe(k, 0);
        DDN_DDe(k, 1) = N(0, b)*M(2, c)*weight;
        ddsum(1) = ddsum(1) + DDN_DDe(k, 1);
        DDN_DDe(k, 2) = N(1, b)*M(1, c)*weight;
        ddsum(2) = ddsum(2) + DDN_DDe(k, 2);
        k++;
      }
    }
    //double sum_2 = pow(sum, 2);
    //double sum_3 = pow(sum, 3);
    double sum_2 = pow(sum, 2);
    double sum_3 = sum_2*sum;
    double inv_sum = 1.0 / sum;
    double inv_sum_2 = 1.0 / sum_2;
    double inv_sum_3 = 1.0 / sum_3;
    // divide through by sum
    for (int k = 0; k<number_of_control_points; k++)
    {
      DDN_DDe(k, 0) = DDN_DDe(k, 0)*inv_sum - 2.0*DN_De(k, 0)*dsum[0] * inv_sum_2
        - r[k] * ddsum[0] * inv_sum_2 + 2.0*r[k] * dsum[0] * dsum[0] * inv_sum_3;
      DDN_DDe(k, 1) = DDN_DDe(k, 1)*inv_sum - 2.0*DN_De(k, 1)*dsum[1] * inv_sum_2
        - r[k] * ddsum[1] * inv_sum_2 + 2.0*r[k] * dsum[1] * dsum[1] * inv_sum_3;
      DDN_DDe(k, 2) = DDN_DDe(k, 2)*inv_sum - DN_De(k, 0)*dsum[1] * inv_sum_2 - DN_De(k, 1)*dsum[0] * inv_sum_2
        - r[k] * ddsum[2] * inv_sum_2 + 2.0*r[k] * dsum[0] * dsum[1] * inv_sum_3;
      DN_De(k, 0) = DN_De(k, 0)*inv_sum - r[k] * dsum[0] * inv_sum_2;
      DN_De(k, 1) = DN_De(k, 1)*inv_sum - r[k] * dsum[1] * inv_sum_2;
      //DDN_DDe(k, 0) = DDN_DDe(k, 0) / sum - 2.0*DN_De(k, 0)*dsum[0] / sum_2
      //  - r[k] * ddsum[0] / sum_2
      //  + 2.0*r[k] * dsum[0] * dsum[0] / sum_3;
      //DDN_DDe(k, 1) = DDN_DDe(k, 1) / sum - 2.0*DN_De(k, 1)*dsum[1] / sum_2
      //  - r[k] * ddsum[1] / sum_2
      //  + 2.0*r[k] * dsum[1] * dsum[1] / sum_3;
      //DDN_DDe(k, 2) = DDN_DDe(k, 2) / sum - DN_De(k, 0)*dsum[1] / sum_2
      //  - DN_De(k, 1)*dsum[0] / sum_2
      //  - r[k] * ddsum[2] / sum_2
      //  + 2.0*r[k] * dsum[0] * dsum[1] / sum_3;
      //DN_De(k, 0) = DN_De(k, 0) / sum - r[k] * dsum[0] / sum_2;
      //DN_De(k, 1) = DN_De(k, 1) / sum - r[k] * dsum[1] / sum_2;
    }
  }

  //  #####################################################################################
  // #######################################################################################
  //
  //  \details    returns the basis fucntions and the first derivative of NURBS basis function w.r.t. u,v
  //              _i,_j are the knot span indices. if unknown, insert 0!
  //
  // ======================================================================================
  //  \param[in]  _i         knotspan index in u-direction
  //  \param[in]  _u         local parameter in u-direction
  //  \param[in]  _j         knotspan index in v-direction
  //  \param[in]  _v         local parameter in v-direction
  //  \param[out] _R         basis func
  //  \param[out] _dR        1st derivatives
  //
  // ======================================================================================
  //  \author     from M.Breitenberger in Carat (12/2009)
  //
  //########################################################################################
  void BrepFace::EvaluateNURBSFunctionsAndDerivative(int span_u, int span_v, double _u, double _v, Matrix& R, std::vector<Matrix>& dR)
  {
    if (span_u == -1) span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, _u);
    if (span_v == -1) span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, _v);

    Matrix N_matrix;
    Matrix M_matrix;
    NurbsUtilities::eval_nonzero_basis_function_with_derivatives(N_matrix, m_knot_vector_u, _u, span_u, m_p, 1);
    NurbsUtilities::eval_nonzero_basis_function_with_derivatives(M_matrix, m_knot_vector_v, _v, span_v, m_q, 1);
    double sum = 0.0;
    double dsum1 = 0.0;
    double dsum2 = 0.0;
    double weight;

    R.resize(m_p + 1, m_q + 1);
    dR.resize(2);
    dR[0].resize(m_p + 1, m_q + 1);
    dR[1].resize(m_p + 1, m_q + 1);

    for (int c = 0; c <= m_q; c++)
    {
      for (int b = 0; b <= m_p; b++)
      {
        // the control point vector is filled up by first going over u, then over v
        int ui = span_u - m_p + b;
        int vi = span_v - m_q + c;
        int m_n_u = m_knot_vector_u.size() - m_p - 1;
        int control_point_index = vi*m_n_u + ui;

		KRATOS_ERROR_IF(control_point_index > m_control_points_ids.size() - 1) << "There is a bug" << std::endl;
        // Evaluate basis function
        weight = mp_model_part->pGetNode(m_control_points_ids[control_point_index])->GetValue(CONTROL_POINT_WEIGHT);
        R(b, c) = N_matrix(0, b)*M_matrix(0, c)*weight;
        sum += R(b, c);

        //First derivatives
        dR[0](b, c) = N_matrix(1, b)*M_matrix(0, c)*weight;
        dsum1 += dR[0](b, c);
        dR[1](b, c) = N_matrix(0, b)*M_matrix(1, c)*weight;
        dsum2 += dR[1](b, c);
      }
    }

    // divide by sum only required in terms of rational basis functions
    double inv_sum = 1.0 / sum;
    // divide through by sum
    for (int c = 0; c <= m_q; c++)
    {
      for (int b = 0; b <= m_p; b++)
      {
        R(b, c) = inv_sum*R(b, c);
        dR[0](b, c) = inv_sum*dR[0](b, c) - R(b, c)*dsum1*inv_sum;
        dR[1](b, c) = inv_sum*dR[1](b, c) - R(b, c)*dsum2*inv_sum;
      }
    }
  }


	///Constructor
	BrepFace::BrepFace(unsigned int brep_id,
		bool is_trimmed,
		bool is_rational,
		TrimmingLoopVector& trimming_loops,
		TrimmingLoopVector& embedded_loops,
		std::vector<EmbeddedPoint>& embedded_points,
		Vector& knot_vector_u,
		Vector& knot_vector_v,
		unsigned int& p,
		unsigned int& q,
		IntVector& control_point_ids,
		Kratos::shared_ptr<ModelPart> model_part)
		  : m_trimming_loops(trimming_loops),
			m_is_trimmed(is_trimmed),
			m_is_rational(is_rational),
			m_embedded_loops(embedded_loops),
			m_embedded_points(embedded_points),
			m_knot_vector_u(knot_vector_u),
			m_knot_vector_v(knot_vector_v),
			m_p(p),
			m_q(q),
			mp_model_part(model_part),
			m_control_points_ids(control_point_ids),
			IndexedObject(brep_id),
			Flags()
	{
	}
	///Destructor
	BrepFace::~BrepFace()
	{}
} // namespace Kratos.