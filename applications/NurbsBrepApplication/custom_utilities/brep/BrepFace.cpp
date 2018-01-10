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

  IntVector BrepFace::GetIntegerVector(const Vector& vector, const int& tolerance)
  {
    //int tolerance = 10e9;
    IntVector new_vector;
    new_vector.resize(vector.size());

    for (unsigned int i = 0; i < vector.size(); i++)
    {
      new_vector[i] = (int)(vector[i] * tolerance);
      //KRATOS_WATCH(new_vector[i])
    }

    return new_vector;
  }

  //IntVector BrepFace::GetIntegerVKnotVector(const int& tolerance)
  //{
  //  //int tolerance = 10e9;
  //  IntVector VKnotVector;
  //  VKnotVector.resize(m_knot_vector_v.size());

  //  for (unsigned int i = 0; i < m_knot_vector_v.size(); i++)
  //  {
  //    VKnotVector[i] = (int)(m_knot_vector_v[i]*tolerance);
  //    KRATOS_WATCH(VKnotVector[i])
  //  }
  //  return VKnotVector;
  //}
  // --------------------------------------------------------------------------
  std::vector<Node<3>::Pointer> BrepFace::GetQuadraturePoints(const int& shapefunction_order)
  {
    int tolerance = 10e7;

    std::vector<Node<3>::Pointer> NodeVector;

    IntVector knot_vector_u = GetIntegerVector(m_knot_vector_u, tolerance);
    IntVector knot_vector_v = GetIntegerVector(m_knot_vector_v, tolerance);

    Vector parameter_span_u = ZeroVector(2);
    Vector parameter_span_v = ZeroVector(2);

    for (unsigned int i = 0; i < knot_vector_u.size() - 1; i++)
    {
      if (abs(knot_vector_u[i + 1] - knot_vector_u[i]) > 1)
      {
        parameter_span_u[0] = m_knot_vector_u[i];
        parameter_span_u[1] = m_knot_vector_u[i + 1];

        for (unsigned int j = 0; j < knot_vector_v.size() - 1; j++)
        {
          if (abs(knot_vector_v[j + 1] - knot_vector_v[j]) > 1)
          {
            parameter_span_v[0] = m_knot_vector_v[j];
            parameter_span_v[1] = m_knot_vector_v[j + 1];

            KnotSpan2d knot_span(0, true, m_p, m_q, parameter_span_u, parameter_span_v);
            std::vector<array_1d<double, 3>> points = knot_span.getIntegrationPointsInParameterDomain();
            //std::cout << "size of points: " << points.size() << std::endl;
            std::vector<Node<3>::Pointer> NodeVectorElement = EnhanceShapeFunctions(points, shapefunction_order);
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

  std::vector<Node<3>::Pointer> BrepFace::GetQuadraturePointsTrimmed(const int& shapefunction_order)
  {
    //std::vector<array_1d<double, 2>> boundary_polygon;
    //for (unsigned int loop_i = 0; loop_i < m_trimming_loops.size(); loop_i++)
    //{
    //  boundary_polygon = m_trimming_loops[loop_i].GetBoundaryPolygon(5);
    //}
    Polygon boundaries(m_trimming_loops);

    std::vector<Node<3>::Pointer> NodeVector;

    int tolerance = 10e6;

    IntVector knot_vector_u = GetIntegerVector(m_knot_vector_u, tolerance);
    IntVector knot_vector_v = GetIntegerVector(m_knot_vector_v, tolerance);

    Vector parameter_span_u = ZeroVector(2);
    Vector parameter_span_v = ZeroVector(2);

    for (unsigned int i = 0; i < knot_vector_u.size() - 1; i++)
    {
      if (abs(knot_vector_u[i + 1] - knot_vector_u[i]) > 1)
      {
        parameter_span_u[0] = m_knot_vector_u[i];
        parameter_span_u[1] = m_knot_vector_u[i + 1];

        for (unsigned int j = 0; j < knot_vector_v.size() - 1; j++)
        {
          if (abs(knot_vector_v[j + 1] - knot_vector_v[j]) > 1)
          {
            parameter_span_v[0] = m_knot_vector_v[j];
            parameter_span_v[1] = m_knot_vector_v[j + 1];

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
            //std::cout << "k: " << std::endl;
            //for (unsigned int k = 0; k < points.size(); k++)
            //{
            //  std::cout << "k: " << k << std::endl;
            //  std::cout << "points: " << points[k][0] << points[k][1] << points[k][2] << std::endl;
            //}
            //std::cout << "size of points: " << points.size() << std::endl;
            std::vector<Node<3>::Pointer> NodeVectorElement = EnhanceShapeFunctions(points, shapefunction_order);
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



  std::vector<Node<3>::Pointer> BrepFace::GetQuadraturePointsEmbedded(const int& shapefunction_order)
  {
    //std::vector<array_1d<double, 2>> boundary_polygon;
    //for (unsigned int loop_i = 0; loop_i < m_trimming_loops.size(); loop_i++)
    //{
    //  boundary_polygon = m_trimming_loops[loop_i].GetBoundaryPolygon(5);
    //}
    Polygon boundaries(m_embedded_loops);

    std::vector<Node<3>::Pointer> NodeVector;

    int tolerance = 10e6;

    IntVector knot_vector_u = GetIntegerVector(m_knot_vector_u, tolerance);
    IntVector knot_vector_v = GetIntegerVector(m_knot_vector_v, tolerance);

    Vector parameter_span_u = ZeroVector(2);
    Vector parameter_span_v = ZeroVector(2);

    for (unsigned int i = 0; i < knot_vector_u.size() - 1; i++)
    {
      if (abs(knot_vector_u[i + 1] - knot_vector_u[i]) > 1)
      {
        parameter_span_u[0] = m_knot_vector_u[i];
        parameter_span_u[1] = m_knot_vector_u[i + 1];

        for (unsigned int j = 0; j < knot_vector_v.size() - 1; j++)
        {
          if (abs(knot_vector_v[j + 1] - knot_vector_v[j]) > 1)
          {
            parameter_span_v[0] = m_knot_vector_v[j];
            parameter_span_v[1] = m_knot_vector_v[j + 1];

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
            //std::cout << "k: " << std::endl;
            //for (unsigned int k = 0; k < points.size(); k++)
            //{
            //  std::cout << "k: " << k << std::endl;
            //  std::cout << "points: " << points[k][0] << points[k][1] << points[k][2] << std::endl;
            //}
            //std::cout << "size of points: " << points.size() << std::endl;
            std::vector<Node<3>::Pointer> NodeVectorElement = EnhanceShapeFunctions(points, shapefunction_order);
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



  std::vector<Node<3>::Pointer> BrepFace::GetQuadraturePointsOfTrimmingCurve(const int& shapefunction_order, const int& trim_index)
  {
    BrepTrimmingCurve trimming_curve = GetTrimmingCurve(trim_index);
    //trimming_curve.PrintData();

    std::vector<double> intersections = trimming_curve.FindIntersections(m_p, m_q, m_knot_vector_u, m_knot_vector_v);

    //int highest_polynomial_order = m_p;
    //if (m_q > m_p)
    //  highest_polynomial_order = m_q;
    int degree = std::max(m_p, m_q);

    std::vector<array_1d<double, 3>> quadrature_points = trimming_curve.GetQuadraturePoints(intersections, degree);

    std::vector<Node<3>::Pointer> NodeVectorElement = EnhanceShapeFunctions(quadrature_points, shapefunction_order);
    for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
    {
      array_1d<double, 2> tangents_basis = trimming_curve.GetBaseVector(quadrature_points[k][2]);
      Vector tangents_basis_vector(2);
      tangents_basis_vector[0] = tangents_basis[0];
      tangents_basis_vector[1] = tangents_basis[1];
      NodeVectorElement[k]->SetValue(TANGENTS_BASIS_VECTOR, tangents_basis_vector);
      //KRATOS_WATCH(tangents_basis_vector)
    }
    return NodeVectorElement;
  }

  void BrepFace::GetLocalParameterOfPointOnTrimmingCurve(const Point& point, const BrepTrimmingCurve& trimming_curve, double& u, double& v)
  {
    bool success = NewtonRaphson(point, u, v);
  }

  std::vector<Node<3>::Pointer> BrepFace::GetQuadraturePointsOfTrimmingCurveWithPoints(const int& shapefunction_order, const int& trim_index, std::vector<Point> intersection_points)
  {
    BrepTrimmingCurve trimming_curve = GetTrimmingCurve(trim_index);
    //trimming_curve.PrintData();

    std::vector<array_1d<double, 2>> intersection_points_2d;
    for (unsigned int i = 0; i < intersection_points.size(); i++)
    {
      //GetLocalParameterOfPointOnTrimmingCurve(intersection_points[i], trimming_curve, )
      double u = 0;
      double v = 0;
      bool success = NewtonRaphson(intersection_points[i], u, v);
      //GetClosestPoint(intersection_points[i], u, v);
      if (!success)
      {
        std::vector<array_1d<double, 2>> trimming_curve_loop = trimming_curve.CreatePolygon(40);
        double distance = 1e10;
        array_1d<double, 2> initial_guess;
        for (auto point_itr = trimming_curve_loop.begin(); point_itr != trimming_curve_loop.end(); ++point_itr)
        {
          Point point_global;
          EvaluateSurfacePoint(point_global, (*point_itr)[0], (*point_itr)[1]);
          double new_distance = sqrt((intersection_points[i][0] - point_global[0]) * (intersection_points[i][0] - point_global[0])
            + (intersection_points[i][1] - point_global[1]) * (intersection_points[i][1] - point_global[1])
            + (intersection_points[i][2] - point_global[2]) * (intersection_points[i][2] - point_global[2]));
          if (new_distance < distance)
          {
            distance = new_distance;
            initial_guess = (*point_itr);
          }
        }
        u = initial_guess[0];
        v = initial_guess[1];
        bool success = NewtonRaphson(intersection_points[i], u, v);
        if (!success) {
          std::cout << "100 iteration points needed." << std::endl;
          trimming_curve_loop = trimming_curve.CreatePolygon(100);
          distance = 1e10;
          for (auto point_itr = trimming_curve_loop.begin(); point_itr != trimming_curve_loop.end(); ++point_itr)
          {
            Point point_global;
            EvaluateSurfacePoint(point_global, (*point_itr)[0], (*point_itr)[1]);
            double new_distance = sqrt((intersection_points[i][0] - point_global[0]) * (intersection_points[i][0] - point_global[0])
              + (intersection_points[i][1] - point_global[1]) * (intersection_points[i][1] - point_global[1])
              + (intersection_points[i][2] - point_global[2]) * (intersection_points[i][2] - point_global[2]));
            if (new_distance < distance)
            {
              distance = new_distance;
              initial_guess = (*point_itr);
            }
          }
          u = initial_guess[0];
          v = initial_guess[1];
          std::cout << "Initial: u=" << u << ", v=" << v << std::endl;
          bool success = NewtonRaphson(intersection_points[i], u, v);
          if (!success)
          {
            u = initial_guess[0];
            v = initial_guess[1];
            //KRATOS_ERROR << "Point not found after 100 pts." << std::endl;
          }
        }
      }
      array_1d<double, 2> point(u, v);
      intersection_points_2d.push_back(point);
    }

    std::vector<double> intersections_slave = trimming_curve.FindIntersectionsWithPoints(intersection_points_2d);
    std::cout << "intersections slave: ";
    for (unsigned int i = 0; i < intersections_slave.size(); i++)
    {
      std::cout << intersections_slave[i] << ", ";
    }
    std::cout << std::endl;
    std::vector<double> intersections_master = trimming_curve.FindIntersections(m_p, m_q, m_knot_vector_u, m_knot_vector_v);
    std::cout << "intersections master: ";
    for (unsigned int i = 0; i < intersections_master.size(); i++)
    {
      std::cout << intersections_master[i] << ", ";
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < intersections_master.size(); i++)
    {
      intersections_slave.push_back(intersections_master[i]);
    }
    std::sort(intersections_slave.begin(), intersections_slave.end());
    std::vector<double> intersections;
    intersections.push_back(intersections_slave[0]);
    for (unsigned int i = 1; i < intersections_slave.size(); i++)
    {
      if (intersections_slave[i - 1] != intersections_slave[i])
        intersections.push_back(intersections_slave[i]);
    }

    int highest_polynomial_order = std::max(m_p, m_q);

    std::vector<array_1d<double, 3>> quadrature_points = trimming_curve.GetQuadraturePoints(intersections, highest_polynomial_order);

    std::vector<Node<3>::Pointer> NodeVectorElement = EnhanceShapeFunctions(quadrature_points, shapefunction_order);
    for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
    {
      array_1d<double, 2> tangents_basis = trimming_curve.GetBaseVector(quadrature_points[k][2]);
      Vector tangents_basis_vector(2);
      tangents_basis_vector[0] = tangents_basis[0];
      tangents_basis_vector[1] = tangents_basis[1];
      NodeVectorElement[k]->SetValue(TANGENTS_BASIS_VECTOR, tangents_basis_vector);
      //KRATOS_WATCH(tangents_basis_vector)
    }
    return NodeVectorElement;
  }

  std::vector<Point> BrepFace::GetIntersectionPoints(const int& trim_index)
  {
    BrepTrimmingCurve trimming_curve = GetTrimmingCurve(trim_index);
    std::vector<double> intersections = trimming_curve.FindIntersections(m_p, m_q, m_knot_vector_u, m_knot_vector_v);
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

  BrepTrimmingCurve BrepFace::GetTrimmingCurve(const int& trim_index)
  {
    for (unsigned int i = 0; i < m_trimming_loops.size(); i++)
    {
      std::vector<BrepTrimmingCurve> trimming_curves = m_trimming_loops[i].GetTrimmingCurves();
      //std::cout << trimming_curves.size() << std::endl;
      for (unsigned int j = 0; j < trimming_curves.size(); j++)
      {
        if (trimming_curves[j].GetIndex() == trim_index)
        {
          //trimming_curves[j].PrintData();
          return trimming_curves[j];
        }
      }
    }
    KRATOS_ERROR << "Brep Trimming Curve with index " << trim_index << " was not found." << std::endl;
  }

  void BrepFace::EnhanceShapeFunctionsSlave(
    std::vector<Node<3>::Pointer>& nodes, const int& trim_index, const int& shapefunction_order)
  {

    BrepTrimmingCurve trimming_curve = GetTrimmingCurve(trim_index);
    //std::vector<Node<3>::Pointer> NodeVector;
    for (unsigned int i = 0; i < nodes.size(); i++)
    {
      Point point(nodes[i]->X(), nodes[i]->Y(), nodes[i]->Z());
      //double u = 0;
      //double v = 0;
      //GetClosestPoint(point, u, v);

      double u = 0;
      double v = 0;

      double parameter_min = 0.0;
      double parameter_max = 0.0;
      double parameter = 0.0;

      bool success = false;// = NewtonRaphson(point, u, v);
      //std::cout << "Very Initial: u=" << u << ", v=" << v << std::endl;
      //GetClosestPoint(point, u, v);
      if (!success)
      {
        std::vector<array_1d<double, 3>> trimming_curve_loop = trimming_curve.CreatePolygonWithParameter(300);
        parameter_min = trimming_curve_loop[0][2];
        parameter_max = trimming_curve_loop[trimming_curve_loop.size()-1][2];
        double distance = 1e10;
        array_1d<double, 3> initial_guess;
        //for (auto point_itr = trimming_curve_loop.begin(); point_itr != trimming_curve_loop.end(); ++point_itr)
        for (int pt_i = 0; pt_i < trimming_curve_loop.size(); ++pt_i)
        {
          Point point_global;
          EvaluateSurfacePoint(point_global, trimming_curve_loop[pt_i][0], trimming_curve_loop[pt_i][1]);

          //std::ofstream file;
          //file.open("slave_points.txt", std::ios_base::app | std::ios_base::out);
          //file << point_global[0] << " " << point_global[1] << " " << point_global[2] << " " << trimming_curve_loop[pt_i][2] << "\n";
          //file.close();

          //std::ofstream file2;
          //file2.open("master_points.txt", std::ios_base::app | std::ios_base::out);
          //file2 << point[0] << " " << point[1] << " " << point[2] << " " << trimming_curve_loop[pt_i][2] << "\n";
          //file2.close();

          double new_distance = sqrt(((point[0] - point_global[0]) * (point[0] - point_global[0]))
            + ((point[1] - point_global[1]) * (point[1] - point_global[1]))
            + ((point[2] - point_global[2]) * (point[2] - point_global[2])));
          if (new_distance < distance)
          {
            distance = new_distance;
            initial_guess = trimming_curve_loop[pt_i];
            parameter = trimming_curve_loop[pt_i][2];
            if (pt_i > 0)
              parameter_min = trimming_curve_loop[pt_i - 1][2];
            if (pt_i < trimming_curve_loop.size() - 1)
              parameter_max = trimming_curve_loop[pt_i + 1][2];
            else
              parameter_max = trimming_curve_loop[trimming_curve_loop.size() - 1][2];
          }
        }
        u = initial_guess[0];
        v = initial_guess[1];
        //std::cout << "Initial: u=" << u << ", v=" << v <<", distance: " << distance << std::endl;
        bool success = NewtonRaphson(point, u, v);
        //std::cout << "nach newton raphson: u=" << u << ", v=" << v << std::endl;
        if (!success) {
          std::cout << "100 iteration points needed." << std::endl;
          //trimming_curve_loop = trimming_curve.CreatePolygonWithParameter(100);
          //distance = 1e10;
          ////for (auto point_itr = trimming_curve_loop.begin(); point_itr != trimming_curve_loop.end(); ++point_itr)
          //for (int pt_i = 0; pt_i < trimming_curve_loop.size(); ++pt_i)
          //{
          //  Point point_global;
          //  EvaluateSurfacePoint(point_global, trimming_curve_loop[pt_i][0], trimming_curve_loop[pt_i][1]);
          //  double new_distance = sqrt((point[0] - point_global[0]) * (point[0] - point_global[0])
          //    + (point[1] - point_global[1]) * (point[1] - point_global[1])
          //    + (point[2] - point_global[2]) * (point[2] - point_global[2]));
          //  if (new_distance < distance)
          //  {
          //    distance = new_distance;
          //    initial_guess = trimming_curve_loop[pt_i];
          //    if (pt_i > 0)
          //      parameter_min = trimming_curve_loop[pt_i - 1][2];
          //    if (pt_i < trimming_curve_loop.size() - 1)
          //      parameter_max = trimming_curve_loop[pt_i + 1][2];
          //    else
          //      parameter_max = trimming_curve_loop[trimming_curve_loop.size() - 1][2];
          //    parameter = trimming_curve_loop[pt_i][2];
          //  }
          //}
          //u = initial_guess[0];
          //v = initial_guess[1];
          //std::cout << "Initial: u=" << u << ", v=" << v << std::endl;
          //bool success = NewtonRaphson(point, u, v);
          //if (!success)
          //  std::cout << "Newton Raphson did not converge for " << std::endl;
            //KRATOS_ERROR << "Point not found after 100 pts." << std::endl;
        }
      }
      //std::cout << "u=" << u << ", v=" << v << std::endl;
      array_1d<double, 2> point2d(u, v);
      Point point3d;
      //Point point3d;
      std::vector<Vector> location;
      //double parameter_2 = parameter;
      //std::cout << "parameter before: " << parameter << std::endl;
      //std::cout << "point2d before: " << point2d[0] << ", " << point2d[1] << std::endl;
      success = trimming_curve.ProjectionNewtonRaphson(parameter, point2d);
      if (!success)
      {
        success = trimming_curve.ProjectionBisection(parameter, point2d, parameter_min, parameter_max);// , parameter_min, parameter_max);
        //if (!success)
        //  parameter = parameter_2;
      }
      //  success = trimming_curve.ProjectionBisection(parameter, point2d, parameter_min, parameter_max);// , parameter_min, parameter_max);
      //std::cout << "parameter afterwards: " << parameter << "success: " << success <<std::endl;
      //std::cout << "point2d afterwards: " << point2d[0] << ", " << point2d[1] << std::endl;
      trimming_curve.GetCurveDerivatives(location, 0, parameter);
      trimming_curve.EvaluateCurvePoint(point3d, parameter);
      double error_distance = sqrt((point3d[0] - location[0][0]) * (point3d[0] - location[0][0])
        + (point3d[1] - location[0][1]) * (point3d[1] - location[0][1]));
      //std::cout << "error_distance: " << error_distance << std::endl;
      //if (error_distance > 0.00001)
      //{
      //  std::vector<Vector> DN_De;
      //  trimming_curve.EvaluateRationalCurveDerivativesPrint(DN_De, 3, parameter);
      //  Matrix DN_De2;
      //  trimming_curve.EvaluateCurveDerivativesPrint(DN_De2, 3, parameter);
      //}
      //  KRATOS_ERROR << "Big error" << std::endl;
      //std::cout << "point2d afterwards: " << point3d[0] << ", " << point3d[1] << std::endl;

      EvaluateShapeFunctionsSlaveNode(point3d[0], point3d[1], shapefunction_order, nodes[i]);

      array_1d<double, 2> tangents_basis = trimming_curve.GetBaseVector(parameter);
      Vector tangents_basis_vector(2);
      tangents_basis_vector[0] = tangents_basis[0];
      tangents_basis_vector[1] = tangents_basis[1];
      nodes[i]->SetValue(TANGENTS_BASIS_VECTOR_SLAVE, tangents_basis_vector);
    }
  }

  std::vector<Node<3>::Pointer> BrepFace::EnhanceShapeFunctions(
    std::vector<array_1d<double, 3>>& points, const int& shapefunction_order)
  {
    std::vector<Node<3>::Pointer> NodeVector;
    for (unsigned int i = 0; i < points.size(); i++)
    {
      Node<3>::Pointer node = EvaluateNode(points[i][0], points[i][1], shapefunction_order);
      node->SetValue(INTEGRATION_WEIGHT, points[i][2]);
      //KRATOS_WATCH(points[i][2])
      NodeVector.push_back(node);
    }
    return NodeVector;
  }

  void BrepFace::EvaluateShapeFunctionsSlaveNode(const double& u, const double& v, const int& shapefunction_order, Node<3>::Pointer node)
  {
    //Point new_point(0, 0, 0);

    int span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, u);
    int span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, v);

    Vector new_point = ZeroVector(3);

    Vector N = ZeroVector((m_q + 1)*(m_p + 1));

    Vector ControlPointIDs = ZeroVector((m_q + 1)*(m_p + 1));
    Vector local_parameter(2);
    local_parameter[0] = u;
    local_parameter[1] = v;

    Matrix R;
    EvaluateNURBSFunctions(span_u, span_v, u, v, R);

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

        if (shapefunction_order > -1)
          N(k) = R(b, c);
        //N(c + (m_q + 1)*b) = R(b, c);

        //ControlPointIDs(c + (m_q + 1)*b) = m_control_points_ids[control_point_index];
        ControlPointIDs(k) = m_control_points_ids[control_point_index];


        new_point[0] += R(b, c) * m_model_part.GetNode(m_control_points_ids[control_point_index]).X();
        new_point[1] += R(b, c) * m_model_part.GetNode(m_control_points_ids[control_point_index]).Y();
        new_point[2] += R(b, c) * m_model_part.GetNode(m_control_points_ids[control_point_index]).Z();
        //std::cout << "c + (m_q + 1)*b" << c + (m_q + 1)*b << std::endl;
        //std::cout << "k" << k << std::endl;

        k++;
      }
    }

    node->SetValue(LOCATION_SLAVE, new_point);

    if (shapefunction_order > -1)
      node->SetValue(SHAPE_FUNCTION_SLAVE, N);
    //KRATOS_WATCH(ControlPointIDs)
    //KRATOS_WATCH(N)

    if (shapefunction_order > 0)
    {
      Matrix DN_De;
      Matrix DDN_DDe;
      EvaluateNURBSFunctionsDerivatives(span_u, span_v, u, v, DN_De, DDN_DDe);
      //KRATOS_WATCH(DN_De)
      //KRATOS_WATCH(DDN_DDe)
      //  Matrix R = ZeroMatrix(4, 2);
      node->SetValue(SHAPE_FUNCTION_DERIVATIVES_SLAVE, DN_De);

      if (shapefunction_order > 1)
      {
        node->SetValue(SHAPE_FUNCTION_SECOND_DERIVATIVES_SLAVE, DDN_DDe);
      }
    }
    node->SetValue(FACE_BREP_ID_SLAVE, this->Id());
    node->SetValue(LOCAL_PARAMETERS_SLAVE, local_parameter);
    node->SetValue(CONTROL_POINT_IDS_SLAVE, ControlPointIDs);

    //KRATOS_WATCH(local_parameter)
    //KRATOS_WATCH(rSurfacePoint->X())
    //KRATOS_WATCH(rSurfacePoint->Y())
    //KRATOS_WATCH(rSurfacePoint->Z())
    //return rSurfacePoint;
  }

  Node<3>::Pointer BrepFace::EvaluateNode(double u, double v, const int& shapefunction_order)
  {
    Point new_point(0, 0, 0);

    int span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, u);
    int span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, v);

    Vector N = ZeroVector((m_q + 1)*(m_p + 1));

    Vector ControlPointIDs = ZeroVector((m_q + 1)*(m_p + 1));
    Vector local_parameter(2);
    local_parameter[0] = u;
    local_parameter[1] = v;

    Matrix R;
    EvaluateNURBSFunctions(span_u, span_v, u, v, R);
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

        if (shapefunction_order > -1)
          N(k) = R(b, c);
        //N(c + (m_q + 1)*b) = R(b, c);

        ControlPointIDs(k) = m_control_points_ids[control_point_index];
        //ControlPointIDs(c + (m_q + 1)*b) = m_control_points_ids[control_point_index];

        new_point[0] += R(b, c) * m_model_part.GetNode(m_control_points_ids[control_point_index]).X();
        new_point[1] += R(b, c) * m_model_part.GetNode(m_control_points_ids[control_point_index]).Y();
        new_point[2] += R(b, c) * m_model_part.GetNode(m_control_points_ids[control_point_index]).Z();

        //std::cout << "c + (m_q + 1)*b" << c + (m_q + 1)*b << std::endl;
        //std::cout << "k" << k << std::endl;

        k++;
      }
    }

    Node<3>::Pointer rSurfacePoint = Node<3>::Pointer(new Node<3>(0, new_point[0], new_point[1], new_point[2]));

    if (shapefunction_order > -1)
      rSurfacePoint->SetValue(SHAPE_FUNCTION_VALUES, N);
    //KRATOS_WATCH(ControlPointIDs)
    //  KRATOS_WATCH(N)

      if (shapefunction_order > 0)
      {
        Matrix DN_De;
        Matrix DDN_DDe;
        EvaluateNURBSFunctionsDerivatives(span_u, span_v, u, v, DN_De, DDN_DDe);
        //KRATOS_WATCH(DN_De)
        //  KRATOS_WATCH(DDN_DDe)
        //  Matrix R = ZeroMatrix(4, 2);
        rSurfacePoint->SetValue(SHAPE_FUNCTION_DERIVATIVES, DN_De);

        if (shapefunction_order > 1)
        {
          rSurfacePoint->SetValue(SHAPE_FUNCTION_SECOND_DERIVATIVES, DDN_DDe);
        }
      }

    rSurfacePoint->SetValue(FACE_BREP_ID, this->Id());
    rSurfacePoint->SetValue(LOCAL_PARAMETERS, local_parameter);
    rSurfacePoint->SetValue(CONTROL_POINT_IDS, ControlPointIDs);

    //KRATOS_WATCH(local_parameter)
    //  KRATOS_WATCH(rSurfacePoint->X())
    //  KRATOS_WATCH(rSurfacePoint->Y())
    //  KRATOS_WATCH(rSurfacePoint->Z())
    return rSurfacePoint;
  }

  void BrepFace::EnhanceNode(Node<3>::Pointer& node, const double& u, const double& v, const int& shapefunction_order)
  {
    Point new_point(0, 0, 0);

    int span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, u);
    int span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, v);

    Vector N = ZeroVector((m_q + 1)*(m_p + 1));

    Vector ControlPointIDs = ZeroVector((m_q + 1)*(m_p + 1));
    Vector local_parameter(2);
    local_parameter[0] = u;
    local_parameter[1] = v;

    Matrix R;
    EvaluateNURBSFunctions(span_u, span_v, u, v, R);
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

        if (shapefunction_order > -1)
          N(k) = R(b, c);

        ControlPointIDs(k) = m_control_points_ids[control_point_index];

        new_point[0] += R(b, c) * m_model_part.GetNode(m_control_points_ids[control_point_index]).X();
        new_point[1] += R(b, c) * m_model_part.GetNode(m_control_points_ids[control_point_index]).Y();
        new_point[2] += R(b, c) * m_model_part.GetNode(m_control_points_ids[control_point_index]).Z();

        k++;
      }
    }
    node->X() = new_point[0];
    node->Y() = new_point[1];
    node->Z() = new_point[2];

    if (shapefunction_order > -1)
      node->SetValue(SHAPE_FUNCTION_VALUES, N);
    if (shapefunction_order > 0)
    {
      Matrix DN_De;
      Matrix DDN_DDe;
      EvaluateNURBSFunctionsDerivatives(span_u, span_v, u, v, DN_De, DDN_DDe);
      node->SetValue(SHAPE_FUNCTION_DERIVATIVES, DN_De);
      if (shapefunction_order > 1)
        node->SetValue(SHAPE_FUNCTION_SECOND_DERIVATIVES, DDN_DDe);
    }

    node->SetValue(FACE_BREP_ID, this->Id());
    node->SetValue(LOCAL_PARAMETERS, local_parameter);
    node->SetValue(CONTROL_POINT_IDS, ControlPointIDs);
  }

  bool BrepFace::NewtonRaphson(const Point& point, double& u, double& v)
  {
    double norm_delta_u = 100000000;
    //unsigned int k = 0;
    unsigned int max_itr = 20;

    for (int i = 0; i<max_itr; ++i)// (norm_delta_u > 1e-8)
    {
      // newton_raphson_point is evaluated
      Point newton_raphson_point;
      EvaluateSurfacePoint(newton_raphson_point, u, v);

      Vector difference = ZeroVector(3); // Distance between current Q_k and P
                                         // The distance between Q (on the CAD surface) and P (on the FE-mesh) is evaluated
      difference(0) = newton_raphson_point[0] - point[0];
      difference(1) = newton_raphson_point[1] - point[1];
      difference(2) = newton_raphson_point[2] - point[2];

      Matrix hessian = ZeroMatrix(2, 2);
      Vector gradient = ZeroVector(2);
      // The distance is used to compute Hessian and gradient
      EvaluateGradientsForClosestPointSearch(difference, hessian, gradient, u, v);

      double det_H = 0;
      Matrix inv_H = ZeroMatrix(2, 2);

      // u and v are updated
      MathUtils<double>::InvertMatrix(hessian, inv_H, det_H);
      Vector delta_u = prod(inv_H, gradient);
      u -= delta_u(0);
      v -= delta_u(1);

      norm_delta_u = norm_2(delta_u);

      //k++;
      if (norm_delta_u > 1e-7)
        return true;
    }
    return false;
  }

  void BrepFace::GetClosestPoint(const Point& point, double& u, double& v)
  {
    double norm_delta_u = 100000000;
    unsigned int k = 0;
    unsigned int max_itr = 40;

    while (norm_delta_u > 1e-5)
    {
      // newton_raphson_point is evaluated
      Point newton_raphson_point;
      EvaluateSurfacePoint(newton_raphson_point, u, v);

      Vector difference = ZeroVector(3); // Distance between current Q_k and P
                                         // The distance between Q (on the CAD surface) and P (on the FE-mesh) is evaluated
      difference(0) = newton_raphson_point[0] - point[0];
      difference(1) = newton_raphson_point[1] - point[1];
      difference(2) = newton_raphson_point[2] - point[2];

      Matrix hessian = ZeroMatrix(2, 2);
      Vector gradient = ZeroVector(2);
      // The distance is used to compute Hessian and gradient
      EvaluateGradientsForClosestPointSearch(difference, hessian, gradient, u, v);

      double det_H = 0;
      Matrix inv_H = ZeroMatrix(2, 2);

      // u and v are updated
      MathUtils<double>::InvertMatrix(hessian, inv_H, det_H);
      Vector delta_u = prod(inv_H, gradient);
      u -= delta_u(0);
      v -= delta_u(1);

      norm_delta_u = norm_2(delta_u);

      k++;
      if (k>max_itr)
        KRATOS_THROW_ERROR(std::runtime_error, "Newton-Raphson to find closest point did not converge in the following number of iterations: ", k - 1);
    }
  }

  void BrepFace::GetProjectPoint(Node<3>::Pointer& node_on_geometry, const Node<3>::Pointer& node_location, const int& shapefunction_order)
  {
    Point point(node_location->X(), node_location->Y(), node_location->Z());

    Vector local_parameter = node_on_geometry->GetValue(LOCAL_PARAMETERS);
    double u = local_parameter(0);
    double v = local_parameter(1);
    GetClosestPoint(point, u, v);

    EnhanceNode(node_on_geometry, u, v, shapefunction_order);
  }

  // --------------------------------------------------------------------------
  void BrepFace::MapNodeNewtonRaphson(const Node<3>::Pointer& node, Node<3>::Pointer& node_on_geometry)
  {
    std::cout << "test hier" << std::endl;
    // Initialize P: point on the mesh
    Vector P = ZeroVector(3);
    P(0) = node->X();
    P(1) = node->Y();
    P(2) = node->Z();
    // Initialize Q_k: point on the CAD surface
    Vector Q_k = ZeroVector(3);
    Q_k(0) = node_on_geometry->X();
    Q_k(1) = node_on_geometry->Y();
    Q_k(2) = node_on_geometry->Z();
    // Initialize what's needed in the Newton-Raphson iteration				
    Vector Q_minus_P = ZeroVector(3); // Distance between current Q_k and P
    Matrix myHessian = ZeroMatrix(2, 2);
    Vector myGradient = ZeroVector(2);
    double det_H = 0;
    Matrix InvH = ZeroMatrix(2, 2);
    Vector local_parameter = node_on_geometry->GetValue(LOCAL_PARAMETERS);
    double u_k = local_parameter(0);
    double v_k = local_parameter(1);
    //Node<3>::Pointer newtonRaphsonPoint;

    double norm_delta_u = 100000000;
    unsigned int k = 0;
    unsigned int max_itr = 20;
    while (norm_delta_u > 1e-8)
    {
      // The distance between Q (on the CAD surface) and P (on the FE-mesh) is evaluated
      Q_minus_P(0) = Q_k(0) - P(0);
      Q_minus_P(1) = Q_k(1) - P(1);
      Q_minus_P(2) = Q_k(2) - P(2);

      // The distance is used to compute Hessian and gradient
      EvaluateGradientsForClosestPointSearch(Q_minus_P, myHessian, myGradient, u_k, v_k);

      // u_k and v_k are updated
      MathUtils<double>::InvertMatrix(myHessian, InvH, det_H);
      Vector delta_u = prod(InvH, myGradient);
      u_k -= delta_u(0);
      v_k -= delta_u(1);

      // Q is updated
      Point point;
      EvaluateSurfacePoint(point, u_k, v_k);
      Q_k(0) = point[0];
      Q_k(1) = point[1];
      Q_k(2) = point[2];

      KRATOS_WATCH(Q_k)

        //Q_k(0) = newtonRaphsonPoint[0];
        //Q_k(1) = newtonRaphsonPoint[1];
        //Q_k(2) = newtonRaphsonPoint[2];
        norm_delta_u = norm_2(delta_u);

      k++;

      if (k>max_itr)
        KRATOS_THROW_ERROR(std::runtime_error, "Newton-Raphson to find closest point did not converge in the following number of iterations: ", k - 1);
    }
    node_on_geometry->X() = Q_k(0);
    node_on_geometry->Y() = Q_k(1);
    node_on_geometry->Z() = Q_k(2);
  }



  bool BrepFace::CheckIfPointIsInside(Vector node_parameters)
  {
    Polygon polygon(m_trimming_loops);
    return polygon.IsInside(node_parameters[0], node_parameters[1]);
  }
  //// Boost is used to check whether point of interest is inside given polygon or not
  //// Type definitions to use boost functionalities
  ////typedef boost::geometry::model::d2::point_xy<double> point_type;
  ////typedef boost::geometry::model::polygon<point_type> polygon_type;
  //// We assume point is inside, check all boundary loops if this is true. If this is not true for a single loop, then point is considered outside of this patch
  ////point_type point(node_parameters[0], node_parameters[1]);
  //bool is_inside = false;
  //// Loop over all boundary loops of current patch
  //for (unsigned int loop_i = 0; loop_i < m_trimming_loops.size(); loop_i++)
  //{
  //  // Initialize necessary variables
  //  //polygon_type poly;
  //  std::vector<array_1d<double, 2>>& boundary_polygon = m_trimming_loops[loop_i].GetBoundaryPolygon(500);

  //  Polygon polygon(boundary_polygon);

  //  is_inside = polygon.IsInside(node_parameters[0], node_parameters[1]);

  //  // Prepare polygon for boost
  //  //for (unsigned int i = 0; i<boundary_polygon.size(); i++)
  //  //  boost::geometry::append(boost::geometry::exterior_ring(poly),
  //  //    boost::geometry::make<point_type>(boundary_polygon[i][0], boundary_polygon[i][1]));
  //  //if (boundary_polygon.size()>0)
  //  //  boost::geometry::append(boost::geometry::exterior_ring(poly),
  //  //    boost::geometry::make<point_type>(boundary_polygon[0][0], boundary_polygon[0][1]));
  //  // Check inside or outside
  //  // is_inside = boost::geometry::within(point, poly);
  //  // Boost does not consider the polygon direction, it always assumes inside as in the interior of the closed polygon.
  //  // If the CAD loop is an inner loop, however, the area which is considered inner is the unbounded side of the closed polygon.
  //  // So we toggle the results from the within search to be correct.
  //  if (!m_trimming_loops[loop_i].IsOuterLoop())
  //    is_inside = !is_inside;

  //  // If a point is considered outside in one loop, it is outside in general, hence we can break the loop
  //  if (!is_inside)
  //    break;
  //}
  //return is_inside;
  //}

  //GEOMETRY FUNCTIONS:

  /**
  * @Author Daniel Baumgaertner
  * @date   December, 2016
  * @brief   returns the cartesian coordinates (global) for a specific point
  * located on the NURBS surface S(u=fixed and v=fixed)
  * Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
  * Algorithm A4.3
  *
  * @param[in]  rSurfacePoint  evaluated point
  * @param[in]  u  local parameter in u-direction
  * @param[in]  v  local parameter in v-direction
  */
  void BrepFace::EvaluateSurfacePoint(Point& rSurfacePoint, const double& u, const double& v)
  {
    //Point new_point(0, 0, 0);
    rSurfacePoint[0] = 0;// (new_point[0], new_point[1], new_point[2]);
    rSurfacePoint[1] = 0;
    rSurfacePoint[2] = 0;

    int span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, u);
    int span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, v);

    Vector ShapeFunctionsN = ZeroVector((m_q + 1)*(m_p + 1));
    Matrix N;
    EvaluateNURBSFunctions(span_u, span_v, u, v, N);

    for (int c = 0; c <= m_q; c++)
    {
      for (int b = 0; b <= m_p; b++)
      {
        // the control point vector is filled up by first going over u, then over v
        int ui = span_u - m_p + b;
        int vi = span_v - m_q + c;
        int m_n_u = m_knot_vector_u.size() - m_p - 1;
        int control_point_index = vi*m_n_u + ui;

        rSurfacePoint[0] += N(b, c) * m_model_part.GetNode(m_control_points_ids[control_point_index]).X();
        rSurfacePoint[1] += N(b, c) * m_model_part.GetNode(m_control_points_ids[control_point_index]).Y();
        rSurfacePoint[2] += N(b, c) * m_model_part.GetNode(m_control_points_ids[control_point_index]).Z();
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

        double cp_x = m_model_part.GetNode(m_control_points_ids[control_point_index]).X();
        double cp_y = m_model_part.GetNode(m_control_points_ids[control_point_index]).Y();
        double cp_z = m_model_part.GetNode(m_control_points_ids[control_point_index]).Z();

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

    R.resize(m_p + 1, m_q + 1);
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
        int m_n_u = m_knot_vector_u.size() - m_p - 1;
        int control_point_index = vi*m_n_u + ui;
        // Evaluate basis function
        R(b, c) = N(b)*M(c)*m_model_part.GetNode(m_control_points_ids[control_point_index]).GetValue(CONTROL_POINT_WEIGHT);
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
    DN_De.resize(number_of_control_points, 2);
    DDN_DDe.resize(number_of_control_points, 3);

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

        // Evaluate basis function
        weight = m_model_part.GetNode(m_control_points_ids[control_point_index]).GetValue(CONTROL_POINT_WEIGHT);

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

        // Evaluate basis function
        weight = m_model_part.GetNode(m_control_points_ids[control_point_index]).GetValue(CONTROL_POINT_WEIGHT);
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
    TrimmingLoopVector& trimming_loops,
    TrimmingLoopVector& embedded_loops,
    Vector& knot_vector_u, Vector& knot_vector_v,
    unsigned int& p, unsigned int& q, IntVector& control_point_ids,
    ModelPart& model_part)
    : m_trimming_loops(trimming_loops),
      m_embedded_loops(embedded_loops),
      m_knot_vector_u(knot_vector_u),
      m_knot_vector_v(knot_vector_v),
      m_p(p),
      m_q(q),
      m_model_part(model_part),
      m_control_points_ids(control_point_ids),
      IndexedObject(brep_id),
      Flags()
  {
  }
  ///Destructor
  BrepFace::~BrepFace()
  {}

} // namespace Kratos.