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


// Project includes
#include "BrepTrimmingCurve.h"
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"


namespace Kratos
{
// --------------------------------------------------------------------------
  std::vector<array_1d<double, 2>> BrepTrimmingCurve::CreatePolygon(unsigned int number_polygon_points)
  {
    std::vector<array_1d<double, 2>> polygon;
    polygon.resize(number_polygon_points);

    // Variables needed
    unsigned int counter = 0;

    double u_min = m_knot_vector_u[0];
    double u_max = m_knot_vector_u[m_knot_vector_u.size() - 1];
    double delta_u = (u_max - u_min) / (number_polygon_points-1);

    // Add points of edge to polygon
    double u_i = u_min - delta_u;
    for (unsigned int i = 0; i<number_polygon_points; i++)
    {
      u_i += delta_u;
      Point<3> curve_point;

      EvaluateCurvePoint(curve_point, u_i);

      polygon[counter][0] = curve_point.X();
      polygon[counter][1] = curve_point.Y();

      counter++;
    }
    return polygon;
  }

  std::vector<array_1d<double, 3>> BrepTrimmingCurve::CreatePolygonWithParameter(unsigned int number_polygon_points)
  {
    std::vector<array_1d<double, 3>> polygon;
    polygon.resize(number_polygon_points);

    // Variables needed
    unsigned int counter = 0;

    //KRATOS_WATCH(m_knot_vector_u)

    double u_min = m_knot_vector_u(0);
    double u_max = m_knot_vector_u(m_knot_vector_u.size() - 1);
    double delta_u = (u_max - u_min) / number_polygon_points;

    //std::cout << "u_min: " << u_min << std::endl;
    // Add points of edge to polygon
    double u_i = u_min - delta_u;
    for (unsigned int i = 0; i<number_polygon_points; i++)
    {
      //std::cout << "i: " << i << std::endl;
      u_i += delta_u;
      Point<3> curve_point;

      //std::cout << "u_i: " << u_i << std::endl;
      EvaluateCurvePoint(curve_point, u_i);

      polygon[counter][0] = curve_point.X();
      polygon[counter][1] = curve_point.Y();
      polygon[counter][2] = u_i;

      counter++;
    }

    //std::cout << "polygon finished" << std::endl;
    return polygon;
  }

  void BrepTrimmingCurve::EvaluateCurvePoint(Point<3>& rCurvePoint, double parameter_u)
  {
    const unsigned int R_Dim = 3;

    //Vector tmpHomCtrlPt;
    Vector N;
    Vector resulting_point(R_Dim+1);
    Vector homPoi(R_Dim+1);
    Matrix homCtrlPts;

    unsigned int span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, parameter_u);
    NurbsUtilities::eval_nonzero_basis_function(N, m_knot_vector_u, parameter_u, span_u, m_p);
    homCtrlPts.resize(R_Dim + 1, m_p + 1);
    for (unsigned int i = 0; i <= m_p; i++)
    {
      int control_point_index = span_u - m_p + i;
      //std::cout << control_point_index << std::endl;
      //std::cout << m_control_points.size() << std::endl;
      // tmpHomCtrlPt = m_control_points[control_point_index].getWeight();
      // tmpHomCtrlPt = Ctrl_Pt[span-m_p+i]->get_Ctrl_Pt_Coo_w();
      Vector tmpHomCtrlPt(R_Dim + 1);

      double u = m_control_points[control_point_index][0];
      double v = m_control_points[control_point_index][1];
      double w = m_control_points[control_point_index][2];
      double weight = m_control_points[control_point_index][3];

      tmpHomCtrlPt[0] = u*weight;
      tmpHomCtrlPt[1] = v*weight;
      tmpHomCtrlPt[2] = w*weight;
      tmpHomCtrlPt[3] = weight;

      for (unsigned int j = 0; j <= R_Dim; j++)
      {
        homCtrlPts(j, i) = tmpHomCtrlPt(j);
      }
    }
    homPoi = prod(homCtrlPts, N);
    resulting_point = (1 / homPoi(R_Dim))*homPoi;
    //std::cout << "test jolo2.1.3.6" << std::endl;

    rCurvePoint[0] = resulting_point[0];
    rCurvePoint[1] = resulting_point[1];
    rCurvePoint[2] = resulting_point[2];

    if (std::abs(resulting_point(R_Dim) - 1.00) > 0.000001)
    {
      KRATOS_THROW_ERROR(std::logic_error, "NURBS 1D: evalutation curve point failed!!!", "");
    }
  }

  std::vector<array_1d<double, 3>> BrepTrimmingCurve::GetQuadraturePoints(std::vector<double> span, double polynomial_order_p)
  {
    std::vector<array_1d<double, 3>> quadrature_points;
    for (unsigned int i = 0; i < span.size()-1; i++)
    {
      Vector parameter_u(2);
      parameter_u[0] = span[i];
      parameter_u[1] = span[i + 1];
      KnotSpan1d knot_span(0, polynomial_order_p, parameter_u);

      std::vector<array_1d<double, 2>> integration_points = knot_span.getIntegrationPointsInParameterDomain();
      for (unsigned int j = 0; j < integration_points.size(); j++)
      {
        Point<3> integration_point;
        EvaluateCurvePoint(integration_point, integration_points[j][0]);

        array_1d<double, 3> quadrature_point;
        quadrature_point[0] = integration_point[0];
        quadrature_point[1] = integration_point[1]; 
        quadrature_point[2] = integration_points[j][1];

        quadrature_points.push_back(quadrature_point);
      }
    }
    return quadrature_points;
  }

  std::vector<double> BrepTrimmingCurve::FindIntersections(const int& p, const int& q, const Vector& knot_vector_u, const Vector& knot_vector_v)
  {
    //KRATOS_WATCH(m_active_range)
    //KRATOS_WATCH(m_knot_vector_u)

    std::vector<array_1d<double, 3>> trim_polygon = CreatePolygonWithParameter(100);
    std::vector<double> intersections;

    std::cout << "trim_polygon: " << trim_polygon.size() << std::endl;

    intersections.push_back(m_active_range[0]);

    int span_u = NurbsUtilities::find_knot_span(p, knot_vector_u, trim_polygon[0][0]);
    int span_v = NurbsUtilities::find_knot_span(q, knot_vector_v, trim_polygon[0][1]);

    for (unsigned int i = 0; i < trim_polygon.size(); i++)
    {
      //std::cout << "polygon point: " << i << std::endl;

      int new_span_u = NurbsUtilities::find_knot_span(p, knot_vector_u, trim_polygon[i][0]);
      if (span_u < new_span_u)
      {
        //std::cout << "change in span u: " << knot_vector_u[new_span_u] << std::endl;
        int intersection_base = 1;
        intersections.push_back(EvaluateIntersection(trim_polygon[i - 1][2], intersection_base, knot_vector_u[new_span_u]));
        span_u = new_span_u;
      }
      if (span_u > new_span_u)
      {
        //std::cout << "change in span u: " << knot_vector_u[new_span_u+1] << std::endl;
        int intersection_base = 1;
        intersections.push_back(EvaluateIntersection(trim_polygon[i - 1][2], intersection_base, knot_vector_u[new_span_u+1]));
        span_u = new_span_u;
      }

      int new_span_v = NurbsUtilities::find_knot_span(q, knot_vector_v, trim_polygon[i][1]);
      if (span_v < new_span_v)
      {
        //std::cout << "change in span v: " << knot_vector_v[new_span_v]  << std::endl;
        int intersection_base = 2;
        intersections.push_back(EvaluateIntersection(trim_polygon[i-1][2], intersection_base, knot_vector_v[new_span_v]));
        span_v = new_span_v;
      }
      if (span_v > new_span_v)
      {
        //std::cout << "change in span v: " << knot_vector_v[new_span_v+1] << std::endl;
        int intersection_base = 2;
        intersections.push_back(EvaluateIntersection(trim_polygon[i - 1][2], intersection_base, knot_vector_v[new_span_v+1]));
        span_v = new_span_v;
      }
    }
    //std::cout << "test jolo2.9" << std::endl;
    intersections.push_back(m_active_range[1]);
    std::vector<double> full_intersections;
    full_intersections.push_back(intersections[0]);
    for (unsigned int j = 1; j < intersections.size(); j++)
    {
      if (intersections[j - 1] != intersections[j])
      {
        full_intersections.push_back(intersections[j]);
        //KRATOS_WATCH(intersections[j])
      }
    }
    //KRATOS_WATCH(full_intersections)
    return full_intersections;
  }

  std::vector<double> BrepTrimmingCurve::FindIntersectionsWithPoints(std::vector<Point<2>> intersection_points)
  {
    std::vector<double> intersections;
    for (unsigned int i = 0; i < intersection_points.size(); i++)
    {
      double parameter = 0;
      GetClosestPoint(intersection_points[i], parameter);
      intersections.push_back(parameter);
    }
    return intersections;
  }

  double BrepTrimmingCurve::EvaluateIntersection(double initial_u, int intersection_base, const double& coordinate_base)
  {
    double local_parameter_u;
    bool converged;
    EvaluateLocalParameter(local_parameter_u, converged, intersection_base, coordinate_base, initial_u, 20, 1e-8);
    return local_parameter_u;
  }

  //  #####################################################################################
  // #######################################################################################
  //#
  //#                  ++++++++++++++++++++++++++++++++++++++++++++++
  //#                  +++  Nurbs1DBasis::eval_Local_Parameter_Ex +++
  //#                  ++++++++++++++++++++++++++++++++++++++++++++++
  //#
  ///   \details     returns the local parameter based one component of the position vector
  ///                  (a*e11, b*e22, c*e33)
  ///                  nonlinear equation is solved by means of a Newton Raphson scheme
  ///
  /// ======================================================================================
  ///   \param[in]   _parameter      local parameter
  ///   \param[in]   _converged      convergence flag
  ///   \param[in]   _baseVec      base vector (1->e11, 2->e22, 3->e33)
  ///   \param[in]   _baseComp        coordinate of the corresponding base vector(a, b, c)
  ///   \param[in]  _uinit        Newton Raphson: inital guess
  ///   \param[in]  _itmax        Newton Raphson: maximal number of iterations
  ///   \param[in]  _iteps        Newton Raphson: tolerance
  ///
  /// ======================================================================================
  ///  \author     A.Widhammer (11/2011)
  //
  //########################################################################################
  void BrepTrimmingCurve::EvaluateLocalParameter(double& parameter, bool& converged, int baseVec,
    double baseComp, double uinit, double itmax, double iteps)
  {
    std::cout << "Evaluate Local Parameter, u_initial: " << uinit << std::endl;
    // local parameters
    int i;
    double u_n;
    double u_np1;
    double dynR;
    double refR;
    double relR;
    double J;
    //array_1d<double, 2> bounds;
    Point<3> point_1, point_2;
    //Point<3> poi_max;
    Matrix dCdu;

    double Tolerance = 1e-9;

    if (fabs(baseComp) < 1e-9)
      baseComp = 0.0;

    // check span
    //this->min_Max_Knot(bounds);
    this->EvaluateCurvePoint(point_1, m_knot_vector_u(0));
    this->EvaluateCurvePoint(point_2, m_knot_vector_u(m_knot_vector_u.size()-1));

    double point_min = point_1(baseVec - 1);
    double parameter_min = m_knot_vector_u(0);
    double point_max = point_2(baseVec - 1);
    double parameter_max = m_knot_vector_u(m_knot_vector_u.size() - 1);

    if (point_min > point_max)
    {
      point_min = point_2(baseVec - 1);
      parameter_min = m_knot_vector_u(m_knot_vector_u.size() - 1);
      point_max = point_1(baseVec - 1);
      parameter_max = m_knot_vector_u(0);
    }

    if (baseComp<point_min)
    {
      parameter = parameter_min;
      converged = false;
    }
    else if (baseComp>point_max)
    {
      parameter = parameter_max;
      converged = false;
    }
    else
    {
      u_n = uinit;
      // START: Newton-Raphson
      for (i = 0; i <= itmax; i++)
      {
        // establish residuum
        Point<3> point;
        this->EvaluateCurvePoint(point, u_n);
        array_1d<double, 2> base_vectors = this->GetBaseVector(u_n);//EvaluateCurveDerivatives(dCdu, 1, u_n);
        dynR = point(baseVec - 1) - baseComp;
        KRATOS_WATCH(dynR)
        J = base_vectors(baseVec - 1);
        KRATOS_WATCH(J)
        // check convergence - establish reference residuum
        if (i == 0)
        {
          refR = fabs(dynR);
          if (refR < 10 * Tolerance)
          {
            converged = true;
            parameter = u_n;
            break;
          }
        }
        relR = fabs(dynR);
        if ((relR < iteps) && (i != 0))
        {
          converged = true;
          parameter = u_n;
          break;
        }
        // compute new iteration step
        u_np1 = u_n - dynR / J;
        // update solution
        if (u_np1 < point_min)
        {
          u_n = parameter_min;
        }
        else if (u_np1 > point_max)
        {
          u_n = parameter_max;
        }
        else
        {
          u_n = u_np1;
        }
      }
      // END: Newton-Raphson
      if (i == itmax)
      {
        parameter = 0.00;
        converged = false;
      }
    }
  }

  void BrepTrimmingCurve::GetClosestPoint(const Point<2>& closest_point, double& parameter)
  {
    //std::cout << "Evaluate Local Parameter, u_initial: " << uinit << std::endl;
    //// local parameters
    //int i;
    //double u_n;
    //double u_np1;
    //double dynR;
    //double refR;
    //double relR;
    //double J;
    ////array_1d<double, 2> bounds;
    //Point<3> point_1, point_2;
    ////Point<3> poi_max;
    //Matrix dCdu;

    double ModelTolerance = 1e-3;

    //if (fabs(baseComp) < 1e-9)
    //  baseComp = 0.0;

    //// check span
    ////this->min_Max_Knot(bounds);
    //this->EvaluateCurvePoint(point_1, m_knot_vector_u(0));
    //this->EvaluateCurvePoint(point_2, m_knot_vector_u(m_knot_vector_u.size() - 1));

    //double point_min = point_1(baseVec - 1);
    //double parameter_min = m_knot_vector_u(0);
    //double point_max = point_2(baseVec - 1);
    //double parameter_max = m_knot_vector_u(m_knot_vector_u.size() - 1);

    //if (point_min > point_max)
    //{
    //  point_min = point_2(baseVec - 1);
    //  parameter_min = m_knot_vector_u(m_knot_vector_u.size() - 1);
    //  point_max = point_1(baseVec - 1);
    //  parameter_max = m_knot_vector_u(0);
    //}

    //if (baseComp<point_min)
    //{
    //  parameter = parameter_min;
    //  converged = false;
    //}
    //else if (baseComp>point_max)
    //{
    //  parameter = parameter_max;
    //  converged = false;
    //}
    //else
    //{
      //u_n = param;
    Vector Q = ZeroVector(2);
    int itmax = 20;
      // START: Newton-Raphson
      for (unsigned int i = 0; i <= itmax; i++)
      {
        // establish residuum
        Point<3> point;
        this->EvaluateCurvePoint(point, parameter);
        array_1d<double, 2> base_vector = this->GetBaseVector(parameter);//EvaluateCurveDerivatives(dCdu, 1, u_n);
        Q(0) = point[0] - closest_point[0];
        Q(1) = point[1] - closest_point[1];
        KRATOS_WATCH(Q)
        //J = base_vectors(baseVec - 1);
        KRATOS_WATCH(base_vector)
          // check convergence - establish reference residuum
          //if (i == 0)
          //{
            double ref_distance = sqrt(Q(0)*Q(0)+Q(1)*Q(1));
            if (ref_distance < 10 * ModelTolerance)
            {
              //converged = true;
              //parameter = u_n;
              break;
            }
          //}
        //relR = fabs(dynR);
        //if ((relR < iteps) && (i != 0))
        //{
        //  converged = true;
        //  parameter = u_n;
        //  break;
        //}
        // compute new iteration step
        parameter = parameter - Q(0)*base_vector(0) - Q(1)*base_vector(1);
        // update solution
      //  if (u_np1 < point_min)
      //  {
      //    u_n = parameter_min;
      //  }
      //  else if (u_np1 > point_max)
      //  {
      //    u_n = parameter_max;
      //  }
      //  else
      //  {
      //    u_n = u_np1;
      //  }
      //}
      // END: Newton-Raphson
      if (i == itmax)
      {
        parameter = 0.00;
        KRATOS_THROW_ERROR(std::logic_error, "BrepTrimmingCurve: Newton Raphson failed!", "");
        //converged = false;
      }
    }
  }

  array_1d<double, 2> BrepTrimmingCurve::GetBaseVector(const int& u)
  {
    int span = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, u);

    Matrix N;
    NurbsUtilities::eval_nonzero_basis_function_with_derivatives(N, m_knot_vector_u, u, span, m_p, 3);

    //this->eval_Derivative_NonzeroBasis_Fct(N, U_Vec, _u, i, P_Deg, 1);
    double sum = 0.0;
    double dsum = 0.0;
    Vector R(m_p + 1);
    Vector dR(m_p + 1);
    double weight;

    for (int b = 0; b <= m_p; b++)
    {
      weight = m_control_points[span - m_p + b][3];
      R(b) = N(0, b)*weight;
      sum += R(b);

      // derivatives
      dR(b) = N(1, b)*weight;
      dsum += dR(b);
    }

    // divide by sum only required in terms of rational basis functions
    if (fabs(sum - weight)> 0.00000001)
    {
      double inv_sum = 1 / sum;
      double inv_sum_2 = 1 / pow(sum, 2);
      // divide through by sum
      for (int b = 0; b <= m_p; b++)
      {
        dR(b) = inv_sum*dR(b) - R(b)*dsum*inv_sum_2;
        R(b) = inv_sum*R(b);
      }
    }

    array_1d<double, 2> g;
    g.clear();


    for (int b = 0; b <= m_p; b++)
    {

      g(0) += dR(b)*m_control_points[span - m_p + b][0];
      g(1) += dR(b)*m_control_points[span - m_p + b][1];
      //g(2) += dR(b)*Ctrl_Pt[i-P_Deg+b]->get_Coo_e33();

    }
    return g;
  }

  void BrepTrimmingCurve::EvaluateCurveDerivatives(Matrix& DN_De, const int& order, const int& u)
  {
    //DN_De = ZeroMatrix(m_control_points.size(), order);

    //for (unsigned int k = 0; k <= order; k++)
    //{
    //  for (unsigned int n_cp_itr = 0; n_cp_itr <= DN_De.size1(); n_cp_itr++)
    //  {
    //    DN_De(n_cp_itr, k) = m_control_points[n_cp_itr][k];
    //  }
    //  for (unsigned int i = 1; i <= k; i++)
    //  {
    //    for (unsigned int n_cp_itr = 0; n_cp_itr <= DN_De.size1(); n_cp_itr++)
    //    {
    //      DN_De(n_cp_itr, k) = DN_De(n_cp_itr, k) - NurbsUtilities::binom(k,i)*DN_De(n_cp_itr,k-i);
    //    }
    //  }
    //}


    DN_De.resize(4, m_control_points.size());
    matrix<double> N;                        //B-Spline basis functions
    int span = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, u);
    NurbsUtilities::eval_nonzero_basis_function_with_derivatives(N, m_knot_vector_u, u, span, m_p, 3);
    
    double sum = 0;    //sum of weights times basis functions
    double dsum = 0;   //sum of weights times 1st derivative of basis functions
    double ddsum = 0;  //sum of weights times 2nd derivative of basis functions
    double dddsum = 0; //sum of weights times 3rd derivative of basis functions

    for (int i = 0; i <= m_p; i++)
    {
      double weight = m_control_points[i + span - m_p][3];

      DN_De(0, i) = N(0, i)*weight;
      sum += DN_De(0, i);
      DN_De(1, i) = N(1, i)*weight;
      dsum += DN_De(1, i);
      DN_De(2, i) = N(2, i)*weight;
      ddsum += DN_De(2, i);
      DN_De(3, i) = N(3, i)*weight;
      dddsum += DN_De(3, i);
    }

    //get derivatives
    for (int i = 0; i <= m_p; i++)
    {
      //3rd derivative
      DN_De(3, i) = DN_De(3, i) / sum - 3 * DN_De(2, i)*dsum / pow(sum, 2) + 4 * DN_De(1, i)*pow(dsum, 2) / pow(sum, 3)
        - (3 * DN_De(1, i)*ddsum + DN_De(0, i)*dddsum) / pow(sum, 2) + 2 * DN_De(0, i)*ddsum*dsum / pow(sum, 3)
        + (2 * DN_De(1, i)*pow(dsum, 2) + 4 * DN_De(0, i)*dsum*ddsum) / pow(sum, 3)
        - 6 * DN_De(0, i)*pow(dsum, 3) / pow(sum, 4);


      //DN_De(3,i) = DN_De(3,i)/sum - 3*DN_De(2,i)*dsum/pow(sum,2) - 3*DN_De(1,i)*ddsum/pow(sum,2) + 4*DN_De(1,i)*dsum*ddsum/pow(sum,3)
      //         - DN_De(0,i)*dddsum/pow(sum,2) + 2*DN_De(1,i)*pow(dsum,2)/pow(sum,3) + 6*DN_De(0,i)*dsum*ddsum/pow(sum,3) - 6*pow(dsum,3)/pow(sum,4);
      //2nd derivative
      DN_De(2, i) = DN_De(2, i) / sum - 2 * DN_De(1, i)*dsum / pow(sum, 2) - DN_De(0, i)*ddsum / pow(sum, 2) + 2 * DN_De(0, i)*pow(dsum, 2) / pow(sum, 3);

      //1st derivative
      DN_De(1, i) = DN_De(1, i) / sum - DN_De(0, i)*dsum / pow(sum, 2);

      //basis functions
      DN_De(0, i) = DN_De(0, i) / sum;
    }
  }


  unsigned int& BrepTrimmingCurve::GetIndex()
  {
    return m_trim_index;
  }
    //BrepTrimmingCurve::BrepTrimmingCurve(BrepTrimmingCurve const& rOther)
    //{}
  /// Print object's data.
  void BrepTrimmingCurve::PrintData()
  {
      KRATOS_WATCH(m_knot_vector_u)
      KRATOS_WATCH(m_control_points.size())
      KRATOS_WATCH(m_control_points[0])
      KRATOS_WATCH(m_control_points[1])
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
    //KRATOS_WATCH(knot_vector_u)
    //KRATOS_WATCH(m_knot_vector_u)
    //KRATOS_WATCH(m_control_points.size())
    //  KRATOS_WATCH(m_control_points[0])
    //  KRATOS_WATCH(m_control_points[1])
  }
//Destructor
BrepTrimmingCurve::~BrepTrimmingCurve()
{}

}  // namespace Kratos.

