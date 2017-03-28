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
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"


namespace Kratos
{
  // --------------------------------------------------------------------------
  Vector& BrepFace::GetUKnotVector()
  {
    return m_knot_vector_u;
  }
  Vector& BrepFace::GetVKnotVector()
  {
    return m_knot_vector_u;
  }
  // --------------------------------------------------------------------------
  void BrepFace::MapNodeNewtonRaphson(const Node<3>::Pointer& node, Node<3>::Pointer& node_on_geometry, ModelPart& model_part)
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
    double v_k = local_parameter(0);
    //Node<3>::Pointer newtonRaphsonPoint;

    double norm_deltau = 100000000;
    unsigned int k = 0;
    unsigned int max_itr = 20;
    while (norm_deltau > 1e-8)
    {
      // The distance between Q (on the CAD surface) and P (on the FE-mesh) is evaluated
      Q_minus_P(0) = Q_k(0) - P(0);
      Q_minus_P(1) = Q_k(1) - P(1);
      Q_minus_P(2) = Q_k(2) - P(2);

      // The distance is used to compute Hessian and gradient
      EvaluateGradientsForClosestPointSearch(Q_minus_P, myHessian, myGradient, u_k, v_k, model_part);

      // u_k and v_k are updated
      MathUtils<double>::InvertMatrix(myHessian, InvH, det_H);
      Vector deltau = prod(InvH, myGradient);
      u_k -= deltau(0);
      v_k -= deltau(1);

      // Q is updated
      EvaluateSurfacePoint(node_on_geometry, u_k, v_k, model_part);
      Q_k(0) = node_on_geometry->X();
      Q_k(1) = node_on_geometry->Y();
      Q_k(2) = node_on_geometry->Z();

      KRATOS_WATCH(Q_k)

      //Q_k(0) = newtonRaphsonPoint[0];
      //Q_k(1) = newtonRaphsonPoint[1];
      //Q_k(2) = newtonRaphsonPoint[2];
      norm_deltau = norm_2(deltau);

      k++;

      if (k>max_itr)
        KRATOS_THROW_ERROR(std::runtime_error, "Newton-Raphson to find closest point did not converge in the following number of iterations: ", k - 1);
    }
    // Update nearest point
    //local_parameter[0] = u_k;
    //local_parameter[1] = v_k;
    //node_on_geometry->X() = Q_k(0);
    //node_on_geometry->Y() = Q_k(1);
    //node_on_geometry->Z() = Q_k(2);
    //node_on_geometry->SetValue(LOCAL_PARAMETERS, local_parameter);

    //// Compute and store span of each parameter to avoid redundant computations later
    //IntVector knot_span_nearest_point = m_patches[patch_itr_of_nearest_point].GetSurface().GetKnotSpan(u_of_nearest_point, v_of_nearest_point);

    //// Set flag to mark control point as relevant for mapping
    //int span_u_of_np = knot_span_nearest_point[0];
    //int span_v_of_np = knot_span_nearest_point[1];
    //m_patches[patch_itr_of_nearest_point].GetSurface().FlagControlPointsForMapping(span_u_of_np, span_v_of_np, u_of_nearest_point, v_of_nearest_point);
  }

  bool BrepFace::CheckIfPointIsInside(Vector node_parameters)
  {
    // Boost is used to check whether point of interest is inside given polygon or not
    // Type definitions to use boost functionalities
    typedef boost::geometry::model::d2::point_xy<double> point_type;
    typedef boost::geometry::model::polygon<point_type> polygon_type;
    // We assume point is inside, check all boundary loops if this is true. If this is not true for a single loop, then point is considered outside of this patch
    point_type point(node_parameters[0], node_parameters[1]);
    bool is_inside = true;
    // Loop over all boundary loops of current patch
    for (unsigned int loop_i = 0; loop_i < m_trimming_loops.size(); loop_i++)
    {
      // Initialize necessary variables
      polygon_type poly;
      std::vector<array_1d<double, 2>>& boundary_polygon = m_trimming_loops[loop_i].GetBoundaryPolygon();

      // Prepare polygon for boost
      for (unsigned int i = 0; i<boundary_polygon.size(); i++)
        boost::geometry::append(boost::geometry::exterior_ring(poly),
          boost::geometry::make<point_type>(boundary_polygon[i][0], boundary_polygon[i][1]));
      if (boundary_polygon.size()>0)
        boost::geometry::append(boost::geometry::exterior_ring(poly),
          boost::geometry::make<point_type>(boundary_polygon[0][0], boundary_polygon[0][1]));
      // Check inside or outside
      is_inside = boost::geometry::within(point, poly);
      // Boost does not consider the polygon direction, it always assumes inside as in the interior of the closed polygon.
      // If the CAD loop is an inner loop, however, the area which is considered inner is the unbounded side of the closed polygon.
      // So we toggle the results from the within search to be correct.
      if (!m_trimming_loops[loop_i].IsOuterLoop())
        is_inside = !is_inside;

      // If a point is considered outside in one loop, it is outside in general, hence we can break the loop
      if (!is_inside)
        break;
    }
    return is_inside;
  }

  //GEOMETRY FUNCTIONS:
  //  #####################################################################################
  // #######################################################################################
  ///
  ///  \details    returns the cartesian coordinates (global) for a specific point 
  ///              located on the NURBS surface S(u=fixed and v=fixed)
  ///              Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
  ///              Algorithm A4.3
  ///
  /// ======================================================================================
  ///  \param[in]  rSurfacePoint   	evaluated point
  ///  \param[in]  _uoi    			local parameter in u-direction
  ///  \param[in]  _voi    			local parameter in v-direction
  ///
  /// ======================================================================================
  ///  \author     Daniel Baumgärtner (12/2016)
  //
  //########################################################################################
  void BrepFace::EvaluateSurfacePoint(Node<3>::Pointer& rSurfacePoint, double u, double v, ModelPart& model_part)
  {
    Point<3> new_point(0, 0, 0);
    //rSurfacePoint->X() = 0;
    //rSurfacePoint->Y() = 0;
    //rSurfacePoint->Z() = 0;

    int span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, u);
    int span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, v);

    Vector ShapeFunctionsN = ZeroVector((m_q + 1)*(m_p + 1));
    //std::vector<int> ControlPointIDs;
    //ControlPointIDs.resize((m_q + 1)*(m_p + 1));
    Vector local_parameter(2);
    local_parameter[0] = u;
    local_parameter[1] = v;

    for (int c = 0; c <= m_q; c++)
    {
      for (int b = 0; b <= m_p; b++)
      {
        // the control point vector is filled up by first going over u, then over v
        int ui = span_u - m_p + b;
        int vi = span_v - m_q + c;
        int m_n_u = m_knot_vector_u.size() - m_p - 1;
        int control_point_index = vi*m_n_u + ui;

        Matrix R;
        EvaluateNURBSFunctions(span_u, span_v, u, v, R, model_part);

        ShapeFunctionsN(b + (m_q + 1)*c) = R(b, c);
        //ControlPointIDs.push_back(m_control_points_ids[control_point_index]);
        //ControlPointIDs[b + (m_q + 1)*c] = m_control_points_ids[control_point_index];

        new_point[0] += R(b, c) * model_part.GetNode(m_control_points_ids[control_point_index]).X();
        new_point[1] += R(b, c) * model_part.GetNode(m_control_points_ids[control_point_index]).Y();
        new_point[2] += R(b, c) * model_part.GetNode(m_control_points_ids[control_point_index]).Z();
      }
    }
    rSurfacePoint->X() = new_point[0];// (new_point[0], new_point[1], new_point[2]);
    rSurfacePoint->Y() = new_point[1];
    rSurfacePoint->Z() = new_point[2];

    KRATOS_WATCH(new_point)

    rSurfacePoint->SetValue(LOCAL_PARAMETERS, local_parameter);
    rSurfacePoint->SetValue(SHAPE_FUNCTION_VALUES, ShapeFunctionsN);
    //rSurfacePoint->SetValue(CONTROL_POINT_IDS, ControlPointIDs);
  }
  // #######################################################################################
  ///
  ///  \details    evaluate Hessian and Gradiend modifying the input objects
  ///
  /// ======================================================================================
  ///  \param[in]  QminP    	 	Distance Vector
  ///  \param[in]  H		     	Hessian reference	
  ///  \param[in]  Gradient    	Gradient reference
  ///  \param[in]  v    			parameter
  ///  \param[in]  u 				parameter 
  ///
  /// ======================================================================================
  ///  \author     Giovanni Filomeno (1/2017) && Massimo Sferza (1/2017)
  //
  //########################################################################################	

  void BrepFace::EvaluateGradientsForClosestPointSearch(Vector QminP, Matrix& Hessian, Vector& Gradient, double& u, double& v, ModelPart& model_part)
  {
    // The derivatives of the basis functions are evaluated
    Matrix dR;
    Matrix ddR;
    EvaluateNURBSFunctionsDerivatives(-1, -1, u, v, dR, ddR, model_part);

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

        double cp_x = model_part.GetNode(m_control_points_ids[control_point_index]).X();
        double cp_y = model_part.GetNode(m_control_points_ids[control_point_index]).Y();
        double cp_z = model_part.GetNode(m_control_points_ids[control_point_index]).Z();

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
  ///  \details    returns the basis functions of NURBS basis function w.r.t. u,v
  ///              span_u, span_v are the knot span indices. if unknown, insert 0!
  ///
  /// ======================================================================================
  ///  \param[in]  span_u     knotspan index in u-direction
  ///  \param[in]  span_v     knotspan index in v-direction
  ///  \param[in]  _u         local parameter in u-direction
  ///  \param[in]  _v         local parameter in v-direction
  ///  \param[out] R         basis func
  ///
  /// ======================================================================================
  ///  \author     Daniel Baumgärtner (12/2016)
  //
  //########################################################################################
  void BrepFace::EvaluateNURBSFunctions(int span_u, int span_v, double _u, double _v, Matrix& R, ModelPart& model_part)
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
        R(b, c) = N(b)*M(c)*model_part.GetNode(m_control_points_ids[control_point_index]).GetValue(CONTROL_POINT_WEIGHT);
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


  // #######################################################################################
  ///
  ///  \details    returns the first and second derivative of NURBS basis function w.r.t. u,v
  ///              span_u,span_v are the knot span indices. if unknown, insert 0!
  ///
  /// ======================================================================================
  ///  \param[in]  span_u     knotspan index in u-direction
  ///  \param[in]  _u     	local parameter in u-direction
  ///  \param[in]  span_v 	knotspan index in v-direction
  ///  \param[in]  _v         local parameter in v-direction
  ///  \param[out] _dR        1st derivatives
  ///  \param[out] _ddR       2nd derivatives
  ///
  /// ======================================================================================
  ///  \author     from M.Breitenberger in Carat (12/2009)
  //
  //########################################################################################

  void BrepFace::EvaluateNURBSFunctionsDerivatives(int span_u, int span_v, double _u, double _v, Matrix& _dR, Matrix& _ddR, ModelPart& model_part)
  {
    if (span_u == -1) span_u = NurbsUtilities::find_knot_span(m_p, m_knot_vector_u, _u);
    if (span_v == -1) span_v = NurbsUtilities::find_knot_span(m_q, m_knot_vector_v, _v);

    int ne = (m_p + 1)*(m_q + 1); // Control Points per element
    Matrix N;              // Basisfunc at _u
    Matrix M;              // Basisfunc at _v
    NurbsUtilities::eval_nonzero_basis_function_with_derivatives(N, m_knot_vector_u, _u, span_u, m_p, 2);
    NurbsUtilities::eval_nonzero_basis_function_with_derivatives(M, m_knot_vector_v, _v, span_v, m_q, 2);

    vector<double> r(ne);
    r.clear();
    _dR.resize(ne, 2);
    _ddR.resize(ne, 3);
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
        weight = model_part.GetNode(m_control_points_ids[control_point_index]).GetValue(CONTROL_POINT_WEIGHT);

        r[k] = N(0, b)*M(0, c)*weight;
        sum += r[k];
        //First derivatives
        _dR(k, 0) = N(1, b)*M(0, c)*weight;
        dsum[0] += _dR(k, 0);
        _dR(k, 1) = N(0, b)*M(1, c)*weight;
        dsum(1) += _dR(k, 1);
        //Second derivatives  1-du^2, 2-dv^2, 3-dudv
        _ddR(k, 0) = N(2, b)*M(0, c)*weight;
        ddsum(0) = ddsum(0) + _ddR(k, 0);
        _ddR(k, 1) = N(0, b)*M(2, c)*weight;
        ddsum(1) = ddsum(1) + _ddR(k, 1);
        _ddR(k, 2) = N(1, b)*M(1, c)*weight;
        ddsum(2) = ddsum(2) + _ddR(k, 2);
        k++;
      }
    }
    double sum_2 = pow(sum, 2);
    double sum_3 = pow(sum, 3);
    // divide through by sum
    for (int k = 0; k<ne; k++)
    {

      _ddR(k, 0) = _ddR(k, 0) / sum - 2.0*_dR(k, 0)*dsum[0] / sum_2
        - r[k] * ddsum[0] / sum_2 + 2.0*r[k] * dsum[0] * dsum[0] / sum_3;
      _ddR(k, 1) = _ddR(k, 1) / sum - 2.0*_dR(k, 1)*dsum[1] / sum_2
        - r[k] * ddsum[1] / sum_2 + 2.0*r[k] * dsum[1] * dsum[1] / sum_3;
      _ddR(k, 2) = _ddR(k, 2) / sum - _dR(k, 0)*dsum[1] / sum_2 - _dR(k, 1)*dsum[0] / sum_2
        - r[k] * ddsum[2] / sum_2 + 2.0*r[k] * dsum[0] * dsum[1] / sum_3;
      _dR(k, 0) = _dR(k, 0) / sum - r[k] * dsum[0] / sum_2;
      _dR(k, 1) = _dR(k, 1) / sum - r[k] * dsum[1] / sum_2;
    }
  }

  //  #####################################################################################
  // #######################################################################################
  ///
  ///  \details    returns the basis fucntions and the first derivative of NURBS basis function w.r.t. u,v
  ///              _i,_j are the knot span indices. if unknown, insert 0!
  ///
  /// ======================================================================================
  ///  \param[in]  _i         knotspan index in u-direction
  ///  \param[in]  _u         local parameter in u-direction
  ///  \param[in]  _j         knotspan index in v-direction
  ///  \param[in]  _v         local parameter in v-direction
  ///  \param[out] _R         basis func
  ///  \param[out] _dR        1st derivatives
  ///
  /// ======================================================================================
  ///  \author     from M.Breitenberger in Carat (12/2009)
  //
  //########################################################################################
  void BrepFace::EvaluateNURBSFunctionsAndDerivative(int span_u, int span_v, double _u, double _v, Matrix& R, std::vector<Matrix>& dR, ModelPart& model_part)
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
        weight = model_part.GetNode(m_control_points_ids[control_point_index]).GetValue(CONTROL_POINT_WEIGHT);
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
    Vector& knot_vector_u, Vector& knot_vector_v,
    unsigned int& p, unsigned int& q, IntVector& control_point_ids)
    : m_trimming_loops(trimming_loops),
    m_knot_vector_u(knot_vector_u),
    m_knot_vector_v(knot_vector_v),
    m_p(p),
    m_q(q),
    m_control_points_ids(control_point_ids),
    IndexedObject(brep_id),
    Flags()
  {
  }
///Destructor
BrepFace::~BrepFace()
{}

}  // namespace Kratos.

