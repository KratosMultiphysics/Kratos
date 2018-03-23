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
#include "KnotSpan2dNIntegrate.h"
//#include "KnotSpan2d.h"


namespace Kratos
{
  // --------------------------------------------------------------------------
  std::vector<array_1d<double, 3>> KnotSpan2dNIntegrate::getIntegrationPointsInFullGaussianDomain()
  {
    std::vector<array_1d<double, 3>> IntegrationPoints;
    //IntegrationPoints.resize((m_p + 1)*(m_q + 1));

    int number_of_nodes_u = (m_p + 1);
    int number_of_nodes_v = (m_q + 1);

    Vector coords_dir_gu = ZeroVector(m_p + 1);
    Vector weights_gu = ZeroVector(m_p + 1);
    Vector coords_dir_gv = ZeroVector(m_q + 1);
    Vector weights_gv = ZeroVector(m_q + 1);

    IntegrationUtilities::getGaussPointLocations(number_of_nodes_u, coords_dir_gu, weights_gu);
    IntegrationUtilities::getGaussPointLocations(number_of_nodes_v, coords_dir_gv, weights_gv);

    for (unsigned int i = 0; i<number_of_nodes_u; i++)
    {
      for (unsigned int j = 0; j<number_of_nodes_v; j++)
      {
        array_1d<double, 3> point;

        point[0] = coords_dir_gu(i);
        point[1] = coords_dir_gv(j);
        point[2] = weights_gu(i)*weights_gv(j);

        //KRATOS_WATCH(point)

        IntegrationPoints.push_back(point);// [i*number_of_nodes_v + j] = point;
      }
    }
    return IntegrationPoints;
  }

  /**
  * @Author T.Oberbichler
  * @date   March, 2017
  * @brief   returns location of gauss points for triangles in gaussian space from -1 to 1.
  *
  * @param [in] triangle in parameter space. Matrix with size 3x2.
  *
  * @return std::vector<array_1d<double, 3>> IntegrationPoints
  */
  std::vector<array_1d<double, 3>> KnotSpan2dNIntegrate::pointsByTriangle(const Matrix &triangle)
  {
    std::vector<array_1d<double, 3>> IntegrationPoints;

    int degree = m_p;
    if (m_q > m_p)
      degree = m_q;

    Matrix coords_dir_gu;// = ZeroVector(number_of_nodes_u);
    Vector weights_gu;// = ZeroVector(number_of_nodes_u);
                      //Vector coords_dir_gv = ZeroVector(number_of_nodes_v);
                      //Vector weights_gv = ZeroVector(number_of_nodes_v);

    IntegrationUtilities::getGaussPointLocationsTriangle(degree, coords_dir_gu, weights_gu);

    Matrix dn = ZeroMatrix(2, 3);
    dn(0, 0) = -1;
    dn(0, 1) = 1;
    dn(0, 2) = 0;
    dn(1, 0) = -1;
    dn(1, 1) = 0;
    dn(1, 2) = 1;
    //KRATOS_WATCH(dn)
    Matrix Jacobian = prod(dn, triangle);

    double detJacobian = Jacobian(0, 0)*Jacobian(1, 1) - Jacobian(1, 0)*Jacobian(0, 1);
    //std::cout << "detJacobian: " << detJacobian << std::endl;

    for (unsigned int i = 0; i<coords_dir_gu.size1(); i++)
    {
      Vector n = ZeroVector(3);
      n[0] = 1 - coords_dir_gu(i,0) - coords_dir_gu(i, 1);
      n[1] = coords_dir_gu(i, 0);
      n[2] = coords_dir_gu(i, 1);

      //KRATOS_WATCH(coords_dir_gu(i, 0))
      //KRATOS_WATCH(coords_dir_gu(i, 1))

      Vector position = prod(n, triangle);
      //KRATOS_WATCH(position)

      array_1d<double, 3> point;

      point[0] = position(0);
      point[1] = position(1);
      point[2] = 0.5 * detJacobian * weights_gu(i);

      //KRATOS_WATCH(point)

      IntegrationPoints.push_back(point);// [i*number_of_nodes_v + j] = point;
    }
    return IntegrationPoints;
  }

  /**
  * @Author T.Oberbichler
  * @date   March, 2017
  * @brief   returns location of gauss points for quads in gaussian space from -1 to 1.
  *
  *
  * @param [in] quad in parameter space. Matrix with size 4x2.
  *
  * @return std::vector<array_1d<double, 3>> IntegrationPoints
  */
  std::vector<array_1d<double, 3>> KnotSpan2dNIntegrate::pointsByQuad(const Matrix &quad)
  {
    std::vector<array_1d<double, 3>> IntegrationPoints;
    //IntegrationPoints.resize((m_p + 1)*(m_q + 1));

    int number_of_nodes_u = (m_p + 1);
    int number_of_nodes_v = (m_q + 1);

    Vector coords_dir_gu, coords_dir_gv;// = ZeroVector(number_of_nodes_u);
    Vector weights_gu, weights_gv;// = ZeroVector(number_of_nodes_u);
    //Vector coords_dir_gv = ZeroVector(number_of_nodes_v);
    //Vector weights_gv = ZeroVector(number_of_nodes_v);

    IntegrationUtilities::getGaussPointLocations(number_of_nodes_u, coords_dir_gu, weights_gu);
    IntegrationUtilities::getGaussPointLocations(number_of_nodes_v, coords_dir_gv, weights_gv);

    for (unsigned int i = 0; i<number_of_nodes_u; i++)
    {
      for (unsigned int j = 0; j<number_of_nodes_v; j++)
      {
        Vector n = ZeroVector(4);
        n[0] = 0.25*(1 - coords_dir_gu(i))*(1 - coords_dir_gv(j));
        n[1] = 0.25*(1 + coords_dir_gu(i))*(1 - coords_dir_gv(j));
        n[2] = 0.25*(1 + coords_dir_gu(i))*(1 + coords_dir_gv(j));
        n[3] = 0.25*(1 - coords_dir_gu(i))*(1 + coords_dir_gv(j));

        Vector position = prod(n, quad);

        Matrix dn = ZeroMatrix(2, 4);
        dn(0,0) =  0.25 * (coords_dir_gv(j) - 1), 0.25 * ( 1 - coords_dir_gv(j)), 0.25 * (1 + coords_dir_gv(j)), 0.25 * (-1 - coords_dir_gv(j)),
              0.25 * (coords_dir_gu(i) - 1), 0.25 * (-1 - coords_dir_gu(i)), 0.25 * (1 + coords_dir_gu(i)), 0.25 * ( 1 - coords_dir_gu(i));

        Matrix Jacobian = prod(dn, quad);

        double detJacobian = Jacobian(0,0)*Jacobian(1, 1) - Jacobian(1, 0)*Jacobian(0, 1);

        array_1d<double, 3> point;

        point[0] = position(0);
        point[1] = position(1);
        point[2] = detJacobian*weights_gu(i)*weights_gv(j);

        //KRATOS_WATCH(point)

        IntegrationPoints.push_back(point);// [i*number_of_nodes_v + j] = point;
      }
    }
    return IntegrationPoints;
  }


  std::vector<array_1d<double, 3>> KnotSpan2dNIntegrate::getIntegrationPointsInParameterDomain()
  {
    std::vector<Matrix> faces = m_polygon.Triangulate();
    //std::cout << "faces: " << faces.size() << std::endl;
    std::vector<array_1d<double, 3>> IntegrationPoints, IntegrationPointsTriangle;

    int degree = m_p;
    if (m_q > m_p)
      degree = m_q;

    for (const auto face : faces) {
      switch (face.size1()) {
      case 3:
        IntegrationPointsTriangle = pointsByTriangle(face);
        for (unsigned int i = 0; i < IntegrationPointsTriangle.size(); i++)
          IntegrationPoints.push_back(IntegrationPointsTriangle[i]);
        break;
      case 4:
        IntegrationPointsTriangle = pointsByQuad(face);
        for (unsigned int i = 0; i < IntegrationPointsTriangle.size(); i++)
          IntegrationPoints.push_back(IntegrationPointsTriangle[i]);
        break;
      default:
        KRATOS_THROW_ERROR(std::runtime_error, "KnotSpan2dNIntegrate::getIntegrationPointsInParameterDomain(): Invalid face type", std::endl);
      }
    }
    std::ofstream file;
    file.open("triangles.txt", std::ios_base::app | std::ios_base::out);
    for (unsigned int i = 0; i < faces.size(); i++)
    {
      file << faces[i](0, 0) << " " << faces[i](1, 0) << " "  << faces[i](2, 0) << " " << "\n";
      file << faces[i](0, 1) << " " << faces[i](1, 1) << " "  << faces[i](2, 1) << " " << "\n";
    }
    file.close();

    return IntegrationPoints;
  }


  // --------------------------------------------------------------------------
  ///Constructor
KnotSpan2dNIntegrate::KnotSpan2dNIntegrate(
  unsigned int knot_span_2d_id,
  bool is_untrimmed,
  int p, int q,
  Vector parameter_span_u,
  Vector parameter_span_v,
  Polygon polygon)
    : KnotSpan2d(knot_span_2d_id, is_untrimmed, p, q, parameter_span_u, parameter_span_v),
      m_polygon(polygon)
  {
  }

  ///Destructor
  KnotSpan2dNIntegrate::~KnotSpan2dNIntegrate()
  {
  }

}  // namespace Kratos.