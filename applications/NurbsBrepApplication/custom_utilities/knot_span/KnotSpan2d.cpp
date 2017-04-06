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
#include "KnotSpan2d.h"


namespace Kratos
{
  // --------------------------------------------------------------------------
  std::vector<array_1d<double, 3>> KnotSpan2d::getIntegrationPointsInFullGaussianDomain()
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

  std::vector<array_1d<double, 3>> KnotSpan2d::getIntegrationPointsInParameterDomain()
  {
    std::vector<array_1d<double, 3>> IntegrationPoints, IntegrationPointsInGaussianDomain;

    IntegrationPointsInGaussianDomain = this->getIntegrationPointsInFullGaussianDomain();

    double u1 = m_parameter_span_u[0];
    double u2 = m_parameter_span_u[1];
    double v1 = m_parameter_span_v[0];
    double v2 = m_parameter_span_v[1];
    double du = u2 - u1;
    double dv = v2 - v1;

    //KRATOS_WATCH(u1)
    //KRATOS_WATCH(u2)
    //KRATOS_WATCH(v1)
    //KRATOS_WATCH(v2)
    //KRATOS_WATCH(du)
    //KRATOS_WATCH(dv)

    double mapping = (u1 - u2)*(v1 - v2)*0.25;

    if (m_is_untrimmed)
    {
      for (unsigned int i = 0; i < IntegrationPointsInGaussianDomain.size(); i++)
      {
        array_1d<double, 3> point;

        point[0] = (u2 + u1 + IntegrationPointsInGaussianDomain[i][0] * du)*0.5;
        point[1] = (v2 + v1 + IntegrationPointsInGaussianDomain[i][1] * dv)*0.5;
        point[2] = mapping*IntegrationPointsInGaussianDomain[i][2];

        //KRATOS_WATCH(point)

        IntegrationPoints.push_back(point);// [i*number_of_nodes_v + j] = point;
      }
    }
    return IntegrationPoints;
  }


  // --------------------------------------------------------------------------
  ///Constructor
  KnotSpan2d::KnotSpan2d(unsigned int knot_span_2d_id,
    bool is_untrimmed,
    int p, int q,
    Vector parameter_span_u,
    Vector parameter_span_v)
    : m_p(p),
      m_q(q),
      m_is_untrimmed(is_untrimmed),
      m_parameter_span_u(parameter_span_u),
      m_parameter_span_v(parameter_span_v),
      IndexedObject(knot_span_1d_id),
      Flags()
  {
  }

  ///Destructor
  KnotSpan2d::~KnotSpan2d()
  {
  }

}  // namespace Kratos.