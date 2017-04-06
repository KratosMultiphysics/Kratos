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
#include "KnotSpan1d.h"


namespace Kratos
{
  // --------------------------------------------------------------------------
  std::vector<array_1d<double, 3>> KnotSpan1d::getIntegrationPointsInFullGaussianDomain()
  {
    std::vector<array_1d<double, 3>> IntegrationPoints;

    int number_of_nodes_u = (m_p + 1);

    Vector coords_dir_gu = ZeroVector(m_p + 1);
    Vector weights_gu = ZeroVector(m_p + 1);

    IntegrationUtilities::getGaussPointLocations(number_of_nodes_u, coords_dir_gu, weights_gu);

    for (unsigned int i = 0; i<number_of_nodes_u; i++)
    {
        array_1d<double, 2> point;

        point[0] = coords_dir_gu(i);
        point[1] = weights_gu(i);

        //KRATOS_WATCH(point)

        IntegrationPoints.push_back(point);// [i*number_of_nodes_v + j] = point;
      }
    }
    return IntegrationPoints;
  }

  std::vector<array_1d<double, 3>> KnotSpan1d::getIntegrationPointsInParameterDomain()
  {
    std::vector<array_1d<double, 3>> IntegrationPoints, IntegrationPointsInGaussianDomain;

    IntegrationPointsInGaussianDomain = this->getIntegrationPointsInFullGaussianDomain();

    double u1 = m_parameter_u[0];
    double u2 = m_parameter_u[1];
    double du = u2 - u1;

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
  KnotSpan1d::KnotSpan1d(unsigned int knot_span_1d_id,
    int p,
    Vector parameter_span_u)
    : m_p(p),
      m_parameter_span_u(parameter_span_u),
      IndexedObject(knot_span_1d_id),
      Flags()
  {
  }

  ///Destructor
  KnotSpan1d::~KnotSpan1d()
  {
  }

}  // namespace Kratos.