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
  std::vector<array_1d<double, 2>> KnotSpan1d::getIntegrationPointsInFullGaussianDomain()
  {
    std::vector<array_1d<double, 2>> IntegrationPoints;

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
    return IntegrationPoints;
  }

  std::vector<array_1d<double, 2>> KnotSpan1d::getIntegrationPointsInParameterDomain()
  {
    std::vector<array_1d<double, 2>> IntegrationPoints, IntegrationPointsInGaussianDomain;

    IntegrationPointsInGaussianDomain = this->getIntegrationPointsInFullGaussianDomain();

    double u1 = m_parameter_u[0];
    double u2 = m_parameter_u[1];
    double du = u2 - u1;

    //KRATOS_WATCH(u1)
    //  KRATOS_WATCH(u2)
    //  KRATOS_WATCH(du)

    double mapping = abs(u2 - u1)*0.5;

    for (unsigned int i = 0; i < IntegrationPointsInGaussianDomain.size(); i++)
    {
      array_1d<double, 3> point;

      point[0] = u1 + du*(IntegrationPointsInGaussianDomain[i][0] + 1)*0.5; //Parameter
      point[1] = mapping*IntegrationPointsInGaussianDomain[i][1];

      IntegrationPoints.push_back(point);
    }
    return IntegrationPoints;
  }


  // --------------------------------------------------------------------------
  ///Constructor
  KnotSpan1d::KnotSpan1d(unsigned int knot_span_1d_id,
    int p,
    Vector parameter_u)
    : m_p(p), m_parameter_u(parameter_u),
      IndexedObject(knot_span_1d_id),
      Flags()
  {
  }

  ///Destructor
  KnotSpan1d::~KnotSpan1d()
  {
  }

}  // namespace Kratos.