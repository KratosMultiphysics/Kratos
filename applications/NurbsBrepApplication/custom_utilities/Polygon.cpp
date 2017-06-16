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
#include "Polygon.h"


namespace Kratos
{
  //bool Polygon::GetOrientation() {
  //  long i1, i2;
  //  double area = 0;
  //  std::vector<PointXYType> const& points = m_polygon.outer();
  //  int number_of_points = boost::geometry::num_points(m_polygon);
  //  for (i1 = 0; i1<number_of_points; i1++) {
  //    i2 = i1 + 1;
  //    if (i2 == number_of_points) i2 = 0;
  //    area += points[i1].x() * points[i2].y() - points[i1].y() * points[i2].x();
  //  }
  //  if (area>0) return false;
  //  if (area<0) return true;
  //  return 0;
  //}

  //void Polygon::Invert() {
  //  int i;
  //  std::vector<PointXYType> const& points = m_polygon.outer();
  //  int number_of_points = boost::geometry::num_points(m_polygon);
  //  std::vector<PointXYType> inverse_points(number_of_points);
  //  //inverse_points.resize(m_polygon.size());
  //  for (i = 0; i<number_of_points; i++) {
  //    inverse_points[i] = points[number_of_points - i - 1];
  //  }
  //  m_polygon.outer = inverse_points;
  //}

  /**
  * @date   Mai, 2017
  * @brief   Basic function to triangulate polygon. Produces minimal number of triangles.
  *          Function from empire-multiphysics. http://empire-multiphysics.com/projects/empire
  * @param [in] polygon, polygon to triangulate space.
  * @param [out] triangles, list of the outcoming triangles.
  */
  bool Polygon::Triangulate_OPT(const PolygonType& polygon, std::vector<Matrix>& triangles)
  {
    PointXYType p1, p2, p3, p4;
    int bestvertex;
    double weight, minweight, d1, d2;
    Diagonal diagonal, newdiagonal;
    std::list<Diagonal> diagonals;
    bool ret = true;

    int n = boost::geometry::num_points(polygon)-1;
    //std::cout << "number of points: " << n << std::endl;
    std::vector<PointXYType> const& points = polygon.outer();
    matrix<DPState> dpstates(n, n);

    //std::cout << "check 2" << std::endl;
    //init states and visibility
    for (unsigned int i = 0; i<(n - 1); i++) {
      p1 = points[i];
      for (unsigned int j = i + 1; j<n; j++) {
        dpstates(j, i).visible = true;
        dpstates(j, i).weight = 0;
        dpstates(j, i).bestvertex = -1;
        if (j != (i + 1)) {
          p2 = points[j];

          //visibility check
          if (i == 0) p3 = points[n - 1];
          else p3 = points[i - 1];
          if (i == (n - 1)) p4 = points[0];
          else p4 = points[i + 1];
          if (!InCone(p3, p1, p4, p2)) {
            dpstates(j, i).visible = false;
            continue;
          }

          if (j == 0) p3 = points[n - 1];
          else p3 = points[j - 1];
          if (j == (n - 1)) p4 = points[0];
          else p4 = points[j + 1];
          if (!InCone(p3, p2, p4, p1)) {
            dpstates(j, i).visible = false;
            continue;
          }

          for (unsigned int k = 0; k<n; k++) {
            p3 = points[k];
            if (k == (n - 1)) p4 = points[0];
            else p4 = points[k + 1];
            if (Intersects(p1, p2, p3, p4)) {
              dpstates(j, i).visible = false;
              break;
            }
          }
        }
      }
    }
    //std::cout << "check 3" << std::endl;
    dpstates(n - 1, 0).visible = true;
    dpstates(n - 1, 0).weight = 0;
    dpstates(n - 1, 0).bestvertex = -1;

    for (unsigned int gap = 2; gap<n; gap++) {
      for (unsigned int i = 0; i<(n - gap); i++) {
        int j = i + gap;
        if (!dpstates(j, i).visible) continue;
        bestvertex = -1;
        for (unsigned int k = (i + 1); k<j; k++) {
          if (!dpstates(k, i).visible) continue;
          if (!dpstates(j, k).visible) continue;

          if (k <= (i + 1)) d1 = 0;
          else d1 = Distance(points[i], points[k]);
          if (j <= (k + 1)) d2 = 0;
          else d2 = Distance(points[k], points[j]);

          weight = dpstates(k, i).weight + dpstates(j, k).weight + d1 + d2;

          if ((bestvertex == -1) || (weight<minweight)) {
            bestvertex = k;
            minweight = weight;
          }
        }
        if (bestvertex == -1) {
          return false;
        }

        dpstates(j, i).bestvertex = bestvertex;
        dpstates(j, i).weight = minweight;
      }
    }

    newdiagonal.index1 = 0;
    newdiagonal.index2 = n - 1;
    diagonals.push_back(newdiagonal);
    while (!diagonals.empty()) {
      diagonal = *(diagonals.begin());
      diagonals.pop_front();
      bestvertex = dpstates(diagonal.index2, diagonal.index1).bestvertex;
      if (bestvertex == -1) {
        ret = false;
        break;
      }
      Matrix triangle(3, 2);
      triangle(0, 0) = points[diagonal.index1].x();
      triangle(0, 1) = points[diagonal.index1].y();
      triangle(1, 0) = points[bestvertex].x();
      triangle(1, 1) = points[bestvertex].y();
      triangle(2, 0) = points[diagonal.index2].x();
      triangle(2, 1) = points[diagonal.index2].y();
      if (abs(GetAreaOfTriangle(triangle))>1e-9)
        triangles.push_back(triangle);
      else
      {
        std::cout << "triangle with zero area" << GetAreaOfTriangle(triangle) << std::endl;
        //KRATOS_WATCH(triangle)
      }
      if (bestvertex > (diagonal.index1 + 1)) {
        newdiagonal.index1 = diagonal.index1;
        newdiagonal.index2 = bestvertex;
        diagonals.push_back(newdiagonal);
      }
      if (diagonal.index2 > (bestvertex + 1)) {
        newdiagonal.index1 = bestvertex;
        newdiagonal.index2 = diagonal.index2;
        diagonals.push_back(newdiagonal);
      }
    }
    //std::cout << "check 5" << std::endl;

    //for (unsigned int i = 1; i<n; i++) {
    //  delete[] dpstates[i];
    //}
    //delete[] dpstates;
    //triangles;// = triangles;
    return true;
  }

  //bool Polygon::Triangulate_OPT(std::vector<Matrix>& triangles)
  //{
  //  //std::vector<Matrix> triangles;
  //  std::cout << "check 1" << std::endl;
  //  //DPState **dpstates;
  //  array_1d<double, 2> p1, p2, p3, p4;
  //  long bestvertex;
  //  double weight, minweight, d1, d2;
  //  Diagonal diagonal, newdiagonal;
  //  std::list<Diagonal> diagonals;
  //  bool ret = true;
  //
  //  int n = m_polygon.size();
  //  matrix<DPState> dpstates(n, n);
  //  //dpstates = new DPState *[n];
  //  //for (unsigned int i = 1; i<n; i++) {
  //  //  dpstates[i] = new DPState[i];
  //  //}
  //
  //  std::cout << "check 2" << std::endl;
  //  //init states and visibility
  //  for (unsigned int i = 0; i<(n - 1); i++) {
  //    p1 = m_polygon[i];
  //    for (unsigned int j = i + 1; j<n; j++) {
  //      dpstates(j,i).visible = true;
  //      dpstates(j,i).weight = 0;
  //      dpstates(j,i).bestvertex = -1;
  //      if (j != (i + 1)) {
  //        p2 = m_polygon[j];
  //
  //        //visibility check
  //        if (i == 0) p3 = m_polygon[n - 1];
  //        else p3 = m_polygon[i - 1];
  //        if (i == (n - 1)) p4 = m_polygon[0];
  //        else p4 = m_polygon[i + 1];
  //        if (!InCone(p3, p1, p4, p2)) {
  //          dpstates(j,i).visible = false;
  //          continue;
  //        }
  //
  //        if (j == 0) p3 = m_polygon[n - 1];
  //        else p3 = m_polygon[j - 1];
  //        if (j == (n - 1)) p4 = m_polygon[0];
  //        else p4 = m_polygon[j + 1];
  //        if (!InCone(p3, p2, p4, p1)) {
  //          dpstates(j,i).visible = false;
  //          continue;
  //        }
  //
  //        for (unsigned int k = 0; k<n; k++) {
  //          p3 = m_polygon[k];
  //          if (k == (n - 1)) p4 = m_polygon[0];
  //          else p4 = m_polygon[k + 1];
  //          if (Intersects(p1, p2, p3, p4)) {
  //            dpstates(j,i).visible = false;
  //            break;
  //          }
  //        }
  //      }
  //    }
  //  }
  //  std::cout << "check 3" << std::endl;
  //  dpstates(n-1,0).visible = true;
  //  dpstates(n-1,0).weight = 0;
  //  dpstates(n-1,0).bestvertex = -1;
  //
  //  for (unsigned int gap = 2; gap<n; gap++) {
  //    for (unsigned int i = 0; i<(n - gap); i++) {
  //      int j = i + gap;
  //      if (!dpstates(j,i).visible) continue;
  //      bestvertex = -1;
  //      for (unsigned int k = (i + 1); k<j; k++) {
  //        if (!dpstates(k,i).visible) continue;
  //        if (!dpstates(j,k).visible) continue;
  //
  //        if (k <= (i + 1)) d1 = 0;
  //        else d1 = Distance(m_polygon[i], m_polygon[k]);
  //        if (j <= (k + 1)) d2 = 0;
  //        else d2 = Distance(m_polygon[k], m_polygon[j]);
  //
  //        weight = dpstates(k,i).weight + dpstates(j,k).weight + d1 + d2;
  //
  //        if ((bestvertex == -1) || (weight<minweight)) {
  //          bestvertex = k;
  //          minweight = weight;
  //        }
  //      }
  //      if (bestvertex == -1) {
  //        //for (unsigned int i = 1; i<n; i++) {
  //        //  delete[] dpstates[i];
  //        //}
  //        //delete[] dpstates;
  //
  //        return false;
  //      }
  //
  //      dpstates(j,i).bestvertex = bestvertex;
  //      dpstates(j,i).weight = minweight;
  //    }
  //  }
  //
  //  newdiagonal.index1 = 0;
  //  newdiagonal.index2 = n - 1;
  //  diagonals.push_back(newdiagonal);
  //  while (!diagonals.empty()) {
  //    diagonal = *(diagonals.begin());
  //    diagonals.pop_front();
  //    bestvertex = dpstates(diagonal.index2,diagonal.index1).bestvertex;
  //    if (bestvertex == -1) {
  //      ret = false;
  //      break;
  //    }
  //    Matrix triangle(3, 2);
  //    triangle(0, 0) = m_polygon[diagonal.index1][0];
  //    triangle(0, 1) = m_polygon[diagonal.index1][1];
  //    triangle(1, 0) = m_polygon[bestvertex][0];
  //    triangle(1, 1) = m_polygon[bestvertex][1];
  //    triangle(2, 0) = m_polygon[diagonal.index2][0];
  //    triangle(2, 1) = m_polygon[diagonal.index2][1];
  //    if (GetAreaOfTriangle(triangle)>1e-9)
  //      triangles.push_back(triangle);
  //    if (bestvertex > (diagonal.index1 + 1)) {
  //      newdiagonal.index1 = diagonal.index1;
  //      newdiagonal.index2 = bestvertex;
  //      diagonals.push_back(newdiagonal);
  //    }
  //    if (diagonal.index2 > (bestvertex + 1)) {
  //      newdiagonal.index1 = bestvertex;
  //      newdiagonal.index2 = diagonal.index2;
  //      diagonals.push_back(newdiagonal);
  //    }
  //  }
  //  std::cout << "check 5" << std::endl;
  //
  //  //for (unsigned int i = 1; i<n; i++) {
  //  //  delete[] dpstates[i];
  //  //}
  //  //delete[] dpstates;
  //  //triangles;// = triangles;
  //  return true;
  //}

  double Polygon::GetAreaOfTriangle(const Matrix& triangle)
  {
    return abs((triangle(0, 0)*(triangle(1, 1) - triangle(2, 1))
      + triangle(1, 0)*(triangle(2, 1) - triangle(0, 1))
      + triangle(2, 0)*(triangle(0, 1) - triangle(1, 1))) / 2);
  }

  bool Polygon::IsInside(const double& u, const double& v)
  {
    PointXYType point(u, v);
    bool is_inside = false;
    for (auto polygon = m_polygon_list.begin(); polygon != m_polygon_list.end(); ++polygon)
    {
      is_inside = boost::geometry::within(point, (*polygon));
      if (is_inside)
        return is_inside;
    }
    return is_inside;
  }

  void Polygon::Reverse(const int& index)
  {
    boost::geometry::reverse(m_polygon_list[index]);
  }

  std::vector<Matrix> Polygon::Triangulate()
  {
    std::vector<Matrix> triangles;
    //std::cout << "m_polygon: " << m_polygon_list.size() << std::endl;
    //for (auto polygon = m_polygon_list.begin(); polygon != m_polygon_list.end(); ++polygon)
    for (unsigned int i = 0; i<m_polygon_list.size(); i++)
    {
      std::cout << boost::geometry::wkt<PolygonType>(m_polygon_list[i]) << std::endl;
      if (boost::geometry::area(m_polygon_list[i]) > 0)
      {
        Reverse(i);
        //std::cout << "reverse polygon: "<< boost::geometry::wkt<PolygonType>(m_polygon_list[i]) << std::endl;
        //std::cout << "reverse polygon are: " << boost::geometry::area(m_polygon_list[i]) << std::endl;
      }

      if (boost::geometry::num_points(m_polygon_list[i]) > 0)
      {
        bool success = Triangulate_OPT((m_polygon_list[i]), triangles);
        if (!success)
          for (std::vector<PolygonType>::const_iterator it = m_polygon_list.begin(); it != m_polygon_list.end(); ++it)
          {
            std::cout << "Nicht Triangulierbares Polygon: " << boost::geometry::wkt<PolygonType>(*it) << std::endl;
          }
          //KRATOS_THROW_ERROR(std::runtime_error, "Polygon::Triangulate: Triangulation failed.", success);
      }
      else
        KRATOS_THROW_ERROR(std::runtime_error, "Polygon::Triangulate: No points in polygon.", std::endl);
    }
    //std::cout << triangles.size() << std::endl;
    //for (unsigned int t = 0; t < triangles.size(); t++)
    //{
    //  //KRATOS_WATCH(triangles[t])
    //}
    return triangles;
  }

  bool Polygon::IsFullKnotSpan()
  {
    return m_is_full_knot_span;
  }

  Polygon Polygon::clipByKnotSpan(const Vector& parameter_span_u, const Vector& parameter_span_v) {
    /// 1.find the knot span which the current element located in.
    //      from minSpanu to maxSpanu in U-direction, and from minSpanV to max SpanV in V-direction
    std::vector<PointXYType> points;
    points.push_back(PointXYType(parameter_span_u[0], parameter_span_v[0]));
    points.push_back(PointXYType(parameter_span_u[1], parameter_span_v[0]));
    points.push_back(PointXYType(parameter_span_u[1], parameter_span_v[1]));
    points.push_back(PointXYType(parameter_span_u[0], parameter_span_v[1]));

    PolygonType polygon;
    boost::geometry::assign_points(polygon, points);
    boost::geometry::correct(polygon);

    PolygonVectorType polygon_vector;
    for (PolygonVectorType::const_iterator it = m_polygon_list.begin(); it != m_polygon_list.end(); ++it)
    {
      clipPolygon((*it), polygon, polygon_vector);
    }
    Polygon new_polygon(polygon_vector);
    //std::cout << "area full knot span: " << abs(boost::geometry::area(polygon)) << ", area of clipped polygon: " << new_polygon.GetArea() << std::endl;
    if (abs(abs(boost::geometry::area(polygon)) - new_polygon.GetArea()) < 1e-7)
      new_polygon.m_is_full_knot_span = true;
    return new_polygon;
  }

  void Polygon::clipPolygon(const PolygonType& polygon_1, const PolygonType& polygon_2, PolygonVectorType& polygon_vector)
  {
    std::deque<PolygonType> output;
    boost::geometry::intersection(polygon_1, polygon_2, output);

    for (std::deque<PolygonType>::const_iterator it = output.begin(); it != output.end(); ++it)
    {
      //std::cout << "Geclipptes Polygon: " << boost::geometry::wkt<PolygonType>(*it) << std::endl;
      polygon_vector.push_back(*it);
    }
  }

  //std::vector<array_1d<double, 2>> Polygon::clipPolygon(std::vector<array_1d<double, 2>> polygon_1, std::vector<array_1d<double, 2>> polygon_2)
  //{
  //  // Type definitions to use boost functionalities
  //  //using namespace boost::geometry;
  //
  //  //typedef model::d2::point_xy<double> PointXYType;
  //  //typedef model::polygon<PointXYType> local_polygon_type;
  //
  //  PolygonType boost_polygon_1, boost_polygon_2;
  //
  //  for (unsigned int i = 0; i < polygon_1.size(); i++)
  //    boost::geometry::append(boost::geometry::exterior_ring(boost_polygon_1),
  //      boost::geometry::make<PointXYType>(polygon_1[i][0], polygon_1[i][1]));
  //  if (polygon_1.size() > 0)
  //    boost::geometry::append(boost::geometry::exterior_ring(boost_polygon_1),
  //      boost::geometry::make<PointXYType>(polygon_1[0][0], polygon_1[0][1]));
  //  boost::geometry::correct(boost_polygon_1);
  //
  //  for (unsigned int i = 0; i < polygon_2.size(); i++)
  //    boost::geometry::append(boost::geometry::exterior_ring(boost_polygon_2),
  //      boost::geometry::make<PointXYType>(polygon_2[i][0], polygon_2[i][1]));
  //  if (polygon_2.size() > 0)
  //    boost::geometry::append(boost::geometry::exterior_ring(boost_polygon_2),
  //      boost::geometry::make<PointXYType>(polygon_2[0][0], polygon_2[0][1]));
  //  boost::geometry::correct(boost_polygon_2);
  //
  //  std::deque<PolygonType> output;
  //  boost::geometry::intersection(boost_polygon_1, boost_polygon_2, output);
  //
  //  std::vector<array_1d<double, 2>> boundary_polygon;
  //  for (std::deque<PolygonType>::const_iterator it = output.begin(); it != output.end(); ++it)
  //  {
  //    std::cout << boost::geometry::wkt<PolygonType>(*it) << std::endl;
  //    std::vector<PointXYType> const& points = (*it).outer();
  //
  //    for (unsigned int i = 0; i < points.size(); i++)
  //    {
  //      array_1d<double, 2> point;
  //      point[0] = points[i].x();
  //      point[1] = points[i].y();
  //
  //      boundary_polygon.push_back(point);
  //    }
  //  }
  //  return boundary_polygon;
  //}

  bool Polygon::InCone(array_1d<double, 2> &p1, array_1d<double, 2> &p2,
    array_1d<double, 2> &p3, array_1d<double, 2> &p) 
  {
    if (IsConvex(p1, p2, p3)) {
      if (!IsConvex(p1, p2, p)) return false;
      if (!IsConvex(p2, p3, p)) return false;
      return true;
    }
    else {
      if (IsConvex(p1, p2, p)) return true;
      if (IsConvex(p2, p3, p)) return true;
      return false;
    }
  }

  bool Polygon::InCone(PointXYType &p1, PointXYType &p2,
    PointXYType &p3, PointXYType &p)
  {
    if (IsConvex(p1, p2, p3)) {
      if (!IsConvex(p1, p2, p)) return false;
      if (!IsConvex(p2, p3, p)) return false;
      return true;
    }
    else {
      if (IsConvex(p1, p2, p)) return true;
      if (IsConvex(p2, p3, p)) return true;
      return false;
    }
  }

  bool Polygon::IsConvex(
    const array_1d<double, 2>& p1, const array_1d<double, 2>& p2, 
    const array_1d<double, 2>& p3) 
  {
    double tmp;
    tmp = (p3[1] - p1[1])*(p2[0] - p1[0]) - (p3[0] - p1[0])*(p2[1] - p1[1]);
    if (tmp>0) return true;
    else return false;
  }

  bool Polygon::IsConvex(
    const PointXYType& p1, const PointXYType& p2,
    const PointXYType& p3)
  {
    double tmp;
    tmp = (p3.y() - p1.y())*(p2.x() - p1.x()) - (p3.x() - p1.x())*(p2.y() - p1.y());
    if (tmp>0) return true;
    else return false;
  }

  double Polygon::Distance(array_1d<double, 2> point_1, array_1d<double, 2> point_2)
  {
    return sqrt(point_1[0] * point_2[0] + point_1[1] * point_2[1]);
  }

  double Polygon::Distance(PointXYType point_1, PointXYType point_2)
  {
    return sqrt(point_1.x() * point_2.x() + point_1.y() * point_2.y());
  }

  //checks if two lines intersect
  bool Polygon::Intersects(array_1d<double, 2> &p11, array_1d<double, 2> &p12,
    array_1d<double, 2> &p21, array_1d<double, 2> &p22)
  {
    if ((p11[0] == p21[0]) && (p11[1] == p21[1])) return false;
    if ((p11[0] == p22[0]) && (p11[1] == p22[1])) return false;
    if ((p12[0] == p21[0]) && (p12[1] == p21[1])) return false;
    if ((p12[0] == p22[0]) && (p12[1] == p22[1])) return false;

    array_1d<double, 2> v1ort, v2ort, v;
    double dot11, dot12, dot21, dot22;

    v1ort[0] = p12[1] - p11[1];
    v1ort[1] = p11[0] - p12[0];

    v2ort[0] = p22[1] - p21[1];
    v2ort[1] = p21[0] - p22[0];

    v[0] = p21[0] - p11[0];
    v[1] = p21[1] - p11[1];
    dot21 = v[0] * v1ort[0] + v[1] * v1ort[1];
    v[0] = p22[0] - p11[0];
    v[1] = p22[1] - p11[1];
    dot22 = v[0] * v1ort[0] + v[1] * v1ort[1];

    v[0] = p11[0] - p21[0];
    v[1] = p11[1] - p21[1];
    dot11 = v[0] * v2ort[0] + v[1] * v2ort[1];
    v[0] = p12[0] - p21[0];
    v[1] = p12[1] - p21[1];
    dot12 = v[0] * v2ort[0] + v[1] * v2ort[1];

    if (dot11*dot12>0) return false;
    if (dot21*dot22>0) return false;

    return true;
  }

  bool Polygon::Intersects(PointXYType &p11, PointXYType &p12,
    PointXYType &p21, PointXYType &p22)
  {
    if ((p11.x() == p21.x()) && (p11.y() == p21.y())) return false;
    if ((p11.x() == p22.x()) && (p11.y() == p22.y())) return false;
    if ((p12.x() == p21.x()) && (p12.y() == p21.y())) return false;
    if ((p12.x() == p22.x()) && (p12.y() == p22.y())) return false;

    PointXYType v1ort, v2ort, v;
    double dot11, dot12, dot21, dot22;

    v1ort.x(p12.y() - p11.y());
    v1ort.y(p11.x() - p12.x());

    v2ort.x(p22.y() - p21.y());
    v2ort.y(p21.x() - p22.x());

    v.x(p21.x() - p11.x());
    v.y(p21.y() - p11.y());
    dot21 = v.x() * v1ort.x() + v.y() * v1ort.y();
    v.x(p22.x() - p11.x());
    v.y(p22.y() - p11.y());
    dot22 = v.x() * v1ort.x() + v.y() * v1ort.y();

    v.x(p11.x() - p21.x());
    v.y(p11.y() - p21.y());
    dot11 = v.x() * v2ort.x() + v.y() * v2ort.y();
    v.x(p12.x() - p21.x());
    v.y(p12.y() - p21.y());
    dot12 = v.x() * v2ort.x() + v.y() * v2ort.y();

    if (dot11*dot12>0) return false;
    if (dot21*dot22>0) return false;

    return true;
  }


  double Polygon::GetArea()
  {
    double area = 0.0;

    for (PolygonVectorType::const_iterator it = m_polygon_list.begin(); it != m_polygon_list.end(); ++it)
    {
      area += abs(boost::geometry::area(*it));
    }
    return area;
  }

  ///Constructor
  Polygon::Polygon(PolygonVectorType polygon_list)
    : m_polygon_list(polygon_list)
  {
  }

  Polygon::Polygon(PolygonType polygon)
  {
    PolygonVectorType polygon_list;
    std::cout << boost::geometry::wkt<PolygonType>(polygon) << std::endl;
    polygon_list.push_back(polygon);
      m_polygon_list = polygon_list;
  }

  Polygon::Polygon(std::vector<array_1d<double, 2>> polygon)
  {
    PolygonVectorType polygon_list;
    PolygonType boost_polygon;

    for (unsigned int i = 0; i < polygon.size(); i++)
      boost::geometry::append(boost::geometry::exterior_ring(boost_polygon),
        boost::geometry::make<PointXYType>(polygon[i][0], polygon[i][1]));
    //if (polygon.size() > 0)
    //  boost::geometry::append(boost::geometry::exterior_ring(boost_polygon_1),
    //    boost::geometry::make<PointXYType>(polygon[0][0], polygon[0][1]));
    boost::geometry::correct(boost_polygon);

    std::cout << boost::geometry::wkt<PolygonType>(boost_polygon) << std::endl;

    polygon_list.push_back(boost_polygon);
    m_polygon_list = polygon_list;
  }

  Polygon::Polygon(std::vector<BrepBoundaryLoop>& boundary_loops)
  {
    PolygonVectorType polygon_outer_list;
    PolygonVectorType polygon_inner_list;
    PolygonVectorType polygon_list;
    std::vector<array_1d<double, 2>> boundary_polygon;
    std::vector<PointXYType> points;
    for (unsigned int loop_i = 0; loop_i < boundary_loops.size(); loop_i++)
    {
      boundary_polygon = boundary_loops[loop_i].GetBoundaryPolygon(5);

      for (unsigned int i = 0; i < boundary_polygon.size(); i++)
        points.push_back(PointXYType(boundary_polygon[i][0], boundary_polygon[i][1]));

      PolygonType polygon;
      boost::geometry::assign_points(polygon, points);
      boost::geometry::correct(polygon);

      if (boundary_loops[loop_i].IsOuterLoop())
        polygon_outer_list.push_back(polygon);
      else
        polygon_inner_list.push_back(polygon);
    }
    for (unsigned int i = 0; i < polygon_inner_list.size(); i++)
    {
      for (unsigned int j = 0; j < polygon_outer_list.size(); j++)
      {
        std::deque<PolygonType> output;
        boost::geometry::difference(polygon_outer_list[j], polygon_inner_list[i], output);
        if (output.size() == 1)
          polygon_outer_list[j] = output[0];
        else
          KRATOS_THROW_ERROR(std::runtime_error, "Error in boundary loop definition.", std::endl);
      }
    }
    m_polygon_list = polygon_outer_list;
  }



  ///Destructor
  Polygon::~Polygon()
  {
  }

}  // namespace Kratos.
