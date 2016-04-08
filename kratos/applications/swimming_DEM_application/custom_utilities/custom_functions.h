//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2011-6-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_CUSTOM_FUNCTIONS)
#define KRATOS_CUSTOM_FUNCTIONS

// /* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

// System includes
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"
#include "processes/find_elements_neighbours_process.h"
#include "processes/find_nodal_neighbours_process.h"

//Database includes
#include "custom_utilities/discrete_particle_configure.h"
#include "discrete_particle_configure.h"
#include "includes/define.h"
#include "../../DEM_application/custom_elements/discrete_element.h"
#include "custom_elements/spheric_swimming_particle.h"
#include "includes/define.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "../../DEM_application/custom_elements/spheric_particle.h"
#include "../swimming_DEM_application.h"
#include "../../../kratos/utilities/geometry_utilities.h"

namespace Kratos
{
template <std::size_t TDim>
class CustomFunctionsCalculator
{
public:

typedef ModelPart::ElementsContainerType::iterator  ElementIterator;
typedef ModelPart::NodesContainerType::iterator     NodeIterator;
typedef ModelPart::NodesContainerType               NodesArrayType;

KRATOS_CLASS_POINTER_DEFINITION(CustomFunctionsCalculator);

CustomFunctionsCalculator(): mPressuresFilled(false), mFirstGradientRecovery(true){}
/// Calculator

virtual ~CustomFunctionsCalculator(){}

/// Default calculator

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void CalculatePressureGradient(ModelPart& r_model_part)
{
    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        noalias(inode->FastGetSolutionStepValue(PRESSURE_GRADIENT)) = ZeroVector(3);
    }

    array_1d <double, 3> grad = ZeroVector(3); // its dimension is always 3
    array_1d <double, TDim + 1 > elemental_pressures;
    array_1d <double, TDim + 1 > N; // shape functions vector
    boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;

    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++){
        // computing the shape function derivatives
        Geometry<Node<3> >& geom = ielem->GetGeometry();
        double Volume;

        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

        // getting the pressure gradients;

        for (unsigned int i = 0; i < TDim + 1; ++i){
            elemental_pressures[i] = geom[i].FastGetSolutionStepValue(PRESSURE);
        }


        array_1d <double, TDim> grad_aux = prod(trans(DN_DX), elemental_pressures); // its dimension may be 2

        for (unsigned int i = 0; i < TDim; ++i){
            grad[i] = grad_aux[i];
        }

        double nodal_area = Volume / static_cast<double>(TDim + 1);
        grad *= nodal_area;

        for (unsigned int i = 0; i < TDim + 1; ++i){
            geom[i].FastGetSolutionStepValue(PRESSURE_GRADIENT) += grad;
        }
    }

    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        inode->FastGetSolutionStepValue(PRESSURE_GRADIENT) /= inode->FastGetSolutionStepValue(NODAL_AREA);
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

// This function is inspired in 2005 Zhang et al. "A new finite element gradient recovery method: superconvergence property"
//
// At a given node i, it is basically performing a least squares fit of a second degree polynomial, p(x), where x are the
// cartesian coordinates with i as their origin, of the nodal values of the pressure, P, for all nodes belonging to the
// elements concurrent to i (but not i). The nodal value of grad(P) at i then set to grad(p(0)).
//
// The algorithm has been simplified, so that if a node does not form an invertible matrix with its concurring elements'
// nodes, the current pressure_gradient is preserved (That's why we are still calling 'CalculatePressureGradient' at the
// beginning). This typically only happens near edges and corners in good-quality meshes.

void RecoverSuperconvergentPressureGradient(ModelPart& r_model_part)
{
    if (mFirstGradientRecovery){
        CalculateVectorsForGradientRecovery(r_model_part);
        mFirstGradientRecovery = false;
    }
    const double h_inv = 1.0 / CalculateTheMaximumEdgeLength(r_model_part); // we use it as a scaling parameter to improve stability

    // Getting nodal estimation from shape functions
    CalculatePressureGradient(r_model_part);

    // Solving least squares problem (Chen, 2016)
    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        array_1d <double, 3>& pressure_grad = inode->FastGetSolutionStepValue(PRESSURE_GRADIENT);
        const Vector& nodal_weights = inode->FastGetSolutionStepValue(NODAL_WEIGHTS);
        double sum = 0.0;
        for (unsigned int i = 0; i < nodal_weights.size(); ++i){
            sum += fabs(nodal_weights[i]);
        }
        if (sum > 0.0){
            pressure_grad = ZeroVector(3);
            WeakPointerVector<Element>& neigh_elems = inode->GetValue(NEIGHBOUR_ELEMENTS);
            unsigned int i = 0;

            for (unsigned int i_el = 0; i_el < neigh_elems.size(); ++i_el){
                Geometry<Node<3> >& geom = neigh_elems[i_el].GetGeometry();

                unsigned int jj = 0; // index of the node in geom corresponding to neighbour neigh_elems[i_el]
                if (geom[jj].Id() == inode->Id()){ // avoiding node inode
                    jj++;
                }

                for (unsigned int j = 0; j < TDim; ++j){
                    for (unsigned int d = 0; d < TDim; ++d){
                        pressure_grad[d] += nodal_weights[3 * (i + j) + d] * geom[jj].FastGetSolutionStepValue(PRESSURE);
                    }
                    jj++;
                }

                i += TDim;
            }

            pressure_grad *= h_inv; // Scaling back to physical units
        }
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void CalculateVectorsForGradientRecovery(ModelPart& r_model_part)
{
    using namespace boost::numeric::ublas;

    FindNodalNeighboursProcess neighbour_finder = FindNodalNeighboursProcess(r_model_part);
    neighbour_finder.Execute();

    vector<unsigned int> nodes_partition;
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_model_part.Nodes().size(), nodes_partition);

    // obtaining largest h to normalize the relative coordinates coordinates
    const double h = CalculateTheMaximumEdgeLength(r_model_part);
    const double h_inv = 1.0 / h; // we use it as a scaling parameter to improve stability
    //#pragma omp parallel for
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
        NodesArrayType& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();
        NodeIterator node_begin = pNodes.ptr_begin() + nodes_partition[k];
        NodeIterator node_end   = pNodes.ptr_begin() + nodes_partition[k + 1];
        int count = 0;
        for (ModelPart::NodesContainerType::iterator inode = node_begin; inode != node_end; ++inode){
            WeakPointerVector<Element>& neigh_elems = inode->GetValue(NEIGHBOUR_ELEMENTS);

            unsigned int n_nodal_neighs = TDim * neigh_elems.size(); // All nodes other than i_node. We obviously repeat nodes here (room for optimizatin!)
            unsigned int n_poly_terms;

            if (TDim == 3) {
                n_poly_terms = 10;
            }
            else { // TDim == 2
                n_poly_terms = 6;
                KRATOS_THROW_ERROR(std::runtime_error,"Gradient recovery not implemented yet in 2D!)","");
            }

            boost::numeric::ublas::matrix<double>A(n_nodal_neighs, n_poly_terms);
            unsigned int i = 0;
            for (unsigned int i_el = 0; i_el < neigh_elems.size(); ++i_el){
                Geometry<Node<3> >& geom = neigh_elems[i_el].GetGeometry();

                unsigned int jj = 0; // index for the position of the node j in geom
                for (unsigned int j = 0; j < TDim; ++j){
                     A(i + j, 0) = 1.0;

                     if (geom[jj].Id() == inode->Id()){ // avoiding node inode
                         jj++;
                     }

                     if (TDim == 3){
                        for (unsigned int d = 1; d < 10; ++d){
                            if (d < 4){
                                A(i + j, d) = (geom[jj][d - 1] - inode->Coordinates()[d - 1]) * h_inv;
                            }
                            else if (d == 4){
                                A(i + j, d) = (geom[jj][0] - inode->Coordinates()[0]) * (geom[jj][1] - inode->Coordinates()[1]) * h_inv * h_inv;
                            }
                            else if (d == 5){
                                A(i + j, d) = (geom[jj][0] - inode->Coordinates()[0]) * (geom[jj][2] - inode->Coordinates()[2]) * h_inv * h_inv;
                            }
                            else if (d == 6){
                                A(i + j, d) = (geom[jj][1] - inode->Coordinates()[1]) * (geom[jj][2] - inode->Coordinates()[2]) * h_inv * h_inv;
                            }
                            else {
                                double difference = (geom[jj][d - 7] - inode->Coordinates()[d - 7]) * h_inv;
                                A(i + j, d) = difference * difference;
                            }
                        }
                     }
                     else {
                         KRATOS_THROW_ERROR(std::runtime_error,"Gradient recovery not implemented yet in 2D!)","");
                     }
                     jj++;
                 }
                 i += TDim;
            }

            boost::numeric::ublas::matrix<double>AtransA(n_poly_terms, n_poly_terms);
            AtransA = prod(trans(A), A);
            Vector& nodal_weights = inode->FastGetSolutionStepValue(NODAL_WEIGHTS);
            nodal_weights.resize(3 * n_nodal_neighs);

            if (fabs(determinant< boost::numeric::ublas::matrix<double> >(AtransA)) < 0.0000001){ // if non-invertible we keep the regular approximation
                for (int i = 0; i < (int)nodal_weights.size(); ++i){
                    nodal_weights[i] = 0.0;
                }
                count++;
            }
            else {
                boost::numeric::ublas::matrix<double>AtransAinvAtrans(n_poly_terms, n_nodal_neighs);
                AtransAinvAtrans = prod(Inverse(AtransA), trans(A));

                // storing rows 1, 2 and 3 in a nodal vector variable
                for (unsigned int i = 0; i < n_nodal_neighs; ++i){
                    for (unsigned int d = 0; d < TDim; ++d){
                        nodal_weights(3 * i + d) = AtransAinvAtrans(d + 1, i);
                    }
                }
            }
        }
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double CalculateTheMaximumEdgeLength(ModelPart& r_model_part)
{
    double max_distance_yet = 0.0;

    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++){
        Geometry<Node<3> >& geom = ielem->GetGeometry();
        unsigned int n_nodes = static_cast<unsigned int>(TDim + 1);

        for (unsigned int k = 1; k < n_nodes - 1; ++k){
            for (unsigned int i = k; i < n_nodes; ++i){
                double distance_2 = 0.0;
                for (unsigned int d = 1; d < 3; ++d){
                    double delta_i = geom[k - 1][d] - geom[i][d];
                    distance_2 += delta_i * delta_i;
                }

                max_distance_yet = max_distance_yet > distance_2 ? max_distance_yet : distance_2;
            }
        }
    }

    return(std::sqrt(max_distance_yet));
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

template<class matrix_T>
double determinant(boost::numeric::ublas::matrix_expression<matrix_T> const& mat_r)
{
  double det = 1.0;

  matrix_T mLu(mat_r() );
  boost::numeric::ublas::permutation_matrix<std::size_t> pivots(mat_r().size1() );

  int is_singular = lu_factorize(mLu, pivots);

  if (!is_singular)
  {
    for (std::size_t i=0; i < pivots.size(); ++i)
    {
      if (pivots(i) != i)
        det *= -1.0;

      det *= mLu(i,i);
    }
  }
  else
    det = 0.0;

  return det;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

bool InvertMatrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse)
{
    double m[16];

    for (unsigned int i = 0; i < 4; ++i){
        for (unsigned int j = 0; j < 4; ++j){
            m[4 * i + j] = input(i, j);
        }
    }

    double inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] -
             m[5]  * m[11] * m[14] -
             m[9]  * m[6]  * m[15] +
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] -
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] +
              m[4]  * m[11] * m[14] +
              m[8]  * m[6]  * m[15] -
              m[8]  * m[7]  * m[14] -
              m[12] * m[6]  * m[11] +
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] -
             m[4]  * m[11] * m[13] -
             m[8]  * m[5] * m[15] +
             m[8]  * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] +
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] -
               m[8]  * m[6] * m[13] -
               m[12] * m[5] * m[10] +
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] +
              m[1]  * m[11] * m[14] +
              m[9]  * m[2] * m[15] -
              m[9]  * m[3] * m[14] -
              m[13] * m[2] * m[11] +
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] -
             m[0]  * m[11] * m[14] -
             m[8]  * m[2] * m[15] +
             m[8]  * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] +
              m[0]  * m[11] * m[13] +
              m[8]  * m[1] * m[15] -
              m[8]  * m[3] * m[13] -
              m[12] * m[1] * m[11] +
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] -
              m[0]  * m[10] * m[13] -
              m[8]  * m[1] * m[14] +
              m[8]  * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] -
             m[1]  * m[7] * m[14] -
             m[5]  * m[2] * m[15] +
             m[5]  * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] +
              m[0]  * m[7] * m[14] +
              m[4]  * m[2] * m[15] -
              m[4]  * m[3] * m[14] -
              m[12] * m[2] * m[7] +
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] -
              m[0]  * m[7] * m[13] -
              m[4]  * m[1] * m[15] +
              m[4]  * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] +
               m[0]  * m[6] * m[13] +
               m[4]  * m[1] * m[14] -
               m[4]  * m[2] * m[13] -
               m[12] * m[1] * m[6] +
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
              m[1] * m[7] * m[10] +
              m[5] * m[2] * m[11] -
              m[5] * m[3] * m[10] -
              m[9] * m[2] * m[7] +
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
               m[0] * m[7] * m[9] +
               m[4] * m[1] * m[11] -
               m[4] * m[3] * m[9] -
               m[8] * m[1] * m[7] +
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;
    double invOut[16];

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    for (unsigned int i = 1; i < 16; ++i){
        inverse(i / 4, i % 4) = invOut[i];
    }
    return true;
}

double CalcDeterminant(const boost::numeric::ublas::matrix<double>& m)
{
  assert(m.size1() == m.size2() && "Can only calculate the determinant of square matrices");
  switch(m.size1())
  {
    case 1:
    {
      return m(0,0);
    }
    case 2:
    {
      const double a = m(0,0);
      const double b = m(0,1);
      const double c = m(1,0);
      const double d = m(1,1);
      const double determinant = (a * d) - (b * c);
      return determinant;
    }
    case 3:
    {
      assert(m.size1() == 3 && m.size2() == 3 && "Only for 3x3 matrices");
      const double a = m(0,0);
      const double b = m(0,1);
      const double c = m(0,2);
      const double d = m(1,0);
      const double e = m(1,1);
      const double f = m(1,2);
      const double g = m(2,0);
      const double h = m(2,1);
      const double k = m(2,2);
      const double determinant
        = (a * ((e*k) - (f*h)))
        - (b * ((k*d) - (f*g)))
        + (c * ((d*h) - (e*g)));
      return determinant;
    }
    default:
      assert(!"Should not get here: unsupported matrix size");
      throw std::runtime_error("Unsupported matrix size");
  }
}

///Chop returns a std::vector of sub-matrices
//[ A at [0]   B at [1] ]
//[ C at [2]   D at [4] ]
const std::vector<boost::numeric::ublas::matrix<double> > Chop(
  const boost::numeric::ublas::matrix<double>& m)
{
  using boost::numeric::ublas::range;
  using boost::numeric::ublas::matrix;
  using boost::numeric::ublas::matrix_range;
  std::vector<matrix<double> > v;
  v.reserve(4);
  const int midy = m.size1() / 2;
  const int midx = m.size2() / 2;
  const matrix_range<const matrix<double> > top_left(    m,range(0   ,midy     ),range(0   ,midx     ));
  const matrix_range<const matrix<double> > bottom_left( m,range(midy,m.size1()),range(0   ,midx     ));
  const matrix_range<const matrix<double> > top_right(   m,range(0   ,midy     ),range(midx,m.size2()));
  const matrix_range<const matrix<double> > bottom_right(m,range(midy,m.size1()),range(midx,m.size2()));
  v.push_back(matrix<double>(top_left));
  v.push_back(matrix<double>(top_right));
  v.push_back(matrix<double>(bottom_left));
  v.push_back(matrix<double>(bottom_right));
  return v;
}

const boost::numeric::ublas::matrix<double> CreateMatrix(
  const std::size_t n_rows,
  const std::size_t n_cols,
  const std::vector<double>& v)
{
  assert(n_rows * n_cols == v.size());
  boost::numeric::ublas::matrix<double> m(n_rows,n_cols);
  for (std::size_t row = 0; row!=n_rows; ++row)
  {
    for (std::size_t col = 0; col!=n_cols; ++col)
    {
      m(row,col) = v[ (col * n_rows) + row];
    }
  }
  return m;
}

const boost::numeric::ublas::matrix<double> CreateRandomMatrix(const std::size_t n_rows, const std::size_t n_cols)
{
  boost::numeric::ublas::matrix<double> m(n_rows,n_cols);
  for (std::size_t row=0; row!=n_rows; ++row)
  {
    for (std::size_t col=0; col!=n_cols; ++col)
    {
      m(row,col) = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
    }
  }
  return m;
}

///Unchop merges the 4 std::vector of sub-matrices produced by Chop
const boost::numeric::ublas::matrix<double> Unchop(
  const std::vector<boost::numeric::ublas::matrix<double> >& v)
{
  //Chop returns a std::vector of sub-matrices
  //[ A at [0]   B at [1] ]
  //[ C at [2]   D at [4] ]
  using boost::numeric::ublas::range;
  using boost::numeric::ublas::matrix;
  using boost::numeric::ublas::matrix_range;
  assert(v.size() == 4);
  assert(v[0].size1() == v[1].size1());
  assert(v[2].size1() == v[3].size1());
  assert(v[0].size2() == v[2].size2());
  assert(v[1].size2() == v[3].size2());
  boost::numeric::ublas::matrix<double> m(v[0].size1() + v[2].size1(),v[0].size2() + v[1].size2());
  for (int quadrant=0; quadrant!=4; ++quadrant)
  {
    const boost::numeric::ublas::matrix<double>& w = v[quadrant];
    const std::size_t n_rows = v[quadrant].size1();
    const std::size_t n_cols = v[quadrant].size2();
    const int offset_x = quadrant % 2 ? v[0].size2() : 0;
    const int offset_y = quadrant / 2 ? v[0].size1() : 0;
    for (std::size_t row=0; row!=n_rows; ++row)
    {
      for (std::size_t col=0; col!=n_cols; ++col)
      {
        m(offset_y + row, offset_x + col) = w(row,col);
      }
    }
  }

  assert(v[0].size1() + v[2].size1() == m.size1());
  assert(v[1].size1() + v[3].size1() == m.size1());
  assert(v[0].size2() + v[1].size2() == m.size2());
  assert(v[2].size2() + v[3].size2() == m.size2());

  return m;
}

const boost::numeric::ublas::matrix<double> Inverse(
  const boost::numeric::ublas::matrix<double>& m)
{
  assert(m.size1() == m.size2() && "Can only calculate the inverse of square matrices");

  switch(m.size1())
  {
    case 1:
    {
      assert(m.size1() == 1 && m.size2() == 1 && "Only for 1x1 matrices");
      const double determinant = CalcDeterminant(m);
      assert(determinant != 0.0);
      assert(m(0,0) != 0.0 && "Cannot take the inverse of matrix [0]");
      boost::numeric::ublas::matrix<double> n(1,1);
      n(0,0) =  1.0 / determinant;
      return n;
    }
    case 2:
    {
      assert(m.size1() == 2 && m.size2() == 2 && "Only for 2x2 matrices");
      const double determinant = CalcDeterminant(m);
      assert(determinant != 0.0);
      const double a = m(0,0);
      const double b = m(0,1);
      const double c = m(1,0);
      const double d = m(1,1);
      boost::numeric::ublas::matrix<double> n(2,2);
      n(0,0) =  d / determinant;
      n(0,1) = -b / determinant;
      n(1,0) = -c / determinant;
      n(1,1) =  a / determinant;
      return n;
    }
    case 3:
    {
      assert(m.size1() == 3 && m.size2() == 3 && "Only for 3x3 matrices");
      const double determinant = CalcDeterminant(m);
      assert(determinant != 0.0);
      const double a = m(0,0);
      const double b = m(0,1);
      const double c = m(0,2);
      const double d = m(1,0);
      const double e = m(1,1);
      const double f = m(1,2);
      const double g = m(2,0);
      const double h = m(2,1);
      const double k = m(2,2);
      boost::numeric::ublas::matrix<double> n(3,3);
      const double new_a =  ((e*k)-(f*h)) / determinant;
      const double new_b = -((d*k)-(f*g)) / determinant;
      const double new_c =  ((d*h)-(e*g)) / determinant;
      const double new_d = -((b*k)-(c*h)) / determinant;
      const double new_e =  ((a*k)-(c*g)) / determinant;
      const double new_f = -((a*h)-(b*g)) / determinant;
      const double new_g =  ((b*f)-(c*e)) / determinant;
      const double new_h = -((a*f)-(c*d)) / determinant;
      const double new_k =  ((a*e)-(b*d)) / determinant;
      n(0,0) = new_a;
      n(1,0) = new_b;
      n(2,0) = new_c;
      n(0,1) = new_d;
      n(1,1) = new_e;
      n(2,1) = new_f;
      n(0,2) = new_g;
      n(1,2) = new_h;
      n(2,2) = new_k;
      return n;
    }
    default:
    {
      //Use blockwise inversion
      //Matrix::Chop returns a std::vector
      //[ A at [0]   B at [1] ]
      //[ C at [2]   D at [4] ]
      const std::vector<boost::numeric::ublas::matrix<double> > v = Chop(m);
      const boost::numeric::ublas::matrix<double>& a = v[0];
      assert(a.size1() == a.size2());
      const boost::numeric::ublas::matrix<double>  a_inv = Inverse(a);
      const boost::numeric::ublas::matrix<double>& b = v[1];
      const boost::numeric::ublas::matrix<double>& c = v[2];
      const boost::numeric::ublas::matrix<double>& d = v[3];
      const boost::numeric::ublas::matrix<double> term
        = d
        - prod(
            boost::numeric::ublas::matrix<double>(prod(c,a_inv)),
            b
          );
      const boost::numeric::ublas::matrix<double> term_inv = Inverse(term);
      const boost::numeric::ublas::matrix<double> new_a
        = a_inv
        + boost::numeric::ublas::matrix<double>(prod(
            boost::numeric::ublas::matrix<double>(prod(
              boost::numeric::ublas::matrix<double>(prod(
                boost::numeric::ublas::matrix<double>(prod(
                  a_inv,
                  b)),
                term_inv)),
             c)),
            a_inv));

      const boost::numeric::ublas::matrix<double> new_b
        =
        - boost::numeric::ublas::matrix<double>(prod(
            boost::numeric::ublas::matrix<double>(prod(
              a_inv,
              b)),
            term_inv));

      const boost::numeric::ublas::matrix<double> new_c
        =
        - boost::numeric::ublas::matrix<double>(prod(
            boost::numeric::ublas::matrix<double>(prod(
              term_inv,
              c)),
            a_inv));

      const boost::numeric::ublas::matrix<double> new_d = term_inv;
      std::vector<boost::numeric::ublas::matrix<double> > w;
      w.push_back(new_a);
      w.push_back(new_b);
      w.push_back(new_c);
      w.push_back(new_d);
      const boost::numeric::ublas::matrix<double> result = Unchop(w);
      return result;
    }
  }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void CalculateVelocityLaplacianRate(ModelPart& r_model_part)
{
    double delta_t_inv = 1.0 / r_model_part.GetProcessInfo()[DELTA_TIME];
    vector<unsigned int> nodes_partition;
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_model_part.Nodes().size(), nodes_partition);

    #pragma omp parallel for
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
        NodesArrayType& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();
        NodeIterator node_begin = pNodes.ptr_begin() + nodes_partition[k];
        NodeIterator node_end   = pNodes.ptr_begin() + nodes_partition[k + 1];

        for (ModelPart::NodesContainerType::iterator inode = node_begin; inode != node_end; ++inode){
            array_1d <double, 3>& laplacian_rate = inode->FastGetSolutionStepValue(VELOCITY_LAPLACIAN_RATE);
            array_1d <double, 3>& laplacian      = inode->FastGetSolutionStepValue(VELOCITY_LAPLACIAN);
            array_1d <double, 3>& old_laplacian  = inode->FastGetSolutionStepValue(VELOCITY_LAPLACIAN, 1);
            noalias(laplacian_rate) = delta_t_inv * (laplacian - old_laplacian);
        }
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

bool AssessStationarity(ModelPart& r_model_part, const double& tol)
{
    if (!mPressuresFilled){
        PerformFirstStepComputations(r_model_part);

        return(false);
    }

    else {
        double max_pressure_change_rate = 0.0; // measure of stationarity
        double mean_celerity = 0.0;            // used to adimensionalize the time step

        // filling up mPressures and calculating the mean velocities and the maximum nodal pressure change

        unsigned int i = 0;

        for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
            array_1d<double, 3> velocity = inode->FastGetSolutionStepValue(VELOCITY);
            mean_celerity += sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);

            double aux    = mPressures[i];
            mPressures[i] = inode->FastGetSolutionStepValue(PRESSURE);
            aux           = fabs(aux - mPressures[i]);

            if (aux > max_pressure_change_rate){
                max_pressure_change_rate = aux;
            }

            ++i;
        }

        mean_celerity /= i;
        double delta_t = r_model_part.GetProcessInfo()[TIME] - mLastMeasurementTime;

        if (delta_t > 0.0){
            // calculating coefficients for adimensionalization of the pressure change rate
            double pressure_variation;
            CalculateVariationWithingVector(mPressures, pressure_variation);
            double char_length         = pow(mTotalVolume, 1/3); // characteristic length of the model. Should be improved: a hydraulic radius or such
            double time_adim_coeff     = mean_celerity / char_length;
            double pressure_adim_coeff = 0.5 * (pressure_variation + mLastPressureVariation);
            mLastPressureVariation     = pressure_variation;

            if (pressure_adim_coeff == 0.0 || time_adim_coeff == 0.0){ // unlikely
                std::cout << "Uniform problem: stationarity check being performed with dimensional values...! " << "\n";

                if (max_pressure_change_rate <= tol){ // go with the absolute value
                    return true;
                }
            }

            max_pressure_change_rate /= time_adim_coeff * delta_t * pressure_adim_coeff ;
        }

        else {
            KRATOS_THROW_ERROR(std::runtime_error,"Trying to calculate pressure variations between two coincident time steps! (null time variation since last recorded time)","");
        }

        std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << "\n";
        std::cout << "The stationarity condition tolerance is " << "\n";
        KRATOS_WATCH(tol)
        std::cout << "The stationarity residual is now " << "\n";
        KRATOS_WATCH(max_pressure_change_rate)
        std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << "\n";

        if (max_pressure_change_rate <= tol){
            return true;
        }

        else {
            return false;
        }
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double CalculateDomainVolume(ModelPart& r_fluid_model_part)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_fluid_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

    double added_volume = 0.0;

    #pragma omp parallel for reduction(+ : added_volume)
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

        for (ElementIterator it = GetElementPartitionBegin(r_fluid_model_part, k); it != GetElementPartitionEnd(r_fluid_model_part, k); ++it){
            added_volume += CalculateElementalVolume(it->GetGeometry());
        }
    }

    return added_volume;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// this function assumes linear elements are used

void CalculateTotalHydrodynamicForceOnParticles(ModelPart& r_dem_model_part, array_1d <double, 3>& force)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_dem_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

    std::vector<array_1d <double, 3> > added_force_vect;
    added_force_vect.resize(OpenMPUtils::GetNumThreads());

    for (unsigned int k = 0; k < added_force_vect.size(); ++k){
        added_force_vect[k] = ZeroVector(3);
    }

    #pragma omp parallel for
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

        for (ElementIterator it = GetElementPartitionBegin(r_dem_model_part, k); it != GetElementPartitionEnd(r_dem_model_part, k); ++it){
            Geometry< Node<3> >& geom = it->GetGeometry();
            array_1d <double, 3> element_force;

            if (geom[0].SolutionStepsDataHas(HYDRODYNAMIC_FORCE)){
                element_force = geom[0].FastGetSolutionStepValue(HYDRODYNAMIC_FORCE);
            }

            else {
                element_force = ZeroVector(3);
            }

            added_force_vect[k] += element_force;
        }
    }

    force = added_force_vect[0];

    for (unsigned int k = 1; k < added_force_vect.size(); ++k){
        force += added_force_vect[k];
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// this function assumes linear elements are used

void CalculateTotalHydrodynamicForceOnFluid(ModelPart& r_fluid_model_part, array_1d <double, 3>& force)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_fluid_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

    std::vector<array_1d <double, 3> > added_force_vect;
    added_force_vect.resize(OpenMPUtils::GetNumThreads());

    for (unsigned int k = 0; k < added_force_vect.size(); ++k){
        added_force_vect[k] = ZeroVector(3);
    }

    #pragma omp parallel for
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

        for (ElementIterator it = GetElementPartitionBegin(r_fluid_model_part, k); it != GetElementPartitionEnd(r_fluid_model_part, k); ++it){
            Geometry< Node<3> >& geom = it->GetGeometry();
            double element_volume;
            array_1d <double, 3> element_force;

            if (geom[0].SolutionStepsDataHas(HYDRODYNAMIC_REACTION) && geom[0].SolutionStepsDataHas(FLUID_FRACTION)){
                element_force  = CalculateVectorIntegralOfLinearInterpolationPerUnitFluidMass(geom, HYDRODYNAMIC_REACTION, element_volume);
            }

            else {
                element_force = ZeroVector(3);
            }

            added_force_vect[k] += element_force;
        }
    }

    force = added_force_vect[0];

    for (unsigned int k = 1; k < added_force_vect.size(); ++k){
        force += added_force_vect[k];
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// this function assumes linear elements are used

double CalculateGlobalFluidVolume(ModelPart& r_fluid_model_part)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_fluid_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

    double added_fluid_volume = 0.0;

    #pragma omp parallel for reduction(+ : added_fluid_volume)
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

        for (ElementIterator it = GetElementPartitionBegin(r_fluid_model_part, k); it != GetElementPartitionEnd(r_fluid_model_part, k); ++it){
            Geometry< Node<3> >& geom = it->GetGeometry();
            double element_volume;
            double element_fluid_volume;

            if (geom[0].SolutionStepsDataHas(FLUID_FRACTION)){
                element_fluid_volume = CalculateScalarIntegralOfLinearInterpolation(geom, FLUID_FRACTION, element_volume);
            }

            else {
                element_fluid_volume = CalculateElementalVolume(geom);
            }

            added_fluid_volume += element_fluid_volume;
        }
    }

    return added_fluid_volume;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

private:

bool mPressuresFilled;
bool mFirstGradientRecovery;
double mLastMeasurementTime;
double mLastPressureVariation;
double mTotalVolume;
std::vector<double> mPressures;
std::vector<vector<double> > mFirstRowsOfB;

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

inline double CalculateVol(const double x0, const double y0, const double z0,
                           const double x1, const double y1, const double z1,
                           const double x2, const double y2, const double z2,
                           const double x3, const double y3, const double z3)
{
    double x10 = x1 - x0;
    double y10 = y1 - y0;
    double z10 = z1 - z0;

    double x20 = x2 - x0;
    double y20 = y2 - y0;
    double z20 = z2 - z0;

    double x30 = x3 - x0;
    double y30 = y3 - y0;
    double z30 = z3 - z0;

    double detJ = x10 * y20 * z30 - x10 * y30 * z20 +
                  y10 * z20 * x30 - y10 * x20 * z30 +
                  z10 * x20 * y30 - z10 * y20 * x30;

    return  detJ * 0.1666666666666666666666667;
}

//***************************************************************************************************************
//***************************************************************************************************************

double CalculateElementalVolume(const Geometry<Node <3> >& geom)
{
    double x0 = geom[0].X();
    double y0 = geom[0].Y();
    double z0 = geom[0].Z();
    double x1 = geom[1].X();
    double y1 = geom[1].Y();
    double z1 = geom[1].Z();
    double x2 = geom[2].X();
    double y2 = geom[2].Y();
    double z2 = geom[2].Z();
    double x3 = geom[3].X();
    double y3 = geom[3].Y();
    double z3 = geom[3].Z();

    double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

    if (vol == 0.0){
        KRATOS_THROW_ERROR(std::logic_error, "element with zero area found with the current geometry ", geom);
    }

    return vol;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double CalculateScalarIntegralOfLinearInterpolation(const Geometry<Node < 3 > >& geom, const Variable<double>& r_var, double& vol)
{
    array_1d<double, 4> N;
    double x0 = geom[0].X();
    double y0 = geom[0].Y();
    double z0 = geom[0].Z();
    double x1 = geom[1].X();
    double y1 = geom[1].Y();
    double z1 = geom[1].Z();
    double x2 = geom[2].X();
    double y2 = geom[2].Y();
    double z2 = geom[2].Z();
    double x3 = geom[3].X();
    double y3 = geom[3].Y();
    double z3 = geom[3].Z();

    double xc = 0.25 * (x0 + x1 + x2 + x3);
    double yc = 0.25 * (y0 + y1 + y2 + y3);
    double zc = 0.25 * (z0 + z1 + z2 + z3);

    vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

    if (vol == 0.0){
        KRATOS_THROW_ERROR(std::logic_error, "Element with zero area found. Its geometry is given by", geom);
    }

    N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc);
    N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc);
    N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc);
    N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc);

    double value_at_gauss_point = N[0] * geom[0].FastGetSolutionStepValue(r_var);

    for (unsigned int i = 1; i != 4; ++i){
        value_at_gauss_point += N[i] * geom[i].FastGetSolutionStepValue(r_var, 0);
    }

    return value_at_gauss_point;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

array_1d <double, 3> CalculateVectorIntegralOfLinearInterpolation(const Geometry<Node < 3 > >& geom, const Variable<array_1d <double, 3> >& r_var, double& vol)
{
    array_1d<double, 4> N;
    double x0 = geom[0].X();
    double y0 = geom[0].Y();
    double z0 = geom[0].Z();
    double x1 = geom[1].X();
    double y1 = geom[1].Y();
    double z1 = geom[1].Z();
    double x2 = geom[2].X();
    double y2 = geom[2].Y();
    double z2 = geom[2].Z();
    double x3 = geom[3].X();
    double y3 = geom[3].Y();
    double z3 = geom[3].Z();

    double xc = 0.25 * (x0 + x1 + x2 + x3);
    double yc = 0.25 * (y0 + y1 + y2 + y3);
    double zc = 0.25 * (z0 + z1 + z2 + z3);

    vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

    if (vol == 0.0){
        KRATOS_THROW_ERROR(std::logic_error, "Element with zero area found. Its geometry is given by", geom);
    }

    N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc);
    N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc);
    N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc);
    N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc);

    array_1d <double, 3> value_at_gauss_point = N[0] * geom[0].FastGetSolutionStepValue(r_var);

    for (unsigned int i = 1; i != 4; ++i){
        value_at_gauss_point += N[i] * geom[i].FastGetSolutionStepValue(r_var);
    }

    return value_at_gauss_point;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

array_1d <double, 3> CalculateVectorIntegralOfLinearInterpolationPerUnitFluidMass(const Geometry<Node < 3 > >& geom, const Variable<array_1d <double, 3> >& r_var, double& vol)
{
    array_1d<double, 4> N;
    double x0 = geom[0].X();
    double y0 = geom[0].Y();
    double z0 = geom[0].Z();
    double x1 = geom[1].X();
    double y1 = geom[1].Y();
    double z1 = geom[1].Z();
    double x2 = geom[2].X();
    double y2 = geom[2].Y();
    double z2 = geom[2].Z();
    double x3 = geom[3].X();
    double y3 = geom[3].Y();
    double z3 = geom[3].Z();

    double xc = 0.25 * (x0 + x1 + x2 + x3);
    double yc = 0.25 * (y0 + y1 + y2 + y3);
    double zc = 0.25 * (z0 + z1 + z2 + z3);

    vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

    if (vol == 0.0){
        KRATOS_THROW_ERROR(std::logic_error, "Element with zero area found. Its geometry is given by", geom);
    }

    N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc);
    N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc);
    N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc);
    N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc);

    array_1d <double, 3> value_at_gauss_point = N[0] * geom[0].FastGetSolutionStepValue(r_var) * geom[0].FastGetSolutionStepValue(DENSITY) * geom[0].FastGetSolutionStepValue(FLUID_FRACTION);

    for (unsigned int i = 1; i != 4; ++i){
        value_at_gauss_point += N[i] * geom[i].FastGetSolutionStepValue(r_var) * geom[i].FastGetSolutionStepValue(DENSITY) * geom[i].FastGetSolutionStepValue(FLUID_FRACTION);
    }

    return value_at_gauss_point;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************


void PerformFirstStepComputations(ModelPart& r_model_part)
{
    mTotalVolume = CalculateDomainVolume(r_model_part);
    mPressures.resize(r_model_part.Nodes().size());
    mLastMeasurementTime = r_model_part.GetProcessInfo()[TIME];

    unsigned int i = 0;

    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++) {
        mPressures[i] = inode->FastGetSolutionStepValue(PRESSURE);
        ++i;
    }

    mPressuresFilled = true;
    CalculateVariationWithingVector(mPressures, mLastPressureVariation);
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

///@}
///@name Member r_variables
///@{
vector<unsigned int> mElementsPartition;

///@}
///@name Un accessible methods
///@{

inline void CalculateVariationWithingVector(const std::vector<double>& vector, double& variation)
{
    double min = vector[0];
    double max = vector[0];

    for (unsigned int i = 0; i != vector.size(); ++i){
        min = std::min(min, mPressures[i]);
        max = std::max(max, mPressures[i]);
    }

    variation = max - min;
}

vector<unsigned int>& GetElementPartition()
{
    return mElementsPartition;
}

ElementIterator GetElementPartitionBegin(ModelPart& r_model_part, unsigned int k)
{
    return r_model_part.GetCommunicator().LocalMesh().Elements().ptr_begin() + mElementsPartition[k];
}

ElementIterator GetElementPartitionEnd(ModelPart& r_model_part, unsigned int k)
{
    return r_model_part.GetCommunicator().LocalMesh().Elements().ptr_begin() + mElementsPartition[k + 1];
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

}; // Class CustomFunctionsCalculator
} // namespace Kratos.

#endif // KRATOS_CREATE_AND_DESTROY  defined


