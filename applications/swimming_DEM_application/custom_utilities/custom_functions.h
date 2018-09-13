//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas (gcasas@cimmne.upc.edu) $
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
#include "custom_utilities/search/discrete_particle_configure.h"
#include "includes/define.h"
#include "../../DEM_application/custom_elements/discrete_element.h"
#include "custom_elements/spheric_swimming_particle.h"
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

CustomFunctionsCalculator(): mPressuresFilled(false), mFirstGradientRecovery(true), mFirstLaplacianRecovery(true), mSomeCloudsDontWork(false), mCalculatingTheGradient(false), mCalculatingTheLaplacian(false), mFirstTimeAppending(true){}
/// Calculator

virtual ~CustomFunctionsCalculator(){}

/// Default calculator

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void CalculatePressureGradient(ModelPart& r_model_part)
{
    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(PRESSURE_GRADIENT)) = ZeroVector(3);
    }

    array_1d <double, 3> grad = ZeroVector(3); // its dimension is always 3
    array_1d <double, TDim + 1 > elemental_pressures;
    array_1d <double, TDim + 1 > N; // shape functions vector
    BoundedMatrix<double, TDim + 1, TDim> DN_DX;

    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem){
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

    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
        inode->FastGetSolutionStepValue(PRESSURE_GRADIENT) /= inode->FastGetSolutionStepValue(NODAL_AREA);
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// This function assesses the stationarity based on the pressure field varaition.
// Its tolerance applies to the adimensinalised pressure variation between consecutive
// measurements.
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

        for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode){
            const array_1d<double, 3>& velocity = inode->FastGetSolutionStepValue(VELOCITY);
            mean_celerity += SWIMMING_MODULUS_3(velocity);

            const double new_pressure = inode->FastGetSolutionStepValue(PRESSURE);
            double& old_pressure = mPressures[i];
            const double delta_p = std::abs(new_pressure - old_pressure);
            max_pressure_change_rate = std::max(delta_p, max_pressure_change_rate);
            old_pressure = new_pressure;

            ++i;
        }

        mean_celerity /= i;
        const double delta_t = r_model_part.GetProcessInfo()[TIME] - mLastMeasurementTime;

        if (delta_t > 0.0){
            max_pressure_change_rate /= delta_t;

            // calculating coefficients for adimensionalization of the pressure change rate
            const double characteristic_length             = std::pow(mTotalDomainVolume, 1.0 / 3); // characteristic length of the model. Should be improved: a hydraulic radius or such
            const double reciprocal_of_characteristic_time = mean_celerity / characteristic_length;
            const double pressure_spatial_variation = GetRangeWithinVector(mPressures);
            mLastPressureVariation = pressure_spatial_variation;
            const double characteristic_pressure_variation = 0.5 * (pressure_spatial_variation + mLastPressureVariation);

            if (characteristic_pressure_variation == 0.0 || reciprocal_of_characteristic_time == 0.0){ // unlikely
                std::cout << "Uniform problem: stationarity check being performed with dimensional values...! " << "\n";

                if (max_pressure_change_rate <= tol){ // go with the absolute value
                    return true;
                }
            }

            max_pressure_change_rate /= reciprocal_of_characteristic_time * characteristic_pressure_variation ;
        }

        else {
            KRATOS_THROW_ERROR(std::runtime_error, "Trying to calculate pressure variations between two coincident time steps! (null time variation since last recorded time)","");
        }

        std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << "\n";
        std::cout << "The stationarity condition tolerance is " << "\n";
        KRATOS_WATCH(tol)
        std::cout << "The stationarity residual is now " << "\n";
        KRATOS_WATCH(max_pressure_change_rate)
        std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << "\n";

        return max_pressure_change_rate <= tol;
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
double CalculateDomainVolume(ModelPart& r_fluid_model_part)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_fluid_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

    double added_volume = 0.0;

    #pragma omp parallel for reduction(+ : added_volume)
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); ++k){

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
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); ++k){

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

void CalculateTotalHydrodynamicForceOnFluid(ModelPart& r_fluid_model_part, array_1d <double, 3>& instantaneous_force, array_1d <double, 3>& mean_force)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_fluid_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

    std::vector<array_1d <double, 3> > added_force_vect;
    added_force_vect.resize(OpenMPUtils::GetNumThreads());
    std::vector<array_1d <double, 3> > added_mean_force_vect;
    added_mean_force_vect.resize(OpenMPUtils::GetNumThreads());

    for (unsigned int k = 0; k < added_force_vect.size(); ++k){
        added_force_vect[k] = ZeroVector(3);
        added_mean_force_vect[k] = ZeroVector(3);
    }

    #pragma omp parallel for
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); ++k){

        for (ElementIterator it = GetElementPartitionBegin(r_fluid_model_part, k); it != GetElementPartitionEnd(r_fluid_model_part, k); ++it){
            Geometry< Node<3> >& geom = it->GetGeometry();
            double element_volume;
            array_1d <double, 3> element_force;
            array_1d <double, 3> element_mean_force;

            if (geom[0].SolutionStepsDataHas(HYDRODYNAMIC_REACTION) && geom[0].SolutionStepsDataHas(FLUID_FRACTION)){
                element_force  = CalculateVectorIntegralOfLinearInterpolationPerUnitFluidMass(geom, HYDRODYNAMIC_REACTION, element_volume);
            }

            else {
                element_force = ZeroVector(3);
            }

            if (geom[0].SolutionStepsDataHas(MEAN_HYDRODYNAMIC_REACTION) && geom[0].SolutionStepsDataHas(FLUID_FRACTION)){
                element_mean_force  = CalculateVectorIntegralOfLinearInterpolationPerUnitFluidMass(geom, MEAN_HYDRODYNAMIC_REACTION, element_volume);
            }

            else {
                element_mean_force = ZeroVector(3);
            }

            added_force_vect[k] += element_force;
            added_mean_force_vect[k] += element_mean_force;
        }
    }

    instantaneous_force = added_force_vect[0];
    mean_force          = added_force_vect[0];

    for (unsigned int k = 1; k < added_force_vect.size(); ++k){
        instantaneous_force += added_force_vect[k];
        mean_force          += added_mean_force_vect[k];
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
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); ++k){

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
const DenseMatrix<double> Inverse(
  const DenseMatrix<double>& m)
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
      DenseMatrix<double> n(1,1);
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
      DenseMatrix<double> n(2,2);
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
      DenseMatrix<double> n(3,3);
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
      const std::vector<DenseMatrix<double> > v = Chop(m);
      const DenseMatrix<double>& a = v[0];
      assert(a.size1() == a.size2());
      const DenseMatrix<double>  a_inv = Inverse(a);
      const DenseMatrix<double>& b = v[1];
      const DenseMatrix<double>& c = v[2];
      const DenseMatrix<double>& d = v[3];
      const DenseMatrix<double> term
        = d
        - prod(
            DenseMatrix<double>(prod(c,a_inv)),
            b
          );
      const DenseMatrix<double> term_inv = Inverse(term);
      const DenseMatrix<double> new_a
        = a_inv
        + DenseMatrix<double>(prod(
            DenseMatrix<double>(prod(
              DenseMatrix<double>(prod(
                DenseMatrix<double>(prod(
                  a_inv,
                  b)),
                term_inv)),
             c)),
            a_inv));

      const DenseMatrix<double> new_b
        =
        - DenseMatrix<double>(prod(
            DenseMatrix<double>(prod(
              a_inv,
              b)),
            term_inv));

      const DenseMatrix<double> new_c
        =
        - DenseMatrix<double>(prod(
            DenseMatrix<double>(prod(
              term_inv,
              c)),
            a_inv));

      const DenseMatrix<double> new_d = term_inv;
      std::vector<DenseMatrix<double> > w;
      w.push_back(new_a);
      w.push_back(new_b);
      w.push_back(new_c);
      w.push_back(new_d);
      const DenseMatrix<double> result = Unchop(w);
      return result;
    }
  }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
void CopyValuesFromFirstToSecond(ModelPart& r_model_part, const Variable<double>& origin_variable, const Variable<double>& destination_variable)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        ModelPart::NodesContainerType::iterator i_particle = r_model_part.NodesBegin() + i;
        Node<3>::Pointer p_node = *(i_particle.base());
        double& destination_value = p_node->FastGetSolutionStepValue(destination_variable);
        const double& origin_value = p_node->FastGetSolutionStepValue(origin_variable);
        destination_value = origin_value;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
void CopyValuesFromFirstToSecond(ModelPart& r_model_part, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >& origin_variable, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >& destination_variable)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        ModelPart::NodesContainerType::iterator i_particle = r_model_part.NodesBegin() + i;
        Node<3>::Pointer p_node = *(i_particle.base());
        double& destination_value = p_node->FastGetSolutionStepValue(destination_variable);
        const double origin_value = p_node->FastGetSolutionStepValue(origin_variable);
        destination_value = origin_value;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
void CopyValuesFromFirstToSecond(ModelPart& r_model_part, const Variable<array_1d<double, 3>>& origin_variable, const Variable<array_1d<double, 3>>& destination_variable)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        ModelPart::NodesContainerType::iterator i_particle = r_model_part.NodesBegin() + i;
        Node<3>::Pointer p_node = *(i_particle.base());
        array_1d<double, 3>& destination_value = p_node->FastGetSolutionStepValue(destination_variable);
        const array_1d<double, 3>& origin_value = p_node->FastGetSolutionStepValue(origin_variable);
        noalias(destination_value) = origin_value;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
void SetValueOfAllNotes(ModelPart& r_model_part, const double& value, const Variable<double>& destination_variable)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        ModelPart::NodesContainerType::iterator i_particle = r_model_part.NodesBegin() + i;
        Node<3>::Pointer p_node = *(i_particle.base());
        double& destination_value = p_node->FastGetSolutionStepValue(destination_variable);
        destination_value = value;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
void SetValueOfAllNotes(ModelPart& r_model_part, const array_1d<double, 3>& value, const Variable<array_1d<double, 3>>& destination_variable)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        ModelPart::NodesContainerType::iterator i_particle = r_model_part.NodesBegin() + i;
        Node<3>::Pointer p_node = *(i_particle.base());
        array_1d<double, 3>& destination_value = p_node->FastGetSolutionStepValue(destination_variable);
        noalias(destination_value) = value;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

private:

bool mPressuresFilled;
bool mFirstGradientRecovery;
bool mFirstLaplacianRecovery;
bool mSomeCloudsDontWork;
bool mCalculatingTheGradient;
bool mCalculatingTheLaplacian;
bool mFirstTimeAppending;
double mLastMeasurementTime;
double mLastPressureVariation;
double mTotalDomainVolume;
std::vector<double> mPressures;
std::vector<DenseVector<double> > mFirstRowsOfB;

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

inline double CalculateArea(const double x0, const double y0,
                            const double x1, const double y1,
                            const double x2, const double y2)
{
    const double x10 = x1 - x0;
    const double y10 = y1 - y0;

    const double x20 = x2 - x0;
    const double y20 = y2 - y0;

    const double area = 0.5 * std::abs(x10 * y20 - x20 * y10);

    return area;
}
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
    double vol;

    if (TDim == 2){
        double x0 = geom[0].X();
        double y0 = geom[0].Y();
        double x1 = geom[1].X();
        double y1 = geom[1].Y();
        double x2 = geom[2].X();
        double y2 = geom[2].Y();
        vol = CalculateArea(x0, y0, x1, y1, x2, y2);
    }

    else {
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

        vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);
    }

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
    mTotalDomainVolume = CalculateDomainVolume(r_model_part);
    mPressures.resize(r_model_part.Nodes().size());
    mLastMeasurementTime = r_model_part.GetProcessInfo()[TIME];

    unsigned int i = 0;

    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); ++inode) {
        mPressures[i] = inode->FastGetSolutionStepValue(PRESSURE);
        ++i;
    }

    mPressuresFilled = true;
    mLastPressureVariation = GetRangeWithinVector(mPressures);
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

struct IsCloser{
    bool operator()(std::pair<unsigned int, double> const& first_pair, std::pair<unsigned int, double> const& second_pair)
    {
        return(first_pair.second < second_pair.second || (first_pair.second == second_pair.second && first_pair.first < second_pair.first));
    }
};

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

inline int Factorial(const unsigned int n){

    if (n == 0){
        return 1;
    }

    unsigned int k = n;

    for (unsigned int i = n - 1; i > 0; --i){
        k *= i;
    }

    return k;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double CalculateTheMaximumEdgeLength(ModelPart& r_model_part)
{
    double max_distance_yet = 0.0;

    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem){
        Geometry<Node<3> >& geom = ielem->GetGeometry();
        unsigned int n_nodes = static_cast<unsigned int>(TDim + 1);

        for (unsigned int k = 1; k < n_nodes - 1; ++k){
            for (unsigned int i = k; i < n_nodes; ++i){
                array_1d <double, 3> delta_i = geom[k - 1] - geom[i];
                double distance_2 = DEM_INNER_PRODUCT_3(delta_i, delta_i);
                max_distance_yet = max_distance_yet > distance_2 ? max_distance_yet : distance_2;
            }
        }
    }

    return(std::sqrt(max_distance_yet));
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double CalculateTheMinumumEdgeLength(ModelPart& r_model_part)
{
    double min_distance_yet = 0.0;

    bool first_node = true;

    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ++ielem){
        Geometry<Node<3> >& geom = ielem->GetGeometry();

        if (first_node){ // assign the distance (squared) between any two nodes to min_distance_yet
            array_1d <double, 3> delta = geom[0] - geom[1];
            double distance_2 = DEM_INNER_PRODUCT_3(delta, delta);
            min_distance_yet = distance_2;
        }

        unsigned int n_nodes = static_cast<unsigned int>(TDim + 1);

        for (unsigned int k = 1; k < n_nodes - 1; ++k){
            for (unsigned int i = k; i < n_nodes; ++i){
                array_1d <double, 3> delta_i = geom[k - 1] - geom[i];
                double distance_2 = DEM_INNER_PRODUCT_3(delta_i, delta_i);

                min_distance_yet = min_distance_yet < distance_2 ? min_distance_yet : distance_2;
            }
        }
    }

    return(std::sqrt(min_distance_yet));
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

// The following block of functions is used to calculate explicit matrix inverses and was taken from
// Richel BilderBeek's website (http://www.richelbilderbeek.nl/CppUblasMatrixExample6.htm), and it is
// transcribed here with a very minor modification

double CalcDeterminant(const DenseMatrix<double>& m)
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
const std::vector<DenseMatrix<double> > Chop(
  const DenseMatrix<double>& m)
{
  using boost::numeric::ublas::range;
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

///Unchop merges the 4 std::vector of sub-matrices produced by Chop
const DenseMatrix<double> Unchop(
  const std::vector<DenseMatrix<double> >& v)
{
  //Chop returns a std::vector of sub-matrices
  //[ A at [0]   B at [1] ]
  //[ C at [2]   D at [4] ]
  using boost::numeric::ublas::range;
  using boost::numeric::ublas::matrix_range;
  assert(v.size() == 4);
  assert(v[0].size1() == v[1].size1());
  assert(v[2].size1() == v[3].size1());
  assert(v[0].size2() == v[2].size2());
  assert(v[1].size2() == v[3].size2());
  DenseMatrix<double> m(v[0].size1() + v[2].size1(),v[0].size2() + v[1].size2());
  for (int quadrant=0; quadrant!=4; ++quadrant)
  {
    const DenseMatrix<double>& w = v[quadrant];
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

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

///@}
///@name Member r_variables
///@{
DenseVector<unsigned int> mElementsPartition;

///@}
///@name Un accessible methods
///@{

double GetRangeWithinVector(const std::vector<double>& vector)
{
    double min = vector[0];
    double max = vector[0];

    for (unsigned int i = 0; i != vector.size(); ++i){
        min = std::min(min, mPressures[i]);
        max = std::max(max, mPressures[i]);
    }

    return (max - min);
}

DenseVector<unsigned int>& GetElementPartition()
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


