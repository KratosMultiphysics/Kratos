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

// System includes

// Project includes
#include "includes/model_part.h"
#include "utilities/timer.h"

//Database includes
#include "custom_utilities/discrete_particle_configure.h"
#include "discrete_particle_configure.h"
#include "includes/define.h"
#include "../../DEM_application/custom_elements/discrete_element.h"
#include "../../DEM_application/custom_elements/spheric_swimming_particle.h"
#include "includes/define.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "../../DEM_application/custom_elements/spheric_particle.h"
#include "../swimming_DEM_application.h"
#include "../../../kratos/utilities/geometry_utilities.h"

//const double prox_tol = 0.00000000001;

namespace Kratos
{
    void CalculateGeometryData2D(
            Element::GeometryType& geom,
            boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX,
            array_1d<double, 4 > & N,
            double& Area) {
        double x10 = geom[1].X() - geom[0].X();
        double y10 = geom[1].Y() - geom[0].Y();

        double x20 = geom[2].X() - geom[0].X();
        double y20 = geom[2].Y() - geom[0].Y();

        //Jacobian is calculated:
        //  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
        //J=|				|=	|			  |
        //  |dy/dxi  dy/deta|	|y1-y0   y2-y0|


        double detJ = x10 * y20 - y10 * x20;

        DN_DX(0, 0) = -y20 + y10;
        DN_DX(0, 1) = x20 - x10;
        DN_DX(1, 0) = y20;
        DN_DX(1, 1) = -x20;
        DN_DX(2, 0) = -y10;
        DN_DX(2, 1) = x10;

        DN_DX /= detJ;
        N[0] = 0.333333333333333;
        N[1] = 0.333333333333333;
        N[2] = 0.333333333333333;

        Area = 0.5 * detJ;
    }

class CustomFunctionsCalculator
{
  public:

        static const std::size_t space_dim                  = 3; ///WARNING: generalize to 2d.
        typedef DiscreteParticleConfigure<space_dim>        Configure;
        typedef Configure::ContainerType                    ParticlePointerVector;
        typedef ParticlePointerVector::iterator             ParticlePointerIterator;
        typedef Configure::IteratorType                     ParticleIterator;


    KRATOS_CLASS_POINTER_DEFINITION(CustomFunctionsCalculator);

  
    CustomFunctionsCalculator():mPressuresFilled(false){}
    /// Calculator

    virtual ~CustomFunctionsCalculator(){}

    /// Default calculator

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CalculatePressureGradient(ModelPart& r_model_part) {

            for (ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++) {
                inode->FastGetSolutionStepValue(AUX_DOUBLE_VAR) = 0.0;
                noalias(inode->FastGetSolutionStepValue(PRESSURE_GRADIENT)) = ZeroVector(3);
            }
            const std::size_t TDim = 3;
            array_1d <double, TDim + 1 > elemental_pressures;
            array_1d <double, TDim> grad;
            array_1d <double, TDim + 1 > N; //Shape functions vector//
            boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;

            for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++) {
                //compute shape function derivatives
                Geometry< Node < 3 > >& geom = ielem->GetGeometry();
                double Volume;

                if (geom.size() == 4)
                    GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
                else 
                    CalculateGeometryData2D(geom, DN_DX, N, Volume);
                    
                //get the pressure gradients;
                
                for (unsigned int i = 0; i < geom.size(); i++)
                    elemental_pressures[i] = geom[i].FastGetSolutionStepValue(PRESSURE);

                noalias(grad) = prod(trans(DN_DX), elemental_pressures);
                double nodal_area = Volume / static_cast<double>(geom.size());
                grad *= nodal_area;

                for (unsigned int i = 0; i < geom.size(); i++)
                    geom[i].FastGetSolutionStepValue(PRESSURE_GRADIENT) += grad;

                for (unsigned int i = 0; i < geom.size(); i++)
                    geom[i].FastGetSolutionStepValue(AUX_DOUBLE_VAR) += nodal_area;

            }

            for (ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++) {
                inode->FastGetSolutionStepValue(PRESSURE_GRADIENT) /= inode->FastGetSolutionStepValue(AUX_DOUBLE_VAR);
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
          double max_pressure_variation_rate = 0.0;
          array_1d<double, 3> velocity;
          double mean_celerity = 0.0;

          unsigned int i = 0;

          for (ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++) {
              double aux = mPressures[i];
              mPressures[i] = inode->FastGetSolutionStepValue(PRESSURE);
              velocity = inode->FastGetSolutionStepValue(VELOCITY);
              mean_celerity += sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);

              aux = fabs(aux - mPressures[i]);

              if (aux > max_pressure_variation_rate){
                  max_pressure_variation_rate = aux;
                }

              i++;
          }

          mean_celerity /= i;

          double delta_t = r_model_part.GetProcessInfo()[TIME] - mLastMeasurementTime;

          if (delta_t > 0.0){
              double min_pressure;
              double max_pressure;
              CalculateVariationWithingVector(mPressures, min_pressure, max_pressure);
              double pressure_variation = max_pressure - min_pressure;
              double char_length = pow(mTotalVolume, 1/3); // should be improved: a hydraulic radius or such
              double time_adim_coeff = mean_celerity / char_length;
              double pressure_adim_coeff = 0.5 * (pressure_variation + mLastPressureVariation);

              if (pressure_adim_coeff == 0.0 || time_adim_coeff == 0.0){ // unlikely

                  std::cout << "Uniform problem: stationarity check being performed with dimensional values...! " << "\n";

                  if (max_pressure_variation_rate <= tol){ // go with the absolute value
                      return(true);
                    }
                }

              max_pressure_variation_rate /= time_adim_coeff * delta_t * pressure_adim_coeff ;
            }

          else {
              KRATOS_ERROR(std::runtime_error,"Trying to calculate max pressure variation between to coincident time steps! (null time variation since last recorded time)","");
            }

          if (max_pressure_variation_rate <= tol){
              return(true);
            }

          else {
              return(false);
            }

        }

    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************


  private:

    bool mPressuresFilled;
    double mLastMeasurementTime;
    double mLastPressureVariation;
    double mTotalVolume;
    std::vector<double> mPressures;

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void PerformFirstStepComputations(ModelPart& r_model_part)
    {
      CalculateDomainVolume(r_model_part);
      mPressures.resize(r_model_part.Nodes().size());
      mLastMeasurementTime = r_model_part.GetProcessInfo()[TIME];

      unsigned int i = 0;

      for (ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++) {
          mPressures[i] = inode->FastGetSolutionStepValue(PRESSURE);
          i++;
      }

      mPressuresFilled = true;
      double min_pressure;
      double max_pressure;

      CalculateVariationWithingVector(mPressures, min_pressure, max_pressure);
      mLastPressureVariation = max_pressure - min_pressure;
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    inline double CalculateVol(	const double x0, const double y0, const double z0,
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

        return  detJ * 0.1666666666666666666667;

    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    inline void CalculateVariationWithingVector(const std::vector<double>& vector, double& min, double& max)
    {
      double min_aux = vector[0];
      double max_aux = vector[0];

      for (unsigned int i = 0; i != vector.size(); ++i){

          min_aux = std::min(min_aux, mPressures[i]);
          max_aux = std::max(max_aux, mPressures[i]);
        }

      min = min_aux;
      max = max_aux;
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CalculateDomainVolume(ModelPart& r_model_part)
    {
      const int n_elem = r_model_part.Elements().size();
      mTotalVolume = 0.0;

      for (int i = 0; i < n_elem; ++i){
          ModelPart::ElementsContainerType::iterator ielem = r_model_part.ElementsBegin() + i;
          Geometry< Node<3> >& geom = ielem->GetGeometry();
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

          mTotalVolume += fabs(CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3));
        }

    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    }; // Class CustomFunctionsCalculator


} // namespace Kratos.

#endif // KRATOS_CREATE_AND_DESTROY  defined


