//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "spheric_swimming_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application.h"


namespace Kratos
{

      SphericSwimmingParticle::SphericSwimmingParticle(): SphericParticle(){}

      SphericSwimmingParticle::SphericSwimmingParticle(IndexType NewId, GeometryType::Pointer pGeometry): SphericParticle(NewId, pGeometry){}

      SphericSwimmingParticle::SphericSwimmingParticle( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
      : SphericParticle(NewId, pGeometry, pProperties){}

      SphericSwimmingParticle::SphericSwimmingParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : SphericParticle(NewId, ThisNodes){}

      Element::Pointer SphericSwimmingParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
         return Element::Pointer(new SphericSwimmingParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
      }     

      /// Destructor.
      SphericSwimmingParticle::~SphericSwimmingParticle(){}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      void SphericSwimmingParticle::ComputeAdditionalForces(array_1d<double, 3>& contact_force,
                                                            array_1d<double, 3>& contact_moment,
                                                            array_1d<double, 3>& additionally_applied_force,
                                                            array_1d<double, 3>& additionally_applied_moment,
                                                            ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          const array_1d<double,3>& gravity         = rCurrentProcessInfo[GRAVITY];
          const double& fluid_density               = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED);
          array_1d<double,3>& buoyancy              = GetGeometry()(0)->FastGetSolutionStepValue(BUOYANCY);
          array_1d<double,3>& drag_force            = GetGeometry()(0)->FastGetSolutionStepValue(DRAG_FORCE);
          array_1d<double,3> virtual_mass_force;//    = GetGeometry()(0)->FastGetSolutionStepValue(VIRTUAL_FORCE);
          array_1d<double,3> lift_force;//            = GetGeometry()(0)->FastGetSolutionStepValue(LIFT_FORCE);

          ComputeBuoyancy(buoyancy, fluid_density, gravity, rCurrentProcessInfo);
          ComputeDragForce(drag_force, fluid_density, rCurrentProcessInfo);
          ComputeVirtualMassForce(virtual_mass_force, fluid_density, rCurrentProcessInfo);
          ComputeLiftForce(lift_force, fluid_density, rCurrentProcessInfo);

          additionally_applied_force[0] = buoyancy[0] + drag_force[0] + virtual_mass_force[0] + lift_force[0] + mRealMass * gravity[0];
          additionally_applied_force[1] = buoyancy[1] + drag_force[1] + virtual_mass_force[1] + lift_force[1] + mRealMass * gravity[1];
          additionally_applied_force[2] = buoyancy[2] + drag_force[2] + virtual_mass_force[2] + lift_force[2] + mRealMass * gravity[2];

          KRATOS_CATCH( "" )
      }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      void SphericSwimmingParticle::ComputeBuoyancy(array_1d<double, 3>& buoyancy, const double& fluid_density, const array_1d<double,3>& gravity, ProcessInfo& rCurrentProcessInfo)
      {
          // Case of identically null buoyancy

          if (mBuoyancyForceType == 0 || fluid_density < 0.00000000001 || GetGeometry()[0].IsFixed(VELOCITY_X) == true){
              noalias(buoyancy) = ZeroVector(3);
              return;
          }

          // General Case

          const array_1d<double,3>& pressure_grad = GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE_GRAD_PROJECTED);
          const double volume = 1.333333333333333 * M_PI * mRadius * mRadius * mRadius;

          noalias(buoyancy) = - volume * pressure_grad;
      }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      void SphericSwimmingParticle::ComputeDragForce(array_1d<double, 3>& drag_force, const double& fluid_density, ProcessInfo& rCurrentProcessInfo)
      {
          // Case of identically null drag force

          if (mDragForceType == 0 || fluid_density < 0.00000000001 || GetGeometry()[0].IsFixed(VELOCITY_X) == true){
              noalias(drag_force) = ZeroVector(3);
              return;
          }

          // General case

          const double fluid_fraction            = 1 - GetGeometry()(0)->FastGetSolutionStepValue(SOLID_FRACTION_PROJECTED);
          // TEMPORARY!!!: the the fluid velocity is modified here. Instead, it should be properly calculated by the fluid code
          const array_1d<double,3> fluid_vel     = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VEL_PROJECTED) / fluid_fraction;
          const array_1d<double,3>& particle_vel = GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
          const array_1d<double,3>& slip_vel     = fluid_vel - particle_vel;
          const double norm_of_slip_vel          = MathUtils<double>::Norm3(slip_vel);
          double drag_coeff;

          // calculating the 'dimensional' drag coefficient, i.e., the factor by which the slip velocity must be multiplied to yield the drag force

          if (mDragForceType == 1){
              drag_coeff = ComputeConstantDragCoefficient(norm_of_slip_vel, fluid_density, rCurrentProcessInfo);
          }

          else if (mDragForceType == 2){ // formulations of Haider (1989) and Chien (1994)
              drag_coeff = ComputeWeatherfordDragCoefficient(norm_of_slip_vel, fluid_density, rCurrentProcessInfo);
          }

          else if (mDragForceType == 3){ // formulation of Ganser (1993)
              drag_coeff = ComputeGanserDragCoefficient(norm_of_slip_vel, fluid_density, rCurrentProcessInfo);
          }

          else {
              std::cout << "The integer value designating the drag coefficient calculation model" << std::endl;
              std::cout << " (mDragForceType = " << mDragForceType << "), is not supported" << std::endl << std::flush;
              return;
          }

          noalias(drag_force) = drag_coeff * slip_vel;

      }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      void SphericSwimmingParticle::ComputeVirtualMassForce(array_1d<double, 3>& virtual_mass_force, const double& fluid_density, ProcessInfo& rCurrentProcessInfo)
      {
          // Case of identically null virtual mass force

          if (mVirtualMassForceType == 0 || fluid_density < 0.00000000001){
              noalias(virtual_mass_force) = ZeroVector(3);
              return;
          }

          // General case

      }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      void SphericSwimmingParticle::ComputeLiftForce(array_1d<double, 3>& lift_force, const double& fluid_density, ProcessInfo& rCurrentProcessInfo)
      {
          // Case of identically null lift force

          if (mLiftForceType == 0 || fluid_density < 0.00000000001){
              noalias(lift_force) = ZeroVector(3);
              return;
          }

          // General case

      }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      void SphericSwimmingParticle::ComputeReynoldsNumber(double norm_of_slip_vel, double kinematic_viscosity, double& reynolds)
      {
          reynolds = 2 * mRadius * norm_of_slip_vel / kinematic_viscosity;
      }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      void SphericSwimmingParticle::AdditionalMemberDeclarationFirstStep(const ProcessInfo& r_process_info)
      {
          mBuoyancyForceType             = r_process_info[BUOYANCY_FORCE_TYPE];
          mDragForceType                 = r_process_info[DRAG_FORCE_TYPE];
          mVirtualMassForceType          = r_process_info[VIRTUAL_MASS_FORCE_TYPE];
          mLiftForceType                 = r_process_info[LIFT_FORCE_TYPE];
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void SphericSwimmingParticle::AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
      {

          if (rVariable == REYNOLDS_NUMBER){
              const array_1d<double,3> fluid_vel      = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
              const array_1d<double,3>& particle_vel  = GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
              const double fluid_fraction_inv         = 1.0 / (1.0 - GetGeometry()(0)->FastGetSolutionStepValue(SOLID_FRACTION_PROJECTED));
              const array_1d<double,3>& avg_fluid_vel = fluid_fraction_inv * fluid_vel;
              const array_1d<double,3>& slip_vel      = avg_fluid_vel - particle_vel;
              const double norm_of_slip_vel           = MathUtils<double>::Norm3(slip_vel);
              const double kinematic_viscosity        = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
              ComputeReynoldsNumber(norm_of_slip_vel, kinematic_viscosity, Output);
          }

      }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      double SphericSwimmingParticle::ComputeConstantDragCoefficient(const double& norm_of_slip_vel, const double fluid_density, ProcessInfo& rCurrentProcessInfo)
      {
            double  drag_coeff = 0.235 * M_PI * fluid_density * mRadius * mRadius * norm_of_slip_vel;
            return drag_coeff;
      }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      double SphericSwimmingParticle::ComputeWeatherfordDragCoefficient(const double& norm_of_slip_vel, const double fluid_density, ProcessInfo& rCurrentProcessInfo)
      {
            KRATOS_TRY

              const double& particle_density             = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_DENSITY);
              const double kinematic_viscosity           = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
              const double sphericity                    = GetGeometry()(0)->GetSolutionStepValue(PARTICLE_SPHERICITY);
              //const array_1d<double,3>& buoyancy       = GetGeometry()(0)->FastGetSolutionStepValue(BUOYANCY);//S
              const array_1d<double,3>& gravity          = rCurrentProcessInfo[GRAVITY];
              const int non_newtonian_option             = rCurrentProcessInfo[NON_NEWTONIAN_OPTION];
              const int manually_imposed_drag_law_option = rCurrentProcessInfo[MANUALLY_IMPOSED_DRAG_LAW_OPTION];
              const int drag_modifier_type               = rCurrentProcessInfo[DRAG_MODIFIER_TYPE];
              const double gel_strength                  = rCurrentProcessInfo[GEL_STRENGTH];
              const double power_law_n                   = rCurrentProcessInfo[POWER_LAW_N];
              const double power_law_K                   = rCurrentProcessInfo[POWER_LAW_K];
              const double initial_drag_force            = rCurrentProcessInfo[INIT_DRAG_FORCE];
              const double drag_law_slope                = rCurrentProcessInfo[DRAG_LAW_SLOPE];
              const double power_law_tol                 = rCurrentProcessInfo[POWER_LAW_TOLERANCE];

              const double area                          = M_PI * mRadius * mRadius;
              const array_1d<double,3> weight            = mRealMass * gravity;
              const array_1d<double,3> buoyancy          = fluid_density / particle_density * weight; // hydrostatic case!! (only for Weatherford)
              double shahs_term_vel                      = 0.0;
              double beta                                = 0.0;
              double F0                                  = 0.0;
              double regularization_v                    = 0.02 * mRadius;
              double reynolds;
              double drag_coeff;

              if (!non_newtonian_option){ //Newtonian
                  ComputeReynoldsNumber(norm_of_slip_vel, kinematic_viscosity, reynolds);

                  if (!non_newtonian_option && reynolds < 0.01){
                      reynolds = 0.01;
                  }

                  CalculateNewtonianDragCoefficient(non_newtonian_option, reynolds, sphericity, drag_coeff, drag_modifier_type);
                  drag_coeff = 0.5 * fluid_density * area * drag_coeff * norm_of_slip_vel;
              }

              else {
                  shahs_term_vel = CalculateShahsTerm(power_law_n, power_law_K, power_law_tol, fluid_density, particle_density, sphericity, drag_modifier_type);

                  if (!manually_imposed_drag_law_option){
                      F0 = 4.0 * gel_strength * area; //initial value
                      beta = (MathUtils<double>::Norm3(weight) - MathUtils<double>::Norm3(buoyancy) - F0) / shahs_term_vel; //S
                  }

                  else {
                      F0 = initial_drag_force; //initial value
                      beta = drag_law_slope; //slope
                  }

                  if (norm_of_slip_vel >= regularization_v){
                      drag_coeff = (F0 + beta * norm_of_slip_vel) / norm_of_slip_vel;
                  }

                  else {
                      drag_coeff = (F0 + beta * regularization_v) / regularization_v;
                  }

              }

              return drag_coeff;

          KRATOS_CATCH("")
      }
        
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SphericSwimmingParticle::CalculateNewtonianDragCoefficient(int non_newtonian_option,
                                                                    const double reynolds,
                                                                    double sphericity,
                                                                    double& drag_coeff,
                                                                    int drag_modifier_type)
    {
        if (reynolds < 1){
            drag_coeff = 24.0; // Reynolds;
        }

        else {

            if (reynolds > 1000){
                drag_coeff = 0.44;
            }

            else{
                drag_coeff = 24.0 / reynolds * (1.0 + 0.15 * pow(reynolds, 0.687));
            }

        }

        if (!non_newtonian_option){ // Newtonian

            if (drag_coeff > 2.0){
                drag_coeff = 2.0; // watch out!
            }

        }

        if (sphericity < 0.9999){
            drag_coeff = CalculateDragCoeffFromSphericity(reynolds, sphericity, drag_modifier_type);
        }

    }

  //**************************************************************************************************************************************************
  //**************************************************************************************************************************************************

    double SphericSwimmingParticle::CalculateDragCoeffFromSphericity(const double reynolds,
                                                                     const double sphericity,
                                                                     const int drag_modifier_type)
    {
        double cdrag = 0.0;

        if (drag_modifier_type == 1){ // visual-Red Book
            double interpolator = (1 - sphericity) / (1 - 0.806);
            double cdrag_modifier = 1 + 0.97 * interpolator + 0.715 * interpolator * log10(reynolds);

            if (reynolds < 1){
                cdrag_modifier += 0.3 * interpolator * pow(- 1.0 * log10(reynolds), 1.6);
            }

            cdrag = cdrag_modifier * cdrag;
        }

        if (drag_modifier_type == 2){ // Hayder
            cdrag = 24 / reynolds * (1 + exp(2.3288 - 6.4581 * sphericity + 2.4486 * sphericity * sphericity) * pow(reynolds, 0.0964 + 0.5565 * sphericity)) + 73.69 * reynolds * exp(- 5.0748 * sphericity) / (reynolds + 5.378 * exp(6.2122 * sphericity));
        }

        if (drag_modifier_type == 3){ // Chien
            cdrag = 30 / reynolds + 67.289 * exp(- 5.03 * sphericity);
        }

        return cdrag;
    }

  //**************************************************************************************************************************************************
  //**************************************************************************************************************************************************

    double SphericSwimmingParticle::CalculateShahsTerm(double power_law_N,
                                                       double power_law_K,
                                                       double power_law_tol,
                                                       const double& fluid_density,
                                                       const double& particle_density,
                                                       double sphericity,
                                                       int drag_modifier_type)
    {
      if (fabs(power_law_N) < power_law_tol || fabs(power_law_K) < power_law_tol){
          std::cout << "WARNING: Shah's method is being used with Power Law data being zero!!" << std::endl << std::flush;
      }
      double shah_A_i = 1 / (6.9148 * power_law_N * power_law_N - 24.838 * power_law_N + 22.642);
      double shah_B_i = 1 / (-0.5067 * power_law_N * power_law_N + 1.3234 * power_law_N - 0.1744);

      double dimensionless_shah = sqrt(pow(13.08, 2 - power_law_N) * pow(2 * mRadius, power_law_N + 2) * pow(fluid_density, power_law_N) * pow(particle_density - fluid_density, 2 - power_law_N) / (pow(2, 2 * (power_law_N - 1)) * power_law_K * power_law_K));
      double reynolds = pow(dimensionless_shah * shah_A_i, shah_B_i);
      double fi_i = CalculateDragCoeffFromSphericity(reynolds, 1.0, drag_modifier_type) / CalculateDragCoeffFromSphericity(reynolds, sphericity, drag_modifier_type);
      dimensionless_shah = sqrt(pow(fi_i, 2 - power_law_N)) * dimensionless_shah;
      reynolds = pow(dimensionless_shah * shah_A_i, shah_B_i);

      double terminal_vel =  pow(pow(2, power_law_N - 1) * power_law_K * reynolds / (pow(2 * mRadius, power_law_N) * fluid_density), 1 / (2 - power_law_N)) ;

      return terminal_vel;

    }

   //**************************************************************************************************************************************************
   //**************************************************************************************************************************************************

     double SphericSwimmingParticle::ComputeGanserDragCoefficient(const double& norm_of_slip_vel, const double fluid_density, ProcessInfo& rCurrentProcessInfo)
     {
       KRATOS_TRY

       const double kinematic_viscosity         = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
       const double sphericity                  = GetGeometry()(0)->GetSolutionStepValue(PARTICLE_SPHERICITY);
       const int isometric_shape                = 1; // TEMPORARY!! yes (1) or no (0); shold be given as data
       const double surface_area                = 4 * M_PI * mRadius * mRadius; // TEMPORARY!! corresponding to a sphere; should be generalized b taking it as a parameter
       const double surface_area_circular_diam  = sqrt(4.0 * surface_area / M_PI);

       double equiv_reynolds;
       double k_1;
       double k_2;
       double drag_coeff;

       ComputeGanserParameters(isometric_shape, sphericity, surface_area_circular_diam, k_1, k_2);
       ComputeReynoldsNumber(norm_of_slip_vel, kinematic_viscosity, equiv_reynolds);
       equiv_reynolds *= k_1 * k_2;

       // calculating adimensional drag coefficient

       drag_coeff =  k_2 * (24 * (1 + 0.1118 * pow((equiv_reynolds), 0.6567)) / (equiv_reynolds) + 0.4305 / (1 + 3305 / equiv_reynolds));

       // and then the dimensional drag coefficient

       drag_coeff *= 0.5 * fluid_density * surface_area * norm_of_slip_vel;


       return drag_coeff;

     KRATOS_CATCH("")
     }

   //**************************************************************************************************************************************************
   //**************************************************************************************************************************************************

     void SphericSwimmingParticle::ComputeGanserParameters(const int isometric_shape, const double sphericity, const double dn, double& k_1, double& k_2)
     {

       if (isometric_shape){
           k_1 = (3 / (1 + 2 / sqrt(sphericity)));
       }

       else {
           k_1 = (3 / (0.5 * dn / mRadius + 2 / sqrt(sphericity)));
       }

       k_2 = pow(10.0, 1.8148 * pow(- log10(sphericity), 0.5743));

     }

   //**************************************************************************************************************************************************
   //**************************************************************************************************************************************************


}  // namespace Kratos.
