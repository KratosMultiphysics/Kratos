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

      void SphericSwimmingParticle::ComputeAdditionalForces(array_1d<double, 3>& contact_force, array_1d<double, 3>& contact_moment,
                                                            array_1d<double, 3>& additionally_applied_force, array_1d<double, 3>& additionally_applied_moment, ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          const int drag_force_type                 = rCurrentProcessInfo[DRAG_FORCE_TYPE];
          const array_1d<double,3>& gravity         = rCurrentProcessInfo[GRAVITY];
          array_1d<double,3>& drag_force            = GetGeometry()(0)->FastGetSolutionStepValue(DRAG_FORCE);
          array_1d<double,3>& buoyancy              = GetGeometry()(0)->FastGetSolutionStepValue(BUOYANCY);

          /*contact_force[0]  = 0.0;
          contact_force[1]  = 0.0;
          contact_force[2]  = 0.0;
          contact_moment[0] = 0.0;
          contact_moment[1] = 0.0;
          contact_moment[2] = 0.0;
          initial_rotation_moment[0]  = 0.0;
          initial_rotation_moment[1]  = 0.0;
          initial_rotation_moment[2]  = 0.0;
	  
          ComputeNewNeighboursHistoricalData();
          ComputeBallToBallContactForce(contact_force, contact_moment, elastic_force, initial_rotation_moment, rCurrentProcessInfo);
          
          if (mLimitSurfaceOption > 0){
			  for (int surface_num = 0; surface_num < mLimitSurfaceOption; surface_num++){
				  ComputeBallToSurfaceContactForce(contact_force, contact_moment, initial_rotation_moment, surface_num, rCurrentProcessInfo);
              }
          }
          
          if (mLimitCylinderOption > 0){
			  for (int cylinder_num = 0; cylinder_num < mLimitCylinderOption; cylinder_num++){
				  ComputeBallToCylinderContactForce(contact_force, contact_moment, initial_rotation_moment, cylinder_num, rCurrentProcessInfo);
              }
	      }*/

          if (drag_force_type == 1){
              ComputeFluidForcesOnParticle(rCurrentProcessInfo);
          }

          else {
              ComputeWeatherfordFluidForcesOnParticle(rCurrentProcessInfo);
          }
          
          additionally_applied_force[0] = buoyancy[0] + drag_force[0] + mRealMass * gravity[0];
          additionally_applied_force[1] = buoyancy[1] + drag_force[1] + mRealMass * gravity[1];
          additionally_applied_force[2] = buoyancy[2] + drag_force[2] + mRealMass * gravity[2];

          KRATOS_CATCH( "" )
      }


    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      void SphericSwimmingParticle::ComputeFluidForcesOnParticle(ProcessInfo& rCurrentProcessInfo)
      {
            KRATOS_TRY

            array_1d<double,3>& pressure_grad       = GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE_GRAD_PROJECTED);
            array_1d<double,3>& drag_force          = GetGeometry()(0)->FastGetSolutionStepValue(DRAG_FORCE);
            array_1d<double,3>& buoyancy            = GetGeometry()(0)->FastGetSolutionStepValue(BUOYANCY);
            double fluid_fraction                   = 1 - GetGeometry()(0)->FastGetSolutionStepValue(SOLID_FRACTION_PROJECTED);

            if (GetGeometry()[0].IsFixed(VELOCITY_X) == false){
                const double& fluid_density = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED);
                double volume = (1.333333333333333 * M_PI) * (mRadius * mRadius * mRadius);
                const array_1d<double,3>& fluid_vel = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
                const array_1d<double,3>& particle_vel = GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);

                if (fluid_density > 0.0000000000001){
                    drag_force = (0.235 * M_PI) * (fluid_density * mRadius * mRadius) * MathUtils<double>::Norm3(fluid_vel - particle_vel) * (fluid_vel - particle_vel);
		    //KRATOS_WATCH(drag_force) //SALVA
                }

                else {
                    noalias(drag_force)    = ZeroVector(3);
                    noalias(pressure_grad) = ZeroVector(3);
                }

                buoyancy = - pressure_grad * volume;
                //KRATOS_WATCH(buoyancy) //SALVA
            }

            else {
               noalias(drag_force)    = ZeroVector(3);
               noalias(pressure_grad) = ZeroVector(3);
               noalias(buoyancy)      = ZeroVector(3);
            }

            KRATOS_CATCH("")
      }


      double SphericSwimmingParticle::CalculateDragCoeffFromSphericity(const double Reynolds, const double Sphericity, int DragModifierType)
      {
          double cdrag = 0.0;

          if (DragModifierType == 1){ // visual-Red Book
              double interpolator = (1 - Sphericity) / (1 - 0.806);
              double cd_modifier = 1 + 0.97 * interpolator + 0.715 * interpolator * log10(Reynolds);

              if (Reynolds < 1){
                  cd_modifier += 0.3 * interpolator * pow(- 1.0 * log10(Reynolds), 1.6);
              }

              cdrag = cd_modifier * cdrag;
          }

          if (DragModifierType == 2){ // Hayder
              cdrag = 24 / Reynolds * (1 + exp(2.3288 - 6.4581 * Sphericity + 2.4486 * Sphericity * Sphericity) * pow(Reynolds, 0.0964 + 0.5565 * Sphericity)) + 73.69 * Reynolds * exp(- 5.0748 * Sphericity) / (Reynolds + 5.378 * exp(6.2122 * Sphericity));
          }

          if (DragModifierType == 3){ // Chien
              cdrag = 30 / Reynolds + 67.289 * exp(- 5.03 * Sphericity);
          }

          return cdrag;
      }



      void SphericSwimmingParticle::CalculateDragCoefficient(int NonNewtonianOption, const double Reynolds, double Sphericity, double& rDrag_coeff, int DragModifierType)
      {
          if (Reynolds < 1){
              rDrag_coeff = 24.0; ///reynolds;
          }

          else {

              if (Reynolds > 1000){
                  rDrag_coeff = 0.44;
              }

              else{
                  rDrag_coeff = 24.0 / Reynolds * (1.0 + 0.15 * pow(Reynolds, 0.687));
              }
          }

          if (!NonNewtonianOption){ //Newtonian

              if (rDrag_coeff > 2.0){
                  rDrag_coeff = 2.0; ///CUIDADO!!!
              }
          }

          if (Sphericity < 0.9999){
              rDrag_coeff = CalculateDragCoeffFromSphericity(Reynolds, Sphericity, DragModifierType);
          }

      }


      void SphericSwimmingParticle::ComputeReynoldsNumber(int NonNewtonianOption, double rNormOfSlipVel, double FluidDensity, double rViscosity, double& rReynolds)
      {
          rReynolds = 2 * mRadius * FluidDensity * rNormOfSlipVel / rViscosity; //dynamic viscosity was received

          if (!NonNewtonianOption && rReynolds < 0.01){
              rReynolds = 0.01;
          }

      }


      double SphericSwimmingParticle::CalculateShahsTerm(double PowerLawN,double PowerLawK, double PowerLawTol, const double& ParticleDensity, const double& FluidDensity, double Sphericity, int DragModifierType)
      {
          if (fabs(PowerLawN) < PowerLawTol || fabs(PowerLawK) < PowerLawTol){
              std::cout << "WARNING: Shah's method is being used with Power Law data being zero!!"<< std::endl << std::flush;
          }

          double shah_A_i = 1 / (6.9148 * PowerLawN * (PowerLawN - 24.838) + 22.642);
          double shah_B_i = 1 / (- 0.5067 * PowerLawN * (PowerLawN + 1.3234) - 0.1744);

          double dimensionless_shah = sqrt(pow(13.08, 2 - PowerLawN) * pow(2 * mRadius, PowerLawN + 2) * pow(FluidDensity, PowerLawN) * pow(ParticleDensity - FluidDensity, 2 - PowerLawN) / (pow(2, 2 * (PowerLawN - 1)) * PowerLawN * PowerLawN));
          double reynolds = pow(dimensionless_shah * shah_A_i, shah_B_i);
          double fi_i = CalculateDragCoeffFromSphericity(reynolds, 1.0, DragModifierType) / CalculateDragCoeffFromSphericity(reynolds, Sphericity, DragModifierType);
          dimensionless_shah = sqrt(pow(fi_i, 2 - PowerLawN)) * dimensionless_shah;
          reynolds = pow(dimensionless_shah * shah_A_i, shah_B_i);

          double terminal_vel =  pow(pow(2, PowerLawN - 1) * PowerLawK * reynolds / (pow(2 * mRadius, PowerLawN) * FluidDensity), 1 / (2 - PowerLawN)) ;

          return terminal_vel;
      }

      void SphericSwimmingParticle::ComputeWeatherfordFluidForcesOnParticle(ProcessInfo& rCurrentProcessInfo)
      {
            KRATOS_TRY

            array_1d<double,3>& pressure_grad       = GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE_GRAD_PROJECTED);
            array_1d<double,3>& drag_force          = GetGeometry()(0)->FastGetSolutionStepValue(DRAG_FORCE);
            array_1d<double,3>& buoyancy            = GetGeometry()(0)->FastGetSolutionStepValue(BUOYANCY);

            const double& particle_density          = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_DENSITY);
            double viscosity                        = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
            double sphericity                       = GetGeometry()(0)->GetSolutionStepValue(PARTICLE_SPHERICITY);

            int non_newtonian_OPTION                = 0; //rCurrentProcessInfo[NON_NEWTONIAN_OPTION];
            int manually_imposed_drag_law_OPTION    = 0; //rCurrentProcessInfo[MANUALLY_IMPOSED_DRAG_LAW_OPTION];
            int drag_modifier_type                  = rCurrentProcessInfo[DRAG_MODIFIER_TYPE]; //2 for Hayder and 3 for CHIEN
            double gel_strength                     = 1.0; //rCurrentProcessInfo[GEL_STRENGTH];
            double power_law_n                      = 0.02; //rCurrentProcessInfo[POWER_LAW_N];
            double power_law_K                      = 0.02; //rCurrentProcessInfo[POWER_LAW_K];
            double initial_drag_force               = 0.5; //rCurrentProcessInfo[INIT_DRAG_FORCE];
            double drag_law_slope                   = 0.05; //rCurrentProcessInfo[DRAG_LAW_SLOPE];
            const double power_law_tol              = 0.000000001; //rCurrentProcessInfo[POWER_LAW_TOL];
            const array_1d<double,3>& gravity       = rCurrentProcessInfo[GRAVITY];

            if (GetGeometry()[0].IsFixed(VELOCITY_X) == false){
                const array_1d<double,3>& fluid_vel    = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
                const array_1d<double,3>& particle_vel = GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
                const double& fluid_density            = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED);
                double volume                          = (1.333333333333333 * M_PI) * (mRadius * mRadius * mRadius);
                double dynamic_viscosity               = viscosity * fluid_density;

                if (fluid_density > 0.0000000000001){
                    double norm_of_slip_vel            = MathUtils<double>::Norm3(fluid_vel - particle_vel);
                    double area                        = M_PI * mRadius * mRadius;
                    double shahs_term_vel              = 0.0;
                    double beta                        = 0.0;
                    double F0                          = 0.0;
                    double regularization_v            = 0.2 * mRadius;
                    const array_1d<double,3>& peso     = volume * particle_density * gravity;
                    double drag_coeff;
                    double reynolds;

                    if (!non_newtonian_OPTION){ //Newtonian
                        ComputeReynoldsNumber(non_newtonian_OPTION, norm_of_slip_vel, fluid_density, dynamic_viscosity, reynolds);
                        CalculateDragCoefficient(non_newtonian_OPTION, reynolds, sphericity, drag_coeff, drag_modifier_type);
                        drag_force = 0.5 * fluid_density * area * drag_coeff * norm_of_slip_vel * (fluid_vel - particle_vel);
                    }

                    else {
                        shahs_term_vel = CalculateShahsTerm(power_law_n, power_law_K, power_law_tol, fluid_density, particle_density, sphericity, drag_modifier_type);

                        if (!manually_imposed_drag_law_OPTION){
                            F0 = gel_strength * 4.0 * area; //initial value
                            beta = (MathUtils<double>::Norm3(peso - buoyancy) - F0) / shahs_term_vel; //slope
                        }

                        else {
                            F0 = initial_drag_force; //initial value
                            beta = drag_law_slope; //slope
                        }

                        if (norm_of_slip_vel) {

                            if (norm_of_slip_vel >= regularization_v){
                                drag_force = (F0 + beta * norm_of_slip_vel) / norm_of_slip_vel * (fluid_vel - particle_vel);
                            }

                            else {
                                drag_force = (F0 + beta * regularization_v) / regularization_v * (fluid_vel - particle_vel);
                            }

                        }

                    }

                }

                buoyancy = - pressure_grad * volume;

            }

            else {
                noalias(drag_force)    = ZeroVector(3);
                noalias(pressure_grad) = ZeroVector(3);
                noalias(buoyancy)      = ZeroVector(3);
            }

          KRATOS_CATCH("")
      }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

}  // namespace Kratos.
