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
         return DiscreteElement::Pointer(new SphericSwimmingParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
      }

      /// Destructor.
      SphericSwimmingParticle::~SphericSwimmingParticle(){}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************


      void SphericSwimmingParticle::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          const array_1d<double,3>& gravity         = rCurrentProcessInfo[GRAVITY];
          array_1d<double,3>& drag_force            = GetGeometry()(0)->FastGetSolutionStepValue(DRAG_FORCE);
          array_1d<double,3>& buoyancy              = GetGeometry()(0)->FastGetSolutionStepValue(BUOYANCY);
          array_1d<double, 3> contact_force;
          array_1d<double, 3> contact_moment;

          contact_force[0]  = 0.0;
          contact_force[1]  = 0.0;
          contact_force[2]  = 0.0;
          contact_moment[0] = 0.0;
          contact_moment[1] = 0.0;
          contact_moment[2] = 0.0;

          ComputeNewNeighboursHistoricalData();
          ComputeBallToBallContactForce(contact_force, contact_moment, rCurrentProcessInfo);
          ComputeBallToSurfaceContactForce(contact_force, contact_moment, rCurrentProcessInfo);
          ComputeDragForces(rCurrentProcessInfo);

          rRightHandSideVector[0] = contact_force[0] + buoyancy[0] + drag_force[0] + mRealMass * gravity[0];
          rRightHandSideVector[1] = contact_force[1] + buoyancy[1] + drag_force[1] + mRealMass * gravity[1];
          rRightHandSideVector[2] = contact_force[2] + buoyancy[2] + drag_force[2] + mRealMass * gravity[2];
          rRightHandSideVector[3] = contact_moment[0];
          rRightHandSideVector[4] = contact_moment[1];
          rRightHandSideVector[5] = contact_moment[2];

          KRATOS_CATCH( "" )
      }


    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      void SphericSwimmingParticle::ComputeDragForces(ProcessInfo& rCurrentProcessInfo)
      {
            KRATOS_TRY

            array_1d<double,3>& pressure_grad       = GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE_GRAD_PROJECTED);
            array_1d<double,3>& drag_force          = GetGeometry()(0)->FastGetSolutionStepValue(DRAG_FORCE);
            array_1d<double,3>& buoyancy            = GetGeometry()(0)->FastGetSolutionStepValue(BUOYANCY);

            if (GetGeometry()[0].IsFixed(VELOCITY_X) == false){
                const double& fluid_density = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED);
                double volume = (1.333333333333333 * M_PI) * (mRadius * mRadius * mRadius);
                const array_1d<double,3>& fluid_vel = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
                const array_1d<double,3>& particle_vel = GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);

                if (fluid_density > 0.0000000000001){
                    drag_force = (0.235 * M_PI) * (fluid_density * mRadius * mRadius) * MathUtils<double>::Norm3(fluid_vel - particle_vel) * (fluid_vel - particle_vel);
                }

                else {
                    noalias(drag_force)    = ZeroVector(3);
                    noalias(pressure_grad) = ZeroVector(3);
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
