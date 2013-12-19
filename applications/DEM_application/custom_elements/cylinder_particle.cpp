//
//   Project Name:        Kratos
//   Last Modified by:    $Author: M.Santasusana $
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
#include "cylinder_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application.h"


namespace Kratos
{
     // using namespace GeometryFunctions;

      CylinderParticle::CylinderParticle()
      : SphericParticle(){mInitializedVariablesFlag = 0;}

      CylinderParticle::CylinderParticle(IndexType NewId, GeometryType::Pointer pGeometry)
      : SphericParticle(NewId, pGeometry){mInitializedVariablesFlag = 0;}

      CylinderParticle::CylinderParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : SphericParticle(NewId, pGeometry, pProperties){mInitializedVariablesFlag = 0;}

      CylinderParticle::CylinderParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : SphericParticle(NewId, ThisNodes){mInitializedVariablesFlag = 0;}

      Element::Pointer CylinderParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
           return SphericParticle::Pointer(new CylinderParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));

      }

      /// Destructor.
      CylinderParticle::~CylinderParticle(){}

      void CylinderParticle::Initialize()
      {
          KRATOS_TRY

          mDimension                = 2;
          mRadius                   = GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
          mYoung                    = GetGeometry()(0)->FastGetSolutionStepValue(YOUNG_MODULUS);         
          mPoisson                  = GetGeometry()(0)->FastGetSolutionStepValue(POISSON_RATIO);
          mTgOfFrictionAngle        = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_FRICTION);
          mLnOfRestitCoeff          = GetGeometry()(0)->FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);
          double& density           = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_DENSITY);
          //double& mass              = GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS);
          double& sqrt_of_mass      = GetGeometry()(0)->FastGetSolutionStepValue(SQRT_OF_MASS);
          double& moment_of_inertia = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);

          double mass               = M_PI * density * mRadius * mRadius * 1.0;
          sqrt_of_mass              = sqrt(mass);
          moment_of_inertia         = 0.5 * mass * mRadius * mRadius;
          mRealMass                 = mass;          
          mSqrtOfRealMass           = sqrt_of_mass;
          mMomentOfInertia          = moment_of_inertia;

          CustomInitialize();

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          array_1d<double, 3> contact_force;
          array_1d<double, 3> contact_moment;
          array_1d<double, 3> additionally_applied_force;
          array_1d<double, 3> additionally_applied_moment;
          array_1d<double, 3> initial_rotation_moment;     
          array_1d<double, 3>& elastic_force = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);

          contact_force.clear();
          contact_moment.clear();
          additionally_applied_force.clear();
          additionally_applied_moment.clear();
          initial_rotation_moment.clear();
          elastic_force.clear();

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
          }
          
          ComputeAdditionalForces(contact_force, contact_moment, additionally_applied_force, additionally_applied_moment, rCurrentProcessInfo);

          rRightHandSideVector[0] = contact_force[0]  + additionally_applied_force[0];
          rRightHandSideVector[1] = contact_force[1]  + additionally_applied_force[1];
          rRightHandSideVector[2] = contact_force[2]  + additionally_applied_force[2];
          //rRightHandSideVector[2] = contact_force[2]  + additionally_applied_force[2];
          rRightHandSideVector[3] = contact_moment[0] + additionally_applied_moment[0];
          rRightHandSideVector[4] = contact_moment[1] + additionally_applied_moment[0];
          rRightHandSideVector[5] = contact_moment[2] + additionally_applied_moment[0];
          //rRightHandSideVector[5] = contact_moment[2] + additionally_applied_moment[0];

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
      void CylinderParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo){}
      void CylinderParticle::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo){}
      void CylinderParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
      void CylinderParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

}  // namespace Kratos.

