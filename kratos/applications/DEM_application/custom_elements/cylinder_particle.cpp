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
      : SphericParticle(){/*mInitializedVariablesFlag = 0;*/}

      CylinderParticle::CylinderParticle(IndexType NewId, GeometryType::Pointer pGeometry)
      : SphericParticle(NewId, pGeometry){/*mInitializedVariablesFlag = 0;*/}

      CylinderParticle::CylinderParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : SphericParticle(NewId, pGeometry, pProperties){/*mInitializedVariablesFlag = 0;*/}

      CylinderParticle::CylinderParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : SphericParticle(NewId, ThisNodes){/*mInitializedVariablesFlag = 0;*/}

      Element::Pointer CylinderParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
           return SphericParticle::Pointer(new CylinderParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));

      }

      /// Destructor.
      CylinderParticle::~CylinderParticle(){}

      void CylinderParticle::Initialize()
      {
          KRATOS_TRY

          mRadius                   = GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
          double density            = GetDensity();          
          double& mass              = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS);
          mass                      = KRATOS_M_PI_3 * density * GetRadius() * GetRadius() * 1.0;
          mRealMass                 = mass;

          //if (mRotationOption){
          if (this->Is(DEMFlags::HAS_ROTATION) ){
              double& moment_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
              moment_of_inertia = 0.5 * mass * GetRadius() * GetRadius();   
              GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = moment_of_inertia;
          }                                                                        

          CustomInitialize();

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderParticle::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double dt, const array_1d<double,3>& gravity, int search_control)
      {
          KRATOS_TRY

          array_1d<double, 3> additionally_applied_force;
          array_1d<double, 3> additionally_applied_moment;
          array_1d<double, 3> initial_rotation_moment;     
          array_1d<double, 3>& elastic_force       = this->GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);
          array_1d<double, 3>& contact_force       = this->GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCES);
          array_1d<double, 3>& rigid_element_force = this->GetGeometry()[0].FastGetSolutionStepValue(RIGID_ELEMENT_FORCE);

          mContactForce.clear();
          mContactMoment.clear();
          additionally_applied_force.clear();
          additionally_applied_moment.clear();
          initial_rotation_moment.clear();
          elastic_force.clear();
          contact_force.clear();
          rigid_element_force.clear();

          bool multi_stage_RHS = false;

          ComputeBallToBallContactForce(/*contact_force, contact_moment, */elastic_force, contact_force, initial_rotation_moment, rCurrentProcessInfo, dt, multi_stage_RHS );

          //Cfeng,RigidFace
          
          if (mFemOldNeighbourIds.size() > 0)
          {
            ComputeBallToRigidFaceContactForce(/*contact_force, contact_moment,*/ elastic_force, contact_force, initial_rotation_moment, rigid_element_force, rCurrentProcessInfo, dt, search_control);
          }
              
          ComputeAdditionalForces(additionally_applied_force, additionally_applied_moment, rCurrentProcessInfo, gravity);

          
          rRightHandSideVector[0] = mContactForce[0]  + additionally_applied_force[0];
          rRightHandSideVector[1] = mContactForce[1]  + additionally_applied_force[1];
          rRightHandSideVector[2] = mContactForce[2]  + additionally_applied_force[2];
          rRightHandSideVector[3] = mContactMoment[0] + additionally_applied_moment[0];
          rRightHandSideVector[4] = mContactMoment[1] + additionally_applied_moment[1];
          rRightHandSideVector[5] = mContactMoment[2] + additionally_applied_moment[2];
          
          array_1d<double,3>& total_forces  = this->GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
          array_1d<double,3>& total_moment = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);

          for (int i = 0; i < 3; i++){
              total_forces[i] = rRightHandSideVector[i];
              total_moment[i] = rRightHandSideVector[3 + i];
          }	  

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

