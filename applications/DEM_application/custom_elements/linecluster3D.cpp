//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva $
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
#include "cluster3D.h"
#include "linecluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"


namespace Kratos {
        
    // using namespace GeometryFunctions;

    LineCluster3D::LineCluster3D() : Cluster3D() {}
            
      
    LineCluster3D::LineCluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    LineCluster3D::LineCluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    LineCluster3D::LineCluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer LineCluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new LineCluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    LineCluster3D::~LineCluster3D() {}

      
    void LineCluster3D::Initialize() {
        
        int number_of_spheres=6;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        mListOfRadii[0]=0.1;
        mListOfRadii[1]=0.1;
        mListOfRadii[2]=0.1;   
        mListOfRadii[3]=0.1;
        mListOfRadii[4]=0.1;
        mListOfRadii[5]=0.1;
        
        mListOfCoordinates[0][0] =-0.35; mListOfCoordinates[0][1] = 0.0; mListOfCoordinates[0][2] = 0.0;
        mListOfCoordinates[1][0] =-0.25; mListOfCoordinates[1][1] = 0.0; mListOfCoordinates[1][2] = 0.0;
        mListOfCoordinates[2][0] =-0.10; mListOfCoordinates[2][1] = 0.0; mListOfCoordinates[2][2] = 0.0;
        mListOfCoordinates[3][0] = 0.10; mListOfCoordinates[3][1] = 0.0; mListOfCoordinates[3][2] = 0.0; 
        mListOfCoordinates[4][0] = 0.25; mListOfCoordinates[4][1] = 0.0; mListOfCoordinates[4][2] = 0.0;  
        mListOfCoordinates[5][0] = 0.35; mListOfCoordinates[5][1] = 0.0; mListOfCoordinates[5][2] = 0.0;  
                        
        mVelocity                  = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        mTotalForces               = GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
        double sqrt_of_mass        = GetGeometry()[0].FastGetSolutionStepValue(SQRT_OF_MASS);
        mSqrtOfRealMass            = sqrt_of_mass;
        mAngularVelocity           = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY); 
        mParticleMoment            = GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);
        mPrincipalMomentsOfInertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        mEulerAngles               = GetGeometry()[0].FastGetSolutionStepValue(EULER_ANGLES);                    
    
    } 
        
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {
        
        KRATOS_TRY

        ElementalDofList.resize(0);

//      for (unsigned int i = 0; i < GetGeometry().size(); i++){
//          ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
//          ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
//          if (GetGeometry().WorkingSpaceDimension() == 3) {
//              ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
//          }
// 
//          ElementalDofList.push_back(GetGeometry()[i].pGetDof(ANGULAR_VELOCITY_X));
//          ElementalDofList.push_back(GetGeometry()[i].pGetDof(ANGULAR_VELOCITY_Y));
//          if (GetGeometry().WorkingSpaceDimension() == 3) {
//              ElementalDofList.push_back(GetGeometry()[i].pGetDof(ANGULAR_VELOCITY_Z));
//          }
// 
//      }

        KRATOS_CATCH("")
        
    }


    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void LineCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void LineCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    

}  // namespace Kratos.

