//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "capsulecluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    CapsuleCluster3D::CapsuleCluster3D() : Cluster3D() {}
      
    CapsuleCluster3D::CapsuleCluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
    CapsuleCluster3D::CapsuleCluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

    CapsuleCluster3D::CapsuleCluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}
      
    Element::Pointer CapsuleCluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new CapsuleCluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }      

    // Destructor
    CapsuleCluster3D::~CapsuleCluster3D() {}

    void CapsuleCluster3D::CustomInitialize() {
        
        int number_of_spheres = 5;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        mListOfRadii[0]= 0.35 * cl;
        mListOfRadii[1]= 0.35 * cl;
        mListOfRadii[2]= 0.35 * cl;   
        mListOfRadii[3]= 0.35 * cl;
        mListOfRadii[4]= 0.35 * cl;
        
        mListOfCoordinates[0][0] = 0.0000 * cl; mListOfCoordinates[0][1] = 0.0; mListOfCoordinates[0][2] = 0.0;
        mListOfCoordinates[1][0] = 0.2325 * cl; mListOfCoordinates[1][1] = 0.0; mListOfCoordinates[1][2] = 0.0;
        mListOfCoordinates[2][0] = 0.4650 * cl; mListOfCoordinates[2][1] = 0.0; mListOfCoordinates[2][2] = 0.0;
        mListOfCoordinates[3][0] = 0.6975 * cl; mListOfCoordinates[3][1] = 0.0; mListOfCoordinates[3][2] = 0.0; 
        mListOfCoordinates[4][0] = 0.9300 * cl; mListOfCoordinates[4][1] = 0.0; mListOfCoordinates[4][2] = 0.0;  
        
        double particle_density = this->SlowGetDensity();
        double cluster_volume = 4.0 * KRATOS_M_PI_3 * 0.35 * 0.35 * 0.35 * cl * cl * cl + KRATOS_M_PI * 0.35 * 0.35 * 0.93 * cl * cl * cl;
        double cluster_mass = particle_density * cluster_volume;
        
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = cluster_mass;
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = 0.5 * cluster_mass * 0.35 * 0.35 * cl * cl;
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = 0.08333333333 * cluster_mass * (3.0 * 0.35 * 0.35 + 4.0 * 1.0 * 1.0) * cl * cl;
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = 0.08333333333 * cluster_mass * (3.0 * 0.35 * 0.35 + 4.0 * 1.0 * 1.0) * cl * cl;
         
        array_1d<double, 3> base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);     
    }     
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CapsuleCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CapsuleCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CapsuleCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CapsuleCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CapsuleCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CapsuleCluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CapsuleCluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void CapsuleCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY
        KRATOS_CATCH("")
    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CapsuleCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CapsuleCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void CapsuleCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double CapsuleCluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

