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
#include "cuboidcluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    CuboidCluster3D::CuboidCluster3D() : Cluster3D() {}
            
      
    CuboidCluster3D::CuboidCluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    CuboidCluster3D::CuboidCluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    CuboidCluster3D::CuboidCluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer CuboidCluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new CuboidCluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    CuboidCluster3D::~CuboidCluster3D() {}
      
    
    void CuboidCluster3D::CustomInitialize(ProcessInfo& r_process_info) {
        
        int number_of_spheres_in_first_sector = 9;
        int number_of_passes = 9;
        int number_of_spheres_per_pass = 4 * number_of_spheres_in_first_sector - 4;
        mListOfRadii.resize(number_of_spheres_per_pass * number_of_passes);
        mListOfCoordinates.resize(number_of_spheres_per_pass * number_of_passes);
        mListOfSphericParticles.resize(number_of_spheres_per_pass * number_of_passes);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        //double sinusoidal_factor_1 = Globals::Pi *
        //PROPERLY CALCULATE THE DIVISIONS
        
        double u, v;
        double a = 0.5, b = 0.5, c = 0.5;
        
        for (int j = 0; j < number_of_passes; j++) {
        
            for (int i = 0; i < number_of_spheres_per_pass; i++) {
                
//                u = j * 2.0 * Globals::Pi / (number_of_passes - 1);
//                v = Globals::Pi * i * 0.25 / (number_of_spheres_in_first_sector - 1);
                u = -Globals::Pi + Globals::Pi * i * 0.5 / (number_of_spheres_in_first_sector - 1);
                v = -0.5 * Globals::Pi + j * Globals::Pi / (number_of_passes - 1);
            
                mListOfRadii[i + number_of_spheres_per_pass * j]= 0.1 * cl;
                
                mListOfCoordinates[i + number_of_spheres_per_pass * j][0] = cl * a * cos(v) * cos(u) / 
                    pow(pow(sin(v), 6) * (pow(sin(u), 6) + pow(cos(u), 6)) + pow(cos(v), 6), 0.166666667); 
                mListOfCoordinates[i + number_of_spheres_per_pass * j][1] = cl * b * cos(v) * sin(u) / 
                    pow(pow(sin(v), 6) * (pow(sin(u), 6) + pow(cos(u), 6)) + pow(cos(v), 6), 0.166666667); 
                mListOfCoordinates[i + number_of_spheres_per_pass * j][2] = cl * c * sin(v) / 
                    pow(pow(sin(v), 6) * (pow(sin(u), 6) + pow(cos(u), 6)) + pow(cos(v), 6), 0.166666667);
//                mListOfCoordinates[i + number_of_spheres_per_pass * j][0] =
//                    
//                    cl * 0.5 * sin(u) * cos(v) / 
//                    pow(pow(sin(u), 6) * (pow(sin(v), 6) + pow(cos(v), 6)) + pow(cos(u), 6), 0.166666667); 
//                
//                mListOfCoordinates[i + number_of_spheres_per_pass * j][1] = 
//                        
//                    cl * 0.375 * sin(u) * sin(v) / 
//                    pow(pow(sin(u), 6) * (pow(sin(v), 6) + pow(cos(v), 6)) + pow(cos(u), 6), 0.166666667); 
//                
//                mListOfCoordinates[i + number_of_spheres_per_pass * j][2] = 
//                        
//                    cl * 0.25 * cos(u) / 
//                    pow(pow(sin(u), 6) * (pow(sin(v), 6) + pow(cos(v), 6)) + pow(cos(u), 6), 0.166666667);

                //double y = pow(8, 0.33333333333333);
                //u = -0.5 * Globals::Pi + j * Globals::Pi / (number_of_passes - 1)
                        
                //v = -Globals::Pi + Globals::Pi * i * 0.5 / (number_of_spheres_in_first_sector - 1)
            
            }
        }
        
        double particle_density = this->SlowGetDensity(); /////////////////////////////USE FAST
         
        double cluster_volume = 0.3333333333 * 4.0 * Globals::Pi * 0.5 * 0.375 * 0.25 * cl * cl * cl; ////APPROXIMATE VOLUME, CALCULATE MORE EXACTLY
        
        double cluster_mass = particle_density * cluster_volume;
        
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = cluster_mass;
        
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = 0.2 * cluster_mass * cl * cl * (0.375 * 0.375 + 0.25 * 0.25);
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = 0.2 * cluster_mass * cl * cl * (0.5 * 0.5 + 0.25 * 0.25);
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = 0.2 * cluster_mass * cl * cl * (0.5 * 0.5 + 0.375 * 0.375);
         
        array_1d<double, 3> base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = cluster_mass;
        GetGeometry()[0].FastGetSolutionStepValue(CLUSTER_VOLUME) = cluster_volume;
        GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MATERIAL) = this->SlowGetParticleMaterial();
        
        Quaternion<double>& Orientation = GetGeometry()[0].FastGetSolutionStepValue(ORIENTATION);
        Orientation.normalize();

        array_1d<double, 3> angular_velocity = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        
        array_1d<double, 3> angular_momentum;
        double LocalTensor[3][3];
        double GlobalTensor[3][3];
        GeometryFunctions::ConstructLocalTensor(base_principal_moments_of_inertia, LocalTensor);
        GeometryFunctions::QuaternionTensorLocal2Global(Orientation, LocalTensor, GlobalTensor);                   
        GeometryFunctions::ProductMatrix3X3Vector3X1(GlobalTensor, angular_velocity, angular_momentum);
        noalias(this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_MOMENTUM)) = angular_momentum;
        
        array_1d<double, 3> local_angular_velocity;
        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, angular_velocity, local_angular_velocity);
        noalias(this->GetGeometry()[0].FastGetSolutionStepValue(LOCAL_ANGULAR_VELOCITY)) = local_angular_velocity;  
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CuboidCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CuboidCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CuboidCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CuboidCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CuboidCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CuboidCluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CuboidCluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void CuboidCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CuboidCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CuboidCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void CuboidCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double CuboidCluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}
    int CuboidCluster3D::SlowGetParticleMaterial()                                  { return GetProperties()[PARTICLE_MATERIAL];}

}  // namespace Kratos.

