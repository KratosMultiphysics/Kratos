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
#include "pillcluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    PillCluster3D::PillCluster3D() : Cluster3D() {}
            
      
    PillCluster3D::PillCluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    PillCluster3D::PillCluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    PillCluster3D::PillCluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer PillCluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new PillCluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    PillCluster3D::~PillCluster3D() {}
      
    
    void PillCluster3D::CustomInitialize() {

//        int number_of_spheres_in_first_sector = 6;
//        int number_of_passes = 5;
//        int number_of_spheres_per_pass = 4 * number_of_spheres_in_first_sector - 4;
//        mListOfRadii.resize(number_of_spheres_per_pass * number_of_passes);
//        mListOfCoordinates.resize(number_of_spheres_per_pass * number_of_passes);
//        mListOfSphericParticles.resize(number_of_spheres_per_pass * number_of_passes);
        
        int total_number_of_particles = 20 + 16 + 12 + 6 + 1;
        
        mListOfRadii.resize(total_number_of_particles);
        mListOfCoordinates.resize(total_number_of_particles);
        mListOfSphericParticles.resize(total_number_of_particles);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        double sinusoidal_factor = 0.5 * KRATOS_M_PI / (20.0/4.0 - 1);
             
        for (int i = 0; i < 20; i++) {           
            mListOfCoordinates[i][0] = 0.45 * cl * cos(sinusoidal_factor * i); ////APPROXIMATE COORDINATES, CALCULATE MORE EXACTLY
            mListOfCoordinates[i][1] = 0.45 * cl * sin(sinusoidal_factor * i);
            mListOfCoordinates[i][2] = 0.0;
        }
        
        sinusoidal_factor = 0.5 * KRATOS_M_PI / (16.0/4.0 - 1);
        
        for (int i = 0; i < 16; i++) {           
            mListOfCoordinates[20 + i][0] = 0.35 * cl * cos(sinusoidal_factor * i);
            mListOfCoordinates[20 + i][1] = 0.35 * cl * sin(sinusoidal_factor * i);
            mListOfCoordinates[20 + i][2] = 0.0;
        }
        
        sinusoidal_factor = 0.5 * KRATOS_M_PI / (12.0/4.0 - 1);
        
        for (int i = 0; i < 12; i++) {           
            mListOfCoordinates[36 + i][0] = 0.25 * cl * cos(sinusoidal_factor * i); 
            mListOfCoordinates[36 + i][1] = 0.25 * cl * sin(sinusoidal_factor * i);
            mListOfCoordinates[36 + i][2] = 0.0;
        }
        
        sinusoidal_factor = 0.333333333333333333 * KRATOS_M_PI;
        
        for (int i = 0; i <  6; i++) {           
            mListOfCoordinates[48 + i][0] = 0.15 * cl * cos(sinusoidal_factor * i); 
            mListOfCoordinates[48 + i][1] = 0.15 * cl * sin(sinusoidal_factor * i);
            mListOfCoordinates[48 + i][2] = 0.0;
        }
        
        for (int i = 0; i <  1; i++) {           
            mListOfCoordinates[54 + i][0] = 0.0; 
            mListOfCoordinates[54 + i][1] = 0.0;
            mListOfCoordinates[54 + i][2] = 0.0;
        }
        
        for (int i = 0; i < 55; i++) {           
            mListOfRadii[i]= 0.1 * cl;
        }        
                
//        double sinusoidal_factor = 0.5 * KRATOS_M_PI / (number_of_spheres_in_first_sector - 1);
//        
//        for (int j = 0; j < number_of_passes; j++) {
//        
//            for (int i = 0; i < number_of_spheres_per_pass; i++) {
//            
//                mListOfRadii[i + number_of_spheres_per_pass * j]= 0.1 * cl;
//            
//                mListOfCoordinates[i + number_of_spheres_per_pass * j][0] = (0.05 + 0.1 * j) * cl * cos(sinusoidal_factor * i); ////APPROXIMATE COORDINATES, CALCULATE MORE EXACTLY
//                mListOfCoordinates[i + number_of_spheres_per_pass * j][1] = (0.05 + 0.1 * j) * cl * sin(sinusoidal_factor * i);
//                mListOfCoordinates[i + number_of_spheres_per_pass * j][2] = 0.0;
//            
//            }
//        }
        
        //double particle_density = this->SlowGetDensity(); /////////////////////////////USE FAST
         
        //double cluster_volume = KRATOS_M_PI * 0.5 * 0.5 * cl * cl * 2.0 * 0.1 * cl; ////APPROXIMATE VOLUME, CALCULATE MORE EXACTLY
        
        //double cluster_mass = particle_density * cluster_volume;
        
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * KRATOS_M_PI * 0.5 * 0.5 * cl * cl * 2.0 * 0.1 * cl;
        
        array_1d<double,3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = 0.25 * cluster_mass * 0.5 * 0.5 * cl * cl;
        base_principal_moments_of_inertia[1] = 0.25 * cluster_mass * 0.5 * 0.5 * cl * cl;
        base_principal_moments_of_inertia[2] = 0.50 * cluster_mass * 0.5 * 0.5 * cl * cl;
         
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void PillCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void PillCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void PillCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void PillCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void PillCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void PillCluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void PillCluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void PillCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void PillCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void PillCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void PillCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double PillCluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

