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
#include "ringcluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    RingCluster3D::RingCluster3D() : Cluster3D() {}
            
      
    RingCluster3D::RingCluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    RingCluster3D::RingCluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    RingCluster3D::RingCluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer RingCluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new RingCluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    RingCluster3D::~RingCluster3D() {}
      
    
    void RingCluster3D::CustomInitialize() {
        
        //total number of spheres in the circumference is 37 + 2 * 36 + 35 = 144
        //each of the four (1/4)circles divided into 36 sectors
        int number_of_spheres_in_first_sector = 6;
        int number_of_passes = 1;
        int number_of_spheres_per_pass = 4 * number_of_spheres_in_first_sector - 4;
        mListOfRadii.resize(number_of_spheres_per_pass * number_of_passes);
        mListOfCoordinates.resize(number_of_spheres_per_pass * number_of_passes);
        mListOfSphericParticles.resize(number_of_spheres_per_pass * number_of_passes);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        double sinusoidal_factor = 0.5 * KRATOS_M_PI / (number_of_spheres_in_first_sector - 1);
        
        for (int j = 0; j < number_of_passes; j++) {
        
            for (int i = 0; i < number_of_spheres_per_pass; i++) {
            
                mListOfRadii[i + number_of_spheres_per_pass * j]= 0.1 * cl;
            
                mListOfCoordinates[i + number_of_spheres_per_pass * j][0] = (0.5 + 0.1 * j) * cl * cos(sinusoidal_factor * i); ////APPROXIMATE COORDINATES, CALCULATE MORE EXACTLY
                mListOfCoordinates[i + number_of_spheres_per_pass * j][1] = (0.5 + 0.1 * j) * cl * sin(sinusoidal_factor * i);
                mListOfCoordinates[i + number_of_spheres_per_pass * j][2] = 0.0;
            
            }
        }
        
        //double particle_density = this->SlowGetDensity(); /////////////////////////////USE FAST
         
        //double cluster_volume = KRATOS_M_PI * 0.1 * 0.1 * cl * cl * 2.0 * KRATOS_M_PI * 0.5 * cl; ////APPROXIMATE VOLUME, CALCULATE MORE EXACTLY
        
        //double cluster_mass = particle_density * cluster_volume;
        
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * KRATOS_M_PI * 0.1 * 0.1 * cl * cl * 2.0 * KRATOS_M_PI * 0.5 * cl;
        
        array_1d<double,3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = 0.5 * cluster_mass * 0.5 * 0.5 * cl * cl;
        base_principal_moments_of_inertia[1] = 0.5 * cluster_mass * 0.5 * 0.5 * cl * cl;
        base_principal_moments_of_inertia[2] =       cluster_mass * 0.5 * 0.5 * cl * cl;
         
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void RingCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void RingCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void RingCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void RingCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void RingCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void RingCluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void RingCluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void RingCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void RingCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void RingCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void RingCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double RingCluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

