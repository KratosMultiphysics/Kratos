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
#include "ellipsoidcluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    EllipsoidCluster3D::EllipsoidCluster3D() : Cluster3D() {}
            
      
    EllipsoidCluster3D::EllipsoidCluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    EllipsoidCluster3D::EllipsoidCluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    EllipsoidCluster3D::EllipsoidCluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer EllipsoidCluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new EllipsoidCluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    EllipsoidCluster3D::~EllipsoidCluster3D() {}
      
    
    void EllipsoidCluster3D::CustomInitialize() {
        
        int number_of_spheres_in_first_sector = 7;
        int number_of_passes = 9;
        int number_of_spheres_per_pass = 4 * number_of_spheres_in_first_sector - 4;
        mListOfRadii.resize(number_of_spheres_per_pass * number_of_passes);
        mListOfCoordinates.resize(number_of_spheres_per_pass * number_of_passes);
        mListOfSphericParticles.resize(number_of_spheres_per_pass * number_of_passes);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
                
        //double sinusoidal_factor_1 = KRATOS_M_PI *
        //PROPERLY CALCULATE THE DIVISIONS
        
        for (int j = 0; j < number_of_passes; j++) {
        
            for (int i = 0; i < number_of_spheres_per_pass; i++) {
            
                mListOfRadii[i + number_of_spheres_per_pass * j]= 0.1 * cl;
            
                mListOfCoordinates[i + number_of_spheres_per_pass * j][0] = 
                        
                    cl * 0.5 * cos(-0.5 * KRATOS_M_PI + j * KRATOS_M_PI / (number_of_passes - 1))
                        * cos(-KRATOS_M_PI + KRATOS_M_PI * i * 0.5 / (number_of_spheres_in_first_sector - 1)); 
                
                mListOfCoordinates[i + number_of_spheres_per_pass * j][1] = 
                        
                    cl * 0.375 * cos(-0.5 * KRATOS_M_PI + j * KRATOS_M_PI / (number_of_passes - 1))
                        * sin(-KRATOS_M_PI + KRATOS_M_PI * i * 0.5 / (number_of_spheres_in_first_sector - 1)); 
                
                mListOfCoordinates[i + number_of_spheres_per_pass * j][2] = 
                        
                    cl * 0.25 * sin(-0.5 * KRATOS_M_PI + j * KRATOS_M_PI / (number_of_passes - 1));
                
                //u = -0.5 * KRATOS_M_PI + j * KRATOS_M_PI / (number_of_passes - 1)
                        
                //v = -KRATOS_M_PI + KRATOS_M_PI * i * 0.5 / (number_of_spheres_in_first_sector - 1)
            
            }
        }
        
        //double particle_density = this->SlowGetDensity(); /////////////////////////////USE FAST
         
        //double cluster_volume = 0.3333333333 * 4.0 * KRATOS_M_PI * 0.5 * 0.375 * 0.25 * cl * cl * cl; ////APPROXIMATE VOLUME, CALCULATE MORE EXACTLY
        
        //double cluster_mass = particle_density * cluster_volume;
        
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.3333333333 * 4.0 * KRATOS_M_PI * 0.5 * 0.375 * 0.25 * cl * cl * cl;
        
        array_1d<double,3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA); 
        
        base_principal_moments_of_inertia[0] = 0.2 * cluster_mass * cl * cl * (0.375 * 0.375 + 0.25 * 0.25);
        base_principal_moments_of_inertia[1] = 0.2 * cluster_mass * cl * cl * (0.5 * 0.5 + 0.25 * 0.25);
        base_principal_moments_of_inertia[2] = 0.2 * cluster_mass * cl * cl * (0.5 * 0.5 + 0.375 * 0.375);
         
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void EllipsoidCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void EllipsoidCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void EllipsoidCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void EllipsoidCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void EllipsoidCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void EllipsoidCluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void EllipsoidCluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void EllipsoidCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void EllipsoidCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void EllipsoidCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void EllipsoidCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double EllipsoidCluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

