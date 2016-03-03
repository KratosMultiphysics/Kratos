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
#include "cubecluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    CubeCluster3D::CubeCluster3D() : Cluster3D() {}
            
      
    CubeCluster3D::CubeCluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    CubeCluster3D::CubeCluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    CubeCluster3D::CubeCluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer CubeCluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new CubeCluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    CubeCluster3D::~CubeCluster3D() {}
      
    
    void CubeCluster3D::CustomInitialize() {
        
        int number_of_spheres=26;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        for (int i = 0; i < 26; i++){ mListOfRadii[i]= 0.4 * cl; }
        
        mListOfCoordinates[ 0][0] =  0.5 * cl; mListOfCoordinates[ 0][1] =  0.5 * cl; mListOfCoordinates[ 0][2] =  0.5 * cl;
        mListOfCoordinates[ 1][0] =  0.5 * cl; mListOfCoordinates[ 1][1] =  0.5 * cl; mListOfCoordinates[ 1][2] =  0.0 * cl;
        mListOfCoordinates[ 2][0] =  0.5 * cl; mListOfCoordinates[ 2][1] =  0.5 * cl; mListOfCoordinates[ 2][2] = -0.5 * cl;
        mListOfCoordinates[ 3][0] =  0.5 * cl; mListOfCoordinates[ 3][1] =  0.0 * cl; mListOfCoordinates[ 3][2] =  0.5 * cl; 
        mListOfCoordinates[ 4][0] =  0.5 * cl; mListOfCoordinates[ 4][1] =  0.0 * cl; mListOfCoordinates[ 4][2] = -0.5 * cl;  
        mListOfCoordinates[ 5][0] =  0.5 * cl; mListOfCoordinates[ 5][1] = -0.5 * cl; mListOfCoordinates[ 5][2] =  0.5 * cl;
        mListOfCoordinates[ 6][0] =  0.5 * cl; mListOfCoordinates[ 6][1] = -0.5 * cl; mListOfCoordinates[ 6][2] =  0.0 * cl;  
        mListOfCoordinates[ 7][0] =  0.5 * cl; mListOfCoordinates[ 7][1] = -0.5 * cl; mListOfCoordinates[ 7][2] = -0.5 * cl;
        mListOfCoordinates[ 8][0] =  0.0 * cl; mListOfCoordinates[ 8][1] =  0.5 * cl; mListOfCoordinates[ 8][2] =  0.5 * cl;  
        mListOfCoordinates[ 9][0] =  0.0 * cl; mListOfCoordinates[ 9][1] =  0.5 * cl; mListOfCoordinates[ 9][2] =  0.0 * cl;
        mListOfCoordinates[10][0] =  0.0 * cl; mListOfCoordinates[10][1] =  0.5 * cl; mListOfCoordinates[10][2] = -0.5 * cl;
        mListOfCoordinates[11][0] =  0.0 * cl; mListOfCoordinates[11][1] =  0.0 * cl; mListOfCoordinates[11][2] =  0.5 * cl;
        mListOfCoordinates[12][0] =  0.0 * cl; mListOfCoordinates[12][1] =  0.0 * cl; mListOfCoordinates[12][2] = -0.5 * cl;
        mListOfCoordinates[13][0] =  0.0 * cl; mListOfCoordinates[13][1] = -0.5 * cl; mListOfCoordinates[13][2] =  0.5 * cl; 
        mListOfCoordinates[14][0] =  0.0 * cl; mListOfCoordinates[14][1] = -0.5 * cl; mListOfCoordinates[14][2] =  0.0 * cl;  
        mListOfCoordinates[15][0] =  0.0 * cl; mListOfCoordinates[15][1] = -0.5 * cl; mListOfCoordinates[15][2] = -0.5 * cl;
        mListOfCoordinates[16][0] = -0.5 * cl; mListOfCoordinates[16][1] =  0.5 * cl; mListOfCoordinates[16][2] =  0.5 * cl;  
        mListOfCoordinates[17][0] = -0.5 * cl; mListOfCoordinates[17][1] =  0.5 * cl; mListOfCoordinates[17][2] =  0.0 * cl;
        mListOfCoordinates[18][0] = -0.5 * cl; mListOfCoordinates[18][1] =  0.5 * cl; mListOfCoordinates[18][2] = -0.5 * cl;  
        mListOfCoordinates[19][0] = -0.5 * cl; mListOfCoordinates[19][1] =  0.0 * cl; mListOfCoordinates[19][2] =  0.5 * cl;
        mListOfCoordinates[20][0] = -0.5 * cl; mListOfCoordinates[20][1] =  0.0 * cl; mListOfCoordinates[20][2] = -0.5 * cl;
        mListOfCoordinates[21][0] = -0.5 * cl; mListOfCoordinates[21][1] = -0.5 * cl; mListOfCoordinates[21][2] =  0.5 * cl;
        mListOfCoordinates[22][0] = -0.5 * cl; mListOfCoordinates[22][1] = -0.5 * cl; mListOfCoordinates[22][2] =  0.0 * cl;
        mListOfCoordinates[23][0] = -0.5 * cl; mListOfCoordinates[23][1] = -0.5 * cl; mListOfCoordinates[23][2] = -0.5 * cl; 
        mListOfCoordinates[24][0] =  0.5 * cl; mListOfCoordinates[24][1] =  0.0 * cl; mListOfCoordinates[24][2] =  0.0 * cl;
        mListOfCoordinates[25][0] = -0.5 * cl; mListOfCoordinates[25][1] =  0.0 * cl; mListOfCoordinates[25][2] =  0.0 * cl; 
        
        double particle_density = this->SlowGetDensity();
        
        double cluster_volume = 1.0 * 1.0 * 1.0 * cl * cl * cl;
        
        double cluster_mass = particle_density * cluster_volume;
        
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = cluster_mass;
        
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = 0.16666666666 * cluster_mass * 1.0 * 1.0 * cl * cl;
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = 0.16666666666 * cluster_mass * 1.0 * 1.0 * cl * cl;
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = 0.16666666666 * cluster_mass * 1.0 * 1.0 * cl * cl;
         
        array_1d<double, 3> base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);     
  
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CubeCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CubeCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CubeCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CubeCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CubeCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CubeCluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CubeCluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void CubeCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CubeCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CubeCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void CubeCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double CubeCluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

