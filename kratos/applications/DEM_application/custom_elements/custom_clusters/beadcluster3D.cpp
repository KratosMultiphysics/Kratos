//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva Latorre
//   Date:                $Date: April 2016

// System includes
#include <string>
#include <iostream>

// Project includes
#include "beadcluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    BeadCluster3D::BeadCluster3D() : Cluster3D()  {}
                  
    BeadCluster3D::BeadCluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
            
    BeadCluster3D::BeadCluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}
      
    BeadCluster3D::BeadCluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}
      
    Element::Pointer BeadCluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {          
        return Cluster3D::Pointer(new BeadCluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));        
    }      

    // Destructor
    BeadCluster3D::~BeadCluster3D() {}      
    
    void BeadCluster3D::CustomInitialize() {
        
        int number_of_spheres = 5;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        mListOfRadii[0]= 0.1 * cl;
        mListOfRadii[1]= 0.1 * cl;
        mListOfRadii[2]= 0.1 * cl;   
        mListOfRadii[3]= 0.1 * cl;
        mListOfRadii[4]= 0.1 * cl;
        
        mListOfCoordinates[0][0] = -0.4 * cl; mListOfCoordinates[0][1] = 0.0; mListOfCoordinates[0][2] = 0.0;
        mListOfCoordinates[1][0] = -0.2 * cl; mListOfCoordinates[1][1] = 0.0; mListOfCoordinates[1][2] = 0.0;
        mListOfCoordinates[2][0] =  0.0     ; mListOfCoordinates[2][1] = 0.0; mListOfCoordinates[2][2] = 0.0;
        mListOfCoordinates[3][0] =  0.2 * cl; mListOfCoordinates[3][1] = 0.0; mListOfCoordinates[3][2] = 0.0; 
        mListOfCoordinates[4][0] =  0.4 * cl; mListOfCoordinates[4][1] = 0.0; mListOfCoordinates[4][2] = 0.0;  
        
        double particle_density = this->SlowGetDensity();
        
        //TODO: Approximating the bead by the cylinder containing it
        double cluster_volume = KRATOS_M_PI * 0.1 * 0.1 * cl * cl;
        
        double cluster_mass = particle_density * cluster_volume;
        
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = cluster_mass;
        //TODO: Approximating the bead by the cylinder containing it
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = 0.5 * cluster_mass * 0.1 * 0.1 * cl * cl;
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = 0.08333333333 * cluster_mass * (3.0 * 0.1 * 0.1 + 0.7 * 0.7) * cl * cl;
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = 0.08333333333 * cluster_mass * (3.0 * 0.1 * 0.1 + 1.0) * cl * cl;
    }    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void BeadCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void BeadCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void BeadCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void BeadCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void BeadCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void BeadCluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void BeadCluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void BeadCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
        KRATOS_TRY
        KRATOS_CATCH("")
    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void BeadCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void   BeadCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info) {}
    void   BeadCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info) {}
    double BeadCluster3D::SlowGetDensity()                                      { return GetProperties()[PARTICLE_DENSITY];}

} // namespace Kratos

