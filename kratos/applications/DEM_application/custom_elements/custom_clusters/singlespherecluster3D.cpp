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
#include <cmath>

// Project includes
#include "singlespherecluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    SingleSphereCluster3D::SingleSphereCluster3D() : Cluster3D() {}
                  
    SingleSphereCluster3D::SingleSphereCluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
            
    SingleSphereCluster3D::SingleSphereCluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}
      
    SingleSphereCluster3D::SingleSphereCluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}
      
    Element::Pointer SingleSphereCluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {          
        return Cluster3D::Pointer(new SingleSphereCluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));        
    }      

    // Destructor
    SingleSphereCluster3D::~SingleSphereCluster3D() {}      
    
    void SingleSphereCluster3D::CustomInitialize(ProcessInfo& r_process_info) {
        
        int number_of_spheres = 1;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        mListOfRadii[0]= 0.5 * cl;
        
        double radius = mListOfRadii[0];
        
        double excentricity = GetGeometry()[0].FastGetSolutionStepValue(EXCENTRICITY);
        excentricity *= cl;
        
        mListOfCoordinates[0][0] = excentricity;
        mListOfCoordinates[0][1] = 0.0;
        mListOfCoordinates[0][2] = 0.0;
        
        double particle_density = this->SlowGetDensity();
         
        double cluster_volume = 4.0 * KRATOS_M_PI_3 * radius * radius * radius;
        
        double cluster_mass = particle_density * cluster_volume;
        
        double inertia_ball = 0.4 * cluster_mass * radius * radius;
        
        double steiner_excentricity_term = 2.0 * cluster_mass * excentricity * excentricity;
        
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = cluster_mass;        
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = inertia_ball;
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = inertia_ball + steiner_excentricity_term;
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = inertia_ball + steiner_excentricity_term;        
    }    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SingleSphereCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SingleSphereCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SingleSphereCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SingleSphereCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SingleSphereCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SingleSphereCluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SingleSphereCluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void SingleSphereCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
        KRATOS_TRY
        KRATOS_CATCH("")
    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SingleSphereCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void   SingleSphereCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info) {}
    void   SingleSphereCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info) {}
    double SingleSphereCluster3D::SlowGetDensity()                                      { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos

