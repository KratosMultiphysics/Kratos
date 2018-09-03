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
#include "custom_utilities/create_and_destroy.h"

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
    
    void SingleSphereCluster3D::Initialize(ProcessInfo& r_process_info) {
        
        RigidBodyElement3D::Initialize(r_process_info); //Skipping the Cluster3D level intentionally!! We should unifiy this when we code ECCENTRICITY and SCALES for Cluster3D

        int number_of_spheres = 1;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        mListOfRadii[0]= 0.5 * cl;
        
        double radius = mListOfRadii[0];
        
        double excentricity = GetProperties()[EXCENTRICITY];
        excentricity *= radius;
        double min_excentricity = 0.05 * excentricity; //TODO: this is a little bit arbitrary
        double max_excentricity = 1.00 * excentricity; //TODO: this is a little bit arbitrary
        double excentricity_std_deviation = GetProperties()[EXCENTRICITY_STANDARD_DEVIATION];
        std::string excentricity_distribution_type = GetProperties()[EXCENTRICITY_PROBABILITY_DISTRIBUTION];
        if (excentricity_distribution_type == "normal") excentricity = ParticleCreatorDestructor::rand_normal(excentricity, excentricity_std_deviation, max_excentricity, min_excentricity);
        else if (excentricity_distribution_type == "lognormal") excentricity = ParticleCreatorDestructor::rand_lognormal(excentricity, excentricity_std_deviation, max_excentricity, min_excentricity);
        else KRATOS_THROW_ERROR(std::runtime_error, "Unknown probability distribution.", "")
        
        mListOfCoordinates[0][0] = excentricity;
        mListOfCoordinates[0][1] = 0.0;
        mListOfCoordinates[0][2] = 0.0;
        
        double particle_density = this->SlowGetDensity();
         
        double cluster_volume = 4.0 * Globals::Pi / 3.0 * radius * radius * radius;
        
        double cluster_mass = particle_density * cluster_volume;
        
        double inertia_ball = 0.4 * cluster_mass * radius * radius;
        
        double steiner_excentricity_term = 2.0 * cluster_mass * excentricity * excentricity;
        
        array_1d<double,3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = inertia_ball;
        base_principal_moments_of_inertia[1] = inertia_ball + steiner_excentricity_term;
        base_principal_moments_of_inertia[2] = inertia_ball + steiner_excentricity_term;
        
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
    int SingleSphereCluster3D::SlowGetParticleMaterial()                                { return GetProperties()[PARTICLE_MATERIAL];}
    
}  // namespace Kratos

