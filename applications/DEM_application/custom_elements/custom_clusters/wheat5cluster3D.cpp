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
#include "wheat5cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Wheat5Cluster3D::Wheat5Cluster3D() : Cluster3D() {}
            
      
    Wheat5Cluster3D::Wheat5Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Wheat5Cluster3D::Wheat5Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Wheat5Cluster3D::Wheat5Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Wheat5Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Wheat5Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Wheat5Cluster3D::~Wheat5Cluster3D() {}
      
    
    void Wheat5Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 5;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 6.164 (in meters) was the length of the particle in GiD (wheat5_01.gid file)
        // the inverse of that number is 0.162232317, that's the number which is multiplying cl
        // to obtain a maximum axis length of 1.0
        
        cl *= 0.162232317;
                
        mListOfRadii[0]= 1.715 * cl;
        mListOfRadii[1]= 1.595 * cl;
        mListOfRadii[2]= 1.595 * cl;
        mListOfRadii[3]= 1.252 * cl;
        mListOfRadii[4]= 1.252 * cl;
        
        mListOfCoordinates[0][0] =  0.0      ;  mListOfCoordinates[0][1] =  0.0;  mListOfCoordinates[0][2] = 0.0;
        mListOfCoordinates[1][0] =  1.0  * cl;  mListOfCoordinates[1][1] =  0.0;  mListOfCoordinates[1][2] = 0.0;
        mListOfCoordinates[2][0] = -1.0  * cl;  mListOfCoordinates[2][1] =  0.0;  mListOfCoordinates[2][2] = 0.0; 
        mListOfCoordinates[3][0] =  1.83 * cl;  mListOfCoordinates[3][1] =  0.0;  mListOfCoordinates[3][2] = 0.0;
        mListOfCoordinates[4][0] = -1.83 * cl;  mListOfCoordinates[4][1] =  0.0;  mListOfCoordinates[4][2] = 0.0;
        
        double particle_density = this->SlowGetDensity();
        
        double a = 3.120 * cl;
        double b = 1.825 * cl;
        double c = 1.715 * cl;
        
        double cluster_volume = 1.33333333333333333 * KRATOS_M_PI * a * b * c;
        
        double cluster_mass = particle_density * cluster_volume;
        
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = cluster_mass;
        
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = 0.2 * cluster_mass * (b * b + c * c);
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = 0.2 * cluster_mass * (a * a + c * c);
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = 0.2 * cluster_mass * (a * a + b * b);
         
        array_1d<double, 3> base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
            
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Wheat5Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Wheat5Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Wheat5Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Wheat5Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Wheat5Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Wheat5Cluster3D::InitializeSolutionStep(const ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Wheat5Cluster3D::FinalizeSolutionStep(const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Wheat5Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Wheat5Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Wheat5Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void Wheat5Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double Wheat5Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

