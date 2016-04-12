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
#include "corn3cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Corn3Cluster3D::Corn3Cluster3D() : Cluster3D() {}
            
      
    Corn3Cluster3D::Corn3Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Corn3Cluster3D::Corn3Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Corn3Cluster3D::Corn3Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Corn3Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Corn3Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Corn3Cluster3D::~Corn3Cluster3D() {}
      
    
    void Corn3Cluster3D::CustomInitialize(ProcessInfo& r_process_info) {
        
        int number_of_spheres = 3;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 6.53 (in meters) was the length of the particle in GiD (corn3_design_01.gid file)
        // the inverse of that number is 0.15314, that's the number which is multiplying cl
        // to obtain a maximum axis length of 1.0
        
        cl *= 0.15314;
                
        mListOfRadii[0]= 2.340 * cl;
        mListOfRadii[1]= 1.872 * cl;
        mListOfRadii[2]= 0.732 * cl;
        
        mListOfCoordinates[0][0] = 0.0; mListOfCoordinates[0][1] =  0.815 * cl; mListOfCoordinates[0][2] = 0.0;
        mListOfCoordinates[1][0] = 0.0; mListOfCoordinates[1][1] = -0.870 * cl; mListOfCoordinates[1][2] = 0.0;
        mListOfCoordinates[2][0] = 0.0; mListOfCoordinates[2][1] = -2.643 * cl; mListOfCoordinates[2][2] = 0.0; 
        
        //double particle_density = this->SlowGetDensity();
        
        //For the time being, we assume, in order to obtain the inertias, an equivalent sphere with an average diameter
        //'a' is the equivalent radius
        
        double a = 2.742 * cl;
                
        //double cluster_volume = 1.33333333333333333 * KRATOS_M_PI * a * a * a;
                
        //double cluster_mass = particle_density * cluster_volume;
        
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 1.33333333333333333 * KRATOS_M_PI * a * a * a;
        
        double inertia = 0.4 * cluster_mass * a * a;
        
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = base_principal_moments_of_inertia[1] = base_principal_moments_of_inertia[2] = inertia;
        
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Corn3Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Corn3Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Corn3Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Corn3Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Corn3Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Corn3Cluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Corn3Cluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Corn3Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Corn3Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Corn3Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void Corn3Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double Corn3Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

