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
#include "linecluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    LineCluster3D::LineCluster3D() : Cluster3D()  {}
            
      
    LineCluster3D::LineCluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    LineCluster3D::LineCluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    LineCluster3D::LineCluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer LineCluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new LineCluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    LineCluster3D::~LineCluster3D() {}
      
    
    void LineCluster3D::CustomInitialize() {
        
        int number_of_spheres=10;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        mListOfRadii[0]= 0.1 * cl;
        mListOfRadii[1]= 0.1 * cl;
        mListOfRadii[2]= 0.1 * cl;   
        mListOfRadii[3]= 0.1 * cl;
        mListOfRadii[4]= 0.1 * cl;
        mListOfRadii[5]= 0.1 * cl;
        mListOfRadii[6]= 0.1 * cl;
        mListOfRadii[7]= 0.1 * cl;
        mListOfRadii[8]= 0.1 * cl;
        mListOfRadii[9]= 0.1 * cl;
        
        mListOfCoordinates[0][0] = -0.45 * cl; mListOfCoordinates[0][1] = 0.0; mListOfCoordinates[0][2] = 0.0;
        mListOfCoordinates[1][0] = -0.35 * cl; mListOfCoordinates[1][1] = 0.0; mListOfCoordinates[1][2] = 0.0;
        mListOfCoordinates[2][0] = -0.25 * cl; mListOfCoordinates[2][1] = 0.0; mListOfCoordinates[2][2] = 0.0;
        mListOfCoordinates[3][0] = -0.15 * cl; mListOfCoordinates[3][1] = 0.0; mListOfCoordinates[3][2] = 0.0; 
        mListOfCoordinates[4][0] = -0.05 * cl; mListOfCoordinates[4][1] = 0.0; mListOfCoordinates[4][2] = 0.0;  
        mListOfCoordinates[5][0] =  0.05 * cl; mListOfCoordinates[5][1] = 0.0; mListOfCoordinates[5][2] = 0.0;
        mListOfCoordinates[6][0] =  0.15 * cl; mListOfCoordinates[6][1] = 0.0; mListOfCoordinates[6][2] = 0.0;  
        mListOfCoordinates[7][0] =  0.25 * cl; mListOfCoordinates[7][1] = 0.0; mListOfCoordinates[7][2] = 0.0;
        mListOfCoordinates[8][0] =  0.35 * cl; mListOfCoordinates[8][1] = 0.0; mListOfCoordinates[8][2] = 0.0;  
        mListOfCoordinates[9][0] =  0.45 * cl; mListOfCoordinates[9][1] = 0.0; mListOfCoordinates[9][2] = 0.0;
        
        double particle_density = this->SlowGetDensity();
         
        double cluster_volume = 1.33333333333 * KRATOS_M_PI * 0.1 * 0.1 * 0.1 * cl * cl * cl
                                              + KRATOS_M_PI * 0.1 * 0.1 * 0.9 * cl * cl * cl;
        
        double cluster_mass = particle_density * cluster_volume;
        
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = cluster_mass;
        
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = 0.5 * cluster_mass * 0.1 * 0.1 * cl * cl;
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = 0.08333333333 * cluster_mass * (3.0 * 0.1 * 0.1 + 0.7 * 0.7) * cl * cl;
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = 0.08333333333 * cluster_mass * (3.0 * 0.1 * 0.1 + 0.7 * 0.7) * cl * cl;
         
        array_1d<double, 3> base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);     
  
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void LineCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void LineCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void LineCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double LineCluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

