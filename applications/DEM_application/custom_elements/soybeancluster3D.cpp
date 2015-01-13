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
#include "soybeancluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    SoyBeanCluster3D::SoyBeanCluster3D() : Cluster3D() {}
            
      
    SoyBeanCluster3D::SoyBeanCluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    SoyBeanCluster3D::SoyBeanCluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    SoyBeanCluster3D::SoyBeanCluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer SoyBeanCluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new SoyBeanCluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    SoyBeanCluster3D::~SoyBeanCluster3D() {}
      
    
    void SoyBeanCluster3D::CustomInitialize() {
        
        int number_of_spheres = 13;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        /////////////////////PAPER INFO DATA TAKEN FROM
        
        double a, b, c, /* FULL-AXES */ rf, resize_factor; //to convert millimeters to meters and divide by two to get semi-axes
        
        a = 6.55;  b = 5.56;  c = 4.53;  resize_factor = 0.001 * 0.5;  rf = resize_factor;
        
        mListOfCoordinates[ 0][0] =  0.0;  mListOfCoordinates[ 0][1] =  0.0;  mListOfCoordinates[ 0][2] =   0.0;
        mListOfCoordinates[ 1][0] =  3.0;  mListOfCoordinates[ 1][1] =  0.0;  mListOfCoordinates[ 1][2] =   0.0;
        mListOfCoordinates[ 2][0] = -3.0;  mListOfCoordinates[ 2][1] =  0.0;  mListOfCoordinates[ 2][2] =   0.0;
        mListOfCoordinates[ 3][0] =  0.0;  mListOfCoordinates[ 3][1] =  2.0;  mListOfCoordinates[ 3][2] =   0.0;
        mListOfCoordinates[ 4][0] =  0.0;  mListOfCoordinates[ 4][1] = -2.0;  mListOfCoordinates[ 4][2] =   0.0;
        mListOfCoordinates[ 5][0] =  2.2;  mListOfCoordinates[ 5][1] =  1.3;  mListOfCoordinates[ 5][2] =   0.0;
        mListOfCoordinates[ 6][0] =  2.2;  mListOfCoordinates[ 6][1] = -1.3;  mListOfCoordinates[ 6][2] =   0.0;
        mListOfCoordinates[ 7][0] = -2.2;  mListOfCoordinates[ 7][1] =  1.3;  mListOfCoordinates[ 7][2] =   0.0;
        mListOfCoordinates[ 8][0] = -2.2;  mListOfCoordinates[ 8][1] = -1.3;  mListOfCoordinates[ 8][2] =   0.0;
        mListOfCoordinates[ 9][0] =  1.7;  mListOfCoordinates[ 9][1] =  0.0;  mListOfCoordinates[ 9][2] =  0.75;
        mListOfCoordinates[10][0] =  1.7;  mListOfCoordinates[10][1] =  0.0;  mListOfCoordinates[10][2] = -0.75;
        mListOfCoordinates[11][0] = -1.7;  mListOfCoordinates[11][1] =  0.0;  mListOfCoordinates[11][2] =  0.75;
        mListOfCoordinates[12][0] = -1.7;  mListOfCoordinates[12][1] =  0.0;  mListOfCoordinates[12][2] = -0.75;
        
        for (int i = 1; i < 13; i++) { 
            mListOfCoordinates[i][0] *= rf;  mListOfCoordinates[i][1] *= rf;  mListOfCoordinates[i][2] *= rf;
        }
        
        mListOfRadii[0]= 4.53 * rf;
        for (int i = 1; i < 13; i++) { mListOfRadii[i]= 3.53 * rf; }
                        
        double particle_density = this->SlowGetDensity(); /////////////////////////////USE FAST
         
        double cluster_volume = 0.3333333333 * 4.0 * KRATOS_M_PI * a * b * c * rf * rf * rf; ////APPROXIMATE VOLUME, CALCULATE MORE EXACTLY
        
        double cluster_mass = particle_density * cluster_volume;
        
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = cluster_mass;
        
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = 0.2 * cluster_mass * (b * b * rf * rf + c * c * rf * rf);
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = 0.2 * cluster_mass * (a * a * rf * rf + c * c * rf * rf);
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = 0.2 * cluster_mass * (a * a * rf * rf + b * b * rf * rf);
         
        array_1d<double, 3> base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);  
  
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void SoyBeanCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void SoyBeanCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double SoyBeanCluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

