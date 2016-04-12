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
      
    
    void SoyBeanCluster3D::CustomInitialize(ProcessInfo& r_process_info) {
        
        int number_of_spheres = 13;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        /////////////////////PAPER INFO DATA TAKEN FROM *************************************** ADD!!!!!
        
        //Element built using 3ds max, there is nothing in GiD
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // Original maximum distance in the geometry is 6.55 m
        // It coincides with the main axis of the geometry
        // So, first of all, the geometry must be converted into unity length,
        // so we multiply by 1/6.55, obtaining 0.152671756
        // and we finally multiply that by the characteristic length given in the problem type
        
        double a, b, c; // Full ellipsoid axes
        
        cl *= 0.152671756;     
                
        a = 6.55 * 0.5 * cl;  b = 5.56 * 0.5 * cl;  c = 4.53 * 0.5 * cl; //Conversion to semi-axes
        
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
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[0]= 4.53 * cl;
        for (int i = 1; i < 13; i++) { mListOfRadii[i]= 3.53 * cl; }
                        
        //double particle_density = this->SlowGetDensity(); /////////////////////////////USE FAST
         
        //double cluster_volume = 0.3333333333 * 4.0 * KRATOS_M_PI * a * b * c; ////APPROXIMATE VOLUME, CALCULATE MORE EXACTLY
        
        //double cluster_mass = (this->SlowGetDensity()) * 0.3333333333 * 4.0 * KRATOS_M_PI * a * b * c;
        
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.3333333333 * 4.0 * KRATOS_M_PI * a * b * c;
        
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = 0.2 * cluster_mass * (b * b + c * c);
        base_principal_moments_of_inertia[1] = 0.2 * cluster_mass * (a * a + c * c);
        base_principal_moments_of_inertia[2] = 0.2 * cluster_mass * (a * a + b * b);
         
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void SoyBeanCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SoyBeanCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void SoyBeanCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double SoyBeanCluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

