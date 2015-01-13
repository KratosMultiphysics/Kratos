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
#include "cornkernelcluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    CornKernelCluster3D::CornKernelCluster3D() : Cluster3D() {}
            
      
    CornKernelCluster3D::CornKernelCluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    CornKernelCluster3D::CornKernelCluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    CornKernelCluster3D::CornKernelCluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer CornKernelCluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new CornKernelCluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    CornKernelCluster3D::~CornKernelCluster3D() {}
      
    
    void CornKernelCluster3D::CustomInitialize() {
        
        int number_of_spheres = 21;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
                
        for (int i = 0; i < 17; i++) { mListOfRadii[i]= 0.0012156; }
        
        //RADIUS 0.0012156
        mListOfCoordinates[ 0][0] =  0.00011844806399999969; mListOfCoordinates[ 0][1] =  0.00287423149800000035; mListOfCoordinates[ 0][2] = 0.0;
        mListOfCoordinates[ 1][0] = -0.00120891420000000025; mListOfCoordinates[ 1][1] =  0.00276718576200000065; mListOfCoordinates[ 1][2] = 0.0;
        mListOfCoordinates[ 2][0] =  0.00129596331599999885; mListOfCoordinates[ 2][1] =  0.00276718576200000065; mListOfCoordinates[ 2][2] = 0.0;
        mListOfCoordinates[ 3][0] = -0.00126624189600000015; mListOfCoordinates[ 3][1] =  0.00105445398600000030; mListOfCoordinates[ 3][2] = 0.0; 
        mListOfCoordinates[ 4][0] =  0.00168128420399999895; mListOfCoordinates[ 4][1] =  0.00135423310200000160; mListOfCoordinates[ 4][2] = 0.0;  
        mListOfCoordinates[ 5][0] =  0.00137918329200000035; mListOfCoordinates[ 5][1] = -0.00015627145799999956; mListOfCoordinates[ 5][2] = 0.0;
        mListOfCoordinates[ 6][0] = -0.00092971519199999995; mListOfCoordinates[ 6][1] = -0.00078205018199999955; mListOfCoordinates[ 6][2] = 0.0;  
        mListOfCoordinates[ 7][0] =  0.00014921489999999976; mListOfCoordinates[ 7][1] =  0.00085791577800000085; mListOfCoordinates[ 7][2] = 0.0;
        mListOfCoordinates[ 8][0] =  0.00077499362400000040; mListOfCoordinates[ 8][1] = -0.00170992981799999985; mListOfCoordinates[ 8][2] = 0.0;  
        mListOfCoordinates[ 9][0] = -0.00060383714400000026; mListOfCoordinates[ 9][1] = -0.00175962354599999995; mListOfCoordinates[ 9][2] = 0.0;
        mListOfCoordinates[10][0] = -0.00131780764800000050; mListOfCoordinates[10][1] =  0.00193457269799999995; mListOfCoordinates[10][2] = 0.0;
        mListOfCoordinates[11][0] = -0.00109543794000000000; mListOfCoordinates[11][1] =  0.00025407863400000091; mListOfCoordinates[11][2] = 0.0;
        mListOfCoordinates[12][0] = -0.00058999146000000033; mListOfCoordinates[12][1] =  0.00285779658599999990; mListOfCoordinates[12][2] = 0.0;
        mListOfCoordinates[13][0] =  0.00150612840000000080; mListOfCoordinates[13][1] =  0.00208500319800000035; mListOfCoordinates[13][2] = 0.0; 
        mListOfCoordinates[14][0] =  0.00156623982000000050; mListOfCoordinates[14][1] =  0.00065963926199999988; mListOfCoordinates[14][2] = 0.0;  
        mListOfCoordinates[15][0] =  0.00121268256000000020; mListOfCoordinates[15][1] = -0.00091943729399999931; mListOfCoordinates[15][2] = 0.0;
        mListOfCoordinates[16][0] =  0.00069562710000000035; mListOfCoordinates[16][1] =  0.00288153725399999985; mListOfCoordinates[16][2] = 0.0;  
        
        //RADIUS 0.0008104
        for (int i = 17; i < 20; i++) { mListOfRadii[i]= 0.0008104; }
        mListOfCoordinates[17][0] = -0.000382275; mListOfCoordinates[17][1] = -0.00274055; mListOfCoordinates[17][2] = 0.0;
        mListOfCoordinates[18][0] =  0.000492665; mListOfCoordinates[18][1] = -0.00266845; mListOfCoordinates[18][2] = 0.0;  
        mListOfCoordinates[19][0] = -0.000112185; mListOfCoordinates[19][1] = -0.00310410; mListOfCoordinates[19][2] = 0.0;
        
        //RADIUS 0.0003039
        for (int i = 20; i < 21; i++) { mListOfRadii[i]= 0.0003039; }
        mListOfCoordinates[20][0] = -0.000118455; mListOfCoordinates[20][1] = -0.00378590; mListOfCoordinates[20][2] = 0.0;
        
        double particle_density = this->SlowGetDensity();
        
        double a = 0.00531;
        double b = 0.00818;
        double c = 0.00251;
        
        double cluster_volume = a * b * c;
        
        double cluster_mass = particle_density * cluster_volume;
        
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = cluster_mass;
        
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = 0.08333333333 * cluster_mass * (b * b + c * c);
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = 0.08333333333 * cluster_mass * (a * a + c * c);
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = 0.08333333333 * cluster_mass * (a * a + b * b);
         
        array_1d<double, 3> base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
            
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void CornKernelCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void CornKernelCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double CornKernelCluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

