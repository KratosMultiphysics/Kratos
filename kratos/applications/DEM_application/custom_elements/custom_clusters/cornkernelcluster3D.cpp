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
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);        
        
        // Original maximum distance in the geometry in corn_kernel_design_05.gid is 0.00818309 m
        // It coincides with the main axis of the geometry, the almost-symmetrical one
        // So, first of all, the geometry must be converted into unity length,
        // so we multiply by 1/0.00818309, obtaining 122.20322641
        // and we finally multiply that by the characteristic length given in the problem type 
        
        cl *= 122.20322641;
        
        for (int i = 0; i < 17; i++) { mListOfRadii[i]= 0.0012156 * cl; }
        
        //RADIUS 0.0012156
        mListOfCoordinates[ 0][0] =  0.000118448063999999 * cl; mListOfCoordinates[ 0][1] =  0.002874231498000000 * cl; mListOfCoordinates[ 0][2] = 0.0;
        mListOfCoordinates[ 1][0] = -0.001208914200000000 * cl; mListOfCoordinates[ 1][1] =  0.002767185762000000 * cl; mListOfCoordinates[ 1][2] = 0.0;
        mListOfCoordinates[ 2][0] =  0.001295963315999998 * cl; mListOfCoordinates[ 2][1] =  0.002767185762000000 * cl; mListOfCoordinates[ 2][2] = 0.0;
        mListOfCoordinates[ 3][0] = -0.001266241896000000 * cl; mListOfCoordinates[ 3][1] =  0.001054453986000000 * cl; mListOfCoordinates[ 3][2] = 0.0; 
        mListOfCoordinates[ 4][0] =  0.001681284203999998 * cl; mListOfCoordinates[ 4][1] =  0.001354233102000001 * cl; mListOfCoordinates[ 4][2] = 0.0;  
        mListOfCoordinates[ 5][0] =  0.001379183292000000 * cl; mListOfCoordinates[ 5][1] = -0.000156271457999999 * cl; mListOfCoordinates[ 5][2] = 0.0;
        mListOfCoordinates[ 6][0] = -0.000929715191999999 * cl; mListOfCoordinates[ 6][1] = -0.000782050181999999 * cl; mListOfCoordinates[ 6][2] = 0.0;  
        mListOfCoordinates[ 7][0] =  0.000149214899999999 * cl; mListOfCoordinates[ 7][1] =  0.000857915778000000 * cl; mListOfCoordinates[ 7][2] = 0.0;
        mListOfCoordinates[ 8][0] =  0.000774993624000000 * cl; mListOfCoordinates[ 8][1] = -0.001709929817999999 * cl; mListOfCoordinates[ 8][2] = 0.0;  
        mListOfCoordinates[ 9][0] = -0.000603837144000000 * cl; mListOfCoordinates[ 9][1] = -0.001759623545999999 * cl; mListOfCoordinates[ 9][2] = 0.0;
        mListOfCoordinates[10][0] = -0.001317807648000000 * cl; mListOfCoordinates[10][1] =  0.001934572697999999 * cl; mListOfCoordinates[10][2] = 0.0;
        mListOfCoordinates[11][0] = -0.001095437940000000 * cl; mListOfCoordinates[11][1] =  0.000254078634000000 * cl; mListOfCoordinates[11][2] = 0.0;
        mListOfCoordinates[12][0] = -0.000589991460000000 * cl; mListOfCoordinates[12][1] =  0.002857796585999999 * cl; mListOfCoordinates[12][2] = 0.0;
        mListOfCoordinates[13][0] =  0.001506128400000000 * cl; mListOfCoordinates[13][1] =  0.002085003198000000 * cl; mListOfCoordinates[13][2] = 0.0; 
        mListOfCoordinates[14][0] =  0.001566239820000000 * cl; mListOfCoordinates[14][1] =  0.000659639261999999 * cl; mListOfCoordinates[14][2] = 0.0;  
        mListOfCoordinates[15][0] =  0.001212682560000000 * cl; mListOfCoordinates[15][1] = -0.000919437293999999 * cl; mListOfCoordinates[15][2] = 0.0;
        mListOfCoordinates[16][0] =  0.000695627100000000 * cl; mListOfCoordinates[16][1] =  0.002881537253999999 * cl; mListOfCoordinates[16][2] = 0.0;  
        
        //RADIUS 0.0008104
        for (int i = 17; i < 20; i++) { mListOfRadii[i]= 0.0008104 * cl; }
        mListOfCoordinates[17][0] = -0.000382275 * cl; mListOfCoordinates[17][1] = -0.00274055 * cl; mListOfCoordinates[17][2] = 0.0;
        mListOfCoordinates[18][0] =  0.000492665 * cl; mListOfCoordinates[18][1] = -0.00266845 * cl; mListOfCoordinates[18][2] = 0.0;  
        mListOfCoordinates[19][0] = -0.000112185 * cl; mListOfCoordinates[19][1] = -0.00310410 * cl; mListOfCoordinates[19][2] = 0.0;
        
        //RADIUS 0.0003039
        for (int i = 20; i < 21; i++) { mListOfRadii[i]= 0.0003039 * cl; }
        mListOfCoordinates[20][0] = -0.000118455 * cl; mListOfCoordinates[20][1] = -0.00378590 * cl; mListOfCoordinates[20][2] = 0.0;
        
        //double particle_density = this->SlowGetDensity();
        
        double a = 0.00531 * cl;
        double b = 0.00818 * cl;
        double c = 0.00251 * cl;
        
        //double cluster_volume = a * b * c;
        
        //double cluster_mass = particle_density * cluster_volume;
        
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * a * b * c;
        
        array_1d<double,3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = 0.08333333333 * cluster_mass * (b * b + c * c);
        base_principal_moments_of_inertia[1] = 0.08333333333 * cluster_mass * (a * a + c * c);
        base_principal_moments_of_inertia[2] = 0.08333333333 * cluster_mass * (a * a + b * b);
            
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void CornKernelCluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void CornKernelCluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void CornKernelCluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double CornKernelCluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

