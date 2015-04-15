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
#include "ballast2cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Ballast2Cluster3D::Ballast2Cluster3D() : Cluster3D() {}
            
      
    Ballast2Cluster3D::Ballast2Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Ballast2Cluster3D::Ballast2Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Ballast2Cluster3D::Ballast2Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Ballast2Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Ballast2Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Ballast2Cluster3D::~Ballast2Cluster3D() {}
      
    
    void Ballast2Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 28;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 1.241 (in meters) was the medium diameter of the rock in GiD (Rock3_01.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.8058,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.8058;
     
        mListOfCoordinates [ 0][0] =-0.112912; mListOfCoordinates [ 0][1] =-0.141381; mListOfCoordinates [ 0][2] =-0.086692;
        mListOfCoordinates [ 1][0] =-0.263107; mListOfCoordinates [ 1][1] =-0.391062; mListOfCoordinates [ 1][2] =-0.015769;
        mListOfCoordinates [ 2][0] =-0.101771; mListOfCoordinates [ 2][1] =-0.144973; mListOfCoordinates [ 2][2] =-0.037814;
        mListOfCoordinates [ 3][0] =-0.074754; mListOfCoordinates [ 3][1] =-0.070499; mListOfCoordinates [ 3][2] =-0.126463;
        mListOfCoordinates [ 4][0] =-0.089507; mListOfCoordinates [ 4][1] =-0.406942; mListOfCoordinates [ 4][2] =-0.195216;
        mListOfCoordinates [ 5][0] =-0.226924; mListOfCoordinates [ 5][1] =-0.134188; mListOfCoordinates [ 5][2] =-0.059330;
        mListOfCoordinates [ 6][0] =-0.105466; mListOfCoordinates [ 6][1] =-0.288940; mListOfCoordinates [ 6][2] =-0.067559;
        mListOfCoordinates [ 7][0] =-0.090178; mListOfCoordinates [ 7][1] =-0.242862; mListOfCoordinates [ 7][2] =-0.063233;
        mListOfCoordinates [ 8][0] =-0.256418; mListOfCoordinates [ 8][1] =-0.367489; mListOfCoordinates [ 8][2] =-0.006342;
        mListOfCoordinates [ 9][0] =-0.158169; mListOfCoordinates [ 9][1] =-0.418187; mListOfCoordinates [ 9][2] =-0.101108;
        mListOfCoordinates [10][0] =-0.048990; mListOfCoordinates [10][1] =-0.390773; mListOfCoordinates [10][2] =-0.065806;
        mListOfCoordinates [11][0] =-0.354623; mListOfCoordinates [11][1] =-0.320103; mListOfCoordinates [11][2] =-0.091540;
        mListOfCoordinates [12][0] =-0.194634; mListOfCoordinates [12][1] =-0.314234; mListOfCoordinates [12][2] =-0.175218;
        mListOfCoordinates [13][0] =-0.089719; mListOfCoordinates [13][1] =-0.155954; mListOfCoordinates [13][2] =-0.021529;
        mListOfCoordinates [14][0] =-0.112499; mListOfCoordinates [14][1] =-0.133166; mListOfCoordinates [14][2] =-0.125492;
        mListOfCoordinates [15][0] =-0.129450; mListOfCoordinates [15][1] =-0.126432; mListOfCoordinates [15][2] =-0.151270;
        mListOfCoordinates [16][0] =-0.018014; mListOfCoordinates [16][1] =-0.261887; mListOfCoordinates [16][2] =-0.286450;
        mListOfCoordinates [17][0] =-0.179901; mListOfCoordinates [17][1] =-0.392846; mListOfCoordinates [17][2] =-0.030392;
        mListOfCoordinates [18][0] =-0.035194; mListOfCoordinates [18][1] =-0.273314; mListOfCoordinates [18][2] =-0.238892;
        mListOfCoordinates [19][0] =-0.284028; mListOfCoordinates [19][1] =-0.233186; mListOfCoordinates [19][2] =-0.165200;
        mListOfCoordinates [20][0] =-0.186650; mListOfCoordinates [20][1] =-0.272654; mListOfCoordinates [20][2] =-0.158157;
        mListOfCoordinates [21][0] =-0.206529; mListOfCoordinates [21][1] =-0.211570; mListOfCoordinates [21][2] =-0.018967;
        mListOfCoordinates [22][0] =-0.031457; mListOfCoordinates [22][1] =-0.185757; mListOfCoordinates [22][2] =-0.003505;
        mListOfCoordinates [23][0] =-0.033867; mListOfCoordinates [23][1] =-0.128483; mListOfCoordinates [23][2] =-0.137474;
        mListOfCoordinates [24][0] =-0.203659; mListOfCoordinates [24][1] =-0.066346; mListOfCoordinates [24][2] =-0.007059;
        mListOfCoordinates [25][0] =-0.250852; mListOfCoordinates [25][1] =-0.057671; mListOfCoordinates [25][2] =-0.017244;
        mListOfCoordinates [26][0] =-0.284695; mListOfCoordinates [26][1] =-0.174947; mListOfCoordinates [26][2] =-0.123858;
        mListOfCoordinates [27][0] =-0.151848; mListOfCoordinates [27][1] =-0.096740; mListOfCoordinates [27][2] =-0.001015;

        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0] = 0.339049;
        mListOfRadii[ 1] = 0.233528;
        mListOfRadii[ 2] = 0.388827;
        mListOfRadii[ 3] = 0.331505;
        mListOfRadii[ 4] = 0.227643;
        mListOfRadii[ 5] = 0.374958;
        mListOfRadii[ 6] = 0.262937;
        mListOfRadii[ 7] = 0.291307;
        mListOfRadii[ 8] = 0.245471;
        mListOfRadii[ 9] = 0.253552;
        mListOfRadii[10] = 0.265224;
        mListOfRadii[11] = 0.269778;
        mListOfRadii[12] = 0.288632;
        mListOfRadii[13] = 0.380106;
        mListOfRadii[14] = 0.300409;
        mListOfRadii[15] = 0.333144;
        mListOfRadii[16] = 0.110335;
        mListOfRadii[17] = 0.268086;
        mListOfRadii[18] = 0.202354;
        mListOfRadii[19] = 0.293389;
        mListOfRadii[20] = 0.284530;
        mListOfRadii[21] = 0.321070;
        mListOfRadii[22] = 0.365016;
        mListOfRadii[23] = 0.329093;
        mListOfRadii[24] = 0.345479;
        mListOfRadii[25] = 0.319487;
        mListOfRadii[26] = 0.268785;
        mListOfRadii[27] = 0.387152;
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.632108 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.130002031;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.108529343;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.157003505; 
  
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast2Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void Ballast2Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double Ballast2Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

